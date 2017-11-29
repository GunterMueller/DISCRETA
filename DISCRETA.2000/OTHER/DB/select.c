/* ----------------------------------------------------------------------- */
/* SELECT.C                                                                */
/*                                                                         */
/* In dieser Datei sind die Funktionen zum Suchen alle Datenseiten die     */
/* nichtleeren Schnitt mit dem Suchraum haben, sowie zum Durchsuchen der   */
/* Datenseiten nach den Datensaetzen innerhalb der Suchgrenzen enthalten.   */
/*                                                                         */
/* Funktionen                                                              */
/* dbms_search                                                             */
/* next_pages                                                              */
/* search_data                                                             */
/* check_and_out                                                           */
/* search_overlay                                                          */
/* push_to_plist                                                           */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
/* #include <conio.h> */

#include "dbmsdef.h"        /* Definitionen                                */
#include "dbms.h"           /* Benutzerstrukturen/operationen              */
#include "dbmstype.h"       /* Interne Verwaltungsstrukturen               */
#include "dbmsprot.h"       /* Prototypen der internen Funktionen          */

/* ----------------------------------------------------------------------- */
/* Modulinterne Funktionen                                                 */
/* ----------------------------------------------------------------------- */
int next_pages( TABLE *, PLIST *, BOUND * );
int search_data( TABLE *, USER_RIGHTS * );
int check_and_out( DINFO *di, USER_RIGHTS * );
int search_overlay( UCHAR *i_data, ULONG offset );
int push_to_plist( PLIST *, ULONG );
int defaultout( char * );

/* ----------------------------------------------------------------------- */
/* Makros fuer Stackhandling                                                */
/* ----------------------------------------------------------------------- */
#define push( x )   {*(stack+spos)=(x);spos++;}
#define pop         *(stack+spos-1);spos--;
#define stack_empty ((spos>0)?(0):(1))

extern TABLE      table;
BASE_PAGE        *base;
extern INDEX      index_db[];
extern UINT       index_count;

/* ----------------------------------------------------------------------- */
/* Makro fuer Projektion eines Datensatzes auf die i-te Komponente          */
/* ----------------------------------------------------------------------- */
#define projection(d,i) ((base->variable)?((d)+*(((UINT*)(d))+i+1)):((d)+(index_db[i].position)))

char   *databuffer = NULL;
char   *pagebuffer = NULL;
int   (*outfunc)( char * ) = NULL;
BOUND   bound;
ULONG   counter;

/* ----------------------------------------------------------------------- */
/* Name       dbms_select                                                  */
/* Definition int dbms_select( SEARCH_PARM *sp );                          */
/* Prototyp   dbms.h                                                       */
/* Funktion   Hauptfunktion von SEARCH.C. Initialisierung der Tabelle,     */
/*            Pruefung der Zugriffsberechtigung, Einrichten der Seitenpuffer*/
/*            Suchoperationen beginnen immer mit der Wurzelseite. Dann     */
/*            werden alle Wege innerhalb der Suchgrenzen weiterverfolgt    */
/*            bis zu den Datenseiten. Zur Bearbeitung wird ein Stack       */
/*            fuer die Seitenadressen benoetigt.Das Durchsuchen der Daten-   */
/*            seiten erledigt die Funktion search_data                     */
/*            Im Suchparameter kann di Funktion zur Ausgabe der Daten      */
/*            frei bestimmt werden. Somit kann die Ausgabe an die          */
/*            Beduerfnisse des Benutzers angepasst werden.                  */
/* Benutzte Funktionen                                                     */
/*            init_table        (INIT.C)                                   */
/*            memory_alloc      (MEMORY.C)                                 */
/*            memory_free       (MEMORY.C)                                 */
/*            read_page         (FILE_OP.C)                                */
/*            next_pages        (SEARCHH.C)                                */
/*            search_data       (SEARCHH.C)                                */
/* ----------------------------------------------------------------------- */
int dbms_select( SELECT_PARM *sp ) {

   TABLE       *tbl = NULL;
   PLIST        plist;
   USER_RIGHTS  usr;
   int          error = DBMS_ALL_OK;

   if (LogFile!=NULL) fprintf(LogFile,"dbms_select\n");

#ifdef _M11MU_
   clock_t t0;
#endif

   error = init_table( &tbl, sp->tblname );
   if ( error ) {
      printf("-- initialization of table failed\n");
      return error;
   }

   if (LogFile!=NULL) fprintf(LogFile,"init ok\n");

   base = (BASE_PAGE*)(tbl->base);

   if ( tbl->read_access == '-' ) {  /* no search operations allowed */
      error = DBMS_ACCESS_DENIED;
      drop_table( tbl );
      printf("-- access denied\n");
      return error;
   }

#ifdef _M11USERCHECK_                             /* if user check is wanted */
   error = get_user_rights( tbl, id, &usr );         /* the users permission on */
   if ( ! error ) {                                  /* this table is recalled  */
      if ( usr.search < (base->tflag).search_limit){ /* by get_user_rights and  */
         error = DBMS_ACCESS_DENIED;                 /* then tables limitation  */
      }                                              /* is proved               */
   }

   if ( error ) {
      drop_table( tbl );
      return error;
   }

   t0 = clock();                               /* search operations    */
   while ( tbl->usage < 0 ) {                  /* can be started only  */
      if ( clock()-t0 > sp->timeout*CLK_TCK ){ /* if no change ops are */
         drop_table( tbl );                    /* active               */
         return( DBMS_TIMEOUT );
      }
   }
   (tbl->usage)++;
#endif

   bound.mode = sp->mode;     /* set boundary    */
   bound.low  = sp->low;      /* according to sp */
   bound.up   = sp->up;
   counter    = 0;

   if ( sp->outfunc != NULL ) {   /* select output function. if NULL */
      outfunc = sp->outfunc;      /* is set in sp, default output is */
   } else {                       /* used, meaning that no output is */
      outfunc = defaultout;       /* made, records are counted only  */
   }

   /* -------------------------------------------------------------------- */
   /* Seitenpuffer fuer eine Seite, Datensatzpuffer und PLIST-Stack anlegen */
   /* -------------------------------------------------------------------- */
   if ( base->max_data_length == 0 ) {
      databuffer = (char*)memory_alloc(base->page_size);
   } else {
      databuffer = (char*)memory_alloc(base->max_data_length);
   }
   pagebuffer = (char*)memory_alloc(base->page_size);
   plist.max  = PLISTCNT;
   plist.list = (ULONG*)memory_alloc(plist.max * sizeof(ULONG));

   if ( pagebuffer == NULL || plist.list == NULL || databuffer == NULL ) {

      if ( pagebuffer != NULL ) memory_free( pagebuffer );
      if ( databuffer != NULL ) memory_free( databuffer );
      if ( plist.list != NULL ) memory_free( plist.list );
      error = DBMS_MEMORY_ERROR;
      drop_table( tbl );
      printf("-- memory failure\n");
      return error;

   } else {

      if (LogFile!=NULL) fprintf(LogFile,"start reading pages\n");

      tbl->page_buf = (UCHAR*)pagebuffer;
      tbl->data_buf = (UCHAR*)databuffer;

      plist.list[0]  = base->root_page;
      plist.act      = 1;

      /* ----------------------------------------------------------------- */
      /* Hier beginnt nun der eigentliche Durchlauf. Nach der              */
      /* Initialisierung befindet sich die Wurzelseite auf dem Stack. Die  */
      /* Schleife wird solange durchlaufen, wie sich Seiten auf diesem     */
      /* befinden. Wenn noch eine Seite auf dem Stack ist, so wird diese   */
      /* gelesen, und dann, falls sie eine    Indexseite ist nach Nach-    */
      /* folgeseiten innerhalb der Suchgrenzen abgesucht (next_pages) oder */
      /* nach Datensaetzen innerhalb der Grenzen, die dann ausgegeben      */
      /* werden (search_data).                                             */
      /* ----------------------------------------------------------------- */

      while ( ! error && plist.act > 0 ) {

         (plist.act)--;
         tbl->active = *((plist.list) + plist.act);
         if (LogFile!=NULL) fprintf(LogFile,"read page %ld\n",tbl->active);
         error = read_page( &table );

         if ( error ) break;

         if ( ((PTYPE*)(pagebuffer))->type == PTYPEDATA ) {
            if (LogFile!=NULL) fprintf(LogFile,"search data\n");
            error = search_data( tbl, &usr );
         } else {
            if (LogFile!=NULL) fprintf(LogFile,"next pages\n");
            error = next_pages( tbl, &plist, &bound );
         }
      }

      memory_free( plist.list );
      memory_free( pagebuffer );
      memory_free( databuffer );
   }

   sp->count = counter;
   drop_table( &table );

   return( error );
}

/* ----------------------------------------------------------------------- */
/* Name       next_pages                                                   */
/* Definition int next_pages( void );                                      */
/* Prototyp   module internal                                              */
/* Funktion   This function searches on index pages (pagebuffer) for all   */
/*            offsets of following pages lying in the search boundary.     */
/*            Therefore a stack for the inode field indices is used. The   */
/*            inode array on the index page is passed in the following     */
/*            order: root - smaller - greater. All page offsets found are  */
/*            pushed to plist stack by push_to_plist.                      */
/* Benutzte Funktionen                                                     */
/*         push_to_plist (module internal function)                        */
/*         push          (module internal macro)                           */
/*         pop           (module internal macro)                           */
/*         stack_empty   (module internal macro)                           */
/*         projection                                                      */
/* ----------------------------------------------------------------------- */
int next_pages( TABLE *tbl, PLIST *plst, BOUND *bnd ) {

   IPAGE *ipg   = (IPAGE*)(tbl->page_buf + sizeof(PTYPE) );
   INODE *ina   = (INODE*)(tbl->page_buf + sizeof(PTYPE) + sizeof(IPAGE));
   ULONG *list  = (ULONG*)(tbl->page_buf + ipg->olist );
   BPOS  *bpos  = (BPOS *)(tbl->page_buf + ipg->bpos );
   INODE *inode;
   UINT   i;                   /* Hilfsvariablen */
   char   okflag;
   UCHAR  *cvn, *cvl, *cvu;        /* Zeiger auf Vergleichswerte */
   SINT   cmpres;
   UINT  *stack;               /* Feld fuer Knotenstack */
   UINT   spos = 0;            /* Akt. Index fuer Knotenstack */
   int    error = DBMS_ALL_OK;

   if ( ipg->node_c == 0 ) {              /* page contains no nodes, means */
      error=push_to_plist( plst, *list ); /* that there is only one page   */
   } else {                               /* following                     */

      stack  = (UINT*)memory_alloc( ipg->node_c * sizeof(UINT));

      if ( stack == NULL ) {
         error = DBMS_MEMORY_ERROR;
      } else {

         push( 0 );       /* start at root node */

         while( ! stack_empty ) {

            inode = ina + pop;

            i = (bpos + inode->bnum)->comp;

            cvn = tbl->page_buf + (bpos+inode->bnum)->pos;

            /* ----------------------------------------------------------- */
            /* If the lower bound is smaller than the nodes compare value  */
            /* the following smaller node is to be examined too. If the    */
            /* smaller pointer of the active node is extern, push the      */
            /* offset of the page to plist stack, otherwise push the nodes */
            /* field index to stack for further examination. Do the same   */
            /* procedure if lower bound is equal to the nodes value and    */
            /* equality is on smaller side.                                */
            /* ----------------------------------------------------------- */
            okflag = 1;
            if ( ! ((bnd->mode)[i] & LBUL) ) {
               cvl = (UCHAR*)proj( tbl, i, bnd->low );
               cmpres = (index_db[i].cmpfunc)( cvl, cvn, index_db[i].length );
               if ( cmpres > 0 ) {
                  okflag = 0;
               } else if ( cmpres == 0 ) {
                  if ( ! (inode->flag.equal_smaller && ((bnd->mode)[i]&LBEQ))) {
                     okflag = 0;
                  }
               }
            }
            if ( okflag ) {
               if ( inode->flag.smaller_extern ) {   /* markiere zugehoerige Seite */
                  if ( *(list+inode->smaller) > 0 ) {
                     error = push_to_plist( plst, *(list+inode->smaller) );
                     *(list+inode->smaller) = 0;
                  }
               } else {
                  push( inode->smaller );    /* bzw. lege Knoten auf Stack */
               }
            }

            /* ----------------------------------------------------------- */
            /* If the upper bound is greater than the nodes compare value  */
            /* the following greater node is to be examined too. If the    */
            /* greater pointer of the active node is extern, push the      */
            /* offset of the page to plist stack, otherwise push the nodes */
            /* field index to stack for further examination. Do the same   */
            /* procedure if upper bound is equal to the nodes value and    */
            /* equality is on greater side.                                */
            /* ----------------------------------------------------------- */
            okflag = 1;
            if ( ! ((bnd->mode)[i] & UBUL) ) {
               cvu = (UCHAR*)proj( tbl, i, bnd->up );
               cmpres = (index_db[i].cmpfunc)( cvu, cvn, index_db[i].length );
               if ( cmpres < 0 ) {
                  okflag = 0;
               } else if ( cmpres== 0 ) {
                  if ( inode->flag.equal_smaller || ! ((bnd->mode)[i]&UBEQ) ) {
                     okflag = 0;
                  }
               }
            }
            if ( okflag ) {
               if ( inode->flag.greater_extern ) {   /* markiere zugehoerige Seite */
                  if ( *(list+inode->greater) > 0 ) {
                     error = push_to_plist( plst, *(list+inode->greater) );
                     *(list+inode->greater) = 0;
                  }
               } else {
                  push( inode->greater );    /* bzw. lege Knoten auf Stack */
               }
            }
         }

         memory_free( stack );
      }
   }

   return( error );
}

/* ----------------------------------------------------------------------- */
/* Name       search_data                                                  */
/* Definition int search_data( void );                                     */
/* Prototyp   Modulinterrn                                                 */
/* Funktion   Die Funktion sucht auf der im Seitenpuffer der Tabelle       */
/*            stehenden Datenseite nach den Datensaetzen, die innerhalb    */
/*            der von bound gesetzten Grenzen liegen. Wenn ein Datensatz   */
/*            gefunden wurde, der diese Bedingung erfuellt, so gibt es     */
/*            zwei Moeglichkeiten: 1) er ist ein Verweis auf eine Daten-   */
/*            seite, dann wird die Suche mittels search_overlay weiter-    */
/*            gefuehrt, 2) er ist ein 'normaler' Datensatz, dann wird er    */
/*            mittels der Ausgabefunktion der Tabelle ausgegeben.          */
/* ----------------------------------------------------------------------- */
int search_data( TABLE *tbl, USER_RIGHTS *usr ) {

   DPAGE  *dpg = (DPAGE*)(pagebuffer + sizeof(PTYPE));
   DNODE  *dna = (DNODE*)(pagebuffer + sizeof(PTYPE) + sizeof(DPAGE));
   DINFO  *dia = (DINFO*)(pagebuffer + dpg->data);
   DNODE  *dnode;
   UINT    i;
   UCHAR  *cvn, *cvl, *cvu;       /* Zeiger auf Vergleichswerte */
   SINT    cmpres;
   char    okflag;
   UINT   *stack;                 /* Stack fuer Knotendurchlauf */
   UINT    spos = 0;
   int     error = DBMS_ALL_OK;

   if ( dpg->dic > 0 ) {

      if ( dpg->dnc == 0 ) {
         error = check_and_out( dia, usr ); /* no nodes => only one record on page */
      } else {

         stack = (UINT*)memory_alloc(dpg->dnc * sizeof(UINT));
         if ( stack == NULL ) {
            error = DBMS_MEMORY_ERROR;
         } else {

            push( 0 );      /* Wurzel auf Stack */

            while ( ! stack_empty ) {

               dnode = dna + pop;

               i  = dnode->comp;

               okflag = 1;
               if ( ! (bound.mode[i] & LBUL) ) {
                  cvn = proj( tbl, i, tbl->page_buf + (dia+dnode->ds)->position );
                  cvl = (UCHAR*)proj( tbl, i, bound.low );
                  cmpres = (index_db[i].cmpfunc)( cvl, cvn, index_db[i].length);
                  if ( cmpres > 0 ) {
                     okflag = 0;
                  } else if ( cmpres == 0 ) {
                     if ( ! ((dnode->mode & EQ_SMALLER) && (bound.mode[i]&LBEQ)) ) {
                        okflag = 0;
                     }
                  }
               }
               if ( okflag ) {
                  if ( dnode->mode & SMALLER_DS ) {
                     error = check_and_out( dia + dnode->smaller, usr );
                     if ( error ) break;
                  } else {
                     push( dnode->smaller );
                  }
               }

               okflag = 1;
               if ( ! (bound.mode[i] & UBUL) ) {
                  cvn = proj( tbl, i, tbl->page_buf + (dia+dnode->ds)->position );
                  cvu = (UCHAR*)proj( tbl, i, bound.up );
                  cmpres = (index_db[i].cmpfunc)( cvu, cvn, index_db[i].length);
                  if ( cmpres < 0 ) {
                     okflag = 0;
                  } else if ( cmpres == 0 ) {
                     if ( (dnode->mode & EQ_SMALLER) || ! (bound.mode[i]&UBEQ) ) {
                        okflag = 0;
                     }
                  }
               }
               if ( okflag ) {
                  if ( dnode->mode & GREATER_DS ) {
                     error = check_and_out( dia + dnode->greater, usr );
                     if ( error ) break;
                  } else {
                     push( dnode->greater );
                  }
               }
            }

            memory_free( stack );
         }
      }
   }

   return( error );
}

/* ----------------------------------------------------------------------- */
/* Name       check_and_out                                                */
/* Definition int check_and_out( TABLE *table, DINFO *di );                */
/* Prototyp   Modulintern                                                  */
/* Funktion   Die Funktion prueft, ob der gefundene Datensatz (di) in den   */
/*            Suchgrenzen liegt. Desweiteren wird die Berechtigung des     */
/*            Benutzers zum Lesen dieses Satzes geprueft. Wenn es sich um   */
/*            einen Verweis auf eine Overlayseite handelt, so wird die     */
/*            Suche auf dieser Seite fortgesetzt, ansonsten wird der Satz  */
/*            mit der jeweils gesetzten Ausgabefunktion ausgegeben         */
/* Benutzte Funktionen                                                     */
/*            search_overlay (Modulinterne Funktion)                       */
/*            projection     (Modulinternes Makro)                         */
/*            outfunc        (Modulinterner Zeiger)                        */
/* ----------------------------------------------------------------------- */
int check_and_out( DINFO *di, USER_RIGHTS *usr ) { 

   UINT       i;                   /* zu vergleichende Komponente         */
   UINT       len;                 /* Satzlaenge                           */
   char      *cvl;                 /* Zeiger auf i-te Komp. der Untergr.  */
   char      *cvu;                 /* Zeiger auf i-te Komp. der Obergr.   */
   char      *cvd;                 /* Zeiger auf i-te Komp. des Satzes    */
   SINT       cmpres;              /* Vergleichsergebnis                  */
   char       okflag = 1;          /* Flag fuer Ausgabekontrolle           */
   int        error = DBMS_ALL_OK;

   for ( i=0; i<index_count; i++ ) { /* check all index components */

      len = index_db[i].length;

      cvd = projection( pagebuffer + di->position, i );

      if ( ! (bound.mode[i] & LBUL) ) { /* check lower bound */

         cvl = projection( bound.low, i );

         cmpres = (index_db[i].cmpfunc)( cvd, cvl, index_db[i].length );

         if ( cmpres < 0 || ( cmpres == 0 && ! (bound.mode[i] & LBEQ) ) ) {
            okflag = 0;
         }
      }

      if ( ! (bound.mode[i] & UBUL) ) { /* check upper bound */

         cvu = projection( bound.up, i  );

         cmpres = (index_db[i].cmpfunc)( cvd, cvu, index_db[i].length );

         if ( cmpres > 0 || ( cmpres == 0 && ! (bound.mode[i] & UBEQ) ) ) {
            okflag = 0;
         }
      }
   }

   if ( okflag ) {

      if ( di->type & DIOVERLAY ) {

         if (LogFile!=NULL) fprintf(LogFile,"search overlay\n");

         error = search_overlay( (UCHAR*)pagebuffer + (di->position), di->lognum );

      } else {  /* output of record */

         error = outfunc( pagebuffer + (di->position) );
         if ( ! error ) { counter++; }
      }

   }

   return( error );
}

/* ----------------------------------------------------------------------- */
/* Name       search_overlay                                               */
/* Definition int search_overlay(TABLE *table,UCHAR *i_data,ULONG offset); */
/* Prototyp   Modulintern                                                  */
/* Funktion   Ausgabe der Daten von Overlay-Seiten. i_data enthaelt die    */
/*            Indexkomponente der Daten auf dieser Overlayseite. Aus ihr   */
/*            und den auf der Seite stehenden Restkomponenten werden die   */
/*            auszugebenden Datensaetze zusammengebaut.                    */
/* Benutze Funktionen                                                      */
/*            memory_alloc                                                 */
/*            memory_free                                                  */
/*            memory_copy                                                  */
/*            read_page                                                    */
/* ----------------------------------------------------------------------- */
int search_overlay( UCHAR *i_data, ULONG offset ) {

   char      *datap = pagebuffer;
   OINFO     *oi0, *oi1;
   UCHAR     *overlay1;
   UCHAR     *overlay2;
   OPAGE     *p0, *p1;
   UCHAR     *r_data;
   UCHAR     *source;
   UCHAR     *dest;
   UINT      *pos;
   unsigned   length;
   UINT       n0, n1;
   int        error = DBMS_ALL_OK;

   if (offset == table.active ) {
      printf("fatal error: bad overlay root position (%ld)\n", offset );
   }

   /* -------------------------------------------------------------------- */
   /* INITIALISIERUNG: Es wird ein Seitenpuffer fuer die Overlay-Seiten    */
   /* eingerichtet. Doppelte Seitengroesse, da Overlaybaum max. Hoehe 2 hat. */
   /* Weiterhin wird ein Puffer fuer die auszugebenden Daten bereitgestellt.*/
   /* Dieser hat halbe Seitengroesse, bei variabler Satzlaenge, sonst die   */
   /* Datensatzlaenge.                                                     */
   /* -------------------------------------------------------------------- */

   overlay1 = (UCHAR*)memory_alloc( base->page_size );
   overlay2 = (UCHAR*)memory_alloc( base->page_size );

   if( overlay1 == NULL || overlay2 == NULL ) {
      if( overlay1 != NULL ) memory_free( overlay1 );
      if( overlay2 != NULL ) memory_free( overlay2 );
      error = DBMS_MEMORY_ERROR;
   } else {

      p0  = (OPAGE*)(overlay1 + sizeof(PTYPE));
      p1  = (OPAGE*)(overlay2 + sizeof(PTYPE));
      oi0 = (OINFO*)(overlay1 + sizeof(PTYPE) + sizeof(OPAGE) );
      oi1 = (OINFO*)(overlay2 + sizeof(PTYPE) + sizeof(OPAGE));
      table.page_buf = overlay1;
      table.active   = offset;

      if ( base->variable ) {

         memory_copy( databuffer, i_data, *(UINT*)i_data );
         pos = (UINT*)databuffer;
         length = pos[0];

      } else {

         memory_copy( databuffer, i_data, index_db[index_count].position );
      }

      if (LogFile!=NULL) fprintf(LogFile,"read root page (overlay) %ld", offset );

      error = read_page( &table );

      if ( ! error ) {

         /* -------------------------------------------------------------- */
         /* Durchlaufe Seite 0 des Overlaybaumes. Falls ein DATA_INFO auf  */
         /* eine weitere Seite verweist, so wird zunaechst diese geholt    */
         /* und dann dort mit der Suche fortgefahren. Bei jeder gefundenen */
         /* Restkomponente wird die Leseberechtigung des Benutzers geprueft */
         /* und falls erlaubt der komplette Datensatz zusammengesetzt und  */
         /* ausgegeben.                                                    */
         /* -------------------------------------------------------------- */

         for( n0 = 0; n0 < p0->oic; n0++ ) {

            if ( oi0[n0].type & DIOVERLAY ) {

               table.page_buf = overlay2;
               table.active   = oi0[n0].lognum;

               error = read_page( &table );

               if ( ! error ) {

                  for( n1 = 0; n1 < p1->oic; n1++ ) {

                        if ( (oi1[n1].type & 7) <= (table.user_rights).search ) {

                        r_data = table.page_buf + oi1[n1].position;
                        source = r_data;

                        if( ! base->variable ) {

                           dest = (UCHAR*)databuffer + index_db[index_count].position;
                           memory_copy( dest, source, oi1[n1].size );

                        }else{

                           dest = (UCHAR*)databuffer + pos[index_count+1];
                           memory_copy( dest, source, oi1[n1].size );
                           *(UINT*)(databuffer) = length + oi1[n1].size;
                        }

                        error = outfunc( databuffer );
                        if ( ! error ) { counter++; }
                     }
                  }
               }

               table.page_buf = overlay1;
               table.active = offset;

            } else {

               /* -------------------------------------------------------- */
               /* Baue Datensatz zum Ausgeben zusammen:                    */
               /* Bei fester Laenge muss nur die Restkomponente in out_buf */
               /* kopiert werden, sonst wird die Index- und Restkomponente */
               /* hineinkopiert                                            */
               /* -------------------------------------------------------- */

               if ( (oi0[n0].type & 7) <= (table.user_rights).search ) {

                  r_data = table.page_buf + oi0[n0].position;
                  source = r_data;

                  if( ! base->variable ) {

                     dest = (UCHAR*)databuffer + index_db[index_count].position;
                     memory_copy( dest, source, oi0[n0].size );

                  } else {

                     dest = (UCHAR*)databuffer + pos[index_count+1];
                     memory_copy( dest, source, oi0[n0].size);
                     *(UINT*)databuffer = length + oi0[n0].size;
                  }

                  error = outfunc( databuffer );

                  if ( ! error ) { counter++; }
               }
            }
            if ( error ) break;
         }
      }
      memory_free( overlay1 );
      memory_free( overlay2 );
   }

   table.page_buf = (UCHAR*)pagebuffer;

   return( error );
}

/* ----------------------------------------------------------------------- */
/* Name       push_to_plist                                                */
/* Definition int push_to_plist( PLIST *plist, ULONG offset );             */
/* Prototyp   Modulintern                                                  */
/* Funktion   Tries to put offset to the page offset stack implemented by  */
/*            PLIST. If there are no more free positions in plist, the     */
/*            list is allocated with PLISTCNT more elements.               */
/*            If push_to_plist fails, DBMS_MEMORY_ERROR is returned.       */
/* ----------------------------------------------------------------------- */
int push_to_plist( PLIST *plst, ULONG offset ) {

   ULONG *h;
   int    error = DBMS_ALL_OK;

   if ( plst->act == plst->max ) {

      h = (ULONG*)memory_alloc( (plst->max + PLISTCNT) * sizeof(ULONG));

      if ( h == NULL ) {
         error = DBMS_MEMORY_ERROR;
      } else {
         memory_copy( h, plst->list, (plst->max) * sizeof(ULONG) );
         memory_free( plst->list );
         plst->list = h;
         (plst->max) += PLISTCNT;
         *((plst->list) + plst->act ) = offset;
         (plst->act)++;
      }
   } else {
      *((plst->list) + plst->act ) = offset;
      (plst->act)++;
   }

   return( error );
}

/* ----------------------------------------------------------------------- */
/* Default Ausgabefunktion: keine Ausgabe                                  */
/* ----------------------------------------------------------------------- */
int defaultout( char *d ) { return( DBMS_ALL_OK ); }

/* ----------------------------------------------------------------------- */
/* Name       dbms_select_f                                                */
/* Definition int dbms_select_f( char *defname, int msg )                  */
/* Funktion   Einlesen der Suchdefinitionsdatei, Ausfuellen der             */
/*            SEARCH_PARM Struktur und Aufruf der Suchfunktion             */
/*            dbms_select                                                  */
/* ----------------------------------------------------------------------- */

/* Globale Variablen fuer dbms_select_f */
FILE     *dbmsoutfile = NULL;
unsigned  dbmsoutreclen = 0;
int       dbmsfout( char *data );

int dbms_select_f( char *defname, int msg ) {

   SELECT_PARM  sp;
   TDF          tdf;
   FILE        *deffile, *tdffile;
   char         line[1024];
   int          i = -1, n, error = 0, j;
   unsigned     l, lc = 0, s = 0, e = 0, p;
   char         icheck[INDEXC];
   char        *ilow[INDEXC], *iup[INDEXC];
   unsigned     ilowsize[INDEXC], iupsize[INDEXC];

   /* --- reset boundary definition fields -------------------------------- */
   for ( n=0; n<INDEXC; n++ ) {
      icheck[n] = 0;  
      ilow[n] = NULL; 
      ilowsize[n] = 0;        
      iup[n] = NULL; 
      iupsize[n] = 0;
   }
   sp.low = NULL;
   sp.up = NULL; 
   sp.tblname = NULL;
   dbmsoutfile = NULL; 
   sp.outfunc = NULL;

   /* --- open search defintion file -------------------------------------- */
   deffile = fopen( defname, "r" );      
   if ( deffile == NULL ) {              

      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Fehler: kann %s nicht oeffnen\n", defname ); break;
         default: printf("Error: Unable to open %s\n", defname );
      }

      return( DBMS_SDF_READ_ERROR );
   }

   switch ( msg ) {
      case 0:  break;
      case 1:  printf("Lese %s ein ...\n", defname ); break;
      default: printf("Reading %s ...\n", defname );
   }

   /* --- read first line of search definition file: from <tblname> ------- */
   if ( fgets( line, 1024, deffile ) == NULL ) {

      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Fehler: Zeile 1 - Lesefehler\n"); break;
         default: printf("Error: Line 1 - unable to read\n");
      }
      fclose( deffile );
      return( DBMS_SDF_READ_ERROR );

   }

   l = strlen(line); 
   if ( line[l-1] == '\n' ) line[l-1] = '\0';

   /* --- look if 'from' is the first word of line 1 ---------------------- */
   n = strlen("from ");

   if ( strncmp( line, "from ", n ) != 0 ) {

      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Syntaxfehler: Zeile 1 - 'from' fehlt \n"); break;
         default: printf("Syntax Error: Line 1 - 'from' missing\n");
      }

      fclose( deffile );
      return( DBMS_SDF_SYNTAX_ERROR );

   }

   /* --- skip all blanks behind 'from' to get position of <tblname> ------ */
   while ( line[n] == ' ' ) n++;
   if ( line[n] == '\0' ) {
      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Syntaxfehler: Zeile 1 - Tabellenname fehlt\n"); break;
         default: printf("Syntax Error: Line 1 - tablename missing\n");
      }
      fclose( deffile );
      return( DBMS_SDF_SYNTAX_ERROR );
   }

   /* --- if tablename ends with '.tbl', remove this extension ------------ */
   l = strlen( line+n );
   if ( l > strlen(TBLEXT) ) {
      if ( strcmp( line+n+l-strlen(TBLEXT), TBLEXT ) == 0 ) {
         line[n+l-strlen(TBLEXT)] = '\0';
      }
   }

   /* --- get length of tablename ----------------------------------------- */
   l = strlen(line+n) + 1 + strlen(TDFEXT);
   if ( getenv("DBMS_TDF_PATH") != NULL ) { 
      l += strlen(getenv("DBMS_TDF_PATH")); 
   }

   /* --- allocate string for tablename ----------------------------------- */
   sp.tblname = (char*)memory_alloc(l);
   if ( sp.tblname == NULL ) {
      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Fehler: zuwenig Speicher fuer Tabellenname\n"); break;
         default: printf("Error: no memory for tablename\n");
      }
      fclose( deffile );
      return( DBMS_MEMORY_ERROR );
   }

   /* --- open tdf, first look in current dir, then in DBMS_TDF_PATH ------ */
   makefname( sp.tblname, NULL, line+n, TDFEXT );
   tdffile = fopen(sp.tblname, "rb" );
   if ( tdffile == NULL && getenv("DBMS_TDF_PATH") != NULL ) {
      makefname( sp.tblname, getenv("DBMS_TDF_PATH"), line+n, TDFEXT );
      tdffile = fopen( sp.tblname, "rb" );
   }

   if ( tdffile == NULL ) {

      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Fehler: %s nicht gefunden\n", sp.tblname ); break;
         default: printf("Error: %s not found\n", sp.tblname );
      }

      fclose( deffile );
      memory_free( sp.tblname );
      return( DBMS_TDF_READ_ERROR );

   }

   /* --- read table defintion file --------------------------------------- */
   if ( fread( &tdf, sizeof(TDF), 1, tdffile ) != 1 ) {
      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Fehler: kann %s nicht lesen\n", sp.tblname ); break;
         default: printf("Error: unable to read %s\n", sp.tblname );
      }
      fclose( tdffile );
      fclose( deffile );
      memory_free( sp.tblname );
      return( DBMS_TDF_READ_ERROR );
   }

   /* --- close table definition file and reset tablename in sp ----------- */
   fclose( tdffile );
   strcpy( sp.tblname, line+n );

   /* --- read line 2 of search definition file: to <outfile> ------------- */
   if ( fgets( line, 1024, deffile ) == NULL ) {

      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Fehler: Zeile 2 - Lesefehler\n"); break;
         default: printf("Error: Line 2 - read error\n");
      }
      fclose( deffile );
      memory_free( sp.tblname );
      return( DBMS_SDF_READ_ERROR );

   }

   /* --- get length of line read ----------------------------------------- */
   l = strlen(line); 
   if (line[l-1]=='\n') line[l-1] = '\0';

   /* --- check 'to' in front of 2nd line --------------------------------- */
   n = strlen("to ");
   if ( strncmp( line, "to ", n ) != 0 ) {
      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Syntaxfehler: Zeile 2 - 'to' fehlt\n"); break;
         default: printf("Syntax Error: Line 2 - 'to' missing\n");
      }
      fclose( deffile );
      memory_free( sp.tblname );
      return( DBMS_SDF_SYNTAX_ERROR );
   }

   /* --- skip blanks behind 'to', to get position of <outfile> ----------- */
   while ( line[n] == ' ' ) n++;
   if ( line[n] == '\0' ) {
      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Fehler: Zeile 2 - Ausgabedatei fehlt\n"); break;
         default: printf("Error: Line 2 - output file missing\n");
      }
      fclose( deffile );
      memory_free( sp.tblname );
      return( DBMS_SDF_SYNTAX_ERROR );
   }

   /* --- open outputfile ------------------------------------------------- */
   dbmsoutfile = fopen( line+n, "wb" );
   if ( dbmsoutfile == NULL ) {
      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Fehler: kann Ausgabedatei %s nicht oeffnen\n", line+n ); break;
         default: printf("Error: unable to open output file %s\n", line+n );
      }
      fclose( deffile );
      memory_free( sp.tblname );
      return( DBMS_OUTPUT_ERROR );
   }

   /* --- set record length and output function --------------------------- */
   dbmsoutreclen = tdf.reclen;
   sp.outfunc    = dbmsfout;

   lc = 2;
   while ( fgets(line,1024,deffile) != NULL ) {

      lc++;

      l = strlen(line); 
      if ( line[l-1] == '\n' ) line[l-1] = '\0';

      /* --- if line format is: where <indexname> ------------------------- */
      p = strlen("where ");
      if ( strncmp( line, "where ", p ) == 0 ) {

         /* --- skip blanks and set '\0' behind indexname ----------------- */
         while ( line[p] == ' ' ) p++;
         n = p; while ( line[n]>32 ) n++; line[n] = '\0';

         /* --- search for indexname in table definition ------------------ */
         i = -1;
         for ( j=0; j<tdf.icnt; j++ ) {
            if ( strncmp( line+p, tdf.name[j], 10 ) == 0 ) { 
               i = j; 
               break;
            }
         }

         /* --- if indexname not found ------------------------------------ */
         if ( i == -1 ) {
            switch ( msg ) {
               case 0:  break;
               case 1:  printf("Fehler: Zeile %u - kein Index %s vorhanden\n", lc, line+p ); break;
               default: printf("Error: Line %u - no Index %s\n", lc, line+p );
            }
            error = DBMS_SDF_INDEX_ERROR;
            break;
         }

         /* --- if index was alread set ----------------------------------- */
         if ( icheck[i] > 0 ) {
            switch ( msg ) {
               case 0:  break;
               case 1:  printf("Fehler: Zeile %u - Index %s schon gesetzt\n", lc, line+p ); break;
               default: printf("Error: Line %u - Index %s already set\n", lc, line+p );
            }
            error = DBMS_SDF_INDEX_ERROR;
            break;
         }

      } else {

         /* --- no valid sign found --------------------------------------- */
         if ( line[0] != '>' && line[0] != '<' && line[0] != '=' ) break;

         /* --- if no index was set by where ... -------------------------- */
         if ( i == -1 ) {
            switch ( msg ) {
               case 0:  break;
               case 1:  printf("Fehler: Zeile %u - kein Index gesetzt\n", lc ); break;
               default: printf("Error: Line %u - no Index set\n", lc );
            }
            error = DBMS_SDF_INDEX_ERROR;
            break;
         }

         /* --- read low boundary for current index ----------------------- */
         if ( line[0] == '>' ) {

            if ( icheck[i] & 1 ) {
               switch ( msg ) {
                  case 0:  break;
                  case 1:  printf("Fehler: Zeile %u - Untergrenze fuer %s schon gesetzt\n", lc, tdf.name[i] ); break;
                  default: printf("Error: Line %u - low bound for %s already set\n", lc, tdf.name[i] );
               }
               error = DBMS_SDF_INDEX_ERROR;
               break;
            }

            icheck[i] |= 1;
            n = 2; if ( line[1] == '=' ) { icheck[i] |= 2; n++; }
            s = 0; e = tdf.elem[i];
            ilow[i] = (char*)read_txt( tdf.type[i], line+n, &e, &s );
            ilowsize[i] = s;
            printf("low: %ld\n", *(long*)(ilow[i]) );
            if ( ilow[i] == NULL ) {
               switch ( msg ) {
                  case 0:  break;
                  case 1:  printf("Fehler beim Einlesen der Untergrenze von %s\n", tdf.name[i] ); break;
                  default: printf("Error reading low bound of %s\n", tdf.name[i] );
               }
               error = DBMS_MEMORY_ERROR;
               break;
            }

         /* --- read upper boundary for current index --------------------- */
         } else if ( line[0] == '<' ) {

            if ( icheck[i] & 4 ) {
               switch ( msg ) {
                  case 0:  break;
                  case 1:  printf("Fehler: Zeile %u - Obergrenze fuer %s schon gesetzt\n", lc, tdf.name[i] ); break;
                  default: printf("Error: Line %u - upper bound for %s already set\n", lc, tdf.name[i] );
               }
               error = DBMS_SDF_INDEX_ERROR;
               break;
            }

            icheck[i] |= 4;
            n = 2; if ( line[1] == '=' ) { icheck[i] |= 8; n++; }
            s = 0; e = tdf.elem[i];
            iup[i] = (char*)read_txt( tdf.type[i], line+n, &e, &s );
            iupsize[i] = s;
            if ( iup[i] == NULL ) {
               switch ( msg ) {
                  case 0:  break;
                  case 1:  printf("Fehler: kein Speicher fuer Obergrenze von %s\n", tdf.name[i] ); break;
                  default: printf("Error: no memory for upper bound of %s\n", tdf.name[i] );
               }
               error = DBMS_MEMORY_ERROR;
               break;
            }

         /* --- read equal boundary --------------------------------------- */
         } else if ( line[0] == '=' ) {

            if ( icheck[i] & 1 || icheck[i] & 4 ) {
               switch ( msg ) {
                  case 0:  break;
                  case 1:  printf("Fehler: Zeile %u - Grenzen fuer %s schon gesetzt\n", lc, tdf.name[i] ); break;
                  default: printf("Error: Line %u - bounds for %s already set\n", lc, tdf.name[i] );
               }
               error = DBMS_SDF_INDEX_ERROR;
               break;
            }

            icheck[i] |= 1+2+4+8;
            n = 2; if ( line[1] == '=' ) n++;
            s = 0; e = tdf.elem[i];
            ilow[i] = (char*)read_txt( tdf.type[i], line+n, &e, &s );
            ilowsize[i] = s;
            if ( ilow[i] == NULL ) {
               switch ( msg ) {
                  case 0:  break;
                  case 1:  printf("Fehler: kein Speicher fuer Grenzen von %s\n", tdf.name[i] ); break;
                  default: printf("Error: no memory for bounds of %s\n", tdf.name[i] );
               }
               error = DBMS_MEMORY_ERROR;
               break;
            }
            iup[i] = (char*)memory_alloc(s);
            if ( iup[i] == NULL ) {
               switch ( msg ) {
                  case 0:  break;
                  case 1:  printf("Fehler: kein Speicher fuer Grenzen von %s\n", tdf.name[i] ); break;
                  default: printf("Error: no memory for bounds of %s\n", tdf.name[i] );
               }
               error = DBMS_MEMORY_ERROR;
               break;
            }

            memory_copy( iup[i], ilow[i], s );
            iupsize[i] = s;
         }
      }

   }

   /* --- close search definition file ------------------------------------ */
   fclose( deffile );

   /* --- if error occured reading index boundaries ----------------------- */
   if ( error ) {
      memory_free( sp.tblname );
      for ( i=0; i<INDEXC; i++ ) {
         if ( ilow[i] != NULL ) memory_free( ilow[i] );
         if ( iup[i]  != NULL ) memory_free( iup[i] );
      }
      return( error );
   }

   /* --- write message: reading of search definition file finished ------- */
   switch ( msg ) {
      case 0:  break;
      case 1:  printf("Lesen von %s abgeschlossen\n", defname ); break;
      default: printf("Reading of %s ready\n", defname );
   }

   s = 0; l = 0;

   for ( i=0; i<tdf.icnt; i++ ) {

      sp.mode[i] = 0;

      if ( icheck[i] & 1 ) {
         sp.mode[i] |= LOW_GREATER;
         if ( icheck[i] & 2 ) sp.mode[i] |= LOW_EQUAL;
         l += ilowsize[i];
      } else { sp.mode[i] |= LOW_UNLIMITED; }

      if ( icheck[i] & 4 ) {
         sp.mode[i] |= UP_SMALLER;
         if ( icheck[i] & 8 ) sp.mode[i] |= UP_EQUAL;
         s += iupsize[i];
      } else { sp.mode[i] |= UP_UNLIMITED; }
   }

   if ( tdf.reclen == 0 ) {
      l = l + (tdf.icnt+2)*sizeof(int);
      s = s + (tdf.icnt+2)*sizeof(int);
   } else {
      l = tdf.reclen;
      s = tdf.reclen;
   }

   sp.low = (char*)memory_alloc(l);
   if ( sp.low == NULL ) {
      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Fehler: kein Speicher fuer Untergrenzen\n"); break;
         default: printf("Error: no memory for low boundary\n");
      }
      memory_free( sp.tblname );
      for ( i=0; i<INDEXC; i++ ) {
         if ( ilow[i] != NULL ) memory_free( ilow[i] );
         if ( iup[i]  != NULL ) memory_free( iup[i] );
      }
      return( DBMS_MEMORY_ERROR );
   }

   if ( tdf.reclen == 0 ) {
      *(unsigned*)(sp.low) = l;
      *((unsigned*)(sp.low)+1) = (tdf.icnt+2)*sizeof(int);
   }

   sp.up  = (char*)memory_alloc(s);
   if ( sp.up == NULL ) {
      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Fehler: kein Speicher fuer Obergrenzen\n"); break;
         default: printf("Error: no memory for upper boundary\n");
      }
      memory_free( sp.tblname );
      memory_free( sp.low );
      for ( i=0; i<INDEXC; i++ ) {
         if ( ilow[i] != NULL ) memory_free( ilow[i] );
         if ( iup[i]  != NULL ) memory_free( iup[i] );
      }
      return( DBMS_MEMORY_ERROR );
   }

   if ( tdf.reclen == 0 ) {
      *(unsigned*)(sp.up) = s;
      *((unsigned*)(sp.up)+1) = (tdf.icnt+2)*sizeof(int);
   }

   if ( tdf.reclen == 0 ) {

      for ( i=0; i<tdf.icnt; i++ ) {

         p = *((unsigned*)(sp.low) + i + 1);
         if ( ilowsize[i] > 0 ) memory_copy( sp.low + p, ilow[i], ilowsize[i] );
         *((unsigned*)(sp.low)+i+2) = p + ilowsize[i];

         p = *((unsigned*)(sp.up) + i + 1);
         if ( iupsize[i] > 0 ) memory_copy( sp.up + p, iup[i], iupsize[i] );
         *((unsigned*)(sp.up)+i+2) = p + iupsize[i];
      }

   } else {

      for ( i=0; i<tdf.icnt; i++ ) {
         if ( ilowsize[i] > 0 ) memory_copy( sp.low + tdf.pos[i], ilow[i], ilowsize[i] );
         if ( iupsize[i] > 0 ) memory_copy( sp.up + tdf.pos[i], iup[i], iupsize[i] );
      }
   }

   for ( i=0; i<INDEXC; i++ ) {
      if ( ilow[i] != NULL ) memory_free( ilow[i] );
      if ( iup[i]  != NULL ) memory_free( iup[i] );
   }

   switch ( msg ) {
      case 0:  break;
      case 1:  printf("Starte Suche ...\n"); break;
      default: printf("Start searching ... \n");
   }

   sp.count = 0;
   error = dbms_select( &sp );

   if ( ! error ) {
      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Fertig: %ld Datensaetze gefunden\n", sp.count ); break;
         default: printf("Ready: %ld records found\n", sp.count );
      }
   } else {
      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Fehler: Suche liefert Fehler %d\n", error ); break;
         default: printf("Error: Search failed with error %d\n", error );
      }
   }

   if ( sp.low != NULL ) memory_free( sp.low );
   if ( sp.up  != NULL ) memory_free( sp.up );
   if ( sp.tblname != NULL ) memory_free( sp.tblname );
   if ( dbmsoutfile != NULL ) fclose ( dbmsoutfile );

   return( error );
}

/* ----------------------------------------------------------------------- */
/* Ausgabefunktion fuer Binaerdateien                                        */
/*                                                                         */
/* Schreibt Binaerdate mit den Datensaetzen im DB-Format                     */
/* ----------------------------------------------------------------------- */
int dbmsfout( char *d ) {
   if ( dbmsoutreclen == 0 ) {
      if ( fwrite( d, *(unsigned*)d, 1, dbmsoutfile ) != 1 )
         return( DBMS_OUTPUT_ERROR );
   } else {
      if ( fwrite( d, dbmsoutreclen, 1, dbmsoutfile ) != 1 )
         return( DBMS_OUTPUT_ERROR );
   }
   return( DBMS_ALL_OK );
}

/* ----------------------------------------------------------------------- */
/* Ende von SELECT.C                                                       */
/* ----------------------------------------------------------------------- */
