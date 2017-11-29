/* ----------------------------------------------------------------------- */
/* INSERTH.C                                                               */
/*                                                                         */
/* Dieses Modul enthaelt alle Hilfsfunktionen die zum Einfuegen von Daten-   */
/* saetzen benoetigt werden.                                                 */
/*                                                                         */
/* Funktionen                                                              */
/*                                                                         */
/*                                                                         */
/* Notiz;   Der umstaendliche Zugriff auf table->index_db[ ... ] bzw.      */
/*          list[ ... ] kommt aufgrund von Sideeffects auf page_c in       */
/*          insert,c die durch diese Zugriffsweise vermieden werden        */
/*          => auf keinen Fall aendern.                                    */
/* ----------------------------------------------------------------------- */

/* #include <iostream.h> */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dbms.h"          /* Benutzeroperationen/strukturen               */
#include "dbmsdef.h"       /* Definitionen                                 */
#include "dbmstype.h"      /* Interne Verwaltungsstrukturen                */
#include "dbmsprot.h"      /* Prototypen der internen Funktionen           */

/* ----------------------------------------------------------------------- */
/* Prototypen fuer dateiinterne Funktionen                                  */
/* ----------------------------------------------------------------------- */
int  ins_data( TABLE *table, UINT *in, ULONG *pages_to_del, UINT *del_c );
UINT split_index_page( TABLE *table, DPS *dps, UINT onum, IPS *ips, UINT mode );
UINT i_new_root( TABLE *table, ULONG *offset, IPS *ips );
UINT i_insert_index( TABLE *table, ULONG *offset, IPS *ips, UINT onum );

/* ----------------------------------------------------------------------- */
/* Name       get_next_page                                                */
/* Definition ULONG get_next_page( TABLE *table, UINT *next );             */
/* Prototyp   dbmsprot.h                                                   */
/* Funktion   Diese Funktion wird in insert und delete Operationen be-     */
/*            noetigt, um die Nachfolgeseite einer Indexseite beim Suchen   */
/*            des Blattes zum Einfuegen zu finden. Auf der Seite befindet   */
/*            sich ein INODE-Baum. Dieser wird durchlaufen, bis ein EXIT   */
/*            Knoten gefunden wird. Der Offset der Nachfolgeseite ist dann */
/*            in der Liste olist zu finden.                                */
/* ----------------------------------------------------------------------- */
ULONG get_next_page( TABLE *table, UINT *next ) {

   IPAGE *page  = (IPAGE*)(table->page_buf + sizeof(PTYPE) );
   INODE *inpos = (INODE*)(table->page_buf + sizeof(PTYPE) + sizeof(IPAGE));
   ULONG *list  = (ULONG*)(table->page_buf + page->olist );
   BPOS  *bpos  = (BPOS *)(table->page_buf + page->bpos );
   INODE *node;         /* Aktiver Knoten                                  */
   INDEX *index_db;        /* Index des Vergleichswertes des Knotens          */
   UCHAR *vdata;        /* Vergleichskomponente des einzufuegenden Satzes   */
   UCHAR *vnode;        /* Vergleichsweret beim Knoten                     */
   SINT   cmpres;       /* Ergebnis des Vergleichs                         */

   if( page->node_c == 0 ) {

      *next = 0;

   } else {

      node = inpos;

      while( 1 ) {

         vdata = proj( table, ((bpos + node->bnum)->comp), table->data_buf );
         vnode = table->page_buf + (bpos + node->bnum)->pos;
         index_db = (table->index_db) + (bpos + node->bnum)->comp;

         cmpres = (index_db->cmpfunc)( vdata, vnode, index_db->length );

         if ( cmpres < 0 || ( cmpres == 0 && node->flag.equal_smaller ) ) {

            if( node->flag.smaller_extern ) {
               *next  = node->smaller;            break;
            } else {
               node = inpos + node->smaller;
            }

         } else {

            if( node->flag.greater_extern ) {
               *next  = node->greater;            break;
            } else {
               node = inpos + node->greater;
            }
         }
      }
   }

   return( *(list + *next) );
}

/* ------------------------------------------------------------------------- */
/* Name       insert_data                                                    */
/* Definition int insert_data( TABLE *table, PAGEBUF *pagebuf, UINT *last,   */
/*                       ULONG *pages_to_del, UINT *del_c, ULONG *offset );  */
/* Prototyp   dbmsprot.h                                                     */
/* Funktion   Hauptfunktion zum Einfuegen der Daten. Fuegt auf Datenseite ein  */
/*            und uebernimmt, falls noetig auch die Steuerung der Seiten-      */
/*            spaltung.                                                      */
/* ------------------------------------------------------------------------- */
int insert_data( TABLE *table, PAGEBUF *pagebuf, UINT *last,
                           ULONG *pages_to_del, UINT *del_c, ULONG *offset ) {

   UINT       in_ds = 0;
   UINT       next;
   DPS        dps;
   IPS        ips;
   int        error = DBMS_ALL_OK;

   ips.ipsn = NULL;
   ips.cnt  = 0;

	Log("insert_data(");
	memory_zero(&dps, (int) sizeof(DPS));
	memory_zero(&ips, (int) sizeof(IPS));
   error = ins_data( table, &in_ds, pages_to_del, del_c );

   if ( ! error ) {

      table->active = getFreeFilePosition( table );
      *offset = table->active;
      pagebuf[*last].offset = *offset;

      if ( table->operation != FIO ) {
         error = write_page( table );
         pagebuf[*last].offset = 0;
      }

   } else {

      if ( error == DB_SPLIT_PAGE ) {

         error = split_data_page( table, &dps, in_ds );

         if ( ! error ) {

            if ( *last == 0 ) {

               error = d_new_root( table, offset, &dps );

            } else {

               (*last)--;
               table->page_buf = (UCHAR*)(pagebuf[*last].page);
               next = pagebuf[*last].next;

               error = d_insert_index( table, offset, &dps, next );

               if ( error == DB_SPLIT_PAGE ) {

                  error = split_index_page( table, &dps, next, &ips, 0 );

                  if ( error == DB_SPLIT_CRASH ) {

                     printf("split index page crashed\n");
                  }

                  if ( ! error ) {

                     if ( *last == 0 ) {

                        error = i_new_root( table, offset, &ips );

                     } else {

                        (*last)--;
                        next = pagebuf[*last].next;
                        table->page_buf = (UCHAR*)(pagebuf[*last].page);
                        error = i_insert_index( table, offset, &ips, next );
                     }
                  }
               }

               while ( error == DB_SPLIT_PAGE ) {

                  error = split_index_page( table, NULL, next, &ips, 1 );

                  if ( error == DB_SPLIT_CRASH ) {

                     printf("split index page crashed\n");
                  }

                  if ( ! error ) {

                     if ( *last == 0 ) {
                        error = i_new_root( table, offset, &ips );
                     } else {
                        (*last)--;
                        next = pagebuf[*last].next;
                        table->page_buf = (UCHAR*)(pagebuf[*last].page);
                        error = i_insert_index( table, offset, &ips, next );
                     }
                  }
               }
            }
         }
      }
   }

	Log(")insert_data\n");
   return( error );
}

/* ------------------------------------------------------------------------- */
/* Name       ins_data                                                       */
/* Definition int ins_data( TABLE *table, UINT *in, ULONG *pages_to_del,     */
/*                                   UINT *del_c );                          */
/* Prototyp   dbmsprot.h                                                     */
/* Funktion   Einfuegen eines Datensatzes auf einer Datenseite bzw. falls     */
/*            noetig auf einer Overlayseite.                                  */
/* ------------------------------------------------------------------------- */
int ins_data( TABLE *table, UINT *in, ULONG *pages_to_del, UINT *del_c ) {

   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   DPAGE     *p  = (DPAGE*)(table->page_buf + sizeof(PTYPE));
   DNODE     *dn = (DNODE*)(table->page_buf + sizeof(PTYPE) + sizeof(DPAGE));
   DINFO     *di = (DINFO*)(table->page_buf + p->data);
   UINT       length;
   UINT       I[INDEXC];
   UINT      *list;
   UINT       lc = 0;
   UINT       n,i;
   UCHAR     *v1;
   UCHAR     *v2;
   UINT       index_id = 1;
   UCHAR     *source;
   UCHAR     *dest;
   SINT       res;
   UINT       x,z;
   UINT       x_next;
   int        error = DBMS_ALL_OK;

	Log("ins_data(");
   if ( base->variable ) {
      length = *(UINT*)(table->data_buf);
   } else {
      length = table->data_length;
   }

   if ( p->dic == 0 ) {
	Log("ins_data: p->dic == 0\n");

      /* ----------------------------------------------------------------- */
      /* Es befindet sich kein Datensatz auf der Seite. Positioniere die   */
      /* Seitenparameter free und data auf die Position nach DPAGE bzw.    */
      /* auf die Position Seitenende - Datensatzlaenge - Groesse eines DINFO */
      /* Initialisiere dann DINFO fuer den neuen Datensatz und kopiere ihn */
      /* nach dem DINFO auf die Seite. Schreibe die Seite.                 */
      /* ----------------------------------------------------------------- */

      p->dic  = 1;
      p->dnc  = 0;
      p->free = sizeof(PTYPE) + sizeof(DPAGE);
      p->data = base->page_size - length - sizeof(DINFO);
	{
	char str[256];

	sprintf(str, "sizeof(PTYPE) = %ld\n", (long) sizeof(PTYPE)); Log(str);
	sprintf(str, "sizeof(DPAGE) = %ld\n", (long) sizeof(DPAGE)); Log(str);
	sprintf(str, "sizeof(DNODE) = %ld\n", (long) sizeof(DNODE)); Log(str);
	sprintf(str, "sizeof(DINFO) = %ld\n", (long) sizeof(DINFO)); Log(str);
	sprintf(str, "base->page_size = %ld\n", (long) base->page_size); Log(str);
	sprintf(str, "length = %ld\n", (long) length); Log(str);
	sprintf(str, "p->data = %ld\n", (long) p->data); Log(str);
	}

      di = (DINFO*)(table->page_buf + p->data);       /* Init. DINFO */
 	Log("ins_data: di[0].position\n");
     di[0].position = p->data + sizeof(DINFO);
 	Log("ins_data: di[0].type\n");
      di[0].type     = table->ditype;
/*    di[0].lognum   = base->lognum; */

      memory_copy( table->page_buf + di[0].position, table->data_buf, length ); /* Kopiere */

   } else {

      if( p->dnc == 0 ) {
	Log("ins_data: p->dnc == 0\n");

         /* -------------------------------------------------------------- */
         /* Es befindet sich nur ein Datensatz auf der Seite.              */
         /* -------------------------------------------------------------- */

         getIndexOrder( table, I );

         for( n = 0; n < table->index_c; n++ ) {

            i = I[n];
            v1 = proj( table, i, table->data_buf);
            v2 = proj( table, i, table->page_buf + di[0].position );

            res = (table->index_db[i].cmpfunc)( v1, v2, table->index_db[i].length );

            if ( res != 0 ) {
               index_id = 0; break;
            }
         }

         if ( index_id ) {

            if ( di[0].type & DIOVERLAY ) {
               if (LogFile!=NULL) fprintf(LogFile,"insert overlay\n");
               error = insert_overlay( table, 0, pages_to_del, del_c );
            } else {
	       if (LogFile!=NULL) fprintf(LogFile,"create new overlay\n");
               error = create_new_overlay( table, 0 );
            }

         } else {

            dn[0].ds    = 0;
            dn[0].comp  = i;
            dn[0].total = base->page_size - p->data - sizeof(DINFO) + length;
            dn[0].mode  = SMALLER_DS | GREATER_DS;

            if ( res > 0 ) {
               dn[0].smaller = 0; dn[0].greater = 1;
               dn[0].mode |= EQ_SMALLER;
            } else {
               dn[0].smaller = 1; dn[0].greater = 0;
            }

            source = table->page_buf + p->data;
            dest   = source - sizeof(DINFO) - length;
            memory_copy( dest, source, sizeof(DINFO) );

            p->data -= sizeof(DINFO) + length;
            p->free += sizeof(DNODE);

            di = (DINFO*)(table->page_buf + p->data);
            di[1].position = p->data + (p->dic+1)*sizeof(DINFO);
            di[1].type     = table->ditype;
/*          di[1].lognum   = base->lognum; */

            memory_copy( table->page_buf + di[1].position, table->data_buf, length );

            p->dic = 2; p->dnc = 1;
         }

      } else {
	Log("ins_data: mehrere Datensaetze auf der Seite\n");

         /* -------------------------------------------------------------- */
         /* Es befinden sich mehrere Datensaetze auf der Seite.            */
         /* -------------------------------------------------------------- */

         list=(UINT*)memory_alloc( p->dic * sizeof(UINT));

         if ( list == NULL ) {
            error = DBMS_MEMORY_ERROR;
         } else {

            x = 0;

            while ( 1 ) {   /* suche Datensatz bei dem eingefuegt wuerde */

               list[lc] = x; lc++; i = dn[x].comp;

               v1 = proj( table, i, table->data_buf);
               v2 = proj( table, i, table->page_buf + di[dn[x].ds].position );
               res = table->index_db[i].cmpfunc( v1, v2, table->index_db[i].length );

               switch ( res ) {

                  case  1: x_next = _GREATER; break;

                  case  0:  if ( dn[x].mode & EQ_SMALLER ) x_next = _SMALLER;
                            else x_next = _GREATER;
                            break;

                  case -1:  x_next = _SMALLER; break;
               }

               if ( x_next == _SMALLER ) {
                  if( dn[x].mode & SMALLER_DS ) { z = dn[x].smaller; *in = z; break; }
                  else   x = dn[x].smaller;
               } else {
                  if( dn[x].mode & GREATER_DS ) { z = dn[x].greater; *in = z; break; }
                  else x = dn[x].greater;
               }
            }

            getIndexOrder( table, I );           /* suche trennende Indexkomponente */

            for ( n = 0; n < table->index_c; n++ ) {
               i = I[n];
               v1 = proj( table, i, table->data_buf);
               v2 = proj( table, i, table->page_buf + di[z].position );

               if( (res=(table->index_db[i].cmpfunc)( v1, v2, table->index_db[i].length )) != 0 )
               {
                  index_id = 0; break;
               }
            }

            if ( index_id ) {
	Log("ins_data: if (index_id) \n");

               /* -------------------------------------------------------- */
               /* Der neue Datensatz stimmt mit dem gefundenen auf dem     */
               /* Index ueberein, ist der gefundene Datensatz ein TID, so  */
               /* fuege auf der zugehoerigen Overlay-Seite ein, ansonsten  */
               /* konstruiere neue Overlay-Seite                           */
               /* -------------------------------------------------------- */

               if ( di[z].type & DIOVERLAY ) {  /* versuche auf Overlayseite einzufuegen */

                  error = insert_overlay( table, z, pages_to_del, del_c );

               } else {      /* Konstruiere neue Overlayseite */

                  length = compSize( table, table->index_c, table->page_buf + di[z].position);
                  error = create_new_overlay( table, z );

                  if ( ! error ) {          /* Korregiere Teilbaumgroessen */
                     for( n = 0; n < lc; n++ ) dn[ list[n] ].total -= length;
                  }
               }

            } else {
	Log("ins_data: if (index_id) else \n");

               if ( p->data - p->free < length + sizeof(DINFO) + sizeof(DNODE) ) {

                  error = DB_SPLIT_PAGE;   /* zuwenig Platz => Seite wird gespalten */

               } else {

                  if ( x_next == _SMALLER ) {

                     dn[x].mode &= NOT_SMALLER_DS;
                     dn[x].smaller = p->dnc;
                  } else {
                     dn[x].mode &= NOT_GREATER_DS;
                     dn[x].greater = p->dnc;
                  }

                  if( ! base->variable ) dn[p->dnc].total = table->data_length;
                  else   dn[p->dnc].total = *(UINT*)(table->page_buf + di[z].position);
                  dn[p->dnc].total += length;
                  dn[p->dnc].ds    = z;
                  dn[p->dnc].comp  = i;   /* Init. DNODE */
                  dn[p->dnc].mode  = SMALLER_DS | GREATER_DS;

                  if ( res == 1 ) {
                     dn[p->dnc].smaller = z; dn[p->dnc].greater = p->dic;
                     dn[p->dnc].mode |= EQ_SMALLER;
                  } else {
                     dn[p->dnc].smaller = p->dic; dn[p->dnc].greater = z;
                  }

                  source = table->page_buf + p->data;       /* Verschiebe Datenblock */    
                  dest   = source - sizeof(DINFO) - length;
                  memory_copy( dest, source, p->dic*sizeof(DINFO) );

                  p->data -= sizeof(DINFO) + length;      /* Setze Seitenparameter */
                  p->free += sizeof(DNODE);

                  di = (DINFO*)(table->page_buf + p->data);
                  di[p->dic].position = p->data + (p->dic+1) * sizeof(DINFO);
                  di[p->dic].type     = table->ditype;
/*                di[p->dic].lognum   = base->lognum; */

                  memory_copy(table->page_buf+di[p->dic].position,table->data_buf,length);

                  (p->dic)++;
                  (p->dnc)++;  /* Seitenparamter setzen und Teilbaumgroessen aendern */

                  for ( n = 0; n < lc; n++ ) {
                     i = list[n];
                     dn[ i ].total += length;
                  }
               }
            }
            memory_free( list );
         }
      }
   }

	Log(")ins_data\n");
   return( error );
}

/* ----------------------------------------------------------------------- */
/* Ende von INSERTH.C                                                      */
/* ----------------------------------------------------------------------- */
