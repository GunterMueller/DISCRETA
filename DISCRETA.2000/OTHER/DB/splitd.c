/* ----------------------------------------------------------------------- */
/* SPLITD.C                                                                */
/*                                                                         */
/* Funktion split_data_page und Hilfsfunktion get_tree_to_extrakt zum      */
/* Spalten von ueberlaufenden Datenseiten. Die Seiten werden so zerlegt,    */
/* dass die die beiden neu entstehenden Seiten zwischen 1/3 und 2/3 der     */
/* Datengroesse der urspruenglichen Seite enthalten. Der neu einzufuegende     */
/* Datensatz wird ebenfalls gespeichert. Die neuen Seiten werden           */
/* geschrieben. Alle wichtigen Informationen zum Spalten werden in der     */
/* Struktur dps zurueckgegeben.                                             */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>

#include "dbmsdef.h"
#include "dbmstype.h"
#include "dbmsprot.h"

/* ----------------------------------------------------------------------- */
/* Makros fuer das Stackhandling beim Suchen des zu entfernenden Teilbaums  */
/* ----------------------------------------------------------------------- */
#define push( x )   { stack[spos] = (x); spos++; }
#define pop( x )    { if( spos > 0 ){ spos--; (x) = stack[spos]; } }
#define stack_empty ( ( spos > 0 )?(0):(1))


/* ----------------------------------------------------------------------- */
/* Modulintern Funktionen                                                  */
/* ----------------------------------------------------------------------- */
UINT get_tree_to_extract( DNODE *dn, UINT *lowlist, UINT *uplist );
void get_bs_bc( TABLE *table, DPS *dps, UINT *bs, UINT *bc );
void put_bounds( TABLE *table, UINT bc, UINT bs, DPS *dps, UINT onum );
void put_nodes( TABLE *table, DPS *dps, UINT in_page, UINT ext_page );


/* ----------------------------------------------------------------------- */
/* Name                                                                    */
/* Definition                                                              */
/* Prototyp                                                                */
/* Funktion                                                                */
/* ----------------------------------------------------------------------- */
int split_data_page( TABLE *table, DPS *dps, UINT in ) {

   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   UCHAR     *pagebuffer = table->page_buf;
   ULONG      active = table->active;
   DPAGE     *p    = (DPAGE*)(pagebuffer + sizeof(PTYPE));
   DNODE     *dn   = (DNODE*)(pagebuffer + sizeof(PTYPE) + sizeof(DPAGE));
   DINFO     *di   = (DINFO*)(pagebuffer + p->data);
   UCHAR     *buffer;
   UCHAR     *d;
   UINT      *stack;
   UINT       spos = 0, x = 0, i,n,s;
   UINT       index_c = table->index_c;
   UINT       uplist[INDEXC];
   UINT       lowlist[INDEXC];
   DINFO      di_new;
   int        error = DBMS_ALL_OK;

	Log("split_data_page(");
   di_new.position = 0;              /* Init DINFO for new */
   di_new.lognum   = base->lognum;   /* record to insert   */
   di_new.type     = table->ditype;

   for ( i=0; i<index_c; i++ ) {      /* Init boundary position */
      uplist[i]  = 0xFFFF;            /* list with 0xFFFF       */
      lowlist[i] = 0xFFFF;
   }

   dps->up  = NULL; 
   dps->low = NULL; 
   dps->upos[ index_c ] = 0;
   dps->lpos[ index_c ] = 0;
   dps->cnt = 0;

   buffer = (UCHAR*)memory_alloc( base->page_size );   /* Alloc. Speicher fuer Buffer und Stack */
   stack  = (UINT*)memory_alloc( p->dnc * sizeof(UINT) );

   if ( buffer == NULL || stack == NULL ) {
      if ( buffer != NULL ) memory_free( buffer );
      if ( stack  != NULL ) memory_free( stack );
      error = DBMS_MEMORY_ERROR;
   } else {

      /* ----------------------------------------------------------------- */
      /* Suche die Wurzel des zu extrahierenden Teilbaumes                 */
      /* ----------------------------------------------------------------- */
      x = get_tree_to_extract( dn, lowlist, uplist );

      push( x );  /* Wurzel des zu extrahierenden Teilbaumes auf Stack */

      for ( i=0; i<index_c; i++ ) { /* ermittle Groesse der Pufferspeicher low und up */

         if ( lowlist[i] != 0xFFFF ) {

            n = ( dn + lowlist[i] )->ds;
            d = pagebuffer + (di+n)->position;
            *( (dps->lpos) + index_c ) += compSize( table, i, d );
         }

         if ( uplist[i] != 0xFFFF ) {

            n = ( dn + uplist[i] )->ds;
            d = pagebuffer + (di+n)->position;
            *( (dps->upos) + index_c ) += compSize( table, i, d );
         }
      }

      /* ----------------------------------------------------------------- */
      /* Hole Speicher fuer die Puffer low und up, wie oben berechnet.     */
      /* low oder up koennen auch leer sein, da nicht immer eine Ober-     */
      /* oder Untergrenze vorhanden sein muessen                            */
      /* ----------------------------------------------------------------- */

      if ( *( (dps->upos) + index_c ) > 0 ) {
         dps->up  = (UCHAR*)memory_alloc( *( (dps->upos) + index_c ) );
         if ( dps->up == NULL ) {
            error = DBMS_MEMORY_ERROR;
         }
      }

      if ( *( (dps->lpos) + index_c ) > 0 ) {
         dps->low = (UCHAR*)memory_alloc( *((dps->lpos) + index_c ) );
         if ( dps->low == NULL ) {
            error = DBMS_MEMORY_ERROR;
         }
      }

      if ( error == DBMS_MEMORY_ERROR ) {
         if ( dps->up  != NULL ) memory_free( dps->up );
         if ( dps->low != NULL ) memory_free( dps->low );
      } else {

         /* -------------------------------------------------------------- */
         /* Trage nun die Grenzen in low und up ein. Dazu werden jeweils   */
         /* die lpos bzw. upos Werte auf die Anfangspositionen der Daten   */
         /* gesetzt und dann werden in Daten in die Puffer an diese        */
         /* Positionen kopiert. 0xFFFF ist eine Markierung, dass bzgl.     */
         /* dieser Indexkomp. keine Begrenzung nach oben bzw. unten        */
         /* existiert. Wenn dies der Fall ist, setze no_?? Flags. Wenn die */
         /* Grenzen eingetragen sind, wird noch ein Flag gesetzt, wo die   */
         /* letzte eingetragene Grenze ist. (bei welcher Komponente)       */
         /* -------------------------------------------------------------- */

         *( dps->lpos ) = 0; 
         *( dps->upos ) = 0;

         for ( i=0; i<index_c; i++ ) {

            dps->flag[i].last = 0;

            if ( lowlist[i] != 0xFFFF ) {     /* trage untere Grenze ein */

               dps->flag[i].no_low = 0;
               d = pagebuffer + ( di + (dn + lowlist[i])->ds )->position;
               s = compSize( table, i, d );
               dps->flag[i].l_eq_s = ((dn+lowlist[i])->mode & EQ_SMALLER) ? 1 : 0;
               memory_copy( dps->low + dps->lpos[i], proj( table, i, d ), s );
               dps->lpos[i+1] = dps->lpos[i] + s;
               (dps->cnt)++;

            } else {

               dps->flag[i].no_low = 1;
               dps->lpos[i+1] = dps->lpos[i];
            }

            if ( uplist[i] != 0xFFFF ) {      /* trage obere Grenze ein */

               dps->flag[i].no_up = 0;
               d = pagebuffer + ( di + ( dn + uplist[i] )->ds )->position;
               s = compSize( table, i, d );
               dps->flag[i].u_eq_s = ((dn+uplist[i])->mode & EQ_SMALLER) ? 1 : 0;
               memory_copy( dps->up + dps->upos[i], proj( table, i, d ), s );
               dps->upos[i+1] = dps->upos[i] + s;
               (dps->cnt)++;

            } else {
               dps->flag[i].no_up = 1;
               dps->upos[i+1] = dps->upos[i];
            }

         }

         i = index_c;        /* Ermittle Komponente mit letztem Eintrag */

         while ( i > 0 )   {
            i--;
            if ( ! ( dps->flag[i].no_low && dps->flag[i].no_up ) ) {
               dps->flag[i].last = 1; break;
            }
         }

         /* -------------------------------------------------------------- */
         /* Initialisiere den buffer als leere Datenseite und kopiere die  */
         /* Datensaetze des zu extrahierenden Teilbaumes auf diese Seite.   */
         /* Falls der neu einzufuegende Datensatz hierbei auftaucht,        */
         /* kopiere auch diesen.                                           */
         /* -------------------------------------------------------------- */

         table->page_buf = buffer;
         newPage( table, DP );

         while ( ! stack_empty )   {

            pop( x );

            if ( ( dn + x )->mode & SMALLER_DS ) {

               d = pagebuffer + ( di + ( dn + x )->smaller )->position;

               copyData( table, d, di + ( dn+x )->smaller );

               if( ( dn + x )->smaller == in ) {

                  d = table->data_buf;
                  copyData( table, d, &di_new );
               }

               ( di + (dn+x)->smaller )->position = 0;

            } else {
               push( dn[x].smaller );
            }

            if ( dn[x].mode & GREATER_DS ) {

               d = pagebuffer + ( di + ( dn + x )->greater )->position;

               copyData( table, d, di + (dn+x)->greater );

               if( dn[x].greater == in ) {

                  d = table->data_buf;
                  copyData( table, d, &di_new );
               }
               ( di + (dn + x)->greater )->position = 0;
            } else {
               push( dn[x].greater );
            }
         }

         /* ----------------------------------------------------------- */
         /* hole freie Seitenposition fuer die int. Seite in dps und     */
         /* schreibe diese Seite                                        */
         /* ----------------------------------------------------------- */
         table->active = getFreeFilePosition( table );
         dps->in_page  = table->active;
         error         = write_page( table );

         if ( ! error ) {

            newPage( table, DP );

            for ( x = 0; x < p->dic; x++ ) {

               if ( (di+x)->position > 0 )   {

                  d = (UCHAR*)pagebuffer + (di+x)->position;

                  copyData( table, d, di + x );

                  if ( x == in ) {
                     copyData( table, table->data_buf, &di_new );
                  }
               }
            }

            /* ----------------------------------------------------------- */
            /* hole freie Seitenposition fuer die ext. Seite in dps und     */
            /* schreibe diese Seite                                        */
            /* ----------------------------------------------------------- */
            table->active = getFreeFilePosition( table );
            dps->ext_page = table->active;
            error         = write_page( table );
         }

         if ( error ) {              /* if function failed, free */
            memory_free( dps->low ); /* memory of dps low and up */
            memory_free( dps->up );
         }
      }
      memory_free( stack );       /* free memory of the stack   */
      memory_free( buffer );      /* and additional page buffer */
   }

   table->page_buf = (UCHAR*)pagebuffer;  /* set page buffer and active */
   table->active   = active;      /* to previous values         */

	Log(")split_data_page\n");
   return( error );
}

/* ----------------------------------------------------------------------- */
/* Name       get_tree_to_extract                                          */
/* Definition UINT get_tree_to_extract(DNODE *dn,UINT *lowlist,            */
/*                                                          UINT *uplist)  */
/* Prototyp   Modulintern                                                  */
/* Funktion   Ermittle den Teilbaum, der auf die neue Seite kommt.         */
/* ----------------------------------------------------------------------- */
UINT get_tree_to_extract( DNODE *dn, UINT *lowlist, UINT *uplist ) {

   float  ratio;
   DNODE *node = dn;
   UINT   x = 0;
   UINT   i;

	Log("get_tree_to_extract(");
   ratio = 1.0;

   while ( ratio > 0.6666 ) {

      i = node->comp;

      if ( node->mode & SMALLER_DS ) { /* suche Teilbaum zum Extrahieren */

         lowlist[i] = x;
         x = node->greater;

      } else {

         if ( node->mode & GREATER_DS )   {

            uplist[i] = x;
            x = node->smaller;

         } else {

            if ( (dn + (node->smaller))->total > (dn + (node->greater))->total ) {
               uplist[i] = x;
               x = node->smaller;
            } else {
               lowlist[i] = x;
               x = node->greater;
            }
         }
      }

      node = dn + x;
      ratio = (float)(node->total) / (float)(dn->total);
   }

	Log(")get_tree_to_extract\n");
   return( x );
}

/* ----------------------------------------------------------------------- */
/* Name       d_new_root                                                   */
/* Definition int d_new_root( TABLE *table, ULONG *offset, DPS *dps )      */
/* Prototyp   dbmsprot.h                                                   */
/* Funktion   Anlegen einer neuen Wurzelseite nach Datenseitensplit.       */
/* ----------------------------------------------------------------------- */
int d_new_root( TABLE *table, ULONG *offset, DPS *dps ) {

   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   IPAGE     *ip = (IPAGE*)(table->page_buf + sizeof(PTYPE));
   INODE     *in = (INODE*)(table->page_buf + sizeof(PTYPE) + sizeof(IPAGE));
   BPOS      *bpos;
   ULONG     *ol;
   UINT       i,s, n, pos = base->page_size;
   UINT       R = table->index_c;
   int        error = DB_ALL_OK;

	Log("d_new_root(");
   (base->level)++;      /* Baumhoehe wird um 1 heraufgesetzt */
   newPage( table, IP );  /* init. Seitenpuffer als leere Indexseite */

   ip->bpos = base->page_size - dps->lpos[ R ] - dps->upos[ R ] - dps->cnt * sizeof(BPOS);
   bpos = (BPOS*)(table->page_buf + ip->bpos);

   ip->olist = ip->bpos - 2 * sizeof(ULONG);    /* trage Seitenoffsets ein */
   ol = (ULONG*)(table->page_buf + ip->olist);

   ol[0] = dps->in_page;
   ol[1] = dps->ext_page;
   ip->o_c    = 2;
   ip->b_c    = 0;
   ip->node_c = 0;

   in[ip->node_c].cnt = dps->cnt;

   for ( n = 1; n < dps->cnt; n++ ) in[ ip->node_c + n ].cnt = 0;

   for ( i = 0; i < table->index_c; i++ ) {
      if ( ! dps->flag[i].no_low ) {
         s = dps->lpos[i+1] - dps->lpos[i];
         pos -= s;
         memory_copy( table->page_buf + pos, dps->low + dps->lpos[i], s );

         bpos[ ip->b_c ].counter = 1;
         bpos[ ip->b_c ].pos     = pos;
         bpos[ ip->b_c ].comp    = i;
         bpos[ ip->b_c ].size    = s;

         in[ ip->node_c ].bnum    = ip->b_c;
         in[ ip->node_c ].smaller = 1;
         in[ ip->node_c ].greater = ip->node_c + 1;
         in[ ip->node_c ].flag.smaller_extern = 1;
         in[ ip->node_c ].flag.greater_extern = 0;
         in[ ip->node_c ].flag.equal_smaller  = dps->flag[i].l_eq_s;
         in[ ip->node_c ].flag.low = 1;

         if( dps->flag[i].last && dps->flag[i].no_up ) {
            in[ ip->node_c ].greater = 0;
            in[ ip->node_c ].flag.greater_extern = 1;
         }

         (ip->b_c)++;
         (ip->node_c)++;
      }

      if( ! dps->flag[i].no_up ) {
         s = dps->upos[i+1] - dps->upos[i];
         pos -= s;
         memory_copy( table->page_buf + pos, dps->up + dps->upos[i], s );

         bpos[ ip->b_c ].counter = 1;
         bpos[ ip->b_c ].pos     = pos;
         bpos[ ip->b_c ].comp    = i;
         bpos[ ip->b_c ].size    = s;

         in[ ip->node_c ].bnum    = ip->b_c;
         in[ ip->node_c ].greater = 1;
         in[ ip->node_c ].smaller = ip->node_c + 1;
         in[ ip->node_c ].flag.smaller_extern = 0;
         in[ ip->node_c ].flag.greater_extern = 1;
         in[ ip->node_c ].flag.equal_smaller  = dps->flag[i].u_eq_s;
         in[ ip->node_c ].flag.low = 0;

         if( dps->flag[i].last ) {
            in[ ip->node_c ].smaller = 0;
            in[ ip->node_c ].flag.smaller_extern = 1;
         }

         (ip->b_c)++;
         (ip->node_c)++;
      }
   }

   ip->free  = ip->olist - sizeof(PTYPE) - sizeof(IPAGE) - ip->node_c * sizeof(INODE);

   table->active = getFreeFilePosition( table );
   *offset = table->active;
   error = write_page( table );

	Log(")d_new_root\n");
   return( error );
}


/* ----------------------------------------------------------------------- */
/* Name       d_insert_index                                               */
/* Definition int d_insert_index( TABLE *table, ULONG *offset, DPS *dps,   */
/*                                                           UINT onum )   */
/* Prototyp   dbmsprot.h                                                   */
/* Funktion   Einfuegen des Spaltungsindex aus dps auf der Indexseite.      */
/* ----------------------------------------------------------------------- */
int d_insert_index( TABLE *table, ULONG *offset, DPS *dps, UINT onum ) {

   IPAGE *ip = (IPAGE*)(table->page_buf + sizeof(PTYPE));
   UINT   bs = 0;
   UINT   bc = 0;
   UINT   needed_space = 0;
   int    error = DB_ALL_OK;

	Log("d_insert_index(");
   get_bs_bc( table, dps, &bs, &bc );

   needed_space = dps->cnt * sizeof(INODE) + sizeof(ULONG) + bc * sizeof(BPOS) + bs;

   if( needed_space > ip->free ) {
      error = DB_SPLIT_PAGE;
   } else {          /* fuege die neuen Knoten ein */
      put_bounds( table, bc, bs, dps, onum );
      put_nodes( table, dps, ip->o_c - 1, onum );

      ip->free = ip->olist - sizeof(PTYPE) - sizeof(IPAGE) - ip->node_c * sizeof(INODE);

      table->active = getFreeFilePosition( table );
      *offset = table->active;
      error = write_page( table );
   }

	Log(")d_insert_index\n");
   return( error );
}

/* ----------------------------------------------------------------------- */
/* Name       get_bs_bc                                                    */
/* Definition void get_bs_bc( TABLE *table, DPS *dps, UINT *bs, UINT *bc );*/
/* Funktion   Die Funktion ermittelt, ob eine von den neu einzufuegenden    */
/*            Grenzwerten schon welche auf der Seite vorhanden sind. Falls */
/*            ja, wird fuer die jeweilige Grenze das _exists Flag gesetzt   */
/*            und in bpl bzw. bpu die neue Nummer der Grenze vermerkt.     */
/*            Somit kann die Anzahl bc der tatsaechlich neu einzufuegenden   */
/*            Grenzen und deren Groesse bs ermittelt werden.                 */
/* ----------------------------------------------------------------------- */
void get_bs_bc( TABLE *table, DPS *dps, UINT *bs, UINT *bc ) {     

   IPAGE *ip = (IPAGE*)(table->page_buf + sizeof(PTYPE));
   BPOS  *bp = (BPOS*)(table->page_buf + ip->bpos);
   UCHAR *v1, *v2;
   UINT   i,n;
   UINT   R = table->index_c;
	int no_low, no_up, low_exists, up_exists;

	Log("get_bs_bc(");
   *bs = dps->lpos[ R ] + dps->upos[ R ];
   *bc = dps->cnt;

   for ( i = 0; i < R; i++ ) { /* belege flags mit 0 vor */
      dps->flag[i].low_exists = 0; dps->flag[i].up_exists  = 0;
   }
	{
	char s[256];
	
	sprintf(s, "R = %ld ip->b_c = %ld\n", 
		(long) R, (long) ip->b_c );
	Log(s);
	}

   for ( n = 0; n < ip->b_c; n++ ) {      /* durchlaufe vorhandene Grenzen */

      i = bp[n].comp;

	{
	char s[256];
	
	sprintf(s, "n = %ld i = %ld\n", 
		(long) n, (long) i );
	Log(s);
	}
	Log("vor low\n");
	no_low = dps->flag[i].no_low;
	low_exists = dps->flag[i].low_exists;
	{
	char s[256];
	
	sprintf(s, "no_low = %ld low_exists = %ld\n", 
		(long) no_low, (long) low_exists );
	Log(s);
	}
	
#if 0
      if ( ! dps->flag[i].no_low && ! dps->flag[i].low_exists ) { /* wenn Untergrenze und noch nicht ex. */
#endif
      if ( ! no_low && ! low_exists ) { /* wenn Untergrenze und noch nicht ex. */

	{
	char s[256];
	
	sprintf(s, "bp[n].pos = %ld\n", 
		(long) bp[n].pos );
	Log(s);
	}
	{
	char s[256];
	
	sprintf(s, "dps->lpos[i] = %ld\n", 
		(long) dps->lpos[i] );
	Log(s);
	}
         v1 = table->page_buf + bp[n].pos;  /* vergleiche die beiden Grenzen */
	if (dps->low == NULL)
		Log("dps->low == NULL");
         v2 = dps->low + dps->lpos[i];         

         if( table->index_db[i].cmpfunc( v1, v2, table->index_db[i].length ) == 0 ) {
            dps->flag[i].low_exists = 1; 
            dps->bpl[i] = n;
            *bs -= dps->lpos[i+1] - dps->lpos[i]; 
            (*bc)--;
         }
      }

	Log("vor up\n");
      if( ! dps->flag[i].no_up && dps->flag[i].up_exists ) { /* ebenso mit der Obergrenze */

         v1 = table->page_buf + bp[n].pos;
         v2 = dps->up + dps->upos[i];         

         if( table->index_db[i].cmpfunc( v1, v2, table->index_db[i].length ) == 0 ) {
            dps->flag[i].up_exists = 1; 
            dps->bpu[i] = n;
            *bs -= dps->upos[i+1] - dps->upos[i]; 
            (*bc)--;
         }
      }
	Log("nach up\n");
   }

	Log(")get_bs_bc\n");
   return;
}


/* ----------------------------------------------------------------------- */
/* Name       put_bounds                                                   */
/* Definition void put_bounds( TABLE *table, UINT bc, UINT bs,             */
/*                                                DPS *dps, UINT onum );   */
/* Prototyp   dbmsprot.h                                                   */
/* Funktion                                                                */
/* ----------------------------------------------------------------------- */
void put_bounds( TABLE *table, UINT bc, UINT bs, DPS *dps, UINT onum ) {

   IPAGE *ip = (IPAGE*)(table->page_buf + sizeof(PTYPE));
   BPOS  *bp;
   UCHAR *dest, *source;
   UINT   bstart = ip->bpos + ip->b_c * sizeof(BPOS);
   ULONG *olist;
   UINT   size;
   UINT   i;

	Log("put_bounds(");
   /*---< trage geaenderten und neuen Seitenoffset ein >-----*/

   source = table->page_buf + ip->olist;
   dest   = source - sizeof(ULONG) - bs - bc * sizeof(BPOS);
   size   = ip->o_c * sizeof(ULONG);
   memory_copy( dest, source, size );

   ip->olist = ip->olist - sizeof(ULONG) - bs - bc * sizeof(BPOS);
   olist = (ULONG*)(table->page_buf + ip->olist );
   olist[ onum ]    = dps->ext_page;
   olist[ ip->o_c ] = dps->in_page;
   (ip->o_c)++;

   /*---< trage neue Grenzen ein >-----------------------------*/

   source = table->page_buf + ip->bpos;
   dest   = source - bs - bc * sizeof(BPOS);
   size   = ip->b_c * sizeof(BPOS);
   memory_copy( dest, source, size );

   ip->bpos = ip->bpos - bs - bc * sizeof(BPOS);
   bp = (BPOS*)(table->page_buf + ip->bpos);

   for ( i = 0; i < table->index_c; i++ ) {

      if ( ! dps->flag[i].no_low && ! dps->flag[i].low_exists ) {

         size    = dps->lpos[i+1] - dps->lpos[i];
         bstart -= size;
         source  = dps->low + dps->lpos[i];
         dest    = table->page_buf + bstart;

         memory_copy( dest, source, size );

         (bp + ip->b_c)->pos     = bstart;
         (bp + ip->b_c)->comp    = i;
         (bp + ip->b_c)->size    = size;
         (bp + ip->b_c)->counter = 0;

         dps->bpl[i] = ip->b_c;
         (ip->b_c)++;
      }

      if ( ! dps->flag[i].no_up && ! dps->flag[i].up_exists ) {

         size   = dps->upos[i+1] - dps->upos[i];
         bstart -= size;
         source = dps->up + dps->upos[i];
         dest   = table->page_buf + bstart;

         memory_copy( dest, source, size );

         (bp + ip->b_c)->pos     = bstart;
         (bp + ip->b_c)->comp    = i;
         (bp + ip->b_c)->size    = size;
         (bp + ip->b_c)->counter = 0;

         dps->bpu[i] = ip->b_c;
         (ip->b_c)++;
      }
   }

	Log(")put_bounds\n");
   return;
}

/* ----------------------------------------------------------------------- */
/* Name       put_nodes                                                    */
/* Definition void put_nodes( TABLE *table, DPS *dps,                      */
/*                                UINT in_page, UINT ext_page );           */
/* Prototyp   Modulintern                                                  */
/* Funktion   Schreibe die Knoten aus dps von Seite in_page                */
/* ----------------------------------------------------------------------- */
void put_nodes( TABLE *table, DPS *dps, UINT in_page, UINT ext_page ) {

   IPAGE *ip = (IPAGE*)(table->page_buf + sizeof(PTYPE));
   INODE *in = (INODE*)(table->page_buf + sizeof(PTYPE) + sizeof(IPAGE));
   BPOS  *bp = (BPOS *)(table->page_buf + ip->bpos);
   UINT   n, i;

	Log("put_nodes(");
   for ( n=0; n<ip->node_c; n++ ) {     /* durchlaufe active-Feld */

      if ( (in+n)->flag.smaller_extern && (in+n)->smaller == ext_page ) {
         (in+n)->smaller             = ip->node_c;
         (in+n)->flag.smaller_extern = 0;
      }

      if ( (in+n)->flag.greater_extern && (in+n)->greater == ext_page ) {
         (in+n)->greater             = ip->node_c;
         (in+n)->flag.greater_extern = 0;
      }
   }

   (in+(ip->node_c))->cnt = dps->cnt;
   for ( n=1; n<dps->cnt; n++ ) (in+(ip->node_c)+n)->cnt = 0;

   for ( i=0; i<table->index_c; i++ ) {

      if ( ! dps->flag[i].no_low ) {

         ((bp + dps->bpl[i])->counter)++;

         (in + ip->node_c)->bnum    = dps->bpl[i];
         (in + ip->node_c)->greater = ip->node_c + 1;
         (in + ip->node_c)->smaller = ext_page;

         (in + ip->node_c)->flag.smaller_extern = 1;
         (in + ip->node_c)->flag.greater_extern = 0;
         (in + ip->node_c)->flag.equal_smaller  = dps->flag[i].l_eq_s;
         (in + ip->node_c)->flag.low            = 1;

         if ( dps->flag[i].last && dps->flag[i].no_up ) {
            (in + ip->node_c)->greater = in_page;
            (in + ip->node_c)->flag.greater_extern = 1;
         }

         (ip->node_c)++;
      }

      if ( ! dps->flag[i].no_up ) {

         ((bp + dps->bpu[i])->counter)++;

         (in + ip->node_c)->bnum    = dps->bpu[i];
         (in + ip->node_c)->greater = ext_page;
         (in + ip->node_c)->smaller = ip->node_c + 1;

         (in + ip->node_c)->flag.smaller_extern = 0;
         (in + ip->node_c)->flag.greater_extern = 1;
         (in + ip->node_c)->flag.equal_smaller  = dps->flag[i].u_eq_s;
         (in + ip->node_c)->flag.low            = 0;

         if ( dps->flag[i].last ) {
            (in + ip->node_c)->smaller = in_page;
            (in + ip->node_c)->flag.smaller_extern = 1;
         }
         (ip->node_c)++;
      }
   }

	Log(")put_nodes\n");
   return;
}

/* ----------------------------------------------------------------------- */
/* Ende von SPLITD.C                                                       */
/* ----------------------------------------------------------------------- */
