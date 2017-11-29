/* ----------------------------------------------------------------------- */
/* SPLITIH.C                                                               */
/*                                                                         */
/* Hilfsfunktionen zur Indexseitenspaltung                                 */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>

#include "dbms.h"
#include "dbmsdef.h"
#include "dbmstype.h"
#include "dbmsprot.h"

/* ----------------------------------------------------------------------- */
/* Name       i_new_root                                                   */
/* Definition UINT i_new_root( TABLE *table, ULONG *offset, IPS *ips );    */
/* Prototyp   dbmsprot.h                                                   */
/* Funktion   Erzeugt eine neue Wurzel nach Spaltung einer Indexseite.     */
/* ----------------------------------------------------------------------- */
UINT i_new_root( TABLE *table, ULONG *offset, IPS *ips ) {

   BASE_PAGE *base   = (BASE_PAGE*)(table->base);
   UCHAR     *old_p  = (UCHAR*)(table->page_buf);
   IPAGE     *ip_old = (IPAGE*)(ips->old_page + sizeof(PTYPE));
   INODE     *in_old = (INODE*)(ips->old_page + sizeof(PTYPE) + sizeof(IPAGE));
   BPOS      *bp_old = (BPOS*)( ips->old_page + ip_old->bpos );
   IPAGE     *ip_new;
   INODE     *in_new;
   BPOS      *bp_new;
   ULONG     *ol_new;
   UCHAR     *buffer;
   UINT      *bp_used;
   UINT       b_s = 0;
   UINT       b_c = 0;
   UINT       b,k,n;
   UINT       pos;
   UINT       result = DB_ALL_OK;

   (base->level)++;

   buffer  = (UCHAR*)memory_alloc( base->page_size );
   bp_used = (UINT*)memory_alloc( ip_old->b_c * sizeof(UINT));

   if ( buffer == NULL || bp_used == NULL ) {                                     
      if( buffer  != NULL ) memory_free( buffer );
      if( bp_used != NULL ) memory_free( bp_used );

      result = DB_MEMORY_ERROR;
   } else {
      table->page_buf = buffer;
      newPage( table, IP );
      ip_new = (IPAGE*)(buffer + sizeof(PTYPE));

      /* ----------------------------------------------------------------- */
      /* Ermittle Groesse und Anzahl der benutzten Grenzen                 */
      /* ----------------------------------------------------------------- */

      for ( b = 0; b < ip_old->b_c; b++ ) bp_used[b] = 0;

      for ( k = 0; k < ips->cnt; k++ ) {
         n = ips->ipsn[k].node; b = in_old[n].bnum;

         if( bp_used[ b ] == 0 ){   b_s += bp_old[ b ].size;   b_c++; }

         bp_used[ b ] = 1;   
      }

      /* ----------------------------------------------------------------- */
      /* Speichere die benoetigten Grenzen                                 */
      /* ----------------------------------------------------------------- */

      ip_new->bpos = base->page_size - b_s - b_c * sizeof(BPOS);
      ip_new->b_c  = b_c;

      bp_new = (BPOS*)(buffer + ip_new->bpos);

      pos = base->page_size;
      b_c = 0;

      for ( b = 0; b < ip_old->b_c; b++ ) {
         if ( bp_used[b] > 0 ) {
            pos -= bp_old[b].size;

            bp_used[b] = b_c;

            bp_new[ b_c ].size    = bp_old[b].size;
            bp_new[ b_c ].comp    = bp_old[b].comp;
            bp_new[ b_c ].counter = 0;
            bp_new[ b_c ].pos     = pos;

            b_c++;

            memory_copy( buffer+pos, old_p+bp_old[b].pos, bp_old[b].size );
         }
      }

      /* ----------------------------------------------------------------- */
      /* Speichere Positionen der Nachfolgeseiten                          */
      /* ----------------------------------------------------------------- */

      ip_new->olist = ip_new->bpos - 2 * sizeof(ULONG);
      ip_new->o_c = 2;
      ol_new = (ULONG*)(buffer + ip_new->olist );
      ol_new[0] = ips->ext_page;
      ol_new[1] = ips->int_page;

      /* ----------------------------------------------------------------- */
      /* Speichere neue Knoten                                             */
      /* ----------------------------------------------------------------- */

      in_new = (INODE*)(buffer + sizeof(PTYPE) + sizeof(IPAGE));
      ip_new->node_c = ips->cnt;

      for ( k = 0; k < ips->cnt; k++ ) {
         n = ips->ipsn[k].node; b = in_old[n].bnum;

         in_new[k].cnt = 0;
         in_new[k].flag.equal_smaller = in_old[n].flag.equal_smaller;
         in_new[k].bnum = bp_used[ b ];
         (bp_new[ bp_used[b] ].counter)++;

         if ( k < ips->cnt - 1 ) {
            if ( ips->ipsn[k].flag.next_node_smaller ) {
               in_new[k].flag.smaller_extern = 0;
               in_new[k].smaller = k+1;
               in_new[k].flag.greater_extern = 1;
               in_new[k].flag.low = 0;

               if ( ips->ipsn[k].flag.other_side_extern ) {
                  if ( ips->ipsn[k].flag.other_side_next ) {
                     in_new[k].flag.greater_extern = 0;
                     in_new[k].greater = ips->ipsn[k].other_next;
                  } else { in_new[k].greater = 0;   }
               } else { in_new[k].greater = 1; }   
            } else {
               in_new[k].flag.smaller_extern = 1;
               in_new[k].flag.greater_extern = 0;
               in_new[k].greater = k+1;
               in_new[k].flag.low = 1;

               if ( ips->ipsn[k].flag.other_side_extern ) {
                  if ( ips->ipsn[k].flag.other_side_next ) {
                     in_new[k].flag.smaller_extern = 0;
                     in_new[k].smaller = ips->ipsn[k].other_next;
                  } else { in_new[k].smaller = 0;   }
               } else { in_new[k].smaller = 1; }   
            }
         } else {
            in_new[k].flag.smaller_extern = 1;
            in_new[k].flag.greater_extern = 1;

            if ( ips->last_int_s == 1 ) {
               in_new[k].smaller  = 1;
               in_new[k].greater  = 0;
               in_new[k].flag.low = 0;
            } else {
               in_new[k].smaller  = 0;
               in_new[k].greater  = 1;
               in_new[k].flag.low = 1;
            }
         }
      }

      in_new[0].cnt = ips->cnt;

      ip_new->free  = ip_new->olist - sizeof(PTYPE) - sizeof(IPAGE) - ip_new->node_c * sizeof(INODE);

      table->active = getFreeFilePosition( table );
      *offset = table->active;
      result = write_page( table );

      memory_free( bp_used );
      memory_free( buffer );
   }

   memory_free( ips->ipsn );
   table->page_buf = old_p;

   return( result );
}

/* ----------------------------------------------------------------------- */
/* Name       i_insert_index                                               */
/* Definition UINT i_insert_index( TABLE *table, ULONG *offset,            */
/*                                                IPS *ips, UINT onum );   */
/* Prototyp   dbmsprot.h                                                   */
/* Funktion   Erzeugt eine neue Wurzel nach Spaltung einer Indexseite.     */
/* ----------------------------------------------------------------------- */
UINT i_insert_index( TABLE *table, ULONG *offset, IPS *ips, UINT onum ) {

   IPAGE  *ip     = (IPAGE*)(table->page_buf + sizeof(PTYPE));
   INODE  *in     = (INODE*)(table->page_buf + sizeof(PTYPE) + sizeof(IPAGE));
   BPOS   *bp     = (BPOS*)(table->page_buf + ip->bpos);
   IPAGE  *ip_old = (IPAGE*)(ips->old_page + sizeof(PTYPE));
   INODE  *in_old = (INODE*)(ips->old_page + sizeof(PTYPE) + sizeof(IPAGE));
   BPOS   *bp_old = (BPOS*)( ips->old_page + ip_old->bpos);
   UINT    b_c = 0;
   UINT    b_s = 0;
   UINT   *bp_used;
   UCHAR  *bex;
   UCHAR  *b1, *b2, *source, *dest;
   ULONG  *list;
   UINT    n,newnum,i,pos,k,b;
   UINT    needed_space;
   UINT    result = DB_ALL_OK;

   bp_used = (UINT*)memory_alloc(  ip_old->b_c * sizeof(UINT) );
   bex     = (UCHAR*)memory_alloc( ip_old->b_c );

   if( bp_used == NULL || bex == NULL ) {
      if( bp_used != NULL ) memory_free( bp_used );
      if( bex     != NULL ) memory_free( bex );
      result = DB_MEMORY_ERROR;
   } else {
      for( b = 0; b < ip_old->b_c; b++ ) {
         bp_used[b] = 0;
         bex[b]     = 0;
      }

      for ( k = 0; k < ips->cnt; k++ ) {
         n = ips->ipsn[k].node; b = in_old[n].bnum;

         if ( bp_used[b] == 0 ) {
            b_c++;
            b_s += bp_old[b].size;
         }

         bp_used[b] = 1;
      }

      for ( b = 0; b < ip_old->b_c; b++ ) {
         if ( bp_used[b] > 0 ) {
            for ( k = 0; k < ip->b_c; k++ ) {
               if ( bp[k].comp == bp_old[b].comp ) {
                  b1 = (UCHAR*)(table->page_buf + bp[k].pos);
                  b2 = ips->old_page   + bp_old[b].pos;
                  i  = bp[k].comp;

                  if ( table->index_db[i].cmpfunc( b1,b2, table->index_db[i].length ) == 0 ) {
                     bp_used[b] = k;
                     bex[b] = 1;
                     b_c--;
                     b_s -= bp[k].size;
                     break;
                  }
               }
            }
         }
      }

      needed_space = b_c * sizeof(BPOS) + b_s + sizeof(ULONG) + ips->cnt * sizeof(INODE);

      if( needed_space > ip->free ) {
         result = DB_SPLIT_PAGE;
      } else {
         /* -------------------------------------------------------------- */
         /* Trage neue Nachfolgeseite ein und aendere Offset der           */
         /* gespaltenen Seite.                                             */
         /* -------------------------------------------------------------- */

         source = (UCHAR*)(table->page_buf + ip->olist);
         dest   = source - sizeof(ULONG) - b_s - b_c * sizeof(BPOS);
         memory_copy( dest, source, ip->o_c * sizeof(ULONG) );
         ip->olist = ip->olist - sizeof(ULONG) - b_s - b_c * sizeof(BPOS);
         list = (ULONG*)(table->page_buf + ip->olist);
         list[ onum ]    = ips->ext_page;
         list[ ip->o_c ] = ips->int_page;
         (ip->o_c)++;

         /* -------------------------------------------------------------- */
         /* Trage neue, noch nicht vorhandene Grenzen ein.                 */
         /* -------------------------------------------------------------- */

         pos = ip->bpos + ip->b_c * sizeof(BPOS);
         source = (UCHAR*)(table->page_buf + ip->bpos);
         dest   = source - b_c * sizeof(BPOS) - b_s;
         memory_copy( dest, source, ip->b_c * sizeof(BPOS) );
         ip->bpos = ip->bpos - b_c * sizeof(BPOS) - b_s;
         bp = (BPOS*)(table->page_buf + ip->bpos);

         for ( b = 0; b < ip_old->b_c; b++ ) {
            if ( bex[b] == 0 && bp_used[b] > 0 ) {
               pos -= bp_old[ b ].size;

               bp[ ip->b_c ].comp    = bp_old[ b ].comp;
               bp[ ip->b_c ].size    = bp_old[ b ].size;
               bp[ ip->b_c ].counter = 0;
               bp[ ip->b_c ].pos     = pos;

               memory_copy( table->page_buf+pos, ips->old_page+bp_old[b].pos, bp_old[b].size );

               bp_used[b] = ip->b_c;
               (ip->b_c)++;
            }
         }

         /* -------------------------------------------------------------- */
         /* Aendere Nachfolgerzeiger aller Knoten, die vorher auf die       */
         /* gespaltene Seite zeigten.                                      */
         /* -------------------------------------------------------------- */

         for ( n = 0; n < ip->node_c; n++ ) {
            if ( in[n].flag.smaller_extern && in[n].smaller == onum ) {
               in[n].flag.smaller_extern = 0;
               in[n].smaller = ip->node_c;
            }

            if ( in[n].flag.greater_extern && in[n].greater == onum ) {
               in[n].flag.greater_extern = 0;
               in[n].greater = ip->node_c;
            }
         }

         /* -------------------------------------------------------------- */
         /* Trage die neuen Knoten ein.                                    */
         /* -------------------------------------------------------------- */

         for ( k = 0; k < ips->cnt; k++ ) {
            newnum = ip->node_c + k; n = ips->ipsn[k].node; b = in_old[n].bnum;

            in[ newnum ].bnum = bp_used[b];
            (bp[ bp_used[b] ].counter)++;
            in[ newnum ].cnt = 0;
            in[ newnum ].flag.equal_smaller = in_old[ n ].flag.equal_smaller;

            if ( k < ips->cnt - 1 ) {
               if ( ips->ipsn[k].flag.next_node_smaller ) {
                  in[ newnum ].flag.smaller_extern = 0;
                  in[ newnum ].smaller = newnum + 1;
                  in[ newnum ].flag.greater_extern = 1;
                  in[ newnum ].flag.low = 0;

                  if ( ips->ipsn[k].flag.other_side_extern ) {
                     if ( ips->ipsn[k].flag.other_side_next ) {
                        in[newnum].flag.greater_extern = 0;
                        in[newnum].greater = ip->node_c + ips->ipsn[k].other_next;
                     } else { in[newnum].greater = onum; }
                  } else { in[newnum].greater = ip->o_c-1; }
               } else {
                  in[ newnum ].flag.smaller_extern = 1;
                  in[ newnum ].flag.greater_extern = 0;
                  in[ newnum ].greater = newnum + 1;
                  in[ newnum ].flag.low = 1;

                  if( ips->ipsn[k].flag.other_side_extern ) {
                     if( ips->ipsn[k].flag.other_side_next ) {
                        in[newnum].flag.smaller_extern = 0;
                        in[newnum].smaller = ip->node_c + ips->ipsn[k].other_next;
                     } else { in[newnum].smaller = onum;   }
                  } else   { in[newnum].smaller = ip->o_c-1; }
               }
            } else {
               in[ newnum ].flag.smaller_extern = 1;
               in[ newnum ].flag.greater_extern = 1;

               if( ips->last_int_s == 1 )
               {
                  in[ newnum ].smaller = ip->o_c - 1;
                  in[ newnum ].greater = onum;
                  in[ newnum ].flag.low = 0;
               }
               else
               {
                  in[ newnum ].smaller = onum;
                  in[ newnum ].greater = ip->o_c - 1;
                  in[ newnum ].flag.low = 1;
               }
            }
         }       

         in[ ip->node_c ].cnt = ips->cnt;
         ip->node_c += ips->cnt;

         ip->free = ip->olist - sizeof(PTYPE) - sizeof(IPAGE) - ip->node_c * sizeof(INODE);

         *offset = getFreeFilePosition( table );
         table->active = *offset;
         result = write_page( table );

         memory_free( ips->ipsn );
      }

      memory_free( bp_used );
      memory_free( bex );
   }

   return( result );
}

/* ----------------------------------------------------------------------- */
/* Ende von SPLITIH.C                                                      */
/* ----------------------------------------------------------------------- */
