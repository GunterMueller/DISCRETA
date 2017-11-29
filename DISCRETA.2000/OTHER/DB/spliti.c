/* ----------------------------------------------------------------------- */
/* SPLITI.C                                                                */
/*                                                                         */
/* Diese Datei und die Datei SPLITIH.C enthalten saemtliche Funktionen      */
/* zur Handhabung der Spaltung von Indexseiten.                            */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
/* #include <conio.h> */

#include "dbms.h"
#include "dbmsdef.h"
#include "dbmstype.h"
#include "dbmsprot.h"


/* ----------------------------------------------------------------------- */
/* Makros fuer die Stackverwaltung                                          */
/* ----------------------------------------------------------------------- */
#define push( x )   {stack[spos]=(x);spos++;}
#define pop( x )    {if(spos>0){spos--;(x)=stack[spos];}}
#define stack_empty ((spos>0)?(0):(1))

/* ----------------------------------------------------------------------- */
/* Modulinterne Funktionen                                                 */
/* ----------------------------------------------------------------------- */
UINT init_page_1( TABLE *table, UINT root, UCHAR *buffer );
UINT init_page_2( TABLE *table, UINT root, UCHAR *buffer );
UINT i_insert_index( TABLE *table, ULONG *offset, IPS *ips, UINT onum );

extern UINT _stop_;

/* ----------------------------------------------------------------------- */
/* Name      split_index_page                                              */
/* Defintion UINT split_index_page( TABLE *table, DPS *dps, UINT onum,     */
/*                                                  IPS *ips, UINT mode )  */
/* Funktion  Spalten von Indexseiten                                       */
/* ----------------------------------------------------------------------- */
UINT split_index_page( TABLE *table, DPS *dps, UINT onum, IPS *ips, UINT mode ) {

   BASE_PAGE *base  = (BASE_PAGE*)(table->base);
   UCHAR     *old_p = table->page_buf;
   IPAGE     *ip    = (IPAGE*)(old_p + sizeof(PTYPE));
   INODE     *in    = (INODE*)(old_p + sizeof(PTYPE) + sizeof(IPAGE));
   IPAGE     *ip_new;
   ULONG     *list;
   UCHAR     *buffer;
   UINT      *stack;
   UINT       spos = 0;
   UINT       x,a, root, n, node;
   UINT       smaller, greater;
   UINT       n_cnt_s, n_cnt_g;
   UINT       splitted_p_in = 0;
   UINT       newnum = 0;
   ULONG      splitted_p = *((ULONG*)(old_p + ip->olist)+onum );
   IPSN      *ipsn;
   UINT       ipsc = 0;
   UINT       last_int_s;
   ULONG      ext_page = 0, int_page = 0;
   UINT       result = DB_ALL_OK;

   stack  = (UINT*)memory_alloc( ip->node_c * sizeof(UINT));
   ipsn   = (IPSN*)memory_alloc( ip->node_c * sizeof(IPSN));
   buffer = (UCHAR*)memory_alloc( base->page_size );

   if ( stack == NULL || buffer == NULL || ipsn == NULL ) {
      if( stack  != NULL ) memory_free( stack );
      if( buffer != NULL ) memory_free( buffer );
      if( ipsn   != NULL ) memory_free( ipsn );
      return( DB_MEMORY_ERROR );
   }

   root = 0;
   ipsc = 0;

   while ( 1 ) {

      for ( n = 0; n < in[root].cnt; n++ ) {

         ipsn[ ipsc + n ].node = root + n;
         ipsn[ ipsc + n ].flag.next_node_smaller = 1;
         ipsn[ ipsc + n ].flag.other_side_extern = 1;
         ipsn[ ipsc + n ].flag.other_side_next   = 0;

         if ( in[root+n].flag.smaller_extern || in[root+n].smaller != root+n+1) {
            ipsn[ ipsc + n ].flag.next_node_smaller = 0;
         }
      }

      spos = 0; n_cnt_s = 0; n_cnt_g = 0;

      a = root + in[root].cnt - 1;

      smaller = in[ a ].smaller; greater = in[ a ].greater;

      if ( ! in[ a ].flag.smaller_extern ) {

         push( smaller );

         while ( ! stack_empty ) {

            pop( x );   x = x + in[x].cnt - 1;
            if( in[ x ].flag.greater_extern ) n_cnt_s++;
            else                              push( in[ x ].greater );

            if( in[ x ].flag.smaller_extern ) n_cnt_s++;
            else                                push( in[ x ].smaller );
         }
      }

      if ( ! in[ a ].flag.greater_extern ) {

         push( greater );

         while ( ! stack_empty ) {

            pop( x ); x = x + in[x].cnt - 1;

            if( in[ x ].flag.greater_extern ) n_cnt_g++;
            else                      push( in[ x ].greater );

            if( in[ x ].flag.smaller_extern ) n_cnt_g++;
            else                              push( in[ x ].smaller );
         }
      }

      if ( n_cnt_s > n_cnt_g ) {
         last_int_s = 1;

         if ( (double)n_cnt_s <= 2.0 * (double)(ip->o_c) / 3.0 ) {
            ipsc += in[root].cnt;
            ipsn[ ipsc-1 ].flag.next_node_smaller = 1;
            root = smaller;
            break;
         } else {

            a = ipsc + in[root].cnt - 1;
            ipsn[ a ].flag.next_node_smaller = 1;

            for ( n = 0; n < in[root].cnt; n++ ) {

               if ( ipsn[ ipsc + n ].flag.next_node_smaller ) {

                  if ( (! in[root+n].flag.greater_extern)
                                    && in[root+n].greater == smaller ) {

                     ipsn[ ipsc + n ].flag.other_side_next = 1;
                     ipsn[ ipsc + n ].other_next = ipsc + in[root].cnt;
                  }

               } else {

                  if ( (! in[root+n].flag.smaller_extern)
                                   && in[root+n].smaller == smaller ) {

                     ipsn[ ipsc + n ].flag.other_side_next = 1;
                     ipsn[ ipsc + n ].other_next = ipsc + in[root].cnt;
                  }
               }
            }
            ipsc += in[root].cnt;
            root = smaller;
         }

      } else {

         last_int_s = 0;

         a = ipsc + in[root].cnt - 1;
         ipsn[ a ].flag.next_node_smaller = 0;

         if ( (double)n_cnt_g <= 2.0 * (double)(ip->o_c) / 3.0 ) {
            ipsc += in[root].cnt;
            ipsn[ ipsc-1 ].flag.next_node_smaller = 0;
            root = greater;
            break;
         } else {

            for ( n = 0; n < in[root].cnt; n++ ) {

               if ( ipsn[ ipsc + n ].flag.next_node_smaller ) {

                  if ( (! in[root+n].flag.greater_extern)
                                   && in[root+n].greater == greater ) {

                     ipsn[ ipsc + n ].flag.other_side_next = 1;
                     ipsn[ ipsc + n ].other_next = ipsc + in[root].cnt;
                  }

               } else {

                  if ( (! in[root+n].flag.smaller_extern)
                                   && in[root+n].smaller == greater ) {

                     ipsn[ ipsc + n ].flag.other_side_next = 1;
                     ipsn[ ipsc + n ].other_next = ipsc + in[root].cnt;
                  }
               }
            }

            ipsc += in[root].cnt;
            root = greater;
         }
      }

   } /* end while */

   for ( n = 0; n < ipsc; n++ ) {

      node = ipsn[n].node;

      if ( ipsn[ n ].flag.next_node_smaller ) {
         if ( ! in[node].flag.greater_extern && in[ node ].greater == root )
            ipsn[n].flag.other_side_extern = 0;
      } else {
         if ( ! in[node].flag.smaller_extern && in[ node ].smaller == root )
            ipsn[n].flag.other_side_extern = 0;
      }
   }

   /* -------------------------------------------------------------------- */
   /* Hier ist die Ermittlung der zu uebertragenden Knoten abgeschlossen   */
   /* -------------------------------------------------------------------- */

   if ( (result = init_page_1( table, root, buffer ) ) == DB_ALL_OK ) {

      table->page_buf = buffer;

      ip_new = (IPAGE*)(buffer + sizeof(PTYPE));
      list = (ULONG*)(buffer + ip_new->olist);

      for ( n = 0; n < ip_new->o_c; n++ ) {
         if( list[n] == splitted_p ) {
            splitted_p_in = 1; newnum = n; break;
         }
      }

      if( splitted_p_in ) {

         if( mode == 0 ) {
            result = d_insert_index( table, &int_page, dps, newnum );
         } else {
            result = i_insert_index( table, &int_page, ips, newnum );
            if( result == DB_SPLIT_PAGE ) result = DB_SPLIT_CRASH;
         }
      } else {
         table->active = getFreeFilePosition( table );
         int_page      = table->active;
         result        = write_page( table );
      }

      table->page_buf = old_p;
   }

   if ( result == DB_ALL_OK ) {

      if ( (result = init_page_2( table, root, buffer ) ) == DB_ALL_OK ) {

         table->page_buf = buffer;

         if ( splitted_p_in ) {
            splitted_p_in = 0;
         } else {

            ip_new = (IPAGE*)(buffer + sizeof(PTYPE));
            list = (ULONG*)(buffer + ip_new->olist);

            for ( n = 0; n < ip_new->o_c; n++ ) {

               if ( list[n] == splitted_p ) {
                  splitted_p_in = 1; newnum = n; break;
               }
            }
         }

         if ( splitted_p_in ) {

            if( mode == 0 ) {
               result = d_insert_index( table, &ext_page, dps, newnum );
            } else {
               result = i_insert_index( table, &ext_page, ips, newnum );
               if( result == DB_SPLIT_PAGE ) result = DB_SPLIT_CRASH;
            }
         } else {
            table->active = getFreeFilePosition( table );
            ext_page      = table->active;
            result        = write_page( table );
         }
         table->page_buf = old_p;
      }
   }

   ips->old_page   = old_p;
   ips->ipsn       = ipsn;
   ips->cnt        = ipsc;
   ips->ext_page   = ext_page;
   ips->int_page   = int_page;
   ips->last_int_s = last_int_s;

   memory_free( buffer );
   memory_free( stack );

   return( result );
}

/* ----------------------------------------------------------------------- */
/* Name       init_page_1                                                  */
/* Definition UINT init_page_1( TABLE *table, UINT root, UCHAR *buffer )   */
/* Funktion   Erzeuge 1. Seite des Spaltergebnisses.                       */
/* ----------------------------------------------------------------------- */
UINT init_page_1( TABLE *table, UINT root, UCHAR *buffer ) {

   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   UCHAR     *old_p = table->page_buf;
   IPAGE     *ip_old = (IPAGE*)(old_p + sizeof(PTYPE));
   INODE     *in_old = (INODE*)(old_p + sizeof(PTYPE) + sizeof(IPAGE));
   BPOS      *bp_old = (BPOS*)(old_p + ip_old->bpos);
   ULONG     *list_old = (ULONG*)(old_p + ip_old->olist);
   IPAGE     *ip_1  = (IPAGE*)(buffer + sizeof(PTYPE));
   INODE     *in_1  = (INODE*)(buffer + sizeof(PTYPE) + sizeof(IPAGE));
   BPOS      *bp_1;
   ULONG     *list;
   UINT       b_c = 0, o_c = 0, b_s = 0;
   UINT      *nodes_used, *o_used, *bp_used;
   UINT       pos;
   UINT       n, b, smaller, greater, node;
   UINT      *stack;
   UINT       spos = 0;
   UINT       result = DB_ALL_OK;

   bp_used    = (UINT*)memory_alloc( ip_old->b_c * sizeof(UINT) );
   o_used     = (UINT*)memory_alloc( ip_old->o_c * sizeof(UINT) );
   stack      = (UINT*)memory_alloc( ip_old->node_c * sizeof(UINT));
   nodes_used = (UINT*)memory_alloc( ip_old->node_c * sizeof(UINT));

   if( bp_used == NULL || o_used == NULL || stack == NULL || nodes_used == NULL ) {
      if( bp_used    != NULL ) memory_free( bp_used );
      if( o_used     != NULL ) memory_free( o_used );
      if( stack      != NULL ) memory_free( stack );
      if( nodes_used != NULL ) memory_free( nodes_used );

      result = DB_MEMORY_ERROR;

   } else {

      table->page_buf = buffer; newPage( table, IP );

      for( b = 0; b < ip_old->b_c; b++ )    bp_used[b] = 0;
      for( n = 0; n < ip_old->o_c; n++ )    o_used[n]  = 0;
      for( n = 0; n < ip_old->node_c; n++ ) nodes_used[n] = 0;

      push( root );

      while ( ! stack_empty ) {

         pop( node );

         nodes_used[node] = 1;

         b = in_old[node].bnum;
         if ( bp_used[b] == 0 ) { b_c++; b_s += bp_old[b].size; }
         bp_used[ b ] = 1;

         smaller = in_old[node].smaller;
         greater = in_old[node].greater;

         if( in_old[node].flag.smaller_extern ) o_used[ smaller ] = 1;
         else                                   push( smaller );

         if( in_old[node].flag.greater_extern ) o_used[ greater ] = 1;
         else                                   push( greater );
      }

      /* ----------------------------------------------------------------- */
      /* Trage einzufuegende Grenzen ein                                   */
      /* ----------------------------------------------------------------- */

      pos        = base->page_size;
      ip_1->bpos = base->page_size - b_s - b_c * sizeof(BPOS);
      ip_1->b_c  = b_c;
      bp_1       = (BPOS*)(buffer + ip_1->bpos);

      b_c = 0;

      for ( b = 0; b < ip_old->b_c; b++ ) {

         if ( bp_used[b] > 0 ) {
            pos -= bp_old[b].size;

            bp_1[b_c].pos     = pos;
            bp_1[b_c].counter = 0;
            bp_1[b_c].size    = bp_old[b].size;
            bp_1[b_c].comp    = bp_old[b].comp;

            bp_used[b] = b_c;

            memory_copy( buffer+pos, old_p+bp_old[b].pos, bp_old[b].size );

            b_c++;
         }
      }

      /* ----------------------------------------------------------------- */
      /* Trage benoetigte Nachfolgeseitenpositionen ein.                   */
      /* ----------------------------------------------------------------- */

      o_c = 0;

      for ( n = 0; n < ip_old->o_c; n++ ) { if( o_used[n] > 0 ) o_c++; }

      ip_1->olist = ip_1->bpos - o_c * sizeof(ULONG);
      ip_1->o_c = o_c;

      list = (ULONG*)(buffer + ip_1->olist);

      o_c = 0;

      for ( n = 0; n < ip_old->o_c; n++ ) {
         if ( o_used[n] > 0 ) {
            list[o_c] = list_old[n];
            o_used[n] = o_c;
            o_c++;
         }
      }

      /* ----------------------------------------------------------------- */
      /* Trage Knoten des extrahierten Teilgraphen ein                     */
      /* ----------------------------------------------------------------- */

      node = 0;

      for ( n = 0; n < ip_old->node_c; n++ ) {

         if ( nodes_used[n] == 1 ) {
            nodes_used[n] = node;

            b = in_old[n].bnum; b = bp_used[b];
            bp_1[b].counter += 1;

            in_1[node].flag    = in_old[n].flag;
            in_1[node].bnum    = b;
            in_1[node].cnt     = in_old[n].cnt;
            in_1[node].smaller = in_old[n].smaller;
            in_1[node].greater = in_old[n].greater;

            node++;
         }
      }

      ip_1->node_c = node;

      /* ----------------------------------------------------------------- */
      /* Berichtige Nachfolgernummern der eingefuegten Knoten              */
      /* ----------------------------------------------------------------- */

      for ( n = 0; n < ip_1->node_c; n++ ) {

         smaller = in_1[n].smaller; greater = in_1[n].greater;

         if( in_1[n].flag.smaller_extern )
            in_1[n].smaller = o_used[ smaller ];
         else
            in_1[n].smaller = nodes_used[ smaller ];

         if( in_1[n].flag.greater_extern )
            in_1[n].greater = o_used[ greater ];
         else
            in_1[n].greater = nodes_used[ greater ];
      }

      memory_free( o_used );
      memory_free( bp_used );
      memory_free( stack );
      memory_free( nodes_used );
   }

   table->page_buf = old_p;

   return( result );
}

/* ----------------------------------------------------------------------- */
/* Name       init_page_2                                                  */
/* Definition UINT init_page_2( TABLE *table, UINT root, UCHAR *buffer )   */
/* Funktion   Erzeuge 2. Seite des Spaltergebnisses.                       */
/* ----------------------------------------------------------------------- */
UINT init_page_2( TABLE *table, UINT root, UCHAR *buffer ) {

   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   UCHAR     *old_p    = table->page_buf;
   IPAGE     *ip_old   = (IPAGE*)(old_p + sizeof(PTYPE));
   INODE     *in_old   = (INODE*)(old_p + sizeof(PTYPE) + sizeof(IPAGE));
   BPOS      *bp_old   = (BPOS*)( old_p + ip_old->bpos);
   ULONG     *list_old = (ULONG*)(old_p + ip_old->olist);
   IPAGE     *ip_2     = (IPAGE*)(buffer + sizeof(PTYPE));
   INODE     *in_2     = (INODE*)(buffer + sizeof(PTYPE) + sizeof(IPAGE));
   BPOS      *bp_2;
   ULONG     *list;
   UINT     *bp_used, *o_used, *nodes_used;
   UINT      b_c = 0, b_s = 0, o_c = 0;
   UINT      n, b, r_top, r_next, r_next_extern,smaller,greater,node;
   UINT      pos;
   UINT     *stack;
   UINT      spos = 0;
   UINT      result = DB_ALL_OK;

   r_top = 0;

   while ( 1 ) {

      node = r_top + in_old[r_top].cnt - 1;

      smaller = in_old[ node ].smaller;
      greater = in_old[ node ].greater;

      if ( ( ! in_old[ node ].flag.smaller_extern ) && smaller == root ) {
         r_next = greater;
         r_next_extern = (in_old[ node ].flag.greater_extern)?(1):(0);
         break;
      }

      if ( ( ! in_old[ node ].flag.greater_extern ) && greater == root ) {
         r_next = smaller;
         r_next_extern = (in_old[ node ].flag.smaller_extern)?(1):(0);
         break;
      }

      if ( in_old[r_top].cnt == 0 ) {
         printf("Fatal error splitting index page (init_page_2)\n");
         printf("Process terminated\n");
         exit(-1);
      }

      r_top += in_old[r_top].cnt;
   }

   bp_used    = (UINT*)memory_alloc( ip_old->b_c * sizeof(UINT) );
   o_used     = (UINT*)memory_alloc( ip_old->o_c * sizeof(UINT) );
   nodes_used = (UINT*)memory_alloc( ip_old->node_c * sizeof(UINT));
   stack      = (UINT*)memory_alloc( ip_old->node_c * sizeof(UINT));

   if ( bp_used == NULL || o_used == NULL || nodes_used == NULL || stack == NULL ) {

      if( bp_used    != NULL ) memory_free( bp_used );
      if( o_used     != NULL ) memory_free( o_used );
      if( nodes_used != NULL ) memory_free( nodes_used );
      if( stack      != NULL ) memory_free( stack );

      result = DB_MEMORY_ERROR;

   } else {

      table->page_buf = buffer; newPage( table, IP );

      for( n = 0; n < ip_old->b_c; n++ )    bp_used[n]    = 0;
      for( n = 0; n < ip_old->o_c; n++ )    o_used[n]     = 0;
      for( n = 0; n < ip_old->node_c; n++ ) nodes_used[n] = 0;

      push( 0 );

      while ( ! stack_empty ) {

         pop( node );

         if ( node == r_top ) {
            if( r_next_extern ) o_used[r_next] = 1;
            else                push( r_next );
         } else {
            nodes_used[node] = 1;

            smaller = in_old[node].smaller; greater = in_old[node].greater;

            b = in_old[node].bnum;
            if( bp_used[b] == 0 ) {
               b_c++; b_s += bp_old[b].size;
            }
            bp_used[b] = 1;

            if( in_old[node].flag.smaller_extern ) o_used[ smaller ] = 1;
            else                              push( smaller );

            if( in_old[node].flag.greater_extern ) o_used[ greater ] = 1;
            else                         push( greater );
         }
      }

      /* ----------------------------------------------------------------- */
      /* Trage benoetigte Grenzen ein                                      */
      /* ----------------------------------------------------------------- */

      pos        = base->page_size;
      ip_2->bpos = base->page_size - b_s - b_c * sizeof(BPOS);
      ip_2->b_c  = b_c;
      bp_2       = (BPOS*)(buffer + ip_2->bpos);

      b_c = 0;

      for ( b = 0; b < ip_old->b_c; b++ ) {

         if( bp_used[b] > 0 ) {
            pos -= bp_old[b].size;

            bp_2[b_c].pos     = pos;
            bp_2[b_c].counter = 0;
            bp_2[b_c].size    = bp_old[b].size;
            bp_2[b_c].comp    = bp_old[b].comp;

            bp_used[b] = b_c;

            memory_copy( buffer+pos, old_p+bp_old[b].pos, bp_old[b].size );

            b_c++;
         }
      }

      /* ----------------------------------------------------------------- */
      /* Trage Positionen der Nachfolgeseiten ein                          */
      /* ----------------------------------------------------------------- */

      o_c = 0;

      for( n = 0; n < ip_old->o_c; n++ ) { if( o_used[n] > 0 ) o_c++; }

      ip_2->olist = ip_2->bpos - o_c * sizeof(ULONG);
      ip_2->o_c   = o_c;
      list        = (ULONG*)(buffer + ip_2->olist);

      o_c = 0;

      for ( n = 0; n < ip_old->o_c; n++ ) {
         if ( o_used[n] > 0 ) {
            list[o_c] = list_old[n];
            o_used[n] = o_c;
            o_c++;
         }
      }

      /* ----------------------------------------------------------------- */
      /* Trage Knoten des verbleibenden Graphen (ohne extrahierten G.) ein */
      /* ----------------------------------------------------------------- */

      node = 0;

      for ( n = 0; n < ip_old->node_c; n++ ) {
         if ( nodes_used[n] == 1 ) {
            nodes_used[n] = node;

            b = in_old[n].bnum; b = bp_used[b];
            bp_2[b].counter += 1;

            in_2[node].flag    = in_old[n].flag;
            in_2[node].bnum    = b;
            in_2[node].cnt     = in_old[n].cnt;
            in_2[node].smaller = in_old[n].smaller;
            in_2[node].greater = in_old[n].greater;

            node++;
         }
      }

      ip_2->node_c = node;

      /* ----------------------------------------------------------------- */
      /* Berichtige Nachfolgernummern der eingefuegten Knoten              */
      /* ----------------------------------------------------------------- */

      for ( n = 0; n < ip_2->node_c; n++ ) {
         smaller = in_2[n].smaller; greater = in_2[n].greater;

         if ( in_2[n].flag.smaller_extern ) {
            in_2[n].smaller = o_used[ smaller ];
         } else {
            if ( smaller == r_top ) {
               if ( r_next_extern ) {
                  in_2[n].flag.smaller_extern = 1;
                  in_2[n].smaller = o_used[r_next];
               } else {
                  in_2[n].smaller = nodes_used[r_next];
               }
            } else {
               in_2[n].smaller = nodes_used[ smaller ];
            }
         }

         if( in_2[n].flag.greater_extern ) {
            in_2[n].greater = o_used[ greater ];
         } else {
            if( greater == r_top ) {
               if( r_next_extern ) {
                  in_2[n].flag.greater_extern = 1;
                  in_2[n].greater = o_used[r_next];
               } else {
                  in_2[n].greater = nodes_used[r_next];
               }
            } else {
               in_2[n].greater = nodes_used[ greater ];
            }
         }
      }

      memory_free( o_used );
      memory_free( nodes_used );
      memory_free( bp_used );
      memory_free( stack );
   }

   table->page_buf = old_p;

   return( result );
}

/* ----------------------------------------------------------------------- */
/* Ende von SPLITI.C                                                       */
/* ----------------------------------------------------------------------- */
