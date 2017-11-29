/* ----------------------------------------------------------------------- */
/* DELETEH.C                                                               */
/*                                                                         */
/* Hilfsfunktionen zum Loeschen von Daten. Werden von dbms_delete und      */
/* dbms_update                                                             */
/* benoetigt. Die Prototypen der Funktionen sind alle in dbmsprot.h        */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dbms.h"
#include "dbmsdef.h"
#include "dbmstype.h"
#include "dbmsprot.h"

#define R (table->index_c)

UINT delete_overlay( TABLE *table, UINT num, ULONG *pages_to_del, UINT *del_cnt, USER_RIGHTS *user );
UINT remove_dinfo( TABLE *table, UCHAR *datap, UCHAR *help, UINT num );
UINT remove_oinfo( UINT oinum, UCHAR *p );
UINT I_id( TABLE *table, UINT dinum, UCHAR *p );
UINT R_id( TABLE *table, UINT oinum, UCHAR *p );

/* ----------------------------------------------------------------------- */
/* Name:       delete_overlay                                              */
/* Definition: UINT delete_overlay( TABLE *table, UINT num,                */
/*                                   ULONG *pages_to_del, UINT *del_cnt ); */
/* Prototyp:   dbmsprot.h                                                  */
/* Funktion:   Loescht Datensatz von Overlayseite. Die Datenseite auf der  */
/*             sich der Verweis auf den zu betrachtenden Overlaybaum       */
/*             befindet muss in table->page_buffer stehen. num ist der     */
/*             Feldindex des Verweises. In pages_to_del werden alle nach   */
/*             erfolgreicher Operation zu loeschenden Seitenoffsets        */
/*             eingetragen. del_cnt wird dazu entsprechend erhoeht.        */
/* Ergebnis:   DB_ALL_OK, falls fehlerfrei ausgefuehrt, sonst Fehlermeldung.*/
/*             DB_MEMORY_ERROR                                             */
/*                Ursache: kein Speicher fuer Overlay-Puffer               */
/*                Woher:   delete_overlay                                  */
/*             DB_READ_ERROR                                               */
/*                Ursache: Seite kann nicht gelesen werden                 */
/*                Woher:   read_page                                       */
/*             DB_DATA_NOT_FOUND                                           */
/*                Ursache: der zu loeschende Datensatz wurde nicht gefunden*/
/*                Woher:   delete_overlay                                  */
/*             DB_WRITE_ERROR                                              */
/*                Ursache: Seite kann nicht geschrieben werden             */
/*                Woher:   write_page                                      */
/* Benutzte Funktionen:                                                    */
/*                memory_alloc           (MEMORY.H)                        */
/*                read_page              (FILE_OP.C)                       */
/*                overlay_identical      (DELETEH.C)                       */
/*                delete_tid             (DELETEH.C)                       */
/*                delete_odata           (DELETEH.C)                       */
/*                getFreeFilePosition    (BASEROUT.C)                      */
/*                write_page             (FILE_OP.C)                       */
/* ----------------------------------------------------------------------- */
UINT delete_overlay( TABLE *table, UINT num, ULONG *pages_to_del, UINT *del_cnt, USER_RIGHTS *user )
{
   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   UCHAR     *old_p = (UCHAR*)(table->page_buf);
   ULONG      old_a = table->active;
   DPAGE     *dp = (DPAGE*)(old_p + sizeof(PTYPE));
   DINFO     *di = (DINFO*)(old_p + dp->data );
   UCHAR     *overlay_0, *overlay_1;
   OINFO     *oi0, *oi1;
   OPAGE     *op0, *op1;
   UINT       n0, n1;
   UINT       found = 0;
   UINT       result = DB_ALL_OK;

   overlay_0 = (UCHAR*)memory_alloc( base->page_size );
   overlay_1 = (UCHAR*)memory_alloc( base->page_size );

   if( overlay_0 == NULL || overlay_1 == NULL )
   {
      if( overlay_0 != NULL ) memory_free( overlay_0 );
      if( overlay_1 != NULL ) memory_free( overlay_1 );
      result = DB_MEMORY_ERROR;
   }
   else
   {
      op0 = (OPAGE*)(overlay_0 + sizeof(PTYPE));
      op1 = (OPAGE*)(overlay_1 + sizeof(PTYPE));
      oi0 = (OINFO*)(overlay_0 + sizeof(PTYPE) + sizeof(OPAGE));
      oi1 = (OINFO*)(overlay_1 + sizeof(PTYPE) + sizeof(OPAGE));

      table->active   = di[num].lognum;
      table->page_buf = overlay_0;

      if( (result = read_page( table )) == DB_ALL_OK )
      {
         for( n0 = 0; n0 < op0->oic; n0++ )
         {
            if ( oi0[ n0 ].type & DIOVERLAY )
            {
               table->active   = oi0[ n0 ].lognum;
               table->page_buf = overlay_1;

               if( (result = read_page( table )) == DB_ALL_OK )
               {
                  for( n1 = 0; n1 < op1->oic; n1++ )
                  {
                     if ( R_id(table,n1,overlay_1)&&user->change >= ((oi1[n1].type & 56)>>3) )
                     {
                        table->ditype = oi1[n1].type;

                        pages_to_del[ *del_cnt ] = oi0[n0].lognum; (*del_cnt)++;
                        pages_to_del[ *del_cnt ] = di[num].lognum; (*del_cnt)++;

                        if( op1->oic > 1 )
                        {
                           remove_oinfo( n1, overlay_1 );

                           table->active = getFreeFilePosition( table );
                           if( (result = write_page( table )) == DB_ALL_OK )
                           {
                              oi0[ n0 ].lognum = table->active;
                              table->page_buf = overlay_0;
                              table->active = getFreeFilePosition( table );
                              result = write_page( table );
                              di[num].lognum = table->active;
                           }
                        }
                        else
                        {
                           if( op0->oic > 1 )
                           {
                              table->page_buf = overlay_0;
                              remove_oinfo( n0, overlay_0 );
                              table->active = getFreeFilePosition( table );
                              result = write_page( table );
                              di[num].lognum = table->active;
                           }
                           else
                           {
                              remove_dinfo( table, old_p, overlay_0, num );
                              result = DB_ALL_OK;
                           }
                        }

                        found = 1; break;
                     }
                  }
               }

               table->active   = di[num].lognum;
               table->page_buf = overlay_0;
            }
            else
            {
               if( R_id( table, n0, overlay_0 ) && user->change >= ((oi0[n0].type & 56)>>3) )
               {
                  table->ditype = oi0[n0].type;

                  pages_to_del[ *del_cnt ] = di[num].lognum; (*del_cnt)++;

                  if( op0->oic > 1 )
                  {
                     remove_oinfo( n0, overlay_0 );

                     table->active = getFreeFilePosition( table );
                     result = write_page( table );
                     di[ num ].lognum = table->active;
                  }
                  else
                  {
                     result = remove_dinfo( table, old_p, overlay_0, num );
                  }

                  found = 1; break;
               }
            }

            if( found ) break;
         }
      }

      memory_free( overlay_0 );
      memory_free( overlay_1 );
   }

   if( result == DB_ALL_OK && ! found ) result = DB_DATA_NOT_FOUND;

   table->page_buf = old_p;
   table->active   = old_a;

   return( result );
}

/* ----------------------------------------------------------------------- */
/* Name:       delete_data                                                 */
/* Definition: UINT delete_data( TABLE *table, UINT num );                 */
/* Prototyp:   dbmsprot.h                                                  */
/* Funktion:   Loescht Datensatz von Datenseite. Die Datenseite muss im    */
/*             Seitenpuffer der Tabelle stehen. num ist der Feldindex des  */
/*             zugehoerigen DATA_INFOs. Wenn der Datensatz nicht der Letzte*/
/*             im Feld der DATA_INFOs ist, wird zunaechst, der Datenblock  */
/*             oberhalb des zu loeschenden DS um dessen Groesse nach unten */
/*             kopiert. Dann die position-Werte aller DATA_INFOs nach num  */
/*             um diese Groesse geaendert. Dann wird der DATA_INFO-Block   */
/*             nach num um die Groesse eines DATA_INFOs nach oben kopiert. */
/*             Damit ist der Datensatz geloescht. Nun werden noch die PAGE-*/
/*             Parameter free_part, data_part und data_count geaendert.    */
/* Benutzte Funktionen: - memory_copy            (MEMORY.H)                */
/* ----------------------------------------------------------------------- */
UINT delete_data( TABLE *table, ULONG *offset, ULONG *pages_to_del, UINT *del_c, USER_RIGHTS *user )
{
   UCHAR     *old_p = table->page_buf;
   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   DPAGE     *dp = (DPAGE*)(old_p + sizeof(PTYPE));
   DNODE     *dn = (DNODE*)(old_p + sizeof(PTYPE) + sizeof(DPAGE));
   DINFO     *di = (DINFO*)(old_p + dp->data);
   UCHAR     *v1, *v2;
   UCHAR     *pNew;   
   UINT       rs1, rs2;
   UINT       i,x,z,x_next;
   UINT       result = DB_DATA_NOT_FOUND;

   if( dp->dic > 0 )
   {
      if( dp->dnc == 0 )
      {
         if( I_id( table, 0, old_p ) && user->change >= ((di[0].type & 56)>>3) )
         {
            if( di[0].type & DIOVERLAY )
            {
               result = delete_overlay( table, 0, pages_to_del, del_c, user );
            }
            else
            {
               v1 = proj( table, R, (UCHAR*)(table->data_buf) );
               v2 = proj( table, R, (UCHAR*)(table->page_buf + di[0].position) );
               rs1 = compSize( table, R, (UCHAR*)(table->data_buf) );
               rs2 = compSize( table, R, (UCHAR*)(table->page_buf + di[0].position) );

               if( rs1 == rs2 && memcmp( v1, v2, rs1 ) == 0 )
               {
                  table->ditype = di[0].type;
                  newPage( table, DP );
                  result = DB_ALL_OK;
               }
            }

            if( result == DB_ALL_OK )             /* schreibe Seite */       
            {
               table->active = getFreeFilePosition( table );
               result = write_page( table );
               *offset = table->active;
            }
         }
      }
      else
      {
         /*-----------------------------------------------------------------------------------
         Es befinden sich mehrere Datensaetze auf der Seite.
         -----------------------------------------------------------------------------------*/

         x = 0; 

         while( 1 )
         {
            i  = dn[x].comp;
            v1 = proj( table, i, table->data_buf);
            v2 = proj( table, i, table->page_buf + di[ dn[x].ds ].position );     

            switch( table->index_db[i].cmpfunc( v1, v2, table->index_db[i].length ) )
            {
               case  1:  x_next = _GREATER; break;

               case  0:  if( dn[x].mode & EQ_SMALLER ) x_next = _SMALLER;
                     else x_next = _GREATER;
                     break;

               case -1:  x_next = _SMALLER; break;
            }

            if( x_next == _SMALLER )
            {
               if( dn[x].mode & SMALLER_DS ) { z = dn[x].smaller; break; }
               else x = dn[x].smaller;
            }
            else
            {
               if( dn[x].mode & GREATER_DS ) { z = dn[x].greater; break; }
               else x = dn[x].greater;
            }
         }

         if( I_id( table, z, old_p ) && user->change >= ((di[z].type&56)>>3) )
         {
            if( di[z].type & DIOVERLAY )     /* versuche auf Overlayseite einzufuegen */
            {
               result = delete_overlay( table, z, pages_to_del, del_c, user );
            }
            else       /* Konstruiere neue Overlayseite */
            {
               v1 = proj( table, R, table->data_buf );
               v2 = proj( table, R, table->page_buf + di[z].position );
               rs1 = compSize( table, R, table->data_buf );
               rs2 = compSize( table, R, table->page_buf + di[z].position );

               if( rs1 == rs2 && memcmp( v1, v2, rs1 ) == 0 )
               {
                  table->ditype = di[z].type;

                  if( (pNew = (UCHAR*)memory_alloc( base->page_size )) == NULL )
                  {
                     result = DB_MEMORY_ERROR;
                  }
                  else
                  {
                     result = remove_dinfo( table, table->page_buf, pNew, z );                                           
                     memory_free( pNew );                                               
                     result = DB_ALL_OK;
                  }
               }
            }

            if( result == DB_ALL_OK )
            {
               table->active = getFreeFilePosition( table );
               *offset = table->active;
               result = write_page( table );
            }
         }
      }
   }

   return( result );
}



UINT remove_dinfo( TABLE *table, UCHAR *datap, UCHAR *help, UINT num ) {

   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   DPAGE     *dp = (DPAGE*)(datap + sizeof(PTYPE));
   DINFO     *di = (DINFO*)(datap + dp->data);
   UCHAR     *data;
   UINT       n;

   table->page_buf = help;
   newPage( table, DP );

   for( n = 0; n < dp->dic; n++ )
   {
      if( n != num )
      {
         data = datap + di[n].position;
         copyData( table, data, di + n );
      }
   }

   memory_copy( datap, help, base->page_size );
   table->page_buf = datap;

   return( DB_ALL_OK );
}

UINT remove_oinfo( UINT oinum, UCHAR *p ) {

   OPAGE *op = (OPAGE*)(p + sizeof(PTYPE));
   OINFO *oi = (OINFO*)(p + sizeof(PTYPE) + sizeof(OPAGE));
   UCHAR *source, *dest;
   UINT   n;

   source = p + op->data;
   dest   = source + oi[ oinum ].size;

   memory_copy( dest, source, oi[ oinum ].position - op->data );
   op->data += oi[ oinum ].size;

   for( n = oinum; n < op->oic; n++ ) oi[ n ].position += oi[oinum].size;

   source = p + sizeof(PTYPE) + sizeof(OPAGE) + (oinum+1) * sizeof(OINFO);
   dest   = source - sizeof(OINFO);
   memory_copy( dest, source, (op->oic - oinum - 1 ) * sizeof(OINFO) );
   (op->oic)--;
   op->free -= sizeof(OINFO);

   return( DB_ALL_OK );
}

UINT I_id( TABLE *table, UINT dinum, UCHAR *p ) {

   DPAGE *dp = (DPAGE*)(p + sizeof(PTYPE));
   DINFO *di = (DINFO*)(p + dp->data);
   UCHAR *v1, *v2;
   UINT   i;

   for( i = 0; i < table->index_c; i++ )
   {
      v1 = proj( table, i, table->data_buf);
      v2 = proj( table, i, table->page_buf + di[dinum].position );

      if( table->index_db[i].cmpfunc( v1, v2, table->index_db[i].length ) != 0 )
      {
         return( 0 );
      }
   }

   return( 1 );
}

UINT R_id( TABLE *table, UINT oinum, UCHAR *p ) {

   OINFO *oi = (OINFO*)(p + sizeof(PTYPE) + sizeof(OPAGE));
   UCHAR *v1, *v2;   
   UINT   rs1, rs2;

   v1  = proj( table, R, table->data_buf );
   v2  = p + oi[ oinum ].position;
   rs1 = compSize( table, R, table->data_buf );
   rs2 = oi[ oinum ].size;

   if( rs1 == rs2 )
   {
      if( memcmp( v1, v2, rs1 ) == 0 )
      {
         return( 1 );
      }
   }

   return( 0 );
}

/* ----------------------------------------------------------------------- */
/* Ende von DELETEH.C                                                      */
/* ----------------------------------------------------------------------- */
