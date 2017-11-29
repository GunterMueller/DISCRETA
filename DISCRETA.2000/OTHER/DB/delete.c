/* ----------------------------------------------------------------------- */
/* DELETE.C                                                                */
/*                                                                         */
/* Funktionen zum Loeschen in Tabellen                                      */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dbms.h"
#include "dbmsdef.h"
#include "dbmstype.h"
#include "dbmsprot.h"

/* ----------------------------------------------------------------------- */
/* Name       dbms_delete                                                  */
/* Definition int dbms_delete( DELETE_PARM *dp );                          */
/* Prototyp   dbms.h                                                       */
/* Funktion   Die Funktion loescht den Datensatz data in der Tabelle       */
/*            filename in dp. Hierzu wird die Tabelle geoeffnet und        */
/*            initialisiert. Dann werden die Seiten von der Wurzel bis     */
/*            zur Datenseite in den Puffer gelesen und der Datensatz       */
/*            gesucht. Kann er nicht gefunden werden, wird abgebrochen,    */
/*            sonst wird er geloescht.                                     */
/* Benutzte Funktionen                                                     */
/*            strcpy                                                       */
/*            init_table       (INIT.C)                                    */
/*            memory_alloc     (MEMORY.H)                                  */
/*            read_page        (FILE_OP.C)                                 */
/*            get_next_page    (INSERTH.C)                                 */
/*            delete_data      (DELETEH.C)                                 */
/*            save_base_page   (BASEROUT.C)                                */
/*            memory_free      (MEMORY.H)                                  */
/*            drop_table       (INIT.C)                                    */
/* ----------------------------------------------------------------------- */
int dbms_delete( DELETE_PARM *dp ) {

   TABLE      *table = NULL;
   BASE_PAGE  *base;
   PAGEBUF    *pagebuf;
   ULONG      *del_pages;
   UINT        page_c;
   UINT        last;
   UINT        n;
   ULONG       offset = 0;
   ULONG       *pages_to_del;
   UINT        del_cnt = 0;
   IPAGE      *ipage;
   ULONG      *list;
   int         error = DBMS_ALL_OK;

   USER_RIGHTS user = *(USER_RIGHTS*)(0);

   /* -------------------------------------------------------------------- */
   /* Initialisiere Tabelle und pruefe die Zugriffsrechte des Benutzers.   */
   /* Lege den Seitenpuffer an und eine Liste fuer Offset der nach erfolg- */
   /* reicher Operation zu loeschenden Seiten. Bei einem Fehler in diesem  */
   /* Bereich wird mit der jeweiligen Fehlernummer abgebrochen.            */
   /* -------------------------------------------------------------------- */

   error = init_table( &table, dp->name );

   if ( ! error ) {

      base = (BASE_PAGE*)(table->base);

      if ( ! error ) {

         table->data_buf      = (UCHAR*)(dp->data);

         pagebuf      = (PAGEBUF*)memory_alloc( base->level * sizeof(PAGEBUF) );  
         pages_to_del = (ULONG  *)memory_alloc( (base->level+2) * sizeof(ULONG) );

         if( pagebuf == NULL || pages_to_del == NULL ) {
            if( pagebuf      != NULL ) memory_free( pagebuf );
            if( pages_to_del != NULL ) memory_free( pages_to_del );
            error = DBMS_MEMORY_ERROR;
         } else {

            /* ----------------------------------------------------------- */
            /* Lese alle Seiten, von der Wurzel bis zur Datenseite in den  */
            /* Puffer ein. Dazu wird fuer jede Seite ein Pufferspeicher    */
            /* angelegt. Geht dies nicht, so wird mit Speicherfehler abge- */
            /* brochen. Es werden jeweils die Filepositionen der Seiten in */
            /* pages_to_del vermerkt, da diese nach erfolgreicher Operation*/
            /* als geloescht ein- getragen werden muessen. Ausserdem wird  */
            /* zu jeder Seite vermerkt, bei welchem Eintrag in olist die   */
            /* Seite verlassen wurde, wegen Aenderung beim Speichern.      */
            /* ----------------------------------------------------------- */

            table->active = base->root_page;
            last = base->level;

            for( page_c = 0; page_c < base->level; page_c++ ) {

               if( (pagebuf[page_c].page = (char*)memory_alloc(base->page_size)) == NULL ) {
                  error = DBMS_MEMORY_ERROR; break;
               }

               table->page_buf = (UCHAR*)(pagebuf[page_c].page);

               error = read_page( table );

               if ( error ) break;

               pages_to_del[del_cnt] = table->active; del_cnt++;

               if( page_c < base->level - 1 ) {
                  table->active = get_next_page( table, &(pagebuf[page_c].next) );
               }

               if ( table->active == 0 ) {
                  error = DBMS_IPG_DEFECT;
                  break;
               }
            }

            /* ----------------------------------------------------------- */
            /* Wenn alle Seiten von der Wurzel bis zur Datenseite erfolg-  */
            /* reich in den Puffer eingelesen werden konnten, wird nun mit */
            /* dem Loeschen begonnen. Dazu wird mit delete_data versucht    */
            /* den Datensatz auf der Datenseite bzw. evtl. auf einer       */
            /* Overlay-Seite zu loeschen. Wird hierbei festgestellt,        */
            /* dass ein identischer Datensatz nicht existiert, so wird     */
            /* mit Fehler abgebrochen.                                     */
            /* ----------------------------------------------------------- */

            if ( ! error ) {
               last--;

               error=delete_data(table,&offset,pages_to_del,&del_cnt,&user);

               /* -------------------------------------------------------- */
               /* Schreibe Seitenpuffer. Aendere dabei jeweils den Eintrag */
               /* in der Nachfolgerliste mit dem Offset der nachfolgenden  */
               /* Seite (wie beim Einlesen vermerkt).                      */
               /* -------------------------------------------------------- */

               if ( ! error ) {

                  while( last > 0 )
                  {
                     last--;
                     table->page_buf = (UCHAR*)(pagebuf[ last ].page);
                     ipage = (IPAGE*)(table->page_buf + sizeof(PTYPE));
                     list  = (ULONG*)(table->page_buf + ipage->olist);
                     list[ pagebuf[ last ].next ] = offset;
                     offset = getFreeFilePosition( table );
                     table->active = offset;
                     error = write_page( table );
                     if ( error ) break;
                  }
               }
            }

            while( page_c > 0 ) {  /* Gebe Seiten im Puffer frei */
               page_c--;
               memory_free( pagebuf[page_c].page );
            }
         }

         /* -------------------------------------------------------------- */
         /* Wenn der Datensatz erfolgreich geloescht wurde und die Seiten  */
         /* bis zur Wurzel geschrieben werden konnten, wird die Wurzel der */
         /* Tabelle auf die Neue gesetzt, die Datensatzanzahl um eins      */
         /* verringert, die zu loeschenden Seiten eingetragen und          */
         /* anschliessend die Basisseite geschrieben.                      */
         /* -------------------------------------------------------------- */

         if ( ! error ) {
            base->root_page = offset;
            (base->data_c)--;

            del_pages = (ULONG*)(table->base + sizeof(BASE_PAGE) );

            for( n = 0; n < del_cnt; n++ )
            {
               del_pages[ base->del_c ] = pages_to_del[n]; 
               (base->del_c)++;
            }
            error = write_base_page( table );
         }

         memory_free( pagebuf );   
         memory_free( pages_to_del );
      }
      drop_table( table );
   }

   return( error );
}

typedef struct _binrec {

   int       error;
   unsigned  reclen;
   unsigned  bsize;
   unsigned  size;
   char     *buffer;
   FILE     *file;

} BINREC;

int readrecord( BINREC *rec );

/* ----------------------------------------------------------------------- */
/* Name       dbms_remove_f                                                */
/* Definition int dbms_remove_f( char *tblname, char *in, int msg )        */
/* Funktion   Loeschen aller Datensaetze einer (binaeren) Datei im            */
/*            DB-Format                                                    */
/* ----------------------------------------------------------------------- */
int dbms_remove_f( char *tblname, char *in, int msg ) {

   DELETE_PARM dp;
   BINREC binrec;
   long   count = 0;
   int    error = DBMS_ALL_OK;

   binrec.error  = 0;
   binrec.reclen = 0;
   binrec.bsize  = 0;
   binrec.size   = 0;
   binrec.buffer = NULL;
   binrec.file   = fopen( in, "rb" );

   if ( binrec.file == NULL ) {
      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Fehler: kann %s nicht oeffnen\n", in ); break;
         default: printf("Error: unable to open %s\n", in );
      }
   } else {

      while ( readrecord( &binrec ) ) {
         count++;
         dp.name = tblname;
         dp.data = binrec.buffer;
         error = dbms_delete( &dp );
         printf("\rdeleted %ld", count );
         if ( error ) { break; }
      }
      printf("\n");
      if ( binrec.buffer != NULL ) memory_free( binrec.buffer );
      fclose( binrec.file );
   }

   return( error );
}


/* ----------------------------------------------------------------------- */
/* Name       readrecord                                                   */
/* Definition int readrecord( BINREC *rec );                               */
/* Funktion   Einlesen der zu loeschenden Datensaetze aus einer Datei im     */
/*            DB-Format (binaer)                                            */
/* ----------------------------------------------------------------------- */
int readrecord( BINREC *rec ) {

   if ( rec->reclen > 0 ) {

      if ( rec->bsize == 0 ) {
         rec->bsize = rec->reclen;
         rec->buffer = (char*)memory_alloc(rec->reclen);
         if ( rec->buffer == NULL ) { rec->error=1; return(0); }
      }
      if ( fread(rec->buffer,rec->reclen,1,rec->file)!=1) { rec->error=2; return(0); }

   } else {

      rec->size = 0;

      if ( fread( &(rec->size), sizeof(int), 1, rec->file ) != 1 ) {
         rec->error = 3; return( 0 );
      } else {
         if ( rec->size > rec->bsize ) {
            if ( rec->buffer != NULL ) { memory_free(rec->buffer); rec->buffer = NULL; }
            rec->buffer = (char*)memory_alloc(rec->size);
            if ( rec->buffer == NULL ) {
               rec->error = 1; return( 0 );
            }
         }
         *(unsigned*)(rec->buffer) = rec->size;
         if ( fread( (rec->buffer)+sizeof(int), (rec->size)-sizeof(int), 1, rec->file ) != 1 ) {
            rec->error = 2; return(0);
         }
      }
   }
   rec->error = 0;
   return( 1 );
}

/* ----------------------------------------------------------------------- */
/* Ende von DELETE.C                                                       */
/* ----------------------------------------------------------------------- */
