/* ----------------------------------------------------------------------- */
/* FIO.C                                                                   */
/*                                                                         */
/* Funktion zum schnellen Einfuegen von Listen von Datensaetzen              */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>

#include "dbms.h"               /* Benutzerstrukturen/operationen          */
#include "dbmsdef.h"            /* Definitionen                            */
#include "dbmstype.h"           /* Interne Verwaltungsstrukturen           */
#include "dbmsprot.h"           /* Prototypen der internen Funktionen      */

inline int readc( FIO_PARM *fp );
char *(*GetMolBuff)(void);

/*------------------------------------------------------------------------ */
/* Name       dbms_fio                                                     */
/* Definition int dbms_fio( FIO_PARM *fiop );                              */
/* Prototyp   dbms.h                                                       */
/* Funktion                                                                */
/* ----------------------------------------------------------------------- */
int dbms_fio( FIO_PARM *fio ) {

   TABLE       *table = NULL;
   BASE_PAGE   *base;
   INDEX       *index_db;
   PAGEBUF      pagebuf[20];
   ULONG       *del_pages;
   static UINT  page_c = 0;
   UINT         actp = 0;
   UINT         last;
   UINT         n,i;
   UCHAR        spe;
   SINT         num;
   ULONG        offset = 0;
   ULONG        pages_to_del[20];
   UINT         del_cnt = 0;
   IPAGE       *ipage;
   ULONG       *list;
   int          error = DBMS_ALL_OK;

   error = init_table( &table, fio->name );   /* Tabelle initialisieren     */

   if ( ! error ) {

//      reclen = table->data_length;            /* Datensatzlaenge             */

      base = (BASE_PAGE*)(table->base);

      /* ------------------------------------------------------------------ */
      /* Pruefe Zugriffsberechtigung                                         */
      /* ------------------------------------------------------------------ */
      if ( table->change_access == '-' ) {
	 error = DBMS_ACCESS_DENIED;
/* !!! Wird derzeit nicht benutzt !!!
      } else if ( (table->user_rights).insert < table->change_access ) {
	 error = DBMS_ACCESS_DENIED;
*/
      }

      if ( ! error ) {

	 /* --------------------------------------------------------------- */
	 /* Hole Speicher fuer die Wurzelseite der Tabelle                   */
	 /* --------------------------------------------------------------- */
	 pagebuf[0].page = (char*)memory_alloc(base->page_size);
         if ( pagebuf[0].page == NULL ) {
            error = DBMS_MEMORY_ERROR;
         } else {

            table->operation = FIO;                     /* Aktive Operation */

            for ( n=0; n<20; n++ ) { pagebuf[n].offset = 0;   }
            del_pages = (ULONG*)(table->base + sizeof(BASE_PAGE) );
            page_c = 0;

            while ( ! error ) {
	       if (LogFile!=NULL) fprintf(LogFile,"\n");
               error = (fio->readb)( fio );
               (fio->count)++;
               if ( error ) {
                  if (LogFile!=NULL) fprintf(LogFile,"error (%u) reading record %lu\n", fio->count );
                  break;
               } else {
                  if (LogFile!=NULL) fprintf(LogFile,"insert record %lu\n", fio->count );
               }

	       del_cnt = 0;

	       table->data_buf = (UCHAR*)(fio->data);
               table->ditype   = (fio->dataaccess) & 63;

               /* --------------------------------------------------------- */
               /* Autozaehlerindizes bearbeiten                              */
               /* --------------------------------------------------------- */
               for ( i=0; i<table->index_c; i++ ) {

                  index_db = &(table->index_db[i]);
                  spe = index_db->special;
                  num = index_db->number;

                  if ( spe == AUTOINCREMENT || spe == AUTODECREMENT ) {

		     if ( base->variable ) {
			n = ((UINT*)(table->data_buf))[i+1];
		     } else {
			n = index_db->position;
		     }

		     *(long*)(table->data_buf+n) = (base->counter)[num];

		     if ( spe == AUTODECREMENT ) {
			((base->counter)[num])--;
		     } else {
			((base->counter)[num])++;
		     }
		  }
	       }

	       if ( base->variable ) {
		  if ( *((UINT*)(table->data_buf)) > base->max_data_length ) {
		     base->max_data_length = *((UINT*)(table->data_buf));
                  }
	       }

               /* -------------------------------------------------------- */
               /* Read pages into buffers                                  */
               /* ---------------------------------------------------------*/

               table->active = base->root_page;

               if (LogFile!=NULL) fprintf(LogFile,"read pages\n");

               for ( actp = 0; actp < base->level; actp++ ) {

                  if ( actp > page_c ) {
                     pagebuf[actp].page = (char*)memory_alloc(base->page_size);
                     if ( pagebuf[actp].page == NULL ) {
                        error = DBMS_MEMORY_ERROR;
                        break;
                     } else {
                        page_c++;
                     }
                     pagebuf[actp].offset = 0;
                  }

		  table->page_buf = (UCHAR*)(pagebuf[actp].page);

                  if ( table->active != pagebuf[actp].offset ) {

                     if ( pagebuf[actp].offset != 0 ) {
                        offset = table->active;
                        table->active = pagebuf[actp].offset;
                        error = write_page( table );
                        pagebuf[actp].offset = 0;
                        table->active = offset;
                     }

                     if ( ! error ) {
			error = read_page( table );
                     }

                     if ( error ) break;
                  } else {
                     if (LogFile!=NULL) fprintf(LogFile,"page in buffer %ld\n", table->active );
                  }

                  pages_to_del[del_cnt] = table->active;
                  del_cnt++;

                  if ( actp < base->level - 1 ) {
                     table->active = get_next_page( table, &(pagebuf[actp].next) );
                  }

                  if ( table->active == 0 ) {
                     error = DBMS_IPG_DEFECT;
                     break;
                  }
               }

               offset = 0;

               if ( ! error )   {

                  last = base->level - 1;

                  if (LogFile!=NULL) fprintf(LogFile,"FIO insert data\n");

                  error = insert_data( table, pagebuf, &last, pages_to_del, &del_cnt, &offset );

                  /* ----------------------------------------------------- */
                  /* Schreibe Seitenpuffer. Aendere dabei jeweils den      */
                  /* Eintrag in der Nachfolgerliste mit dem Offset der     */
                  /* nachfolgenden Seite (wie beim Einlesen vermerkt).     */
                  /* ----------------------------------------------------- */
                  if ( ! error ) {

                     while ( last > 0 ) {

                        last--;
			table->page_buf = (UCHAR*)(pagebuf[last].page);
                        ipage = (IPAGE*)(table->page_buf + sizeof(PTYPE));
                        list  = (ULONG*)(table->page_buf + ipage->olist);
                        *(list + pagebuf[last].next ) = offset;
                        offset = getFreeFilePosition( table );
                        pagebuf[last].offset = offset;
                     }

                     base->root_page = offset;
                     base->lognum++;
                     (base->data_c)++;

                     for( n = 0; n < del_cnt; n++ ) {
                        del_pages[base->del_c] = pages_to_del[n];
                        (base->del_c)++;
                     }
		  }
               }

               /* -------------------------------------------------------- */
               /* Behandlung der FIO-Spezialfaelle                          */
               /* -------------------------------------------------------- */
               if ( error == DBMS_DATA_IDENTICAL ) {
                  if (LogFile!=NULL) fprintf(LogFile,"Identical %lu\n", fio->count );
                  if ( fio->special & NO_BREAK_ON_IDENTICALS ) {
                     error = DBMS_ALL_OK;
                  }
               }

            }

            if ( error == DBMS_FIO_END ) { error = DBMS_ALL_OK; }

            for ( n=0; n<=page_c; n++ ) {

               if ( pagebuf[n].offset != 0 && ! error ) {
		  table->page_buf = (UCHAR*)(pagebuf[n].page);
                  table->active   = pagebuf[n].offset;
                  error = write_page( table );
               }
            }

	    for ( n=0; n<=page_c; n++ ) {
               memory_free( pagebuf[n].page );
            }
         }
      }

      if ( ! error ) { error = write_base_page( table ); }

      drop_table( table );
   }

   return( error );
}

/* ----------------------------------------------------------------------- */
/* Name       dbms_insert_f                                                */
/* Definition int dbms_insert_f( char *tblname, char *in, int msg );       */
/* Funktion   Einfuegen einer Liste von Datensaetzen aus der Datei 'in'      */
/*            in die Tabelle tblname. msg gibt an, ob Textmeldungen        */
/*            ausgegeben werden sollen oder nicht.                         */
/* Ergebnis   0 OK, sonst Fehlermeldung                                    */
/* ----------------------------------------------------------------------- */

int dbms_insert_c( char *tblname, char *(*ReadFn)(void))
{
  FIO_PARM fp;
  int error = DBMS_ALL_OK;

  GetMolBuff = ReadFn;
  LogFile = NULL; // fopen("log.txt","w");

  fp.name       = tblname;
  fp.data       = NULL;
  fp.readb      = readc;
  fp.dataaccess = 0;
  fp.count      = 0;
  fp.special    = NO_BREAK_ON_IDENTICALS;

  error = dbms_fio( &fp );

  return( error );
}

/* ----------------------------------------------------------------------- */
/* Name       readf                                                        */
/* Definition int readf( FIO_PARM *fp );                                   */
/* Funktion   Einlesen der Datensaetze von der Datei.                       */
/*            Wenn der letzte Datensatz gelesen wurde, wird als Ergebnis   */
/*            DBMS_FIO_END zurueckgegeben                                   */
/* ----------------------------------------------------------------------- */

inline int readc( FIO_PARM *fp )
{
   if ((fp->data=(GetMolBuff)())==NULL)
	return( DBMS_FIO_END );
   else
	return( DBMS_ALL_OK );
}

/* ----------------------------------------------------------------------- */
/* Ende von FIO.C                                                          */
/* ----------------------------------------------------------------------- */
