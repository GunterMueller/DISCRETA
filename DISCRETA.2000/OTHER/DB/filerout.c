/* ----------------------------------------------------------------------- */
/*                                                                         */
/* FILEROUT.C                                                              */
/*                                                                         */
/* Funktionen fuer die Dateioperationen. Es wurden ausschliesslich Standard  */
/* C-Funktionen benutzt, damit die Portabilitaet gewaehrleistet ist.         */
/* Bei einigen Systemen koennte die Funktion access (in file_exists) nicht  */
/* vorhanden sein, so dass diese nachgebildet werden muss.                   */
/*                                                                         */
/* Benutzte Funktionen:                                                    */
/*                                                                         */
/* file_exists                                                             */
/* create_table                                                            */
/* open_table                                                              */
/* close_table                                                             */
/* delete_table                                                            */
/* read_page                                                               */
/* write_page                                                              */
/* read_base_page                                                          */
/* write_base_page                                                         */
/*                                                                         */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <unistd.h>
/* Aenderung wolfb: statt #include <io.h> */

#include "dbmsdef.h"            /* Definitionen                            */
#include "dbmstype.h"           /* Interne Verwaltungsstrukturen           */


/* ----------------------------------------------------------------------- */
/* Name       file_exists                                                  */
/* Definition int file_exists( char *filename );                           */
/* Prototyp   dbmsprot.h                                                   */
/* Funktion   Wenn Datei mit Namen filename existiert, ist das Ergebnis 1  */
/*            sonst 0.                                                     */
/* Benutzt    access (IO.H)                                                */
/* ----------------------------------------------------------------------- */
int file_exists( char *filename ) {

   if ( access( filename, 0 ) == -1 ) return( 0 );
   return( 1 );
}

/* ----------------------------------------------------------------------- */
/* Name       create_table                                                 */
/* Definition int create_table( TABLE *table );                            */
/* Prototyp   dbmsprot.h                                                   */
/* Funktion   Erzeugt Datei mit Namen table->filename. Wenn dies nicht     */
/*            gelingt, so wird DBMS_CREATE_ERROR zurueckgegeben, sonst      */
/*            DBMS_ALL_OK.                                                 */
/* Benutzt    fopen (STDIO.H)                                              */
/* ----------------------------------------------------------------------- */
int create_table( TABLE *table ) {

   if ( (table->fileptr = fopen( table->filename, "wb" )) == NULL ) {
      return( DBMS_CREATE_ERROR );
   }

   return( DBMS_ALL_OK );
}

/* ----------------------------------------------------------------------- */
/* Name       open_table                                                   */
/* Definition int open_table( TABLE *table );                              */
/* Prototyp   dbmsprot.h                                                   */
/* Funktion   oeffnet die Tabelle mit Namen table->filename zur Bearbeitung */
/*            Wenn dies fehlschlaegt, so ist das Ergebnis DBMS_OPEN_ERROR,  */
/*            andernfalls DBMS_ALL_OK.                                     */
/* Benutzt    fopen (STDIO.H)                                              */
/* ----------------------------------------------------------------------- */
int open_table( TABLE *table ) {

   if ( (table->fileptr = fopen( table->filename, "r+b" )) == NULL ) {
      return( DBMS_OPEN_ERROR );
   }
   return( DBMS_ALL_OK );
}

/* ------------------------------------------------------------------------ */
/* Name       close_table                                                   */
/* Definition int close_table( TABLE *table );                              */
/* Prototyp   dbmsprot.h                                                    */
/* Funktion   Schliesst eine mit open_table geoeffnete Tabellendatei wieder.  */
/*            Im Fehlerfall wird DBMS_CLOSE_ERROR, sonst DBMS_ALL_OK        */
/*            zurueckgegeben.                                                */
/* Benutzt    fclose      (STDIO.H)                                         */
/* ------------------------------------------------------------------------ */
int close_table( TABLE *table ) {

   if ( fclose( table->fileptr ) != 0 ) return( DBMS_CLOSE_ERROR );
   table->fileptr = NULL;
   return( DBMS_ALL_OK );
}

/* ------------------------------------------------------------------------ */
/* Name       delete_table                                                  */
/* Definition int delete_table( TABLE *table );                             */
/* Prototyp   dbmsprot.h                                                    */
/* Funktion   Loescht die Tabelle table->filename. Falls dies fehlschlug,    */
/*            wird DBMS_DELETE_ERROR sonst DBMS_ALL_OK zurueckgegeben.       */
/* Benutzt    remove (STDIO.H>                                              */
/* ------------------------------------------------------------------------ */
int delete_table( TABLE *table ) {

   if ( remove( table->filename ) != 0 ) {
      return( DBMS_DELETE_ERROR );
   }
   return( DBMS_ALL_OK );
}

/* ------------------------------------------------------------------------ */
/* Name       read_page                                                     */
/* Definition int read_page( TABLE *table );                                */
/* Prototype  dbmsprot.h                                                    */
/* Funktion   Liest eine Seite von Position table->active (relativ zum      */
/*            Dateianfang) aus der Datei der geoeffneten Tabelle. Die Groesse  */
/*            der Seite ist base->page_size. Die Seite wird in den Seiten-  */
/*            speicher table->page_buffer geladen. Bei einem Lesefehler     */
/*            ist das Ergebnis DBMS_READ_ERROR, sonst DBMS_ALL_OK.          */
/* Benutzt    fseek (STDIO.h)                                               */
/*            fread (STDIO.H)                                               */
/* ------------------------------------------------------------------------ */
int read_page( TABLE *table ) {

   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   int        error = DBMS_ALL_OK;

	Log("read_page(");
   if ( fseek( table->fileptr, table->active, FILE_BEGIN ) != 0 ) {
      if (LogFile!=NULL) fprintf(LogFile,"read page error fseek %ld\n", table->active );
      error = DBMS_READ_ERROR;
   } else {
      if ( fread( table->page_buf, base->page_size, 1, table->fileptr )!=1 ) {
	 if (LogFile!=NULL) fprintf(LogFile,"read page error fread %ld\n", base->page_size );
         error = DBMS_READ_ERROR;
      }
   }
	Log(")\n");
   return( error );
}

/* ------------------------------------------------------------------------ */
/* Name       write_page                                                    */
/* Definition int write_page( TABLE *table );                               */
/* Prototyp   dbmsprot.h                                                    */
/* Funktion   Schreibt die Seite im Seitenspeicher auf die Position         */
/*            table->active (rel. zum Dateianfang). Im Fehlerfall ist das   */
/*            Ergebniss DBMS_WRITE_ERROR, sonst DBMS_ALL_OK.                */
/* Benutzt    fseek   (STDIO.H)                                             */
/*            fwrite  (STDIO.H)                                             */
/*            fflush  (STDIO.H)                                             */
/* ------------------------------------------------------------------------ */
int write_page( TABLE *table ) {

   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   int        error = DBMS_ALL_OK;

	Log("write_page(");
   if ( fseek( table->fileptr, table->active, FILE_BEGIN ) != 0 ) {
      error = DBMS_WRITE_ERROR;
   } else {
      if ( fwrite( table->page_buf, base->page_size, 1, table->fileptr )!=1 ){
         error = DBMS_WRITE_ERROR;
      }
   }
   fflush( table->fileptr );

	Log(")\n");
   return( error );
}

/* ------------------------------------------------------------------------ */
/* Name       read_base_page                                                */
/* Definition int read_base_page( TABLE *table );                           */
/* Prototyp   dbmsprot.h                                                    */
/* Funktion   Liest die Basisseite der Tabelle ein. Diese Seite befindet    */
/*            auf Pos. 0 der Datei. Im Fehlerfall wird DBMS_READ_ERROR,     */
/*            sonst DBMS_ALL_OK zurueckgegeben.                              */
/* Benutzt    fseek  (STDIO.h)                                              */
/*            fread  (STDIO.H)                                              */
/* ------------------------------------------------------------------------ */
int read_base_page( TABLE *table ) {

   int error = DBMS_ALL_OK;

	Log("read_base_page(");
   if ( fseek( table->fileptr, 0, FILE_BEGIN ) != 0 ) {
      error = DBMS_READ_ERROR;
   } else {
      if ( fread( table->base, BASE_PAGE_SIZE, 1, table->fileptr ) != 1 ) {
         error = DBMS_READ_ERROR;
      }
   }

	Log(")\n");
   return( error );
}

/* ------------------------------------------------------------------------ */
/* Name:      write_base_page                                               */
/* Definition int write_base_page( TABLE *table );                          */
/* Prototyp   dbmsprot.h                                                    */
/* Funktion   Schreibt die Basisseite der Tabelle (Pos. 0). Im Fehlerfall   */
/*            ist das Ergebnis DBMS_WRITE_ERROR, sonst DBMS_ALL_OK.         */
/* Benutzt    fseek    (STDIO.H)                                            */
/*            fwrite   (STDIO.H)                                            */
/*            fflush   (STDIO.h)                                            */
/* ------------------------------------------------------------------------ */
int write_base_page( TABLE *table ) {

   int error = DB_ALL_OK;

	Log("write_base_page(");
   if ( fseek( table->fileptr, 0, FILE_BEGIN ) != 0 ) {
      error = DBMS_WRITE_ERROR;
   } else {
      if ( fwrite( table->base, BASE_PAGE_SIZE, 1, table->fileptr ) != 1 ) {
         error = DBMS_WRITE_ERROR;
      }
   }
   fflush( table->fileptr );

	Log(")\n");
   return( error );
}

/* ----------------------------------------------------------------------- */
/* Ende von FILEROUT.C                                                     */
/* ----------------------------------------------------------------------- */
