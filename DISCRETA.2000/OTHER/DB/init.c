/* ----------------------------------------------------------------------- */
/* INIT.C                                                                  */
/*                                                                         */
/* Funktionen zum Initialisieren und Schliessen einer Tabelle               */
/* Dieses Modul wird fuer alle Operationen ausser DBCREATE benoetigt          */
/*                                                                         */
/*                                                                         */
/* Funktionen:                                                             */
/* init_table                                                              */
/* drop_table                                                              */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dbms.h"        /* Benutzer Funktionen/Strukturen                 */
#include "dbmsdef.h"     /* Definitionen                                   */
#include "dbmstype.h"    /* Interne Verwaltungsstrukturen                  */
#include "dbmsprot.h"    /* Prototypen der internen Funktionen             */

extern INDEX    index_db[];
extern UINT     index_count;
extern UINT     data_length;
extern TRIGHTS  table_rights;
extern CHAR     read_access;
extern CHAR     change_access;

TABLE        table;

/* ----------------------------------------------------------------------- */
/* Name       init_table                                                   */
/* Definition int init_table( TABLE *table, char *name );                  */
/* Prototype  dbmsprot.h                                                   */
/* Funktion   oeffnet eine Tabelle zur Bearbeitung. Der Basisseitenpuffer   */
/*            wird allociert und die Datei zum Lesen und Schreiben ge-     */
/*            oeffnet. Ausserdem wird die TABLE-Struktur table gemaess den     */
/*            Werten der Basisseite initialisiert.                         */
/*            Die Funktion bekommt hierzu die Struktur table uebergeben.    */
/*            Hierin muss der Filename der Tabelle angegeben sein.          */
/* Benutzt    open_table                                                   */
/*            read_base_page                                               */
/*            get_definition                                               */
/* ----------------------------------------------------------------------- */
int init_table( TABLE **tbl, char *name ) {

   unsigned   l1 = 0, l2 = 0;
   char      *tblpath, *tdfpath;
   int        error = DBMS_ALL_OK;

   if (LogFile!=NULL) {
	fprintf(LogFile,"INIT '%s'\n", name);
	fflush(LogFile);
	}

   tblpath = getenv("DB_TBL_PATH"); tdfpath = getenv("DB_TDF_PATH");
   if ( tblpath != NULL ) { l1 = strlen(tblpath); if ( tblpath[l1-1] != DIRSEP ) l1++; }
   if ( tdfpath != NULL ) { l2 = strlen(tdfpath); if ( tdfpath[l2-1] != DIRSEP ) l2++;	}
   if ( l1 < l2 ) l1 = l2;
   if ( strlen(TBLEXT) > strlen(TDFEXT) ) l1+=strlen(TBLEXT); else l1+=strlen(TDFEXT);
   l1 += strlen(name) + 1; /* + 1 for '\0' */

   table.filename = (char*)memory_alloc( l1 );

   if ( table.filename == NULL ) {
      if (LogFile!=NULL) {
	fprintf(LogFile,"- Speicherfehler (Filename)\n" );
	fflush(LogFile);
	}
      error = DBMS_MEMORY_ERROR;
   } else {
      makefname( table.filename, NULL, name, TDFEXT );
      if (LogFile!=NULL) {
	fprintf(LogFile,"* get_definition [1] %s\n", table.filename );
	fflush(LogFile);
	}
      error = get_definition( table.filename );
      if ( error ) {
         makefname( table.filename, tdfpath, name, TDFEXT );
         if (LogFile!=NULL) {
	fprintf(LogFile,"* get_definition [2] %s\n", table.filename );
  	fflush(LogFile);
	}
       error = get_definition( table.filename );
      }
   }

   if ( error ) {
      if (LogFile!=NULL) {
	fprintf(LogFile,"- TDF kann nicht geoeffnet werden [%d]\n",error );
	fflush(LogFile);
	}
   } else {

      if (LogFile!=NULL) {
	fprintf(LogFile,"+ TDF gelesen\n");
	fflush(LogFile);
	}

      table.user_rights.insert = 0;
      table.user_rights.search = 0;
      table.user_rights.change = 0;

      makefname( table.filename, NULL, name, TBLEXT );
      if (LogFile!=NULL) {
	fprintf(LogFile,"* open_table [1] %s\n", table.filename );
	fflush(LogFile);
	}
      error = open_table( &table );
      if ( error ) {
         makefname( table.filename, tblpath, name, TBLEXT );
         if (LogFile!=NULL) fprintf(LogFile,"* open_table [2] %s\n", table.filename );
         error = open_table( &table );
      }

      if ( error ) {
         if (LogFile!=NULL) {
	fprintf(LogFile,"- TBL kann nicht geoeffnet werden [%d]\n", error );
	fflush(LogFile);
	}
      } else {

         if (LogFile!=NULL) {
            fprintf(LogFile,"- TBL geoeffnet\n" );
            fprintf(LogFile,"* read_base_page\n" );
	fflush(LogFile);
         }

         error = read_base_page( &table );

         if ( error ) {

            if (LogFile!=NULL) {
	fprintf(LogFile,"- read_base_page [%d]\n", error );
	fflush(LogFile);
	}

         } else {

            if (LogFile!=NULL) {
	fprintf(LogFile,"+ read_base_page ok\n" );
	fflush(LogFile);
	}

            if ( table.read_access == '-' && table.change_access == '-' ) {

               error = DBMS_TBL_LOCKED;
               if (LogFile!=NULL) {
	fprintf(LogFile,"- kein Zugriff\n" );
	fflush(LogFile);
	}

            } else {

               table.data_length = data_length;
               table.page_buf    = NULL;
               table.active      = 0;
               table.data_buf    = NULL;
               table.counter     = 0;
               table.read_access = read_access;
               table.change_access = change_access;

               table.index_db   = index_db;
               table.index_c = index_count;
            }
         }
         if ( error ) {
            if (LogFile!=NULL) {
	fprintf(LogFile,"schliesse Tabelle\n");
	fflush(LogFile);
	}
            close_table( &table );
         } else {
	    if (LogFile!=NULL) {
		fprintf(LogFile,"INIT OK\n");
		fflush(LogFile);
		}
            (*tbl) = &table;
         }
      }
   }

   if ( table.filename != NULL ) memory_free( table.filename );

   return( error );
}

/* ----------------------------------------------------------------------- */
/* Name       drop_table                                                   */
/* Definition void drop_table( TABLE *table );                             */
/* Prototyp   dbmsprot.h                                                   */
/* Funktion   Beendet die Bearbeitung einer Tabelle                        */
/* Benutzt    close_table  (FILEROUT.C)                                    */
/* ----------------------------------------------------------------------- */
void drop_table( TABLE *table ) {
   close_table( table );
   return;
}

/* ----------------------------------------------------------------------- */
/* Ende von INIT.C                                                         */
/* ----------------------------------------------------------------------- */
