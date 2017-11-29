/* ----------------------------------------------------------------------- */
/* DBINSERT.C                                                              */
/*                                                                         */
/* Einfuegen einer Liste von Datensaetzen aus einer Datei in eine            */
/* bestehende Tabelle                                                      */
/*                                                                         */
/* Funktionen                                                              */
/* printInfo                                                               */
/* printErrorMsg                                                           */
/* main                                                                    */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "dbms.h"       /* Benutzerstrukturen/operationen                  */
#include "dbmserr.h"     /* Fehlercodes der Operationen                     */

FILE *LogFile; /* Aenderung wolfb: zusaetzliche Zeile! */
FILE *log;

/* ----------------------------------------------------------------------- */
/*                                                                         */
/* printInfo -  Ausgabe der Aufrufinformation                              */
/*                                                                         */
/* Deutsch/Englisch je nach gewaehlter Sprache.                             */
/*                                                                         */
/* ----------------------------------------------------------------------- */
void printInfo( int language ) {

   if ( language==1 ) {

      printf("--------------------------------------------------------\n");
      printf("| DBMS    Version 1.22 - Juni 1993    Juergen Mueller  |\n");
      printf("|------------------------------------------------------|\n");
      printf("| DBINSERT: Einfuegen in bestehende Tabelle            |\n");
      printf("|                                                      |\n");
      printf("| Aufruf: dbinsert from:<dat1> to:<dat2>               |\n");
      printf("|------------------------------------------------------|\n");
      printf("| dat1 Datei mit Daten im binaeren DB-Format            |\n");
      printf("| dat2 Name der Tabelle in die eingefuegt wird         |\n");
      printf("--------------------------------------------------------\n");

   } else {

      printf("--------------------------------------------------------\n");
      printf("| DBMS    Version 1.22 - June 1993    Juergen Mueller  |\n");
      printf("|------------------------------------------------------|\n");
      printf("| DBINSERT: Insert record into existing table          |\n");
      printf("|                                                      |\n");
      printf("| Syntax: dbinsert from:<dat1> to:<dat2>               |\n");
      printf("|------------------------------------------------------|\n");
      printf("| dat1 File with record in binary DB-format            |\n");
      printf("| dat2 Name of table                                   |\n");
      printf("--------------------------------------------------------\n");

   }

   return;
}

/* ----------------------------------------------------------------------- */
/*                                                                         */
/* printErrorMsg - Ausgabe der Fehlerinterpretation nach Tabellenerzeugung */
/*                                                                         */
/* Deutsch/Englisch je nach gewaehlter Sprache.                             */
/*                                                                         */
/* ----------------------------------------------------------------------- */
void printErrorMsg( int error, int language ) {

   if ( language==1 ) {  /* Deutsch */

      switch ( error ) {
         case DBMS_ALL_OK:
            printf("++ Tabelle erzeugt\n");
            break;
         case DBMS_DEF_OPEN_ERROR:
            printf("-- Definitionsdatei nicht gefunden\n");
            break;
         case DBMS_DEF_READ_ERROR:
            printf("-- kann Definitionsdatei nicht lesen\n");
            break;
         case DBMS_DEF_WRONG_ACC:
            printf("-- falsche Zugriffsdefinition\n");
            break;
         case DBMS_DEF_NO_INDEX:
            printf("-- kein Index definiert\n");
            break;
         case DBMS_MEMORY_ERROR:
            printf("-- Speicherfehler\n");
            break;
         case DBMS_FNAME_NOT_UNIQUE:
            printf("-- Datei gleichen Namens existiert schon\n");
            break;
         case DBMS_TDF_WRITE_ERROR:
            printf("-- kann TDF nicht erzeugen\n");
            break;
         case DBMS_CREATE_ERROR:
            printf("-- kann Tabelle nicht erzeugen\n");
            break;
         case DBMS_WRITE_ERROR:
            printf("-- kann Seiten nicht schreiben\n");
            break;
         default:
            printf("-- unbekannte Fehlernummer (%d)\n");
      }

   } else { /* Default Englisch */

      switch ( error ) {
         case DBMS_ALL_OK:
            printf("++ Table created\n");
            break;
         case DBMS_DEF_OPEN_ERROR:
            printf("-- unable to open definition file\n");
            break;
         case DBMS_DEF_READ_ERROR:
            printf("-- unable to read definition file\n");
            break;
         case DBMS_DEF_WRONG_ACC:
            printf("-- incorrect access definition\n");
            break;
         case DBMS_DEF_NO_INDEX:
            printf("-- no index defined\n");
            break;
         case DBMS_MEMORY_ERROR:
            printf("-- not enough memory available\n");
            break;
         case DBMS_FNAME_NOT_UNIQUE:
            printf("-- file with that name already exists\n");
            break;
         case DBMS_TDF_WRITE_ERROR:
            printf("-- unable to create TDF\n");
            break;
         case DBMS_CREATE_ERROR:
            printf("-- unable to create table file\n");
            break;
         case DBMS_WRITE_ERROR:
            printf("-- unable to write pages\n");
            break;
         default:
            printf("-- unknown return value (%d)\n");
      }
   }
}


/* ----------------------------------------------------------------------- */
/*                                                                         */
/* main - Hauptfunktion zum Einfuegen von Datensaetzen                       */
/*                                                                         */
/* Aufruf mit dbinsert from: <dat1> to: <dat2>                             */
/* Durch Setzen der Umgebungsvaribale DBMS_LANGUAGE kann man zwischen      */
/* der deutschen und englischen Sprachumgebung waehlen.                     */
/*                                                                         */
/* ----------------------------------------------------------------------- */
int main( int argc, char *argv[] ) {

   float  t0,t1;
   char  *h;
   int    error = 0;
   int    language = 0;

   log = NULL;
   LogFile = NULL;
  LogFile = fopen( "dbinsert.log", "w" );
    if ( LogFile!=NULL ) { fprintf(LogFile,"** DBINSERT\n"); }

   if ( ( h=getenv("DB_LANGUAGE") ) != NULL ) {
      if ( strcmp(h,"GERMAN")==0 || strcmp(h,"DEUTSCH")==0 )
         language = 1;
   }

   if ( argc < 3 ) {
      printInfo(language);
      if (LogFile !=NULL) { fprintf(LogFile,"- Zuwenig Parameter\n"); }
   } else {
      printf("** DBMS - DBINSERT (c) Juergen Mueller 1993\n");
      t0 = clock();
      if ( strncmp( argv[1], "from:", 5 ) == 0 ) {
         if ( strncmp( argv[2], "to:", 3 ) == 0 ) {
            if (LogFile!=NULL) {
               fprintf(LogFile,"von  %s\n", argv[1]+5 ); fflush(LogFile);
               fprintf(LogFile,"nach %s\n", argv[2]+3 ); fflush(LogFile);
            }
            error = dbms_insert_f( argv[2]+3, argv[1]+5, 1 );
         } else {
            if (LogFile!=NULL) {
	fprintf(LogFile,"- Tabelle nicht angegeben\n"); fflush(LogFile);
	}
         }
      } else if ( strncmp( argv[1], "to:", 3 ) == 0 ) {
         if ( strncmp( argv[2], "from:", 5 ) == 0 ) {
            if (LogFile!=NULL) {
               fprintf(LogFile,"von  %s\n", argv[2]+5 ); fflush(LogFile);
               fprintf(LogFile,"nach %s\n", argv[1]+3 ); fflush(LogFile);
            }
            error = dbms_insert_f( argv[1]+3, argv[2]+5, 1 );
         } else {
            if (LogFile!=NULL) {
	fprintf(LogFile,"- Quelle nicht angegeben\n"); fflush(LogFile);
	}
         }
      } else {
         printInfo(language);
      }
      if ( error ) {
         if (LogFile !=NULL) {
	fprintf(LogFile,"- Fehlernummer %d\n", error ); fflush(LogFile);
	}
         printErrorMsg( error, language );
      } else {
         t1 = clock();
         if ( language == 1 ) {
            printf("Benoetigte Zeit: %.0f s\n", (t1-t0) / CLOCKS_PER_SEC );
         } else {
            printf("Time required: %.0f s\n", (t1-t0) /  CLOCKS_PER_SEC );
         }
      }
   }
  if (LogFile)
	fclose(LogFile);
}

/* ----------------------------------------------------------------------- */
/* Ende von DBINSERT.C                                                     */
/* ----------------------------------------------------------------------- */
