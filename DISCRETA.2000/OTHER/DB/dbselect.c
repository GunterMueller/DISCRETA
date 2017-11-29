/* ----------------------------------------------------------------------- */
/* DBSELECT.C                                                              */
/*                                                                         */
/* Suchfunktion der Datenbank. Die Eingabe der Suchgrenzen erfolgt ueber    */
/* eine Textdatei. Dort kann die zu durchsuchende Tabelle, die Zieldatei   */
/* fuer die gefundenen Datensaetze und die Grenzen des Suchraumes angegeben  */
/* werden.                                                                 */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dbms.h"       /* Benutzerstrukturen/operationen                  */
#include "dbmserr.h"     /* Fehlercodes der Benutzeroperationen             */

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
      printf("| DBSELECT: Suchen in einer Tabelle                    |\n");
      printf("|                                                      |\n");
      printf("| Aufruf: dbselect <Suchdefinitionsdatei>              |\n");
      printf("|------------------------------------------------------|\n");
      printf("| Aufbau einer Definitionsdatei                        |\n");
      printf("|                                                      |\n");
      printf("| Zeile 1: from <Tabellenname>                         |\n");
      printf("| Zeile 2: to <Zieldatei fuer Ergebnis>                 |\n");
      printf("| Zeile 3 ff.:                                         |\n");
      printf("|          where <Indexname>                           |\n");
      printf("|          >[=]  <Wert der Untergrenze>                |\n");
      printf("|          <[=]  <Wert der Obergrenze>                 |\n");
      printf("--------------------------------------------------------\n");

   } else {

      printf("--------------------------------------------------------\n");
      printf("| DBMS    Version 1.22 - June 1993    Juergen Mueller  |\n");
      printf("|------------------------------------------------------|\n");
      printf("| DBSELECT: Searching in tables                        |\n");
      printf("|                                                      |\n");
      printf("| Syntax: dbselect <search definition file>            |\n");
      printf("|------------------------------------------------------|\n");
      printf("| Structure of Definitionsfiles                        |\n");
      printf("|                                                      |\n");
      printf("| Line 1: from <table name>                            |\n");
      printf("| Line 2: to <file for search result>                  |\n");
      printf("| Line 3: etc.                                         |\n");
      printf("|         where <name of index>                        |\n");
      printf("|          >[=]  <value of lower bound>                |\n");
      printf("|          <[=]  <value of upper bound>                |\n");
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
/* main - Hauptfunktion zum Durchsuchen von Tabellen                       */
/*                                                                         */
/* Aufruf mit dbselect <definion file>                                     */
/* Durch Setzen der Umgebungsvaribale DBMS_LANGUAGE kann man zwischen      */
/* der deutschen und englischen Sprachumgebung waehlen.                     */
/*                                                                         */
/* ----------------------------------------------------------------------- */
void main( int argc, char *argv[] ) {

   char  *h;
   int   error = 0, language = 100, msg;

   log = fopen( "dbselect.log", "w" );

   if ( (h=getenv( "DBMS_LANGUAGE" )) != NULL ) {
      if ( strcmp( h, "GERMAN") == 0 || strcmp( h, "DEUTSCH" ) == 0 ) {
         language = 1;
      }
   }

   if ( argc < 2 ) {

      printInfo( language );

   } else {

      msg = language;
      if ( argc > 2 && strcmp( argv[2], "-nomsg" ) == 0 ) msg = 0;
      error = dbms_select_f( argv[1], msg );
   }

   if ( error ) {
      printErrorMsg( error, language );
   }

   exit ( error );
}

/* ----------------------------------------------------------------------- */
/* Ende von DBSELECT.C                                                     */
/* ----------------------------------------------------------------------- */
