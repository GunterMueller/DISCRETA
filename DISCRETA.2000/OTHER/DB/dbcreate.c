/* ----------------------------------------------------------------------- */
/*                                                                         */
/* DBCREATE.C - Erzeugen einer Tabelle von der Kommandozeile               */
/*                                                                         */
/* (c) J. Mueller Juni 1993  -  Version 1.22                                */
/*                                                                         */
/* Funktionen:                                                             */
/* printInfo        Ausgabe der Information zum Aufruf der Operation       */
/* printErrorMsg    Ausgabe des Operationsergebnisses als Textmeldung      */
/* main             Hauptfunktion fuer Tabellenerzeugung                    */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dbms.h"       /* Strukturen/Prototypen fuer Benutzeroperationen   */
#include "dbmserr.h"     /* Fehlercodes der Benutzeroperationen             */

FILE *LogFile;

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
      printf("| DBCREATE: Erzeugen einer neuen DBMS-Tabelle          |\n");
      printf("|                                                      |\n");
      printf("| Aufruf: dbcreate <Definitionsdatei> [-nomsg]         |\n");
      printf("|------------------------------------------------------|\n");
      printf("| Aufbau einer Definitionsdatei                        |\n");
      printf("|                                                      |\n");
      printf("| Zeile 1: Dateiname der neuen Tabelle                 |\n");
      printf("| Zeile 2: <acc> <psize> <reclen> <infopos> <infosize> |\n");
      printf("| Zeile 3 ff.:                                         |\n");
      printf("|          'Index' <type> <name> <elem> <pos> <flags>  |\n");
      printf("--------------------------------------------------------\n");

   } else {

      printf("--------------------------------------------------------\n");
      printf("| DBMS    Version 1.22 - June 1993    Juergen Mueller  |\n");
      printf("|------------------------------------------------------|\n");
      printf("| DBCREATE: Create a new table                         |\n");
      printf("|                                                      |\n");
      printf("| Syntax: dbcreate <definiton file> [-nomsg]           |\n");
      printf("|------------------------------------------------------|\n");
      printf("| Structure of Definitionsfiles                        |\n");
      printf("|                                                      |\n");
      printf("| Line 1: Filename of table to create                  |\n");
      printf("| Line 2: <acc> <psize> <reclen> <restpos> <restsize>  |\n");
      printf("| Line 3: etc.                                         |\n");
      printf("|          'Index' <type> <name> <elem> <pos> <flags>  |\n");
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
/* main - Hauptfunktion zum Erzeugen von Tabellen                          */
/*                                                                         */
/* Aufruf mit dbcreate <definion file> [-nomsg]                            */
/* Durch Setzen der Umgebungsvaribale DBMS_LANGUAGE kann man zwischen      */
/* der deutschen und englischen Sprachumgebung waehlen.                     */
/*                                                                         */
/* ----------------------------------------------------------------------- */
int main( int argc, char *argv[] ) {

   char  *h;
   int   error = 0, language = 0;

   LogFile = NULL;

   if ( ( h=getenv("DB_LANGUAGE") ) != NULL ) {
      if ( strcmp(h,"GERMAN")==0 || strcmp(h,"DEUTSCH")==0 )
         language = 1;
   }

   if ( argc < 2 ) {

      printInfo(language);
      error = DBMS_NO_DEF_FILE;

   } else {

      if ( ! ( argc > 2 && strcmp( argv[2],"-nomsg" )==0 ) ) {
         printf("** DBMS Version 1.22 Juergen Mueller\n");
         printf("   DBCREATE %s\n", argv[1] );
         error = dbmsCreateF( argv[1] );
         printErrorMsg( error, language );
      } else {
         error = dbmsCreateF( argv[1] );
      }

   }

   return( error );
}

/* ----------------------------------------------------------------------- */
/* Ende von DBCREATE.C                                                     */
/* ----------------------------------------------------------------------- */
