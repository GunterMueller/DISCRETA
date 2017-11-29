/* ----------------------------------------------------------------------- */
/*                                                                         */
/* CREATE.C - Funktionen zum Erzeugen von Tabellen                         */
/*                                                                         */
/* (c) J. Mueller Juni 1993  -  Version 1.22                                */
/*                                                                         */
/* Funktionen:                                                             */
/*                                                                         */
/* - dbmsCreate     Erzeugen mittels CREATE_PARM                           */
/* - dbmsCreateF    Erzeugen via Definitionsdatei                          */
/*                                                                         */
/* Aenderungen:                                                             */
/*                                                                         */
/* 10.12.91 Juergen Mueller (jm)                                             */
/* 10.02.92 (jm) set_table_structure durch class_settings ersetzt.         */
/* 28.02.92 (jm) Moeglichkeit geschaffen die Seitengroesse beim               */
/*               Erzeugen anzugeben.                                       */
/* 16.03.92 (jm) type-Wahl bei new_page geschaffen.                        */
/* 17.03.92 (jm) memset in new_base_page und new_page entfernt da          */
/*               fehlerhaft. Ursache konnte nicht gefunden werden.         */
/* 18.03.92 (jm) Abfrage der Zugriffsrechte eingebaut, Funktion            */
/*               check_page_size zum Pruefen der Seitengroesse.               */
/*               Funktionen new_base_page und new_page nach BASEROUT.C     */
/* 30.03.92 (jm) Tabellenname wird nicht mehr doppelt, in table und        */
/*               auf der Basiseite gefuehrt. tablename in table ist nur     */
/*               noch Zeiger auf Basisseite.                               */
/* 28.04.92 (jm) Aenderung von TABLE und BASE_PAGE, new_base_page           */
/*               wieder hier.                                              */
/* 07.12.92 (jm) user_rights in TABLE eingebaut. Abfrage auch hier.        */
/*                                                                         */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dbms.h"       /* Strukturen/Prototypen fuer Benutzeroperationen   */
#include "dbmsdef.h"    /* Definitionen                                    */
#include "dbmstype.h"   /* Interne Verwaltungsstrukturen                   */
#include "dbmsprot.h"   /* Prototypen fuer interne Funktionen               */

/* ----------------------------------------------------------------------- */
/*                                                                         */
/* dbmsCreate                                                              */
/*                                                                         */
/* Erzeuge Tabelle gemaess Parameter (CREATE_PARM). Dieser enthaelt alle      */
/* Informationen zur Beschreibung der Tabelle.                             */
/* Es wird sowohl ein TDF (Definitionsdatei binaer) als auch die Tabellen-  */
/* datei selbst erzeugt (TBL). Diese enthaelt die Basisseite und eine       */
/* leere Datenseite.                                                       */
/*                                                                         */
/* Benutzte Funktionen:                                                    */
/*                                                                         */
/* file_exists           filerout                                          */
/* create_table          filerout                                          */
/* write_base_page       filerout                                          */
/* write_page            filerout                                          */
/* close_table           filerout                                          */
/* remove                stdio                                             */
/* new_page              baserout                                          */
/* makefname             baserout                                          */
/* memory_alloc          Makro                                             */
/* memory_free           Makro                                             */
/*                                                                         */
/* Ergebnis:                                                               */
/*                                                                         */
/* 0 OK, sonst groesser 0                                                    */
/*                                                                         */
/* ----------------------------------------------------------------------- */
int dbmsCreate( CREATE_PARM *cp ) {

   TABLE     *table;
   TDF        tdf;
   FILE      *file;
   FILE      *log;
   BASE_PAGE *base;
   UINT       n, m, l1 = 0, l2 = 0;
   CHAR      *tblpath, *tdfpath;
   int        error = DBMS_ALL_OK;

   log = fopen( "dbcreate.log", "w" );

   if ( log != NULL ) {
      fputs("Protokoll - dbmsCreate\n", log);
      fprintf(log,"  Tabelle:          %s\n", cp->filename );
      fprintf(log,"  Indexanzahl:      %u\n", cp->index_count );
      fprintf(log,"  Lesezugriff:      %c\n", cp->read_access );
      fprintf(log,"  Aenderungszugriff: %c\n", cp->change_access );
      fprintf(log,"  Seitengroesse:      %u\n", cp->page_size );
      fprintf(log,"  Satzlaenge:        %u\n", cp->record_length );
      fprintf(log,"  Infoposition:     %u\n", cp->info_pos );
      fprintf(log,"  Infolaenge:        %u\n", cp->info_size );
   }

   /* hole Speicher fuer Tabellenstruktur */
   table = (TABLE*)memory_alloc( sizeof(TABLE) );
   if ( table == NULL ) {
      if ( log != NULL ) {
         fputs("- Kein Speicher fuer die Tabelle\n", log );
         fputs("-- Abbruch\n", log );
         fclose( log );
      }
      return DBMS_MEMORY_ERROR;
   }
   if ( log != NULL ) {
      fprintf(log,"+ Speicher fuer Tabelle ok [%p]\n", table );
   }

   /* setze Zeiger auf die Basisseite der Tabelle */
   base = (BASE_PAGE*)(table->base);

   /* pruefe, ob Pfad fuer Tabellendatei eingestellt ist   */
   /* l1 ist Laenge des Pfades (falls vorhanden, sonst 0) */
   tblpath = getenv("DBMS_TBL_PATH");
   if ( tblpath != NULL ) {
      l1 = strlen(tblpath);
      if ( tblpath[l1-1] != DIRSEP ) l1++;
      if ( log != NULL ) fprintf(log,"  TBL-Pfad: %s\n", tblpath);
   } else {
      if ( log != NULL ) fputs("  Kein TBL-Pfad\n", log );
   }

   /* pruefe, ob Pfad fuer Tabellendefinitionsdatei eingestellt ist */
   /* l2 ist Laenge des Pfades (falls vorhanden, sonst 0) */
   tdfpath = getenv("DBMS_TDF_PATH");
   if ( tdfpath != NULL ) {
      l2 = strlen(tdfpath);
      if ( tdfpath[l2-1] != DIRSEP ) l2++;
      if ( log != NULL ) fprintf(log,"  TDF-Pfad: %s\n", tblpath);
   } else {
      if ( log != NULL ) fputs("  Kein TBL-Pfad\n", log );
   }

   /* Hole Speicher fuer Puffer fuer die Dateinamen      */
   /* dieser muss sowohl den TDF- wie auch den TBL-Pfad */
   /* beinhalten koennen, d.h.                          */
   /* max(l1,l2) + Laenge Tabellenname + Extension + 1  */
   if ( l1 < l2 ) l1 = l2;
   if ( strlen(TBLEXT) > strlen(TDFEXT) )
      l1 += strlen(TBLEXT);
   else
      l1 += strlen(TDFEXT);
   l1 += strlen(cp->filename) + 1; /* + 1 for '\0' */
   if ( log != NULL ) {
      fprintf(log,"  Hole Speicher fuer Dateinamen: %u Bytes\n", l1 );
   }
   table->filename = (CHAR*)memory_alloc( l1 );
   if ( table->filename == NULL ) {
      if ( log != NULL ) {
         fputs("- Kein Speicher fuer Dateinamen\n",log);
         fprintf(log, "  Gebe Speicher fuer Tabelle frei [%p]\n", table );
         fputs("-- Abbruch\n",log);
         fclose( log );
      }
      memory_free( table );
      return( DBMS_MEMORY_ERROR );
   }
   if ( log != NULL ) fprintf(log,"+ Speicher ok [%p]\n", table->filename );

   /* pruefe, ob gleichnamige Tabelle schon existiert */
   /* erzeuge dazu den vollstaendigen Tabellennamen   */
   /* mit Pfad, Name und Erweiterung                 */
   makefname( table->filename, tblpath, cp->filename, TBLEXT );
   if ( file_exists( table->filename ) ) {
      if ( log != NULL ) {
         fprintf(log,"- Tabelle %s existiert bereits\n", table->filename );
         fprintf(log,"  Gebe Speicher fuer Dateinamen frei [%p]\n",
                                                        table->filename );
         fprintf(log, "  Gebe Speicher fuer Tabelle frei [%p]\n", table );
         fputs("-- Abbruch\n", log);
         fclose( log );
      }
      memory_free( table->filename );
      memory_free( table );
      return( DBMS_FNAME_NOT_UNIQUE );
   }

   if ( log != NULL ) fputs("  Pruefe Eindeutigkeit der Indexnamen\n", log );

   /* pruefe, ob die Indexnamen eindeutig vergeben wurden */
   for ( n=0; n<cp->index_count-1; n++ ) {
      if ( log != NULL ) fprintf(log, "  %2u %s\n", n+1, cp->name[n] );
      for ( m=n+1; m<cp->index_count; m++ ) {
         if ( strncmp( cp->name[n], cp->name[m], 10 )==0 ) {
            if ( log != NULL ) {
               fputs("- Indexname schon vergeben\n", log );
               fprintf(log,"  Gebe Speicher fuer Dateinamen frei [%p]\n",
                                                    table->filename );
               fprintf(log, "  Gebe Speicher fuer Tabelle frei [%p]\n",
                                                              table );
               fputs("-- Abbruch\n", log);
               fclose( log );
            }
            memory_free( table->filename );
            memory_free( table );
            return ( DBMS_INAME_NOT_UNIQUE );
         }
      }
   }
   if ( log != NULL ) fputs("+ Indexnamen korrekt\n", log );

   /* initialisiere die Tabellendefinition (TDF) */
   tdf.status        = 0;
   tdf.read_access   = cp->read_access;
   tdf.change_access = cp->change_access;
   tdf.psize         = cp->page_size;
   tdf.reclen        = cp->record_length;
   tdf.icnt          = cp->index_count;
   tdf.rpos          = cp->info_pos;
   tdf.rsize         = cp->info_size;

   for ( n=0; n<cp->index_count; n++ ) {
      tdf.pos[n]   = cp->pos[n];
      tdf.elem[n]  = cp->elem[n];
      tdf.type[n]  = cp->type[n];
      tdf.flags[n] = cp->flags[n];
      strncpy( tdf.name[n], cp->name[n], 10 );
   }

   /* erzeuge Dateinamen fuer die TDF-Datei und lege die  */
   /* Datei an. (Binaerdatei, wird zum Schreiben geoeffnet */
   makefname( table->filename, tdfpath, cp->filename, TDFEXT );
   file = fopen( table->filename, "wb" );

   if (log!=NULL) fprintf(log,"  Erzeuge TDF-Datei %s\n", table->filename );

   /* Erzeugung von TDF schiefgegangen */
   if ( file == NULL ) {
      if ( log != NULL ) {
         fputs("- Erzeugung TDF fehlgeschlagen\n", log );
         fprintf(log,"  Gebe Speicher fuer Dateinamen frei [%p]\n",
                                                  table->filename );
         fprintf(log, "  Gebe Speicher fuer Tabelle frei [%p]\n",
                                                           table );
         fputs("-- Abbruch\n", log);
         fclose( log );
      }
      memory_free( table->filename );
      memory_free( table );
      return( DBMS_TDF_WRITE_ERROR );
   }

   if ( log != NULL ) fputs("+ TDF-Datei erzeugt\n", log);

   /* schreibe die initialisierte TDF-Struktur in die */
   /* neu angelegte TDF-Datei                         */
   if ( fwrite( &tdf, sizeof(TDF), 1, file ) != 1 ) {
      if ( log != NULL ) {
         fputs("- Schreiben von TDF fehlgeschlagen\n", log );
         fputs("  Schliesse und loesche TDF-Datei\n", log );
         fprintf(log,"  Gebe Speicher fuer Dateinamen frei [%p]\n",
                                                  table->filename );
         fprintf(log, "  Gebe Speicher fuer Tabelle frei [%p]\n",
                                                           table );
         fputs("-- Abbruch\n", log);
         fclose( log );
      }
      fclose( file );
      remove( table->filename );
      memory_free( table->filename );
      memory_free( table );
      return( DBMS_TDF_WRITE_ERROR );
   }

   /* schliesse die Datei wieder, damit ist die Tabellen- */
   /* definitionsdatei erfolgreich angelegt worden       */
   fclose( file );

   if ( log != NULL ) fputs("+ TDF geschrieben\n", log);

   /* nun folgt die Erzeugung der Tabellendatei */

   /* setze alle Bytes der Basisseite auf 0 */
   for ( n=0; n<BASE_PAGE_SIZE; n++ ) {
      ((char*)(table->base))[n] = 0;
   }

   /* Trage Copyrightvermerk auf der Basisseite ein */
   strcpy( base->crjm, "MARVIN V1.22 (c) J. Mueller 1994" );

   /* initialisiere die Basiseite */

   /* dann die physikalischen Tabellendaten */
   base->page_size       = cp->page_size;
   base->max_data_length = cp->record_length;
   base->variable        = (cp->record_length>0)?(0):(1);
   base->version         = 0;
   base->file_end        = BASE_PAGE_SIZE + base->page_size;
   base->level           = 1;
   base->root_page       = BASE_PAGE_SIZE;
   base->ovl_root        = 0;
   base->lognum          = 1;
   base->del_c           = 0;
   base->data_c          = 0;

   /* setze alle internen Zaehler auf 0 */
   for ( n=0; n<BASCNT; n++ ) {
      base->counter[n] = 0;
   }

   /* hole Speicher fuer den Seitenpuffer zur */
   /* Erzeugung einer leeren Datenseite      */
   table->page_buf = (unsigned char*)memory_alloc( base->page_size );
   if ( table->page_buf == NULL ) {
      if ( log != NULL ) {
         fprintf(log,"- Kein Speicher fuer Seitenpuffer (%u Bytes)\n",
                                                      base->page_size );
         fputs("  Loesche TDF-Datei\n", log );
         fprintf(log,"  Gebe Speicher fuer Dateinamen frei [%p]\n",
                                                  table->filename );
         fprintf(log, "  Gebe Speicher fuer Tabelle frei [%p]\n", table );
         fputs("-- Abbruch\n", log);
         fclose( log );
      }
      remove( table->filename );
      memory_free( table->filename );
      memory_free( table );
      return( DBMS_MEMORY_ERROR );
   }

   if ( log != NULL ) {
      fprintf(log,"+ Speicher fuer Seitenpuffer ok, %u Bytes [%p]\n",
                                      base->page_size, table->page_buf );
   }

   /* erzeuge Dateinamen der neuen Tabellendatei und */
   /* erzeuge die Tabelle                            */
   makefname( table->filename, tblpath, cp->filename, TBLEXT );
   error = create_table( table );
   if ( error ) {
      if ( log != NULL ) {
         fprintf(log,"- Fehler beim Erzeugen der Tabelle %s\n",
                                                  table->filename );
         fputs("  Loesche TDF-Datei\n", log );
         fprintf(log,"  Gebe Speicher fuer Dateinamen frei [%p]\n",
                                                  table->filename );
         fprintf(log,"  Gebe Speicher fuer Seitenpuffer frei [%p]\n",
                                                  table->page_buf );
         fprintf(log,"  Gebe Speicher fuer Tabelle frei [%p]\n", table );
         fputs("-- Abbruch\n", log);
         fclose( log );
      }
      makefname( table->filename, tblpath, cp->filename, TDFEXT );
      remove( table->filename );
      memory_free( table->filename );
      memory_free( table->page_buf );
      memory_free( table );
      return( error );
   }

   if ( log != NULL ) {
      fprintf(log,"+ Tabellendatei %s angelegt\n", table->filename );
   }

   /* schreibe nun die Basisseite  */
   error = write_base_page( table );
   if ( error ) {
      if ( log != NULL ) {
         fprintf(log,"- Fehler beim Schreiben der Basisseite\n" );
         fputs("  Schliesse und loesche Tabellendatei\n", log );
         fputs("  Loesche TDF-Datei\n", log );
         fprintf(log,"  Gebe Speicher fuer Dateinamen frei [%p]\n",
                                                  table->filename );
         fprintf(log,"  Gebe Speicher fuer Seitenpuffer frei [%p]\n",
                                                  table->page_buf );
         fprintf(log,"  Gebe Speicher fuer Tabelle frei [%p]\n", table );
         fputs("-- Abbruch\n", log);
         fclose( log );
      }
      close_table( table );
      remove( table->filename );
      makefname( table->filename, tblpath, cp->filename, TDFEXT );
      remove( table->filename );
      memory_free( table->filename );
      memory_free( table->page_buf );
      memory_free( table );
      return( error );
   }

   if ( log != NULL ) fputs("+ Basisseite geschrieben\n", log);

   /* initialisiere den Seitenpuffer als leere Datenseite */
   /* und schreibe diese unmittelbar nach der Basiseite   */
   /* in die Tabellendatei, schliesse anschliessend die     */
   /* Datei. Damit ist die Tabelle angelegt worden        */
   newPage( table, DP );
   table->active = BASE_PAGE_SIZE;
   error = write_page( table );
   close_table( table );

   /* falls ein Fehler auftrat, entferne die neu angelegten */
   /* TDF- und TBL-Dateien wieder                           */
   if ( error ) {
      if ( log != NULL ) {
         fprintf(log,"- Fehler beim Schreiben der Datenseite\n" );
         fputs("  Schliesse und loesche Tabellendatei\n", log );
         fputs("  Loesche TDF-Datei\n", log );
         fprintf(log,"  Gebe Speicher fuer Dateinamen frei [%p]\n",
                                                  table->filename );
         fprintf(log,"  Gebe Speicher fuer Seitenpuffer frei [%p]\n",
                                                  table->page_buf );
         fprintf(log,"  Gebe Speicher fuer Tabelle frei [%p]\n", table );
         fputs("-- Abbruch\n", log);
         fclose( log );
      }
      remove( table->filename );
      makefname( table->filename, tblpath, cp->filename, TDFEXT );
      remove( table->filename );
      memory_free( table->filename );
      memory_free( table->page_buf );
      memory_free( table );
      return( error );
   }

   if ( log != NULL ) {
      fputs("+ Datenseite geschrieben\n", log);
      fputs("  Schliesse Tabellendatei\n", log);
      fprintf(log,"  Gebe Speicher fuer Dateinamen frei [%p]\n",
                                                 table->filename );
      fprintf(log,"  Gebe Speicher fuer Seitenpuffer frei [%p]\n",
                                                 table->page_buf );
      fprintf(log,"  Gebe Speicher fuer Tabelle frei [%p]\n", table );
      fputs("++ dbmsCreate erfolgreich beendet\n", log );
      fclose( log );
   }

   memory_free( table->filename );
   memory_free( table->page_buf );
   memory_free( table );

   return( DBMS_ALL_OK );
}

/* ------------------------------------------------------------------------ */
/*                                                                          */
/* dbmsCreateF - Erzeugen von Tabellen via Definitionsdatei                 */
/*                                                                          */
/* Die Definitionsdatei wird zuerst im aktuellen Verzeichnis, und dann,     */
/* falls gesetzt, im Verzeichnis DBMS_DEF_PATH (Umgebungsvariable).         */
/*                                                                          */
/* ------------------------------------------------------------------------ */
int dbmsCreateF( char *defname ) {

   CREATE_PARM *cp = (CREATE_PARM*)memory_alloc(sizeof(CREATE_PARM));
   FILE        *infile = NULL;
   int          n = 0;
   CHAR         line1[512], line2[512];
   int          error;

   /* falls kein Speicher fuer den CREATE_PARM vorhanden war */
   /* breche mit Fehler ab                                  */
   if ( cp == NULL ) {
      return DBMS_MEMORY_ERROR;
   }

   /* erzeuge Dateinamen der Definitionsdatei und suche im aktuellen */
   /* Pfad nach dieser Datei, falls diese dort nicht gefunden werden */
   /* suche (falls angegeben) im DBMS_DEF_PATH, falls gefunden,      */
   /* oeffne die Datei (ASCII) zum Lesen.                             */
   makefname( line1, NULL, defname, DEFEXT );
   infile = fopen( line1, "r" );
   if ( infile == NULL ) {
      if ( getenv( "DB_DEF_PATH" ) != NULL ) {
         makefname( line1, getenv("DB_DEF_PATH"), defname, DEFEXT );
         infile = fopen( line1, "r" );
      }
   }
   if ( infile == NULL ) {
      memory_free( cp );
      return( DBMS_DEF_OPEN_ERROR );
   }

   /* lese die 1. Zeile; diese beinhaltet den Namen der zu erzeugenden */
   /* Tabelle, setze den Namen im CREATE_PARM auf die gelesene Zeile   */
   if ( fgets( line1, 512, infile ) == NULL ) {
      fclose( infile );
      memory_free( cp );
      return( DBMS_DEF_READ_ERROR );
   }
   line1[strlen(line1)-1] = '\0';
   cp->filename = line1;

   /* lese die 2. Zeile; diese hat folgende Gestalt              */
   /* <Zugriff> <Seitengroesse> <Satzlaenge> <Infopos.> <Infogroesse> */
   if ( fgets( line2, 512, infile ) == NULL ) {
      fclose( infile );
      memory_free( cp );
      return( DBMS_DEF_READ_ERROR );
   }

   /* Zugriff: rc wobei r,c aus '0',...'9', fuer c auch '-' */
   /* r Zugriffsbeschraenkung fuer Lesen,                    */
   /* c Zugriffsbeschraenkung fuer Aenderungen                */
   /* '0',...,'9' bedeutet, dass der Benutzer mind. den     */
   /* jeweiligen Wert als Berechtigung haben muss,          */
   /* '-' bei den Aenderungsrechten bedeutet, dass fuer die   */
   /* Tabelle keine Aenderungen zugelassen sind.            */
   if ( ! (line2[0]>='0' && line2[0]<='9') ) {
      fclose( infile );
      memory_free( cp );
      return( DBMS_DEF_WRONG_ACC );
   }
   if ( ! (line2[1]=='-'||(line2[1]>='0'&&line2[1]<='9')) ) {
      fclose( infile );
      memory_free( cp );
      return( DBMS_DEF_WRONG_ACC );
   }
   cp->read_access   = line2[0];
   cp->change_access = line2[1];

   /* lese die Seitengroesse (DOS max. 64K) sowie Satzlaenge */
   /* Position und Groesse des Infoteils                    */
   sscanf( line2+3, "%u %u %u %u", &(cp->page_size),
              &(cp->record_length), &(cp->info_pos), &(cp->info_size) );

   /* lese die Indexdefinitionen                              */
   /* gelesen werden alle Zeilen bis eine Zeile nicht mit 'I' */
   /* oder ';' beginnt oder die maximale Indexzahl (INDEXC)   */
   /* erreicht wurde.                                         */
   /* Zeilen mit ';' am Anfang werden als Kommentarzeilen     */
   /* uebersprungen und nicht ausgewertet                      */
   n = 0;
   while ( fgets( line2, 512, infile ) != NULL ) {

      if ( line2[0] == ';' ) continue;
      if ( line2[0] != 'I' || line2[1] != ' ' ) break;
      if ( n == INDEXC ) break;
      sscanf( line2+2, "%u %s %u %u %u",
                 &(cp->type[n]),
                 cp->name[n],
                 &(cp->elem[n]),
                 &(cp->pos[n]),
                 &(cp->flags[n]) );
      n++;
   }

   /* schliesse die Definitionsdatei, es wurden alle zum */
   /* Erzeugen benoetigten Informationen gelesen         */
   fclose( infile );

   /* es muss mindestens ein Index definiert sein */
   /* sonst erfolgt Abbruch mit Fehler           */
   if ( n==0 ) {
      memory_free( cp );
      return( DBMS_DEF_NO_INDEX );
   }

   /* setze die Indexanzahl in CREATE_PARM und rufe dbmsCreate */
   /* zur Erzeugung der Tabelle auf. Kehre mit dem Ergebnis    */
   /* dieser Funktion zurueck. Falls dieses DBMS_ALL_OK war     */
   /* so konnte die Tabelle erfolgreich erzeugt werden         */
   cp->index_count = n;
   error = dbmsCreate( cp );
   memory_free( cp );
   return ( error );
}

/* ----------------------------------------------------------------------- */
/* Ende von CREATE.C                                                       */
/* ----------------------------------------------------------------------- */






