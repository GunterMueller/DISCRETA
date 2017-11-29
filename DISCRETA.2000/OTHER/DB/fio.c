/* ----------------------------------------------------------------------- */
/* FIO.C                                                                   */
/*                                                                         */
/* Funktion zum schnellen Einfuegen von Listen von Datensaetzen            */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* #include <conio.h> */

#include "dbms.h"               /* Benutzerstrukturen/operationen          */
#include "dbmsdef.h"            /* Definitionen                            */
#include "dbmstype.h"           /* Interne Verwaltungsstrukturen           */
#include "dbmsprot.h"           /* Prototypen der internen Funktionen      */

FILE     *infile = NULL;
	/* in dbms_insert_f() geoeffnet */
uint      reclen = 0;
	/* benoetigt in readf() */
char     *buffer = NULL;
	/* allociert in readf(); 
	 * bei fester Datensatzlaenge wird 
	 * derselbe Puffer immer wieder verwendet, 
	 * ansonsten wird er jedesmal neu allociert */
uint      bsize = 0;
	/* size of buffer;
	 * wenn ungleich Null, so ist buffer allociert. */


int readf( FIO_PARM *fp );
	/* Einlesefunktion;
	 * Funktionspointer nach FIO_PARM->readb */

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
	/* Index der letzten allocierten Seite in pagebuf[] */
   UINT         actp = 0;
   UINT         last;
   UINT         n,i;
   UCHAR        spe;
   SINT         num;
   ULONG        offset = 0;

	/* Liste der Seiten des Zugriffsweges;
	 * werden im Anschluss an den Zugriff freigegeben, 
	 * d.h. in der Base Page notiert in der Liste der 
	 * freigegebenen Seiten. */
   ULONG        pages_to_del[20];
   UINT         del_cnt = 0;
			/* Laenge pages_to_del[] */

   IPAGE       *ipage;
   ULONG       *list;
   int          error = DBMS_ALL_OK;

	error = init_table( &table, fio->name );
		/* Tabelle initialisieren */

	if (LogFile) {
		fprintf(LogFile, "fio.c: nach init_table()\n"); 
		fflush(LogFile);
		}
	if ( ! error ) {

	reclen = table->data_length;
		/* Datensatzlaenge */
	if (LogFile) {
		fprintf(LogFile, "fio.c: reclen = %d\n", reclen); 
		fflush(LogFile);
		}

	base = (BASE_PAGE*)(table->base);

	/* Pruefe Zugriffsberechtigung: */

	if ( table->change_access == '-' ) {
		error = DBMS_ACCESS_DENIED;
		/* !!! Wird derzeit nicht benutzt !!!
		} else if ( (table->user_rights).insert < table->change_access ) {
		error = DBMS_ACCESS_DENIED;
		*/
		}

	/* eine Klammerebene offen */

	if ( ! error ) {

	/* Hole Speicher fuer die Wurzelseite der Tabelle: */
	pagebuf[0].page = (char *) memory_alloc(base->page_size);

	if ( pagebuf[0].page == NULL ) {
		error = DBMS_MEMORY_ERROR;
		}
	else {

		table->operation = FIO;
		/* Aktive Operation */

		for ( n = 0; n < 20; n++ ) { 
			pagebuf[n].offset = 0;
			}
		
		del_pages = (ULONG*)(table->base + sizeof(BASE_PAGE) );
			/* geloeschte Seiten werden in der base_page gespeichert */

		page_c = 0;
			/* 0-te page (pagebuf[0].page) ist allociert */

	/* 3 Klammerebenen offen */


	/* 
	 * Beginn der EINLESESCHLEIFE: 
	 */

	while ( ! error ) {
		if (LogFile != NULL)
			fprintf(LogFile,"\n");

		/* naechsten Datensatz einlesen (mit readf()): */
		error = (fio->readb)( fio );

		(fio->count)++;
		if ( error ) {
			if (LogFile != NULL) {
				fprintf(LogFile, "error (%d) reading record %lu\n", 
					error, fio->count );
				fflush(LogFile);
				}
			break; /* break while */
			}
		else {
			if (LogFile != NULL)
				fprintf(LogFile, "insert record %lu\n", fio->count );
			}

		del_cnt = 0;
			/* bezieht sich auf pages_to_del[] */

		table->data_buf = (UCHAR *) (fio->data);
			/* der eingelesene Datensatz steht in table->data_buf */
		table->ditype   = (fio->dataaccess) & 63;

		/* Autozaehlerindizes bearbeiten: */
		for ( i = 0; i < table->index_c; i++ ) {

			index_db = &(table->index_db[i]);
			spe = index_db->special;
			num = index_db->number;

			if ( spe == AUTOINCREMENT || spe == AUTODECREMENT ) {

				if ( base->variable ) {
					n = ((UINT *)(table->data_buf))[i+1];
					} 
				else {
					n = index_db->position;
					}

				*(long *)(table->data_buf+n) = (base->counter)[num];

				if ( spe == AUTODECREMENT ) {
					((base->counter)[num])--;
					} 
				else {
					((base->counter)[num])++;
					}
				}
			} /* Ende Autozaehlerindizes */


		if ( base->variable ) {
			/* Datensatzlaenge holen, 
			 * base->max_data_length setzen: */
			
			if ( *((UINT *)(table->data_buf)) > base->max_data_length ) {
				base->max_data_length = *((UINT *)(table->data_buf));
				}
			}

		/* 
		 * Read pages into buffers: 
		 */

		table->active = base->root_page;
			/* jeder Zugriff beginnt bei der Wurzelseite;
			 * base->root_page ist der Seitenoffset */

		if (LogFile != NULL)
			fprintf(LogFile,"read pages\n");

		for ( actp = 0; actp < base->level; actp++ ) {

			/* evtl Seite allocieren: */
			if ( actp > page_c ) { /* Seite noch nicht allociert ? */
				pagebuf[actp].page = 
					(char*) memory_alloc(base->page_size);
				if ( pagebuf[actp].page == NULL ) {
					error = DBMS_MEMORY_ERROR;
					break;
					}
				else {
					page_c++;
					}
				pagebuf[actp].offset = 0;
					/* setze offset auf 0, 
					 * damit diese Seite nicht 
					 * gleich abgespeichert wird
					 * (offset 0 kommt fuer Seiten nicht vor, 
					 * er wird fuer die base_page benuetzt). */
				}

			table->page_buf = (UCHAR *)(pagebuf[actp].page);

			if ( table->active != pagebuf[actp].offset ) {
				/* es wurde in eine andere Seite verzweigt 
				 * als beim vorigen Zugriffsweg;
				 * die alte Seite wird gerettet */

				/* war es wirklich eine bereits gueltige Seite ? */
				if ( pagebuf[actp].offset != 0 ) {
					
					/* alte Seite schreiben: */
					
					offset = table->active; /* merken */
					table->active = pagebuf[actp].offset;
					error = write_page( table );
						/* alte Seite sichern */
					pagebuf[actp].offset = 0;
					table->active = offset; /* restaurieren */
					}

				if ( ! error ) {
					/* Seite wird jetzt gelesen: */
					error = read_page( table );
					}

				if ( error ) 
					break;
				} 
			else { /* active == pagebuf[].offset */
				/* Seite muss nicht erneut gelesen werden; 
				 * sie ist bereits aus dem vorigen 
				 * Zugriff vorhanden. */
				if (LogFile != NULL)
					fprintf(LogFile,"page in buffer %ld\n", 
						table->active );
				}

			/* Seite als nach dem Gebrauch freizugeben markieren 
			 * (der vollstaendige Zugriffsweg 
			 * wird spaeter neu abgespeichert -
			 * ob die Seiten geaendert wurden oder nicht. */
			
			pages_to_del[del_cnt] = table->active;
			del_cnt++;
			
			/* die Seiten koennen deswegen noch nicht 
			 * 'richtig' (also in der base_page Liste) 
			 * freigegeben werden, weil sie sonst 
			 * bereits fuer die neuen Seiten beim 
			 * Abspeichern des Zugriffsweges 
			 * wiederverwendet wuerden. */

			if ( actp < base->level - 1 ) {
				/* noch nicht bei Datenseite angekommen;
				 * Indexseite wird durchlaufen, 
				 * EXIT Knoten gesucht;
				 * pagebuf[actp].next bekommt den Index 
				 * ins olist[] Array der Indexseite, 
				 * wo die Nachfolgerseite zu finden ist;
				 * table->active wird auf die Nachfolgerseite gesetzt 
				 * (siehe: inserth.c get_next_page()) */
				table->active = 
					get_next_page( table, &(pagebuf[actp].next) );
				}

			if ( table->active == 0 ) {
				error = DBMS_IPG_DEFECT;
				break;
				}

			} /* end for actp */

		/* Der Zugriffsweg samt Datenseite am Ende 
		 * steht nun komplett im pagebuf[] Array;
		 * es sind genau base->level Seiten, 
		 * also base->level - 1 Indexseiten. */
		offset = 0;

		if ( ! error )   {

			last = base->level - 1;
				/* index der Datenseite in pagebuf[] */

			if (LogFile!=NULL)
				fprintf(LogFile,"FIO insert data\n");

			/* Einfuegen des Datensatzes: */
			error = insert_data( table, pagebuf, 
				&last, pages_to_del, &del_cnt, &offset );
				/* inserth.c */
				/* offset enthaelt den neuen 
				 * Seitenoffset der Datenpage */

			/* bei FIO wird in insert_data keine Seite geschrieben */


		/* eine Klammerebene if ! error offen */

		/* Schreibe Seitenpuffer. Aendere dabei jeweils den      */
		/* Eintrag in der Nachfolgerliste mit dem Offset der     */
		/* nachfolgenden Seite (wie beim Einlesen vermerkt).     */

		if ( ! error ) {


			/* Zugriffsweg vorbereiten zum Zurueckschreiben:
			 * von hinten her werden neue Seiten allociert;
			 * die olist[] Seitenoffsets werden daraufhin geaendert: */

			while ( last > 0 ) {

				last--;
					/* vorletzte Seite und frueher
					 * (also keine Datenseite mehr) */

				/* offset enthaelt nun den Seitenoffset 
				 * der Seite last + 1 
				 * (beim ersten mal also den Offset der Datenseite 
				 * aus insert_data()) */

				table->page_buf = (UCHAR*)(pagebuf[last].page);
				ipage = (IPAGE *)(table->page_buf + sizeof(PTYPE));
				list  = (ULONG *)(table->page_buf + ipage->olist);
				*(list + pagebuf[last].next ) = offset;
					/* olist[] setzen */

				/* diese Seite wird allociert: */
				offset = getFreeFilePosition( table );
				pagebuf[last].offset = offset;
				}

			base->root_page = offset;
				/* 0-te Seite des Zugriffsweges ist die root-page;
				 * auch ihr offset wurde geaendert. */
			base->lognum++;
			(base->data_c)++;

			/* jetzt koennen die bislang nur 
			 * vorlaeufig freigegebenen Seiten 
			 * wirklich als frei markiert werden: */
			for( n = 0; n < del_cnt; n++ ) {
				del_pages[base->del_c] = pages_to_del[n];
				(base->del_c)++;
				}
			} /* if ! error */
		} /* eine Klammerebene if ! error */

		/* Behandlung der FIO-Spezialfaelle: */
		if ( error == DBMS_DATA_IDENTICAL ) {
			if (LogFile != NULL)
				fprintf(LogFile,"Identical %lu\n", fio->count );
			if ( fio->special & NO_BREAK_ON_IDENTICALS ) {
				error = DBMS_ALL_OK;
				}
			}


		/* 
		 * Ende der EINLESESCHLEIFE 
		 */
		} /* end while */


	/* alle Datensaetze sind eingelesen und 
	 * eingefuegt worden 
	 * (oder es ist ein Fehler aufgetreten). */

	if ( error == DBMS_FIO_END ) { 
		error = DBMS_ALL_OK;
		}

	/* ein letztes zurueckschreiben der Seiten 
	 * des Zugriffsweges: */

	for ( n = 0; n <= page_c; n++ ) {

		if ( pagebuf[n].offset != 0 && ! error ) {
			table->page_buf = (UCHAR *)(pagebuf[n].page);
			table->active   = pagebuf[n].offset;
			error = write_page( table );
			}
            	}

	/* Speicher der Seiten freigeben: */

	for ( n = 0; n <= page_c; n++ ) {
		memory_free( pagebuf[n].page );
		}

	/* es muessen noch 3 Klammerebenen geschlossen werden: */

	} 
	}

	if ( ! error ) { 
		error = write_base_page( table ); 
		}

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
int dbms_insert_f( char *tblname, char *in, int msg ) {

   FIO_PARM fp;
   int error = DBMS_ALL_OK;

   infile = fopen( in, "rb" );
   if ( infile == NULL ) {
      switch ( msg ) {
         case 0:  break;
         case 1:  printf("Fehler: kann %s nicht oeffnen\n", in ); break;
         default: printf("Error: unable to open %s\n", in );
      }
      error = DBMS_INFILE_OPEN_ERROR;

   } else {

      fp.name       = tblname;
      fp.data       = NULL;
      fp.readb      = readf;
      fp.dataaccess = 0;
      fp.count      = 0;
      fp.special    = NO_BREAK_ON_IDENTICALS;

      error = dbms_fio( &fp );

      if ( buffer != NULL ) memory_free( buffer );
      fclose( infile );
   }
   return( error );
}

/* ----------------------------------------------------------------------- */
/* Name       readf                                                        */
/* Definition int readf( FIO_PARM *fp );                                   */
/* Funktion   Einlesen der Datensaetze von der Datei.                       */
/*            Wenn der letzte Datensatz gelesen wurde, wird als Ergebnis   */
/*            DBMS_FIO_END zurueckgegeben                                   */
/* ----------------------------------------------------------------------- */
int readf( FIO_PARM *fp ) {
   uint  size = 0;

   if ( reclen > 0 ) {
      if ( bsize == 0 ) {
         buffer = (char*)memory_alloc(reclen);
         if ( buffer == NULL ) {
            return( DBMS_MEMORY_ERROR );
         }
      }
      bsize = reclen;
      if ( fread( buffer, reclen, 1, infile ) != 1 ) {
         return( DBMS_FIO_END );
      }
      fp->data = buffer;
   } else {
      if ( fread( &size, sizeof(int), 1, infile ) != 1 ) {
         return( DBMS_FIO_END );
      } else {
         if ( size > bsize ) {
            if ( buffer != NULL ) { memory_free(buffer); buffer = NULL; }
            buffer = (char*)memory_alloc(size);
            if ( buffer == NULL ) {
               return( DBMS_MEMORY_ERROR );
            }
         }
         bsize = size;
         *(uint*)buffer = size;
         if ( fread( buffer+sizeof(int), size-sizeof(int), 1, infile ) != 1 ) {
            return( DBMS_FIO_END );
         }
         fp->data = buffer;
      }
   }
   return( DBMS_ALL_OK );
}

/* ----------------------------------------------------------------------- */
/* Ende von FIO.C                                                          */
/* ----------------------------------------------------------------------- */
