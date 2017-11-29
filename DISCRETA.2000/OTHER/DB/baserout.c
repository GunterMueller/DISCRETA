/* ----------------------------------------------------------------------- */
/*                                                                         */
/* BASEROUT.C                                                              */
/*                                                                         */
/* Grundlegende Funktionen, die fuer alle wichtigen Operationen benoetigt    */
/* werden.                                                                 */
/*                                                                         */
/* Funktionen:                                                             */
/*                                                                         */
/* newPage   - Initialisieren von Seiten                                   */
/* getFreeFilePosition                                                     */
/*                                                                         */
/* makefname - Erzeugung der benoetigten Dateinamen                         */
/*                                                                         */
/* ----------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "dbms.h"             /* Benutzerstrukturen/operationen            */
#include "dbmsdef.h"          /* Definitionen                              */
#include "dbmstype.h"         /* Interne Verwaltungsstrukturen             */
#include "dbmsprot.h"         /* Funktionsprototypen (intern)              */

/* ----------------------------------------------------------------------- */
/* Name       newPage                                                      */
/* Definition void newPage( TABLE *table, UINT type )                      */
/* Prototyp   dbmsprot.h                                                   */
/* Funktion   Initialisiert eine leere Datenseite der Tabelle table.       */
/*            Hierzu muss der Seitenpuffer der Tabelle allociert sein.      */
/*            Der Parameter type kann die Werte                            */
/*               DP fuer Datenseite,                                       */
/*               IP fuer Indexseite und                                    */
/*               OP fuer Overlayseite                                      */
/*            annehemen.                                                   */
/* ----------------------------------------------------------------------- */
void newPage( TABLE *table, UINT type ) {

   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   PTYPE     *pt = (PTYPE*)(table->page_buf);
   DPAGE     *dp = (DPAGE*)(table->page_buf + sizeof(PTYPE));
   IPAGE     *ip = (IPAGE*)(table->page_buf + sizeof(PTYPE));
   OPAGE     *op = (OPAGE*)(table->page_buf + sizeof(PTYPE));

   switch ( type ) {

      case DP: {
         pt->type = DP;
         dp->dic = dp->dnc = 0;
         dp->free = sizeof(PTYPE) + sizeof(DPAGE);
         dp->data = base->page_size;
      } break;

      case IP: {
         pt->type = IP;
         ip->node_c = 0;
         ip->olist = base->page_size;
         ip->o_c   = 0;
         ip->b_c   = 0;
         ip->bpos  = base->page_size;
         ip->free  = base->page_size - sizeof(PTYPE) - sizeof(IPAGE);
      } break;

      case OP: {
         pt->type = OP;
         op->oic = 0;
         op->free = sizeof(PTYPE) + sizeof(OPAGE);
         op->data = base->page_size;
      } break;

   }

   return;
}

#if ! defined(__DBCREATE__)

/* ----------------------------------------------------------------------- */
/* Name:       getFreeFilePosition                                         */
/* Definition: ULONG getFreeFilePosition( TABLE *table );                  */
/* Prototyp:   dbmsprot.h                                                  */
/* Funktion:   Sucht eine freie Seitenadresse in der Tabellendatei zum     */
/*             Speichern einer Seite. Hierzu wird zunaechst die geloeschten  */
/*             Seiten benutzt, wenn keine mehr vorhanden ist, so wird am   */
/*             Dateiende angehaengt.                                        */
/* Ergebnis:   freie Seitenadresse                                         */
/* ----------------------------------------------------------------------- */
ULONG getFreeFilePosition( TABLE *table ) {

   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   ULONG     *del_pages = (ULONG*)(table->base + sizeof(BASE_PAGE));
   ULONG      result;

   if ( base->del_c > 0 ) {
      base->del_c--;
      result = del_pages[base->del_c];
   } else {
      result = base->file_end;
      base->file_end += base->page_size;
   }

   return( result );
}

/* ----------------------------------------------------------------------- */
/* Name:       getIndexOrder                                               */
/* Definition: void getIndexOrder( TABLE *table, UINT *I );                */
/* Prototyp:   dbmsprot.h                                                  */
/* Funktion:   Die Funktion ermittelt eine Permutation von 1, ... ,        */
/*             table->index_c - 1 fuer die Auswahl der Reihenfolge der      */
/*             Indexkomponenten beim Erzeugen neuer Knoten auf Daten- bzw. */
/*             Indexseiten. I muss mindestens mit der Laenge table->index_c */
/*             * sizeof(UINT) allociert sein und enthaelt nach Aufruf die   */
/*             ermittelte Reihenfolge.                                     */
/* ----------------------------------------------------------------------- */
void getIndexOrder( TABLE *table, UINT *I ) {

   UINT   len = table->index_c;
   UINT   x;
   UINT   h;
   UINT   n;
   time_t timer;

   srand( (unsigned)time( &timer ) );

   for ( n=0;n <table->index_c; n++ ) {
      I[n] = n;
   }

   for ( n = 0; n < table->index_c; n++ ) {
      x = (rand() % len) + n;
      h    = I[n]; I[n] = I[x]; I[x] = h;
      len--;
   }

   return;
}

/* ----------------------------------------------------------------------- */
/* Name:       compSize                                                    */
/* Definition: UINT compSize( TABLE *table, UINT comp, UCHAR *data );      */
/* Prototyp:   dbmsprot.h                                                  */
/* Funktion:   Ermittelt die Groesse der Komponente comp des Datensatzes     */
/*             data in Bytes. Kann auch fuer die Restkomponente verwendet   */
/*             werden.                                                     */
/* Ergebnis:   Groesse der Komponente in Bytes.                            */
/* ----------------------------------------------------------------------- */
/* Berechnet die Groesse in bytes der Komponente 'comp';
 * Im Falle base->variable (d.h. variable Datensatzlaenge), 
 * wird die Groesse aus der Differenz des i und i+1 ten 
 * Schluessel-Zeigers berechnet (an Position i+1 und i+2).
 * zur Erinnerung: an Position 0 ist die gesamte L"ange des 
 * Datensatzes eingetragen, dann folgen die Pointer fuer die 
 * Schluesseloffsets, dann die Schluessel selber 
 * und schliesslich der Datensatz selbst:
 * 
 * UINT pos[0] : gesamte L"ange des Datensatzes 
 *        (samt Schluessel und Pointer)
 * UINT pos[1...table->index_c] : die Offsets der Schluessel
 * UINT pos[table->index_c + 1] : der Beginn des Datensatzes 
 *        (sogenannter INFO-Teil)
 * hier kommen die Schluessel
 * hier kommt der INFO-Teil
 * 
 * Im Falle comp = table->index_c wird die L"ange des INFO-Teiles 
 *       zur"uckgegeben.
 *
 * Im Falle fester Datensatzgroesse werden die Informationen 
 * aus table->index_db[] und table->data_length berechnet. */

UINT compSize( TABLE *table, UINT comp, UCHAR *data ) {

   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   UINT      *pos = (UINT*)data;
   UINT       result = 0;

   if ( base->variable ) {
      if( comp == table->index_c )
		result = pos[0] - pos[comp+1];
      else
		result = pos[comp+2] - pos[comp+1];
   } else {
      if( comp == table->index_c ) 
         result = table->data_length - 
		table->index_db[comp].position;
      else 
         result = table->index_db[comp+1].position - 
		table->index_db[comp].position;
   }

   return( result );
}

/* ----------------------------------------------------------------------- */
/* Name:       copyData                                                    */
/* Definition: void copyData( TABLE *table, UCHAR *data, UINT length,      */
/*                            DINFO *dinfo );                              */
/* Prototyp:   dbmsprot.h                                                  */
/* Funktion:   Kopiert einen Datensatz data samt seines DINFO dinfo auf    */
/*             die aktuelle Seite im Seitenpuffer der Tabelle. length ist  */
/*             die Groesse von data in Bytes.                                */
/* ----------------------------------------------------------------------- */
void copyData( TABLE *table, UCHAR *data, DINFO *dinfo ) {

   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   DPAGE     *p  = (DPAGE*)(table->page_buf + sizeof(PTYPE));
   DNODE     *dn = (DNODE*)(table->page_buf + sizeof(PTYPE) + sizeof(DPAGE));
   DINFO     *di = (DINFO*)(table->page_buf + p->data);
   UINT       I[INDEXC];     /* Feld fuer Indexreihenfolge                 */
   UINT       n,i;           /* Laufvariablen                              */
   UCHAR     *v1, *v2;       /* Zeiger auf Vergleichswerte                 */
   UCHAR     *source, *dest; /* Zeiger auf Speicherbereiche zum Umkopieren */
   SINT       res;           /* Ergebnis der Vergleichsfunktionsaufrufe    */
   UINT       x,z;           /* Hilfsvariablen                             */
   UINT       x_next;        /* _SMALLER / _GREATER je nach Nachfolger     */
                             /* beim Suchen                                */
   UINT       length;        /* Groesse des Datensatzes in Bytes           */

   if ( base->variable ) {
      length = *(UINT*)data;
   } else {
      if( dinfo->type & DIOVERLAY )
         length = table->index_db[ table->index_c ].position;
      else
         length = table->data_length;
   }

   if ( p->dic == 0 ) {

      /* ----------------------------------------------------------------- */
      /* Es befindet noch kein Datensatz auf der Seite. Lege DATA_INFO fuer */
      /* neuen Datensatz an.  di_start wird auf Seitengroesse - length -     */
      /* Groesse eines DATA_INFOs gesetzt und auf diese Position, wird auch  */
      /* das erste DATA_INFO gelegt. nach diesem wird der Datensatz auf    */
      /* die Seite kopiert.                                                */
      /* ----------------------------------------------------------------- */

      p->dic  = 1;
      p->dnc  = 0;
      p->free = sizeof(PTYPE) + sizeof(DPAGE);
      p->data = base->page_size - length - sizeof(DINFO);

      di = (DINFO*)(table->page_buf + p->data);

      di[0].position = p->data + sizeof(DINFO);
      di[0].type     = dinfo->type;
      di[0].lognum   = dinfo->lognum;

      memory_copy( table->page_buf + di[0].position, data, length );

   } else {

      if ( p->dnc == 0 ) {

         /* -------------------------------------------------------------- */
         /* Es befindet sich nur ein Datensatz auf der Seite.   Dazu wird  */
         /* eine Reihenfolge der Indexkomponenten ermittelt und die beiden */
         /* Datensaetze anhand dieser Reihenfolge bzgl. der jeweiligen Kom- */
         /* ponenten verglichen um einen Trennwert zu ermittlen.           */
         /* -------------------------------------------------------------- */

         getIndexOrder( table, I );               /* Reihenfolge ermitteln */

         x = 1;

         for ( n=0;n<table->index_c;n++ ) { /* vergleiche komponentenweise */

            i  = I[n];
            v1 = proj( table, i, data);
            v2 = proj( table, i, table->page_buf + di[0].position );   

            res = (table->index_db[i].cmpfunc)( v1, v2, table->index_db[i].length );
            if ( res != 0 ) break;
         }

         dn[0].ds    = 0;                          /* Lege neuen Knoten an */
         dn[0].comp  = i;
         dn[0].total = base->page_size - p->data - sizeof(DINFO) + length;
         dn[0].mode  = SMALLER_DS | GREATER_DS;

                                      /* Vergleichswert ist Komponente     */
         if ( res == 1 ) {            /* des schon vorhandenen Datensatzes */
            dn[0].smaller = 0; 
            dn[0].greater = 1;       /* falls diese der kleinere Wert ist, */
            dn[0].mode |= EQ_SMALLER;  /* setze equal_smaller, sonst nicht */
         } else {
            dn[0].smaller = 1; 
            dn[0].greater = 0;
         }

         source = table->page_buf + p->data;         /* schaffe Platz fuer */
         dest   = source - sizeof(DINFO) - length;   /* neues DINFO       */
         memory_copy( dest, source, sizeof(DINFO) );

         p->data -= sizeof(DINFO) + length;
         p->free += sizeof(DNODE);

         di = (DINFO*)(table->page_buf + p->data);       /* lege DINFO an */
         di[1].position = p->data + (p->dic+1)*sizeof(DINFO);
         di[1].type     = dinfo->type;
         di[1].lognum   = dinfo->lognum;

         memory_copy( table->page_buf + di[1].position, data, length );

         p->dic = 2;
         p->dnc = 1;

      } else {

         /* -------------------------------------------------------------- */
         /* Es befinden sich mehrere Datensaetze auf der Seite.            */
         /* Suche Blattknoten bei dem der neue Datensatz eingefuegt werden  */
         /* muss. Rest wie bei nur einem Datensatz.                        */
         /* -------------------------------------------------------------- */
         x = 0;

         while ( 1 ) {
            dn[ x ].total += length;          /* erhoehe total um die Groesse */
            i  = dn[x].comp;                  /* des neuen DS              */
            v1 = proj( table, i, data);
            v2 = proj( table, i, table->page_buf + di[dn[x].ds].position );   

            res = table->index_db[i].cmpfunc(v1,v2,table->index_db[i].length);

            switch ( res ) {
               case  1:    x_next = _GREATER; break;
               case  0:  if( dn[x].mode & EQ_SMALLER ) x_next = _SMALLER;
                     else x_next = _GREATER;
                     break;
               case -1:  x_next = _SMALLER; break;
            }

            if ( x_next == _SMALLER ) {
               if( dn[x].mode & SMALLER_DS ) { z = dn[x].smaller; break; }
               else   x = dn[x].smaller;
            } else {
               if( dn[x].mode & GREATER_DS ) { z = dn[x].greater; break; }
               else x = dn[x].greater;
            }
         }

         getIndexOrder( table, I );                    /* Indexreihenfolge */

         for ( n = 0; n < table->index_c; n++ ) {       /* suche Trennwert */
            i = I[n];
            v1 = proj( table, i, data);
            v2 = proj( table, i, table->page_buf + di[z].position );   
            res=(table->index_db[i].cmpfunc)( v1, v2, table->index_db[i].length );
            if ( res != 0 ) break;
         }

         if ( x_next == _SMALLER ) {             /* fuege neuen DNODE ein   */
            dn[x].mode &= NOT_SMALLER_DS;        /* als kleinerer NF von x */
            dn[x].smaller = p->dnc;
         } else {
            dn[x].mode &= NOT_GREATER_DS;        /* bzw. groesserer NF       */
            dn[x].greater = p->dnc;
         }

         if( ! base->variable )
            dn[p->dnc].total = table->data_length + length;
         else   
            dn[p->dnc].total = *(UINT*)(table->page_buf + di[z].position) + length;

         dn[p->dnc].ds    = z;                     /* erzeuge neuen Knoten */
         dn[p->dnc].comp  = i;
         dn[p->dnc].mode = SMALLER_DS | GREATER_DS;

         if ( res == 1 ) { /* setze equal-Flag nach Wert des vorhandenen DS */
            dn[p->dnc].smaller = z; dn[p->dnc].greater = p->dic;
            dn[p->dnc].mode |= EQ_SMALLER;
         } else {
            dn[p->dnc].smaller = p->dic; dn[p->dnc].greater = z;
         }

         source = table->page_buf + p->data;     /* schaffe Platz fuer DINFO */
         dest   = source - sizeof(DINFO) - length;
         memory_copy( dest, source, p->dic*sizeof(DINFO) );

         p->data -= sizeof(DINFO) + length;
         p->free += sizeof(DNODE);

         di = (DINFO*)(table->page_buf + p->data);         /* erzeuge DINFO */

         di[p->dic].position = p->data + (p->dic+1)*sizeof(DINFO);
         di[p->dic].type     = dinfo->type;
         di[p->dic].lognum   = dinfo->lognum;

         memory_copy(table->page_buf+di[p->dic].position,data,length);

         (p->dic)++; (p->dnc)++;
      }
   }

   return;
}

#endif

/* ----------------------------------------------------------------------- */
/* Name       makefname                                                    */
/* Definition void makefname( char*fname, char *path, char *filename,      */
/*                            char *extension )                            */
/* Prototyp   dbmsprot.h                                                   */
/* Funktion   Erzeugt aus dem angegebenen Pfad, dem Dateinamen und der     */
/*            Extension einen vollstaendigen Dateinamen und schreibt diesen */
/*            nach fname. fname muss gross genug sein, um den kompletten     */
/*            String aufnehmen zu koennen!                                  */
/*            Wird fuer Pfad oder Extenension NULL angegeben, so wird       */
/*            der Name ohne den jeweiligen Teilstring gebildet.            */
/*            Die Verzeichnistrennungssymbole (DIRSEP) werden automatisch  */
/*            ergaenzt.                                                     */
/* ----------------------------------------------------------------------- */
void makefname( char *fname, char *path, char *filename, char *extension ) {

   int l;

   if ( filename[0] == '%' ) {
      strcpy( fname, filename+1 );
   } else {
      if ( path != NULL ) {
         strcpy( fname, path );
         l = strlen(fname);
         if ( fname[l-1] != DIRSEP ) {
            fname[l]   = DIRSEP;
            fname[l+1] = '\0';
         }
         strcat( fname, filename );
      } else {
         strcpy( fname, filename );
      }
   }
   if ( extension != NULL ) {
	l = strlen(filename);
	if ( l > strlen(extension) ) {
		if ( strcmp( filename + l - strlen(extension), extension ) != 0 ) {
			strcat( fname, extension );
			}
		}
	else {
		strcat( fname, extension );
		}
	}
   return;
}

#undef DEBUG

void Log(char *s)
{
#ifdef DEBUG
	if (LogFile) {
		fprintf(LogFile, "%s", s);
		fflush(LogFile);
		}
	printf("%s", s);
	fflush(stdout);
#endif
}

void *Memmove(void *s1, const void *s2, int len)
{
#if 0
#ifdef DEBUG
	Log("Memmove(");
#endif
#endif
	memmove(s1, s2, len);
#if 0
#ifdef DEBUG
	Log(")\n");
#endif
#endif
	return s1;
}

void Memory_zero(void *p, int l)
{
	char *pp = (char *) p;
	int i;

	for (i = 0; i < l; i++) 
		pp[i] = (char) 0;
}


/* ----------------------------------------------------------------------- */
/* Ende von BASEROUT.C                                                     */
/* ----------------------------------------------------------------------- */
