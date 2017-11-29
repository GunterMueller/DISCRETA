/* ----------------------------------------------------------------------- */
/* INSERTO.C                                                               */
/*                                                                         */
/* Diese Datei enthaelt die Funktionen zum Erzeugen neuer Overlayseiten     */
/* und Einfuegen auf diesen. Derzeit gibt eine Beschraenkung fuer die Hoehe    */
/* des Overlaybaumes auf die Hoehe 2. Als Seite 0 wird die Wurzelseite      */
/* dieses Baumes bezeichnet, Seite 1 ist eine Seite der naechsttieferen    */
/* Ebene.                                                                  */
/*                                                                         */
/* Funktionen                                                              */
/* create_new_overlay                                                      */
/* insert_overlay                                                          */
/* ----------------------------------------------------------------------- */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dbms.h"            /* Benutzerstrukturen/operationen             */
#include "dbmsdef.h"         /* Definitionen                               */
#include "dbmstype.h"        /* Interne Verwaltungsstrukturen              */
#include "dbmsprot.h"        /* Prototypen interner Funktionen             */

/* ----------------------------------------------------------------------- */
/* Name       create_new_overlay                                           */
/* Definition UINT create_new_overlay( TABLE *table, UCHAR *overlay,       */
/*                                                            UINT num );  */
/* Prototyp   dbmsprot.h                                                   */
/* Funktion   Die Funktion erzeugt eine neue, leere Overlayseite (in       */
/*            overlay). Dann werden die Restkomponenten des Datensatzes    */
/*            auf der Datenseite und des einzufuegenden Datensatzes         */
/*            ermittelt und auf der Overlayseite abgelegt. Der ursprueng-   */
/*            liche Datensatz der Datenseite wird in ein TID umgewandelt,  */
/*            d.h. die lognum wird als Verweis auf eine Overlayseite       */
/*            benutzt und die Daten um die Restkomponente verkuerzt. Dazu   */
/*            wird der Datenblock von di_start bis zum Ende der verkuerzten */
/*            TID-Daten um die Laenge der Restkomponente nach 'hinten'      */
/*            kopiert. Anschliessend werden die Datenpositionen der        */
/*            DATA_INFO's nach dem TID berichtigt.                         */
/* Parameter  table   betroffene Tabelle                                   */
/*            overlay Zeiger auf Puffer fuer die neue Overlayseite         */
/*                    (muss allociert sein)                                */
/*            num     Feldindex des betroffenen Datensatzes auf der        */
/*                    Datenseite.                                          */
/* Ergebnis   Wenn die Operation erfolgreich ausgefuehrt wurde: DB_ALL_OK, */
/*            sonst Fehlermeldung:                                         */
/*            DB_WRITE_ERROR                                               */
/*            Ursache: die neu konstruierte Overlay-Seite konnte nicht     */
/*            geschrieben werden.                                          */
/*            Woher:   write_page                                          */
/*            DB_DATA_IDENTICAL                                            */
/*            Ursache: Der einzufuegende und der vorhandene Datensatz       */
/*            sind identisch.                                              */
/*            Woher:   create_new_overlay                                  */
/* Benutzte Funktionen                                                     */
/*            comp_size              (INSERTH.C)                           */
/*            proj - Makro                                                 */
/*            put_R_comp             (INSERT2.C)                           */
/*            get_free_file_position (BASEROUT.C)                          */
/*            write_page             (FILE_OP.C)                           */
/*            memory_copy            (MEMORY.C)                            */
/* ----------------------------------------------------------------------- */
int create_new_overlay( TABLE *table, UINT num ) {

   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   UCHAR     *overlay;
   UCHAR     *data_p = table->page_buf;                    
   ULONG      old_a  = table->active;                      
   DPAGE     *p      = (DPAGE*)(table->page_buf+sizeof(PTYPE));
   DINFO     *di     = (DINFO*)(data_p + p->data);
   DINFO      di_new;                         
   UINT       R_size_new;                     
   UINT       R_size_dp;                      
   UINT       new_size;                       
   UCHAR     *source;                         
   UINT       n;
   UINT       R = table->index_c;         
   int        error = DBMS_ALL_OK;

   R_size_new = compSize( table, R, table->data_buf );
   R_size_dp  = compSize( table, R, data_p + di[num].position );
   source = proj( table, R, table->page_buf + di[num].position );

   if ( base->variable ) {

      if ( R_size_new == R_size_dp ) {

         if ( memcmp( source, proj(table,R,table->data_buf), R_size_dp ) == 0 ) {
            error = DBMS_DATA_IDENTICAL;
         }
      }

   } else {

      if ( memcmp( source, proj(table,R,table->data_buf), R_size_new ) == 0 ) {
         error = DBMS_DATA_IDENTICAL;
      }
   }

   if ( ! error ) {

      /* ----------------------------------------------------------------- */
      /* overlay wird als leere Overlayseite initialisiert, ein DATA_INFO  */
      /* info_new fuer die Restkomponente des einzufuegenden Datensatzes   */
      /* angelegt. Dann werden die Restkomponente des neuen und des Daten- */
      /* satzes num auf der Datenseite auf die Overlayseite kopiert und    */
      /* diese geschrieben.                                                */
      /* ----------------------------------------------------------------- */

      overlay = (UCHAR*)memory_alloc( base->page_size);

      if( overlay == NULL ) {

         error = DBMS_MEMORY_ERROR;

      } else {

         table->page_buf = overlay;       /* initialisiere Seite */
         newPage( table, OP );

         di_new.type    = table->ditype;    /* lege DATA_INFO an */
         di_new.lognum  = base->lognum;

         put_R_comp( overlay, proj( table, R, table->data_buf ), R_size_new, &di_new );
         put_R_comp( overlay, proj( table, R, data_p + di[num].position ), R_size_dp , di+num );

         table->active = getFreeFilePosition( table );  /* schreibe Seite */

         error = write_page( table );

         if ( ! error ) {

            /* ----------------------------------------------------------- */
            /* Wandle DATA_INFO auf der Datenseite in TID um               */
            /* Nun wird der urspruengliche Datensatz auf der Datenseite in  */
            /* einen Verweis auf eine Overlayseite (TID) umgewandelt. Dazu */
            /* wird lognum als offset fuer die Overlayseite verwendet. Das */
            /* tid-Flag wird gesetzt und die Flags fuer die Zugriffs auf 0 */
            /* gesetzt um den Zugriff auf die Daten zu gewaehrleisten.     */
            /* Dann wird die Restkomonente des Datensatzes geloescht. Die   */
            /* Berechnungen hierzu (new_size, size, source, dest) sind in  */
            /* der Dokumentation zu finden.                                */
            /* ----------------------------------------------------------- */

            di[num].lognum = table->active;
            di[num].type   = DIOVERLAY;

            if ( base->variable ) { /* new_size = neue Groesse des DS auf der Datenseite */

               new_size = *(UINT*)(data_p + di[num].position) - R_size_dp;
               *(UINT*)(data_p + di[num].position) = new_size;

            } else {

               new_size = table->data_length - R_size_dp;

            }

            source = data_p + p->data;             
            memory_copy( source + R_size_dp, source, di[num].position + new_size - p->data );

            p->data += R_size_dp;                             /* Berichtige Positionen */
            di = (DINFO*)(data_p + p->data);
            for ( n = num; n < p->dic; n++ ) di[n].position += R_size_dp;
         }
         memory_free( overlay );
      }
   }

   table->page_buf = data_p;
   table->active   = old_a;

   return( error );
}

/* ----------------------------------------------------------------------- */
/* Name       insert_overlay                                               */
/* Definition UINT insert_overlay( TABLE *table, UINT num,                 */
/*                                 ULONG *pages_to_del, UINT *del_cnt );   */
/* Prototyp   insert.h                                                     */
/* Funktion   Die Funktion fuegt Datensaetze auf existierenden Overlay-    */
/*            seiten ein. Hierzu wird zuerst geprueft, ob ein Datensatz     */
/*            gleicher Restkomponente schon existiert. Dazu wird der       */
/*            Overlay-Baum (max. Hoehe 2) durchlaufen. Ist eine solche      */
/*            vorhanden, so muss nicht eingefuegt werden. Wenn nicht, so     */
/*            wird versucht, auf Seite 0 einzufuegen, geht dies nicht, wird */
/*            eine Seite der Stufe 1 gesucht, auf die der Datensatz passt. */
/*            Kann auch eine solche nicht gefunden werden, versuche, ob    */
/*            noch ein DATA_INFO als TID auf Seite 0 passt. Geht dies      */
/*            nicht, dann wandle, falls ex. ein DATA_INFO in TID um.       */
/*            Sind bereits alle DATA_INFO's TID, so kann der neue Daten-   */
/*            satz nicht mehr eingefuegt werden. Normalerweise sollte die   */
/*            Konzeption fuer alle Overlaydaten ausreichen. Tritt dieser    */
/*            Fehler auf, ist vermutlich die Datenmodellierung nicht ganz  */
/*            optimal. Ausweichmoeglichkeit ist die Seitengroesse der Tabelle */
/*            zu erhoehen.                                                  */
/* Parameter  table        betroffene Tabelle                              */
/*            num          Feldindex des betroffenen TID auf der Datenseite*/
/*            pages_to_del Feld fuer die Seiten, die nach erfolgreicher     */
/*            del_cnt      Operation geloescht werden und deren Anzahl      */
/* Ergebnis   Wenn die Operation erfolgreich ausgefuehrt wurde: DB_ALL_OK  */
/*            sonst Fehlermeldung:                                         */
/*            DB_MEMORY_ERROR                                              */
/*            Ursache: Kein Speicherplatz fuer Overlay-Puffer              */
/*            Woher:   insert_overlay                                      */
/*            DB_READ_ERROR                                                */
/*            Ursache: Seite 0 oder 1 kann nicht gelesen werden.           */
/*            Woher:   read_page                                           */
/*            DB_WRITE_ERROR                                               */
/*            Ursache: die neu konstruierte Overlay-Seite konnte nicht     */
/*            geschrieben werden.                                          */
/*            Woher:   write_page                                          */
/*            DB_DATA_IDENTICAL                                            */
/*            Ursache: der Datensatz existiert schon. (kein Fehler)        */
/*            Woher:   insert_overlay                                      */
/*            DB_TOO_MUCH_OVERLAY                                          */
/*            Ursache: kein Platz, fuer neue Overlay-Daten                 */
/*            Woher:   insert_overlay                                      */
/* Benutzte Funktionen                                                     */
/*            R_size      (INSERTH.C)                                      */
/*            R_comp      (INSERTH.C)                                      */
/*            put_R_comp  (INSERT2.C)                                      */
/*            write_page  (FILE_OP.C)                                      */
/*            read_page   (FILE_OP.C)                                      */
/*            memory_copy (MEMORY.C)                                       */
/*            memcmp      (Standardfunktion)                               */
/* ----------------------------------------------------------------------- */
int insert_overlay( TABLE *table, UINT num, ULONG *pages_to_del, UINT *del_cnt ) {

   BASE_PAGE *base = (BASE_PAGE*)(table->base);
   UCHAR     *old_p = table->page_buf;
   ULONG      old_a = table->active;
   DPAGE     *dp = (DPAGE*)(old_p + sizeof(PTYPE));
   DINFO     *di = (DINFO*)(old_p + dp->data);
   UCHAR     *overlay1;
   UCHAR     *overlay2;
   UCHAR     *R_comp_new;
   UINT       R_size_new;
   ULONG      op0 = di[num].lognum;
   OINFO     *oi0, *oi1;
   OPAGE     *p0, *p1;
   UINT       n0, n1;              /* Laufvariablen zum Durchlaufen der S. */
   UCHAR     *data;                /* Zeiger auf Daten eines DATA_INFO's   */
   UINT       size;                /* Groesse der Restk.                   */
   DINFO      di_new;              /* DATA_INFO fuer einzufuegenden DS     */
   UINT       op1_free = 0;        /* Flags fuer Seite mit Platz gefunden  */
   UINT       di_free;             /* Feldindex des zug. DATA_INFO's       */
   UINT       all_tid = 1;         /* Flag fuer alle DATA_INFO's auf 0 TID */
   UCHAR     *source;              /* Anfangsadresse des umzukopierenden   */
   UCHAR     *dest;                /* Speicherblocks und Zieladresse sowie */
   UINT       blk_size;            /* Groesse des Blocks in Bytes          */
   int        error = DBMS_ALL_OK;

   if (op0 == table->active ) {
      printf("fatal error: bad overlay root position (%ld)\n",op0);
      return DBMS_READ_ERROR;
   }

   /* ---< Initialisierung >---------------------------------------------- */
   /* Bestimme die Restkomponente des einzufuegenden Datensatzes und deren */
   /* Groesse. Lege DATA_INFO fuer den neuen Datensatz an und belege Speicher */
   /* fuer den overlay-Puffer. Dieser wird mit doppelter Seitengroesse belegt,*/
   /* da der Overlay-Baum max. Hoehe 2 haben kann. Setze dann PAGE2-Zeiger  */
   /* und DATA_INFO-Zeiger auf die beiden Seiten.                          */
   /* -------------------------------------------------------------------- */

   R_comp_new = proj( table, table->index_c, table->data_buf );      /* Restk. ermitteln */
   R_size_new = compSize( table, table->index_c, table->data_buf );
   size       = R_size_new + sizeof(OINFO);

   di_new.type    = table->ditype;    /* DATA_INFO anlegen */
   di_new.lognum  = base->lognum;

   /* -------------------------------------------------------------------- */
   /* Lege Pufferspeicher fuer die beiden Seiten an. Falls dies nicht       */
   /* moeglich ist, breche mit Speicherfehler ab                            */
   /* -------------------------------------------------------------------- */
   overlay1 = (UCHAR*)memory_alloc( base->page_size );
   overlay2 = (UCHAR*)memory_alloc( base->page_size );
   if ( overlay1 == NULL || overlay2 == NULL ) {
      if ( overlay1 != NULL ) memory_free( overlay1 );
      if ( overlay2 != NULL ) memory_free( overlay2 );
      if (LogFile!=NULL) fprintf(LogFile,"no memory for overlay page buffers\n");
      error = DBMS_MEMORY_ERROR;
      return error;
   }

   /* -------------------------------------------------------------------- */
   /* Setze Zeiger auf die OPAGE-Struktur und das OINFO-Feld               */
   /* p0  Overlay-Seitenstruktur Seite 0                                   */
   /* p1  Overlay-Seitenstruktur Seite 1                                   */
   /* oi0 Feld mit den OINFO's der Seite 0                                 */
   /* oi1 Feld mit den OINFO's der Seite 1                                 */
   /* -------------------------------------------------------------------- */
   p0  = (OPAGE*)(overlay1 + sizeof(PTYPE) );
   p1  = (OPAGE*)(overlay2 + sizeof(PTYPE) );
   oi0 = (OINFO*)(overlay1 + sizeof(PTYPE) + sizeof(OPAGE));
   oi1 = (OINFO*)(overlay2 + sizeof(PTYPE) + sizeof(OPAGE));

   /* ---< Lese Overlayseite 0 in den Puffer ein >--------------------- */
   table->active = op0;
   table->page_buf = overlay1;
   error = read_page( table );

   if ( error ) {
      if (LogFile!=NULL) fprintf(LogFile,"error reading overlay root page (%ld)\n", table->active );
   } else {

      if (LogFile!=NULL) fprintf(LogFile,"reading overlay root page (%ld) ok\n", table->active );

      /* ----------------------------------------------------------------- */
      /* Durchlaufe alle OINFOs der Seite 0, und suche identische          */
      /* Restkomponente. Wird dabei ein tid gefunden,  so wird die         */
      /* zugehoerige Seite geladen und auf dieser weitergesucht. Ausserdem   */
      /* wird geprueft, ob auf der Seite Platz zum Einfuegen des neuen     */
      /* Datensatz waere. Ist dies der Fall, so wird o1_free auf 1         */
      /* gesetzt und der Feldindex des DATA_INFO's zu dieser Seite in      */
      /* di_free vermerkt.                                                 */
      /* ----------------------------------------------------------------- */

      if (LogFile!=NULL) fprintf(LogFile,"searching for identical record\n");

      for ( n0 = 0; n0 < p0->oic; n0++ ) {
         error = DBMS_DATA_IDENTICAL;

         if ( oi0[n0].type & DIOVERLAY ) {

            /* ----------------------------------------------------------- */
            /* Weitere Seite unterhalb des aktuellen OINFO's               */
            /* Lade diese Seite und suche dort weiter                      */
            /* Bei Lesefehler erfolgt sofortiger Abbruch                   */
            /* ----------------------------------------------------------- */
            table->page_buf = overlay2;
            table->active   = oi0[n0].lognum;
            error = read_page( table );
            if ( error ) {
               if (LogFile!=NULL) fprintf(LogFile,"error reading page on level 1 (%ld)\n", table->active );
               break;
            }

            if (LogFile!=NULL) fprintf(LogFile,"searching on page level 1 (%ld)\n", table->active );

            /* ----------------------------------------------------------- */
            /* Wenn genug Platz auf diese Seite zur Aufnahme des neuen     */
            /* Satzes ist, so setze Flag fuer Einfuegen auf Stufe 1 moeglich  */
            /* (op1_free) und vermerke die Nr. des OINFO's auf Seite 0 an  */
            /* dem die Seite haengt (di_free)                               */
            /* ----------------------------------------------------------- */
            if ( (!op1_free) && (p1->data - p1->free > size) ) {
               op1_free = 1;
               di_free = n0;
            }

            /* ----------------------------------------------------------- */
            /* Durchsuche die Seite der Stufe 1                            */
            /* Wenn identische R-Komp. gefunden wurde => Abbruch           */
            /* ----------------------------------------------------------- */
            for ( n1 = 0; n1 < p1->oic; n1++ ) {
               error = DBMS_DATA_IDENTICAL;
               data = table->page_buf + oi1[n1].position;
               if ( base->variable ) {
                  if ( oi1[n1].size == R_size_new ) {
                     if ( memcmp( data, R_comp_new, R_size_new ) != 0 ) {
                        error = DBMS_ALL_OK;
                     }
                  } else {
                     error = DBMS_ALL_OK;
                  }
               } else {
                  if ( memcmp( data, R_comp_new, R_size_new ) != 0 ) {
                     error = DBMS_ALL_OK;
                  }
               }
               if ( error == DBMS_DATA_IDENTICAL ) {
                  if (LogFile!=NULL) fprintf(LogFile,"identical found (%u)\n", n1 );
                  break;
               }
            }

            /* ----------------------------------------------------------- */
            /* Kehre zur Seite 0 zurueck                                    */
            /* ----------------------------------------------------------- */
            table->active   = op0;
            table->page_buf = overlay1;
            if ( error ) break;

         } else {

            /* ----------------------------------------------------------- */
            /* Vergleiche ob R-Komp. zum OINFO identisch zum der des       */
            /* neuen Datensatzes ist, falls ja => Abbruch                  */
            /* ----------------------------------------------------------- */
            all_tid = 0;
            data = table->page_buf + oi0[n0].position;
            if ( base->variable ) {
               if ( oi0[n0].size == R_size_new ) {
                  if ( memcmp( data, R_comp_new, R_size_new ) != 0 ) {
                     error = DBMS_ALL_OK;
                  }
               } else {
                  error = DBMS_ALL_OK;
               }
            } else {
               if ( memcmp( data, R_comp_new, R_size_new ) != 0 ) {
                 error = DBMS_ALL_OK;
               }
            }
            if ( error == DBMS_DATA_IDENTICAL ) {
               if (LogFile!=NULL) fprintf(LogFile,"identical found (%u)\n", n0 );
               break;
            }
         }
      }
   }

   /* ---< Restkomponente einfuegen ? >------------------------------------ */
   /* Wenn alle Operationen bisher erfolgreich abgelaufen sind und keine   */
   /* Restkomponente gefunden wurde, die mit der der einzufuegenden         */
   /* Datensatzes identisch ist, so muss neu eingefuegt werden. Dazu wird    */
   /* zuerst versucht, die neuen Daten auf der Seite 0 des Overlaybaumes   */
   /* unterzubrigen. Geht dies nicht, so wird nachgesehen, ob beim         */
   /* Suchen eine Seite der Ebene 1 ausgemacht wurde, auf die eingefuegt   */
   /* werden kann. Ist dies nicht der Fall, so muss eine neue Overlay-      */
   /* seite der Ebene 1 erzeugt werden. Hierfuer wird geprueft, ob ein     */
   /* DATA_INFO als Tid auf die Seite passt, wenn nicht wird, falls noch    */
   /* moeglich, ein Datensatz in ein TID umgewandelt. Geht alles nicht,     */
   /* so ist der Overlay-Baum voll => Fehlermeldung. siehe oben.           */
   /* -------------------------------------------------------------------- */

   if ( ! error ) {

      if (LogFile!=NULL) fprintf(LogFile,"no identical record found\n");

      if ( p0->data - p0->free > size ) {

         /* -------------------------------------------------------------- */
         /* Wenn auf der Overlayseite 0 ausreichend Platz vorhanden ist,   */
         /* so wird auf dieser eingefuegt. Dann wird die Seite geschrieben  */
         /* -------------------------------------------------------------- */
         if (LogFile!=NULL) {
            fprintf(LogFile,"store info part on root page\n");
            fprintf(LogFile,"data=%u ", p0->data );
            fprintf(LogFile,"free=%u ", p0->free );
            fprintf(LogFile,"size=%u\n", size );
         }
         put_R_comp( table->page_buf, R_comp_new, R_size_new, &di_new );

      } else {

         if ( op1_free ) {

            /* ----------------------------------------------------------- */
            /* Auf einer Seite der Stufe 1 ist Platz zum Einfuegen          */
            /* Lade die Seite in den Overlay-Puffer auf die 2.Position     */
            /* und fuege die Restkomponente des neuen Datensatz ein. An-    */
            /* schliessend schreibe die Seite, aendere Seitenverweis im      */
            /* zugehoerigem DATA_INFO der Overseite 0 auf den neuen         */
            /* Seitenoffset. Schreibe dann die Seite 0.                    */
            /* ----------------------------------------------------------- */
            table->page_buf = overlay2;
            table->active = oi0[di_free].lognum;
            if (LogFile!=NULL) fprintf(LogFile,"store info part on page level 1 (%ld)\n", table->active );
            error = read_page( table );

            if ( error ) {
               if (LogFile!=NULL) fprintf(LogFile,"error reading page\n");
            } else {
               put_R_comp( table->page_buf, R_comp_new, R_size_new, &di_new );
               table->active = getFreeFilePosition( table );
               error = write_page( table );
               if ( error ) {
                  if (LogFile!=NULL) fprintf(LogFile,"error writing page (%ld)\n", table->active );
               } else {
                  if (LogFile!=NULL) fprintf(LogFile,"writing page (%ld) ok\n", table->active );
                  pages_to_del[*del_cnt] = oi0[di_free].lognum; (*del_cnt)++;
                  oi0[di_free].lognum = table->active;
               }
            }

         } else {

            /* ----------------------------------------------------------- */
            /* Wenn all_tid 1 ist, heisst dies, dass alle Datensaetze      */
            /* auf der Seite schon tid's sind. Da keine Seite mit Platz    */
            /* zum Einfuegen gefunden wurde, ist somit das Einfuegen des    */
            /* neuen Datensatz nicht mehr moeglich. Es wird mit Fehler-    */
            /* meldung abgebrochen. Abhilfe sollte dann bessere Daten-     */
            /* modellierung oder groessere Seitengroesse bieten.               */
            /* ----------------------------------------------------------- */
            if ( all_tid ) {
               if (LogFile!=NULL) fprintf(LogFile,"error: no more overlay data possible\n");
               error = DBMS_TOO_MUCH_OVERLAY;
            } else {

               /* ------------------------------------------------------- */
               /* Lege neue Overlayseite an und fuege dort die R-Komp. des */
               /* neuen Satzes ein.                                       */
               /* ------------------------------------------------------- */
               if (LogFile!=NULL) fprintf(LogFile,"store data on new overlay page\n");
               table->page_buf = overlay2;
               table->active   = getFreeFilePosition( table );
               newPage( table, OP );
               put_R_comp(table->page_buf,R_comp_new,R_size_new,&di_new);

               /* ------------------------------------------------------- */
               /* Wenn noch Platz fuer ein TID auf Seite 0 ist, so wird    */
               /* dies dort angelegt und die neue Seite angehaengt.        */
               /* Sonst muss ein Datensatz der Seite in ein TID umgewandelt*/
               /* werden                                                  */
               /* ------------------------------------------------------- */
               if ( p0->data - p0->free > sizeof(OINFO) ) {
                  if (LogFile!=NULL) fprintf(LogFile,"store OINFO on overlay root\n");
                  error = write_page( table );
                  if ( error ) {
                     if (LogFile!=NULL) fprintf(LogFile,"error writing new level 1 page (%ld)\n", table->active );
                  } else {
                     if (LogFile!=NULL) fprintf(LogFile,"new level 1 page written (%ld)\n", table->active );
                     oi0[p0->oic].lognum = table->active;
                     oi0[p0->oic].position = 0;
                     oi0[p0->oic].type = DIOVERLAY;
                     (p0->oic)++;
                     p0->free += sizeof(OINFO);
                  }

               } else {

                  /* ---------------------------------------------------- */
                  /* Suche Datensatz auf Overlayseite 0, der in TID       */
                  /* umgewandelt werden kann. Wenn dieser gefunden wurde, */
                  /* werden dessen Daten auf die neue Overlayseite kopiert*/
                  /* und dann auf der Seite 0 geloescht.                   */
                  /* ---------------------------------------------------- */

                  for ( n0=0; n0 < p0->oic; n0++ ) {
                     if ( ! (oi0[n0].type & DIOVERLAY) ) break;
                  }
                  if ( n0 == p0->oic ) {
                     if (LogFile!=NULL) fprintf(LogFile,"fatal error: no OINFO to convert\n");
                  } else {
                     if (LogFile!=NULL) fprintf(LogFile,"convert record %u to OINFO\n", n0 );

                     data = overlay1 + oi0[n0].position;
                     size = oi0[n0].size;
                     di_new.type    = oi0[n0].type;
                     di_new.lognum  = oi0[n0].lognum;
                     put_R_comp( table->page_buf, data, size, &di_new );
                     error = write_page( table );

                     if ( error ) {
                     } else {
                        oi0[n0].lognum = table->active;
                        oi0[n0].type   = DIOVERLAY;
                        source   = overlay1 + p0->data;
                        dest     = source + size;
                        blk_size = oi0[n0].position - p0->data;
                        p0->data += size;

                        if ( blk_size > 0 ) {
                           memory_copy( dest, source, blk_size );
                           for ( n1 = n0+1; n1 < p0->oic; n1++)
                              if (oi0[n1].position > 0)
                                 oi0[n1].position += size;
                        }
                     }
                  }
               }
            }
         }
      }

      /* ----------------------------------------------------------------- */
      /* Setze Seite 0 des Overlaybaumes aktiv, hole freie Position und    */
      /* schreibe die Seite.                                               */
      /* ----------------------------------------------------------------- */
      if ( ! error ) {
         table->page_buf = overlay1;
         table->active = getFreeFilePosition( table );
         error = write_page( table );
         if (LogFile!=NULL) {
            if ( error ) {
               fprintf(LogFile,"error writing overlay root page (%ld)\n", table->active );
            } else {
	       fprintf(LogFile,"overlay root page written (%ld)\n", table->active );
            }
         }
      }

      /* ----------------------------------------------------------------- */
      /* Wenn die Operation bis hierhin erfolgreich ausgefuehrt wurde, so   */
      /* wird im tid der Datenseite, das zu dem geaenderten Overlaybaum    */
      /* gehoert der Offset der zugehoerigen Overlay-Seite geaendert und die  */
      /* bisherige Seite 0 des Overlay-Baumes als zu loeschen eingetragen   */
      /* ----------------------------------------------------------------- */
      if ( ! error ) {
         di[num].lognum  = table->active;
         pages_to_del[*del_cnt] = op0;
         (*del_cnt)++;
      }
   }

   /* -------------------------------------------------------------------- */
   /* Gebe Speicher fuer die Seitenpuffer wieder frei und setze den         */
   /* Seitenpuffer und den Offset der aktiven Seite wieder zurueck auf      */
   /* die Eingangswerte.                                                   */
   /* -------------------------------------------------------------------- */
   memory_free( overlay1 );
   memory_free( overlay2 );
   table->page_buf = old_p;
   table->active = old_a;

   return( error );
}

/* ----------------------------------------------------------------------- */
/* Name       put_R_comp                                                   */
/* Definition void put_R_comp( UCHAR *overlay, UCHAR *r_comp,              */
/*                                         UINT r_size, DINFO *di )        */
/* Funktion   Speichere den Infoteil des Datensatzes auf der Overlay-      */
/*            Seite.                                                       */
/* ----------------------------------------------------------------------- */
void put_R_comp( UCHAR *overlay, UCHAR *r_comp, UINT r_size, DINFO *di ) {

   OPAGE *p  = (OPAGE*)(overlay + sizeof(PTYPE));
   OINFO *oi = (OINFO*)(overlay + sizeof(PTYPE) + sizeof(OPAGE));

   /* -------------------------------------------------------------------- */
   /* Setze den Zeiger auf den Beginn des Datenblocks um die Groesse der     */
   /* einzufuegenden Komponente nach 'vorne' und kopiere die neue Komp. auf */
   /* diese Position fuelle dann das OINFO mit den noetigen Daten aus        */
   /* -------------------------------------------------------------------- */
   p->data -= r_size;
   memory_copy( overlay + p->data, r_comp, r_size );
   oi[p->oic].position = p->data;
   oi[p->oic].lognum   = di->lognum;
   oi[p->oic].size     = r_size;
   oi[p->oic].type     = di->type;

   /* -------------------------------------------------------------------- */
   /* Setze den Zeiger auf den freien Platz der Seite um ein OINFO weiter  */
   /* nach 'hinten' underhoehe die Anzahl der OINFO's um eins               */
   /* -------------------------------------------------------------------- */
   p->free += sizeof(OINFO);
   (p->oic)++;

   return;
}

/* ----------------------------------------------------------------------- */
/* Ende von INSERTO.C                                                      */
/* ----------------------------------------------------------------------- */
