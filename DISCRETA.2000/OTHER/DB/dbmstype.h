/*
	dbmstype.h

   MODUL:        DB_TYPES.H

   BESCHREIBUNG: Headerfile mit Typdeklarationen

   INHALT:       TRIGHTS   - Zugriffsrechte fuer eine Tabellenklasse
           INDEX     - Beschreibung einer Indexkomponente
           DINFO     - Verwaltungsblock fuer einen Datensatz
           BOUND     - Darstellung der Begrenzung fuer Suchvorgaenge
           TABLE     - Zentrale Verwaltungsstruktur fuer eine Tabelle
           BASE_PAGE - Informationen der Basisseite
           PTYPE     - Art und Status einer Seite
           DPAGE     - Verwaltungsstruktur einer Datenseite
           IPAGE     - Verwaltungsstruktur einer Indexseite
           BPOS      - Verwaltung der Grenzen auf einer Indexseite
           DNODE     - Knoten auf einer Datenseite
           INODE     - Knoten auf einer Indexseite
           PAGEBUF   - Verwaltung des Seitenpuffer beim Einfuegen/Aendern
           PLIST     - Seitenoffset-Liste beim Suchen
           DPS       - Ergebnis von Datenseitenspaltung

   VERWEIS:      Datenbankprojekt / Typdeklarationen

   AENDERUNG:    28.04.91 Juergen Mueller
	zusaetzliche Kommentare von Anton Betten 1995

*/


typedef struct _trights   /* access limitation for tables */
{
      UINT insert_denied :1; /* deny insert operations               */
      UINT insert_limit  :3; /* insert only if permission >= limit   */
      UINT search_denied :1; /* deny insert operations               */
      UINT search_limit  :3; /* insert only if permission >= limit   */
      UINT change_denied :1; /* deny insert operations               */
      UINT change_limit  :3; /* insert only if permission >= limit   */
      UINT reserved      :4;
}
TRIGHTS;

typedef struct _user_rights { /* user permission                */

      UINT insert   :3;          /* for insert operations          */
      UINT search   :3;          /* for search operations          */
      UINT change   :3;          /* for change / delete operations */
      UINT reserved :7;

} USER_RIGHTS;

/*
 * Das Format der Datensaetze:
 * 
 * UINT pos[0] : gesamte L"ange des Datensatzes 
 *        (samt Schluessel und Pointer)
 * UINT pos[1...table->index_c] : die Offsets der Schluessel
 * UINT pos[table->index_c + 1] : der Beginn des Datensatzes 
 *        (sogenannter INFO-Teil)
 * hier kommen die Schluessel
 * hier kommt der INFO-Teil
 * 
 */

/*
   Typ: INDEX

   Enthaelt alle Informationen fuer vergleichbare Datensatzkomponenten.
   length ist die Anzahl der internen Komponenten einer Datensatzkomponente. 
   z.B. fuer Strings der Laenge 10 Zeichen, 10. Wenn die Anzahl variabel ist,
   so ist length 0. position gibt den Offset der Datensatzkomponente relativ 
   zum Datensatzanfang an. Bei Tabellen mit variabler Satzlaenge ist position
   0. cmpfunc ist ein Zeiger auf die Funktion zum Vergleich von Datensaetzen 
   bzgl. der Komponente ( siehe USER_M.C ). 
   Die Nummer ist nur fuer Verwaltungszwecke relevant. */

typedef struct _index_db {

      UCHAR  special;	/* AUTOINCREMENT, AUTODECREMENT */
      SINT   number;
      UINT   length;
	/* Anzahl der internen Komponenten des Vgl.Wertes  */

      UINT   position;
	/* Position des Vgl.wertes rel. zum Satzanfang     */

      SINT (*cmpfunc)(void *,void *,UINT);
	/* zugehoerige Vergleichsfkt.*/
	/* realisiert in tdefini.c:
	 * int ubyte_cmp( void *a, void *b, unsigned );
	 * usw. 
	 * liefert -1 genau dann wenn a < b,
	 *          0                 a = b, 
	 *          1                 a > b.
	 */

} INDEX;



/*   BOUND  -  Structure required for search operations */

typedef struct _bound {

      char  *mode;   /* mode of handling the boundary */
      char  *low;    /* pointer to lower boundary     */
      char  *up;     /* pointer to upper boundary     */

} BOUND;

#define LBEQ   1
#define UBEQ   2
#define LBUL   4
#define UBUL   8



   /* TABLE  */
   /* Verwaltungsstruktur der Tabellen fuer saemtliche Operationen */

typedef struct _table {

      int           usage;

      UINT          index_c;     /* Anzahl der Indexkomponenten  */
      INDEX        *index_db;    /* Feld mit Indexdefinitionen   */

      UINT          data_length;
		/* Satzlaenge, (0 hei"st variable L"ange) */
		/* in init.c init_table(): aus tdf.reclen */

      CHAR         *filename;
		/* Name der Tabellendatei  */
      FILE         *fileptr;
		/* Filestream fuer Tabellendatei */

      UCHAR         base[BASE_PAGE_SIZE];
		/* Puffer fuer Basisseite */
		/* mit Offset sizeof(BASE_PAGE) folgt del_pages, 
		 * eine Liste freizugebender Seiten (fio.c) 
		 * die Anzahl der Eintraege steht in base->del_c */
		/* getFreeFilePosition() in baserout.c 
		 * vergibt neue Seiten zuerst aus dieser Liste; 
		 * nur wenn diese leer ist, wird die neue Seite 
		 * an base->file_end vergeben und 
		 * base->file_end wird um base->page_size erhoeht. */

      UCHAR        *page_buf;
		/* Puffer fuer andere Seiten  */

      ULONG         active;
		/* Seitenadresse fuer read/write, 
		 * read_page, write_page schreiben an dieser Position */
		/* Offset in der Datenbank in bytes */

      UCHAR        *data_buf;
		/* Puffer fuer Datensatz */
		/* hier hinein wird von fio.c gelesen */
		/* ins_data() in inserth.c schreibt 
		 * den Datensatz von hier aus auf die Datenseite. */

      char          ditype;        /* Zugriffsrechte fuer Datensatz */
      ULONG         counter;       /* Temporaerer Zaehler (intern)   */
      int           operation;     /* Derzeit aktive Operation     */
      CHAR          read_access;   /* Zugriffsbeschraenkung Lesen   */
      CHAR          change_access; /* Zugriffsbeschraenkung Aendern  */
      USER_RIGHTS   user_rights;   /* Zugriffsrechte des Benutzers */

} TABLE;


/* BASE_PAGE   */
/* Struktur zur Speicherung der permanenten Tabellendaten */


typedef struct _base_page {

      CHAR     crjm[36];
		/* Copyright  */
      CHAR     variable;
		/* Flag fuer Satzlaenge variabel */
      UINT     max_data_length;
		/* Groesse des laengsten Datensatzes */
		/* wird beim Einlesen (fio.c) 
		 * auf das Maximum aller Datensatzl"angen gesetzt. */
      ULONG    version;
		/* Aktuelle Versionsnummer der Tabelle */
      UINT     defnum;
		/* number of table's definition */
      CHAR     info[TBLINF];
		/* Infotext */

      ULONG    file_end;
		/* Adresse des Dateiendes */
		/* erhoeht nur (?) durch getFreeFilePosition(), 
		 * wenn keine geloeschten Seiten mehr vorhanden sind. */

      UINT     level;
		/* Hoehe des Baumes */
		/* wird in splitd.c d_new_root() erhoeht. */

      ULONG    root_page;
		/* Adresse der Wurzelseite */
		/* an diesem Wert fuer table->active 
		 * beginnt jeder Zugriffsweg (fio.c ) */
		/* wird veraendert beim Abspeichern 
		 * jedes Zugriffsweges (fio.c ) */
		
      long     ovl_root;
      ULONG    lognum;
		/* Naechste logische Satznummer */

      UINT     del_c;
		/* Anzahl gloeschter Seiten 
		 * (in table->base[sizeof(BASE_PAGE)]) */
	
      UINT     page_size;
		/* Seitengroesse in Bytes  */
		/* so viel bytes werden von read_page, write_page 
		 * gelesen bzw. geschrieben. */ 

      ULONG    data_c;
		/* Anzahl Datensaetze in der Tabelle */
      long     counter[BASCNT];
		/* Interne Zaehler */

} BASE_PAGE;


/*   PTYPE  -  Structure at the beginning of all pages */

typedef struct _ptype {

      unsigned type     :2;   /* see definitions below */
      unsigned reserved :14;

} PTYPE;

#define PTYPEDATA     1
#define PTYPEOVERLAY  2
#define PTYPEINDEX    3


/* DPAGE  -  Structure behind PTYPE at each data page  */
/* DNODE  -  Structure for tree nodes on data pages */
/* DINFO  -  Structure for each record on data pages  */

/* Der Aufbau einer Datenseite:
 * 
 * ganz vorne genau eine DPAGE Struktur.
 * 
 * Dann kommen die DNODE Strukturen, welche einen binaeren Baum 
 * bilden.
 * 
 * deren Ende ist der Beginn des freien Speicherbereichs (free Pointer) 
 * 
 * Am oberen Ende (von oben nach unten (vorne) wachsend 
 * die DINFOs und die Datensaetze selber. 
 * Der data Pointer zeigt auf den Beginn der DNODE Strukturen.
 * 
 * Der Bereich zwischen free und data pointer ist frei.
 * 
 * Jeder Datensatz hat ein DINFO, wo in position sein Offset 
 * vermerkt ist. 
 * Die DINFO Strukturen werden beim Einfuegen 
 * eines neuen Datensatzes nach vorne kopiert 
 * und wachsen dann 'am Ende';
 * Neue Datensaetze werden vor den bisherigen Datensaetzen 
 * angelegt. Datensatz Null ist und bleibt somit 
 * immer am oberen Ende der Datenpage. */

typedef struct _dpage {

      UINT   dic;        /* Count of DINFO's = Count of records */
      UINT   dnc;        /* Count of DNODE's                    */

      UINT   data;       /* Position of first DINFO             */
	/* Beginn der DINFOs, am oberen Ende der Seite, 
	 * waechst nach vorne (unten). 
	 * Nach den DINFOs kommen die Datensaetze selber 
	 * (variable L"ange) */

      UINT   free;       /* Position where free area starts     */
	/* Beginn freier Speicher auf der Datenseite, 
	 * nach den DNODEs */

} DPAGE;


typedef struct _dnode {

      UINT ds;      
	/* Nr. des Datensatz, der den Vergleichswert enthaelt    */
	/* = Index in das DINFO - Array */
 
      UINT comp;
	/* Komponente des Vergleichswertes */
	/* der i-te Schluessel ist also durch 
	 * 
	 * proj (table, i, table->page_buf + di[ dn[x].ds ].position)
	 * 
	 * erhaeltlich
	 * (di ist Pointer auf die DINFO Strukturen, 
	 * dn zeigt auf die DNODE Strukturen). */

      UINT smaller; /* Kleinerer Nachfolger */
      UINT greater; /* Groesserer Nachfoler */

      UINT total;
	/* Groesse der Datensaetze im Baum mit Knoten als Wurzel */
	/* = Summe der Datensatzlaengen aller Datensaetze 
	 * im Teilbaum; wird benoetigt um beim Seiten - Splitting 
	 * die 1/3 - 2/3 (oder besser) Teilung zu finden. */

      char mode;
	/* Flags zum Knoten s.u. */

} DNODE;

/* fuer mode: */

#define SMALLER_DS      1
#define NOT_SMALLER_DS  254
#define GREATER_DS      2
#define NOT_GREATER_DS  253
#define EQ_SMALLER      4


typedef struct _dinfo {

      UINT   position;  /* Position of the record on the data page   */
      long   lognum;    /* Logical record number / offset of overlay */
      char   type;      /* bit 0-2 read 3-6 change access, 7 overlay */

} DINFO;

#define DIOVERLAY 128



/* IPAGE (Verwaltungsstruktur fuer Indexseiten) */

/* Aufbau einer Indexseite:
 * 
 * Zuerst genau eine IPAGE Struktur.
 * 
 * Dann kommen INODE Strukturen, 
 * IPAGE->node_c ist die Anzahl der INODE Strukturen.
 * Die INODE-Strukturen bilden einen 
 * bin"aren Baum, die Wurzel ist der 0-te INODE.
 * mit dem Ende dieses Feldes beginnt der freie Speicherplatz.
 * 
 * hier ist also freier Speicherbereich. 
 * (IPAGE->free ist NICHT der Offset, sondern die LAENGE 
 * des freien Speicherbereichs)
 * 
 * Von oben her nach vorne wachsend folgt dann:
 *
 * mit Offset olist:
 * ULONGS, welche Offsets von anderen Index bzw. Daten-Seiten angeben.
 * es handelt sich hierbei um die Seiten-Offsets in der 
 * Datendatei, d.h. Werte, die nach table->active gehen, 
 * bevor load_page() gerufen wird.
 * IPAGE->o_c ist die Anzahl der Nachfolgerseiten (L"ange von olist[])
 * 
 * unmittelbar darauf:
 * mit Offset bpos:
 * Das Feld der BPOS Strukturen; hier sind die Schranken 
 * eingetragen; die Werte der Schranken selber stehen hier 
 * nicht, da Schranken unterschiedliche Wertebereiche haben 
 * z.B. int- Variablen und char Arrays;
 * IPAGE->b_c ist die Anzahl der BPOS Strukturen.
 * Die BPOS Struktur:
 *   counter ist die Anzahl der Verweise auf diese Schranke (?)
 *   pos ist der Offset in der Seite der Schranke (s.u.)
 *   comp ist der Index (table->index_db[]) der Schranke
 *   size ist die Schrankengroesse in bytes. 
 * 
 * unmittelbar darauf:
 * (kein Offset in der IPAGE Struktur hierf"ur vorgesehen, 
 * berechnet sich aus IPAGE->bpos + IPAGE->b_c * sizeof(BPOS) )
 * die Schranken selber; 
 * Offsets auf die jeweilige Schranke sind 
 * in der BPOS Struktur zu finden (pos)
 * Die Schranken werden von ober herab angelegt 
 * (Wachstum nach vorne bzw. unten). 
 * Ganz oben die erste Schranke.
 */

typedef struct _ipage {

      UINT node_c;  /* Anzahl der INODE's                                  */
      UINT olist;   /* Position des 1. Eintrags der Nachfolgeroffset-Liste */
      UINT o_c;     /* Anzahl der Nachfolgeseiten                          */
      UINT b_c;     /* Anzahl der Grenzen                                  */
      UINT bpos;    /* Position des 1. Eintrages der Liste der Grenzen     */
      UINT free;    /* Groesse des freien Platzes auf der Seite            */

} IPAGE;


typedef struct _inf {
      UINT smaller_extern :1;
		/* wenn 1: Kleiner Nachfolger 
		 * ist auf (externer) Nachfolgerseite */
      UINT greater_extern :1;
		/* wenn 1: Groesserer Nachfolger 
		 * ist auf (externer) Nachfolgerseite     */
		/* falls 1, so handelt es sich um einen sog. EXIT-Knoten;
		 * d.h. der Nachfolger ist eine andere Datenseite;
		 * bei INODE->smaller bzw. INODE->greater handelt es sich 
		 * dann um einen Index in das olist[] - Array; 
		 * dort ist dann der Seitenoffset (table->active) 
		 * der Nachfolgerseite auszulesen. 
		 * Das Durchlaufen eines INODE- Baumes geschieht 
		 * z.Bsp. in inserth.c - get_next_page() */

      UINT equal_smaller  :1;
		/* Bei Gleichheit weiter mit kleinerem Nachfolger    */
		/* wenn equal_smaller 1: 
		 * bei Gleichheit gehe in den smaller Zweig. */
		/* aus dps->flag[i].l_eq_s */
	
      UINT low            :1;
		/* wenn 1: Knoten ist Untergrenze 
		 * des extrahierten Bereiches */
      UINT reserved       :12;
} INF;


typedef struct _inode {
      UINT bnum;
	/* Nummer des Vergleichswertes in der BPOS-Liste */
	/* also ein index in das BPOS Array */

      UINT smaller;    /* Kleinerer Nachfolger */
      UINT greater;    /* Groesserer Nachfolger */
	/* hier sind f"ur die beiden Nachfolger Indizes in das 
	 * INDOE Array zu finden. 
	 * wenn das entsprechende smaller_extern bzw. greater_extern
	 * Flag gesetzt ist, so handelt es sich um 
	 * einen EXIT-Knoten; smaller bzw. greater 
	 * sind dann ein Index in das olist[] Array. */

      UINT cnt;        /* Anzahl der Knoten fuer Grenzen des extr. Ber. */
      INF  flag;
} INODE;


typedef struct _bpos {
      UINT counter;  /* Anzahl der Verweise auf die jeweilige Grenze */
      UINT pos;      /* Position der Daten der Grenze relativ zum Seitenanfang */
      UINT comp;     /* Komponente der Grenze */
      UINT size;     /* Groesse der Grenze in Bytes */
} BPOS;


/*   OPAGE  -  Structure on every overlay page behind PTYPE */
/*                                                                          */
/*   OINFO  -  Structure to every record part on overlay pages              */
/*           Positioned behind OPAGE                                        */

typedef struct _opage {

      UINT  oic;    /* Count of overlay data on page    */
      UINT  free;   /* Begin of free area               */
      UINT  data;   /* Begin of data area               */
      ULONG next;   /* Offset of next page (0 if none ) */

} OPAGE;

typedef struct _oinfo {

      UINT  position; /* Position of record part on overlay page */
      ULONG lognum;   /* Logical record number                   */
      UINT  size;     /* Size of record part in bytes            */
      char  type;     /* Flags of the record                     */

} OINFO;



/* PAGEBUF  -  Structure required in insert, delete and update operations */

typedef struct _pagebuf {

      char  *page;
		/* pointer to page stored at this buffer position */
		/* malloc allocated memory */
      UINT   next;
		/* Index in olist[] der Indexseite (in page), 
		 * wo die Nachfolgerseite zu finden ist.
		 * wird in fio.c von der Funktion 
		 * get_next_page() gesetzt. */
		/* wird benoetigt, um beim schreiben der evtl. 
		 * geaenderten Seiten die neuen 
		 * Nachfolgeroffsets eintragen zu koennen 
		 * (siehe fio.c). */
		
		/* Position in olist array where page was left */

      ULONG  offset;
		/* offset in der Datendatei (Beginn der Seite in bytes) */
		/* Offset of page at this buffer position */
		/* wird beim Abspeichern der Seiten in fio.c 
		 * mittels getFreeFilePosition() auf 
		 * neue Seitenadressen gesetzt 
		 * (der Zugriffsweg wird rueckwaerts wieder 
		 * abgespeichert). */

} PAGEBUF;



/* PLIST  -  Structure for page offset 'stack' in search function */

typedef struct _plist {

      ULONG *list;
      UINT   act;
      UINT   max;

} PLIST;

#define PLISTCNT  50



/*   DPS (Ergebnis des Datenseitenspaltens) */

typedef struct _dpsf {
      UINT u_eq_s     :1;
		/* bei Knoten in up ist = beim kleineren Nachfolger  */
      UINT l_eq_s     :1;
		/* bei Knoten in low ist = beim kleineren Nachfolger */
		/* EQ_SMALLER flag aus DNODE ->mode */

      UINT no_low     :1;      /* keine untere Grenze bzgl. dieser Komponente       */
      UINT no_up      :1;      /* keine oberer Grenze bzgl. dieser Komponente       */

      UINT low_exists :1;      /* untere Grenze existiert schon (d_insert_index)    */
      UINT up_exists  :1;      /* obere Grenze existiert schon (d_insert_index)     */
		/* wenn 1, dann steht in DPS->bpl[] bzw. bpu[]
		 * die Nummer der bereits existierenden Grenze (BPOS Index) */

      UINT last       :1;
		/* 1 bei der letzten Komponente, fuer die mindestens 
		 * eine Schranke existiert */

      UINT reserved   :9;
} DPSF;

typedef struct _dps {
      UINT   cnt; /* # Grenzen */

      DPSF   flag[INDEXC];

      UINT   upos[INDEXC+1];
      UINT   lpos[INDEXC+1];
	/* Offsets der Grenzen in up, low */
	/* upos[index_c] = gesamte L"ange der Grenzen in 
	 * up, low */

      UINT   bpl[INDEXC];
      UINT   bpu[INDEXC];
	/* hier kommen in put_bounds() die 
	 * BPOS Indizes rein, welche f"ur die neu 
	 * angelegten Schranken zustaendig sind. */

      UCHAR *up;
      UCHAR *low;
	/* Speicher der Grenzen, 
	 * Grenzen direkt aufeinanderfolgend, 
	 * Offsets (Beginn) wie in upos[], lpos[] */

      ULONG  in_page;
	/* Adresse (Offset) der ersten gespalteten Seite */
      ULONG  ext_page;
	/* zweite gespaltete Seite */
} DPS;


typedef struct _ipsnf {
      UINT next_node_smaller :1;
      UINT other_side_extern :1;
      UINT other_side_next   :1;
      UINT reserved          :13;
} IPSNF;

typedef struct _ipsn {
      UINT   node;
      IPSNF  flag;
      UINT   other_next;
} IPSN;

typedef struct _ips {
      UCHAR  *old_page;
      UINT    cnt;
      IPSN   *ipsn;
      UINT    last_int_s;
      ULONG   ext_page;
      ULONG   int_page;
} IPS;


typedef struct _tdf {

   char     status;
   char     read_access;
   char     change_access;
   unsigned psize;
   unsigned reclen;
   unsigned icnt; /* # indices */
   unsigned rpos;
   unsigned rsize;
   unsigned pos[INDEXC];
   unsigned elem[INDEXC];
   unsigned type[INDEXC];
   unsigned flags[INDEXC];
   char     name[INDEXC][10];

} TDF;

/* zur Initialisierung siehe 
 * tdefini.c : get_definition() */
