/* bt.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef DB_TRUE

#include <DISCRETA/db.h>

#ifdef SYSTEMMAC
#include <console.h>
#include <time.h>
#include <unix.h>
#endif

INT bt_debug = FALSE;

#undef DEBUG_BT_CREATE
#undef DEBUG_SEARCH_PAGE

/* 
 * BAYERTREE
 */

#define BTREEMAXPAGESIZE 32
#define BTREEHALFPAGESIZE  16

/* Dateiformat:
 * In Block 0 sind AllocRec/NextFreeRec/RootRec gesetzt.
 * Block 1..AllocRec sind Datenpages.
 * Die freien Bloecke sind ueber NextFreeRec verkettet.
 * Der letzte freie Block hat NIL als Nachfolger.
 * Dateigroesse = (AllocRec + 1) * sizeof(PageTyp) */

typedef struct itemtyp {
	KEYTYPE Key;
	DATATYPE Data;
	INT4 Childs; /* Anzahl der Nachfolger ueber Ref */
	INT4 Ref;
} ItemTyp;

typedef struct pagetyp {
	INT4 AllocRec;
	INT4 NextFreeRec;
	INT4 RootRec;

	INT4 NumItems;
	ItemTyp Item[BTREEMAXPAGESIZE + 1];
/* Item[0]           enthaelt keine Daten, 
 *                   nur Ref/Childs ist verwendet.
 * Item[1..NumItems] fuer Daten und 
 *                   Ref/Childs verwendet. */
} PageTyp;

typedef struct buffer {
	INT4 PageNum;
	INT4 unused;
	PageTyp Page;
	long align;
} Buffer;

static INT WRITE1(FILE *f, PageTyp *page);
static INT READ1(FILE *f, PageTyp *page);
static INT WriteInfo(BAYERTREE_OP p);
/* Schreibt die Variablen 'AllocRec', 'FreeRec' und
 * 'Root' in die 0te-Seite der Datenbank */
static INT AllocateRec(BAYERTREE_OP p, INT *x);
static INT ReleaseRec(BAYERTREE_OP p, INT x);
static INT LoadPage(BAYERTREE_OP p, Buffer *BF, INT x);
static INT SavePage(BAYERTREE_OP p, Buffer *BF);
static INT SearchPage(BAYERTREE_OP p, Buffer *buffer, void *pSearchKey, 
	DATATYPE *pSearchData, INT *cur, INT *x, INT *Found);
static INT SearchBtree(BAYERTREE_OP p, INT page, INT *idx, INT *f_found);
static INT page_i_th(BAYERTREE_OP p, INT l, Buffer *buffer, 
	INT *cur, INT *i, INT *found);
static INT Update(BAYERTREE_OP p, INT Node, INT *Rise, ItemTyp *RisenItem, 
	INT *RisenNeighbourChilds);
static void bt_item_copy(ItemTyp *a, ItemTyp *b);
static INT Split(BAYERTREE_OP p, Buffer *BF, ItemTyp *Item, 
	INT x, INT *RisenNeighbourChilds);
static INT Del(BAYERTREE_OP p, INT Node, INT *Underflow);
static INT FindGreatest(BAYERTREE_OP p, INT Node1, INT *Underflow, 
	Buffer *DKBF, INT x);
static INT Compensate(BAYERTREE_OP p, INT Precedent, 
	INT Node, INT Path, INT *Underflow);
static void page_print(BAYERTREE_OP p, Buffer *BF);
static void fpage_print(BAYERTREE_OP p, Buffer *BF, FILE *fp);
// static void item_print(BAYERTREE_OP p, ItemTyp *item);

/* check_page_datref */

INT BAYERTREE_OB::field_name(INT i, INT j, BYTE *str)
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "f_duplicatekeys"; break;
	case 1: s = "nb_bt_key"; break;
	case 2: s = "bt_key"; break;
	case 3: s = "filename"; break;
	case 4: s = "f_open"; break;
	case 5: s = "stream"; break;
	case 6: s = "buf_idx"; break;
	case 7: s = "Root"; break;
	case 8: s = "FreeRec"; break;
	case 9: s = "AllocRec"; break;
	case 10: s = "btree_idx"; break;
	default:
		return error("BT::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

INT BAYERTREE_OB::init(BYTE *file_name, INT f_duplicatekeys, INT btree_idx)
{
	INT erg = OK;
	
	erg += m_il(11);
	c_obj_k(BAYERTREE);
	s_f_duplicatekeys()->m_i(f_duplicatekeys);
	s_nb_bt_key()->m_i(0);
	s_bt_key()->m_il(VECTOR_OVERSIZE);
	s_filename()->init(file_name);
	s_f_open()->m_i(0);
	s_stream()->m_i(0);
	s_buf_idx()->m_i(0);
	s_Root()->m_i(0);
	s_FreeRec()->m_i(0);
	s_AllocRec()->m_i(0);
	s_btree_idx()->m_i(btree_idx);
	return erg;
}

INT BAYERTREE_OB::add_key_INT4(INT field1, INT field2)
{
	BT_KEY_OB bt_key;
	INT erg = OK;
	
	erg += bt_key.init_INT4(field1, field2);
	erg += s_bt_key()->append_element(s_nb_bt_key(), &bt_key);
	return erg;
}

INT BAYERTREE_OB::add_key_INT2(INT field1, INT field2)
{
	BT_KEY_OB bt_key;
	INT erg = OK;
	
	erg += bt_key.init_INT2(field1, field2);
	erg += s_bt_key()->append_element(s_nb_bt_key(), &bt_key);
	return erg;
}

INT BAYERTREE_OB::add_key_STRING(INT output_size, 
	INT field1, INT field2)
{
	BT_KEY_OB bt_key;
	INT erg = OK;
	
	erg += bt_key.init_STRING(output_size, field1, field2);
	erg += s_bt_key()->append_element(s_nb_bt_key(), &bt_key);
	return erg;
}

INT BAYERTREE_OB::add_key_INT4_VEC(INT field1, INT field2, 
	INT vec_fst, INT vec_len)
{
	BT_KEY_OB bt_key;
	INT erg = OK;
	
	erg += bt_key.init_INT4_VEC(field1, field2, vec_fst, vec_len);
	erg += s_bt_key()->append_element(s_nb_bt_key(), &bt_key);
	return erg;
}

INT BAYERTREE_OB::add_key_INT2_VEC(INT field1, INT field2, 
	INT vec_fst, INT vec_len)
{
	BT_KEY_OB bt_key;
	INT erg = OK;
	
	erg += bt_key.init_INT2_VEC(field1, field2, vec_fst, vec_len);
	erg += s_bt_key()->append_element(s_nb_bt_key(), &bt_key);
	return erg;
}

INT BAYERTREE_OB::sprint(BYTE *s)
{
	BYTE str[512];
	
	sprintf(str, "BAYERTREE: %s "
		"nb_bt_key = %ld Root = %ld FreeRec = %ld AllocRec = %ld", 
		s_filename_s(), s_nb_bt_key_i(), 
		s_Root_i(), s_FreeRec_i(), s_AllocRec_i());
	if (strlen(s) + strlen(str) < 200)
		strcat(s, str);
	return OK;
}

#define MAX_ROOT_BUF 10

INT f_RootBF_free[MAX_ROOT_BUF];
Buffer *RootBF = NIL;
Buffer *tmpBF = NIL;
	/* used in: bt_open(), WriteInfo(), AllocateRec(), 
	 * ReleaseRec() */

INT database_init(void)
{
	INT size, i;
	
	/* printf("size ItemTyp = %ld\n", (INT)sizeof(ItemTyp));
	printf("size PageTyp = %ld\n", (INT)sizeof(PageTyp)); */
	size = MAX_ROOT_BUF * sizeof(Buffer);
	/* printf("size for root buf = %ld\n", size); */
	RootBF = (Buffer *) my_malloc(size, "bt.C: database_init()");
	if (RootBF == NIL) {
		printf("no memory for RootBF\n");
		return ERROR;
		}
	for (i = 0; i < MAX_ROOT_BUF; i++)
		f_RootBF_free[i] = TRUE;
	f_RootBF_free[0] = FALSE;
	tmpBF = RootBF;
	bt_key_init();
	return (OK);
}

void database_exit(void)
{
	bt_key_exit();
	if (RootBF) {
		my_free(RootBF);
		RootBF = NIL;
		}
}

INT root_buf_alloc(void)
/* returns -1 if no buffer free */
{
	INT i;
	
	for (i = 0; i < MAX_ROOT_BUF; i++) {
		if (f_RootBF_free[i]) {
			f_RootBF_free[i] = FALSE;
			return i;
			}
		}
	return -1;
}

void root_buf_free(INT i)
{
	if (i < 0 || i >= MAX_ROOT_BUF) {
		printf("root_buf_free()|i illegal");
		return;
		}
	f_RootBF_free[i] = TRUE;
}

/* OS-Interface: */

FILE *OPEN1(BYTE *name)
{
	FILE *f = fopen(name, "r+b");
	
	if (f == NIL) {
		printf("OPEN1()|error in fopen");
		}
	return(f);
}

FILE *CREATE1(BYTE *name)
{
	FILE *f = fopen(name, "w+b");
	
	if (f == NIL) {
		printf("CREATE1()|error in fopen");
		}
	return(f);
}

INT CLOSE1(FILE *f)
{
	if (fclose(f) != 0) {
		printf("CLOSE1()|error in fclose");
		return ERROR;
		}
	return OK;
}

static INT WRITE1(FILE *f, PageTyp *page)
{
	INT size = (INT)sizeof(PageTyp);
	
	if ((INT)fwrite((BYTE *)page, 1, size, f) != size) {
		return error("WRITE1() error in fwrite()");
		}
	return OK;
}

static INT READ1(FILE *f, PageTyp *page)
{
	INT size = (INT)sizeof(PageTyp);
	
	if ((INT) fread((BYTE *)page, 1, size, f) != size) {
		return error("READ1() error in fread()");
		}
	return OK;
}

INT SEEK1(FILE *f, INT l)
{
	INT off;
	
	off = l * (INT)sizeof(PageTyp);
	if (f_seek_set(f, off) != 0) {
		return error("SEEK1() error in f_seek_set()");
		}
	return OK;
}

INT BAYERTREE_OB::create()
{
	FILE *f;
	INT root_idx;
	Buffer *Root_BF;
	
	if (s_f_open_i() == TRUE) {
		return error("BT::create() already open");
		}
#ifdef DEBUG_BT_CREATE
	printf("in BT::create()\n"); fflush(stdout);
#endif
	f = CREATE1(s_filename_s());
	if (f == NIL) {
		return error("BT::create() error in CREATE1()");
		}
	s_stream()->m_i((INT)f);
	root_idx = root_buf_alloc();
	if (root_idx == -1) {
		return error("BT::create() no free root buffer");
		}
	Root_BF = RootBF + root_idx;
	fill_char(Root_BF, (INT)sizeof(Buffer), 0);
	s_buf_idx()->m_i(root_idx);
	s_FreeRec()->m_i(0);
	s_AllocRec()->m_i(0);
	s_Root()->m_i(0);
	s_f_open()->m_i(TRUE);
	if (WriteInfo(this) != OK) {
		return error("BT::create() error in WriteInfo()");
		}
#ifdef DEBUG_BT_CREATE
	printf("in BT::create() nach WriteInfo\n");
	fflush(stdout);
#endif
	return OK;
}

INT BAYERTREE_OB::open()
{
	FILE *f;
	Buffer *BF = tmpBF;
	INT root_idx;
	Buffer *Root_BF;
	
	if (s_f_open_i()) {
		return error("BT::open() already open");
		}
	f = OPEN1(s_filename_s());
	if (f == NIL) {
		return error("BT::open() error in OPEN1()");
		}
	s_stream()->m_i((INT)f);
	root_idx = root_buf_alloc();
	if (root_idx == -1) {
		return error("BT::open() no free root buffer");
		}
	Root_BF = RootBF + root_idx;
	fill_char(Root_BF, (INT)sizeof(Buffer), 0);
	s_buf_idx()->m_i(root_idx);

	SEEK1(s_stream_i(), 0);
	READ1(s_stream_i(), &BF->Page);
	s_FreeRec()->m_i(BF->Page.NextFreeRec);
	s_AllocRec()->m_i(BF->Page.AllocRec);
	s_Root()->m_i(BF->Page.RootRec);
#if 0
	{
	BYTE str[256];
	
	sprintf(str, "FreeRec = %ld AllocRec = %ld Root = %ld", 
		s_FreeRec_i(), s_AllocRec_i(), s_Root_i());
	Srfs("BT::open", str);
	}
#endif
	if (s_Root_i() != 0) {
		SEEK1(s_stream_i(), s_Root_i());
		READ1(s_stream_i(), &Root_BF->Page);
		Root_BF->PageNum = s_Root_i();
		}
	else
		Root_BF->PageNum = 0;
	s_f_open()->m_i(TRUE);
#ifdef DEBUG_BT_CREATE
	printf("in BT::open() btree opened\n");
	fflush(stdout);
#endif
	return OK;
}

INT BAYERTREE_OB::close()
{
	if (s_f_open_i() == FALSE) {
		return error("BT::close() not open");
		}
	CLOSE1(s_stream_i());
	root_buf_free(s_buf_idx_i());
	s_f_open()->m_i(FALSE);
#ifdef DEBUG_BT_CREATE
	printf("in BT::close() btree closed\n");
	fflush(stdout);
#endif
	return OK;
}

INT BAYERTREE_OB::len(INT *len)
{
	INT l;
	INT j, pagelen;
	Buffer *Root_BF;
	
	if (s_f_open_i() == FALSE) {
		return error("BT::len() not open");
		}
	Root_BF = RootBF + s_buf_idx_i();
	l = 0;
	pagelen = Root_BF->Page.NumItems;
	/* sprintf(str, "pagelen = %ld", pagelen);
	Srfs("BT::len", str); */
	for (j = 0; j <= pagelen; j++) { /* mit nulltem Index ! */
		l += Root_BF->Page.Item[j].Childs;
		}
	l += pagelen;
	*len = l;
	return OK;
}

/*
 * interne Routinen: 
 */

static INT WriteInfo(BAYERTREE_OP p)
/* Schreibt die Variablen 'AllocRec', 'FreeRec' und 'Root' 
 * als 'AllocRec', 'NextFreeRec' und 'RootRec' 
 * in die 0-te Seite der Datenbank. */
{
	INT size, erg = OK;
	Buffer *BF = tmpBF;
	
	if (p->s_f_open_i() == FALSE) {
		return error("WriteInfo() not open");
		}
	size = sizeof(Buffer);
	fill_char((BYTE *)BF, size, 0);
	BF->Page.AllocRec = p->s_AllocRec_i();
	BF->Page.NextFreeRec = p->s_FreeRec_i();
	BF->Page.RootRec = p->s_Root_i();
	/* printf("WriteInfo()|AllocRec = %ld FreeRec = %ld Root = %ld\n", 
		p->s_AllocRec_i(), p->s_FreeRec_i(), p->s_Root_i()); */
	erg += SEEK1(p->s_stream_i(), 0);
	erg += WRITE1(p->s_stream_i(), &BF->Page);
	return erg;
}

static INT AllocateRec(BAYERTREE_OP p, INT *x)
/* Ein freier Record der Datanbank wird ermittelt
 * --
 * INT *x - Gibt Nummer eines freien Records an */
{
	Buffer *BF = tmpBF;
	INT erg = OK, size;
	
	if (p->s_f_open_i() == FALSE) {
		return error("AllocateRec() not open");
		}
	if (p->s_FreeRec_i() == 0) {
		p->s_AllocRec()->inc();
		*x = p->s_AllocRec_i();
		/* printf("AllocateRec()|p->FreeRec == 0|x = %ld\n", *x); */
		if (WriteInfo(p) != OK) {
			return error("AllocateRec() error in WriteInfo()");
			}
		return OK;
		}
	else {
		size = (INT)sizeof(Buffer);
		fill_char((BYTE *)BF, size, 0);
		*x = p->s_FreeRec_i();
		/* printf("AllocateRec()|x = p->FreeRec = %ld\n", *x); */
		if (SEEK1(p->s_stream_i(), *x) != OK) {
			return error("AllocateRec() error in SEEK1()");
			}
		if (READ1(p->s_stream_i(), &BF->Page) != OK) {
			return error("AllocateRec() error in READ1()");
			}
		p->s_FreeRec()->m_i(BF->Page.NextFreeRec);
		if (WriteInfo(p) != OK) {
			return error("AllocateRec() error in WriteInfo()");
			}
		}
	return erg;
}

static INT ReleaseRec(BAYERTREE_OP p, INT x)
/* Gibt einen Record wieder frei
 * Der Block wird an den Anfang der Frei-Liste eingefuegt.
 * --
 * INT x - Nummer des freizugebenden Records */
{
	INT size;
	Buffer *BF = tmpBF;
	INT erg = OK;
	
	if (p->s_f_open_i() == FALSE) {
		return error("ReleaseRec() not open");
		}
	size = (INT)sizeof(Buffer);
	fill_char((BYTE *)BF, size, 0);
	BF->Page.NextFreeRec = p->s_FreeRec_i();
	if (SEEK1(p->s_stream_i(), x) != OK) {
		return error("ReleaseRec() error in SEEK1()");
		}
	if (WRITE1(p->s_stream_i(), &BF->Page) != OK) {
		return error("ReleaseRec() error in WRITE1()");
		}
	p->s_FreeRec()->m_i(x);
	if (WriteInfo(p) != OK) {
		return error("ReleaseRec() error in WriteInfo()");
		}
	return erg;
}

static INT LoadPage(BAYERTREE_OP p, Buffer *BF, INT x)
/* Laedt eine Seite in den Speicher. Soll die Wurzel ge-
 * laden werden, so wird nur der Puffer kopiert
 * --
 * Buffer *BF - Puffer enthaelt nach Aufruf die Seite
 * INT x      - zu ladende Seite */
{
	Buffer *Root_BF;

	if (p->s_f_open_i() == FALSE) {
		return error("LoadPage() not open");
		}
	Root_BF = RootBF + p->s_buf_idx_i();
	if (x == p->s_Root_i())
		*BF = *Root_BF;
	else {
		BF->PageNum = x;
		if (SEEK1(p->s_stream_i(), x) != OK) {
			return error("LoadPage() error in SEEK1()");
			}
		if (READ1(p->s_stream_i(), &BF->Page) != OK) {
			return error("LoadPage() error in READ1()");
			}
		}
	return OK;
}

static INT SavePage(BAYERTREE_OP p, Buffer *BF)
/* Eine Seite wird auf den Hintergrundspeicher geschrieben
 * --
 * Buffer *BF - Zu speichernde Seite */
{
	Buffer *Root_BF;

	if (p->s_f_open_i() == FALSE) {
		return error("SavePage() not open");
		}
	Root_BF = RootBF + p->s_buf_idx_i();
	if (BF->PageNum == p->s_Root_i())
		*Root_BF = *BF;
	if (SEEK1(p->s_stream_i(), BF->PageNum) != OK) {
		return error("SavePage() error in SEEK1()");
		}
	if (WRITE1(p->s_stream_i(), &BF->Page) != OK) {
		return error("SavePage() error in WRITE1()");
		}
	fflush(p->s_stream_i());
	return OK;
}

static INT SearchPage(BAYERTREE_OP p, 
	Buffer *buffer, void *pSearchKey, 
	DATATYPE *pSearchData, 
	INT *cur, INT *x, INT *Found)
/* Fuehrt binaere Suche innerhalb einer Seite aus.
 * --
 * Buffer *BF            - zu untersuchende Seite.
 * KEYTYPE  *SearchKey   - Zu suchender Schluessel.
 * DATATYPE *pSearchData - optional, nur datref verwendet. 
 *                         Bei gleichen Schluesseln werden 
 *                         zusaetzlich die datref's verglichen.
 * INT *x                - Gibt Position an, falls Suche 
 *                         erfolgreich. Ansonsten naechst 
 *                         kleinerers Element.
 *                         x kann also auch null werden.
 * INT *Found            - Gibt an, ob Schluessel gefunden.
 * 
 * Es wird aufgerufen:
 *   WORD (*cmp_func)(void *key1, void *key2, INT *res);
 * Key1 stammt aus der Datenbank, key2 ist der durchgeschleifte, 
 * zu suchende Schluessel. 
 * Ergebnis:
 *   res < 0:  key1 < key2
 *   res == 0: key1 == key2
 *   res > 0:  key1 > key2
 * Bei gleichen Schluesseln werden intern noch die datref's 
 * verglichen, sofern pSearchData gesetzt ist. 
 * Bei gleichen Schluesseln (und evtl. gleichen datref's) 
 * wird der letzte Eintrag gesucht und zurueckgegeben. 
 * *cur wird erhoeht um die Childs Zaehler auf dieser Seite 
 * (incl. 0-tem Eintrag) bis unmittelbar vor dem gefundenen 
 * Datensatz, und dann + 1, so dass cur die aktuelle 
 * Datensatznummer enthaelt, wenn es vorher im Suchbaum 
 * schon mitgefuehrt wurde. */
{
	ItemTyp *item = NIL;
	INT searchdatref;
	INT childs;
	INT r, l, i;
	INT res;
	
	if (pSearchData != NIL) {
		searchdatref = pSearchData->datref;
		}
	l = 0;
	r = buffer->Page.NumItems + 1;
	*Found = FALSE;
	while (l + 1 < r) {
		*x = ((l + r) >> 1);
		item = &buffer->Page.Item[*x];
#ifdef DEBUG_SEARCH_PAGE
		printf("SearchPage()|l = %ld *x = %ld r = %ld\n", l, *x, r);
		item_print(p, item);
#endif
		if (bt_key_cmp(item->Key.c, (BYTE *) pSearchKey, 
			p->s_bt_key(), p->s_nb_bt_key_i(), &res) != OK) {
			return error("SearchPage() error in bt_key_cmp()");
			}
		if (res == 0) {
			/* wenn pSearchData gesetzt, 
			 * dann bei gleichen Schluesseln 
			 * Suche nach datref */
			if (pSearchData != NIL) {
				if (item->Data.datref > searchdatref)
					res = 1;
				else if (item->Data.datref < searchdatref)
					res = -1;
				}
			}
		if (res == 0) {
			*Found = TRUE;
			l = *x;
			}
		else {
			if (res > 0) /* Page.Item[].Key > *pSearchKey */
				r = *x;
			else {
				if (*Found) {
					return error("SearchPage() not ascending");
					}
				l = *x;
				}
			}
		}
	*x = l;
	item = buffer->Page.Item;
	for (i = 0; i < *x; i++) {
		childs = item[i].Childs;
		*cur += childs;
		*cur += 1;
		}
	return OK;
}

/*static INT *FKpFound;
static INT FKx;
static INT FKidx;*/
static void *FKpSearchKey;
static DATATYPE *FKpData;
static Buffer *FKBF;

INT BAYERTREE_OB::search(void *pSearchKey, 
	DATATYPE *pData, INT *idx, INT *f_found)
/* void pointer pSearchKey hier.
 * pSearchKey wird nicht benutzt, nur an (*cmp_func)() 
 * durchgeschleift. Insbesondere kann pSearchKey 
 * auf laengere Daten als KEYTYPE zeigen. 
 * Anwendung: ein evtl. langer Suchstring, 
 * der nicht in einen KEYTYPE passen wuerde. 
 * Es wird zusatzlich mit pData->datref 
 * gesucht, sofern pData != NIL.
 * idx enthaelt Nummer des gefundenen bzw des naechst 
 * kleineren Datensatzes. 
 * idx kann also auch -1 werden (nicht gefunden - 
 * vor dem 0-ten Datensatz) */
{
	Buffer BF;
	
	if (s_f_open_i() == FALSE) {
		return error("BT::search() not open");
		}
	FKpSearchKey = pSearchKey;
	FKpData = pData;
	FKBF = &BF;
	/*FGBF1 = &BF1;*/
	*idx = -1L;
	if (SearchBtree(this, s_Root_i(), idx, f_found) != OK) {
		return error("BT::search() error in SearchBtree()");
		}
	return OK;
}

static INT SearchBtree(BAYERTREE_OP p, 
	INT page, INT *idx, INT *f_found)
/* Sucht einen Schluessel in der Datenbank
 * --
 * KEYTYPE  FKpSearchKey - Zu suchender Schluessel
 * DATATYPE *FKpData     - Zum Schluessel gehoerende Daten
 * INT      *idx         - Enthaelt Nummer des gefundenen bzw des naechst 
 *                         kleineren Datensatzes
 * INT      *f_found     - Gibt an, ob Daten gefunden wurden */
{
	INT x, f_found1, idx1;
	
	if (page == 0)
		*f_found = FALSE;
	else {
		if (LoadPage(p, FKBF, page) != OK) {
			return error("SearchBtree() error in LoadPage()");
			}
		if (SearchPage(p, FKBF, FKpSearchKey, NIL, idx, 
			&x, f_found) != OK) {
			return error("SearchBtree() error in SearchPage()");
			}
		f_found1 = *f_found;
		idx1 = *idx;
		if (*f_found) {
			if (FKpData)
				*FKpData = FKBF->Page.Item[x].Data;
			}
		if (SearchBtree(p, 
			FKBF->Page.Item[x].Ref, 
			idx, f_found) != OK) {
			return error("SearchBtree() error in SearchBtree()");
			}
		if (!(*f_found) && f_found1) {
			if (FKpData) {
				if (LoadPage(p, FKBF, page) != OK) {
					return error("SearchBtree() error in LoadPage()");
					}
				*FKpData = FKBF->Page.Item[x].Data;
				}
			*f_found = TRUE;
			*idx = idx1;
			}
#if 0
		if (*f_found) {
			if (FKpData)
				*FKpData = FKBF->Page.Item[x].Data;
			}
		else {
			if (SearchBtree(p, 
				FKBF->Page.Item[x].Ref, 
				idx, f_found) != OK) {
				return error("SearchBtree() error in SearchBtree()");
				}
			}
#endif
		}
	return OK;
}

static INT page_i_th(BAYERTREE_OP p, 
	INT l, Buffer *buffer, 
	INT *cur, INT *i, INT *found)
{
	INT childs;
	INT page_len, j;
	
	page_len = buffer->Page.NumItems;
	for (j = 0; j <= page_len; j++) {
		childs = buffer->Page.Item[j].Childs;
		if (*cur + childs > l) {
			*i = j;
			*found = FALSE;
			return OK;
			}
		if (*cur + childs == l) {
			if (j == page_len) {
				return error("page_i_th() j == page_len");
				}
			/* gefunden: (in j + 1) */
			*i = j + 1;
			*found = TRUE;
			return OK;
			}
		/* naechster Zweig: */
		*cur += childs + 1;
		}
	return error("page_i_th() not found");
}

INT BAYERTREE_OB::ith(INT l, KEYTYPE *key, DATATYPE *data)
/* key muss auf einen ganzen KEYTYPE zeigen, 
 * da der gesamte keycarrier mittels struct-
 * zuweisung kopiert wird. */
{
	INT cur, page, ref;
	INT i, found, f_debug, ret;
	FILE *fp;
	Buffer buffer;
	BYTE str[256];

	if (FALSE) {
		f_debug = TRUE;
		fp = fopen("ith.log", "a");
		fprintf(fp, "BT::ith(%ld):\n", l);
		}
	else {
		f_debug = FALSE;
		}
	if (s_f_open_i() == FALSE) {
		return error("BT::ith() not open");
		}
	page = s_Root_i();
	cur = 0;
	while (TRUE) {
		if (LoadPage(this, &buffer, page) != OK) {
			Srfs("BT::ith", "LoadPage");
			return ERROR;
			}
		if (f_debug) {
			fprintf(fp, "page = %ld\n", page);
			fpage_print(this, &buffer, fp);
			}
		if (page_i_th(this, l, 
			&buffer, &cur, &i, &found) != OK) {
			Srff("BT::ith", "page_i_th");
			return ERROR;
			}
		if (f_debug) {
			fprintf(fp, "page_i_th: cur = %ld found = %ld i = %ld\n", 
				cur, found, i);
			fflush(fp);
			}
		if (found) {
			*key = buffer.Page.Item[i].Key;
			*data = buffer.Page.Item[i].Data;
			ret = OK;
			break;
			}
		ref = buffer.Page.Item[i].Ref;
		if (ref == 0) {
			Srfs("BT::ith", "ref == 0");
			sprintf(str, "l = %ld rootpage = %ld "
				"page = %ld i = %ld cur = %ld", 
				l, s_Root_i(), page, i, cur);
			Srfs("BT::ith", str);
			ret = ERROR;
			break;
			}
		page = ref;
		}
	if (f_debug) {
		fprintf(fp, "end of BT::ith(%ld):\n", l);
		fclose(fp);
		}
	return ret;
}

INT BAYERTREE_OB::ith0(INT l, KEYTYPE *key, DATATYPE *data)
{
	INT erg = ERROR;
	INT f_open;
	
	f_open = s_f_open_i();
	if (!f_open) {
		if (open() != OK) {
			Srfs("BT::ith0", "BT::open");
			return ERROR;
			}
		}
	if (ith(l, key, data) != OK) {
		Srfs("BT::ith0", "BT::ith");
		goto l_exit;
		}
	erg = OK;
l_exit:
	if (!f_open) {
		if (close() != OK) {
			Srfs("BT::ith0", "BT::close");
			return ERROR;
			}
		}
	return erg;
}

static KEYTYPE *IKpKey;
static DATATYPE *IKpData;
static Buffer *IKBF;
static INT IKFound;
static INT IKRisen;
static INT f_keyadded;

INT BAYERTREE_OB::insert_key(KEYTYPE *pKey, DATATYPE *pData)
/* Fuegt einen Schluessel in die Datenbank ein. Wenn
 * der Schluessel schon existiert, werden nur die Daten
 * ('Data') aktualisiert
 * --
 * KEYTYPE Key   - Einzufuegender Schluessel
 * DATATYPE Data - Die zum Schluessel gehoerenden Daten */
{
	INT RootSplit;
	ItemTyp RootItem;
	Buffer *NewRoot = NIL;
	Buffer *BF1 = NIL;
	Buffer *Root_BF;
	INT NewNeighbourChilds, size, new_page_num;
	INT erg = ERROR;
	
	if (s_f_open_i() == FALSE) {
		Srfs("BT::insert_key", "!open");
		return ERROR;
		}
	size = (INT)sizeof(Buffer);
	NewRoot = (Buffer *) my_malloc(size, "BT::insert_key");
	BF1 = (Buffer *) my_malloc(size, "BT::insert_key");
	if (NewRoot == NIL || BF1 == NIL) {
		Srfs("BT::insert_key", "no memory for Buffer");
		goto l_exit;
		}
	fill_char((BYTE *)NewRoot, size, 0);
	fill_char((BYTE *)BF1, size, 0);
	f_keyadded = FALSE;
	IKpKey = pKey;
	IKpData = pData;
	RootSplit = FALSE;
	IKBF = BF1;
	if (Update(this, s_Root_i(), 
		&RootSplit, &RootItem, 
		&NewNeighbourChilds) != OK) {
		Srff("BT::insert_key", "Update");
		goto l_exit;
		}
	if (RootSplit) {
		if (AllocateRec(this, 
			&new_page_num) != OK) {
			Srff("BT::insert_key", "AllocateRec");
			goto l_exit;
			}
		NewRoot->PageNum = new_page_num;
		/*printf("InsertKey()|RootSplit "
		"- NewRoot->PageNum = %ld\n", 
			NewRoot->PageNum);*/

		NewRoot->Page.NumItems = 1;
		NewRoot->Page.Item[0].Ref = s_Root_i();
		NewRoot->Page.Item[0].Childs = NewNeighbourChilds;
		NewRoot->Page.Item[1] = RootItem;
		if (SavePage(this, NewRoot) != OK) {
			Srff("BT::insert_key", "SavePage");
			goto l_exit;
			}
		Root_BF = RootBF + s_buf_idx_i();
		*Root_BF = *NewRoot;
		s_Root()->m_i(NewRoot->PageNum);
		if (WriteInfo(this) != OK) {
			Srff("BT::insert_key", "WriteInfo");
			goto l_exit;
			}
		}
	erg = OK;
l_exit:
	if (NewRoot) {
		my_free(NewRoot);
		NewRoot = NIL;
		}
	if (BF1) {
		my_free(BF1);
		BF1 = NIL;
		}
	IKBF = NIL;
	return erg;
}

INT BAYERTREE_OB::add(KEYTYPE *pKey, DATATYPE *pData)
{
	INT ret = ERROR;
	INT f_open;
	
	f_open = s_f_open_i();
	if (!f_open) {
		if (open() != OK) {
			Srfs("BT::add", "BT::open");
			return ERROR;
			}
		}
	if (insert_key(pKey, pData) != OK) {
		Srfs("BT::add", "BT::insert_key");
		goto l_exit;
		}
	ret = OK;
l_exit:
	if (!f_open) {
		if (close() != OK) {
			Srfs("BT::add", "BT::close");
			return ERROR;
			}
		}
	return ret;
}

static INT Update(BAYERTREE_OP p, 
	INT Node, INT *Rise, ItemTyp *RisenItem, 
	INT *RisenNeighbourChilds)
/* Einfuegen in den Zweig Node.
 * RisenNeighbourChilds nur gesetzt, wenn Rise TRUE ist. */
{
	INT x;
	INT idx, z;

	if (Node == 0) {
		/* Auf Blattebene angekommen wird von update()
		 * eingefuegt als wenn ein RisenItem eingefuegt 
		 * werden muesste. Deswegen wird hier diese
		 * Situation vorbereitet. */
		*Rise = TRUE;
		*RisenNeighbourChilds = 0;
		/* Nachbar darf ebenfalls keine Nachfolger haben */
		RisenItem->Key = *IKpKey;
		RisenItem->Data = *IKpData;
		RisenItem->Ref = 0;
		RisenItem->Childs = 0;
		f_keyadded = TRUE;
		return OK;
		}
	if (LoadPage(p, IKBF, Node) != OK) {
		Srff("Update", "LoadPage");
		return ERROR;
		}
	if (SearchPage(p, IKBF, (void *)IKpKey, IKpData, &idx, 
		&x, &IKFound) != OK) {
		Srff("Update", "SearchPage");
		return ERROR;
		}
	if (IKFound && !p->s_f_duplicatekeys_i()) {
		IKBF->Page.Item[x].Data = *IKpData;
		if (SavePage(p, IKBF) != OK) {
			Srff("Update", "SavePage");
			return ERROR;
			}
		/* keine doppelten Schluessel erlaubt */
		/* Rise bleibt unveraendert FALSE */
		return OK;
		}
	/* Einfuegen in den Zweig von x: */
	IKRisen = FALSE;
	if (Update(p, IKBF->Page.Item[x].Ref, 
		&IKRisen, RisenItem,
		RisenNeighbourChilds) != OK) {
		Srff("Update", "Update");
		printf("Node = %ld x = %ld "
		"IKBF->Page.Item[x].Ref = %ld IKRisen = %ld\n", 
			Node, x, (INT) IKBF->Page.Item[x].Ref, IKRisen);
		return ERROR;
		}
	if (LoadPage(p, IKBF, Node) != OK) {
		Srff("Update", "LoadPage");
		return ERROR;
		}
	/* Neuladen der Seite, da der Buffer in der Rekursion
	 * benutzt wird. */
	if (IKRisen) {
		/* RisenItem muss nach x eingefuegt werden: */
		IKBF->Page.Item[x].Childs = *RisenNeighbourChilds;
		/* Nach Seiten-Split hat der linke Nachbar weniger 
		 * Nachfolger */
		if (IKBF->Page.NumItems < BTREEMAXPAGESIZE) {
			/* Einfuegen auf dieser Seite: */
			IKBF->Page.NumItems++;
			for (z = IKBF->Page.NumItems - 1; 
				z >= x + 1; z--) {
				/* IKBF->Page.Item[z + 1] = 
					IKBF->Page.Item[z]; */
				bt_item_copy(
					&IKBF->Page.Item[z], 
					&IKBF->Page.Item[z + 1]);
				}
			/* IKBF->Page.Item[x + 1] = *RisenItem; */
			bt_item_copy(RisenItem, &IKBF->Page.Item[x + 1]);
			*Rise = FALSE;
			}
		else { /* Seite voll: splitting */
			/*printf("Update()| Seite voll: "
			"calling Split()| x = %ld\n", x);*/
			*RisenNeighbourChilds = 0; /* redundant */
			if (Split(p, IKBF, RisenItem, 
				x, RisenNeighbourChilds) != OK) {
				Srff("Update", "Split");
				return ERROR;
				}
			*Rise = TRUE;
			}
		}
	else { /* IKRisen == FALSE */
		if (f_keyadded)
			IKBF->Page.Item[x].Childs++;
		}
	if (SavePage(p, IKBF) != OK) {
		Srff("Update", "SavePage");
		printf("IKRisen = %ld x = %ld "
			"f_keyadded = %ld RisenNeighbourChilds = %ld\n", 
			IKRisen, x, f_keyadded, *RisenNeighbourChilds);
		printf("IKFound = %ld\n", IKFound);
		return ERROR;
		}
	return OK;
}

static void bt_item_copy(ItemTyp *a, ItemTyp *b)
{
	INT i, len;
	BYTE *pca = (BYTE *)a;
	BYTE *pcb = (BYTE *)b;
	
	len = sizeof(ItemTyp);
	for (i = 0; i < len; i++)
		pcb[i] = pca[i];
}

static INT Split(BAYERTREE_OP p, 
	Buffer *BF, ItemTyp *Item, 
	INT x, INT *RisenNeighbourChilds)
/* Fuegt Item in volle Seite BF nach position x ein.
 * Die uebervolle Seite wird zerlegt, die 2. Haelfte in 
 * eine neue Seite kopiert. Das mittlere Element wird 
 * angehoben und bekommt die neue Seite als Nachfolger.
 * Der alte Nachfolger des mittleren Elements wird in die 
 * neue Seite an 0ter Position eingehaengt.
 * RisenNeighbourChilds wird als Summe der Datensatze der 
 * verkleinerten Seite berechnet und muss von der auf-
 * rufenden Funktion links neben Item eingetragen werden.
 * Die neu generierte Seite wird abgespeichert; die ver-
 * kleinerte muss von der rufenden Funktion gesichert 
 * werden. */
{
	ItemTyp SplitItem;
	Buffer *SplitBF = NIL;
	INT sum1, sum2;
	INT z;
	INT size, new_page_num;
	INT erg = ERROR;
	
	size = (INT)sizeof(Buffer);
	SplitBF = (Buffer *) my_malloc(size, "BT::Split");
	if (SplitBF == NIL) {
		Srfs("Split", "no memory for Buffer");
		return ERROR;
		}
	fill_char((BYTE *)SplitBF, size, 0);
	if (AllocateRec(p, &new_page_num) != OK) {
		Srff("Split", "AllocateRec");
		return ERROR;
		}
	SplitBF->PageNum = new_page_num;
	/*printf("Split()|nach AllocateRec()|"
	"SplitBF->PageNum = %ld\n", SplitBF->PageNum);*/
	if (x < BTREEHALFPAGESIZE) {
		SplitItem = BF->Page.Item[BTREEHALFPAGESIZE];
		for (z = BTREEHALFPAGESIZE - 1; z >= x + 1; z--)
			BF->Page.Item[z + 1] = BF->Page.Item[z];
		BF->Page.Item[x + 1] = *Item;
		}
	else {
		if (x > BTREEHALFPAGESIZE) {
			SplitItem = BF->Page.Item[BTREEHALFPAGESIZE + 1];
			for (z = BTREEHALFPAGESIZE + 2L; z <= x; z++)
				BF->Page.Item[z - 1] = BF->Page.Item[z];
			BF->Page.Item[x] = *Item;
			}
		else
			SplitItem = *Item;
		}
	SplitBF->Page.Item[0].Ref = SplitItem.Ref;
	SplitBF->Page.Item[0].Childs = sum2 = SplitItem.Childs;
	sum1 = BF->Page.Item[0].Childs;
	for (z = 1L; z <= BTREEHALFPAGESIZE; z++) {
		SplitBF->Page.Item[z] = 
			BF->Page.Item[BTREEHALFPAGESIZE + z];
		sum2 += SplitBF->Page.Item[z].Childs;
		sum1 += BF->Page.Item[z].Childs;
		}
	sum1 += BTREEHALFPAGESIZE;
	sum2 += BTREEHALFPAGESIZE;
	BF->Page.NumItems = BTREEHALFPAGESIZE;
	SplitBF->Page.NumItems = BTREEHALFPAGESIZE;
	SplitItem.Ref = SplitBF->PageNum;
	SplitItem.Childs = sum2;
	*Item = SplitItem;
	if (SavePage(p, SplitBF) != OK) {
		Srff("Split", "SavePage");
		return ERROR;
		}
	*RisenNeighbourChilds = sum1;
	erg = OK;
	if (SplitBF) {
		my_free(SplitBF);
		SplitBF = NIL;
		}
	return erg;
}

static INT DKFound;
static INT DKidx, DKcur;
/*static Buffer *DKBF;*/
static INT f_keydeleted;

INT BAYERTREE_OB::del(INT idx)
/* Loescht einen Datensatz in der Datenbank
 * --
 * DelKey - Zu loeschender Schluessel */
{
	INT z, z2;
	INT Underflow;
	/* Buffer BF; */
	Buffer *Root_BF;
	
	if (s_f_open_i() == FALSE) {
		Srfs("BT::del", "!open");
		return ERROR;
		}
	DKidx = idx;
	DKcur = 0L;
	/* DKBF = &BF; */
	f_keydeleted = FALSE;
	if (Del(this, s_Root_i(), &Underflow) != OK) {
		Srff("BT::del", "Del");
		return ERROR;
		}
	Root_BF = RootBF + s_buf_idx_i();
	if (Underflow && Root_BF->Page.NumItems == 0) {
		z = s_Root_i();
		z2 = Root_BF->Page.Item[0].Ref;
		if (z2 == 0) { /* leere Datenbank */
			if (WriteInfo(this) != OK) { /* unnoetig ? */
				Srff("BT::del", "WriteInfo");
				return ERROR;
				}
			return OK;
			}
		if (LoadPage(this, Root_BF, z2) != OK) {
			Srff("BT::del", "LoadPage");
			return ERROR;
			}
		if (ReleaseRec(this, z) != OK) {
			Srff("BT::del", "ReleaseRec");
			return ERROR;
			}
		s_Root()->m_i(Root_BF->PageNum);
		if (WriteInfo(this) != OK) {
			Srff("BT::del", "WriteInfo");
			return ERROR;
			}
		}
	return OK;
}

static INT Del(BAYERTREE_OP p, 
	INT Node, INT *Underflow)
{
	INT x, y, z;
	INT size;
	Buffer *DKBF = NIL;
	INT erg = ERROR;

	if (Node == 0) {
		*Underflow = FALSE;
		return OK;
		}
	size = (INT)sizeof(Buffer);
	DKBF = (Buffer *) my_malloc(size, "BT::Del");
	if (DKBF == NIL) {
		Srfs("Del", "no memory for Buffer");
		return ERROR;
		}
	fill_char((BYTE *)DKBF, size, 0);
	if (LoadPage(p, DKBF, Node) != OK) {
		Srff("Del", "LoadPage");
		return ERROR;
		}
/*	SearchPage(p, DKBF, DKpDelKey, NIL, &idx, 
		&x, &DKFound);*/
	if (page_i_th(p, DKidx, DKBF, 
		&DKcur, &x, &DKFound) != OK) {
		Srff("Del", "page_i_th");
		return ERROR;
		}
	if (DKFound) {
		y = DKBF->Page.Item[x - 1].Ref;
		if (y == 0L) {
			DKBF->Page.NumItems--;
			*Underflow = 
				DKBF->Page.NumItems < BTREEHALFPAGESIZE;
			for (z = x; z <= DKBF->Page.NumItems; z++)
				DKBF->Page.Item[z] = DKBF->Page.Item[z + 1];
			if (SavePage(p, DKBF) != OK) {
				Srff("Del", "SavePage");
				return ERROR;
				}
			f_keydeleted = TRUE;
			}
		else {
			if (FindGreatest(p, y, 
				Underflow, DKBF, x) != OK) {
				Srff("Del", "FindGreatest");
				return ERROR;
				}
			if (f_keydeleted) {
				DKBF->Page.Item[x - 1].Childs--;
				if (SavePage(p, DKBF) != OK) {
					Srff("Del", "SavePage");
					return ERROR;
					}
				}
			if (*Underflow) {
				if (Compensate(p, Node, 
					y, x - 1, Underflow) != OK) {
					Srff("Del", "Compensate");
					return ERROR;
					}
				}
			}
		}
	else {
		y = DKBF->Page.Item[x].Ref;
		if (Del(p, y, Underflow) != OK) {
			Srff("Del", "Del");
			return ERROR;
			}
		if (f_keydeleted) {
			if (LoadPage(p, DKBF, Node) != OK) {
				/* kann entfallen, da DKBF jetzt lokal ! */
				Srff("Del", "LoadPage");
				return ERROR;
				}
			DKBF->Page.Item[x].Childs--;
			if (SavePage(p, DKBF) != OK) {
				Srff("Del", "SavePage");
				return ERROR;
				}
			}
		if (*Underflow) {
			if (Compensate(p, Node, y, x, Underflow) != OK) {
				Srff("Del", "Compensate");
				return ERROR;
				}
			}
		}
	erg = OK;
	if (DKBF) {
		my_free(DKBF);
		DKBF = NIL;
		}
	return erg;
}

static INT FindGreatest(BAYERTREE_OP p, 
	INT Node1, INT *Underflow, 
	Buffer *DKBF, INT x)
{
	INT Node2;
	INT NumBF;
	INT size;
	Buffer *buf = NIL; /* vormals FGBF1 */
	Buffer *Root_BF;
	INT erg = ERROR;

	size = (INT)sizeof(Buffer);
	buf = (Buffer *) my_malloc(size, "BT::FindGreatest");
	if (buf == NIL) {
		Srfs("FindGreatest", "no memory for Buffer");
		return ERROR;
		}
	fill_char((BYTE *)buf, size, 0);
	if (LoadPage(p, buf, Node1) != OK) {
		Srff("FindGreatest", "LoadPage");
		return ERROR;
		}
	NumBF = buf->Page.NumItems;
	Node2 = buf->Page.Item[NumBF].Ref;
	if (Node2 != 0) {
		if (FindGreatest(p, Node2, 
			Underflow, DKBF, x) != OK) {
			Srff("FindGreatest", "FindGreatest");
			return ERROR;
			}
		if (f_keydeleted) {
			if (LoadPage(p, buf, Node1) != OK) {
				/* kann entfallen, da Buffer lokal gemacht ! */
				Srff("FindGreatest", "LoadPage");
				return ERROR;
				}
			buf->Page.Item[NumBF].Childs--;
			if (SavePage(p, buf) != OK) {
				Srff("FindGreatest", "SavePage");
				return ERROR;
				}
			}
		if (*Underflow) {
			if (Compensate(p, Node1, Node2, 
				NumBF, Underflow) != OK) {
				Srff("FindGreatest", "Compensate");
				return ERROR;
				}
			}
		}
	else {
		DKBF->Page.Item[x].Key = buf->Page.Item[NumBF].Key;
		DKBF->Page.Item[x].Data = buf->Page.Item[NumBF].Data;
		f_keydeleted = TRUE;
		if (DKBF->PageNum == p->s_Root_i()) {
			Root_BF = RootBF + p->s_buf_idx_i();
			*Root_BF = *DKBF;
			}
		NumBF--;
		*Underflow = (NumBF < BTREEHALFPAGESIZE);
		buf->Page.NumItems = NumBF;
		if (SavePage(p, buf) != OK) {
			Srff("FindGreatest", "SavePage");
			return ERROR;
			}
		}
	erg = OK;
	if (buf) {
		my_free(buf);
		buf = NIL;
		}
	return erg;
}

static INT Compensate(BAYERTREE_OP p, INT Precedent, 
	INT Node, INT Path, INT *Underflow)
{
	INT Neighbour;
	INT NumBF2, NumBF3;
	INT x, z;
	INT sum;
	INT size;
	Buffer *BF1 = NIL;
	Buffer *BF2 = NIL;
	Buffer *BF3 = NIL;
	INT erg = ERROR;
	
	size = (INT)sizeof(Buffer);
	BF1 = (Buffer *) my_malloc(size, "BT::Compensate");
	BF2 = (Buffer *) my_malloc(size, "BT::Compensate");
	BF3 = (Buffer *) my_malloc(size, "BT::Compensate");
	if (BF1 == NIL || BF2 == NIL || BF3 == NIL) {
		Srfs("Compensate", "no memory for Buffer");
		return ERROR;
		}
	fill_char((BYTE *)BF1, size, 0);
	fill_char((BYTE *)BF2, size, 0);
	fill_char((BYTE *)BF3, size, 0);
	if (LoadPage(p, BF1, Node) != OK) {
		Srff("Compensate", "LoadPage");
		return ERROR;
		}
	if (LoadPage(p, BF3, Precedent) != OK) {
		Srff("Compensate", "LoadPage");
		return ERROR;
		}
	NumBF3 = BF3->Page.NumItems;
	if (Path < NumBF3) {
		/* Blatt nicht rechts aussen, dh. ein rechter
		 * Nachbar existiert. */
		Neighbour = BF3->Page.Item[Path + 1].Ref;
		if (LoadPage(p, BF2, Neighbour) != OK) {
			Srff("Compensate", "LoadPage");
			return ERROR;
			}
		NumBF2 = BF2->Page.NumItems;
		x = (NumBF2 + 1 - BTREEHALFPAGESIZE) / 2;
		BF1->Page.Item[BTREEHALFPAGESIZE].Key = 
			BF3->Page.Item[Path + 1].Key;
		BF1->Page.Item[BTREEHALFPAGESIZE].Data = 
			BF3->Page.Item[Path + 1].Data;
		BF1->Page.Item[BTREEHALFPAGESIZE].Ref = 
			BF2->Page.Item[0].Ref;
		BF1->Page.Item[BTREEHALFPAGESIZE].Childs = 
			sum = BF2->Page.Item[0].Childs;
		if (x > 0) {
			/*printf("Fall I x = %ld\n", x);*/
			for (z = 1; z <= x - 1; z++) {
				BF1->Page.Item[BTREEHALFPAGESIZE + z] = 
					BF2->Page.Item[z];
				sum += BF2->Page.Item[z].Childs;
				}
			sum += x;
			BF3->Page.Item[Path].Childs += sum;
			BF3->Page.Item[Path + 1].Childs -= sum;
			BF3->Page.Item[Path + 1].Key = BF2->Page.Item[x].Key;
			BF3->Page.Item[Path + 1].Data = BF2->Page.Item[x].Data;
			BF2->Page.Item[0].Ref = BF2->Page.Item[x].Ref;
			BF2->Page.Item[0].Childs = BF2->Page.Item[x].Childs;
			NumBF2 = NumBF2 - x;
			for (z = 1; z <= NumBF2; z++)
				BF2->Page.Item[z] = BF2->Page.Item[x + z];
			BF2->Page.NumItems = NumBF2;
			BF1->Page.NumItems = BTREEHALFPAGESIZE + x - 1;
			if (SavePage(p, BF1) != OK) {
				Srff("Compensate", "SavePage");
				return ERROR;
				}
			if (SavePage(p, BF2) != OK) {
				Srff("Compensate", "SavePage");
				return ERROR;
				}
			if (SavePage(p, BF3) != OK) {
				Srff("Compensate", "SavePage");
				return ERROR;
				}
			*Underflow = FALSE;
			}
		else {
			/*printf("Fall II x = %ld\n", x);*/
			BF3->Page.Item[Path].Childs += 
				BF3->Page.Item[Path + 1].Childs + 1;
			for (z = 1; z <= BTREEHALFPAGESIZE; z++)
				BF1->Page.Item[BTREEHALFPAGESIZE + z] = 
					BF2->Page.Item[z];
			for (z = Path + 1; z <= NumBF3 - 1; z++)
				BF3->Page.Item[z] = BF3->Page.Item[z + 1];
			BF1->Page.NumItems = BTREEMAXPAGESIZE;
			BF3->Page.NumItems = NumBF3 - 1L;
			*Underflow = (NumBF3 <= BTREEHALFPAGESIZE);
			if (SavePage(p, BF1) != OK) {
				Srff("Compensate", "SavePage");
				return ERROR;
				}
			if (SavePage(p, BF3) != OK) {
				Srff("Compensate", "SavePage");
				return ERROR;
				}
			if (ReleaseRec(p, Neighbour) != OK) {
				Srff("Compensate", "ReleaseRec");
				return ERROR;
				}
			}
		}
	else {
		/* Blatt rechts aussen; Nachbar links */
		Neighbour = BF3->Page.Item[Path - 1].Ref;
		if (LoadPage(p, BF2, Neighbour) != OK) {
			Srff("Compensate", "LoadPage");
			return ERROR;
			}
		NumBF2 = BF2->Page.NumItems;
		x = (NumBF2 + 1 - BTREEHALFPAGESIZE) / 2;
		if (x > 0) {
			/*printf("Fall III x = %ld\n", x);*/
			for (z = BTREEHALFPAGESIZE - 1L; z >= 1L; z--)
				BF1->Page.Item[z + x] = BF1->Page.Item[z];
			BF1->Page.Item[x].Key = BF3->Page.Item[Path].Key;
			BF1->Page.Item[x].Data = BF3->Page.Item[Path].Data;
			BF1->Page.Item[x].Ref = BF1->Page.Item[0].Ref;
			BF1->Page.Item[x].Childs = BF1->Page.Item[0].Childs;
			NumBF2 = NumBF2 - x;
			BF1->Page.Item[0].Ref = BF2->Page.Item[NumBF2 + 1].Ref;
			BF1->Page.Item[0].Childs = sum =
				BF2->Page.Item[NumBF2 + 1].Childs;
			for (z = x - 1; z >= 1; z--) {
				BF1->Page.Item[z] = BF2->Page.Item[NumBF2 + 1 + z];
				sum += BF2->Page.Item[NumBF2 + 1 + z].Childs;
				}
			sum += x;
			BF3->Page.Item[Path].Key = BF2->Page.Item[NumBF2 + 1].Key;
			BF3->Page.Item[Path].Data = BF2->Page.Item[NumBF2 + 1].Data;
			BF3->Page.Item[Path].Childs += sum;
			BF3->Page.Item[Path - 1].Childs -= sum;
			BF2->Page.NumItems = NumBF2;
			BF1->Page.NumItems = x + BTREEHALFPAGESIZE - 1L;
			if (SavePage(p, BF1) != OK) {
				Srff("Compensate", "SavePage");
				return ERROR;
				}
			if (SavePage(p, BF2) != OK) {
				Srff("Compensate", "SavePage");
				return ERROR;
				}
			if (SavePage(p, BF3) != OK) {
				Srff("Compensate", "SavePage");
				return ERROR;
				}
			*Underflow = FALSE;
			}
		else {
			/*printf("Fall IV x = %ld\n", x);*/
			BF2->Page.Item[NumBF2 + 1].Key = BF3->Page.Item[Path].Key;
			BF2->Page.Item[NumBF2 + 1].Data = BF3->Page.Item[Path].Data;
			BF2->Page.Item[NumBF2 + 1].Ref = BF1->Page.Item[0].Ref;
			BF2->Page.Item[NumBF2 + 1].Childs = BF1->Page.Item[0].Childs;
			for (z = 1; z <= BTREEHALFPAGESIZE - 1; z++)
				BF2->Page.Item[NumBF2 + 1 + z] = BF1->Page.Item[z];
			BF2->Page.NumItems = BTREEMAXPAGESIZE;
			BF3->Page.NumItems = NumBF3 - 1;
			BF3->Page.Item[Path - 1].Childs += 
				BF3->Page.Item[Path].Childs + 1;
			*Underflow = (NumBF3 <= BTREEHALFPAGESIZE);
			if (SavePage(p, BF2) != OK) {
				Srff("Compensate", "SavePage");
				return ERROR;
				}
			if (SavePage(p, BF3) != OK) {
				Srff("Compensate", "SavePage");
				return ERROR;
				}
			if (ReleaseRec(p, Node) != OK) {
				Srff("Compensate", "ReleaseRec");
				return ERROR;
				}
			}
		}
	erg = OK;
	if (BF1) {
		my_free(BF1);
		BF1 = NIL;
		}
	if (BF2) {
		my_free(BF2);
		BF2 = NIL;
		}
	if (BF3) {
		my_free(BF3);
		BF3 = NIL;
		}
	return erg;
}

INT BAYERTREE_OB::print()
/* Gibt die Datenbank auf den Bildschirm aus. */
{
	INT i, len1;
	KEYTYPE key;
	DATATYPE data;
	BYTE str[1024];
	
	if (s_f_open_i() == FALSE) {
		Srfs("BT::print", "!open");
		return ERROR;
		}
	if (len(&len1) != OK) {
		Srff("BT::print", "BT::len");
		return ERROR;
		}
	printf("Root = %ld\n", s_Root_i());
	printf("FreeRec = %ld\n", s_FreeRec_i());
	printf("AllocRec = %ld\n", s_AllocRec_i());
	printf("len = %ld\n", len1);
	for (i = 0; i < len1; i++) {
		if (ith(i, &key, &data) != OK) {
			Srff("BT::print", "bt_ith");
			return ERROR;
			}
		sprintf(str, "%ld: datref: %ld size: %ld ", 
		i, (INT) data.datref, (INT) data.data_size);
		bt_key_sprint(key.c, s_bt_key(), 
			s_nb_bt_key_i(), str);
		printf("%s\n", str);
		}
	return OK;
}

INT BAYERTREE_OB::print_pages()
{
	INT f_open;
	
	f_open = s_f_open_i();
	if (!f_open) {
		if (open() != OK) {
			Srff("BT::print_pages", "BT::open");
			return ERROR;
			}
		}
	
	print_page(s_Root_i());
	
	if (!f_open) {
		if (close() != OK) {
			Srff("BT::print_pages", "BT::close");
			return ERROR;
			}
		}
	return OK;
}

INT BAYERTREE_OB::print_page(INT x)
{
	INT y;
	Buffer BF;

	if (x == 0)
		return OK;
	LoadPage(this, &BF, x);
	printf("page %ld:\n", x);
	page_print(this, &BF);
	for (y = 0; y <= BF.Page.NumItems; y++) {
		print_page(BF.Page.Item[y].Ref);
		}
	return(TRUE);
}

static void page_print(BAYERTREE_OP p, Buffer *BF)
{
	INT i, len, childs, ref, datref, data_size;
	BYTE str[1024];
	ItemTyp *item = NIL;
	
	len = BF->Page.NumItems;
	printf("BF->Page.NumItems = %ld\n", len);
	for (i = 0; i <= len; i++) {
		item = &BF->Page.Item[i];
		childs = item->Childs;
		ref = item->Ref;
		printf("item %ld: Childs=%ld Ref=%ld", 
			i, childs, ref);
		if (i != 0) {
			datref = item->Data.datref;
			data_size = item->Data.data_size;
			sprintf(str, " (%ld/%ld): ", 
				datref, data_size);
			bt_key_sprint(item->Key.c, 
			p->s_bt_key(), p->s_nb_bt_key_i(), str);
			printf("%s\n", str);
			}
		else {
			printf("\n");
			}
		}
	printf("\n");
}

static void fpage_print(BAYERTREE_OP p, Buffer *BF, FILE *fp)
{
	INT i, len, childs, ref, datref, data_size;
	BYTE str[1024];
	ItemTyp *item = NIL;
	
	len = BF->Page.NumItems;
	fprintf(fp, "BF->Page.NumItems = %ld\n", len);
	for (i = 0; i <= len; i++) {
		item = &BF->Page.Item[i];
		childs = item->Childs;
		ref = item->Ref;
		fprintf(fp, "item %ld: Childs=%ld Ref=%ld", 
			i, childs, ref);
		if (i != 0) {
			datref = item->Data.datref;
			data_size = item->Data.data_size;
			sprintf(str, " (%ld/%ld): ", 
				datref, data_size);
			/* bt_key_sprint(item->Key.c, 
			p->s_bt_key(), p->s_nb_bt_key_i(), str); */
			fprintf(fp, "%s\n", str);
			}
		else {
			fprintf(fp, "\n");
			}
		}
	fprintf(fp, "\n");
}

#if 0
static void item_print(BAYERTREE_OP p, ItemTyp *item)
{
	BYTE str[1024];
	INT childs, ref, datref, data_size;

	childs = item->Childs;
	ref = item->Ref;
	datref = item->Data.datref;
	data_size = item->Data.data_size;
	sprintf(str, "Childs=%ld Ref=%ld (%ld/%ld): ", 
		childs, ref, datref, data_size);
	bt_key_sprint(item->Key.c, p->s_bt_key(), 
		p->s_nb_bt_key_i(), str);
	printf("%s\n", str);
}
#endif

#endif /* DB_TRUE */




