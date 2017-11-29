/* db.C */

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

#undef DEBUG_DB_PUT_FILE_SIZE
#undef DEBUG_DB_ADD_DATA_DB
// #undef DEBUG_DB_ADD_OP
#undef DEBUG_DB_GET_OP

/* 
 * DATABASE
 */

#define DB_SIZEOF_HEADER 16
#define DB_POS_FILESIZE 4

INT DATABASE_OB::field_name(INT i, INT j, BYTE *str)
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "nb_btree"; break;
	case 1: s = "btree"; break;
	case 2: s = "filename"; break;
	case 3: s = "sym_type"; break;
	case 4: s = "f_open"; break;
	case 5: s = "stream"; break;
	case 6: s = "file_size"; break;
	default:
		return error("DB::"
		"field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

INT DATABASE_OB::init(BYTE *file_name, INT sym_type)
{
	INT erg = OK;
	
	erg += m_il(7);
	c_obj_k(DATABASE);
	s_nb_btree()->m_i(0);
	erg += s_btree()->m_il(VECTOR_OVERSIZE);
	erg += s_filename()->init(file_name);
	s_sym_type()->m_i(sym_type);
	s_f_open()->m_i(0);
	s_file_size()->m_i(0);
	return erg;
}

INT DATABASE_OB::sprint(BYTE *s)
{
	BYTE str[512];
	
	sprintf(str, "DATABASE: %s "
		"file_size = %ld sym_type = %ld", 
		s_filename_s(), s_file_size_i(), 
		s_sym_type_i());
	if (strlen(s) + strlen(str) < 200)
		strcat(s, str);
	return OK;
}

INT DATABASE_OB::add_btree(BYTE *file_name, 
	INT f_duplicatekeys, INT btree_idx)
{
	BAYERTREE_OB btree;
	INT erg = OK;
	
	erg += btree.init(file_name, 
		f_duplicatekeys, btree_idx);
	erg += s_btree()->append_element(s_nb_btree(), &btree);
	return erg;
}

INT DATABASE_OB::put_file_size()
{
	FILE *f;
	INT4 l, l1;
	
	if (s_f_open_i() == FALSE) {
		Srfs("DB::put_file_size", "!open");
		return ERROR;
		}
	f = s_stream_i();
	if (f_seek_set(f, DB_POS_FILESIZE) != 0) {
		Srff("DB::put_file_size", "f_seek_set");
		return ERROR;
		}
	l = l1 = s_file_size_i();
	block_swap_bytes((SCHAR *)&l, sizeof(INT4), 1);
	if (fwrite((BYTE *)&l, sizeof(INT4), 1, f) != 1) {
		Srff("DB::put_file_size", "fwrite");
		return ERROR;
		}
#ifdef DEBUG_DB_PUT_FILE_SIZE
	printf("nach put_file_size %ld\n", (long)l1);
	printf("testing file_size:"); fflush(stdout);
	get_file_size();
	if (s_file_size_i() != l1) {
		printf("read file size of %ld; error\n", 
		(long)s_file_size_i()); fflush(stdout);
		return ERROR;
		}
	printf("OK\n"); fflush(stdout);
#endif
	return OK;
}

INT DATABASE_OB::get_file_size()
{
	FILE *f;
	INT4 l;
	
	if (s_f_open_i() == FALSE) {
		Srfs("DB::get_file_size", "!open");
		return ERROR;
		}
	f = s_stream_i();
	if (f_seek_set(f, DB_POS_FILESIZE) != 0) {
		Srff("DB::get_file_size", "f_seek_set");
		return ERROR;
		}
	if (fread((BYTE *)&l, sizeof(INT4), 1, f) != 1) {
		Srff("DB::get_file_size", "fread");
		return ERROR;
		}
	block_swap_bytes((SCHAR *)&l, sizeof(INT4), 1);
	s_file_size()->m_i(l);
	return OK;
}

INT DATABASE_OB::delete_files()
{
	BAYERTREE_OP bt;
	INT len, i;
	BYTE cmd[1024];

	sprintf(cmd, "rm %s", s_filename_s());
	call_system(cmd);
	len = s_nb_btree_i();
	for (i = 0; i < len; i++) {
		bt = s_btree_i(i);
		sprintf(cmd, "rm %s", bt->s_filename_s());
		call_system(cmd);
		}
	return OK;
}

INT DATABASE_OB::create()
{
	BAYERTREE_OP bt;
	INT len, i;
	BYTE buf[DB_SIZEOF_HEADER];
	FILE *f;
	
	if (s_f_open_i()) {
		Srfs("DB::create", "already open");
		return ERROR;
		}
	f = CREATE1(s_filename_s());
	if (f == NIL) {
		Srff("DB::create", "CREATE1");
		return ERROR;
		}
	s_stream()->m_i((INT)f);
	len = s_nb_btree_i();
	for (i = 0; i < len; i++) {
		bt = s_btree_i(i);
		if (bt->create() != OK) {
			Srff("DB::create", "BT::create");
			return ERROR;
			}
		}
	s_f_open()->m_i(TRUE);
	s_file_size()->m_i(DB_SIZEOF_HEADER);
	for (i = 0; i < DB_SIZEOF_HEADER; i++)
		buf[i] = 0;
	if (f_seek_set(f, 0) != 0) {
		Srff("DB::create", "f_seek_set");
		return ERROR;
		}
	if (fwrite(buf, DB_SIZEOF_HEADER, 1, f) != 1) {
		Srff("DB::create", "fwrite");
		return ERROR;
		}
	if (put_file_size() != OK) {
		Srff("DB::create", "DB::put_file_size");
		return ERROR;
		}
	return OK;
}

INT DATABASE_OB::open()
{
	BAYERTREE_OP bt;
	INT i, len;
	
	if (open_DB() != OK) {
		Srff("DB::open", "DB::open_DB");
		return ERROR;
		}
	len = s_nb_btree_i();
	for (i = 0; i < len; i++) {
		bt = s_btree_i(i);
		if (bt->open() != OK) {
			Srff("DB::open", "BT::open");
			return ERROR;
			}
		}
	return OK;
}

INT DATABASE_OB::open_DB()
{
	FILE *f;
	
	if (s_f_open_i()) {
		Srfs("DB::open_DB", "already open");
		return ERROR;
		}
	f = OPEN1(s_filename_s());
	if (f == NIL) {
		Srff("DB::open_DB", "OPEN1");
		return ERROR;
		}
	s_stream()->m_i((INT)f);
	s_f_open()->m_i(TRUE);
	if (get_file_size() != OK) {
		Srff("DB::open_DB", "DB::get_file_size");
		return ERROR;
		}
	return OK;
}

INT DATABASE_OB::close()
{
	BAYERTREE_OP bt;
	INT i, len;
	
	if (close_DB() != OK) {
		Srff("DB::close", "DB::close_DB");
		return ERROR;
		}
	len = s_nb_btree_i();
	for (i = 0; i < len; i++) {
		bt = s_btree_i(i);
		if (bt->close() != OK) {
			Srff("DB::close", "BT::close");
			return ERROR;
			}
		}
	return OK;
}

INT DATABASE_OB::close_DB()
{
	if (!s_f_open_i()) {
		Srfs("DB::close_DB", "!open");
		return ERROR;
		}
	if (put_file_size() != OK) {
		Srff("DB::close_DB", "DB::put_file_size");
		return ERROR;
		}
	CLOSE1(s_stream_i());
	s_stream()->m_i(0);
	s_file_size()->m_i(0);
	s_f_open()->m_i(FALSE);
	return OK;
}

#define MAGIC_SYNC 762873656L

static void user2total(INT user, INT *total, INT *pad);

INT DATABASE_OB::free_data_DB(INT datref, INT size)
{
	INT f_open, total;
	INT ret = ERROR;
	INT4 header[8];
	
	f_open = s_f_open_i();
	if (!f_open) {
		if (open_DB() != OK) {
			Srfs("DB::free_data_DB", 
			"DB::open_DB");
			return ERROR;
			}
		}
	
	if (f_seek_set(s_stream_i(), datref) != 0) {
		Srff("DB::free_data_DB", "f_seek_set");
		goto l_exit;
		}
	total = 8 * (INT)sizeof(INT4);
	if (fread(header, (INT)sizeof(BYTE), 
		total, s_stream_i()) != total) {
		Srff("DB::free_data_DB", "fread");
		return ERROR;
		}
	block_swap_bytes((SCHAR *)header, 
		sizeof(INT4), 8);
	if (header[0] != MAGIC_SYNC) {
		printf("DB::free_data_DB()|"
		"header: no MAGIC_SYNC\n");
		goto l_header2;
		}
	if (!header[1]) {
		printf("DB::free_data_DB()|"
		"header: data is not used\n");
		goto l_header2;
		}
	if (header[2] != size) {
		printf("DB::free_data_DB()|"
		"header: header[2] != size\n");
		goto l_header2;
		}
l_header2:
	if (header[4] != MAGIC_SYNC) {
		printf("DB::free_data_DB()|"
		"header2: no MAGIC_SYNC\n");
		goto l_exit;
		}
	if (!header[5]) {
		printf("DB::free_data_DB()|"
		"header2: data is not used\n");
		goto l_exit;
		}
	if (header[6] != size) {
		printf("DB::free_data_DB()|"
		"header2: header[6] != size\n");
		goto l_exit;
		}
	if (header[7] != header[3]) {
		printf("DB::free_data_DB()|"
		"header2: header[7] != header[3]\n");
		goto l_exit;
		}
	header[1] = FALSE;
	header[5] = FALSE;
	block_swap_bytes((SCHAR *)header, sizeof(INT4), 8);
	if (f_seek_set(s_stream_i(), datref) != 0) {
		Srff("DB::free_data_DB", "f_seek_set");
		goto l_exit;
		}
	if (fwrite(header, (INT)sizeof(BYTE), 
		total, s_stream_i()) != total) {
		Srff("DB::free_data_DB", "fwrite");
		return ERROR;
		}

	ret = OK;
l_exit:
	if (!f_open) {
		if (close_DB() != OK) {
			Srfs("DB::free_data_DB", 
			"close_DB");
			return ERROR;
			}
		}
	return ret;
}

INT DATABASE_OB::add_data_DB(void *d, INT size, INT *datref)
{
	INT f_open, total, pad;
	BYTE *data2 = NIL;
	BYTE *pc, *pc0;
	INT i;
	INT4 *pi;
	INT4 header[4];
	INT4 new_header[4];
		/* 0: SYNC
		 * 1: f_used
		 * 2: length user data
		 * 3: total length (header inclusive), 
		 *    a multiple of 16, 
		 *    one unused full 16 byte block garanteed.
		 */
	INT old_file_size;
	INT ret = ERROR;
	
#ifdef DEBUG_DB_ADD_DATA_DB
	printf("in DB::add_data_DB\n");
#endif
	f_open = s_f_open_i();
	if (!f_open) {
		if (open_DB() != OK) {
			Srfs("DB::add_data_DB", "DB::open_DB");
			return ERROR;
			}
		}
	user2total(size, &total, &pad);
	data2 = (BYTE *) my_malloc(total, "DB::add_data_DB");
	if (data2 == NIL) {
		Srff("DB::add_data_DB", "no memory for data2");
		goto l_exit;
		}
	header[0] = MAGIC_SYNC;
	header[1] = TRUE;
	header[2] = size;
	header[3] = total;
	block_swap_bytes((SCHAR *)header, sizeof(INT4), 4);
	pi = (INT4 *)data2;
	pi[0] = header[0];
	pi[1] = header[1];
	pi[2] = header[2];
	pi[3] = header[3];
	pi[4] = header[0];
	pi[5] = header[1];
	pi[6] = header[2];
	pi[7] = header[3];
	block_swap_bytes((SCHAR *)header, sizeof(INT4), 4);
		/* header zurueckswappen, 
		 * es folgt noch ein Test. */
	pc = (BYTE *)(pi + 8);
	pc0 = (BYTE *)d;
	/* printf("size = %ld pad = %ld total = %ld\n", size, pad, total); */
	for (i = 0; i < size; i++)
		pc[i] = pc0[i];
	for (i = 0; i < pad; i++)
		pc[size + i] = 0;
	old_file_size = s_file_size_i();
	if (f_seek_set(s_stream_i(), s_file_size_i()) != 0) {
		Srff("DB::add_data_DB", "f_seek_set");
		goto l_exit;
		}
	if (fwrite(data2, (INT)sizeof(BYTE), 
		total, s_stream_i()) != total) {
		Srfs("DB::add_data_DB", "error in fwrite() - data2");
		goto l_exit;
		}
	fflush(s_stream_i());
	*datref = s_file_size_i();
	s_file_size()->m_i(s_file_size_i() + total);
	if (put_file_size() != OK) {
		Srff("DB::add_data_DB", "DB::put_file_size");
		goto l_exit;
		}
	fflush(s_stream_i());
	
	if (f_seek_set(s_stream_i(), 
		old_file_size) != 0) {
		Srff("DB::add_data_DB", "f_seek_set");
		goto l_exit;
		}
	if (fread((BYTE *)new_header, 
		(INT)sizeof(INT4), 4, s_stream_i()) != 4) {
		Srfs("DB::add_data_DB", "error in fread");
		goto l_exit;
		}
	block_swap_bytes((SCHAR *)new_header, sizeof(INT4), 4);
	if (header[0] != new_header[0]) {
		printf("header[0] != new_header[0]\n");
		}
	if (header[1] != new_header[1]) {
		printf("header[1] != new_header[1]\n");
		}
	if (header[2] != new_header[2]) {
		printf("header[2] != new_header[2]\n");
		}
	if (header[3] != new_header[3]) {
		printf("header[3] != new_header[3]\n");
		}
	ret = OK;
l_exit:
	if (!f_open) {
		if (close_DB() != OK) {
			Srfs("DB::add_data_DB", "DB::close_DB");
			return ERROR;
			}
		}
	if (data2) {
		my_free(data2);
		data2 = NIL;
		}
#ifdef DEBUG_DB_ADD_DATA_DB
	printf("nach DB::add_data_DB\n");
#endif
	return ret;
}

static void user2total(INT user, INT *total, INT *pad)
{
	INT r, r1;
	
	r = user % 16;
	if (r != 0)
		r1 = 16 - r;
	else
		r1 = 0;
	*pad = r1 + 16;
	*total = sizeof(INT4) * 4 /* header */ + 
		sizeof(INT4) * 4 /* header 2 */ + 
		user + r1 + 16;
}

INT DATABASE_OB::ith(INT btree_idx, 
	INT l, KEYTYPE *key, DATATYPE *d)
/* key muss auf einen ganzen KEYTYPE zeigen, 
 * da der gesamte keycarrier mittels struct-
 * zuweisung kopiert wird. */
{
	BAYERTREE_OP bt = NIL;
	INT f_open;
	INT ret = ERROR;
	
	bt = s_btree_i(btree_idx);
	f_open = s_f_open_i();
	if (!f_open) {
		if (open_DB() != OK) {
			Srfs("DB::ith", "DB::open_DB");
			return ERROR;
			}
		if (bt->open() != OK) {
			Srff("DB::ith", "BT::open");
			return ERROR;
			}
		}
	
	if (bt->ith(l, key, d) != OK) {
		Srff("DB::ith", "BT::ith");
		goto l_exit;
		}

	ret = OK;
l_exit:
	if (!f_open) {
		if (close_DB() != OK) {
			Srfs("DB::ith", "DB::close_DB");
			return ERROR;
			}
		if (bt->close() != OK) {
			Srff("DB::ith", "BT::close");
			return ERROR;
			}
		}
	return ret;
}

INT DATABASE_OB::add_op(SYM_OP v)
{
	return add_op1(v, FALSE, FALSE);
}

INT DATABASE_OB::add_op1(SYM_OP v, INT f_v, INT f_vv)
{
	MEM_OP mem = (MEM_OP) callocobject("DATABASE::add_op() mem");
	BAYERTREE_OP bt;
	INT size, i, j, len, datref;
	BYTE *pc;
	KEYTYPE *key_type = NIL;
	DATATYPE data_type;
	INT ret = ERROR;
	
	if (f_v) {
		printf("in DB::add_op()\n");
		fflush(stdout);
		}
	key_type = (KEYTYPE *) my_malloc(sizeof(KEYTYPE), "DB::add_op");
	if (key_type == NIL) {
		Srfs("DB::add_op", "no memory for key_type");
		return ERROR;
		}
	if (v->s_obj_k() != s_sym_type_i()) {
		Srfs("DB::add_op", "v->s_obj_k() != s_sym_type_i()");
		return ERROR;
		}
	mem->init(0, NIL);
	if (f_v) {
		printf("in DB::add_op(): vor v->write_mem()\n");
		fflush(stdout);
		}
	v->write_mem(mem, 0 /* debug_depth */);
	if (f_v) {
		printf("in DB::add_op(): nach v->write_mem()\n");
		fflush(stdout);
		}
	mem->compress(f_v);
	size = mem->s_used_length_i();
	pc = mem->ob_self.ob_charpointer;
	if (f_v) {
		printf("in DB::add_op() nach v->write_mem()\n");
		fflush(stdout);
		}
	if (add_data_DB((void *)pc, size, &datref) != OK) {
		Srff("DB::add_op", "DB::add_data_DB");
		goto l_exit;
		}
	if (f_v) {
		printf("in DB::add_op() nach add_data_DB()\n");
		fflush(stdout);
		}
	data_type.datref = datref;
	data_type.data_size = size;
	len = s_nb_btree_i();
	for (i = 0; i < len; i++) {
		for (j = 0; j < BTREEMAXKEYLEN; j++) {
			key_type->c[j] = 0;
			}
		bt = s_btree_i(i);
		bt_key_fill_in(key_type->c, 
			(VECTOR_OP) v, bt->s_bt_key(), 
			bt->s_nb_bt_key_i());

		if (f_v) {
			printf("in DB::add_op() vor bt->add()\n");
			fflush(stdout);
			}
		if (bt->add(key_type, &data_type) != OK) {
			Srff("DB::add_op", "bt->add");
			goto l_exit;
			}
		if (f_v) {
			printf("in DB::add_op() nach bt->add()\n");
			fflush(stdout);
			}
		}
	
	ret = OK;
l_exit:
	freeall(mem);
	if (f_v) {
		printf("nach freeall(mem)\n");
		fflush(stdout);
		}
	if (key_type) {
		my_free(key_type);
		key_type = NIL;
		}

	if (f_v) {
		printf("nach DB::add_op()\n");
		fflush(stdout);
		}
	return ret;
}

INT DATABASE_OB::del_op(SYM_OP v, INT datref, INT size)
{
	BAYERTREE_OP bt;
	INT i, j, len;
	INT idx, f_found;
	KEYTYPE key_type;
	DATATYPE data_type;
	INT ret = ERROR;
	
	if (v->s_obj_k() != s_sym_type_i()) {
		Srfs("DB::del_op", 
		"v->s_obj_k() != s_sym_type_i()");
		return ERROR;
		}
	data_type.datref = datref;
	data_type.data_size = size;
	len = s_nb_btree_i();
	for (i = 0; i < len; i++) {
		for (j = 0; j < BTREEMAXKEYLEN; j++) {
			key_type.c[j] = 0;
			}
		bt = s_btree_i(i);
		bt_key_fill_in(key_type.c, 
		(VECTOR_OP) v, bt->s_bt_key(), 
		bt->s_nb_bt_key_i());
		if (bt->search(key_type.c, 
			&data_type, &idx, &f_found) != OK) {
			Srfs("DB::del_op", 
			"warning: btree entry not found");
			continue;
			}
		if (bt->del(idx) != OK) {
			Srff("DB::del_op", "bt->del");
			goto l_exit;
			}
		}
	if (free_data_DB(datref, size) != OK) {
		Srff("DB::del_op", "DB::free_data_DB");
		goto l_exit;
		}
	
	ret = OK;
l_exit:

	return ret;
}

INT DATABASE_OB::check_DB()
{
	FILE *f;
	INT4 header[4];
	INT4 header2[4];
	INT4 l;
	INT f_open, size, total, pad;
	INT cur_pos, entry_nr;
	INT ret = ERROR;
	
	f_open = s_f_open_i();
	if (!f_open) {
		if (open_DB() != OK) {
			Srfs("DB::check_DB", 
				"DB::open_DB");
			return ERROR;
			}
		}
	
	f = s_stream_i();
	
	if (f_seek_set(f, DB_POS_FILESIZE) != 0) {
		Srff("DB::check_DB", "f_seek_set");
		goto l_exit;
		}
	if (fread((BYTE *)&l, sizeof(INT4), 1, f) != 1) {
		Srff("DB::check_DB", "fread");
		goto l_exit;
		}
	block_swap_bytes((SCHAR *)&l, sizeof(INT4), 1);
	if (l != s_file_size_i()) {
		Srfs("DB::check_DB", "l != s_file_size_i()");
		goto l_exit;
		}
	cur_pos = DB_SIZEOF_HEADER;
	entry_nr = 0;
	while (TRUE) {
		if (cur_pos >= s_file_size_i())
			break;
		if (f_seek_set(f, cur_pos) != 0) {
			Srff("DB::check_DB", "f_seek_set");
			goto l_exit;
			}
		if (fread((BYTE *)header, (INT)sizeof(INT4), 4, f) != 4) {
			Srff("DB::check_DB", "fread");
			goto l_exit;
			}
		if (fread((BYTE *)header2, (INT)sizeof(INT4), 4, f) != 4) {
			Srff("DB::check_DB", "fread");
			goto l_exit;
			}
		block_swap_bytes((SCHAR *)header, sizeof(INT4), 4);
		block_swap_bytes((SCHAR *)header2, sizeof(INT4), 4);
		if (header[0] != MAGIC_SYNC) {
			printf("entry_nr %ld cur_pos = %ld: "
				"no MAGIC_SYNC in header\n", 
				entry_nr, cur_pos);
			goto l_exit;
			}
		if (header2[0] != MAGIC_SYNC) {
			printf("entry_nr %ld cur_pos = %ld: "
				"no MAGIC_SYNC in header2\n", 
				entry_nr, cur_pos);
			goto l_exit;
			}
		if (header[1] != header2[1]) {
			printf("entry_nr %ld cur_pos = %ld: "
			"header[1] = %ld != header2[1] = %ld\n", 
			entry_nr, cur_pos, 
			(INT) header[1], (INT) header2[1]);
			goto l_exit;
			}
		if (header[2] != header2[2]) {
			printf("entry_nr %ld cur_pos = %ld: "
			"header[2] = %ld != header2[2] = %ld\n", 
			entry_nr, cur_pos, 
			(INT) header[2], (INT) header2[2]);
			goto l_exit;
			}
		size = header[2];
		user2total(size, &total, &pad);
		if (header[3] != header2[3]) {
			printf("entry_nr %ld cur_pos = %ld: "
			"header[3] = %ld != header2[3] = %ld\n", 
			entry_nr, cur_pos, 
			(INT) header[3], (INT) header2[3]);
			goto l_exit;
			}
		if (header[3] != total) {
			printf("entry_nr %ld cur_pos = %ld: "
			"header[3] = %ld != total = %ld\n", 
			entry_nr, cur_pos, 
			(INT) header[3], total);
			goto l_exit;
			}
		printf("%ld at %ld: %ld size = %ld "
			"pad = %ld total = %ld next = %ld\n", 
			entry_nr, cur_pos, (INT) header[1], 
			size, pad, total, cur_pos + total);
		cur_pos += total;
		entry_nr++;
		}

	ret = OK;
l_exit:
	if (!f_open) {
		if (close_DB() != OK) {
			Srfs("DB::check_DB", 
			"DB::close_DB");
			return ERROR;
			}
		}
	return ret;
}

INT DATABASE_OB::get_op(
	DATATYPE *data_type, SYM_OP op)
{
	FILE *f;
	BYTE *pc, *pc1;
	MEM_OP mem = (MEM_OP) callocobject("DATABASE::get_op(): mem");
	BYTE *d = NIL;
	INT4 *header;
	INT4 *header2;
	INT f_open, size, total, pad, i;
	INT ret = ERROR;
	
	size = data_type->data_size;
	user2total(size, &total, &pad);
#ifdef DEBUG_DB_GET_OP
	printf("in DB::get_op()\n"); fflush(stdout);
#endif
	if (mem->alloc(size) != OK) {
		Srff("DB::get_op", "mem->alloc");
		return ERROR;
		}
	pc = mem->ob_self.ob_charpointer;
	if (pc == NIL) {
		Srfs("DB::get_op", 
		"mem->ob_self.ob_charpointer == NIL");
		return ERROR;
		}
	d = (BYTE *) my_malloc(total, "DB::get_op");
	if (d == NIL) {
		printf("DB::get_op()| no memory\n");
		return ERROR;
		}
	op->freeself();
	
	f_open = s_f_open_i();
	if (!f_open) {
		if (open_DB() != OK) {
			Srfs("DB::get_op", 
			"DB::open_DB");
			return ERROR;
			}
		}
	
	f = s_stream_i();
	
	if (f_seek_set(f, 0) != 0) {
		Srff("DB::get_op", "f_seek_set");
		return ERROR;
		}
	if (fread(d, (INT)sizeof(BYTE), 16, f) 
		!= 16) {
		Srff("DB::get_op", "fread - 16");
		return ERROR;
		}
	
	if (f_seek_set(f, data_type->datref) != 0) {
		Srff("DB::get_op", "f_seek_set");
		return ERROR;
		}
	if (fread(d, (INT)sizeof(BYTE), total, f) 
		!= total) {
		Srff("DB::get_op", "fread");
		return ERROR;
		}
	
	header = (INT4 *) d;
	header2 = header + 4;
	pc1 = d + 8 * (INT)sizeof(INT4);
	block_swap_bytes((SCHAR *)header, 
		sizeof(INT4), 4);
	block_swap_bytes((SCHAR *)header2, 
		sizeof(INT4), 4);
	if (header[0] != MAGIC_SYNC) {
		printf("DB::get_op()|header: "
		"no MAGIC_SYNC\n");
		goto l_header2;
		}
	if (!header[1]) {
		printf("DB::get_op()|header: "
		"data is not used\n");
		goto l_header2;
		}
	if (header[2] != size) {
		printf("DB::get_op()|header: "
		"header[2] != size\n");
		goto l_header2;
		}
	if (header[3] != total) {
		printf("DB::get_op()|header: "
		"header[3] != total\n");
		goto l_header2;
		}
l_header2:
	if (header2[0] != MAGIC_SYNC) {
		printf("DB::get_op()|header2: "
		"no MAGIC_SYNC\n");
		goto l_exit;
		}
	if (!header2[1]) {
		printf("DB::get_op()|header2: "
		"data is not used\n");
		goto l_exit;
		}
	if (header2[2] != size) {
		printf("DB::get_op()|header2: "
		"header2[2] != size\n");
		goto l_exit;
		}
	if (header2[3] != total) {
		printf("DB::get_op()|header2: "
		"header2[3] != total\n");
		goto l_exit;
		}
	for (i = 0; i < size; i++)
		pc[i] = pc1[i];
	mem->c_used_length(data_type->data_size);
	mem->decompress(FALSE /* f_verbose */);
	mem->c_cur_pointer(0);
	/* printf("allocated = %ld "
		"used = %ld cur = %ld\n", 
		mem->s_alloc_length_i(), 
		mem->s_used_length_i(), 
		mem->s_cur_pointer_i()); */
#ifdef DEBUG_DB_GET_OP
	printf("in DB::get_op() "
	"vor op->read_mem\n"); fflush(stdout);
#endif
	op->read_mem(mem, 0 /* debug_depth */);
	if (op->s_obj_k() != s_sym_type_i()) {
		Srfs("DB::get_op", 
		"op->s_obj_k() != s_sym_type_i()");
		return ERROR;
		}

	ret = OK;
l_exit:
	if (!f_open) {
		if (close_DB() != OK) {
			Srfs("DB::get_op", 
			"DB::close_DB");
			return ERROR;
			}
		}
	freeall(mem);
	if (d) {
		my_free(d);
		d = NIL;
		}
#ifdef DEBUG_DB_GET_OP
	printf("in DB::get_op() fertig\n");
	fflush(stdout);
#endif
	return ret;
}

INT DATABASE_OB::print(INT btree_idx)
/* Gibt die Datenbank auf den Bildschirm aus. */
{
	INT i, len1;
	INT f_open;
	KEYTYPE key;
	DATATYPE data;
	BAYERTREE_OP bt;
	SYM_OP tmp_op = callocobject("DATABASE::print");
	BYTE str[1024];
	
	f_open = s_f_open_i();
	bt = s_btree_i(btree_idx);
	if (!f_open) {
		if (open_DB() != OK) {
			Srfs("DB::print", "DB::open_DB");
			return ERROR;
			}
		if (bt->open() != OK) {
			Srff("DB::print", "BT::open");
			return ERROR;
			}
		}
	if (bt->len(&len1) != OK) {
		Srff("DB::print", "BT::len");
		return ERROR;
		}
	printf("Root = %ld\n", bt->s_Root_i());
	printf("FreeRec = %ld\n", bt->s_FreeRec_i());
	printf("AllocRec = %ld\n", bt->s_AllocRec_i());
	printf("len = %ld\n", len1);
	for (i = 0; i < len1; i++) {
		if (bt->ith(i, &key, &data) != OK) {
			Srff("DB::print", "BT::ith");
			return ERROR;
			}
		sprintf(str, "%ld: datref: %ld "
		"size: %ld ", i, 
		(INT) data.datref, (INT) data.data_size);
		bt_key_sprint(key.c, bt->s_bt_key(), 
		bt->s_nb_bt_key_i(), str);
		printf("%s\n", str);
		get_op(&data, tmp_op);
		tmp_op->println();
		}
	
	if (!f_open) {
		if (close_DB() != OK) {
			Srfs("DB::print", 
			"DB::close_DB");
			return ERROR;
			}
		if (bt->close() != OK) {
			Srff("DB::print", 
			"BT::close");
			return ERROR;
			}
		}
	freeall(tmp_op);
	return OK;
}

#undef DEBUG_RESTORE

INT DATABASE_OB::restore(FILE *fp_txt)
{
	BYTE restore_fname[1024];
	BYTE cmd[1000];

	sprintf(restore_fname, "%s.tmp", s_filename_s());
	sprintf(cmd, "mv %s %s", s_filename_s(), restore_fname);
	call_system(cmd);
	return restore_from_file(restore_fname, TRUE /* f_delete_afterwards */, fp_txt);
}

INT DATABASE_OB::restore_from_file(BYTE *restore_fname, 
	INT f_delete_afterwards, FILE *fp_txt)
{
	BYTE cmd[1000];
	FILE *f;
	INT4 fsize;
	INT4 header[4];
	INT4 header2[4];
	INT size, total, pad;
	INT cur_pos, entry_nr;
	INT ret = ERROR;

#ifdef DEBUG_RESTORE
	fprintf(fp_txt, "DATABASE::restore_from_file(%s)...\n", restore_fname);
	fflush(fp_txt);
#endif

	if (s_f_open_i())
		return error("DATABASE::restore() error: db already open");

	create();
#ifdef DEBUG_RESTORE
	fprintf(fp_txt, "DATABSE::restore(): database created");
	fflush(fp_txt);
#endif

	f = fopen(restore_fname, "r");
	if (f_seek_set(f, DB_POS_FILESIZE) != 0) {
		Srff("DB::restore", "f_seek_set");
		goto l_exit;
		}
	if (fread((BYTE *)&fsize, sizeof(INT4), 1, f) != 1) {
		Srff("DB::restore", "fread");
		goto l_exit;
		}
	block_swap_bytes((SCHAR *)&fsize, sizeof(INT4), 1);

	cur_pos = DB_SIZEOF_HEADER;
	entry_nr = 0;
	while (TRUE) {
		if (cur_pos >= fsize)
			break;
		if (f_seek_set(f, cur_pos) != 0) {
			Srff("DB::restore", "f_seek_set");
			goto l_exit;
			}
		if (fread((BYTE *)header, (INT)sizeof(INT4), 4, f) != 4) {
			Srff("DB::restore", "fread");
			goto l_exit;
			}
		if (fread((BYTE *)header2, (INT)sizeof(INT4), 4, f) != 4) {
			Srff("DB::restore", "fread");
			goto l_exit;
			}
		block_swap_bytes((SCHAR *)header, sizeof(INT4), 4);
		block_swap_bytes((SCHAR *)header2, sizeof(INT4), 4);
		if (header[0] != MAGIC_SYNC) {
			printf("entry_nr %ld cur_pos = %ld: "
				"no MAGIC_SYNC in header\n", 
				entry_nr, cur_pos);
			goto l_exit;
			}
		if (header2[0] != MAGIC_SYNC) {
			printf("entry_nr %ld cur_pos = %ld: "
				"no MAGIC_SYNC in header2\n", 
				entry_nr, cur_pos);
			goto l_exit;
			}
		if (header[1] != header2[1]) {
			printf("entry_nr %ld cur_pos = %ld: "
			"header[1] = %ld != header2[1] = %ld\n", 
			entry_nr, cur_pos, 
			(INT) header[1], (INT) header2[1]);
			goto l_exit;
			}
		if (header[2] != header2[2]) {
			printf("entry_nr %ld cur_pos = %ld: "
			"header[2] = %ld != header2[2] = %ld\n", 
			entry_nr, cur_pos, 
			(INT) header[2], (INT) header2[2]);
			goto l_exit;
			}
		size = header[2];
		user2total(size, &total, &pad);
		if (header[3] != header2[3]) {
			printf("entry_nr %ld cur_pos = %ld: "
			"header[3] = %ld != header2[3] = %ld\n", 
			entry_nr, cur_pos, 
			(INT) header[3], (INT) header2[3]);
			goto l_exit;
			}
		if (header[3] != total) {
			printf("entry_nr %ld cur_pos = %ld: "
			"header[3] = %ld != total = %ld\n", 
			entry_nr, cur_pos, 
			(INT) header[3], total);
			goto l_exit;
			}
		fprintf(fp_txt, "DATABASE::restore(): %ld at %ld: used=%ld size = %ld "
			"size = %ld total = %ld = %ld next = %ld\n", 
			entry_nr, cur_pos, (INT) header[1], (INT) header[2], 
			size, total, (INT) header[3], cur_pos + total);
		if (!header[1])
			goto loop;
		if (f_seek_set(f, cur_pos + 8 * 4) != 0) {
			Srff("DB::restore", "f_seek_set");
			goto l_exit;
			}
		{
			MEM_OB mem;
			SYM_OB ob;
			BYTE *pc;

			if (mem.alloc(size) != OK) {
				Srff("DB::restore", "mem->alloc");
				return ERROR;
				}
			pc = mem.ob_self.ob_charpointer;
			if (fread(pc, 1, size, f) != size) {
				Srff("DB::restore", "fread");
				goto l_exit;
				}
			mem.c_used_length(size);
			mem.decompress(FALSE /* f_verbose */);
			mem.c_cur_pointer(0);
			/* printf("allocated = %ld used = %ld cur = %ld\n", 
				mem.s_alloc_length_i(), 
				mem.s_used_length_i(), 
				mem.s_cur_pointer_i()); */
#ifdef DEBUG_RESTORE
			fprintf(fp_txt, "in DB::restore() read data into MEM_OB; calling op.read_mem\n");
			fflush(fp_txt);
#endif
			ob.read_mem(&mem, 0 /* debug_depth */);
			if (ob.s_obj_k() != s_sym_type_i()) {
				Srfs("DB::restore", "ob.s_obj_k() != s_sym_type_i()");
				return ERROR;
				}

#ifdef DEBUG_RESTORE
			fprintf(fp_txt, "in DB::restore() calling add_op()...");
			fflush(fp_txt);
#endif
			add_op(&ob);
#ifdef DEBUG_RESTORE
			fprintf(fp_txt, "finished !\n");
			fflush(fp_txt);
#endif
		}

loop:
		cur_pos += total;
		entry_nr++;
		}
	if (f_delete_afterwards) {
		sprintf(cmd, "rm %s", restore_fname);
		call_system(cmd);
		}

	ret = OK;
l_exit:
	close();
	return ret;
}


#endif /* DB_TRUE */




