/* geo_store.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1997
 */


#include <DISCRETA/discreta.h>
#include <DISCRETA/geo.h>
#include <DISCRETA/db.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/perm.h>
#include <DISCRETA/lb.h>
#include <DISCRETA/ladder.h>
#include <DISCRETA/divs.h>


/* 
 * GEO_BY_BASE_BLOCKS
 */

INT GEO_BY_BASE_BLOCKS_OB::field_name(INT i, INT j, BYTE *str)
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "hash"; break;
	case 1: s = "I"; break;
	case 2: s = "BB"; break; // base blocks
	default:
		return error("GEO_BY_BASE_BLOCKS::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

INT GEO_BY_BASE_BLOCKS_OB::init(VECTOR_OP BB, INT id)
{
	m_il(4);
	c_obj_k(GEO_BY_BASE_BLOCKS_KIND);
	s_hash()->init("");
	s_I()->m_ilih(0, 0);
	BB->copy(s_BB());
	s_id()->m_i(id);
	return OK;
}

INT GEO_BY_BASE_BLOCKS_OB::print()
{
	INT i, l;
	
	l = s_BB()->s_li();
	printf("%ld base blocks:", l);
	for (i = 0; i < l; i++) {
		s_BB_i(i)->println();
		}
	printf("\n");
	printf("hash: %s id = %ld\n", s_hash()->s_str(), s_id_i());
	return OK;
}

INT GEO_BY_BASE_BLOCKS_OB::sprint(BYTE *s)
{
	sprintf(s + strlen(s), "GEO_BY_BASE_BLOCKS_OB");
	return OK;
}

INT GEO_BY_BASE_BLOCKS_OB::fprint_incidences(FILE *fp)
{
	INT l, h, i, j, a;
	
	l = s_I()->s_li();
	h = s_I()->s_hi();
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			a = s_I()->s_iji(i, j);
			if (a == 0)
				continue;
			fprintf(fp, "%ld ", i * l + j);
			}
		}
	fprintf(fp, "\n");
	return OK;
}

INT GEO_BY_BASE_BLOCKS_OB::nb_incidences()
{
	INT l, h, i, j, a, n = 0;
	
	l = s_I()->s_li();
	h = s_I()->s_hi();
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			a = s_I()->s_iji(i, j);
			if (a == 0)
				continue;
			n++;
			}
		}
	return n;
}

INT GEO_BY_BASE_BLOCKS_OB::border(INT n)
{
	INT v, b, nb_V, nb_B;
	INT i, j, ii, jj, a;
	MATRIX_OB I;

	v = s_I()->s_hi();
	b = s_I()->s_li();
	nb_B = b / n;
	nb_V = v / n;
	I.m_ilih_n(b + nb_V, v + nb_B);
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			a = s_I()->s_iji(i, j);
			I.m_iji(nb_B + i, nb_V + j, a);
			}
		}
	for (jj = 0; jj < nb_B; jj++) {
		for (j = 0; j < n; j++) {
			I.m_iji(jj, nb_V + jj * n + j, 1);
			}
		}
	for (ii = 0; ii < nb_V; ii++) {
		for (i = 0; i < n; i++) {
			I.m_iji(nb_B + ii * n + i, ii, 1);
			}
		}
	I.swap(s_I());
	return OK;
}

INT GEO_BY_BASE_BLOCKS_OB::span_design(VECTOR_OP gen, INT f_v, INT f_vv)
{
	VECTOR_OB O, OO, tmp;
	VECTOR_OP B;
	MATRIX_OB I;
	PERMUTATION_OB p, q;
	INT v, b;
	INT i, l, ii, ll;
	
	if (f_v) {
		printf("GEO_BY_BASE_BLOCKS_OB::span_design()\n");
		fflush(stdout);
		}
	v = ((PERMUTATION_OP) gen->s_i(0))->s_li();
	if (f_v) {
		printf("a permutation group of degree %ld\n", v);
		fflush(stdout);
		}
	O.m_il(0);
	l = s_BB()->s_li();
	for (i = 0; i < l; i++) {
		B = s_BB_i(i);
		if (f_v) {
			printf("base block %ld: ", i + 1);
			ll = B->s_li();
			printf("[");
			for (ii = 0; ii < ll; ii++) {
				printf("%ld ", B->s_ii(ii) + 1);
				}
			printf("]\n");
			// B->println();
			fflush(stdout);
			}
		
		dc_calc_set_orbit(gen, B, &OO, FALSE /* f_v */, FALSE /* f_vv */);
		if (f_v) {
			printf("found an orbit of length %ld\n", OO.s_li());
			fflush(stdout);
			}
		O.append(&OO, &tmp);
		tmp.swap(&O);
		
		}
	if (f_v) {
		printf("altogether %ld blocks\n", O.s_li());
		fflush(stdout);
		}

	I.build_incidence_matrix(&O, v, FALSE /* f_entries_start_with_1 */);
	v = I.s_hi();
	b = I.s_li();
	I.swap(s_I());
	if (f_v) {
		printf("incidence matrix v=%ld b=%ld\n", v, b);
		if (f_vv) {
			I.print_incidences();
			fflush(stdout);
			}
		}
	return OK;
}
	
INT GEO_BY_BASE_BLOCKS_OB::canonical_form(INT f_tdo, INT f_tdo2, 
	INT f_V, VECTOR_OP V, 
	INT f_B, VECTOR_OP B, 
	INT f_v, INT f_vv)
{
	MATRIX_OP I;
	MATRIX_OB tdo_scheme;
	PERMUTATION_OB p, q;
	LABRA_OB aut1;
	VECTOR_OB aut_gen;
	VECTOR_OB row_decomp, col_decomp;
	INT v, b;
	
	if (f_v) {
		printf("GEO_BY_BASE_BLOCKS_OB::canonical_form()\n");
		fflush(stdout);
		if (f_tdo) {
			printf("using TDO\n");
			fflush(stdout);
			}
		if (f_tdo2) {
			printf("using TDO2\n");
			fflush(stdout);
			}
		}
	I = s_I();
	v = I->s_hi();
	b = I->s_li();
	
	if (f_V) {
		V->copy(&row_decomp);
		}
	else {
		row_decomp.m_il(1);
		row_decomp.m_ii(0, v);
		}
	if (f_B) {
		B->copy(&col_decomp);
		}
	else {
		col_decomp.m_il(1);
		col_decomp.m_ii(0, b);
		}
	
	if (f_tdo || f_tdo2) {
		s_I()->calc_TDO(f_tdo2 /* f_calc_second_tdo */, 
			FALSE /* f_ddp */, 
			FALSE /* f_ddb */, 
			&p, &q, 
			&row_decomp, &col_decomp, 
			NIL /* DDp */,
			NIL /* DDb */,
			&tdo_scheme, 
			f_vv, FALSE /* f_vv */);
		s_I()->apply_perms(&p, &q);
#if 0
		s_I()->print_char_TD(stdout, FALSE /* f_tex */, 
			&row_decomp, &col_decomp);
		print_tdo_scheme(&tdo_scheme, 
			&row_decomp, &col_decomp, FALSE /* f_tex */, stdout);
#endif
		}
	
	I->canonical_form(TRUE /* f_row_perms */, 
		&row_decomp, &col_decomp, 
		FALSE /* f_ddp */, FALSE /* f_ddb */,
		NULL /* DDp */, NULL /* DDb */, 
		&p, &q, 
		&aut1, TRUE /* f_aut_on_canonical_form */, &aut_gen, 
		TRUE /* f_apply_perms_to_canonical_form */, 
		f_v /* f_v */, FALSE /* f_vv */);
	if (f_v) {
		printf("canonical form\n");
		if (f_vv) {
			I->print_incidences();
			fflush(stdout);
			}
		}
	
	return OK;
}

INT GEO_BY_BASE_BLOCKS_OB::calc_hash(INT f_v, INT f_vv)
{
	INT key_len, al_len;
	BYTE *alphabet = NIL;
	BYTE *inc = NIL;
	BYTE *h = NIL;
	static INT primes[] = { 
		2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 
		31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 
		73, 79, 83, 89, 97 } ; // the first 25 primes 
	INT i0, i, j, k, v, b, nb_inc, pr, x, y;
	BYTE c;
	
	key_len = BTREEMAXKEYLEN;
	v = s_I()->s_hi();
	b = s_I()->s_li();
	nb_inc = 0;
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			if (s_I()->s_iji(i, j) != 0)
				nb_inc++;
			}
		}
			
	al_len = MAXIMUM(256, b);
	alphabet = (BYTE *) my_malloc(sizeof(BYTE) * (al_len + 1), "alphabet for GEO_BY_BASE_BLOCKS");
	inc = (BYTE *) my_malloc(sizeof(BYTE) * (nb_inc + 1), "alphabet for GEO_BY_BASE_BLOCKS");
	h = (BYTE *) my_malloc(sizeof(BYTE) * (key_len + 1), "alphabet for GEO_BY_BASE_BLOCKS");
	i0 = 0;
	k = 0;
	while (k < al_len) {
		for (i = 0; i < 10; i++) {
			alphabet[k] = '0' + i;
			k++;
			if (k >= al_len)
				break;
			}
		for (i = 0; i < 26; i++) {
			alphabet[k] = 'a' + i;
			k++;
			if (k >= al_len)
				break;
			}
		for (i = 0; i < 26; i++) {
			alphabet[k] = 'A' + i;
			k++;
			if (k >= al_len)
				break;
			}
		} // while
	alphabet[al_len] = 0;
	if (f_vv) {
		printf("GEO_BY_BASE_BLOCKS: alphabet: %s\n", alphabet);
		}
	
	k = 0;
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			if (s_I()->s_iji(i, j) == 0)
				continue;
			c = alphabet[j];
			inc[k] = c;
			k++;
			}
		}
	inc[nb_inc] = 0;
	if (f_vv) {
		printf("GEO_BY_BASE_BLOCKS: incidences: %s\n", inc);
		}
	
	j = 0;
	for (k = 0; k < key_len; k++) {
		pr = primes[k % 25];
		x = 0;
		for (i = 0; i < nb_inc; i++) {
			y = (INT)inc[j];
			x = (x + y) % 256;
			j += pr;
			if (j >= nb_inc)
				j = 0;
			}
		h[k] = alphabet[x];
		// printf("k=%ld pr = %ld x=%ld h[k]=%c\n", k, pr, x, h[k]);
		}
	h[key_len - 1] = 0;
	if (f_v) {
		printf("GEO_BY_BASE_BLOCKS: hash: %s\n", h);
		}
	s_hash()->init(h);
	
#if 0
	if (f_v) {
		printf("hash key:\n");
		s_hash()->println();
		fflush(stdout);
		}
#endif
	my_free(alphabet);
	my_free(inc);
	my_free(h);
	return OK;
}

#if TEXDOCU
INT GEO_BY_BASE_BLOCKS_OB::Equal(GEO_BY_BASE_BLOCKS_OP G1, INT f_v)
#endif
{
	INT i, j, h, l, a, b;
	MATRIX_OP I, I1;

	h = s_I()->s_hi();
	l = s_I()->s_li();
	if (G1->s_I()->s_hi() != h)
		return error("GEO_BY_BASE_BLOCKS_OB::Equal(): h differs");
	if (G1->s_I()->s_li() != l)
		return error("GEO_BY_BASE_BLOCKS_OB::Equal(): l differs");
	I = s_I();
	I1 = G1->s_I();
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			a = I->s_iji(i, j);
			b = I1->s_iji(i, j);
			if (a != b)
				return FALSE;
			}
		}
	return TRUE;
}

#if TEXDOCU
INT i_geo_by_base_blocks_db(DATABASE_OP *db, BYTE *db_prefix)
#endif
{
	DATABASE_OP db1;
	BAYERTREE_OP bt;
	BYTE fn_db[1024];
	BYTE fn_idx[1024];

	sprintf(fn_db, "%sgeo_by_base_blocks.db", db_prefix);
	sprintf(fn_idx, "%sgeo_by_base_blocks.idx", db_prefix);

	db1 = (DATABASE_OP) callocobject("i_geo_by_base_blocks_db");
	db1->init(fn_db, GEO_BY_BASE_BLOCKS_KIND);
	db1->add_btree(fn_idx, TRUE /* f_duplicatekeys */, 0 /* btree_idx */ );
	bt = db1->s_btree_i(0);
	bt->add_key_STRING(BTREEMAXKEYLEN, 0, 0);

	*db = db1;
	return OK;
}

#if TEXDOCU
void e_geo_by_base_blocks_db(DATABASE_OP db)
#endif
{
	freeall(db);
}

#if TEXDOCU
INT db_geo_by_base_blocks_create_and_close(BYTE *db_prefix)
#else
creates the database, closes it (an empty database). 
#endif
{
	DATABASE_OP db;
	
	if (i_geo_by_base_blocks_db(&db, db_prefix) != OK) {
		return error("db_geo_by_base_blocks_create_and_close() error in i_geo_by_base_blocks_db()");
		}
	if (db->create() != OK) {
		return error("db_geo_by_base_blocks_create_and_close() error in DB::create()");
		}
	if (db->close() != OK) {
		return error("db_geo_by_base_blocks_create_and_close() error in DB::close()");
		}
	e_geo_by_base_blocks_db(db);
	return OK;
}

INT db_geo_by_base_blocks_find_isomorphic_via_hash(DATABASE_OP db, 
	FILE *fp_txt, GEO_BY_BASE_BLOCKS_OP G, INT f_verbose, 
	INT *idx_btree_hash, INT *f_found)
{
	BAYERTREE_OP btree;
	INT btree_idx_hash = 0;
	KEYTYPE key_type, key_type1;
	DATATYPE data_type, data_type1;
	INT idx, f_found1, idx1, res;
	GEO_BY_BASE_BLOCKS_OB G1;

	if (f_verbose) {
		fprintf(fp_txt, "db_geo_by_base_blocks_find_isomorphic_via_hash(");
		fflush(fp_txt);
		}

	btree = db->s_btree_i(btree_idx_hash);
	if (btree->open() != OK) {
		return f_error(fp_txt, "db_geo_by_base_blocks_find_isomorphic_via_hash() "
			"error in BT::open");
		}

	bt_key_fill_in(key_type.c, G, 
		btree->s_bt_key(), btree->s_nb_bt_key_i());
	if (btree->search(&key_type, &data_type, 
		&idx, &f_found1) != OK) {
		return f_error(fp_txt, "db_geo_by_base_blocks_find_isomorphic_via_hash() "
			"error in btree->search()");
		}
	if (f_found1) {
		if (f_verbose) {
			fprintf(fp_txt, "hash found at %ld ", idx);
			fflush(fp_txt);
			}
		db->open_DB();
		idx1 = idx;
		*f_found = FALSE;
		while (idx1 >= 0) {
			if (btree->ith(idx1, &key_type1, &data_type1) != OK) {
				return error("db_geo_by_base_blocks_find_isomorphic_via_hash() "
					"error in btree->ith");
				}
			bt_key_cmp(key_type.c, key_type1.c, 
				btree->s_bt_key(), btree->s_nb_bt_key_i(), &res);
			if (res != 0 && idx1 == idx) {
				db->close();
				return f_error(fp_txt, "db_geo_by_base_blocks_find_"
					"isomorphic_via_hash(): res != 0 && idx1 == idx");
				}
			if (res != 0) { /* not in database; add it */
				*f_found = FALSE;
				break;
				}
			if (db->get_op(&data_type1, &G1) != OK) {
				return f_error(fp_txt, "db_geo_by_base_blocks_find_"
					"isomorphic_via_hash() error in db->get_op(&G1)");
				}
			res = G->Equal(&G1, FALSE /* f_v */);
			if (res) {
				if (f_verbose) {
					fprintf(fp_txt, "isomorphic, skipped");
					fflush(fp_txt);
					}
				*f_found = TRUE;
				*idx_btree_hash = idx1;
				break;
				}
			if (f_verbose) {
				fprintf(fp_txt, ".");
				fflush(fp_txt);
				}
			idx1--;
			}
		/* now f_found is set in any case */
		db->close_DB();
		if (f_verbose) {
			fprintf(fp_txt, ")\n");
			fflush(fp_txt);
			}
		}
	else {
		if (f_verbose) {
			fprintf(fp_txt, "hash value not found)\n");
			fflush(fp_txt);
			}
		*f_found = FALSE;
		}

	btree->close();
	return OK;
}

INT db_geo_by_base_blocks_export(BYTE *db_prefix, BYTE *geo_fname)
{
	FILE *fp;
	INT i, l, j, ll, lll;
	GEO_BY_BASE_BLOCKS_OB geo;
	DATABASE_OP db;
	BAYERTREE_OP bt;
	KEYTYPE key_type;
	DATATYPE data_type;
	MATRIX_OP I;
	VECTOR_OP bb;
	INT v, b, ii, jj, a, nb_inc;

	fp = fopen(geo_fname, "w");
	
	if (i_geo_by_base_blocks_db(&db, db_prefix) != OK) {
		return error("db_geo_by_base_blocks_export() error in i_geo_by_base_blocks_db()");
		}
	bt = db->s_btree_i(0);
	if (db->open() != OK) {
		return error("db_geo_by_base_blocks_export() error in db->open()");
		}
	bt->len(&l);
	printf("db_geo_by_base_blocks_export() database has length %ld\n", l);
	for (i = 0; i < l; i++) {
		if (bt->ith(i, &key_type, &data_type) != OK) {
			return error("db_geo_by_base_blocks_export() "
				"error in btree->ith");
			}
		if (db->get_op(&data_type, &geo) != OK) {
			return error("db_geo_by_base_blocks_export()"
				"error in db->get_op(&geo)");
			}
		ll = geo.s_BB()->s_li();
		// printf("base blocks: ");
		// for (j = 0; j < ll; j++) {
		//	geo.s_BB_i(j)->println();
		//	}
		for (j = 0; j < ll; j++) {
			bb = geo.s_BB_i(j);
			lll = bb->s_li();
			for (jj = 0; jj < lll; jj++) {
				a = bb->s_ii(jj) + 1;
				printf("%ld ", a);
				}
			printf("\n");
			}
		printf("\n");
		
		I = geo.s_I();
		v = I->s_hi();
		b = I->s_li();
		if (i == 0) {
			nb_inc = 0;
			for (ii = 0; ii < v; ii++) {
				for (jj = 0; jj < b; jj++) {
					if (I->s_iji(ii, jj)) {
						nb_inc++;
						}
					}
				}
			fprintf(fp, "%ld %ld %ld\n", v, b, nb_inc);
			}
		for (ii = 0; ii < v; ii++) {
			for (jj = 0; jj < b; jj++) {
				if (I->s_iji(ii, jj)) {
					a = ii * b + jj;
					fprintf(fp, "%ld ", a);
					}
				}
			}
		fprintf(fp, "\n");
		}
	
	fprintf(fp, "-1 %ld geometries\n", l);
	
	if (db->close() != OK) {
		return error("db_geo_by_base_blocks_export() error in db->close()");
		}
	
	fclose(fp);
	return OK;
}

INT db_geo_by_base_blocks_export_id(BYTE *db_prefix, BYTE *id_fname)
{
	FILE *fp;
	INT i, l;
	GEO_BY_BASE_BLOCKS_OB geo;
	DATABASE_OP db;
	BAYERTREE_OP bt;
	KEYTYPE key_type;
	DATATYPE data_type;
	INT id;

	fp = fopen(id_fname, "w");
	
	if (i_geo_by_base_blocks_db(&db, db_prefix) != OK) {
		return error("db_geo_by_base_blocks_export_id() error in i_geo_by_base_blocks_db()");
		}
	bt = db->s_btree_i(0);
	if (db->open() != OK) {
		return error("db_geo_by_base_blocks_export_id() error in db->open()");
		}
	bt->len(&l);
	printf("db_geo_by_base_blocks_export_id() database has length %ld\n", l);
	fprintf(fp, "%ld\n", l);
	for (i = 0; i < l; i++) {
		if (bt->ith(i, &key_type, &data_type) != OK) {
			return error("db_geo_by_base_blocks_export() "
				"error in btree->ith");
			}
		if (db->get_op(&data_type, &geo) != OK) {
			return error("db_geo_by_base_blocks_export()"
				"error in db->get_op(&geo)");
			}
		id = geo.s_id_i();
		fprintf(fp, "%ld\n", id);
		}
	
	fprintf(fp, "-1 %ld geometries\n", l);
	
	if (db->close() != OK) {
		return error("db_geo_by_base_blocks_export() error in db->close()");
		}
	
	fclose(fp);
	return OK;
}

#define BUFSIZE 50000

INT db_geo_by_base_blocks_add(BYTE *db_prefix, INT f_0, INT f_first, 
	BYTE *group_fname, BYTE *base_block_fname, 
	INT f_inc_file, BYTE *inc_file_name, 
	INT f_tdo, INT f_tdo2, 
	INT f_tda, INT f_range, INT first, INT len, 
	INT f_test2design, 
	INT f_V, VECTOR_OP decomp_V, 
	INT f_B, VECTOR_OP decomp_B, 
	INT f_do_not_add, 
	INT f_v, INT f_vv)
{
	BYTE buf[BUFSIZE];
	FILE *in_fp;
	FILE *OUT_fp;
	FILE *TDA_fp;
	BYTE OUT_fname[1024];
	BYTE TDA_fname[1024];
	INT a, i, j;
	VECTOR_OB group_generators;
	LABRA_OB aut, aut1;
	VECTOR_OB BB, B;
	INT v, nb_base_blocks, nb_bb, bb_len = 0;
	INT b, nb_inc, x, y;
	GEO_BY_BASE_BLOCKS_OB geo;
	DATABASE_OP db;
	BAYERTREE_OP bt;
	INT a1, a2, a3, a4;
	INT idx, f_found, nb_input, nb_geo;
	INT design_k, design_lambda, k, l, f_first_design = TRUE;

	if (f_first) {
		db_geo_by_base_blocks_create_and_close(db_prefix);
		}
	
	if (f_inc_file) {
		printf("reading .inc file %s\n", inc_file_name);
		in_fp = fopen(inc_file_name, "r");
		if (in_fp == NIL) {
			printf("inc_file_name = %s\n", inc_file_name);
			return error("can't open inc-file for reading");
			}
		strcpy(OUT_fname, inc_file_name);
		strcpy(TDA_fname, inc_file_name);
		}
	else {
		read_file_of_generators(&group_generators, group_fname);
		printf("read %ld group generators\n", group_generators.s_li());
		aut.init_quick(&group_generators);
		aut.print_group_order();
		in_fp = fopen(base_block_fname, "r");
		if (in_fp == NIL) {
			printf("base_block_fname = %s\n", base_block_fname);
			return error("can't open base-block-file for reading");
			}
		strcpy(OUT_fname, base_block_fname);
		strcpy(TDA_fname, base_block_fname);
		}
	
	
	strcat(OUT_fname, ".reduced");
	if (f_first) {
		OUT_fp = fopen(OUT_fname, "w");
		}
	else {
		OUT_fp = fopen(OUT_fname, "a");
		}

	if (f_tda) {
		strcat(TDA_fname, ".tda");
		if (f_first) {
			TDA_fp = fopen(TDA_fname, "w");
			}
		else {
			TDA_fp = fopen(TDA_fname, "a");
			}
		}
	
	
	if (f_inc_file) {
		fscanf(in_fp, "%ld %ld %ld", &v, &b, &nb_inc);
		printf("inc_file: v=%ld; b=%ld nb_inc=%ld\n", 
			v, b, nb_inc);
		fflush(stdout);
		fprintf(OUT_fp, "%ld %ld %ld\n", v, b, nb_inc);
		fflush(OUT_fp);
		}
	else {
		fscanf(in_fp, "%ld %ld %ld", &v, &nb_base_blocks, &bb_len);
		printf("base-block-file: v=%ld; "
			"reading %ld base_blocks of size %ld each\n", 
			v, nb_base_blocks, bb_len);
		fprintf(OUT_fp, "%ld\n", v);
		fflush(OUT_fp);
		}
	
	if (i_geo_by_base_blocks_db(&db, db_prefix) != OK) {
		return error("db_geo_by_base_blocks_add() error in i_geo_by_base_blocks_db()");
		}
	bt = db->s_btree_i(0);
	if (db->open() != OK) {
		return error("db_geo_by_base_blocks_add() error in db->open()");
		}
	
	// BB.m_il(bb_len);
	nb_input = 0;
	nb_geo = 0;
	while (TRUE) {
		fscanf(in_fp, "%ld", &a);
		if (a == -1)
			break;
		nb_input++;
		if (f_range) {
			if (nb_input < first || nb_input >= first + len) {
				// skip the rest of the input of this geometry:
				if (f_inc_file) {
					for (i = 0; i < nb_inc; i++) {
						fscanf(in_fp, "%ld", &a);
						}
					fgets(buf, BUFSIZE, in_fp);
					}
				else {
					if (nb_base_blocks == -1) {
						nb_bb = a;
						fscanf(in_fp, "%ld", &a);
						}
					else {
						nb_bb = nb_base_blocks;
						}
					for (i = 0; i < nb_bb; i++) {
						for (j = 0; j < bb_len; j++) {
							if (i == 0 && j == 0) {
								;
								}
							else {
								fscanf(in_fp, "%ld", &a);
								}
							}
						}
					}
				continue;
				}
			}
		printf("\n\nreading input geometry no %ld\n", nb_input);
		if (f_inc_file) {
			BB.m_il(0);
			geo.init(&BB, nb_input - 1);
			geo.s_I()->m_ilih_n(b, v);
			for (i = 0; i < nb_inc; i++) {
				if (i != 0) {
					fscanf(in_fp, "%ld", &a);
					}
				if (a == -1)
					return error("db_geo_by_base_blocks_add() a == -1");
				y = a % b;
				x = a / b;
				geo.s_I()->m_iji(x, y, 1);
				}
			fgets(buf, BUFSIZE, in_fp);
			}
		else {
			if (nb_base_blocks == -1) {
				nb_bb = a;
				fscanf(in_fp, "%ld", &a);
				}
			else {
				nb_bb = nb_base_blocks;
				}
			printf("nb_bb = %ld\n", nb_bb);
			fflush(stdout);
			BB.m_il(nb_bb);
			for (i = 0; i < nb_bb; i++) {
				B.m_il(bb_len);
				for (j = 0; j < bb_len; j++) {
					if (i == 0 && j == 0) {
						;
						}
					else {
						fscanf(in_fp, "%ld", &a);
						if (a == -1)
							return error("db_geo_by_base_blocks_add() a == -1");
						}
					if (f_0) {
						a++;
						}
					if (a < 1) {
						printf("a = %ld i = %ld j = %ld\n", a, i, j);
						fflush(stdout);
						return error("elements must start with 1");
						}
					if (a > v) {
						printf("a = %ld\n", a);
						fflush(stdout);
						return error("a > v");
						}
					a--;
					B.m_ii(j, a);
					}
				printf("base block: ");
				B.println();
				B.swap((VECTOR_OP) BB.s_i(i));
				
				}
			geo.init(&BB, nb_input - 1);
			geo.span_design(&group_generators, TRUE /* f_v */, FALSE /* f_vv */);
			}


		if (f_test2design) {
			
			if (!geo.s_I()->test_2_design(&k, &l)) {
				printf("WARNING: NOT A 2-DESIGN !!!\n");
				}
			else {
				printf("a 2-(%ld,%ld,%ld) design\n", geo.s_I()->s_hi(), k, l);
				if (f_first_design) {
					design_k = k;
					design_lambda = l;
					f_first_design = FALSE;
					}
				else {
					if (design_k != k)
						printf("warning: design_k changes\n");
					if (design_lambda != l)
						printf("warning: design_lambda changes\n");
					}
				}
			}
		if (f_v) {
			printf("geo no %ld\n", nb_input);
			geo.s_I()->Print();
			}

		geo.canonical_form(f_tdo, f_tdo2, 
			f_V, decomp_V, f_B, decomp_B, 
			TRUE /* f_v */, FALSE /* f_vv */);
		if (f_v) {
			printf("canonical geo no %ld\n", nb_input);
			geo.s_I()->Print();
			}
		geo.calc_hash(TRUE /* f_v */, FALSE /* f_vv */);
		
		bt->len(&l);
		printf("current db len = %ld\n", l);
		fflush(stdout);
	
		if (db->close() != OK) {
			return error("db_geo_by_base_blocks_add() error in db->close()");
			}
		
		db_geo_by_base_blocks_find_isomorphic_via_hash(db, 
			stdout /* fp_txt */, &geo, TRUE /* f_verbose */, 
			&idx /* idx_btree_hash */, &f_found);
		
		if (db->open() != OK) {
			return error("db_geo_by_base_blocks_add() error in db->open()");
			}
		
		if (f_found) {
			printf("geometry already there, skipping !\n");
			}
		else {
			if (f_do_not_add) {
				printf("could not find this geometry\n");
				}
			else {
				nb_geo++;
				printf("adding new geometry no %ld (nb_input = %ld)\n", l, nb_input - 1);
				printf("this is the %ld-th geometry\n", nb_geo);
				if (f_tda) {
					geo.s_I()->calc_TDA(TDA_fp, l, TRUE /* f_row_perms */, 
						TRUE /* f_has_aut */, &aut1, 
						&a1, &a2, &a3, &a4, 
						TRUE /* f_v */, FALSE /* f_vv */, FALSE /* f_incidences */);
					}
				if (db->add_op(&geo) != OK) {
					return error("db_geo_by_base_blocks_add() error in db->add_op()");
					}
				}
			}

		}
	bt->len(&l);
	fprintf(OUT_fp, "-1 %ld geometries\n\n", nb_geo);
	fclose(OUT_fp);

	fclose(in_fp);
	
	if (f_tda) {
		fclose(TDA_fp);
		}
	
	db->close();
	e_geo_by_base_blocks_db(db);
	
	return OK;
}

