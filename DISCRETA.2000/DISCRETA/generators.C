/* generators.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef GENERATORS_TRUE

#include <DISCRETA/perm.h>
#include <DISCRETA/generators.h>

#if TEXDOCU
INT GENERATORS_OB::field_name(INT i, INT j, BYTE *str)
#endif
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "id"; break;
	case 1: s = "deg"; break;
	case 2: s = "gen"; break;
	case 3: s = "label"; break;
	case 4: s = "label_tex"; break;
	case 5: s = "order"; break;
	case 6: s = "go"; break;
	default:
		return error("GENERATORS::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

#if TEXDOCU
INT GENERATORS_OB::init()
#endif
{
	m_il(7);
	c_obj_k(GENERATORS_KIND);
	s_id()->m_i(-1);
	s_deg()->m_i(0);
	s_gen()->m_il(0);
	s_label()->init("");
	s_label_tex()->init("");
	s_order()->init("");
	s_go()->m_i_i(0);
	return OK;
}

#if TEXDOCU
INT GENERATORS_OB::sprint(BYTE *s)
#endif
{
	char str[1024];

	str[0] = 0;
	sprintf(str + strlen(str), "GENERATORS go=%s for %s\tid=%ld\tdegree=%ld ", 
		s_order_s(), s_label_s(), s_id_i(), s_deg_i());
	if (strlen(s) + strlen(str) < 256)
		strcat(s, str);
	return OK;
}

/* 
 * btree #0: id
 * btree #1: name, id
 * btree #2: go, id
 * btree #3: deg, id
 */

#if TEXDOCU
INT i_db_generators(DATABASE_OP *p_db, BYTE *path)
#endif
/* path without trailing '/' */
{
	DATABASE_OP db;
	BAYERTREE_OP bt;
	BYTE file_db[1024];
	BYTE file_idx0[1024];
	BYTE file_idx1[1024];
	BYTE file_idx2[1024];
	BYTE file_idx3[1024];

	strcpy(file_db, path);
	strcat(file_db, "/generators.db");
	strcpy(file_idx0, path);
	strcat(file_idx0, "/generators0.idx");
	strcpy(file_idx1, path);
	strcat(file_idx1, "/generators1.idx");
	strcpy(file_idx2, path);
	strcat(file_idx2, "/generators2.idx");
	strcpy(file_idx3, path);
	strcat(file_idx3, "/generators3.idx");
	
	db = (DATABASE_OP) callocobject("generators.C: i_generators_dp()");
	db->init(file_db, GENERATORS_KIND);

	db->add_btree(file_idx0, 
		TRUE /* f_duplicatekeys */, 
		0 /* btree_idx */ );
	bt = db->s_btree_i(0);
	bt->add_key_INT4(0 /* id */, 0);

	db->add_btree(file_idx1, 
		TRUE /* f_duplicatekeys */, 
		1 /* btree_idx */ );
	bt = db->s_btree_i(1);
	bt->add_key_STRING(40, 3 /* label */, 0);
	bt->add_key_INT4(0 /* id */, 0);

	db->add_btree(file_idx2, 
		TRUE /* f_duplicatekeys */, 
		2 /* btree_idx */ );
	bt = db->s_btree_i(2);
	bt->add_key_STRING(40, 5 /* order */, 0);
	bt->add_key_INT4(0 /* id */, 0);

	db->add_btree(file_idx3, 
		TRUE /* f_duplicatekeys */, 
		3 /* btree_idx */ );
	bt = db->s_btree_i(3);
	bt->add_key_INT4(1 /* deg */, 0);
	bt->add_key_INT4(0 /* id */, 0);

	*p_db = db;
	return OK;
}

#if TEXDOCU
void e_db_generators(DATABASE_OP db)
#else
Frees the object.
#endif
{
	freeall(db);
}

#if TEXDOCU
INT do_db_generators_create_and_close(BYTE *path)
#else
Creates an empty database and closes it.
#endif
{
	DATABASE_OP db;
	
	if (i_db_generators(&db, path) != OK) {
		return error("do_db_generators_create_and_close() i_db_dp()");
		}
	if (db->create() != OK) {
		return error("do_db_generators_create_and_close() create()");
		}
	if (db->close() != OK) {
		return error("do_db_generators_create_and_close() close()");
		}
	e_db_generators(db);
	return OK;
}

#if TEXDOCU
INT do_db_generators_info(BYTE *path, INT btree_idx)
#else
//TEX
Calls the btree function "print\_pages" for the btree index "btree\_idx".
///TEX
#endif
{
	BAYERTREE_OP btree;
	DATABASE_OP db;
	
	if (i_db_generators(&db, path) != OK) {
		return error("do_db_generators_info() i_db_generators()");
		}
	/* da_print(db, 0); */
	btree = db->s_btree_i(btree_idx);
	btree->print_pages();
	e_db_generators(db);
	return OK;
}

#if TEXDOCU
INT do_db_generators_dump(BYTE *path, INT btree_idx)
#endif
{
	BAYERTREE_OP btree;
	DATABASE_OP db;
	INT i, l;
	KEYTYPE key;
	DATATYPE data;
	GENERATORS_OB g;
	
	if (i_db_generators(&db, path) != OK) {
		return error("do_db_generators_dump() i_db_generators()");
		}
	btree = db->s_btree_i(btree_idx);
	//btree->print_pages();
	
	if (db->open() != OK) {
		return error("db_generators_dump() error in DB::open()");
		}
	if (btree->len(&l) != OK) {
		return error("db_generators_dump() error in BT::len()");
		}
	for (i = 0; i < l; i++) {
		if (btree->ith(i, &key, &data) != OK) {
			return error("db_generators_dump() error in BT::ith()");
			}
		if (db->get_op(&data, &g) != OK) {
			return error("db_generators_dump() error in DB::get_op()");
			}
		printf("%ld: ", i);
		g.println();
		fflush(stdout);
		} // next i
	
	if (db->close() != OK) {
		return error("db_generators_dump() error in DB::close()");
		}
	e_db_generators(db);
	return OK;
}

#if TEXDOCU
INT db_generators_highest_id(DATABASE_OP db, INT *highest_id)
#else
returns the highest id which is used in the database.
#endif
{
	BAYERTREE_OP btree;
	INT btree_idx;
	KEYTYPE key;
	INT4 *pi;
	INT id, len;
	DATATYPE data;

	btree_idx = 0;
	btree = db->s_btree_i(btree_idx);
	if (btree->open() != OK) {
		return error("db_generators_highest_id() error in BT::open()");
		}
	if (btree->len(&len) != OK) {
		return error("db_generators_highest_id() error in BT::len()");
		}
	if (len == 0) {
		btree->close();
		*highest_id = -1;
		return OK;
		}
	if (btree->ith(len - 1, &key, &data) != OK) {
		return error("db_generators_highest_id() error in BT::ith()");
		}
	pi = (INT4 *) &key;
	id = pi[0];
	*highest_id = id;
	if (btree->close() != OK) {
		return error("db_generators_highest_id() error in BT::close()");
		}
	return OK;
}

#define GENERATORS_DEBUG_LOAD_ID
#define GENERATORS_DEBUG_LOAD_ID_VERBOSE


#if TEXDOCU
INT db_generators_load_id(DATABASE_OP db, INT id, GENERATORS_OP p)
#else
Loads the dataset by id.
#endif
{
	BAYERTREE_OP btree;
	INT btree_idx;
	KEYTYPE key, key1;
	INT4 *pi, *pi1;
	INT idx, f_found;
	DATATYPE data;

#ifdef GENERATORS_DEBUG_LOAD_ID
	printf("db_generators_load_id() loading id = %ld\n", id);
#endif
	btree_idx = 0; /* id */
	btree = db->s_btree_i(btree_idx);
	if (btree->open() != OK) {
		return error("db_generators_load_id() error in BT::open()");
		}
	pi = (INT4 *) &key;
	pi1 = (INT4 *) &key1;
	pi[0] = id;
	pi[1] = 0;
	pi[2] = 0;
	pi[3] = 0;
	pi[4] = 0;
	if (btree->search(pi /* pSearchKey */, 
		&data, &idx, &f_found) != OK) {
		return error("db_generators_load_id() error in BT::search()");
		}
#ifdef GENERATORS_DEBUG_LOAD_ID_VERBOSE
	printf("db_generators_load_id(): f_found = %ld %ld "
		"idx = %ld\n", 
		f_found, (INT) pi[0], idx);
#endif
	if (!f_found)
		return error("db_generators_load_id(): id not found");

	idx--;
	if (db->open_DB() != OK) {
		return error("db_generators_load_id() error in DB::open_DB()");
		}
	if (btree->ith(idx, &key1, &data) != OK) {
		return error("db_generators_load_id() error in BT::ith()");
		}
	if (pi1[0] != id) {
		return error("db_generators_load_id() pi1[0] != id (id not found)");
		}
	if (db->get_op(&data, p) != OK) {
		return error("db_generators_load_id() error in DB::get_op()");
		}
	if (db->close_DB() != OK) {
		return error("db_generators_load_id() error in DB::close()");
		}
	if (btree->close() != OK) {
		return error("db_generators_load_id() error in DB::close()");
		}
	return OK;
}

#if TEXDOCU
INT db_generators_add(DATABASE_OP db, VECTOR_OP gen, 
	BYTE *label, BYTE *label_tex)
#endif
{
	LABRA_OB L;
	GENERATORS_OB g;
	INT deg, highest_id;
	SYM_OB go;
	BYTE str[1024], str1[1024];
	INT i, l;

	g.init();
	gen->copy(g.s_gen());
	g.s_label()->init(label);
	g.s_label_tex()->init(label_tex);
	deg = perm_vec_get_degree(gen);
	g.s_deg()->m_i(deg);
	
	db_generators_highest_id(db, &highest_id);
	g.s_id()->m_i(highest_id + 1);

	reduce_generators_labra(gen, &go, FALSE /* f_verbose */, &L);
	str[0] = 0;
	str1[0] = 0;
	go.sprint(str);
	l = strlen(str);
	for (i = l; i < 20; i++) {
		strcat(str1, " ");
		}
	strcat(str1, str);
	printf("$|G| = %s$\\\\\n", str1);
	
	g.s_order()->init(str1);
	go.copy(g.s_go());
	printf("adding generators to database:\n");
	g.println();
	fflush(stdout);

	
	if (db->open_DB() != OK) {
		return error("db_generators_add() error in DB::open_DB()");
		}
	if (db->add_op(&g) != OK) {
		return error("db_generators_add() error in DB::add_op()");
		}
	if (db->close_DB() != OK) {
		return error("db_generators_add() error in DB::close()");
		}
	return OK;
}


#endif /* GENERATORS_TRUE */


