/* db.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef DB_INCLUDED
#define DB_INCLUDED

#ifndef DIVS_INCLUDED
#include <DISCRETA/divs.h>
#endif

#define BTREEMAXKEYLEN 48

typedef struct keycarrier {
	BYTE c[BTREEMAXKEYLEN];
} KEYCARRIER;

typedef struct datatype {
	INT4 datref;
	INT4 data_size;
} DATATYPE;

typedef KEYCARRIER KEYTYPE;

class bt_key_ob : public VECTOR_OB {
public:
	INTEGER_OP s_type() { 
		return((INTEGER_OP)s_i(0)); };
	INT s_type_i() { 
		return(s_type()->s_i()); };
	/* type:
	 * 0 = INT
	 * 1 = STRING
	 * 2 = INT VEC
	 */
	INTEGER_OP s_output_size() { 
		return((INTEGER_OP)s_i(1)); };
	INT s_output_size_i() { 
		return(s_output_size()->s_i()); };
	/* fuer type INT nur moeglich: 1, 2, 4, evtl 8 */

	INTEGER_OP s_int_vec_fst() { 
		return((INTEGER_OP)s_i(2)); };
	INT s_int_vec_fst_i() { 
		return(s_int_vec_fst()->s_i()); };
	INTEGER_OP s_int_vec_len() { 
		return((INTEGER_OP)s_i(3)); };
	INT s_int_vec_len_i() { 
		return(s_int_vec_len()->s_i()); };
	
	INTEGER_OP s_field1() { 
		return((INTEGER_OP)s_i(4)); };
	INT s_field1_i() { 
		return(s_field1()->s_i()); };
	INTEGER_OP s_field2() { 
		return((INTEGER_OP)s_i(5)); };
	INT s_field2_i() { 
		return(s_field2()->s_i()); };
	INTEGER_OP s_f_ascending() { 
		return((INTEGER_OP)s_i(6)); };
	INT s_f_ascending_i() { 
		return(s_f_ascending()->s_i()); };

INT field_name(INT i, INT j, BYTE *str);
INT init(INT type, 
	INT output_size, INT field1, INT field2);
INT init_INT4(INT field1, INT field2);
INT init_INT2(INT field1, INT field2);
INT init_STRING(INT output_size, 
	INT field1, INT field2);
INT init_INT4_VEC(INT field1, INT field2, 
	INT vec_fst, INT vec_len);
INT init_INT2_VEC(INT field1, INT field2, 
	INT vec_fst, INT vec_len);
INT sprint(BYTE *s);
};

class bt_ob : public VECTOR_OB {
public:
	INTEGER_OP s_f_duplicatekeys() { 
		return((INTEGER_OP)s_i(0)); };
	INT s_f_duplicatekeys_i() { 
		return(s_f_duplicatekeys()->s_i()); };

	INTEGER_OP s_nb_bt_key() { 
		return((INTEGER_OP)s_i(1)); };
	INT s_nb_bt_key_i() { 
		return(s_nb_bt_key()->s_i()); };
	VECTOR_OP s_bt_key() { 
		return((VECTOR_OP)s_i(2)); };
	BT_KEY_OP s_bt_key_i(INT i) { 
		return((BT_KEY_OP)s_bt_key()->s_i(i)); };

	STRING_OP s_filename() {
		return((STRING_OP)s_i(3)); };
	BYTE *s_filename_s() {
		return(s_filename()->s_str()); };
	INTEGER_OP s_f_open() { 
		return((INTEGER_OP)s_i(4)); };
	INT s_f_open_i() { 
		return(s_f_open()->s_i()); };
	INTEGER_OP s_stream() { 
		return((INTEGER_OP)s_i(5)); };
	FILE *s_stream_i() { 
		return((FILE *)s_stream()->s_i()); };
	INTEGER_OP s_buf_idx() { 
		return((INTEGER_OP)s_i(6)); };
	INT s_buf_idx_i() { 
		return(s_buf_idx()->s_i()); };
	INTEGER_OP s_Root() { 
		return((INTEGER_OP)s_i(7)); };
	INT s_Root_i() { 
		return(s_Root()->s_i()); };
	INTEGER_OP s_FreeRec() { 
		return((INTEGER_OP)s_i(8)); };
	INT s_FreeRec_i() { 
		return(s_FreeRec()->s_i()); };
	INTEGER_OP s_AllocRec() { 
		return((INTEGER_OP)s_i(9)); };
	INT s_AllocRec_i() { 
		return(s_AllocRec()->s_i()); };
	INTEGER_OP s_btree_idx() { 
		return((INTEGER_OP)s_i(10)); };
	INT s_btree_idx_i() { 
		return(s_btree_idx()->s_i()); };

INT field_name(INT i, INT j, BYTE *str);
INT init(BYTE *file_name, 
	INT f_duplicatekeys, INT btree_idx);
INT add_key_INT4(INT field1, INT field2);
INT add_key_INT2(INT field1, INT field2);
INT add_key_STRING(INT output_size, INT field1, INT field2);
INT add_key_INT4_VEC(INT field1, INT field2, 
	INT vec_fst, INT vec_len);
INT add_key_INT2_VEC(INT field1, INT field2, 
	INT vec_fst, INT vec_len);
INT sprint(BYTE *s);
/* INT get_cmp_func(void **v);
INT get_sprint_func(void **v); */
INT create();
INT open();
INT close();
INT len(INT *len);
INT search(void *pSearchKey, DATATYPE *pData, INT *idx, INT *f_found);
INT ith(INT l, KEYTYPE *key, DATATYPE *data);
INT ith0(INT l, KEYTYPE *key, DATATYPE *data);
INT insert_key(KEYTYPE *pKey, DATATYPE *pData);
INT add(KEYTYPE *pKey, DATATYPE *pData);
INT del(INT idx);
INT print();
INT print_pages();
INT print_page(INT x);
};

/* OS-Interface: */

FILE *OPEN1(BYTE *name);
FILE *CREATE1(BYTE *name);
INT CLOSE1(FILE *f);
INT SEEK1(FILE *f, INT l);

class db_ob : public VECTOR_OB {
public:
	INTEGER_OP s_nb_btree() { 
		return((INTEGER_OP)s_i(0)); };
	INT s_nb_btree_i() { 
		return(s_nb_btree()->s_i()); };
	VECTOR_OP s_btree() { 
		return((VECTOR_OP)s_i(1)); };
	BAYERTREE_OP s_btree_i(INT i) { 
		return((BAYERTREE_OP)s_btree()->s_i(i)); };
	STRING_OP s_filename() {
		return((STRING_OP)s_i(2)); };
	BYTE *s_filename_s() {
		return(s_filename()->s_str()); };
	INTEGER_OP s_sym_type() { 
		return((INTEGER_OP)s_i(3)); };
	INT s_sym_type_i() { 
		return(s_sym_type()->s_i()); };
	INTEGER_OP s_f_open() { 
		return((INTEGER_OP)s_i(4)); };
	INT s_f_open_i() { 
		return(s_f_open()->s_i()); };
	INTEGER_OP s_stream() { 
		return((INTEGER_OP)s_i(5)); };
	FILE *s_stream_i() { 
		return((FILE *)s_stream()->s_i()); };
	INTEGER_OP s_file_size() { 
		return((INTEGER_OP)s_i(6)); };
	INT s_file_size_i() { 
		return(s_file_size()->s_i()); };

INT field_name(INT i, INT j, BYTE *str);
INT init(BYTE *file_name, INT sym_type);
INT sprint(BYTE *s);
INT add_btree(BYTE *file_name, INT f_duplicatekeys, INT btree_idx);
INT put_file_size();
INT get_file_size();
INT delete_files();
INT create();
INT open();
INT open_DB();
INT close();
INT close_DB();
INT free_data_DB(INT datref, INT size);
INT add_data_DB(void *data, INT size, INT *datref);
INT ith(INT btree_idx, INT l, KEYTYPE *key, DATATYPE *data);
INT add_op(SYM_OP v);
INT add_op1(SYM_OP v, INT f_v, INT f_vv);
INT del_op(SYM_OP v, INT datref, INT size);
INT check_DB();
INT get_op(DATATYPE *data_type, SYM_OP op);
INT print(INT btree_idx);
INT restore(FILE *fp_txt);
INT restore_from_file(BYTE *restore_fname, 
	INT f_delete_afterwards, FILE *fp_txt);
};

INT database_init();
void database_exit();
INT root_buf_alloc();
/* returns -1 if no buffer free */
void root_buf_free(INT i);

/* bt_key.C: */
INT bt_key_init(void);
void bt_key_exit(void);
INT lex_cmp_init();
INT bt_lexicographic_cmp(BYTE *p1, BYTE *p2, INT *res);
INT bt_key_int_cmp(BYTE *p1, BYTE *p2, INT *res);
INT bt_key_int2_cmp(BYTE *p1, BYTE *p2, INT *res);
INT bt_sprint_int(BYTE *key, BYTE *str_255);
INT bt_sprint_int2(BYTE *key, BYTE *str_255);
INT bt_sprint_string(BYTE *key, BYTE *str_255);
INT bt_key_sprint(BYTE *key, VECTOR_OP V, INT len, BYTE *str);
INT bt_sprint_INT4(BYTE **p_key, BYTE *str);
INT bt_sprint_INT2(BYTE **p_key, BYTE *str);
INT bt_key_cmp(BYTE *key1, BYTE *key2, VECTOR_OP V, INT len, INT *res);
INT bt_lex_n_cmp(BYTE *p1, BYTE *p2, INT n, INT *res);
INT bt_cmp_INT4(BYTE **p_key1, BYTE **p_key2, INT *res);
INT bt_cmp_INT2(BYTE **p_key1, BYTE **p_key2, INT *res);
INT bt_key_fill_in(BYTE *key, VECTOR_OP theOp, VECTOR_OP V, INT len);
INT bt_key_INT4_fill_in(BYTE **p_key, SYM_OP key_op);
INT bt_key_INT4_fill_in_(BYTE **p_key, INT value);
INT bt_key_INT2_fill_in(BYTE **p_key, SYM_OP key_op);
INT bt_key_INT2_fill_in_(BYTE **p_key, INT value);
INT bt_key_STR_fill_in(BYTE **p_key, INT output_size, SYM_OP key_op);
INT bt_key_change_INT(BYTE *key, VECTOR_OP V, INT idx, INT new_value);
INT bt_get_INT4(BYTE **p_key, INT *p_i);
INT bt_get_INT2(BYTE **p_key, INT *p_i);

/* db1.C: */
INT i_custom_db(DATABASE_OP *db_fg);
void e_custom_db(DATABASE_OP db);
INT c_db_custom_create_and_close();
INT c_db_custom_info();
INT db_custom_import(BYTE *fname, 
	BYTE begin_mark, BYTE end_mark, 
	BYTE eof_mark, BYTE comment_mark);
INT db_custom_HTML();
INT db_custom_exp_tex(INT btree_idx, BYTE *tex_file_name);



#endif /* DB_INCLUDED */

