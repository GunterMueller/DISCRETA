/* divs.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef DIVS_INCLUDED
#define DIVS_INCLUDED

class bitvec_ob : public VECTOR_OB {
public:
	VECTOR_OP s_ints() {
		return((VECTOR_OP) s_i(0)); };
	INTEGER_OP s_ints_i(INT i) {
		return((INTEGER_OP) s_ints()->s_i(i)); };
	INT s_ints_ii(INT i) {
		return(s_ints_i(i)->s_i()); };

INT field_name(INT i, INT j, BYTE *str);
INT init_ascii(BYTE *s);
INT print();
INT sprint(BYTE *s);
};

class string_ob : public SYM_OB {
public:
INT freeself();
INT init(BYTE *str);
BYTE *s_str();
INT c_str(BYTE *str);
/* overwrites an existing string; does not copy str. */
INT copy(STRING_OP b);
INT append(STRING_OP b);
/* concatenates strings: first this then b, result into this. */
/* warning: former routinue string_append(): first b then a, result into b */
INT sprint(BYTE *str);
/* appends to str. writes to maximal strlength of 200. */
INT sscan(BYTE *str);
/* in iof.C: */
INT write_mem(MEM_OP mem);
INT read_mem(MEM_OP mem);
};
INT string_test();

class mem_ob : public SYM_OB {
public:
/*
 * INT - offset - 3 + 0: allocated length in BYTEs
 *              - 3 + 1: used length in BYTES
 *              - 3 + 2: current pointer: 0 <= pointer < used length. 
 */

INT freeself();
INT s_alloc_length_i();
INT s_used_length_i();
INT s_cur_pointer_i();
INT *s_alloc_length_pi();
INT *s_used_length_pi();
INT *s_cur_pointer_pi();
INT c_alloc_length(INT alloc_length);
INT c_used_length(INT used_length);
INT c_cur_pointer(INT cur_pointer);
INT init(INT length, BYTE *d);
INT compress(INT f_verbose);
INT decompress(INT f_verbose);
INT mult_char_c(BYTE c);
INT alloc(INT length);
/* setzt alloc_length auf length + MEM_OVERSIZE, 
 *       used_length auf length, 
 *       cur_pointer auf 0.
 */
INT copy(MEM_OP b);
INT append(INT length, BYTE *data);
INT realloc(INT new_length);
INT sprint(BYTE *str);
/* appends to str. writes to maximal strlength of 200. */
INT write_char(BYTE c);
INT read_char(BYTE *c);
INT write_int(INT i);
INT read_int(INT *i);
INT length_in_entries(INT entry_size);
BYTE *pc_ith(INT i, INT entry_size);
INT insert_at(INT entry_size, BYTE *data, INT i);
INT isd(INT entry_size, 
	BYTE *add_this, INT *f_found, 
	INT (*cmp_func_data)(void *a, void *b, void *data), 
	void *data);
INT sd(INT entry_size, void *search_this, 
	INT *idx, INT *f_found, 
	INT (*cmp_func_data)(void *a, void *b, void *data), 
	void *data);
INT print_data(INT entry_size, 
	INT (*print_data)(void *a, void *data), 
	void *data);
/* in iof.C: */
INT calc_size_on_file();
INT write_mem(MEM_OP mem);
INT read_mem(MEM_OP mem);
INT slurp_file(BYTE *fname, INT f_v, FILE *fp_txt);
INT write_to_file(BYTE *fname, INT f_v, FILE *fp_txt);
};

class conti_ob : public VECTOR_OB {
public:
	INTEGER_OP s_nb_L() { 
		return((INTEGER_OP)s_i(0)); };
	INT s_nb_L_i() { 
		return(s_nb_L()->s_i()); };
	VECTOR_OP s_L() { 
		return((VECTOR_OP)s_i(1)); };
	STRING_OP s_L_i(INT i) { 
		return((STRING_OP)s_L()->s_i(i)); };
	BYTE *s_L_i_str(INT i) { 
		return(s_L_i(i)->s_str()); };
	INTEGER_OP s_nb_C() { 
		return((INTEGER_OP)s_i(2)); };
	INT s_nb_C_i() { 
		return(s_nb_C()->s_i()); };
	VECTOR_OP s_C() { 
		return((VECTOR_OP)s_i(3)); };
	SYM_OP s_C_i(INT i) { 
		return(s_C()->s_i(i)); };
	VECTOR_OP s_CS() { 
		return((VECTOR_OP)s_i(4)); };
	STRING_OP s_CS_i(INT i) { 
		return((STRING_OP)s_CS()->s_i(i)); };
	BYTE *s_CS_i_str(INT i) { 
		return(s_CS_i(i)->s_str()); };

INT field_name(INT i, INT j, BYTE *str);
INT init();
INT sprint(BYTE *s);
INT add_label(BYTE *label);
INT add_object(BYTE *label, SYM_OP op);
INT find_object(BYTE *label, INT *idx, INT *f_found);
};

INT ascii_formatter(VECTOR_OP text, VECTOR_OP in, INT length);
INT ascii_formatter_paragraph(VECTOR_OP text, BYTE *heading, INT indent, 
	VECTOR_OP in, INT length);
int find_wrap_position(BYTE *p, INT max);

#endif /* DIVS_INCLUDED */


