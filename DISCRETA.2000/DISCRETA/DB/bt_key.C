/* bt_key.C */

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

/* 
 * BT_KEY
 */

INT BT_KEY_OB::field_name(INT i, INT j, BYTE *str)
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "type"; break;
	case 1: s = "output_size"; break;
	case 2: s = "int_vec_fst"; break;
	case 3: s = "int_vec_len"; break;
	case 4: s = "field1"; break;
	case 5: s = "field2"; break;
	case 6: s = "f_ascending"; break;
	default:
		return error("BT_KEY::field_name()|"
		"i out of range");
	}
	strcpy(str, s);
	return OK;
}

INT BT_KEY_OB::init(INT type, INT output_size, INT field1, INT field2)
{
	INT erg = OK;
	
	erg += m_il(7);
	c_obj_k(BT_KEY_KIND);
	s_type()->m_i(type);
	s_output_size()->m_i(output_size);
	s_field1()->m_i(field1);
	s_field2()->m_i(field2);
	s_int_vec_fst()->m_i(0);
	s_int_vec_len()->m_i(0);
	s_f_ascending()->m_i(TRUE);
	return erg;
}

INT BT_KEY_OB::init_INT4(INT field1, INT field2)
{
	INT erg = OK;
	
	erg += init(0 /* type */, 4 /* output_size */, field1, field2);
	return erg;
}

INT BT_KEY_OB::init_INT2(INT field1, INT field2)
{
	INT erg = OK;
	
	erg += init(0 /* type */, 2 /* output_size */, field1, field2);
	return erg;
}

INT BT_KEY_OB::init_STRING(INT output_size, INT field1, INT field2)
{
	INT erg = OK;
	
	erg += init(1 /* type */, output_size, field1, field2);
	return erg;
}

INT BT_KEY_OB::init_INT4_VEC(INT field1, INT field2, INT vec_fst, INT vec_len)
{
	INT erg = OK;
	
	erg += init(2 /* type */, 4 /* output_size */, field1, field2);
	s_int_vec_fst()->m_i(vec_fst);
	s_int_vec_len()->m_i(vec_len);
	return erg;
}

INT BT_KEY_OB::init_INT2_VEC(INT field1, INT field2, INT vec_fst, INT vec_len)
{
	INT erg = OK;
	
	erg += init(2 /* type */, 2 /* output_size */, field1, field2);
	s_int_vec_fst()->m_i(vec_fst);
	s_int_vec_len()->m_i(vec_len);
	return erg;
}

INT BT_KEY_OB::sprint(BYTE *s)
{
	BYTE str[512];
	
	sprintf(str, "BT_KEY: type = %ld output_size = %ld "
		"field1 = %ld field2 = %ld", 
		s_type_i(), s_output_size_i(), 
		s_field1_i(), s_field2_i() );
	if (strlen(s) + strlen(str) < 200)
		strcat(s, str);
	return OK;
}

static BYTE *lex_cmp = NIL; /* [256] */

INT bt_key_init(void)
{
	INT size;
	
	size = 256 * sizeof(BYTE);
	lex_cmp = (BYTE *) my_malloc(size, "bt_key_init()");
	if (lex_cmp == NIL) {
		printf("no memory for lex_cmp\n");
		return ERROR;
		}
	lex_cmp_init();
	return (OK);
}

void bt_key_exit(void)
{
	if (lex_cmp) {
		my_free(lex_cmp);
		lex_cmp = NIL;
		}
}

INT lex_cmp_init(void)
{
	INT i;
	UBYTE *ulex = (UBYTE *)lex_cmp;
	
	for (i = 0; i < 256; i++)
		ulex[i] = (UBYTE)i;
#ifdef SYSTEMMAC
	ulex[0x80] = ulex[0x81] = ulex[0xae] = 
	ulex[0xcb] = ulex[0xcc] = (UBYTE)0x61; 
		/* Aoben2Punkte, AobenKringel, 
		 * AE, Agrave, Acirconflex -> a */
	ulex[0x82] = (UBYTE)0x63; /* C cedille -> c */
	ulex[0x83] = (UBYTE)0x65; /* E grave -> e */
	ulex[0x84] = (UBYTE)0x6e; /* N circonflex -> n */
	ulex[0x85] = ulex[0xaf] = ulex[0xcd] = (UBYTE)0x6f;
		/* O 2 Punkte, O durchgestrichen, 
		 * O circonflex -> o */
	ulex[0x86] = (UBYTE)0x75; /* U 2 Punkte -> u */
	ulex[0x87] = ulex[0x88] = ulex[0x89] = 
	ulex[0x8a] = ulex[0x8b] = ulex[0x8c] = 
	ulex[0xbb] = ulex[0xbe] = (UBYTE)0x61; 
		/* a akut, a grave, a dach, 
		 * a 2 Punkte, a circonflex, 
		 * a oben Kringel, a unterstrichen, 
		 * ae -> a */
	ulex[0x8d] = (UBYTE)0x63; /* c cedille -> c */
	ulex[0x8e] = ulex[0x8f] = ulex[0x90] = 
	ulex[0x91] = (UBYTE)0x65;
		/* e akut, e grave, e dach, e 2 punkte -> e */
	ulex[0x92] = ulex[0x93] = ulex[0x94] = 
	ulex[0x95] = (UBYTE)0x69;
		/* i akut, i grave, i dach, i 2 punkte -> i */
	ulex[0x96] = (UBYTE)0x6e;
		/* n circonflex -> n */
	ulex[0x97] = ulex[0x98] = ulex[0x99] = 
	ulex[0x9a] = ulex[0x9b] = ulex[0xbc] = 
	ulex[0xbf] = ulex[0xcf] = (UBYTE)0x6f;
		/* o akut, o grave, o dach, 
		 * o 2 punkte, o circonflex, 
		 * o unterstrichen, o durchgestrichen, 
		 * oe -> o */
	ulex[0x9c] = ulex[0x9d] = ulex[0x9e] = 
	ulex[0x9f] = (UBYTE)0x75;
		/* u akut, u grave, u dach, u 2 punkte -> u */
	ulex[0xa7] = (UBYTE)0x73;
		/* sharp s -> s */
#endif
	for (i = 0x41; i <= 0x5a; i++) {
		ulex[i] = i + 0x20;
		}
	return(OK);
}

INT bt_lexicographic_cmp(BYTE *p1, BYTE *p2, INT *res)
{
	UBYTE c1, c2, d1, d2;
	
	while (TRUE) {
		c1 = *(UBYTE *)p1;
		c2 = *(UBYTE *)p2;
		d1 = lex_cmp[c1];
		d2 = lex_cmp[c2];
		if (d1 < d2) {
			/* Beachte: hier auch Abbruch, 
			 * falls nur String 1 zuende. */
			*res = -1;
			return(OK);
			}
		if (d1 > d2) {
			*res = 1;
			return(OK);
			}
		/* Gleichheit */
		if (c1 == 0) {
			/* Abbruchbedingung: 
			 * beide Strings zuende. */
			*res = 0;
			return(OK);
			}
		p1++;
		p2++;
		}
}

INT bt_key_int_cmp(BYTE *p1, BYTE *p2, INT *res)
/* bt_key_long_cmp */
{
	INT4 *p_l1, *p_l2;
	
	p_l1 = (INT4 *) p1;
	p_l2 = (INT4 *) p2;
	if (*p_l1 < *p_l2) {
		*res = -1;
		return(OK);
		}
	if (*p_l1 > *p_l2) {
		*res = 1;
		return(OK);
		}
	*res = 0;
	return(OK);
}

INT bt_key_int2_cmp(BYTE *p1, BYTE *p2, INT *res)
{
	INT4 *p_l1, *p_l2;
	
	p_l1 = (INT4 *) p1;
	p_l2 = (INT4 *) p2;
	if (*p_l1 < *p_l2) {
		*res = -1;
		return(OK);
		}
	if (*p_l1 > *p_l2) {
		*res = 1;
		return(OK);
		}
	if (p_l1[1] < p_l2[1]) {
		*res = -1;
		return(OK);
		}
	if (p_l1[1] > p_l2[1]) {
		*res = 1;
		return(OK);
		}
	*res = 0;
	return(OK);
}

INT bt_sprint_int(BYTE *key, BYTE *str_255)
{
	INT4 *pi;
	
	pi = (INT4 *) key;
	sprintf(str_255, "%ld", (long) *pi);
	return OK;
}

INT bt_sprint_int2(BYTE *key, BYTE *str_255)
{
	INT4 *pi;
	
	pi = (INT4 *) key;
	sprintf(str_255, "%ld %ld", 
		(long) pi[0], (long) pi[1]);
	return OK;
}

INT bt_sprint_string(BYTE *key, BYTE *str_255)
{
	sprintf(str_255, "%s", key);
	return OK;
}

INT bt_key_sprint(BYTE *key, VECTOR_OP V, INT len, BYTE *str)
/* es wird an str angehaengt */
{
	BT_KEY_OP bt_key;
	BYTE *the_key = key;
	BYTE c;
	INT i, j, l, l1, type, output_size;
	
	for (i = 0; i < len; i++) {
		bt_key = (BT_KEY_OP) V->s_i(i);
		type = bt_key->s_type_i();
		output_size = bt_key->s_output_size_i();
		if (type == 0) { /* INT */
			if (output_size == 4) {
				bt_sprint_INT4(&the_key, str);
				}
			else if (output_size == 2) {
				bt_sprint_INT2(&the_key, str);
				}
			else {
				Srfs("bt_key_sprint", 
				"unknown output size");
				return ERROR;
				}
			}
		else if (type == 1) { /* STRING */
			l = strlen(str);
			for (j = 0; j < output_size; j++) {
				if (the_key[j] == 0) {
					break;
					}
				}
			l1 = j;
			for (j = 0; j < output_size; j++) {
				if (j < l1)
					c = *the_key;
				else
					c = ' ';
				str[l + j] = c;
				the_key++;
				}
			str[l + j] = 0;
			}
		else if (type == 2) { /* INT VEC */
			strcat(str, "(");
			for (j = 0; 
				j < bt_key->s_int_vec_len_i(); 
				j++) {
				if (output_size == 4) {
					bt_sprint_INT4(&the_key, str);
					}
				else if (output_size == 2) {
					bt_sprint_INT2(&the_key, str);
					}
				else {
					Srfs("bt_key_sprint", 
					"unknown output size");
					return ERROR;
					}
				if (j < bt_key->s_int_vec_len_i())
					strcat(str, ", ");
				}
			strcat(str, ")");
			}
		else {
			Srfs("bt_key_sprint", "unknown type");
			return ERROR;
			}
		if (i < len - 1)
			strcat(str, " ");
		}
	return OK;
}

INT bt_sprint_INT4(BYTE **p_key, BYTE *str)
/* es wird an str angehaengt */
{
	INT ii;
	INT4 i;
	BYTE *p;
	
	bt_get_INT4(p_key, &ii);
	i = (INT4) ii;
	block_swap_bytes((SCHAR *) &i, sizeof(INT4), 1);
	p = (BYTE *) &i;
	for (i = 0; i < 4; i++) {
		sprintf(str + strlen(str), "%ld", (INT)p[i]);
		}
	/* sprintf(str + strlen(str), "%ld", (long) i); */
	return OK;
}

INT bt_sprint_INT2(BYTE **p_key, BYTE *str)
/* es wird an str angehaengt */
{
	INT i;
	
	bt_get_INT2(p_key, &i);
	sprintf(str + strlen(str), "%ld", (long) i);
	return OK;
}

INT bt_key_cmp(BYTE *key1, BYTE *key2, 
	VECTOR_OP V, INT len, INT *res)
{
	BT_KEY_OP bt_key;
	BYTE *the_key1 = key1;
	BYTE *the_key2 = key2;
	INT i, j, type, output_size;
	
	for (i = 0; i < len; i++) {
		bt_key = (BT_KEY_OP) V->s_i(i);
		type = bt_key->s_type_i();
		output_size = bt_key->s_output_size_i();
		if (type == 0) { /* INT */
			if (output_size == 4) {
				bt_cmp_INT4(&the_key1, &the_key2, res);
				if (*res)
					return OK;
				}
			else if (output_size == 2) {
				bt_cmp_INT2(&the_key1, &the_key2, res);
				if (*res)
					return OK;
				}
			else {
				Srfs("bt_key_cmp", 
				"unknown output size");
				return ERROR;
				}
			}
		else if (type == 1) { /* STRING */
			bt_lex_n_cmp(the_key1, the_key2, 
				output_size, res);
			if (*res)
				return OK;
			the_key1 += output_size;
			the_key2 += output_size;
			}
		else if (type == 2) { /* INT VEC */
			for (j = 0; 
				j < bt_key->s_int_vec_len_i(); 
				j++) {
				if (output_size == 4) {
					bt_cmp_INT4(&the_key1, &the_key2, res);
					if (*res)
						return OK;
					}
				else if (output_size == 2) {
					bt_cmp_INT2(&the_key1, &the_key2, res);
					if (*res)
						return OK;
					}
				else {
					Srfs("bt_key_cmp", 
					"unknown output size");
					return ERROR;
					}
				}
			}
		else {
			Srfs("bt_key_cmp", "unknown type");
			return ERROR;
			}
		}
	*res = 0;
	return OK;
}

INT bt_lex_n_cmp(BYTE *p1, BYTE *p2, 
	INT n, INT *res)
{
	UBYTE c1, c2, d1, d2;
	INT l = 0;
	
	while (TRUE) {
		c1 = *(UBYTE *) p1;
		c2 = *(UBYTE *) p2;
		d1 = lex_cmp[c1];
		d2 = lex_cmp[c2];
		if (d1 < d2) {
			/* Beachte: hier auch Abbruch, 
			 * falls nur String 1 zuende. */
			*res = -1;
			return(OK);
			}
		if (d1 > d2) {
			*res = 1;
			return(OK);
			}
		/* Gleichheit */
		if (c1 == 0) {
			/* Abbruchbedingung: 
			 * beide Strings zuende. */
			*res = 0;
			return(OK);
			}
		p1++;
		p2++;
		l++;
		if (l == n) {
			*res = 0;
			return(OK);
			}
		}
}

INT bt_cmp_INT4(BYTE **p_key1, BYTE **p_key2, INT *res)
{
	INT4 int4_1, int4_2;
	INT i;
	BYTE c;
	BYTE *pc_1 = (BYTE *) &int4_1;
	BYTE *pc_2 = (BYTE *) &int4_2;
	
	if (sizeof(INT4) != 4) {
		Srfs("bt_cmp_INT4", 
		"sizeof(INT4) != 4");
		return ERROR;
		}
	for (i = 0; i < 4; i++) {
		c = **p_key1;
		(*p_key1)++;
		pc_1[i] = c;
		c = **p_key2;
		(*p_key2)++;
		pc_2[i] = c;
		}
	if (int4_1 < int4_2) {
		*res = -1;
		return(OK);
		}
	if (int4_1 > int4_2) {
		*res = 1;
		return(OK);
		}
	*res = 0;
	return OK;
}

INT bt_cmp_INT2(BYTE **p_key1, BYTE **p_key2, INT *res)
{
	INT2 int2_1, int2_2;
	INT i;
	BYTE c;
	BYTE *pc_1 = (BYTE *) &int2_1;
	BYTE *pc_2 = (BYTE *) &int2_2;
	
	if (sizeof(INT2) != 2) {
		Srfs("bt_cmp_INT2", 
		"sizeof(INT2) != 2");
		return ERROR;
		}
	for (i = 0; i < 2; i++) {
		c = **p_key1;
		(*p_key1)++;
		pc_1[i] = c;
		c = **p_key2;
		(*p_key2)++;
		pc_2[i] = c;
		}
	if (int2_1 < int2_2) {
		*res = -1;
		return(OK);
		}
	if (int2_1 > int2_2) {
		*res = 1;
		return(OK);
		}
	*res = 0;
	return OK;
}

INT bt_key_fill_in(BYTE *key, VECTOR_OP theOp, VECTOR_OP V, INT len)
{
	BT_KEY_OP bt_key;
	SYM_OP key_op, key_op1;
	VECTOR_OP key_V;
	BYTE *the_key = key;
	INT i, j, type, output_size, fst;
	INTEGER_OB null_ob;
	
	null_ob.m_i(0);
	for (i = 0; i < len; i++) {
		bt_key = (BT_KEY_OP) V->s_i(i);
		type = bt_key->s_type_i();
		output_size = bt_key->s_output_size_i();
		key_op = theOp->s_i( bt_key->s_field1_i() );
		if (type == 0) { /* INT */
			if (output_size == 4) {
				bt_key_INT4_fill_in(&the_key, key_op);
				}
			else if (output_size == 2) {
				bt_key_INT2_fill_in(&the_key, key_op);
				}
			else {
				Srfs("bt_key_fill_in", 
					"unknown output size");
				return ERROR;
				}
			}
		else if (type == 1) { /* STRING */
			bt_key_STR_fill_in(&the_key, output_size, key_op);
			}
		else if (type == 2) { /* INT VEC */
			if (key_op->s_obj_k() != VECTOR) {
				Srfs("bt_key_fill_in", 
					"not a vector");
				return ERROR;
				}
			key_V = (VECTOR_OP) key_op;
			fst = bt_key->s_int_vec_fst_i();
			for (j = 0; 
				j < bt_key->s_int_vec_len_i(); 
				j++) {
				if (fst + j < key_V->s_li())
					key_op1 = key_V->s_i(fst + j);
				else 
					key_op1 = &null_ob;
				if (output_size == 4) {
					bt_key_INT4_fill_in(&the_key, key_op1);
					}
				else if (output_size == 2) {
					bt_key_INT2_fill_in(&the_key, key_op1);
					}
				else {
					Srfs("bt_key_fill_in", 
					"unknown output size");
					return ERROR;
					}
				}
			}
		else {
			Srfs("bt_key_fill_in", "unknown type");
			return ERROR;
			}
		}
	return OK;
}

INT bt_key_INT4_fill_in(BYTE **p_key, SYM_OP key_op)
{
	INTEGER_OP key_op_int;
	INT erg = OK;

	if (key_op->s_obj_k() != INTEGER) {
		Srfs("bt_key_INT4_fill_in", 
			"key_op->s_obj_k() != INTEGER");
		return ERROR;
		}
	key_op_int = (INTEGER_OP) key_op;
	erg += bt_key_INT4_fill_in_(p_key, key_op_int->s_i());
	return erg;
}

INT bt_key_INT4_fill_in_(BYTE **p_key, INT value)
{
	INT4 int4;
	INT i;
	BYTE c;
	BYTE *pc = (BYTE *) &int4;
	
	if (sizeof(INT4) != 4) {
		Srfs("bt_key_INT4_fill_in_", 
			"sizeof(INT4) != 4");
		return ERROR;
		}
	int4 = (INT4) value;
	for (i = 0; i < 4; i++) {
		c = pc[i];
		**p_key = c;
		(*p_key)++;
		}
	return OK;
}

INT bt_key_INT2_fill_in(BYTE **p_key, SYM_OP key_op)
{
	INTEGER_OP key_op_int;
	INT erg = OK;

	if (key_op->s_obj_k() != INTEGER) {
		Srfs("bt_key_INT2_fill_in", 
			"key_op->s_obj_k() != INTEGER");
		return ERROR;
		}
	key_op_int = (INTEGER_OP) key_op;
	erg += bt_key_INT2_fill_in_(p_key, key_op_int->s_i());
	return erg;
}

INT bt_key_INT2_fill_in_(BYTE **p_key, INT value)
{
	INT2 int2;
	INT i;
	BYTE c;
	BYTE *pc = (BYTE *) &int2;
	
	if (sizeof(INT2) != 2) {
		Srfs("bt_key_INT2_fill_in_", 
			"sizeof(INT2) != 2");
		return ERROR;
		}
	int2 = (INT2) value;
	for (i = 0; i < 2; i++) {
		c = pc[i];
		**p_key = c;
		(*p_key)++;
		}
	return OK;
}

INT bt_key_STR_fill_in(BYTE **p_key, INT output_size, SYM_OP key_op)
{
	STRING_OP key_op_str;
	INT i, k;
	BYTE c;
	BYTE *pc;
	
	if (key_op->s_obj_k() != STRING) {
		Srfs("bt_key_STR_fill_in", 
			"key_op->s_obj_k() != STRING");
		return ERROR;
		}
	key_op_str = (STRING_OP) key_op;
	pc = key_op_str->s_str();
	k = strlen(pc);
	for (i = 0; i < output_size; i++) {
		if (i < k)
			c = pc[i];
		else
			c = 0;
		**p_key = c;
		(*p_key)++;
		}
	return OK;
}

INT bt_key_change_INT(BYTE *key, 
	VECTOR_OP V, INT idx, INT new_value)
/* idx ten Eintrag in key auf new_value setzen; 
 * es muss sich um ein 
 * INT2 oder INT4 Schluesselfeld handeln
 * (nicht fuer INT VEC Schluessel geeignet). */
{
	BT_KEY_OP bt_key;
	BYTE *the_key = key;
	INT i, j, type, output_size, len;
	
	for (i = 0; i <= idx; i++) {
		bt_key = (BT_KEY_OP) V->s_i(i);
		type = bt_key->s_type_i();
		output_size = bt_key->s_output_size_i();
		if (type == 0) { /* INT */
			if (i == idx) {
				if (output_size == 4) {
					bt_key_INT4_fill_in_(&the_key, new_value);
					}
				else if (output_size == 2) {
					bt_key_INT2_fill_in_(&the_key, new_value);
					}
				else {
					Srfs("bt_key_change_INT", 
					"unknown output size");
					return ERROR;
					}
				return OK;
				}
			else {
				if (output_size == 4) {
					the_key += 4;
					}
				else if (output_size == 2) {
					the_key += 2;
					}
				else {
					Srfs("bt_key_change_INT", 
					"unknown output size");
					return ERROR;
					}
				}
			}
		else if (type == 1) { /* STRING */
			if (i == idx) {
				Srfs("bt_key_change_INT", 
				"type == 1");
				return ERROR;
				}
			the_key += output_size;
			}
		else if (type == 2) { /* INT VEC */
			if (i == idx) {
				Srfs("bt_key_change_INT", 
				"type == 2");
				return ERROR;
				}
			len = bt_key->s_int_vec_len_i();
			for (j = 0; j < len; j++) {
				if (output_size == 4) {
					the_key += 4;
					}
				else if (output_size == 2) {
					the_key += 2;
					}
				else {
					Srfs("bt_key_change_INT", 
					"unknown output size");
					return ERROR;
					}
				}
			}
		else {
			Srfs("bt_key_change_INT", "unknown type");
			return ERROR;
			}
		}
	return OK;
}

INT bt_get_INT4(BYTE **p_key, INT *p_i)
{
	INT4 int4;
	INT i;
	BYTE c;
	BYTE *pc = (BYTE *) &int4;
	
	if (sizeof(INT4) != 4) {
		Srfs("bt_get_INT4", 
			"sizeof(INT4) != 4");
		return ERROR;
		}
	for (i = 0; i < 4; i++) {
		c = **p_key;
		(*p_key)++;
		pc[i] = c;
		}
	*p_i = (INT) int4;
	return OK;
}

INT bt_get_INT2(BYTE **p_key, INT *p_i)
{
	INT2 int2;
	INT i;
	BYTE c;
	BYTE *pc = (BYTE *) &int2;
	
	if (sizeof(INT2) != 2) {
		Srfs("bt_get_INT2", 
			"sizeof(INT2) != 2");
		return ERROR;
		}
	for (i = 0; i < 2; i++) {
		c = **p_key;
		(*p_key)++;
		pc[i] = c;
		}
	*p_i = (INT) int2;
	return OK;
}

#endif /* DB_TRUE */

