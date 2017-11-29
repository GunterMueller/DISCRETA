/* divs.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef DIVS_TRUE

#include <DISCRETA/divs.h>

/* 
 * BITVEC
 */

INT BITVEC_OB::field_name(INT i, INT j, BYTE *str)
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "ints"; break;
	default:
		return error("BITVEC::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

INT BITVEC_OB::init_ascii(BYTE *s)
{
	UINT4 *pi4;
	UINT w;
	INT i, ii, l, nb_ints;
	
	
	m_il(1);
	c_obj_k(BITVEC_KIND);
	l = strlen(s);
	nb_ints = (l >> 5) + 1;
	s_ints()->m_il_n(nb_ints);
	pi4 = (UINT4 *) my_malloc(nb_ints * sizeof(INT4), "BITVEC_OB::init_ascii");
	for (i = 0; i < nb_ints; i++) {
		pi4[i] = 0;
		}
	ii = 0;
	w = 0;
	for (i = 0; i < l; i++) {
		w <<= 1;
		if (s[i] == '1')
			w |= (UINT) 1;
		if ((i + 1) % 32 == 0) {
			pi4[ii] = w;
			ii++;
			w = 0;
			}
		}
	block_swap_bytes((SCHAR *) pi4, sizeof(INT4), nb_ints);
	for (i = 0; i < nb_ints; i++) {
		w = pi4[i];
		s_ints_i(i)->m_i(w);
		}
	my_free(pi4);
	return OK;
}

INT BITVEC_OB::print()
{
	BYTE str[10000];
	
	str[0] = 0;
	sprint(str);
	printf("%s\n", str);
	return OK;
}

INT BITVEC_OB::sprint(BYTE *s)
{
	UINT4 *pi4;
	UINT w;
	INT i, ii, l, nb_ints, bit;
	
	
	nb_ints = s_ints()->s_li();
	l = nb_ints << 5;
	pi4 = (UINT4 *) my_malloc(nb_ints * sizeof(INT4), "BITVEC_OB::print");
	for (i = 0; i < nb_ints; i++) {
		pi4[i] = s_ints_ii(i);
		}
	block_swap_bytes((SCHAR *) pi4, sizeof(INT4), nb_ints);
	ii = 0;
	w = pi4[ii];
	for (i = 0; i < l; i++) {
		bit = w & 1;
		if (bit)
			sprintf(s + strlen(s), "1");
		else
			sprintf(s + strlen(s), "0");
		w >>= 1;
		if ((i + 1) % 32 == 0) {
			w = pi4[++ii];
			}
		}
	my_free(pi4);
	return OK;
}

/* 
 * STRING
 */

INT string_test()
{
	STRING_OB s1, s2;
	
	s1.init("this is a string");
	s1.println();
	s1.copy(&s2);
	s2.println();
	printf("vor s1.freeself():\n");
	s1.freeself();
	printf("nach s1.freeself():\n");
	return OK;
}

INT STRING_OB::freeself()
{
	if (s_obj_k() != EMPTY) {
		if (s_obj_k() != STRING)
			return error("s_freeself(): kind != STRING");
		if (ob_self.ob_charpointer)
			tstrcls(&ob_self.ob_charpointer);
		}
	ob_self.ob_charpointer = NIL;
	c_obj_k(EMPTY);
	return OK;
}

INT STRING_OB::init(BYTE *str)
{
	freeself();
	if (str && strlen(str)) {
		if (!tstropn(&ob_self.ob_charpointer, str))
			return ERROR;
		}
	c_obj_k(STRING);
	return OK;
}

BYTE *STRING_OB::s_str()
{
	if (s_obj_k() != STRING) {
		error("s_str()|not a STRING object");
		return NULL;
		}
	return(tstr(ob_self.ob_charpointer));
}

INT STRING_OB::c_str(BYTE *str)
/* overwrites an existing string; 
 * does not copy str. */
{
	if (s_obj_k() != STRING) {
		error("c_str()|not a STRING object");
		return ERROR;
		}
	ob_self.ob_charpointer = str;
	return OK;
}

INT STRING_OB::copy(STRING_OP b)
{
	INT ret = OK;
	
	if (b == this)
		return OK;
	b->freeself();
	ret += b->init(ob_self.ob_charpointer);
	return ret;
}

INT STRING_OB::append(STRING_OP b)
/* concatenates strings: 
 * first this then b, result into this. */
/* warning: former routinue string_append(): 
 * first b then a, result into b */
{
	INT ret = OK;
	
	if (s_obj_k() != STRING)
		return error("append(): a is not a string");
	if (b->s_obj_k() != STRING)
		return error("append(): b is not a string");
	if (b == this)
		return error("append()|b == this");
	if (!tstrcat(&ob_self.ob_charpointer, 
		b->s_str()))
		ret = ERROR;
	return ret;
}

INT STRING_OB::sprint(BYTE *str)
/* appends to str. 
 * writes to maximal strlength of 200. */
{
	INT i, j;
	BYTE *s;
	
	s = s_str();
	j = strlen(str);
	i = strlen(s);
	if (i + j > 200) {
		strncpy(str + j, s, 200 - j);
		str[200] = 0;
		return OK;
		}
	strcpy(str + j, s);
	return OK;
}

INT STRING_OB::sscan(BYTE *str)
{
	return init(str);
}

/* 
 * MEM
 */

#define MEM_OVERSIZE 32
#define MEM_OVERSIZE1 512

/*
 * INT - offset - 3 + 0: alloc_length
 *              - 3 + 1: used_length
 *              - 3 + 2: cur_pointer
 * 
 * alloc_length: allocated length in BYTEs
 * used_length: used length in BYTES
 * cur_pointer:
 *          0 <= pointer < used length. 
 */

INT MEM_OB::freeself()
{
	BYTE *pc;
	INT *pi;
	
	if (s_obj_k() != EMPTY) {
		if (s_obj_k() != MEM)
			return error("MEM::freeself(): kind != MEM");
		pc = ob_self.ob_charpointer;
		if (pc) {
			pi = (INT *)pc;
			pi -= 3;
			my_free(pi);
			}
		}
	ob_self.ob_charpointer = NIL;
	c_obj_k(EMPTY);
	return OK;
}

INT MEM_OB::s_alloc_length_i()
{
	BYTE *pc;
	INT *pi;
	
	pc = ob_self.ob_charpointer;
	pi = (INT *)pc;
	pi -= 3;
	return pi[0];
}

INT MEM_OB::s_used_length_i()
{
	BYTE *pc;
	INT *pi;
	
	pc = ob_self.ob_charpointer;
	pi = (INT *)pc;
	pi -= 3;
	return pi[1];
}

INT MEM_OB::s_cur_pointer_i()
{
	BYTE *pc;
	INT *pi;
	
	pc = ob_self.ob_charpointer;
	pi = (INT *)pc;
	pi -= 3;
	return pi[2];
}

INT *MEM_OB::s_alloc_length_pi()
{
	BYTE *pc;
	INT *pi;
	
	pc = ob_self.ob_charpointer;
	pi = (INT *)pc;
	pi -= 3;
	return pi;
}

INT *MEM_OB::s_used_length_pi()
{
	BYTE *pc;
	INT *pi;
	
	pc = ob_self.ob_charpointer;
	pi = (INT *)pc;
	pi -= 3;
	return pi + 1;
}

INT *MEM_OB::s_cur_pointer_pi()
{
	BYTE *pc;
	INT *pi;
	
	pc = ob_self.ob_charpointer;
	pi = (INT *)pc;
	pi -= 3;
	return pi + 2;
}

INT MEM_OB::c_alloc_length(INT alloc_length)
{
	BYTE *pc;
	INT *pi;
	
	pc = ob_self.ob_charpointer;
	pi = (INT *)pc;
	pi -= 3;
	pi[0] = alloc_length;
	return OK;
}

INT MEM_OB::c_used_length(INT used_length)
{
	BYTE *pc;
	INT *pi;
	
	pc = ob_self.ob_charpointer;
	pi = (INT *)pc;
	pi -= 3;
	pi[1] = used_length;
	return OK;
}

INT MEM_OB::c_cur_pointer(INT cur_pointer)
{
	BYTE *pc;
	INT *pi;
	
	pc = ob_self.ob_charpointer;
	pi = (INT *)pc;
	pi -= 3;
	pi[2] = cur_pointer;
	return OK;
}

INT MEM_OB::init(INT length, BYTE *d)
{
	INT erg = OK, i;
	BYTE *pc;
	
	erg += alloc(length);
	pc = ob_self.ob_charpointer;
	for (i = 0; i < length; i++) {
		pc[i] = d[i];
		}
	return erg;
}

static INT code(UBYTE *pc, INT l, UBYTE *pc2, UBYTE code_char);
static INT decode(UBYTE *pc2, INT l2, UBYTE *pc, UBYTE code_char);

INT MEM_OB::compress(INT f_verbose)
/* Wolfgang Boessenecker 9/94 */
{
	MEM_OB mem2;
	BYTE *pc, *pc2;
	INT l, l2, l_c;

	pc = ob_self.ob_charpointer;
	l = s_used_length_i();
	if (f_verbose) {
		printf("compressing from %ld to ", l);
		fflush(stdout);
		}
	l_c = mult_char_c((BYTE) 0);
	l2 = l - l_c + ((l + 7) >> 3);
	mem2.alloc(l2); /* sets used_length to l2 */
	pc2 = mem2.ob_self.ob_charpointer;
	code((UBYTE *) pc, l, (UBYTE *) pc2, (UBYTE) 0);
#if 0
	if (l3 != l2) {
		printf("MEM::compress() warning: l2 = %ld != l3 = %ld\n", l2, l3);
		fflush(stdout);
		}
#endif
	mem2.swap(this);
	if (f_verbose) {
		printf("%ld Bytes.\n", l2);
		fflush(stdout);
		}
	return OK;
}

INT MEM_OB::decompress(INT f_verbose)
{
	MEM_OB mem;
	BYTE *pc, *pc2;
	INT l, l2;
	
	pc2 = ob_self.ob_charpointer;
	l2 = s_used_length_i();
	if (f_verbose) {
		printf("decompressing from %ld to ", l2);
		fflush(stdout);
		}
	l = decode((UBYTE *) pc2, l2, NIL, (UBYTE) 0);
	mem.alloc(l);
	pc = mem.ob_self.ob_charpointer;
	decode((UBYTE *) pc2, l2, (UBYTE *) pc, (UBYTE) 0);
	mem.swap(this);
	if (f_verbose) {
		printf("%ld Bytes.\n", l);
		fflush(stdout);
		}
	return OK;
}

static INT code(UBYTE *pc, INT l, UBYTE *pc2, UBYTE code_char)
/* WB 940919 */
{
	UBYTE cc;
	INT pos = 0, pos2 = 0, pos2h = 0, i;

	while (pos < l) {
		pos2++;
		cc = 0;
#if 0
		if ((posf % 100000) == 0) {
			printf("%ld\n", posf);
			fflush(stdout);
			}
#endif
		for (i = 0; i < 8; i++) {
			cc <<= 1;
			if (pos < l) {
				if (pc[pos] == code_char)
					cc = cc | 0X1U;
				else {
					pc2[pos2] = pc[pos];
					pos2++;
					}
				pos++;
				}
			}
		pc2[pos2h] = cc;
		pos2h = pos2;
		}
	return OK;
}

static INT decode(UBYTE *pc2, INT l2, UBYTE *pc, UBYTE code_char)
/* returns length of decompressed data 
 * pc may be NIL */
{
	UBYTE cc = 0;
	INT pos = 0, pos2 = 0, i = 8;
	
	while (TRUE) {
	/* for (; pos2 < l2; ) { */
		if (pos2 >= l2 && i >= 8)
			break;
		if (i == 8) {
			cc = pc2[pos2];
			pos2++;
			i = 0;
			}
		if (cc & (UBYTE) 128U) {
			if (pc) {
				pc[pos] = code_char;
				}
			pos++;
			}
		else {
			if (pos2 < l2) {
				if (pc) {
					pc[pos] = pc2[pos2];
					}
				pos2++;
				pos++;
				}
			}
		cc <<= 1;
		i++;
		}
	return pos;
}

INT MEM_OB::mult_char_c(BYTE c)
{
	BYTE *pc;
	INT i, l = 0, len;
	
	pc = ob_self.ob_charpointer;
	len = s_used_length_i();
	for (i = 0; i < len; i++)
		if (pc[i] == c)
			l++;
	return l;
}

INT MEM_OB::alloc(INT length)
/* setzt alloc_length auf 
 *  length + MEM_OVERSIZE, 
 *       used_length auf length, 
 *       cur_pointer auf 0.
 */
{
	INT size, mem_oversize;
	INT *pi;
	BYTE *pc;
	
	if (length >= MEM_OVERSIZE) {
		mem_oversize = MEM_OVERSIZE1;
		}
	else {
		mem_oversize = MEM_OVERSIZE;
		}
	freeself();
	size = length + mem_oversize + 
		3 * sizeof(INT);
	pi = (INT *) my_malloc(size, "MEM::alloc");
	if (pi == NIL) {
		Srfs("MEM::alloc", "no memory");
		return ERROR;
		}
	pi[0] = length + mem_oversize;
	pi[1] = length;
	pi[2] = 0;
	pc = (BYTE *)(pi + 3);
	ob_self.ob_charpointer = pc;
	c_obj_k(MEM);
	return OK;
}

INT MEM_OB::copy(MEM_OP b)
{
	INT ret = OK;
	INT length;
	
	if (b == this)
		return OK;
	b->freeself();
	length = s_used_length_i();
	ret += b->init(length, ob_self.ob_charpointer);
	b->c_cur_pointer(s_cur_pointer_i());
	return ret;
}

INT MEM_OB::append(INT length, BYTE *d)
{
	BYTE *pc;
	INT i, old_length, new_length;
	INT erg = OK;
	
	old_length = s_used_length_i();
	new_length = old_length + length;
	if (new_length > s_alloc_length_i()) {
		erg += realloc(new_length);
		}
	else {
		c_used_length(new_length);
		}
	pc = ob_self.ob_charpointer;
	for (i = 0; i < length; i++) {
		pc[old_length + i] = d[i];
		}
	return erg;
}

INT MEM_OB::realloc(INT new_length)
{
	INT old_length;
	INT old_cur_pointer;
	INT erg = OK, i;
	BYTE *old_pc, *pc;
	INT *old_pi;
	
	old_pc = ob_self.ob_charpointer;
	old_pi = (INT *)old_pc - 3;
	old_length = s_used_length_i();
	old_cur_pointer = s_cur_pointer_i();
	if (new_length < old_length)
		printf("MEM::realloc()|"
			"warning: new_length < old_length\n");
	ob_self.ob_charpointer = NIL;
	erg += alloc(new_length);
	pc = ob_self.ob_charpointer;
	for (i = 0; 
		i < MINIMUM(old_length, new_length); 
		i++) {
		pc[i] = old_pc[i];
		}
	for (i = old_length; i < new_length; i++) {
		pc[i] = 0;
		}
	my_free(old_pi);
	c_cur_pointer(old_cur_pointer);
#ifdef DEBUG_MEM
	printf("MEM::realloc()|%ld\n", s_used_length_i());
#endif
	return erg;
}

INT MEM_OB::sprint(BYTE *str)
/* appends to str. 
 * writes to maximal strlength of 200. */
{
	BYTE s[256];
	
	sprintf(s, "MEM: alloc_length = %ld "
		"used_length = %ld cur_pointer = %ld", 
		s_alloc_length_i(), s_used_length_i(), 
		s_cur_pointer_i());
	if (strlen(s) + strlen(str) < 200)
		strcat(str, s);
	return OK;
}

INT MEM_OB::write_char(BYTE c)
{	
	return append(1, &c);
}

INT MEM_OB::read_char(BYTE *c)
{
	INT l1, cur_p, used_length;
	INT erg = OK;
	BYTE *cp;
	
	cur_p = s_cur_pointer_i();
	used_length = s_used_length_i();
	l1 = used_length - cur_p;
	if (1 > l1) {
		Srfs("MEM::read_char", "1 > l1");
		return ERROR;
		}
	cp = ob_self.ob_charpointer + cur_p;
	*c = *cp;
	c_cur_pointer(cur_p + 1);
	return erg;
}

INT MEM_OB::write_int(INT i)
{
	INT4 i1 = (INT4) i;
	INT erg = OK;
	
#ifdef DEBUG_MEM
	printf("MEM::write_int()|%ld\n", 
		(INT) i1);
#endif
	block_swap_bytes((SCHAR *)&i1, 
		sizeof(INT4), 1);
	erg += append(sizeof(INT4), 
		(BYTE *) &i1);
	return erg;
}

INT MEM_OB::read_int(INT *i)
{
	INT4 i1;
	INT l1, j, cur_p, used_length;
	INT erg = OK;
	BYTE *cp, *cp1;
	
	cur_p = s_cur_pointer_i();
	used_length = s_used_length_i();
	l1 = used_length - cur_p;
	if ((INT) sizeof(INT4) > l1) {
		return error("MEM::read_int() sizeof(INT4) > l1");
		}
	cp = ob_self.ob_charpointer + cur_p;
	cp1 = (BYTE *) &i1;
	for (j = 0; j < (INT) sizeof(INT4); j++) {
		*cp1 = *cp;
		cp1++;
		cp++;
		}
	/* i1 = *(INT *) (cp + cur_p); */
	c_cur_pointer(cur_p + sizeof(INT4));
	block_swap_bytes((SCHAR *) &i1, 
		sizeof(INT4), 1);
#ifdef DEBUG_MEM
	printf("MEM::read_int()|%ld\n", (long) i1);
#endif
	*i = (INT) i1;
	return erg;
}

INT MEM_OB::length_in_entries(INT entry_size)
{
	INT len;
	
	len = s_used_length_i() / entry_size;
	if (len * entry_size != s_used_length_i()) {
		Srfs("MEM::length_in_entries", 
			"len * entry_size != s_used_length_i");
		return 0;
		}
	return len;
}

BYTE *MEM_OB::pc_ith(INT i, INT entry_size)
{
	BYTE *pc;
	
	pc = ob_self.ob_charpointer;
	return pc + i * entry_size;
}

INT MEM_OB::insert_at(INT entry_size, BYTE *d, INT i)
{
	INT len, new_length;
	INT j;
	BYTE *pc;
	
	len = length_in_entries(entry_size);
	if (i > len) {
		Srfs("MEM::insert_at", "i > len");
		return ERROR;
		}
	if ((len + 1) * entry_size > s_alloc_length_i()) {
		new_length = (len + 1) * entry_size;
		if (realloc(new_length) != OK) {
			Srff("MEM::insert_at", "realloc");
			return ERROR;
			}
		}
	pc = ob_self.ob_charpointer;
	for (j = (len + 1) * entry_size - 1; 
		j >= (i + 1) * entry_size; j--)
		pc[j] = pc[j - entry_size];
	for (j = 0; j < entry_size; j++) {
		pc[i * entry_size + j] = d[j];
		}
	c_used_length((len + 1) * entry_size);
	return OK;
}

INT MEM_OB::isd(INT entry_size, 
	BYTE *add_this, INT *f_found, 
	INT (*cmp_func_data)(void *a, void *b, void *data), 
	void *data)
{
	INT idx;
	
	if (sd(entry_size, add_this, &idx, f_found, 
		cmp_func_data, data) != OK) {
		Srff("MEM_OB::isd", "sd");
		return ERROR;
		}
	if (*f_found) {
		return OK;
		}
	if (insert_at(entry_size, add_this, idx) != OK) {
		Srff("MEM::isd", "insert_at");
		return ERROR;
		}
	return OK;
}

INT MEM_OB::sd(INT entry_size, 
	void *search_this, 
	INT *idx, INT *f_found, 
	INT (*cmp_func_data)(
		void *a, void *b, void *data), 
	void *data)
{
	INT len;
	INT l, r, m, res;
	BYTE *pc;
	
	len = length_in_entries(entry_size);
	pc = ob_self.ob_charpointer;
	if (len == 0) {
		*idx = 0;
		*f_found = FALSE;
		return OK;
		}
	l = 0;
	r = len;
	*f_found = FALSE;
	/* Invariante:
	 * V[i] <= search_this fuer i < l;
	 * V[i] >  search_this fuer i >= r;
	 * r - l ist Laenge des Suchbereichs. 
	 */
	while (l < r) {
		m = (l + r) >> 1;
		/* Bei einem Suchbereich 
		 * gerader Laenge wird das Element 
		 * oberhalb der Mitte genommen ! */
		res = (*cmp_func_data)(
			pc + m * entry_size, search_this, data);
		if (res <= 0) {
			l = m + 1;
			if (res == 0)
				*f_found = TRUE;
			}
		else
			r = m;
		}
	/* Jetzt: l == r; 
	 * f_found ist evtl. schon gesetzt worden */
	*idx = l;
	return OK;
}

INT MEM_OB::print_data(INT entry_size, 
	INT (*print_data)(void *a, void *data), 
	void *data)
{
	INT i, len;
	BYTE *pc;
	
	len = s_used_length_i() / entry_size;
	if (len * entry_size != s_used_length_i()) {
		Srfs("MEM::print_data", 
			"len * entry_size != s_used_length_i");
		return ERROR;
		}
	pc = ob_self.ob_charpointer;
	for (i = 0; i < len; i++) {
		printf("%ld: ", i);
		(*print_data)(pc + i * entry_size, data);
		}
	return OK;
}

INT MEM_OB::slurp_file(BYTE *fname, INT f_v, FILE *fp_txt)
{
	FILE *fp;
	INT fsize;
	BYTE *pc;

	fsize = file_size(fname);
	alloc(fsize);
	pc = ob_self.ob_charpointer;
	fp = fopen(fname, "r");
	if ((INT) fread(pc, 1 /* size */, fsize /* nitems */, fp) != fsize)
		return error("MEM::surp_file() error in fread");
	fclose(fp);
	c_used_length(fsize);
	c_cur_pointer(0);
	if (f_v) {
		fprintf(fp_txt, "MEM::slurp_file()|"
			"read file %s of size %ld\n", 
			fname, file_size(fname));
		fflush(fp_txt);
		}
	return OK;
}

INT MEM_OB::write_to_file(BYTE *fname, INT f_v, FILE *fp_txt)
{
	FILE *fp;
	INT size;
	BYTE *pc;

	size = s_used_length_i();
	pc = ob_self.ob_charpointer;
	
	fp = fopen(fname, "wb");

	fwrite(pc, 1 /* size */, size /* items */, fp);
	
	fclose(fp);
	if (file_size(fname) != size)
		return error("write_op_file(): file_size(fname) != size");
	if (f_v) {
		fprintf(fp_txt, "MEM::write_to_file()|"
			"wrote file %s of size %ld\n", 
			fname, file_size(fname));
		fflush(fp_txt);
		}
	return OK;
}

/*
 * CONTI
 */

INT CONTI_OB::field_name(INT i, INT j, BYTE *str)
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "nb_L"; break;
	case 1: s = "L"; break;
	case 2: s = "nb_C"; break;
	case 3: s = "C"; break;
	case 4: s = "CS"; break;
	default:
		return error("CONTI::field_name()|"
			"i out of range");
	}
	strcpy(str, s);
	return OK;
}

INT CONTI_OB::init()
{
	INT erg = OK;
	
	erg += m_il(5);
	c_obj_k(CONTI);
	s_nb_L()->m_i(0);
	s_nb_C()->m_i(0);
	erg += s_L()->m_il(VECTOR_OVERSIZE);
	erg += s_C()->m_il(VECTOR_OVERSIZE);
	erg += s_CS()->m_il(VECTOR_OVERSIZE);
	return erg;
}

INT CONTI_OB::sprint(BYTE *s)
{
	BYTE str[512];
	
	sprintf(str, "CONTI: nb_L = %ld nb_C = %ld ", 
		s_nb_L_i(), s_nb_C_i());
	if (strlen(s) + strlen(str) < 200)
		strcat(s, str);
	s_L()->sprint(s);
	return OK;
}

INT CONTI_OB::add_label(BYTE *label)
{
	STRING_OB s;
	INT erg = OK;
	
	erg += s.init(label);
	erg += s_L()->append_element(s_nb_L(), &s);
	return erg;
}

INT CONTI_OB::add_object(
	BYTE *label, SYM_OP op)
{
	STRING_OB s;
	INT erg = OK;
	
	erg += s.init(label);
	erg += s_CS()->append_element(s_nb_C(), &s);
	s_nb_C()->dec();
	erg += s_C()->append_element(s_nb_C(), op);
	return erg;
}

INT CONTI_OB::find_object(
	BYTE *label, INT *idx, INT *f_found)
{
	INT l, len;
	
	len = s_nb_C_i();
	for (l = 0; l < len; l++) {
		if (strcmp(s_CS_i_str(l), label) == 0) {
			*idx = l;
			*f_found = TRUE;
			return OK;
			}
		}
	*f_found = FALSE;
	return OK;
}

INT ascii_formatter(VECTOR_OP text, VECTOR_OP in, INT length)
{
	STRING_OB s1;
	STRING_OP s2;
	INT i, l;
	INT l1, l2;

	l = in->s_li();
	s1.init("");
	/* text->m_li(0); */
	for (i = 0; i < l; i++) {
		s2 = (STRING_OP) in->s_i(i);
		if (s2->s_obj_k() == EMPTY) { /* empty object means new line ! */
			text->inc();
			s1.swap(text->s_i(text->s_li() - 1));
			s1.init("");
			continue;
			}
		l2 = strlen(s2->s_str());
		l1 = strlen(s1.s_str());
		if (l1 + l2 > length) { /* flush s1 */
			text->inc();
			s1.swap(text->s_i(text->s_li() - 1));
			s1.init("");
			}
		s1.append(s2);
		}
	if (strlen(s1.s_str())) {
		text->inc();
		s1.swap(text->s_i(text->s_li() - 1));
		}
	return OK;
}

INT ascii_formatter_paragraph(VECTOR_OP text, BYTE *heading, INT indent, 
	VECTOR_OP in, INT length)
{
	STRING_OB s1, ss2, s3;
	STRING_OP s2;
	BYTE *p2;
	INT i, l, wrap;
	INT l1, l2;
	BYTE ind[1024], c;

	for (i = 0; i < indent; i++)
		ind[i] = ' ';
	ind[indent] = 0;
	l = in->s_li();
	s1.init(heading);
	/* text->m_li(0); */
	for (i = 0; i < l; i++) {
		s2 = (STRING_OP) in->s_i(i);
		ss2.init(s2->s_str());
#if 0
		if (s2->s_obj_k() == EMPTY) { /* empty object means new line ! */
			text->inc();
			s1.swap(text->s_i(text->s_li() - 1));
			s1.init("");
			continue;
			}
#endif
		p2 = ss2.s_str();
		l2 = strlen(p2);
		l1 = strlen(s1.s_str());
		if (l1 + l2 > length) { /* flush s1 */
			wrap = find_wrap_position(p2, length - l1);
			c = p2[wrap];
			p2[wrap] = 0;
			s3.init(p2);
			s1.append(&s3);
			p2[wrap] = c;
			while (p2[wrap] == ' ')
				wrap++;
			ss2.init(&p2[wrap]);
			text->inc();
			s1.swap(text->s_i(text->s_li() - 1));
			s1.init(ind);
			}
		s1.append(&ss2);
		}
	if (strlen(s1.s_str())) {
		text->inc();
		s1.swap(text->s_i(text->s_li() - 1));
		}
	return OK;
}

int find_wrap_position(BYTE *p, INT max)
{
	INT i, ii;
	BYTE c;
	
	for (i = max; i > 0; i--) {
		c = p[i - 1];
		if (c == ',' || c == ':' || c == '>' || c == '=') {
		
			/* have to test if it is not an inner ',': 
			 * inside (  ) or <   > */
			if (c == ',' && i - 1 > 0) {
				for (ii = i - 2; ii >= 0; ii--) {
					if (p[ii] == ')' || p[ii] == '>')
						return i;
					if (p[ii] == '<' || p[ii] == '(')
						break;
					}
				if (ii >= 0)
					continue; /* next i */
				}
			
			return i;
			}
		}
	return 0;
}

#endif /* DIVS_TRUE */

