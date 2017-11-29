/* vec.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef VECTORTRUE

#include <DISCRETA/ma.h>
#include <DISCRETA/perm.h>
#include <DISCRETA/lb.h>
#include <DISCRETA/geo.h>
#include <DISCRETA/poly.h>

/*
 * VECTOR
 */

#if TEXDOCU
SYM_OP VECTOR_OB::s_i(INT i)
#endif
{ 
#ifdef VECTOR_SAFE
		if (s_obj_bk() != VECTOR) {
			printobjectkind();
			error("VECTOR access violation: vector obj_bk != VECTOR");
			}
		if (s_s() == NULL)
			error("VECTOR access violation: vector not allocated");
		if (i < 0 || i >= s_li()) {
			printf("i = %ld, length = %ld\n", i, s_li()); fflush(stdout);
			error("VECTOR access violation.");
			}
#endif
		return(s_s() + i);
}

#if TEXDOCU
INT VECTOR_OB::Print()
#else
prints the vector elements each on a separate line.
#endif
{
	INT i, l;

	l = s_li();
	for (i = 0; i < l; i++) {
		s_i(i)->println();
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::latex(FILE *fp)
#else
Prints the vector in latex format to FILE fp.
#endif
{
	INT i, l;
	
	l = s_li();
	fprintf(fp, "\\left(\n");
	for (i = 0; i < l; i++) {
		s_i(i)->latex(fp);
		if (i < l - 1)
			fprintf(fp, ", ");
		}
	fprintf(fp, "\\right)\n");
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::Latex(FILE *fp)
#else
latexs all vector elements each on a line. 
The lines are set into math mode.
They are terminated by double backslash.
#endif
{
	INT i, l;

	l = s_li();
	for (i = 0; i < l; i++) {
		fprintf(fp, "$");
		s_i(i)->latex(fp);
		fprintf(fp, "$ \\\\\n");
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::fprint_GAP(FILE *fp)
#endif
{
	BYTE str[10000];
	INT i, l;

	fprintf(fp, "[");
	l = s_li();
	for (i = 0; i < l; i++) {
		str[0] = 0;
		s_i(i)->sprint(str);
		s_i(i)->fprint_GAP(fp);
		if (i < l - 1) {
			fprintf(fp, ", ");
			if (strlen(str) > 30)
				fprintf(fp, "\n");
			}
		}
	fprintf(fp, "]\n");
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::sprint(BYTE *str)
#endif
{
	sprint_len(s_li(), str);
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::sprint_len(INT len, BYTE *str)
#endif
{
	INT i;
	BYTE str1[1024];
	
	strcat(str, "(");
	for (i = 0; i < len; i++) {
			str1[0] = 0;
			s_i(i)->sprint(str1);
			if (strlen(str) + strlen(str1) < 200)
				strcat(str, str1);
			else
				return OK;
		if (i < s_li() - 1)
			strcat(str, ", ");
		}
	strcat(str, ")");
	return(OK);
}

#if TEXDOCU
INT vec_test()
#endif
{
	VECTOR_OB V1, V2;
	INTEGER_OB V1len, V2len, z;
	INT i, j, idx, f_found;
	
	V1.m_il(VECTOR_OVERSIZE);
	V2.m_il(VECTOR_OVERSIZE);
	V1len.m_i(0);
	V2len.m_i(0);
	j = 7;
	for (i = 0; i < 20; i++) {
		j = (j * j) % 61;
		z.m_i(j);
		V1.append_element(&V1len, &z);
		V2.search(V2len.s_i(), TRUE, &z, &idx, &f_found);
		V2.insert_at(&V2len, idx, &z);
		/* V2.insert_sorted(&V2len, TRUE, &z); */
		}
	V1.println();
	V1.quicksort(V1len.s_i(), TRUE);
	V1.println();
	V2.println();
	if (V1.sym_comp(&V2) != 0) {
		printf("not equal\n");
		return ERROR;
		}
	printf("equal, OK\n");
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::freeself()
#else
frees all vector elements and destroys the vector 
structure. This function is called automatically from the 
freeself function of the SYM\_OB class.
It is directly called if vector objects are destroyed 
for example in a destruction process for automoatic variables 
of a function.
#endif
{
	OBJECTSELF d;
	INT i, erg = OK;
	SYM_OP z;

	if (emptyp())
		return OK;
	if (s_obj_bk() != VECTOR) {
		// printf("VEC::freeself() warning: freeing non-vector object!\n");
		return ((SYM_OP) this)->freeself();
		}
	d = s_obj_s();
	if (d.ob_vector == NIL) {
		c_obj_k(EMPTY);
		/* printf("VEC::freeself(): "
		"d.ob_vector == NIL\n"); fflush(stdout); */
		return OK;
		}
	if (s_l() == NIL) {
		c_obj_k(EMPTY);
		return error("VEC::freeself(): "
		"s_l() == NIL");
		}
	z = s_s();
	if (z == NIL) {
		my_free(d.ob_vector);
		c_obj_k(EMPTY);
		return erg;
		}
	if (s_li() > 0) {
		for (i = 0; i < s_li(); i++, z++) {
			if (! z->emptyp())
				if (z->s_obj_k() != INTEGER)
					erg += z->freeself();
			}
		my_free(s_s());
		}
	freeall(s_l()); 
	my_free(d.ob_vector);
	c_obj_k(EMPTY);
	return erg;
}

#if TEXDOCU
INT VECTOR_OB::freeself_debug(INT print_depth)
#else
frees all vector elements and destroys the vector 
structure. This function is called automatically from the 
freeself function of the SYM\_OB class.
#endif
{
	OBJECTSELF d;
	INT i, erg = OK, pd1 = 0;
	SYM_OP z;

	if (print_depth > 0)
		pd1 = print_depth - 1;
	if (emptyp())
		return OK;
	if (s_obj_bk() != VECTOR) {
		printf("VEC::freeself() warning: freeing non-vector object!\n");
		return ((SYM_OP) this)->freeself();
		}
	d = s_obj_s();
	if (d.ob_vector == NIL) {
		c_obj_k(EMPTY);
		/* printf("VEC::freeself(): "
		"d.ob_vector == NIL\n"); fflush(stdout); */
		return OK;
		}
	if (s_l() == NIL) {
		c_obj_k(EMPTY);
		return error("VEC::freeself(): "
		"s_l() == NIL");
		}
	z = s_s();
	if (z == NIL) {
		my_free(d.ob_vector);
		c_obj_k(EMPTY);
		return erg;
		}
	if (s_li() > 0) {
		for (i = 0; i < s_li(); i++, z++) {
			if (! z->emptyp())
				if (z->s_obj_k() != INTEGER)
					erg += z->freeself_debug(pd1);
			}
		my_free(s_s());
		}
	freeall(s_l()); 
	my_free(d.ob_vector);
	c_obj_k(EMPTY);
	return erg;
}

#if TEXDOCU
INT VECTOR_OB::m_il(INT il)
#else
make\_integerlength\_vector:
allocates space for il elements. 
All elements are set to EMPTY.
#endif
{
	INT erg = OK, i;
	
	if (il < 0) 
		return error("VECTOR::m_il(): negative length");
	if (! emptyp())
		freeself();
	if (il == 0) 
		erg += b_ls(
			(INTEGER_OP) callocobject("VEC::m_il: l"), NULL);
	else {
		erg += b_ls((INTEGER_OP)callocobject("VEC::m_il: l"),
			(SYM_OP) my_malloc(il * sizeof(struct object), "VEC::m_il: s"));
		}
	s_l()->c_obj_k(INTEGER);
	s_l()->m_i(il);
	for (i=0; i < il; i++)
		s_i(i)->c_obj_k(EMPTY);
	if (erg != OK)
		error("VECTOR::m_il(): error during computation");
	return erg;
}

#if TEXDOCU
INT VECTOR_OB::m_il_n(INT il)
#else
make\_integerlength\_null\_vector:
Allocates space for il vector elements.
The elements are initialized to 0 (INTEGER\_OB).
#endif
{
	INT i;
	
	m_il(il);
	for (i = 0; i < s_li(); i++)
		m_ii(i, 0);
	return OK;
}

#if TEXDOCU
struct vector *callocvectorstruct()
#else
Internal function for allocating struct vector.
#endif
{
	struct vector * ergebnis =
	(struct vector *) my_malloc(sizeof(struct vector), "callocvectorstruct()");
	if (ergebnis == NULL) 
		error("callocvectorstruct: no memory");
	return(ergebnis);
}

#if TEXDOCU
INT VECTOR_OB::b_ls(INTEGER_OP length, SYM_OP self)
#else
build\_length\_self: allocates a vectorstruct and puts self and 
length into it.
#endif
{
	OBJECTSELF d;

	if (! emptyp()) 
		freeself();
	d.ob_vector = callocvectorstruct();

	b_ks(VECTOR, d);
	
	c_s(self);
	c_l(length);
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::realloc_z(INT new_len)
#else
Reallocates the vector to the given length.
#endif
{
	INTEGER_OB nl;
	INT erg = OK;
	
	nl.m_i(new_len);
	erg += realloc(&nl);
	return(erg);
}

#if TEXDOCU
INT VECTOR_OB::realloc(INTEGER_OP new_len)
#else
Reallocates the vector to the given length.
The length is specified as an object.
#endif
{
	INT old_len, i;
	VECTOR_OB new_vec;
	INTEGER_OP length;
	SYM_OP mem;
	
	old_len = s_li();
	if (new_len->s_i() <= old_len) {
		for (i = new_len->s_i(); i < old_len; i++) {
			s_i(i)->freeself();
			}
		s_l()->m_i(new_len->s_i());
		return OK;
		}
	length = (INTEGER_OP) my_malloc(sizeof(struct object), "VEC::realloc(): l");
	length->c_obj_k(EMPTY);
	new_len->copy(length);
	mem = (SYM_OP) my_malloc(new_len->s_i() * sizeof(struct object), "VEC::realloc(): s");
	if (mem == NIL)
		return error("VECTOR::realloc(): no memory");
	if (new_vec.b_ls(length, mem) != OK)
		return error("VECTOR::realloc(): error in b_ls_v()");
	for (i = 0; i < new_len->s_i(); i++)
		new_vec.s_i(i)->c_obj_k(EMPTY);
	for (i = 0; i < old_len; i++)
		new_vec.s_i(i)->swap(s_i(i));
	swap(&new_vec);
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::append_i(INT i)
#else
#endif
{
	INT l;
	
	l = s_li();
	inc();
	m_ii(l, i);
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::append_element(INTEGER_OP len, SYM_OP b)
#else
Appends the object b to the vector. 
len is the length of the vector 
(i.e. the number of used elements).
#endif
{
	return append_element_itself(len, b, FALSE);
}

#if TEXDOCU
INT VECTOR_OB::append_element_itself(INTEGER_OP len, SYM_OP b, INT f_itself)
#else
Appends the element, swaps the object into the vector 
if f\_itself is TRUE.
#endif
{
	INT N, new_n;
	
	N = s_li();
	if (len->s_i() >= N) {
		new_n = len->s_i() + VECTOR_OVERSIZE;
		if (realloc_z(new_n) != OK)
			return error("VECTOR::append_element_itself(): "
			"error in v_realloc_z()");
		}
	if (f_itself)
		b->swap(s_i(len->s_i()));
	else
		b->copy(s_i(len->s_i()));
	len->inc();
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::insert_at(INTEGER_OP len, INT i, SYM_OP b)
#else
Inserts object b into the vector at position i.
Len is incremented.
The vector is reallocated if necessary.
#endif
{
	return insert_at_itself(len, i, b, FALSE);
}

#if TEXDOCU
INT VECTOR_OB::insert_at_itself(INTEGER_OP len, INT i, SYM_OP b, INT f_itself)
#endif
{
	INT N, new_n, l;
	
	if (i > len->s_i())
		return error("VECTOR::insert_at(): i > len");
	N = s_li();
	if (len->s_i() >= N) {
		new_n = len->s_i() + VECTOR_OVERSIZE;
		if (realloc_z(new_n) != OK)
			return error("VECTOR::insert_at_itself(): "
			"error in v_realloc_z()");
		}
	for (l = len->s_i(); l > i; l--) {
		s_i(l)->swap(s_i(l - 1));
		}
	if (f_itself)
		b->swap(s_i(i));
	else
		b->copy(s_i(i));
	len->inc();
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::delete_ith(INTEGER_OP len, INT i)
#else
Deletes the element i of the vector. All later elements 
are moved one position down.
#endif
{
	INT l;
	
	if (i < 0)
		return error("VECTOR::delete_ith(): i < 0");
	if (i >= len->s_i())
		return error("VECTOR::delete_ith(): i >= len");
	for (l = i + 1; l < len->s_i(); l++) {
		s_i(l)->swap(s_i(l - 1));
		}
	len->dec();
	s_i(len->s_i())->freeself();
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::v_del_ith2(INT len, INT i)
#endif
{
	INT j;
	
	if (i < 0)
		return error("v_del_ith2()|i < 0");
	s_i(i)->freeself();
	for (j = i; j < len - 1; j++)
		s_i(j)->swap(s_i(j + 1));
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::v_insert_at(INT i, SYM_OP b)
#endif
{
	INT len, l;
	
	len = s_li();
	if (i > len)
		return error("v_insert_at()|i > len");
	inc();
	for (l = len; l > i; l--)
		s_i(l)->swap(s_i(l - 1));
	b->copy(s_i(i));
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::v_realloc(INT len)
#endif
{
	INT old_len = s_li();
	INT erg = OK, i;
	VECTOR_OB p1;

	if (len < old_len)
		return error("v_realloc()|len < old_len");
	erg += p1.m_il(len);
	for (i = 0; i < len; i++)
		p1.s_i(i)->c_obj_k(EMPTY);
	for (i = 0; i < old_len; i++)
		s_i(i)->swap(p1.s_i(i));
	swap(&p1);
	return erg;
}

#if TEXDOCU
INT VECTOR_OB::v_shorten(INT len)
#endif
{
	INT old_len = s_li();
	INT erg = OK, i;
	VECTOR_OB p1;

	if (len < 0)
		return error("v_shorten()|len < 0");
	if (len > old_len)
		return error("v_shorten()|len > old_len");
	erg += p1.m_il(len);
	for (i = 0; i < len; i++)
		p1.s_i(i)->c_obj_k(EMPTY);
	for (i = 0; i < len; i++)
		p1.s_i(i)->swap(s_i(i));
	swap(&p1);
	return erg;
}

#if TEXDOCU
INT VECTOR_OB::v_minus_v_ip(VECTOR_OP q)
#else
this := this $\setminus$ q.
q must be sorted ascendingly.
#endif
{
	INT i, len, new_len, idx, f_found;
	INT erg = OK;

	len = s_li();
	new_len = 0;
	for (i = 0; i < len; i++) {
		erg += q->search(q->s_li(), 
			TRUE /* f_ascending */, 
			s_i(i), &idx, &f_found);
		if (f_found)
			continue;
		if (new_len != i)
			s_i(new_len)->swap(s_i(i));
		new_len++;
		}
	if (new_len < len)
		v_shorten(new_len);
	return erg;
}

#if TEXDOCU
INT VECTOR_OB::v_minus_v(VECTOR_OP q, VECTOR_OP p_minus_q)
#else
p\_minus\_q := this $\setminus$ q.
q must be sorted ascendingly.
#endif
{
	INT i, len, len_pq, idx, f_found;
	INT erg = OK;
	
	len = s_li();
	p_minus_q->m_il(len);
	len_pq = 0;
	for (i = 0; i < len; i++) {
		erg += q->search(q->s_li(), 
			TRUE /* f_ascending */, 
			s_i(i), &idx, &f_found);
		if (f_found)
			continue;
		erg += s_i(i)->copy(
			p_minus_q->s_i(len_pq));
		len_pq++;
		}
	if (len_pq < len)
		p_minus_q->v_shorten(len_pq);
	return erg;
}

#if TEXDOCU
INT VECTOR_OB::do_search(
	INT len, INT f_ascending, SYM_OP v, 
	INT *idx, INT *f_found, 
	INT type, void *data)
#endif
{
	return(VS_search(s_s(), len, 
		f_ascending, v, idx, f_found, type, data));
}

#if TEXDOCU
INT VECTOR_OB::search(INT len, INT f_ascending, SYM_OP v, 
	INT *idx, INT *f_found)
#endif
{
	return(VS_search(s_s(), len, 
		f_ascending, v, idx, f_found, 0, NIL));
}

#if TEXDOCU
INT VS_search(
	SYM_OP p, INT len, 
	INT f_ascending, SYM_OP v, 
	INT *idx, INT *f_found, 
	INT type, void *data)
#else
Binary search in the array p of length len for the object v.
Assumes the array p is sorted, ascendingly if f\_ascendingly is TRUE.
Returns in idx the position of the element which is strictly 
larger than v. This is the position where one usually has to insert 
v to add it to the sorted list. 
This function works also for vectors which have repeated entries.
In this case one adds a new element at the end of the list of 
equal elements.
Note that idx can become equal to len.
The function sets f\_found to TRUE if the element at position 
idx - 1 is equal to v. 
Usually p is the self-component of a vector.
#endif
{
	INT l, r, m, res;
	
	if (v == NIL || idx == NIL || f_found == NIL)
		return error("VS_search(): args nil");
	if (len == 0) {
		*idx = 0;
		*f_found = FALSE;
		return(OK);
		}
	l = 0;
	r = len;
	*f_found = FALSE;
	/* Invariante:
	 * p[i] <= v fuer i < l;
	 * p[i] >  v fuer i >= r;
	 * r - l ist Laenge des Suchbereichs. 
	 */
	while (l < r) {
		m = (l + r) >> 1;
		/* Bei einem Suchbereich gerader Laenge 
		 * wird das Element 
		 * oberhalb der Mitte genommen ! */
		res = do_comp(p + m, v, type, data);
		/* res = (p + m)->sym_comp(v); */
		if (!f_ascending)
			res = - res;
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

#if TEXDOCU
INT VECTOR_OB::insert_sorted(INTEGER_OP len, INT f_ascending, SYM_OP v)
#endif
{
	INT idx, f_found;

	if (search(len->s_i(), f_ascending, 
		v, &idx, &f_found) != OK)
		return error("VECTOR::insert_sorted(): "
		"error in v_search()");
	if (f_found) {
		printf("VECTOR::insert_sorted(): "
		"found at %ld: ", idx - 1);
		v->println();
		return error("VECTOR::insert_sorted(): f_found");
		}
	if (insert_at(len, idx, v) != OK)
		return error("VECTOR::insert_sorted(): "
		"error in v_insert_at()");
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::search_and_insert(INT len, SYM_OP v)
#endif
{
	INT idx, f_found, i;

	if (search(len, TRUE /* f_ascending */, 
		v, &idx, &f_found) != OK)
		return error("VECTOR::search_and_insert(): "
		"error in search()");
	if (f_found) {
		printf("VECTOR::search_and_insert(): "
		"found at %ld: ", idx - 1);
		v->println();
		return ERROR;
		}
	for (i = len; i > idx; i--) {
		s_i(i)->swap(s_i(i - 1));
		}
	v->copy(s_i(idx));
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::search_and_insert_int(INT len, INT i)
#endif
{
	INTEGER_OB int_ob;

	int_ob.m_i(i);
	return search_and_insert(len, &int_ob);
}

#if TEXDOCU
INT VECTOR_OB::set_minus(
	VECTOR_OP B, INTEGER_OP len_A, 
	INTEGER_OP len_B, INT f_ascending)
#endif
{
	SYM_OP arrA;
	struct object tmp_ob; /* we use our own swap */
	INT idx, l, f_found, i, nb_deleted;
	BYTE *to_be_deleted = NIL;
	
	to_be_deleted = (BYTE *) my_malloc(len_A->s_i() * sizeof(BYTE), "VEC::set_minus()");
	for (l = 0; l < len_A->s_i(); l++)
		to_be_deleted[l] = 0;
	nb_deleted = 0;
	for (l = len_B->s_i() - 1; l >= 0; l--) {
		if (search(len_A->s_i(), f_ascending, 
			B->s_i(l), &idx, &f_found) != OK)
			return error("VECTOR::set_minus(): error in v_search()");
		if (!f_found)
			continue;
		idx--;
		to_be_deleted[idx] = 1;
		nb_deleted++;
#if FALSE
		if (v_delete_ith(A, len_A, idx) != OK)
			return error("v_set_minus(): error in v_delete_ith()");
#endif
		}
	arrA = s_s();
	i = 0; /* next element inserted at i */
	for (l = 0; l < len_A->s_i(); l++) {
		if (to_be_deleted[l]) {
			(arrA + l)->freeself();
			continue; /* don't inc i */
			}
		if (l != i) {
			tmp_ob = *(OP)(arrA + i);
			*(OP)(arrA + i) = *(OP)(arrA + l);
			*(OP)(arrA + l) = tmp_ob;
			}
		i++;
		}
	len_A->m_i(len_A->s_i() - nb_deleted);
	my_free(to_be_deleted);
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::set_minus2(SYM_OP B, INTEGER_OP len_A, 
	INT len_B, INT f_ascending)
#endif
/* former name: vs_set_minus(VECTOR_OP A, OP B, ... ) */
{
	INT idx, l, f_found;
	
	for (l = len_B - 1; l >= 0; l--) {
		if (search(len_A->s_i(), f_ascending, 
			B + l, &idx, &f_found) != OK)
			return error("VECTOR::set_minus2(): "
			"error in v_search()");
		if (!f_found)
			continue;
		idx--;
		if (delete_ith(len_A, idx) != OK)
			return error("VECTOR::set_minus2(): "
			"error in v_delete_ith()");
		}
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::subseteq(VECTOR_OP B, INT lenA, INT lenB, 
	INT f_ascending, INT *f_subseteq)
#endif
{
	INT l, idx, f_found;
	
	if (lenB < lenA) {
		*f_subseteq = FALSE;
		return(OK);
		}
	for (l = 0L; l < lenA; l++) {
		if (B->search(lenB, f_ascending, 
			s_i(l), &idx, &f_found) != OK)
			return error("VECTOR::subseteq(): "
			"error in v_search()");
		if (!f_found) {
			*f_subseteq = FALSE;
			return(OK);
			}
		}
	*f_subseteq = TRUE;
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::equal(VECTOR_OP B, INT lenA, INT lenB, 
	INT f_ascending, INT *f_equal)
#else
Sets f\_equal to TRUE ia A and B contain the same elements.
(have the same length and all elements of A are contained in B).
#endif
{
	if (lenA != lenB) {
		*f_equal = FALSE;
		return(OK);
		}
	return(subseteq(B, lenA, lenB, f_ascending, f_equal));
}

#if TEXDOCU
INT VS_subseteq(SYM_OP A, SYM_OP B, 
	INT lenA, INT lenB, 
	INT f_ascending, INT *f_subseteq)
#else
Sets f\_subseteq to TRUE if A is a subset of B.
#endif
{
	INT l, idx, f_found;
	
	if (lenB < lenA) {
		*f_subseteq = FALSE;
		return OK;
		}
	for (l = 0; l < lenA; l++) {
		if (VS_search(B, lenB, f_ascending, 
			A + l, &idx, &f_found, 
			0 /* type */, NIL /* data */) != OK)
			return error("VS_subseteq(): error in VS_search()");
		if (!f_found) {
			*f_subseteq = FALSE;
			return(OK);
			}
		}
	*f_subseteq = TRUE;
	return(OK);
}

static INT q2sort(VECTOR_OP p, 
	INT f_ascending, INT left, INT right);
static INT partition(VECTOR_OP p, 
	INT f_ascending, INT left, INT right, 
	INT *middle); 

#if TEXDOCU
INT VECTOR_OB::quicksort(INT len, INT f_ascending)
#else
Sorts the first len elements in the vector ascendingly / descendingly 
using quicksort.
#endif
{
	if (q2sort(this, f_ascending, 0, len - 1) != OK)
		return error("VECTOR::quicksort: error in q2sort()");
	return(OK);
}

static INT q2sort(VECTOR_OP p, INT f_ascending, INT left, INT right)
{
	INT middle;
	
	if (left < right) {
		if (partition(p, f_ascending, left, right, &middle) != OK)
			return error("q2sort: error in partition()");
		if (q2sort(p, f_ascending, left, middle - 1) != OK)
			return error("q2sort: error in q2sort(left part)");
		if (q2sort(p, f_ascending, middle + 1, right) != OK)
			return error("q2sort: error in q2sort(right part)");
		}
	return(OK);
}

static INT partition(VECTOR_OP p, 
	INT f_ascending, INT left, INT right, INT *middle) 
{
	INT l, r, m, len, m1, res;
	SYM_OP pivot;
	
	/* Pivot Strategie: nimm' mittleres Element: */
	len = right + 1 - left;
	m1 = len >> 1;
	pivot = p->s_i(left);
	if (m1)
		pivot->swap(p->s_i(left + m1));
	l = left;
	r = right;
	while (l < r) {
		while (TRUE) {
			if (l > right)
				break;
			res = p->s_i(l)->sym_comp(pivot);
			if (!f_ascending)
				res *= -1;
			if (res > 0)
				break;
			l++;
			}
		while (TRUE) {
			if (r < left)
				break;
			res = p->s_i(r)->sym_comp(pivot);
			if (!f_ascending)
				res *= -1;
			if (res <= 0)
				break;
			r--;
			}
		if (l < r)
			p->s_i(l)->swap(p->s_i(r));
		}
	m = r;
	if (left != m)
		p->s_i(left)->swap(p->s_i(m));
	*middle = m;
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::nullp()
#else
TRUE if all elements satisfy nullp().
#endif
{
	INT i;
	
	for (i = 0; i < s_li(); i++) {
		if (! s_i(i)->nullp())
			return FALSE;
		}
	return TRUE;
}

#if TEXDOCU
INT VECTOR_OB::einsp()
#else
TRUE if all elements satisfy onep().
#endif
{
	INT i;
	
	for (i = 0; i < s_li(); i++) {
		if (! s_i(i)->einsp())
			return FALSE;
		}
	return TRUE;
}

#if TEXDOCU
INT VECTOR_OB::vectorp()
#endif
{
	if ((s_obj_k() == VECTOR) ||
		(s_obj_k() == WORD) ||
		(s_obj_k() == KRANZ) ||
		(s_obj_k() == COMP))
		return TRUE;
	return FALSE;
}

#if TEXDOCU
INT VECTOR_OB::m_o(SYM_OP ob)
#else
make\_object\_vector:
Builds a vector object containing ob as its single element.
#endif
{ 
	m_il(1);
	return (ob->copy(s_i(0))); 
}

#if TEXDOCU
INT VECTOR_OB::b_o(SYM_OP ob)
#else
build\_object\_vector: the object pointer is used for the self part of the vector.
#endif
{ 
	if (! emptyp())
		freeself();
	b_ls((INTEGER_OP) callocobject("VEC::b_o"), ob);
	s_l()->m_i(1);
	return(OK); 
}

#if TEXDOCU
INT VECTOR_OB::add_apply(VECTOR_OP b)
#else
Componentwise addition of elements: b := this + b.
#endif
{
	INT i;
	VECTOR_OP c;

	if (s_li() > b->s_li()) {
		c = (VECTOR_OP) callocobject("VEC::add_apply");
		copy(c);
		for (i = 0; i < s_li(); i++) {
			if (i < b->s_li()) {
				s_i(i)->add(b->s_i(i), c->s_i(i));
				}
			else break;
			}
		b->freeself();
		*b = *c;
		c->c_obj_k(EMPTY);
		freeall(c);
		}
	else {
		for (i = 0; i < b->s_li(); i++) {
			if (i < s_li())
				s_i(i)->add_apply(b->s_i(i));
			else break;
			}
		}
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::add(VECTOR_OP b, VECTOR_OP c)
#endif
{
	INT i;

	if (s_li() > s_li()) {
		copy(c);
		for (i = 0; i < s_li(); i++) {
			if (i < b->s_li())
				s_i(i)->add(b->s_i(i), c->s_i(i));
			else break;
			}
		}
	else {
		b->copy(c);
		for (i = 0; i < b->s_li(); i++) {
			if (i < s_li())
				s_i(i)->add(b->s_i(i), c->s_i(i));
			else break;
			}
		}
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::addinvers(VECTOR_OP erg)
#else
erg becomes vector holding additive inverses of elements of this.
#endif
{
	INT i;

	copy(erg);
	for (i = 0; i < s_li(); i++)
		s_i(i)->addinvers(erg->s_i(i));
	return(OK);
}


#if TEXDOCU
INT VECTOR_OB::addinvers_apply()
#endif
{
	INT i, erg = OK;
	
	for (i = 0; i < s_li(); i++) 
		erg += s_i(i)->addinvers_apply();
	return erg;
}


#if TEXDOCU
INT VECTOR_OB::addtoallelements(SYM_OP zahl, VECTOR_OP ergebnis)
#endif
/* before: addtoallvectorelements(
 *    OP zahl, OP vector, OP ergebnis) */
{
	INT i;
	
	copy(ergebnis);
	for (i = 0; i < s_li(); i++)
	    zahl->add(ergebnis->s_i(i), ergebnis->s_i(i));
	return OK;
}


#if TEXDOCU
INT VECTOR_OB::copy(VECTOR_OP res)
#endif
{
	INT i;

	if (s_obj_k() == EMPTY) {
		res->freeself();
		return OK;
		}
	res->m_il(s_li());

	for (i = 0; i < s_li(); i++) {
		if (s_i(i)->s_obj_k() == INTEGER)
			* res->s_i(i) = * s_i(i);
		else
			s_i(i)->copy(res->s_i(i));
		}
	res->c_obj_k(s_obj_k());
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::compare(VECTOR_OP b)
#endif
{
	INT i, erg;
	
	for (i = 0; i < s_li(); i++) {
		if (i >= b->s_li())
			return(1);
		erg = s_i(i)->sym_comp(b->s_i(i));
		if (erg != 0)
			return(erg);
		}
	return(0);
}

#if TEXDOCU
INT VECTOR_OB::lastof(SYM_OP b)
#endif
{
	return (s_i(s_li() - 1)->copy(b));
}

#if TEXDOCU
INT VECTOR_OB::length(INTEGER_OP b)
#endif
{
	return (s_l()->copy(b));
}

#if TEXDOCU
INT VECTOR_OB::inc()
#else
Reallocates the vector by one element in the end which will 
be empty. The remaining entries are swaped into the new vector.
#endif
{
	INT i;
	SYM_OP z;


	if (s_li() == 0) {
		z = (SYM_OP) my_malloc(sizeof(struct object), "VEC::inc()");
		z->c_obj_k(EMPTY);
		}
	else {
		i = (s_li() + 1) * (sizeof(struct object));
		z = (SYM_OP) my_malloc(i, "VEC::inc()");
		my_memcpy(z, s_s(), s_li() * sizeof(struct object));
		(z + s_li())->c_obj_k(EMPTY);
		my_free(s_s());
		}

	if (z == NULL)
		error("VECTOR::inc:self == NULL");
	c_s(z);
	s_l()->inc();
	s_i(s_li() - 1)->c_obj_k(EMPTY);
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::dec()
#else
Decrenent the vector by one element. 
The element is freed and the length of the vector is decreased by one 
(the vector is not reallocated).
#endif
{
	if (s_li() == 0) {
		error("vector der laenge 0 in decvector");
		return(ERROR);
		}
	if (! s_i(s_li() - 1)->emptyp())
		s_i(s_li() - 1)->freeself();
	/* freigeben des speicherplatzes 
	 * des letzten vectorelements */
	s_l()->dec();
	/* verkuerzen der laenge um eins */
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::append(SYM_OP b, VECTOR_OP c)
#else
Append the vector b to this and put the result into c.
If b is not a vector, the element b will be added 
as last element of the vector.
#endif
{
	INT i, length;
	INT erg = OK;
	VECTOR_OP b1;

	if (b->s_obj_k() != VECTOR) {
		VECTOR_OP d = (VECTOR_OP) callocobject("VEC::append");
		erg += d->m_o(b);
		erg += append((SYM_OP)d, c);
		erg += freeall(d);
		return erg;
		}
	b1 = (VECTOR_OP)b;
	length = s_li() + b1->s_li();
	erg += c->m_il(length);
	for (i = 0; i < length; i++) {
		if (i < s_li())
			erg += s_i(i)->copy(c->s_i(i));
		else
			erg += b1->s_i(i - s_li())->copy(c->s_i(i));
		}
	return erg;
}

#if TEXDOCU
INT VECTOR_OB::append_in_place(SYM_OP b)
#else
Append the object b to the vector.
#endif
{
	INT l;
	
	l = s_li();
	inc();
	b->copy(s_i(l));
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::max_vector(SYM_OP m)
#else
Search for the maximal element in the vector and copy this element to m.
#endif
{
	INT i;
	SYM_OP zm;
	
	zm = s_i(0);
	for (i = 1; i < s_li(); i++)
		if (s_i(i)->gt(zm))
			zm = s_i(i);
	return (zm->copy(m));
}

#if TEXDOCU
INT VECTOR_OB::sum(INTEGER_OP ergebnis)
#else
berechnet die summe der vectorelemente;
nur fuer INTEGER vectoren.
#endif
{
	INT i;
	
	ergebnis->m_i(0);
	for (i = 0; i < s_li(); i++)
		ergebnis->add(s_i(i), ergebnis);
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::mult_scalar(SYM_OP a, VECTOR_OP c)
#else
scalar multiplication: all elements get multiplied by a.
The result is put into c.
#endif
{
	INT i = 0;
	INT erg = OK;
	
	erg += c->m_il(s_li());
	for (i = 0; i < c->s_li(); i++)
		erg += a->mult(s_i(i), c->s_i(i));
	return erg;
}

#ifdef MATRIXTRUE
#if TEXDOCU
INT VECTOR_OB::mult_matrix(MATRIX_OP b, VECTOR_OP c)
#else
/* before: mult\_vector\_matrix(OP a, OP b, OP c) */
#endif
{
	INT i, j;
	INT erg = OK;
	SYM_OP d;
	
	if (s_obj_k() != VECTOR) return ERROR;
	if (b->s_obj_k() != MATRIX) return ERROR;
	if (this == c) return ERROR;
	if (s_li() != b->s_hi()) return ERROR;
	erg += c->m_il_n(b->s_li());
	d = callocobject("VEC::mult_matrix");
	for (i = 0; i < c->s_li(); i++) {
		for (j = 0; j < s_li(); j++) {
			erg += s_i(j)->mult(b->s_ij(j, i), d);
			erg += d->add_apply(c->s_i(i));
			}
		}
	erg += freeall(d);
	return erg;
}
#endif /* MATRIXTRUE */

#if TEXDOCU
INT VECTOR_OB::mult_vector(VECTOR_OP b, VECTOR_OP c)
#else
componentwise multiplication
#endif
{
	INT i = 0;

	if (s_li() !=  b->s_li()) {
		fprintf(stderr,"size of a: %ld\n", s_li());
		fprintf(stderr,"size of b: %ld\n", b->s_li());
		return error("VECTOR::mult_vector: "
		"different size of vectors ");
		}
	else {
		if (! c->emptyp())
			c->freeself();
		c->m_il(s_li());
		for (i = 0; i < b->s_li(); i++)
			s_i(i)->mult(b->s_i(i), c->s_i(i));
		}
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::scalarproduct(VECTOR_OP b, SYM_OP d)
#else
scalar product of two vectors this and b. The result 
is put into d.
#endif
{
	INT i;
	INT erg = OK;
	SYM_OB c;

	if (s_li() != b->s_li()) {
		error("VECTOR::scalarproduct:different length");
		return(ERROR); 
		}
	for (i = s_li() - 1; i >= 0; i--) {
		erg += s_i(i)->mult(b->s_i(i), &c);
		erg += c.add_apply(d);
		}
	return erg;
}

#if TEXDOCU
INT VECTOR_OB::mult_apply_vector(VECTOR_OP b)
#else
componentwise multiplication: b := this * b.
#endif
{
	INT i = 0;

	if (s_li() !=  b->s_li()) {
		fprintf(stderr,"size of a: %ld\n", s_li());
		fprintf(stderr,"size of b: %ld\n", b->s_li());
		return error("VECTOR::mult_apply_vector: "
		"different size of vectors ");
		}
	else {
		for (i = 0; i < b->s_li(); i++)
			s_i(i)->mult_apply(b->s_i(i));
		};
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::mult_apply(SYM_OP b)
#else
/* before: mult\_apply\_vector(OP a, OP b) */
#endif
{
	INT erg = OK;
	
	switch (b->s_obj_k()) {
		case VECTOR: 
			erg += mult_apply_vector((VECTOR_OP)b); 
			break;
		default:
			erg = error("VECTOR::mult_apply: wrong type");
			break;
		}
	return erg;
}

#if TEXDOCU
INT VECTOR_OB::search_linear(SYM_OP p)
#else
Linear search in the vector for the element p.
Returns -1 if p is not found, otherwise the position where it was found.
#endif
{
	INT i, l, r;

	l = s_li();
	for (i = 0; i < l; i++) {
		r = p->sym_comp(s_i(i));
		if (r == 0)
			return i;
		}
	return -1;
}

#if TEXDOCU
INT VECTOR_OB::join2(VECTOR_OP V1, VECTOR_OP V2)
#else
The vectors V1 and V2 must be of equal length. 
A new vector is created (into the this object) 
of the same length as V1 and V2 holding 
vectors of length 2 with (V1[i], V2[i]) in its i-th position. 
#endif
{
	VECTOR_OB key;
	INT i, l;

	l = V1->s_li();
	if (V2->s_li() != l)
		return error("VEC::join2(): different legth");
	m_il(l);
	for (i = 0; i < l; i++) {
		key.m_il(2);
		V1->s_i(i)->copy(key.s_i(0));
		V2->s_i(i)->copy(key.s_i(1));
		key.swap(s_i(i));
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::multiplicities(VECTOR_OP val, VECTOR_OP mult)
#endif
{
	SYM_OP a;
	INT i, ii, idx, len, l, f_found;
	
	l = 0;
	val->m_il(l);
	mult->m_il(l);
	len = s_li();
	for (i = 0; i < len; i++) {
		a = s_i(i);
		if (val->search(l, TRUE /* f_ascending */, 
			a, &idx, &f_found) != OK)
			return error("VECTOR::multiplicities(): error in search()");
		if (f_found) {
				mult->s_i(idx - 1)->inc();
				}
		else {
			val->inc();
			mult->inc();
			for (ii = l; ii > idx; ii--) {
				val->s_i(ii)->swap(val->s_i(ii - 1));
				mult->s_i(ii)->swap(mult->s_i(ii - 1));
				}
			a->copy(val->s_i(idx));
			mult->m_ii(idx, 1);
			l++;
			}
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::classify(VECTOR_OP val, VECTOR_OP classes)
#endif
{
	SYM_OP a;
	VECTOR_OP c;
	INT i, ii, idx, len, l, f_found;
	
	l = 0;
	val->m_il(l);
	classes->m_il(l);
	len = s_li();
	for (i = 0; i < len; i++) {
		a = s_i(i);
		if (val->search(l, TRUE /* f_ascending */, 
			a, &idx, &f_found) != OK)
			return error("VECTOR::classify(): error in search()");
		if (f_found) {
				c = (VECTOR_OP) classes->s_i(idx - 1);
				c->inc();
				c->m_ii(c->s_li() - 1, i);
				}
		else {
			val->inc();
			classes->inc();
			for (ii = l; ii > idx; ii--) {
				val->s_i(ii)->swap(val->s_i(ii - 1));
				classes->s_i(ii)->swap(classes->s_i(ii - 1));
				}
			a->copy(val->s_i(idx));
			c = (VECTOR_OP) classes->s_i(idx);
			c->m_il(1);
			c->m_ii(0, i);
			l++;
			}
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::classify_and_reorder(VECTOR_OP val, VECTOR_OP classes, 
	PERMUTATION_OP p, VECTOR_OP class_lengths)
#endif
{
	VECTOR_OP pclass;
	INT i, ii, j, l, ll, n, a;
	
	classify(val, classes);
	l = classes->s_li();
	class_lengths->m_il(l);
	n = s_li();
	p->m_il(n);
	i = 0;
	for (ii = 0; ii < l; ii++) {
		pclass = (VECTOR_OP) classes->s_i(ii);
		ll = pclass->s_li();
		for (j = 0; j < ll; j++) {
			a = pclass->s_ii(j);
			p->m_ii(a, ++i);
			}
		class_lengths->m_ii(ii, ll);
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::norm(VECTOR_OP mult)
#endif
{
	VECTOR_OP p;
	INT i, l, l1;
	
	l = s_li();
	mult->m_il(l);
	for (i = 0; i < l; i++) {
		p = (VECTOR_OP) s_i(i);
		if (p->s_obj_k() != VECTOR)
			return error("VEC::norm entry is not a vector");
		l1 = p->s_li();
		mult->m_ii(i, l1);
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::sprint_multiplicities(VECTOR_OP mult, BYTE *str)
#endif
{
	SYM_OP p;
	VECTOR_OP key;
	INT i, j, l, ll, a, m;
	
	l = s_li();
	for (i = 0; i < l; i++) {
		m = mult->s_ii(i);
		p = s_i(i);
		if (p->s_obj_k() == INTEGER) {
			a = s_ii(i);
			sprintf(str + strlen(str), "%ld x %ld", m, a);
			}
		else if (p->s_obj_k() == VECTOR) {
			key = (VECTOR_OP) p;
			ll = key->s_li();
			sprintf(str + strlen(str), "%ld x ", m);
			for (j = 0; j < ll; j++) {
				a = key->s_ii(j);
				sprintf(str + strlen(str), "%ld", a);
				if (j < ll - 1)
					sprintf(str + strlen(str), "/");
				}
			}
		if (i < l - 1)
			strcat(str, ", ");
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::gcd_all_elements(SYM_OP g)
#else
Computes the gcd of all vector elements into g.
#endif
{
	SYM_OB g1, g2, u, v;
	INT i, l;
	
	l = s_li();
	if (l == 0)
		return error("gcd_vector() length zero!");
	if (l == 1) {
		s_i(0)->copy(g);
		return OK;
		}
	s_i(0)->copy(&g1);
	for (i = 1; i < l; i++) {
		bezout(&g1, s_i(i), &u, &v, &g2);
		g2.swap(&g1);
		}
	g1.copy(g);
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::gcd_all_elements_and_divide_out(SYM_OP g, VECTOR_OP V)
#else
Computes the gcd of all vector elements into g. 
V contains the vector after taking out that common factor.
#endif
{
	INT i, l;
	
	l = s_li();
	gcd_all_elements(g);
	V->m_il(l);
	for (i = 0; i < l; i++) {
		s_i(i)->ganzdiv(g, V->s_i(i));
		}
	return OK;
}

#if TEXDOCU
INT v_println(VECTOR_OP V, INT f_numerated)
#endif
{
	INT i, len;

	len = V->s_li();
	for (i = 0; i < len; i++) {
		if (f_numerated)
			printf("%ld: ", i);
		V->s_i(i)->println();
		}
	return OK;
}

#if TEXDOCU
INT v_do_println(VECTOR_OP V, INT f_numerated, INT n_on_a_row, 
	INT type, void *data)
#endif
{
	INT i, len;

	len = V->s_li();
	for (i = 0; i < len; i++) {
		if (f_numerated)
			printf("%ld: ", i);
		do_print(V->s_i(i), type, data);
		if (n_on_a_row > 0) {
			if (i < len - 1)
				printf(", ");
			if (((i + 1) % n_on_a_row) == 0)
				printf("\n");
			}
		else
			printf("\n");
		}
	if (n_on_a_row > 0)
		printf("\n");
	fflush(stdout);
	return OK;
}

#if TEXDOCU
INT v_do_fprintln(VECTOR_OP V, FILE *fp, INT f_numerated, INT n_on_a_row, 
	INT type, void *data)
#endif
{
	INT i, len;
	BYTE str[1024];

	len = V->s_li();
	for (i = 0; i < len; i++) {
		if (f_numerated)
			fprintf(fp, "%ld: ", i);
		str[0] = 0;
		do_sprint(V->s_i(i), str, type, data);
		fprintf(fp, "%s", str);
		if (n_on_a_row > 0) {
			if (i < len - 1)
				fprintf(fp, ", ");
			if (((i + 1) % n_on_a_row) == 0)
				fprintf(fp, "\n");
			}
		else
			fprintf(fp, "\n");
		}
	if (n_on_a_row > 0)
		fprintf(fp, "\n");
	fflush(stdout);
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::add_apply_elementwise(VECTOR_OP V2, INT sign)
#endif
{
	INT i, l;
	SYM_OP pa, pb;
	SYM_OB tmp;

	l = s_li();
	if (V2->s_li() != l)	
		return error("VEC::add_apply_elementwise() different length");
	for (i = 0; i < l; i++) {
		pa = s_i(i);
		pb = V2->s_i(i);
		if (sign < 0)
			pb->addinvers_apply();
		pa->add(pb, &tmp);
		tmp.swap(pa);
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::divide_out(SYM_OP val)
#endif
{
	INT i, l;
	SYM_OP pa;
	SYM_OB tmp, tmp2, tmp3;

	l = s_li();
	for (i = 0; i < l; i++) {
		pa = s_i(i);
		pa->ganzdiv(val, &tmp);
		tmp.mult(val, &tmp2);
		tmp2.addinvers_apply();
		pa->add(&tmp2, &tmp3);
		if (!tmp3.nullp()) {
			printf("division problem !, the vector is:\n");
			Print();
			printf("the divisor: ");
			val->println();
			fflush(stdout);
			return error("VEC::divide_out(): no integer division !");
			}
		tmp.swap(pa);
		// (*pa) /= (*val);
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::multiply_elementwise(SYM_OP val)
#endif
{
	INT i, l;
	SYM_OP pa;

	l = s_li();
	for (i = 0; i < l; i++) {
		pa = s_i(i);
		(*pa) *= (*val);
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::apply_perm(PERMUTATION_OP p)
#else
applies p to the entries of the vector. 
This means that $v[i]$ is mapped to $v[p(i)]$.
#endif
{
	VECTOR_OB I;
	INT l, i, i1;
	
	l = s_li();
	I.m_il(l);
	for (i = 0; i < l; i++) {
		i1 = p->s_ii(i) - 1;
		s_i(i)->copy(I.s_i(i1));
		}
	swap(&I);
	return TRUE;
}

#if TEXDOCU
INT VECTOR_OB::apply_perm_to_vector_of_pairs(PERMUTATION_OP p, INT n)
#else
applies p to the entries of the vector, whose entries are indexed by the pairs of 
elements 1,...,n. 
This means that $v[ij2k(i, j)]$ is mapped to 
$v[ij2k(p(i), p(j))]$ for all $0 \le i < j < n$.
#endif
{
	VECTOR_OB I;
	INT l, i, i1, j, j1, k, k1;
	
	l = s_li();
	I.m_il(l);
#if 0
	printf("VECTOR_OB::apply_perm_to_vector_of_pairs()\n");
	printf("p=");
	p->println();
	printf("vector=");
	println();
#endif
	for (i = 0; i < n; i++) {
		i1 = p->s_ii(i) - 1;
		for (j = i + 1; j < n; j++) {
			j1 = p->s_ii(j) - 1;
			k = ij2k(i, j, n);
			k1 = ij2k(i1, j1, n);
			// printf("i=%ld i1=%ld j=%ld j1=%ld k=%ld k1=%ld v[k]=%ld\n", i, i1, j, j1, k, k1, s_ii(k));
			s_i(k)->copy(I.s_i(k1));
			}
		}
#if 0
	printf("I=\n");
	I.println();
#endif
	swap(&I);
	return TRUE;
}

#if TEXDOCU
INT VECTOR_OB::m_standard_1_n(INT l)
#else
Makes the vector $[1,2,\ldots,l+1]$.
#endif
{
	INT i;
	
	m_il(l);
	for (i = 0; i < l; i++) {
		s_i(i)->m_i_i(i + 1);
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::m_standard_0_nm1(INT l)
#else
Makes the vector $[0,1,\ldots,l]$.
#endif
{
	INT i;
	
	m_il(l);
	for (i = 0; i < l; i++) {
		s_i(i)->m_i_i(i);
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::m_intvect1(INT a)
#else
Makes the vector of length one holding the integer a.
#endif
{
	m_il(1);
	s_i(0)->m_i_i(a);
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::m_i_times_j(INT i, INT j)
#else
Makes the vector $[j,\ldots,j]$ of length $i$.
#endif
{
	INT a;
	
	m_il(i);
	for (a = 0; a < i; a++)
		s_i(a)->m_i_i(j);
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::power_elementwise(INT k, VECTOR_OP q)
#else
Makes q the vector holding the $k$-th powers of the elements of this.
#endif
{
	INT i, l;
	
	l = s_li();
	q->m_il(l);
	for (i = 0; i < l; i++) {
		s_i(i)->power_int(q->s_i(i), k);
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::power_elementwise_apply(INT k)
#else
Raises all elements to the $k$-th power using power\_int\_apply.
#endif
{
	INT i, l;
	
	l = s_li();
	for (i = 0; i < l; i++) {
		s_i(i)->power_int_apply(k);
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::first_subset(INT n)
#else
Makes the first subset of $\{0,1,\ldots,n-1\}$. 
This is always the empty set.
#endif
{
	m_il(0);
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::next_subset(INT n)
#else
Makes the lexicographically next subset of $\{0,1,\ldots,n-1\}$. 
Returns TRUE if there is another subset.
#endif
{
	INT *set;
	INT i, l, a, ret;
	
	l = s_li();
	set = (INT *) my_malloc(n * sizeof(INT), "VEC::next_subset");
	for (i = 0; i < n; i++)
		set[i] = 0;
	for (i = 0; i < l; i++) {
		a = s_ii(i);
		set[a] = 1;
		}
	ret = FALSE;
	for (i = n - 1; i >= 0; i--) {
		if (set[i] == 0) {
			set[i] = 1;
			ret = TRUE;
			break;
			}
		else
			set[i] = 0;
		}
	l = 0;
	for (i = 0; i < n; i++)
		if (set[i])
			l++;
	m_il_n(l);
	l = 0;
	for (i = 0; i < n; i++)
		if (set[i]) {
			m_ii(l, i);
			l++;
			}
	my_free(set);
	return ret;
}

#if TEXDOCU
INT VECTOR_OB::first_subset_ordered(INT n)
#else
#endif
{
	m_il(0);
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::next_subset_ordered(INT n)
#else
#endif
{
	INT *choice;
	INT i, k, a, ret;
	
	k = s_li();
	choice = (INT *) my_malloc(n * sizeof(INT), "VEC::next_subset_ordered");
	for (i = 0; i < k; i++) {
		a = s_ii(i);
		choice[i] = a;
		}
	ret = FALSE;
	while (k < n) {
		if (n_Choose_k_next(choice, n, k)) {
			ret = TRUE;
			m_il(k);
			for (i = 0; i < k; i++) {
				a = choice[i];
				m_ii(i, a);
				}
			break;
			}
		k++;
		if (n_Choose_k_first(choice, n, k)) {
			ret = TRUE;
			m_il(k);
			for (i = 0; i < k; i++) {
				a = choice[i];
				m_ii(i, a);
				}
			break;
			}
		}
	my_free(choice);
	return ret;
}

#if TEXDOCU
INT VECTOR_OB::all_k_subsets(INT n, INT k)
#endif
{
	VECTOR_OB v;
	INT l = 0;
	
	m_il(0);
	v.first_k_subset(n, k);
	do {
		inc();
		v.copy((VECTOR_OP) s_i(l));
		l++;
		} while (v.next_k_subset(n, k));
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::first_k_subset(INT n, INT k)
#endif
{
	if (k <= n) {
		m_standard_0_nm1(k);
		return TRUE;
		}
	return FALSE;
}

#if TEXDOCU
INT VECTOR_OB::next_k_subset(INT n, INT k)
#endif
{
	INT *choice;
	INT i, a;
	
	if (k != s_li()) {
		return error("VECTOR_OB::next_k_subset() k != s_li()");
		}
	choice = (INT *) my_malloc(n * sizeof(INT), "VEC::next_k_subset");
	for (i = 0; i < k; i++) {
		a = s_ii(i);
		choice[i] = a;
		}
	if (n_Choose_k_next(choice, n, k)) {
		for (i = 0; i < k; i++) {
			a = choice[i];
			m_ii(i, a);
			}
		my_free(choice);
		return TRUE;
		}
	my_free(choice);
	return FALSE;
}

#if TEXDOCU
INT VECTOR_OB::is_subset(VECTOR_OP set2)
#endif
{
	INT f_subseteq;
	
	subseteq(set2, s_li(), set2->s_li(), 
		TRUE /* f_ascending */, &f_subseteq);
	return f_subseteq;
}

#if TEXDOCU
INT VECTOR_OB::is_fix_under(PERMUTATION_OP p)
#endif
{
	INT i, j, a, b, l;
	
	l = s_li();
	for (i = 0; i < l; i++) {
		a = s_ii(i);
		b = p->s_ii(a) - 1;
		for (j = 0; j < l; j++) {
			if (s_ii(j) == b)
				break;
			}
		if (j == l)
			return FALSE;
		}
	return TRUE;
}

#if TEXDOCU
INT VECTOR_OB::dec_all_entries()
#endif
{
	INT i, l;
	
	l = s_li();
	for (i = 0; i < l; i++) {
		s_i(i)->dec();
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::inc_all_entries()
#endif
{
	INT i, l;
	
	l = s_li();
	for (i = 0; i < l; i++) {
		s_i(i)->inc();
		}
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::canonicize_map(LABRA_OP G, LABRA_OP aut, 
	PERMUTATION_OP transporter, INT f_v)
#endif
{
	VECTOR_OB values, classes;
	VECTOR_OP c;
	INT n, k, *theX;
	INT a, i, j, l;
	INT f_vv = FALSE;
	LABRA_OB G1, G2;
	SYM_OB ago, go;
	MATRIX_OB TG;
	PERMUTATION_OB q, q1, q2;
	VECTOR_OB gen;
	
	G->copy(&G1);
	n = s_li();
	q1.m_il(n);
	q1.one();
	if (f_v) {
		printf("VEC::canonicize_map(): ");
		println();
		}
	classify(&values, &classes);
	if (f_v) {
		printf("values=");
		values.println();
		printf("classes=");
		classes.println();
		}
	l = values.s_li();
	for (i = l - 1; i >= 0; i--) {
		c = (VECTOR_OP) classes.s_i(i);
		k = c->s_li();
		theX = (INT *) my_malloc(k * sizeof(INT), "VEC::canonicize_map theX");
		for (j = 0; j < k; j++) {
			a = c->s_ii(j);
			theX[j] = a;
			}
		G1.calc_transversal_matrix(&TG, FALSE);
		if (f_v) {
			printf("TG=");
			TG.Print();
			printf("i=%ld: {", i);
			for (j = 0; j < k; j++) {
				printf("%ld ", theX[j]);
				}
			printf("} -> ");
			}
		geo_canonicize_set(n, k, theX, 
			&G1, &TG, f_v, f_vv, 
			TRUE /* f_get_aut_group */, &G2, &ago, &q);

#if 0
		geo_canon_simple(FALSE /* f_maxtest */, &back_to,
			1 /* nrow */, n /* ncol */, k /* nb_X */, 
			FALSE /* f_print_backtrack_points */, theX, 
			&p, &q, 
			FALSE /*f_transposed */, TRUE /* f_get_aut_group */, 
			&G2, f_v, f_vv);
#endif
		
		G2.group_order(&ago);
		G2.generators(&gen);
		if (f_v) {
			printf("{");
			for (j = 0; j < k; j++) {
				printf("%ld ", theX[j]);
				}
			printf("} ago=");
			ago.println();
			// printf("p="); p.println();
			printf("q="); q.println();
			}
		my_free(theX);
		apply_perm(&q);
		q1.mult(&q, &q2);
		q2.swap(&q1);
		if (f_v) {
			printf("permuted vector=");
			println();
			printf("stabilizer:");
			gen.Print();
			}
		G1.init_quick(&gen);
		// G1.freeself();
		// reduce_generators_labra(&gen, &go, FALSE /* f_verbose */, &G1);
		
		// G2.swap(&G1);
		classify(&values, &classes);
		if (f_v) {
			printf("values=");
			values.println();
			printf("classes=");
			classes.println();
			}
		}
	q1.copy(transporter);
	aut->init_quick(&gen);
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::Nabla(LABRA_OP G, POLYNOM_OP res, INT n, INT f_v)
#endif
{
	VECTOR_OB VC, vc;
	LABRA_OB aut;
	PERMUTATION_OB transporter;
	SYM_OB stab_i, stab_j, a, b, c;
	INT go;
	POLYNOM_OB p;
	POLYNOM_OP q;
	MONOM_OP m;
	VECTOR_OP v;
	POLYNOM_OP tmp, res1 = NIL;
	INT f_vv = FALSE;
	
	G->group_order(&a);
	go = a.s_i_i();
	
	copy(&VC);
	VC.canonicize_map(G, &aut, &transporter, f_vv);
	aut.group_order(&stab_i);
	if (f_v) {
		printf("Nabla(): ");
		VC.print();
		printf(" stab order i= ");
		stab_i.println();
		}
	
	nabla(&p, n, f_vv);
	q = &p;
	res->freeself();
	while (q != NIL) {
		m = q->s_mo();
		if (f_v) {
			m->print();
			printf(" -> ");
			}
		m->s_k()->copy(&a);
		v = m->s_s();
		v->copy(&vc);
		vc.canonicize_map(G, &aut, &transporter, f_vv);
		tmp = (POLYNOM_OP) callocobject("VEC::Nabla() tmp");
		tmp->m_skn(&vc, &a, NIL);
		if (f_v) {
			tmp->println();
			}
		if (res1 == NIL) {
			res1 = tmp;
			}
		else
			((LIST_OP) res1)->insert_list((LIST_OP) tmp, add_koeff, comp_monomvector_monomvector);
		q = q->s_n();
		}
	if (f_v) {
		printf("res1=\n");
		res1->Print();
		}
	q = res1;
	while (q != NIL) {
		m = q->s_mo();
		if (f_vv) {
			m->print();
			}
		m->s_k()->copy(&a);
		v = m->s_s();
		v->copy(&vc);
		vc.canonicize_map(G, &aut, &transporter, f_vv);
		aut.group_order(&stab_j);
		if (f_vv) {
			printf(" stab order ");
			stab_j.println();
			}
		a.mult(&stab_j, &b);
		b.ganzdiv_integral(&stab_i, &c);
		c.copy(m->s_k());

		q = q->s_n();
		}
	if (f_v) {
		printf("res1=\n");
		res1->Print();
		}
	res1->swap(res);
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::nabla(POLYNOM_OP res, INT n, INT f_v)
#endif
{
	VECTOR_OB V;
	INT n2, i, j, k, ki, kj, ij, m, m1, m2, m3;
	POLYNOM_OP tmp, res1 = NIL;
	INTEGER_OB c;
	
	n2 = (n * (n - 1)) >> 1;
	// res->b_skn(NIL, NIL, NIL);
	// res->s_s()->m_il_n(n2);
	// ((INTEGER_OP) res->s_k())->m_i(0);
	
	for (i = 0; i < n2; i++) {
		m = s_ii(i);
		if (m < 2)
			continue;
		c.m_i(4 * m * (m - 1));
		copy(&V);
		V.m_ii(i, m - 1);
		tmp = (POLYNOM_OP) callocobject("VEC::nabla() tmp");
		tmp->m_skn(&V, &c, NIL);
		if (res1 == NIL) {
			res1 = tmp;
			}
		else
			((LIST_OP) res1)->insert_list((LIST_OP) tmp, add_koeff, comp_monomvector_monomvector);
		}

	for (i = 0; i < n2; i++) {
		m = s_ii(i);
		if (m < 1)
			continue;
		c.m_i(2 * m);
		copy(&V);
		V.m_ii(i, m - 1);
		tmp = (POLYNOM_OP) callocobject("VEC::nabla() tmp");
		tmp->m_skn(&V, &c, NIL);
		if (res1 == NIL) {
			res1 = tmp;
			}
		else
			((LIST_OP) res1)->insert_list((LIST_OP) tmp, add_koeff, comp_monomvector_monomvector);
		}

	for (k = 0; k < n; k++) {
		for (i = 0; i < n; i++) {
			if (i == k)
				continue;
			ki = ij2k(k, i, n);
			m1 = s_ii(ki);
			if (m1 < 1)
				continue;
			for (j = i + 1; j < n; j++) {
				if (j == k)
					continue;
				kj = ij2k(k, j, n);
				m2 = s_ii(kj);
				if (m2 < 1)
					continue;
				ij = ij2k(i, j, n);
				m3 = s_ii(ij);
				
				c.m_i(2 * m1 * m2);
				copy(&V);
				V.m_ii(kj, m2 - 1);
				tmp = (POLYNOM_OP) callocobject("VEC::nabla() tmp");
				tmp->m_skn(&V, &c, NIL);
				if (res1 == NIL) {
					res1 = tmp;
					}
				else
					((LIST_OP) res1)->insert_list((LIST_OP) tmp, add_koeff, comp_monomvector_monomvector);

				c.m_i(2 * m1 * m2);
				copy(&V);
				V.m_ii(ki, m1 - 1);
				tmp = (POLYNOM_OP) callocobject("VEC::nabla() tmp");
				tmp->m_skn(&V, &c, NIL);
				if (res1 == NIL) {
					res1 = tmp;
					}
				else
					((LIST_OP) res1)->insert_list((LIST_OP) tmp, add_koeff, comp_monomvector_monomvector);

				c.m_i(-2 * m1 * m2);
				copy(&V);
				V.m_ii(ki, m1 - 1);
				V.m_ii(kj, m2 - 1);
				V.m_ii(ij, m3 + 1);
				tmp = (POLYNOM_OP) callocobject("VEC::nabla() tmp");
				tmp->m_skn(&V, &c, NIL);
				if (res1 == NIL) {
					res1 = tmp;
					}
				else
					((LIST_OP) res1)->insert_list((LIST_OP) tmp, add_koeff, comp_monomvector_monomvector);

				}
			}
		}
	res1->swap(res);
	res1->freeself();
	printf("nabla() res=");
	res->println();
	return OK;
}

#if TEXDOCU
INT VECTOR_OB::compare_images_as_unordered_multisets(VECTOR_OP V2)
#endif
{
	VECTOR_OB S1, S2;
	
	copy(&S1);
	V2->copy(&S2);
	S1.quicksort(S1.s_li(), TRUE);
	S2.quicksort(S2.s_li(), TRUE);
	return S1.compare(&S2);
}


#endif /* VECTORTRUE */



