/* ma.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef MATRIXTRUE

#include <DISCRETA/ma.h>

#ifdef PERMTRUE
#include <DISCRETA/perm.h>
#endif
#include <DISCRETA/gfq.h> // for perm_rep_mtx
#include <DISCRETA/part.h>
#include <DISCRETA/in.h> // for stirling_first / stirling_second

/*
 * MA
 */

static struct sym_matrix *callocmatrix(void);

#if TEXDOCU
INT MATRIX_OB::freeself()
#endif
{
	OBJECTSELF d;
	INT k, erg = OK;
	SYM_OP z;

	if (s_obj_k() == EMPTY) {
		//printf("MA::freeself(): s_obj_k() is EMPTY, nothing to free !");
		//fflush(stdout);
		return OK;
		}
	if (s_obj_bk() != MATRIX) {
		printf("warning: MA::freeself(): s_obj_bk() is not MATRIX, calling global free !");
		fflush(stdout);
		((SYM_OP) this)->freeself();
		return OK;
		}
	d = s_obj_s();
	z = s_ij(s_hi() - 1, s_li() - 1);
	k = s_hi() * s_li();
	for (; k > 0; k--, z--) {
		if (! z->emptyp())
			if (z->s_obj_k() != INTEGER)
				erg += z->freeself();
		}
	my_free(s_s());
	freeall(s_l()); 
	freeall(s_h()); 
	my_free(d.ob_matrix);
	c_obj_k(EMPTY);
	return erg;
}

#if TEXDOCU
INT MATRIX_OB::m_lh(INTEGER_OP len, INTEGER_OP height)
#else
make\_length\_height\_matrix: allocates a matrix of dimension height $\times$ length.
height and len will be copied.
#endif
{
	return m_ilih(len->s_i(), height->s_i());
}

#if TEXDOCU
INT MATRIX_OB::m_ilih(INT len, INT height)
#else
make\_length\_height\_matrix: allocates a matrix of dimension height $\times$ length.
length and height are specified as integers.
#endif
{
	INTEGER_OP len1 = (INTEGER_OP)callocobject("MA::m_ilih");
	INTEGER_OP height1 = (INTEGER_OP)callocobject("MA::m_ilih");
	
	len1->m_i(len);
	height1->m_i(height);
	return b_lh(len1, height1);
}

#if TEXDOCU
INT MATRIX_OB::m_lh_n(INTEGER_OP len, INTEGER_OP height)
#else
make\_length\_height\_null\_matrix: allocates a matrix of dimension height $\times$ length.
height and len will be copied.
The matrix is initialized with zero.
#endif
{
	INTEGER_OP len1 = (INTEGER_OP)callocobject("MA::m_lh_n");
	INTEGER_OP height1 = (INTEGER_OP)callocobject("MA::m_lh_n");
	
	len1->m_i(len->s_i());
	height1->m_i(height->s_i());
	return b_lh_n(len1, height1);
}

#if TEXDOCU
INT MATRIX_OB::m_ilih_n(INT len, INT height)
#else
make\_integer\_length\_integer\_height\_null\_matrix: 
allocates a matrix of dimension height $\times$ length.
height and len are specified as integers.
The matrix is initialized with zero.
#endif
{
	INTEGER_OP len1 = (INTEGER_OP)callocobject("MA::m_ilih_n");
	INTEGER_OP height1 = (INTEGER_OP)callocobject("MA::m_ilih_n");
	
	len1->m_i(len);
	height1->m_i(height);
	return b_lh_n(len1, height1);
}

#if TEXDOCU
INT MATRIX_OB::b_lhs(INTEGER_OP len, INTEGER_OP height, SYM_OP self)
#else
bulid\_length\_height\_self\_matrix: initializes a matrix using the 
objects len and height and self. These objects become part of the 
matrix and may not be freed anywhere else.
#endif
{
	OBJECTSELF d;

	d.ob_matrix = callocmatrix();
	if (d.ob_matrix == NIL)
		return ERROR;
	freeself();
	c_obj_k(MATRIX);
	c_obj_s(d);
	c_l(len);
	c_h(height);
	c_s(self);
	return(OK);
}

#if TEXDOCU
INT MATRIX_OB::b_lh(INTEGER_OP l, INTEGER_OP h)
#else
build\_length\_height\_matrix: builds a matrix using the l and h object for 
the dimensions. The self part is allocated.
#endif
/* build_length_height_matrix */
/* height und length werden nicht kopiert */
{
	SYM_OP s;
	INT i, len;
	
	len = l->s_i() * h->s_i();
	s = (SYM_OP) my_malloc(len * sizeof(struct object), "MA::b_lh");
	for (i = 0; i < len; i++)
		(s + i)->c_obj_k(EMPTY);
	return b_lhs(l, h, s);
}

#if TEXDOCU
INT MATRIX_OB::b_lh_n(INTEGER_OP l, INTEGER_OP h)
#else
build\_length\_height\_null\_matrix: build a matrix and initializes all 
entries to zero (integer zero).
#endif
/* mit 0 vorbesetzen */
/* build_length_height_null_matrix */
{
	INT i, erg = OK;
	SYM_OP z;
	
	erg += b_lh(l, h);
	if (erg)
		return erg;
	for (z = s_s(), i = s_hi() * s_li(); i > 0; i--, z++)
		((INTEGER_OP)z)->m_i(0);
	return OK;
}

static struct sym_matrix * callocmatrix(void)
{
	struct sym_matrix *ergebnis;

	ergebnis = (struct sym_matrix *) my_malloc(sizeof(struct sym_matrix), "callocmatrix()");
	if (ergebnis == NULL)
		error("callocmatrix:no mem");
	return(ergebnis);
}

INT MATRIX_OB::sprint(BYTE *str)
/* appends to str. 
 * writes to maximal strlength of 200. */
{
	INT i, j;
	BYTE str1[256];
	
	strcat(str, "(");
	for (i = 0; i < s_hi(); i++) {
		strcat(str, "(");
		for (j = 0; j < s_li(); j++) {
			str1[0] = 0;
			s_ij(i, j)->sprint(str1);
			if (strlen(str) + strlen(str1) < 200)
				strcat(str, str1);
			else
				return OK;
			if (j < s_li() - 1)
				strcat(str, ":");
			}
		strcat(str, ")");
		if (i < s_hi() - 1)
			strcat(str, ", ");
		}
	strcat(str, ")");
	return(OK);
}

INT MATRIX_OB::Print()
{
	fPrint(stdout);
	return OK;
}

INT MATRIX_OB::fPrint(FILE *fp)
{
	INT i, j;
	BYTE str1[256];
	
	if (s_obj_k() == EMPTY) {
		printf("MA::Print(): EMPTY object");
		fflush(stdout);
		return OK;
		}
	fprintf(fp, "(");
	for (i = 0; i < s_hi(); i++) {
		fprintf(fp, "(");
		for (j = 0; j < s_li(); j++) {
			str1[0] = 0;
			s_ij(i, j)->sprint(str1);
			fprintf(fp, "%s", str1);
			if (j < s_li() - 1)
				fprintf(fp, ":");
			}
		fprintf(fp, ")\n");
		}
	fprintf(fp, ")\n");
	return(OK);
}

INT MATRIX_OB::fprint_GAP(FILE *fp)
{
	INT i, j;

	if (s_obj_k() == EMPTY) {
		printf("MA::Print(): EMPTY object");
		fflush(stdout);
		return OK;
		}
	fprintf(fp, "[\n");
	for (i = 0; i < s_hi(); i++) {
		fprintf(fp, "[");
		for (j = 0; j < s_li(); j++) {
			s_ij(i, j)->fprint_GAP(fp);
			if (j < s_li() - 1)
				fprintf(fp, ", ");
			}
		fprintf(fp, "]");
		if (i < s_hi() - 1)
			fprintf(fp, ", \n");
		fprintf(fp, "\n");
		}
	fprintf(fp, "]\n");
	return(OK);
}

#if TEXDOCU
INT MATRIX_OB::calc_column_width(VECTOR_OP col_width, INT f_latex)
#endif
{
	INT i, j, h, l, c;
	BYTE str1[1000];
	
	h = s_hi();
	l = s_li();
	col_width->m_il_n(l);
	for (j = 0; j < l; j++) {
		for (i = 0; i < h; i++) {
			str1[0] = 0;
			if (f_latex) {
				s_ij(i, j)->sprint_latex(str1);
				}
			else {
				s_ij(i, j)->sprint(str1);
				}
			c = col_width->s_ii(j);
			c = MAXIMUM(c, (INT) strlen(str1));
			col_width->m_ii(j, c);
			}
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::fprint_raw(FILE *f)
#else
Prints the matrix to FILE. The colums will be aligned using space characters.
#endif
{
	INT *cl = NIL; /* column length */
	INT i, j;
	BYTE str1[1000];
	BYTE str2[1000];
	INT f_verbose = FALSE;
	
	/* if (s_li() > 500)
		f_verbose = TRUE; */
	cl = (INT *) my_malloc(s_li() * sizeof(INT), "MA::fprint_raw()");
	if (cl == NIL)
		return error("MA::frint_raw(): no memory");
	if (f_verbose) {
		printf("MA::fprint_raw(calculating col width:");
		fflush(stdout);
		}
	for (j = 0; j < s_li(); j++) {
		if (f_verbose) {
			printf(".");
			if (((j + 1) % 25) == 0)
				printf("%ld ", j + 1);
			if (((j + 1) % 50) == 0)
				printf("\n");
			fflush(stdout);
			}
		cl[j] = 0;
		for (i = 0; i < s_hi(); i++) {
			str1[0] = 0;
			s_ij(i, j)->sprint(str1);
			cl[j] = MAX(cl[j], (INT) strlen(str1));
			}
		}
	if (f_verbose) {
		printf(";writing to file:");
		fflush(stdout);
		}
	for (i = 0; i < s_hi(); i++) {
		for (j = 0; j < s_li(); j++) {
			str1[0] = 0;
			s_ij(i, j)->sprint(str1);
			str2[0] = 0;
			strfill(str2, cl[j] - strlen(str1), ' ');
			strcat(str2, str1);
			fprintf(f, "%s ", str2);
			/* if (j < s_li() - 1)
				fprintf(f, " "); */
			}
		fprintf(f, "\n");
		}
	fprintf(f, "\n");
	if (f_verbose) {
		printf(")\n");
		fflush(stdout);
		}
	if (cl) {
		my_free(cl);
		cl = NIL;
		}
	return(OK);
}

#if TEXDOCU
INT MATRIX_OB::sprint_latex(BYTE *str)
#else
Prints the matrix in latex format to str.
#endif
{
	BYTE str0[10000];
	BYTE str1[10000];
	BYTE str2[10000];
	VECTOR_OB col_width;
	INT i, j, l, h;
	
	l = s_li();
	h = s_hi();
	calc_column_width(&col_width, TRUE /* f_latex */);
	sprintf(str0, "\\begin{array}{*{%ld}{c}}\n", l);
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			str1[0] = 0;
			s_ij(i, j)->sprint_latex(str1);
			str2[0] = 0;
			strfill(str2, col_width.s_ii(j) - strlen(str1), ' ');
			strcat(str2, str1);
			sprintf(str0 + strlen(str0), "%s ", str2);
			if (j < l - 1)
				sprintf(str0 + strlen(str0), " & ");
			}
		sprintf(str0 + strlen(str0), "\\\\\n");
		}
	sprintf(str0 + strlen(str0), "\\end{array}\n");
	strcat(str, str0);
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::latex(FILE *fp)
#else
Prints the matrix in latex format to FILE fp.
#endif
{
	VECTOR_OB col_width;
	INT i, j, l, h;
	BYTE str1[10000];
	BYTE str2[10000];
	
	l = s_li();
	h = s_hi();
	calc_column_width(&col_width, TRUE);
	fprintf(fp, "\\begin{array}{*{%ld}{r}}\n", l);
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			str1[0] = 0;
			s_ij(i, j)->sprint_latex(str1);
			str2[0] = 0;
			strfill(str2, col_width.s_ii(j) - strlen(str1), ' ');
			strcat(str2, str1);
			fprintf(fp, "%s ", str2);
			if (j < l - 1)
				fprintf(fp, " & ");
			}
		fprintf(fp, "\\\\\n");
		}
	fprintf(fp, "\\end{array}\n");
	return OK;
#if 0
	
	l = s_li();
	h = s_hi();
	// fprintf(fp, "\\left(\n");
	fprintf(fp, "\\begin{array}{*{%ld}{c}}\n", l);
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			s_ij(i, j)->latex(fp);
			if (j < l - 1)
				fprintf(fp, " & ");
			}
		fprintf(fp, " \\\\\n");
		}
	fprintf(fp, "\\end{array}\n");
	// fprintf(fp, "\\right)\n");
	return OK;
#endif
}

#if TEXDOCU
INT MATRIX_OB::sprint_latex_decomposed_vector(BYTE *str, 
	VECTOR_OP p, VECTOR_OP q)
#else
Prints the matrix in latex format to str.
#endif
{
	PARTITION_OB pp, qq;
	INT i, l;
	
	l = p->s_li();
	pp.m_kli(VECTOR, l);
	for (i = 0; i < l; i++) {
		pp.m_ii(i, p->s_ii(i));
		}
	l = q->s_li();
	qq.m_kli(VECTOR, l);
	for (i = 0; i < l; i++) {
		qq.m_ii(i, q->s_ii(i));
		}
	return sprint_latex_decomposed(str, &pp, &qq);
}

#if TEXDOCU
INT MATRIX_OB::sprint_latex_decomposed(BYTE *str, 
	PARTITION_OP p, PARTITION_OP q)
#else
Prints the matrix in latex format to str.
#endif
{
	BYTE str0[10000];
	BYTE str1[10000];
	BYTE str2[10000];
	VECTOR_OB col_width, w;
	INT i, j, l, h, ll, a, b;
	
	l = s_li();
	h = s_hi();
	calc_column_width(&col_width, TRUE /* f_latex */);
	sprintf(str0, "\\begin{array}{|");
	ll = q->s_li();
	b = 0;
	for (i = 0; i < ll; i++) {
		a = q->s_ii(i);
		sprintf(str0 + strlen(str0), "*{%ld}{c}|", a);
		b += a;
		}
	sprintf(str0 + strlen(str0), "}\n");
	if (b != l) {
		printf("WARNING: MA::sprint_latex_decomposed decomposition q does not match length\n");
		}
	w.m_il_n(h);
	ll = p->s_li();
	b = 0;
	for (i = 0; i < ll; i++) {
		a = p->s_ii(i);
		w.m_ii(b, 1);
		b += a;
		}
	if (b != h) {
		printf("WARNING: MA::sprint_latex_decomposed decomposition p does not match height\n");
		}
	
	for (i = 0; i < h; i++) {
		if (w.s_ii(i))
			sprintf(str0 + strlen(str0), "\\hline\n");
		for (j = 0; j < l; j++) {
			str1[0] = 0;
			s_ij(i, j)->sprint_latex(str1);
			str2[0] = 0;
			strfill(str2, col_width.s_ii(j) - strlen(str1), ' ');
			strcat(str2, str1);
			sprintf(str0 + strlen(str0), "%s ", str2);
			if (j < l - 1)
				sprintf(str0 + strlen(str0), " & ");
			}
		sprintf(str0 + strlen(str0), "\\\\\n");
		}
	sprintf(str0 + strlen(str0), "\\hline\n");
	sprintf(str0 + strlen(str0), "\\end{array}\n");
	strcat(str, str0);
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::latex_lower_tri(FILE *fp)
#else
Prints the matrix in latex format to FILE fp.
#endif
{
	INT i, j, l, h;
	
	l = s_li();
	h = s_hi();
	fprintf(fp, "\\left( \\begin{array}{*{%ld}{c}}\n", l);
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			if (j <= i)
				s_ij(i, j)->latex(fp);
			if (j < l - 1)
				fprintf(fp, " & ");
			}
		fprintf(fp, " \\\\\n");
		}
	fprintf(fp, "\\end{array} \\right)\n");
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::latex_upper_tri(FILE *fp)
#else
Prints the matrix in latex format to FILE fp.
#endif
{
	INT i, j, l, h;
	
	l = s_li();
	h = s_hi();
	fprintf(fp, "\\left( \\begin{array}{*{%ld}{c}}\n", l);
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			if (j >= i)
				s_ij(i, j)->latex(fp);
			if (j < l - 1)
				fprintf(fp, " & ");
			}
		fprintf(fp, " \\\\\n");
		}
	fprintf(fp, "\\end{array} \\right)\n");
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::latex_upper_tri_block_wise(FILE *fp, INT block_size)
#else
Prints the matrix in latex format to FILE fp.
#endif
{
	MATRIX_OB A;
	INT i, j, l, h, ll, hh, ii, jj;
	
	l = s_li();
	h = s_hi();
	if (l != h)
		return error("MATRIX_OB::latex_upper_tri_block_wise() l != h");
	for (i = 0; i < l; i += block_size) {
		hh = MINIMUM(i + block_size, h) - i;
		for (j = i; j < l; j += block_size) {
			ll = MINIMUM(j + block_size, l) - j;
			A.m_ilih(ll, hh);
			for (ii = 0; ii < hh; ii++) {
				for (jj = 0; jj < ll; jj++) {
					s_ij(i + ii, j + jj)->copy(A.s_ij(ii, jj));
					}
				}
			fprintf(fp, "\\[\n");
			fprintf(fp, "(%ld,%ld)=", i, j);
			A.latex(fp);
			fprintf(fp, "\\]\n");
			}
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::integer_print_dense()
#endif
{
	INT i, i0, i1, j, m, n, a, idx, f_found;
	BYTE alphabet[100];
	BYTE extra_alphabet[100];
	VECTOR_OB large_values;
	INTEGER_OB large_values_len;
	
	for (i = 0; i < 10; i++)
		alphabet[i] = '0' + i;
	i0 = i;
	for (i = 0; i < 26; i++) {
		alphabet[i0 + i] = 'a' + i;
		}
	i0 += i;
	alphabet[i0] = 0;
	printf("alphabet:\n%s\n", alphabet);

	for (i = 0; i < 26; i++) {
		extra_alphabet[i] = 'A' + i;
		}
	i1 = i;
	extra_alphabet[i1] = 0;
	printf("extra_alphabet:\n%s\n", extra_alphabet);
	
	m = s_hi();
	n = s_li();
	large_values.m_il(VECTOR_OVERSIZE);
	large_values_len.m_i(0);
	
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			a = s_iji(i, j);
			if (a < 0 || a >= i0) {
				printf("warning: entry out of range\n");
				printf("a=%ld at i=%ld j=%ld\n", 
					a, i, j);
				large_values.search(large_values_len.s_i(), 
					TRUE /* f_ascending */, s_ij(i, j), 
					&idx, &f_found);
				if (!f_found) {
					large_values.insert_at(&large_values_len, idx, s_ij(i, j));
					}
				// return OK;
				}
			}
		}
	if (large_values_len.s_i()) {
		if (large_values_len.s_i() > i1) {
			printf("too many large values, cannot print\n");
			return OK;
			}
		printf("found the following large values:\n");
		for (i = 0; i < large_values_len.s_i(); i++) {
			large_values.s_i(i)->print();
			printf(" becomes %c\n", extra_alphabet[i]);
			}
		}
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			a = s_iji(i, j);
			if (a < 0 || a >= i0) {
				large_values.search(large_values_len.s_i(), 
					TRUE /* f_ascending */, s_ij(i, j), 
					&idx, &f_found);
				if (!f_found)
					return error("large value not found");
				printf("%c", extra_alphabet[idx - 1]);
				}
			else {
				printf("%c", alphabet[a]);
				}
			}
		printf("\\\\\n");
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::nullp()
#else
returns TRUE if the matrix is zero.
#endif
{
	INT i, j;
	for (i = 0; i < s_hi(); i++)
		for (j = 0; j < s_li(); j++)
			if (! s_ij(i, j)->nullp())
				return FALSE;
	return TRUE;
}

#if TEXDOCU
INT MATRIX_OB::einsp()
#else
returns TRUE if the matrix is quadratic and an identity matrix.
#endif
{
	INT i, j;
	
	if (s_hi() != s_li()) return FALSE;
#if 0
	if (norm_projective() != OK)
		return error("MA::einsp()|"
		"error in norm_projective()");
#endif
	for (i = 0; i < s_hi(); i++)
		for (j = 0; j < s_li(); j++)
			if (i == j) {
				if (! s_ij(i, j)->einsp())
					return FALSE;
				}
			else {
				if (! s_ij(i, j)->nullp())
					return FALSE;
				}
	return TRUE;
}

#if TEXDOCU
INT MATRIX_OB::compare(MATRIX_OP b)
#else
returns $-1 \iff this < b$, $0 \iff this = b$,  $1 \iff this > b$.
#endif
{
	INT i, j, erg;
	SYM_OP x, y;
	
	x = s_s();
	y = b->s_s();	
	for (i = 0; i < s_hi(); i++) {
		if (i >= b->s_hi())
			return 1;
		else {
			for (j = 0; j < s_li(); j++) {
				if (j >= b->s_li())
					return 1;
				else {
					erg = x->sym_comp(y);
					if (erg != 0)
						return(erg);
					x++; y++;
					}
				}
			}
		}
	if (b->s_hi() > s_hi()) return(-1);
	return(0); /* matrizen sind gleich */
}

#if TEXDOCU
INT MATRIX_OB::add_apply_matrix(MATRIX_OP b)
#else
Adds the matrices in this and b. The result is written into b.
#endif
{
	SYM_OP c, d;
	INT i, erg = OK;
	
	if ((s_hi() == b->s_hi()) &&
		(s_li() == b->s_li())) {
		i = s_hi() * s_li();
		c = s_s();
		d = b->s_s();
		while (i-- > 0) {
			erg += c->add_apply(d);
			c++;
			d++;
			}
		}
	else {
		return error("MA::add_apply_matrix()|different format");
		
		}
#if 0
	if (b->norm_projective() != OK)
		return error("MA::add_apply_matrix()|error in norm_projective()");
#endif
	return erg;
}

INT MATRIX_OB::add_apply(SYM_OP b)
{
	if (b->emptyp()) 
		return(copy((MATRIX_OP)b));
	switch (b->s_obj_k()) {
		case KRANZTYPUS:
		case MATRIX: 
			return(add_apply_matrix(
				(MATRIX_OP)b));
		default:
			b->printobjectkind();
			error("MA::add_apply:wrong second type");
			return(ERROR);
	}
}

INT MATRIX_OB::add(MATRIX_OP b, MATRIX_OP ergeb)
{
	INT i, j;
	INT len, height;
	SYM_OP z;
	
	len = MAXIMUM(s_li(), b->s_li());
	height = MAXIMUM(s_hi(), b->s_hi());

	if (!ergeb->emptyp())
		ergeb->freeself();
	ergeb->m_ilih(len, height);
	ergeb->c_obj_k(s_obj_k());

	z = ergeb->s_s();
	for (i = 0; i < height; i++)
		for (j = 0; j < len; j++, z++) {
			if ((i < s_hi()) && (i < b->s_hi()) &&
				(j < s_li()) && (j < b->s_li())) { 
				s_ij(i, j)->add(b->s_ij(i, j), z); 
				}
			else if ((i < s_hi()) && (j < s_li()))
				s_ij(i, j)->copy(z); 
			else if ((i < b->s_hi()) && (j < b->s_li()))
				b->s_ij(i, j)->copy(z); 
			}
#if 0
	if (ergeb->norm_projective() != OK)
		return error("MA::add()|"
		"error in norm_projective()");
#endif
	return(OK);
}

INT MATRIX_OB::copy(MATRIX_OP b)
{
	INT k;
	SYM_OP z, w;

	b->m_ilih(s_li(), s_hi());
	b->c_obj_k(s_obj_k());
	z = s_ij(s_hi() - 1, s_li() - 1);
	w = b->s_ij(s_hi() - 1, s_li() - 1);
	k = s_hi() * s_li();
	for (; k > 0; k--, z--, w--) { 
		if (z->s_obj_k() == INTEGER)
				*w = *z;
			else if (z->emptyp())
				*w = *z;
			else
				z->copy(w);
		}
	return(OK);
}

#if TEXDOCU
INT MATRIX_OB::mult_apply_scalar_matrix(SYM_OP a)
#else
Multiplies the matrix in this elementwise by the scalar a.
The result is written into this again.
#endif
/* before: mult_apply_scalar_matrix(
 *    OP a, OP b), a is scalar, b matrix */
{
	SYM_OP z = s_s();
	INT grenze = s_li() * s_hi();
	INT i;
	
	for (i = 0; i < grenze; i++, z++)
		a->mult_apply(z);
	return(OK);
}

#if TEXDOCU
INT MATRIX_OB::mult_scalar(SYM_OP scalar, MATRIX_OP res)
#else
Multiplies the matrix in this elementwise by the scalar a.
The result is written into res.
#endif
/* before: mult_scalar_matrix(
 * OP scalar, OP mat, OP res) 
 * ACHTUNG: mat is now this ! */
/* res != mat notwendig ! */
/* komponentenweise multiplikation */
{
	INT i, j, erg = OK;
	INT height, len;

	if (s_s() == NULL)
		return error("MA::mult_scalar:self == NULL");
	height = s_hi();
	len = s_li();
	erg += res->m_ilih(len, height);
	res->c_obj_k(s_obj_k());
	for (i = 0; i < height; i++)
		for (j = 0; j < len; j++)
			erg += scalar->mult(
				s_ij(i, j), res->s_ij(i, j));
	return erg;
}

#if TEXDOCU
INT MATRIX_OB::mult_matrix(MATRIX_OP b, MATRIX_OP c)
#else
Multiplies the matrices this and b. Writes the result into c.
c must be an object different from this and b. 
#endif
/* c = a * b; c must be different from a and b. */
/* before: mult_matrix_matrix(OP a, OP b, OP c) */
{
	INT i, j, k;
	SYM_OP z;
	
	if (s_li() != b->s_hi())
		return error("MA::mult_matrix: "
		"cannot be multiplied");
	if (c == b || c == this)
		return error("MA::mult_matrix: "
		"c == b || c == this");
	c->m_ilih(b->s_li(), s_hi());
	c->c_obj_k(s_obj_k());
	z = (SYM_OP) callocobject("MA::mult_matrix");
	for (i = 0; i < s_hi(); i++)
		for (j = 0; j < b->s_li(); j++)
			for (k = 0; k < s_li(); k++) { 
				s_ij(i, k)->mult(b->s_ij(k, j), z);
				z->add(c->s_ij(i, j), c->s_ij(i, j));
			}
	freeall(z); 
#if 0
	if (c->norm_projective() != OK)
		return error("MA::mult_matrix()|"
		"error in norm_projective()");
#endif
	return(OK);
}

#if TEXDOCU
INT MATRIX_OB::mult_vector(VECTOR_OP b, VECTOR_OP c)
#else
Multiplies the matrix this and b. Writes the result into c.
c must be an object different from this and b. 
#endif
{
	INT i, k, m, n;
	SYM_OB z, z1;
	
	m = s_hi();
	n = s_li();
	if (n != b->s_li())
		return error("MA::mult_vector: cannot be multiplied");
	if (c == b)
		return error("MA::mult_vector: c == b");
	c->m_il(m);
	for (i = 0; i < m; i++) {
		for (k = 0; k < s_li(); k++) { 
			s_ij(i, k)->mult(b->s_i(k), &z);
			if (k == 0)
				z.swap(&z1);
			else
				z.add_apply(&z1);
			}
		z1.swap(c->s_i(i));
		}
	return(OK);
}

#if TEXDOCU
INT MATRIX_OB::mult_apply_matrix(MATRIX_OP b)
#else
Multiplies the matrices this and b and writes the result into b again.
#endif
/* b = a * b */
/* before: mult_apply_matrix_matrix(OP a, OP b) */
{
	MATRIX_OP c = (MATRIX_OP) callocobject("MA::mult_apply_matrix");
	INT erg = OK;
	
	*c = *b;
	b->c_obj_k(EMPTY);
	erg += mult_matrix(c, b);
	erg += freeall(c);
	return erg;
}

#if TEXDOCU
INT MATRIX_OB::mult(SYM_OP b, MATRIX_OP d)
#else
Multiplies the matrix in this by something in b and writes the result into d.
Something in b can be either a matrix or a scalar. 
#endif
/* before: mult_matrix(OP a, OP b, OP d) */
{
	switch (b->s_obj_k()) {
	case BRUCH:
	case INTEGER:
	case LONGINT: 
		return mult_scalar(b, d);
	case MATRIX: 
		return mult_matrix((MATRIX_OP)b, d);
	/* case VECTOR: 
		return mult_matrix_vector(a,b,d); !!! */
	default:
		{
			b->printobjectkind();
			error("MA::mult:wrong second type"); 
			return(ERROR);
		}
	}
}

#if TEXDOCU
INT MATRIX_OB::mult_apply(SYM_OP b)
#else
Multiplies this and b and writes the result into b.
b can be either a scalar or a matrix.
#endif
/* before: mult_apply_matrix(OP a, OP b) */
{
	switch (b->s_obj_k()) {
	case MATRIX: 
		return(mult_apply_matrix((MATRIX_OP)b));
	default:
		{
			b->printobjectkind();
			error("mult_apply:wrong second type"); 
			return(ERROR);
		}
	}
}

#if TEXDOCU
INT MATRIX_OB::swap_col(INT i, INT j)
#else
swaps the entries in columns $i$ and $j$.
#endif
/* before: change_column_ij(OP a, INT i, INT j) */
{
	INT k;
	
	if (s_obj_k() != MATRIX) 
		error("MA::swap_col: no matrix");
	if ((i < 0) || (i >= s_li())) 
		error("MA::swap_col: wrong i");
	if ((j < 0) || (j >= s_li())) 
		error("MA::swap_col: wrong j");

	for (k = 0; k < s_hi(); k++)
		s_ij(k, i)->swap(s_ij(k, j));
	return(OK);
}

#if TEXDOCU
INT MATRIX_OB::swap_row(INT i, INT j)
#else
swaps the entries in rows $i$ and $j$.
#endif
/* before: change_row_ij(OP a, INT i, INT j) */
{
	INT k;

	if (s_obj_k() != MATRIX) 
		error("MA::swap_row: no matrix");
	if ((i < 0) || (i >= s_hi())) 
		error("MA::swap_row: wrong i");
	if ((j < 0) || (j >= s_hi())) 
		error("MA::swap_row: wrong j");

	for (k = 0; k < s_li(); k++)
		s_ij(i, k)->swap(s_ij(j, k));
	return(OK);
}

#if TEXDOCU
INT MATRIX_OB::zero()
#else
Sets all entries to zero.
#endif
/* l/h muessen schon gesetzt sein. */
{
	INT i, j, erg = OK;
	
	for (i = 0; i < s_hi(); i++) {
		for (j = 0; j < s_li(); j++) {
			erg += s_ij(i, j)->zero();
			}
		}
	return erg;
}

#if TEXDOCU
INT MATRIX_OB::one()
#else
Sets all entries $x_{i,i}$ to one, all other entries to zero.
#endif
/* l/h muessen schon gesetzt sein. */
{
	INT i, j, erg = OK;
	
	for (i = 0; i < s_hi(); i++) {
		for (j = 0; j < s_li(); j++) {
			if (i == j)
				erg += s_ij(i, j)->one();
			else
				erg += s_ij(i, j)->zero();
			}
		}
	return erg;
}

#if TEXDOCU
INT MATRIX_OB::m_one()
#else
Sets all entries $x_{i,i}$ to $-1$, all other entries to zero.
#endif
/* l/h muessen schon gesetzt sein. */
{
	INT i, j, erg = OK;
	
	for (i = 0; i < s_hi(); i++) {
		for (j = 0; j < s_li(); j++) {
			if (i == j)
				erg += s_ij(i, j)->m_one();
			else
				erg += s_ij(i, j)->zero();
			}
		}
	return erg;
}

#if TEXDOCU
INT MATRIX_OB::homo_z(INT z)
#else
Sets all entries $x_{i,i}$ to homo\_z($z$), all other entries to zero.
#endif
/* l/h muessen schon gesetzt sein. */
{
	INT i, j, erg = OK;
	
	for (i = 0; i < s_hi(); i++) {
		for (j = 0; j < s_li(); j++) {
			if (i == j)
				erg += s_ij(i, j)->homo_z(z);
			else
				erg += s_ij(i, j)->zero();
			}
		}
	return erg;
}

#if TEXDOCU
INT MATRIX_OB::norm_projective()
#else
This function is no longer in use!
#endif
/* Normalform der projektiven Matrixgruppen:
 * Bestimme minimales j, so dass 
 * a[n-1][j] != 0 aber a[n-1][k] = 0
 * fuer k = 0, 1, .. j - 1.
 * Dann skalare Multiplikation mit a[n-1][j]^-1 */ 
{
	return error("in MA::norm_projective()");
	return(OK);
#if 0
	SYM_OP a_j, a_j_inv;
	MATRIX_OP q;
	INT j;
	
	if (s_obj_k() != MATRIX)
		return error("MA::norm_projective()|"
		"s_ind_k() != MATRIX");
	if (s_hi() != s_li()) {
		/* printf("MA::norm_projective()|"
		"not quadratic; do not normalize"); */
		return OK;
		}
	a_j_inv = (SYM_OP) callocobject("MA::norm_projective");
	q = (MATRIX_OP) callocobject("MA::norm_projective");
	for (j = 0; j < s_li(); j++) {
		a_j = s_ij(s_hi() - 1, j);
		if (!a_j->nullp())
			break;
		}
	if (j == s_li())
		return error("MA::norm_projective()|singulary matrix");
	if (!a_j->einsp()) {
		if (a_j->invers(a_j_inv) != OK)
			return error("MA::norm_projective()|a_j not invertible");
		swap(q);
		if (q->mult_scalar(a_j_inv, this) != OK)
			return error("MA::norm_projective()|error in MA::mult_scalar()");
		}
	freeall(a_j_inv);
	freeall(q);
	return(OK);
#endif
}

#if TEXDOCU
INT MATRIX_OB::invers(MATRIX_OP b)
#else
Computes the inverse of the matrix this into b.
Gauss-algorithm.
#endif
/* AK 290388 nach stoer (dietmar) */
/* umgewandelt aus pascal */
{
	INT erg = OK;
	INT i, j, k, r;
	/* r ist die selectierte spalte */
	/* r = 0 ... n */
	INT n = s_li() - 1;
	INT singulaer = FALSE;
	SYM_OB hr, hs;
	VECTOR_OB hv, p;

	if (! quadraticp())
		 return error("MA::invers: not quadratic"); 

#if 0
	if (norm_projective() != OK)
		return error("MA::invers()|"
		"error in norm_projective()");
#endif

	erg += p.m_il(n + 1);

	for (j = 0; j <= n; j++)
		p.m_ii(j, j);

	j = -1;
	erg += copy(b);
	while ((j++ < n) && (! singulaer)) {
		/* pivotsuche */
		for (r = j; r <= n; r++) {
			if (! b->s_ij(r, j)->nullp())
				break;
			}

		if (r == n + 1) /* nur nullen in der spalte j */
			singulaer = TRUE;
		else {
			/* zeilentausch */
			if (r > j) {
				for (k = 0; k <= n; k++) 
					erg += b->s_ij(j, k)->
						swap(b->s_ij(r, k));
				erg += p.s_i(j)->swap(p.s_i(r));
				}
			/* transformation */
			erg += b->s_ij(j, j)->invers(&hr);
			for (i = 0; i <= n; i++)
				erg += hr.mult_apply(b->s_ij(i, j));
			erg += hr.copy(b->s_ij(j, j));
			erg += hr.addinvers_apply();
			for (k = 0; k <= n; k++) {
				/* spalte j nicht anwenden: */
				if (k == j)
					k++;
				if (k > n)
					break;
				for (i = 0; i <= n; i++) {
					/* auf zeile j nicht anwenden: */
					if (i==j)
						i++; 
					if (i>n) break;
					erg += b->s_ij(i, j)->mult(
						b->s_ij(j, k), &hs);
					erg += hs.addinvers_apply();
					erg += hs.add_apply(b->s_ij(i, k));
					}
				erg += hr.mult_apply(b->s_ij(j, k));
				}
			} /* end else */
		} /* end while */
	if (erg != OK)
		return error("MA::invers: (1) error in computation");

	erg += hv.m_il(n + 1);
	if (! singulaer) {
		/* spaltentausch */
		for (i = 0; i <= n; i++) {
			for (k = 0; k <= n; k++) {
				erg += b->s_ij(i, k)->copy(hv.s_i(p.s_ii(k)));
				}
			for (k = 0; k <= n; k++) {
				erg += hv.s_i(k)->copy(b->s_ij(i, k));
				}
			}
		}

	if (singulaer)	{ 
		b->freeself();
		error("MA::invers: singulary");
		return(SINGULAER); 
	}
	if (erg != OK)
		return error("MA::invers: (2) "
		"error in computation");
#if 0
	if (b->norm_projective() != OK)
		return error("MA::invers()|"
		"error in b->norm_projective()");
#endif
	return erg;
}

#if TEXDOCU
INT MATRIX_OB::transpose(MATRIX_OP At)
#else
Computes the transposed matrix into $At$.
#endif
{
	INT l, h, i, j;
	
	if (this == At)
		return error("MA::transpose: this == a");
	l = s_li();
	h = s_hi();
	At->m_ilih(h, l);
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			s_ij(i, j)->copy(At->s_ij(j, i));
			}
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::Asup2Ainf(MATRIX_OP Ainf)
#else
Computes the Plesken matrix $A^\wedge$ (Ainf)  from $A^\vee$ (Asup).
Compare Plesken~\cite{Plesken82}.
#endif
{
	INT erg = OK;
	MATRIX_OB C;
	INTEGER_OB op1;
	INT l, h, i, j;
	INT *orbit_size;
	
	l = s_li();
	h = s_hi();
	if (l != h)
		return(ERROR);
	copy(&C);
	orbit_size = (INT *) my_malloc(l * sizeof(INT), "MA::Asup2Ainf");
	for (i = 0; i < l; i++) {
		orbit_size[i] = C.s_iji(i, l - 1);
		if (orbit_size[i] == 0)
			return(ERROR);
		}
	for (j = 0; j < l; j++) {
		op1.m_i(orbit_size[j]);
		for (i = 0; i < l; i++) {
			erg += op1.mult_apply(C.s_ij(i, j));
			}
		}
	for (i = 0; i < l; i++) {
		op1.m_i(orbit_size[i]);
		for (j = 0; j < l; j++) {
			erg += C.s_ij(i, j)->sym_div(&op1, C.s_ij(i, j));
			}
		}
	// erg += C.transpose(Ainf);
	C.swap(Ainf);
	my_free(orbit_size);
	return erg;
}

#if TEXDOCU
INT MATRIX_OB::Ainf2Asup(MATRIX_OP Asup)
#else
Computes the Plesken matrix $A^\vee$ (Asup)  from $A^\wedge$ (Ainf).
#endif
{
	INT erg = OK;
	MATRIX_OB C;
	INTEGER_OB op1;
	INT l, h, i, j;
	INT *orbit_size;
	
	l = s_li();
	h = s_hi();
	if (l != h)
		return(ERROR);
	copy(&C);
	orbit_size = (INT *) my_malloc(l * sizeof(INT), "MA::Asup2Ainf");
	for (i = 0; i < l; i++) {
		orbit_size[i] = C.s_iji(0, i);
		if (orbit_size[i] == 0)
			return(ERROR);
		}
	for (i = 0; i < l; i++) {
		op1.m_i(orbit_size[i]);
		for (j = 0; j < l; j++) {
			erg += op1.mult_apply(C.s_ij(i, j));
			}
		}
	for (j = 0; j < l; j++) {
		op1.m_i(orbit_size[j]);
		for (i = 0; i < l; i++) {
			erg += C.s_ij(i, j)->sym_div(&op1, C.s_ij(i, j));
			}
		}
	// erg += C.transpose(Asup);
	C.swap(Asup);
	my_free(orbit_size);
	return erg;
}

#if TEXDOCU
INT MATRIX_OB::Asup2Acover(MATRIX_OP Acover)
#else
Computes the cover relations of the poset defined by this (=Asup).
Assumes that Asup is upper triangular.
#endif
{
	INT l, h, i, j, k;
	
	l = s_li();
	h = s_hi();
	copy(Acover);
	/* replace the numbers by a 0/1 flag */
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			if (Acover->s_iji(i, j))
				Acover->m_iji(i, j, 1);
			}
		}
	for (i = h - 1; i >= 0; i--) {
		for (j = i + 1; j < l; j++) {
			if (Acover->s_iji(i, j)) {
				/* B_i < B_j in the orbit - 
				 * poset of the lattice */
				for (k = 0; k < i; k++) {
					if (Acover->s_iji(k, i)) {
						/* B_k < B_i < B_j
						 * therefore the entry 
						 * B_k < B_j can be deleted. */
						Acover->m_iji(k, j, 0);
						}
					}
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::Asup2Acover_arbitrary(MATRIX_OP Acover)
#else
Computes the cover relations of the poset defined by this (=Asup).
#endif
{
	INT l, h, i, j, k;
	
	l = s_li();
	h = s_hi();
	if (l != h)
		return error("MATRIX_OB::Asup2Acover_arbitrary() l != h");
	copy(Acover);
	/* replace the numbers by a 0/1 flag */
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			if (Acover->s_iji(i, j))
				Acover->m_iji(i, j, 1);
			}
		}
	for (i = 0; i < h; i++) {
		for (j = 0; j < h; j++) {
			if (j == i)
				continue;
			if (!s_iji(i, j))
				continue;
			for (k = 0; k < h; k++) {
				if (k == i)
					continue;
				if (k == j)
					continue;
				if (!s_iji(j, k))
					continue;
				if (s_iji(i, k)) {
					printf("%ld-%ld-%ld clearing(%ld,%ld)\n", i, j, k, i, k);
					Acover->m_iji(i, k, 0);
					}
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::Acover2nl(VECTOR_OP nl)
#else
Computes the \lq neighbour-list\rq of the poset whose cover relations are given 
in this (=Acover). This list nl is used as input for the lattice-placement 
program \lq vbp \rq.
#endif
{
	INT l, h, i, j, k, len, cur;
	
	l = s_li();
	h = s_hi();
	/* count the non diagonal entries: */
	k = 0;
	for (i = 0; i < h; i++) {
		for (j = i + 1; j < l; j++) {
			if (s_iji(i, j))
				k++;
			}
		}
	/* length nl = h + 1 + k */
	len = h + 1 + k;
	nl->m_il(len);
	nl->m_ii(h, len);
	cur = h + 1;
	for (i = 0; i < h; i++) {
		nl->m_ii(i, cur);
			/* start of neighbour list of i in nl */
		for (j = i + 1; j < l; j++) {
			if (s_iji(i, j)) {
				nl->m_ii(cur, j);
					/* found a neighbour of i */
				cur++;
				}
			}
		}
	return OK;
}

#ifdef PERMTRUE
#if TEXDOCU
INT MATRIX_OB::perm_cols_rows(PERMUTATION_OP p, PERMUTATION_OP q, MATRIX_OP M1)
#else
Permutes the rows by p and the columns by q and writes the result into M1.
#endif
{
	INT l, h, i, j, ii, jj;
	
	l = s_li();
	h = s_hi();
	M1->m_ilih(l, h);
	M1->c_obj_k(s_obj_k());
	for (i = 0; i < h; i++) {
		ii = p->s_ii(i) - 1;
		for (j = 0; j < l; j++) {
			jj = q->s_ii(j) - 1;
			s_ij(i, j)->copy(M1->s_ij(ii, jj));
			}
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::perm_cols(PERMUTATION_OP p, MATRIX_OP M1)
#else
Permutes the columns by q and writes the result into M1.
#endif
{
	INT l, h, i, j, j1;
	
	l = s_li();
	h = s_hi();
	M1->m_ilih(l, h);
	M1->c_obj_k(s_obj_k());
	for (j = 0; j < l; j++) {
		j1 = p->s_ii(j) - 1;
		for (i = 0; i < h; i++) {
			s_ij(i, j)->copy(M1->s_ij(i, j1));
			}
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::perm_rows(PERMUTATION_OP p, MATRIX_OP M1)
#else
Permutes the rows by p and writes the result into M1.
#endif
{
	INT l, h, i, j, i1;
	
	l = s_li();
	h = s_hi();
	M1->m_ilih(l, h);
	M1->c_obj_k(s_obj_k());
	for (i = 0; i < h; i++) {
		i1 = p->s_ii(i) - 1;
		for (j = 0; j < l; j++) {
			s_ij(i, j)->copy(M1->s_ij(i1, j));
			}
		}
	return OK;
}
#endif

#if TEXDOCU
INT MATRIX_OB::determinante(SYM_OP erg)
#else
computes $\det$. The result is written into erg. Gauss-algorithm.
#endif
/* before: det_mat_tri(OP a, OP erg) */
/* determinante einer matrix a mittels triangulation
 * algorithmus 41 in CACM */
{
	INT r, i, j, y, count, sign = 1, n;
	INT debug = 0;
	MATRIX_OP b = (MATRIX_OP) callocobject("MA::determinante");
	SYM_OP product = callocobject("MA::determinante");
	SYM_OP factor = callocobject("MA::determinante");
	SYM_OP temp = callocobject("MA::determinante");
	
	if (s_li() != s_hi())
		return error("MA::determinante(): l != h");
	n = s_li();
	
	s_ij(0, 0)->copy(product);
	product->one();
	
#if 0
	if (norm_projective() != OK)
		return error("MA::determinante()|"
		"error in MA::norm_projective()");
#endif

	copy(b);

	for (r = 0; r < n - 1; r++) {
		if (b->s_ij(r, r)->nullp()) {
			for (count = r + 1; count < n; count++) {
				if (!b->s_ij(count, r)->nullp()) {
					for (y = r; y < n; y++) {
						b->s_ij(count, y)->swap(b->s_ij(r, y));
						}
					if (debug) {
						fprintf(stderr, "MA::determinante: "
						"nach swap b = ");
						b->fprintln(stderr);
						}
					sign = -sign;
					break;
					}
				} /* for count */
			if (count == n) {
				s_ij(0, 0)->copy(erg);
				erg->zero();
				goto l_exit;
				}
			}
		if (debug) {
			b->s_ij(r, r)->fprint(stderr);
			fprintf(stderr," ungleich null\n");
			}
		for (i = r + 1; i < n; i++) {
			if (debug) {
				fprintf(stderr, "MA::determinante: "
				"b vor div = ");
				b->fprintln(stderr); 
				fprintf(stderr, "r = %ld\n", r);
				}
			b->s_ij(i, r)->sym_div(b->s_ij(r, r), factor);
			factor->addinvers_apply();
			if (debug) {
				fprintf(stderr, "MA::determinante: factor= ");
				factor->fprintln(stderr);
				}
			for (j = r; j < n; j++) {
				factor->mult(b->s_ij(r, j), temp);
				temp->add_apply(b->s_ij(i, j));
				}
			}
		if (debug) {
			fprintf(stderr, "MA::determinante: b = ");
			b->fprintln(stderr);
			}
		} /* for r */

	for (i = 0; i < n; i++) {
		if (i == 0) {
			b->s_ij(0, 0)->copy(product);
			}
		else {
			b->s_ij(i, i)->mult_apply(product);
			}
		}

	if (debug) {
		fprintf(stderr, "MA::determinante: b = ");
		b->fprintln(stderr);
		}

	if (sign == -1)
		product->addinvers(erg);
	else
		product->copy(erg);
	
l_exit:
	freeall(product);
	freeall(temp);
	freeall(factor);
	freeall(b);
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::inc_row()
#else
Increase the number of rows. The contents of the original matrix is preserved.
#endif
{
	return realloc_row(s_hi() + 1);
}

#if TEXDOCU
INT MATRIX_OB::realloc_row(INT new_len)
#else
Changes the number of rows of the matrix in this to new\_len.
new\_len can be smaller, equal to or larger than the original height.
#endif
{
	INT old_len, i, j;
	SYM_OP mem, old_mem;
	
	old_len = s_hi();
	if (new_len <= old_len) {
		for (i = new_len; i < old_len; i++) {
			for (j = 0; j < s_li(); j++)
				s_ij(i, j)->freeself();
			}
		s_h()->m_i(new_len);
		return OK;
		}
	old_mem = s_s();
	mem = (SYM_OP) my_malloc(new_len * s_li() * sizeof(struct object), "MA::realloc_row()");
	if (mem == NIL)
		return error("MATRIX::realloc_row(): no memory");
	my_memcpy(mem, old_mem, old_len * s_li() * sizeof(struct object));
	for (i = old_len * s_li(); i < new_len * s_li(); i++)
		(mem + i)->c_obj_k(EMPTY);
	c_s(mem);
	s_h()->m_i(new_len);
	if (old_mem)
		my_free(old_mem);
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::inc_col()
#else
Increase the number of columns. The contents of the original matrix is preserved.
#endif
{
	return realloc_column(s_li() + 1);
}

#if TEXDOCU
INT MATRIX_OB::realloc_column(INT new_len)
#else
Changes the number of columns of the matrix in this to new\_len.
new\_len can be smaller, equal to or larger than the original length.
#endif
{
	MATRIX_OB M;
	INT i, j, j1;
	
	M.m_ilih(new_len, s_hi());
	j1 = MINIMUM(new_len, s_li());
	for (i = 0; i < s_hi(); i++) {
		for (j = 0; j < j1; j++) {
			s_ij(i, j)->swap(M.s_ij(i, j));
			}
		}
	swap(&M);
	return OK;
}

INT MATRIX_OB::insert_row_at(INTEGER_OP len, INT i, VECTOR_OP V, INT V_len)
{
	INT N, new_n, l, j;
	
	if (i > len->s_i())
		return error("MATRIX::insert_row_at(): i > len");
	N = s_hi();
	if (len->s_i() >= N) {
		new_n = len->s_i() + MATRIX_OVERSIZE;
		if (realloc_row(new_n) != OK)
			return error("MATRIX::insert_row_at(): error in realloc_row()");
		}
	/* in l ist jetzt eine EMPTY Zeile, 
	 * diese wird vorgeswaped: */
	for (l = len->s_i(); l > i; l--) {
		for (j = 0; j < s_li(); j++)
			s_ij(l, j)->swap(s_ij(l - 1, j));
		}
	/* jetzt ist in i eine EMPTY Zeile */
	for (j = 0; j < V_len; j++)
		V->s_i(j)->copy(s_ij(i, j));
	len->inc();
	return OK;
}

static INT mtx_cmp_row_vec(MATRIX_OP M, VECTOR_OP V, INT row, INT nb_keys);

INT MATRIX_OB::search_row(
	INT len, INT nb_keys, 
	INT f_ascending, VECTOR_OP V, 
	INT *idx, INT *f_found)
/* V ist Vector mit mindestens nb_keys Eintraegen;
 * die matrix muss mindestens len Zeilen 
 * und nb_keys Spalten haben. */
{
	INT l, r, m, res;
	
	if (V == NIL || idx == NIL || 
		f_found == NIL)
		return error("MATRIX::search_row(): args nil");
	if (len == 0) {
		*idx = 0;
		*f_found = FALSE;
		return(OK);
		}
	l = 0;
	r = len;
	*f_found = FALSE;
	/* Invariante:
	 * M[i] <= V fuer i < l;
	 * M[i] >  V fuer i >= r;
	 * r - l ist Laenge des Suchbereichs. 
	 */
	while (l < r) {
		m = (l + r) >> 1;
		/* Bei einem Suchbereich gerader Laenge 
		 * wird das Element 
		 * oberhalb der Mitte genommen ! */
		res = mtx_cmp_row_vec(this, V, m, nb_keys);
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
	return(OK);
}

static INT mtx_cmp_row_vec(MATRIX_OP M, VECTOR_OP V, INT row, INT nb_keys)
{
	INT j, res;

	for (j = 0; j < nb_keys; j++) {
		res = M->s_ij(row, j)->sym_comp(V->s_i(j));
		if (res)
			return res;
		}
	return 0;
}

INT mtx_test()
{
	VECTOR_OB the_key;
	MATRIX_OB M;
	INTEGER_OB M_len;
	INT i, j, idx, f_found;
	
	M.m_ilih(2, MATRIX_OVERSIZE);
	the_key.m_il(2);
	M_len.m_i(0);
	j = 7; /* hoffentlich primitiv mod 61 (?!) nein ! */
	for (i = 0; i < 20; i++) {
		j = (j * j) % 61;
		((INTEGER_OP)the_key.s_i(0))->m_i(j);
		((INTEGER_OP)the_key.s_i(1))->m_i(i);
		M.search_row(M_len.s_i(), 1 /* nb_keys */, 
			TRUE, &the_key, &idx, &f_found);
		M.insert_row_at(&M_len, idx, &the_key, 2);
		}
	M.println();
	return OK;
}

INT MATRIX_OB::split_column_wise(VECTOR_OP S, INT col_width)
{
	MATRIX_OB M;
	INT m, n, i, j, jj, J, j0, JJ;

	m = s_hi();
	n = s_li();
	JJ = n / col_width;
	if (col_width * JJ < n)
		JJ++;
	S->m_il(JJ);
	for (J = 0; J < JJ; J++) {
		if (J < JJ - 1)
			jj = col_width;
		else
			jj = n - J * col_width;
		M.m_ilih(jj, m);
		j0 = J * col_width;
		for (i = 0; i < m; i++) {
			for (j = 0; j < jj; j++) {
				s_ij(i, j0 + j)->copy(M.s_ij(i, j));
				}
			}
		M.swap(S->s_i(J));
		}
	return OK;
}

INT MATRIX_OB::split_row_wise(VECTOR_OP S, INT nrow)
{
	MATRIX_OB M;
	INT m, n, i, j, ii, I, i0, II;

	m = s_hi();
	n = s_li();
	II = m / nrow;
	if (nrow * II < m)
		II++;
	S->m_il(II);
	for (I = 0; I < II; I++) {
		if (I < II - 1)
			ii = nrow;
		else
			ii = m - I * nrow;
		M.m_ilih(n, ii);
		i0 = I * nrow;
		for (i = 0; i < ii; i++) {
			for (j = 0; j < n; j++) {
				s_ij(i0 + i, j)->copy(M.s_ij(i, j));
				}
			}
		M.swap(S->s_i(I));
		}
	return OK;
}

INT MATRIX_OB::cyclic_block_matrix(VECTOR_OP blocks, INT f_row_wise, INT f_forward)
{
	MATRIX_OP p;
	INT n = 0, N, Nn, n1;
	INT I, J, J1, i0, j0, i, j;
	
	N = blocks->s_li();
	for (i = 0; i < N; i++) {
		p = (MATRIX_OP) blocks->s_i(i);
		if (p->s_obj_k() == MATRIX) {
			n1 = p->s_li();
			if (p->s_hi() != n1)
				return error("MA::cyclic_block_matrix() not symmetrical");
			}
		else {
			n1 = 1;
			}
		if (i == 0) {
			n = n1;
			}
		else {
			if (n1 != n)
				return error("MA::cyclic_block_matrix() matrices of different size");
			}
		}
	Nn = n * N;
	m_ilih(Nn, Nn);
	if (f_row_wise) {
		for (I = 0; I < N; I++) {
			i0 = I * n;
			for (J = 0; J < N; J++) {
				if (f_forward) {
					J1 = J - I;
					if (J1 < 0)
						J1 += N;
					}
				else {
					J1 = J + I;
					if (J1 >= N)
						J1 -= N;
					}
				if (J1 < 0 || J1 >= N)
					return error("MA::cyclic_block_matrix() J1 out of range");
				j0 = J * n;
#if 0
				printf("MA::cyclic_block_matrix() I=%ld J=%ld J1=%ld i0=%ld j0=%ld\n", 
					I, J, J1, i0, j0);
#endif
				p = (MATRIX_OP) blocks->s_i(J1);
				if (p->s_obj_k() != MATRIX) {
					((SYM_OP)p)->copy(s_ij(i0, j0));
					}
				else {
					for (i = 0; i < n; i++) {
						for (j = 0; j < n; j++) {
							p->s_ij(i, j)->copy(s_ij(i0 + i, j0 + j));
							}
						}
					}
				}
			}
		}
	else {
		return error("MA::cyclic_block_matrix() !f_row_wise not yet implemented");
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::perm_rep(PERM_REP *P, PERMUTATION_OP p)
#else
Computes the permutation of a permutation matrix.
Assumes that this is a permutation matrix.
Needs a PERM\_REP P and calls the function perm\_rep\_mtx.
#endif
{
	INT *mtx;
	INT n, q, i, j, a;
	
	q = P->q;
	n = s_hi();
	if (n != s_li())
		return error("MATRIX_OB::perm_rep() not quadratic");
	if (n != P->n) 
		return error("MATRIX_OB::perm_rep() n != P->n");
	mtx = (INT *) my_malloc(n * n * sizeof(INT), "MATRIX_OB::perm_rep(): mtx");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a = s_iji(i, j);
			if (a >= q || a < 0)
				return error("MATRIX_OB::perm_rep() a >= q || a < 0");
			mtx[i * n + j] = a;
			}
		}
	perm_rep_mtx(P, mtx, p);
	my_free(mtx);
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::matrix_of_perm_rep(PERM_REP *P, PERMUTATION_OP p)
#else
Computes a permutation matrix for a given permutation p.
#endif
{
	INT *mtx;
	INT n, q, i, j, a;
	
	q = P->q;
	n = P->n;
	m_ilih_n(n, n);
	
	mtx = (INT *) my_malloc(n * n * sizeof(INT), "MATRIX_OB::matrix_of_perm_rep(): mtx");
	perm_rep_to_mtx(P, p, mtx);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a = mtx[i * n + j];
			if (a >= q || a < 0)
				return error("MATRIX_OB::matrix_of_perm_rep() a >= q || a < 0");
			m_iji(i, j, a);
			}
		}
	my_free(mtx);
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::complement()
#else
For a zero/one matrix this, the function computes the complementary matrix 
which has zeros and ones interchanged.
#endif
{
	INT h, l, i, j, a, b;
	
	h = s_hi();
	l = s_li();
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			a = s_iji(i, j);
			if (a != 0 && a != 1) 
				return error("MA::complement not a 0/1-matrix");
			b = 0;
			if (a == 0)
				b = 1;
			m_iji(i, j, b);
			}
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::stirling_first(INT n, INT f_extended, INT f_signless, INT f_v)
#endif
{
	INT i, j;
	
	if (f_extended) {
		m_ilih_n(n + 1, n + 1);
		for (i = 0; i <= n; i++) {
			for (j = 0; j <= n; j++) {
				::stirling_first(i, j, 
					f_signless, s_ij(i, j), f_v);
				}
			}
		}
	else {
		m_ilih_n(n, n);
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				::stirling_first(i + 1, j + 1, 
					f_signless, s_ij(i, j), f_v);
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::stirling_second(INT n, INT f_extended, INT f_ordered, INT f_v)
#endif
{
	INT i, j;
	
	if (f_extended) {
		m_ilih_n(n + 1, n + 1);
		for (i = 0; i <= n; i++) {
			for (j = 0; j <= n; j++) {
				::stirling_second(i, j, 
					f_ordered, s_ij(i, j), f_v);
				}
			}
		}
	else {
		m_ilih_n(n, n);
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				::stirling_second(i + 1, j + 1, 
					f_ordered, s_ij(i, j), f_v);
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::binomial(INT n, INT f_extended, INT f_inverse, INT f_v)
#endif
{
	INT i, j;
	
	if (f_extended) {
		m_ilih_n(n + 1, n + 1);
		for (i = 0; i <= n; i++) {
			for (j = 0; j <= n; j++) {
				Binomial(j, i, s_ij(i, j));
				if (f_inverse) {
					if (ODD(i + j))
						s_ij(i, j)->addinvers_apply();
					}
				}
			}
		}
	else {
		m_ilih_n(n, n);
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				Binomial(j + 1, i + 1, s_ij(i, j));
				if (f_inverse) {
					if (ODD(i + j))
						s_ij(i, j)->addinvers_apply();
					}
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::rank_p(INT p)
#endif
{
	INT i, j, m, n, a, r;
	G_INT_TYPE *M;
	
	m = s_hi();
	n = s_li();
	M = (G_INT_TYPE *) my_malloc(m * n * sizeof(G_INT_TYPE), "MA::rank_p");
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			a = s_iji(i, j);
			a = g_asr(a, p);
			M[i * n + j] = a;
			}
		}
	r = g_rank(M, m, n, p);
	my_free(M);
	return r;
}

#if TEXDOCU
INT MATRIX_OB::calc_hash_key(INT key_len, BYTE *hash_key, INT f_v)
#endif
{
	INT al_len;
	BYTE *alphabet = NIL;
	BYTE *inc = NIL;
	static INT primes[] = { 
		2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 
		31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 
		73, 79, 83, 89, 97 } ; // the first 25 primes 
	INT i0, i, j, k, v, b, nb_inc, pr, x, y;
	BYTE c;
	INT f_vv = FALSE;
	
	v = s_hi();
	b = s_li();
	nb_inc = 0;
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			if (s_iji(i, j) != 0)
				nb_inc++;
			}
		}
			
	al_len = MAXIMUM(256, b);
	alphabet = (BYTE *) my_malloc(sizeof(BYTE) * (al_len + 1), "alphabet for GEO_BY_BASE_BLOCKS");
	inc = (BYTE *) my_malloc(sizeof(BYTE) * (nb_inc + 1), "alphabet for GEO_BY_BASE_BLOCKS");
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
		printf("alphabet: %s\n", alphabet);
		}
	
	k = 0;
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			if (s_iji(i, j) == 0)
				continue;
			c = alphabet[j];
			inc[k] = c;
			k++;
			}
		}
	inc[nb_inc] = 0;
	if (f_vv) {
		printf("incidences: %s\n", inc);
		}
	
	j = 0;
	for (k = 0; k < key_len; k++) {
		pr = primes[k % 25];
		x = 0;
		for (i = 0; i < nb_inc; i++) {
			y = (INT) inc[j];
			x = (x + y) % 256;
			j += pr;
			if (j >= nb_inc)
				j = 0;
			}
		hash_key[k] = alphabet[x];
		// printf("k=%ld pr = %ld x=%ld h[k]=%c\n", k, pr, x, h[k]);
		}
	hash_key[key_len - 1] = 0;
	if (f_v) {
		printf("MATRIX_OB::calc_hash_key() (len=%ld) hash=%s\n", 
			key_len, hash_key);
		}
	my_free(alphabet);
	my_free(inc);
	return OK;
}


#endif /* MATRIXTRUE */

