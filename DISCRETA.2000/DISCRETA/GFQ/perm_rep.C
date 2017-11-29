/* perm_rep.C
 * 
 * Anton Betten 
 * Jul 30, 1995 
 */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */

#include <DISCRETA/discreta.h>

#ifdef GFQ_TRUE

#include <DISCRETA/perm.h>
#include <DISCRETA/gfq.h>
#include <DISCRETA/ladder.h>

#undef DEBUG_PERM_REP_MTX
#undef DEBUG_PERM_REP_TO_MATRIX
#undef DEBUG_PERM_REP_DETERMINANTE
#undef DEBUG_GL2SL_GO

#define MAX_DEGREE_FOR_PRINT 70

static void print_generators(VECTOR_OP V);

#if TEXDOCU
INT perm_rep_GL_n_q_degree(INT f_affine, 
	INT f_semi, INT f_special, INT q, INT n, 
	INT f_verbose)
#else
Gets the degree of the permutation representation 
by opening a temporary PERM\_REP structure.
#endif
{
	PERM_REP *P;
	INT deg;
	
	P = open_perm_rep(q, f_affine, n, f_verbose);
	deg = P->deg;
	perm_rep_free(P);
	return deg;
}

#if TEXDOCU
INT perm_rep_T_n_q_generators(INT q, INT n, VECTOR_OP V, INT f_verbose)
#else
Gets generators for the normal subgroup of translastions.
#endif
{
	PERM_REP *P;
	
	// f_verbose = TRUE;
	if (f_verbose) {
		printf("perm_rep_T_n_q_generators(): "
			"calling open_perm_rep()\n");
		fflush(stdout);
		}
	P = open_perm_rep(q, TRUE /* f_affine */, n, f_verbose);
	perm_rep_translation_generators(P, V, f_verbose);
	perm_rep_free(P);
	return OK;
}

#if TEXDOCU
INT perm_rep_GL_n_q_generators(INT f_affine, 
	INT f_semi, INT f_special, INT f_translations, INT q, INT n, 
	VECTOR_OP V, INT f_verbose)
#else
Gets generators for $\PSL$, $\PGL$, $\PSSL$, $\PGGL$, 
$\ASL$, $\AGL$, $\ASSL$, $\AGGL$
of dimension $n$ over $GF(q)$.
The generators are written into the vector V.
#endif
{
	PERM_REP *P;
	
	// f_verbose = TRUE;
	if (f_verbose) {
		printf("perm_rep_GL_n_q_generators(): "
			"calling open_perm_rep()\n");
		fflush(stdout);
		}
	P = open_perm_rep(q, f_affine, n, f_verbose);
	if (f_verbose) {
		printf("perm_rep_GL_n_q_generators(): "
			"calling perm_rep_GL_generators()\n");
		fflush(stdout);
		}
	perm_rep_GL_generators(P, f_semi, f_special, f_translations, V, f_verbose);
	if (f_verbose) {
		printf("perm_rep_GL_n_q_generators(): "
			"calling perm_rep_free()\n");
		fflush(stdout);
		}
	perm_rep_free(P);
	if (f_verbose) {
		printf("perm_rep_GL_n_q_generators(): "
			"finished!\n");
		fflush(stdout);
		}
	return OK;
}

#if TEXDOCU
PERM_REP *open_perm_rep(INT q, 
	INT f_affine, INT n, INT f_verbose)
#else
Opens a PERM\_REP structure. 
If f\_affine ist TRUE, the representation is affine, 
otherwise it is projective.
If $q$ is a prime power $q=p^f$ with $f > 1$, 
a ZECH\_DATA structure is opened via zech\_open.
#endif
{
	PERM_REP *P1;
	VECTOR_OB vp, ve;
	
	P1 = (PERM_REP *) my_malloc(sizeof(PERM_REP), "open_perm_rep");
	factor_integer(q, &vp, &ve);
	if (vp.s_li() > 1) {
		error("q not prime(power)");
		return NULL;
		}
	P1->f_affine = f_affine;
	P1->q = q;
	P1->f = ve.s_ii(0);
	P1->p = vp.s_ii(0);
	P1->n = n;
	P1->f_verbose = f_verbose;
	if (P1->f > 1) {
		P1->Z = zech_open(P1->f, P1->p, f_verbose);
		printf("open_perm_rep(): Zech table opened;\n");
		fflush(stdout);
		P1->idx_zero = 0;
		P1->idx_one = P1->Z->Num_inv[0];
		}
	else {
		P1->Z = NULL;
		P1->idx_zero = 0;
		P1->idx_one = 1;
		}
	
	if (f_affine)
		P1->deg = z_affine_degree(P1->q, P1->n);
	else
		P1->deg = z_projective_degree(P1->q, P1->n);
	P1->vec1 = (INT *) my_malloc(n * sizeof(INT), "open_perm_rep");
	P1->vec2 = (INT *) my_malloc(n * sizeof(INT), "open_perm_rep");
	if (f_verbose) {
		printf("open_perm_rep(): PERM_REP opened;\n"
			"q = %ld f_affine = %ld n = %ld degree = %ld\n", 
			P1->q, P1->f_affine, P1->n, P1->deg);
		fflush(stdout);
		if (P1->deg < MAX_DEGREE_FOR_PRINT) {
			perm_rep_print_vectors(P1);
			fflush(stdout);
			}
		}
	return P1;
}

#if TEXDOCU
void perm_rep_free(PERM_REP *P)
#else
Closes the PERM\_REP structure.
#endif
{
	if (P->vec1)
		my_free(P->vec1);
	if (P->vec2)
		my_free(P->vec2);
	if (P->Z) {
		zech_free(P->Z);
		}
	my_free(P);
}

#if TEXDOCU
INT perm_rep_i2vec(PERM_REP *P, INT i, INT *vec)
#else
converts an integer to the corresponding vector in the 
affine / projective representation determined by PERM\_REP.
#endif
{
	if (P->f_affine)
		return z_affine_i2vec(P->q, P->n, i, vec);
	else {
		return z_projective_i2vec(P->q, P->idx_zero, P->idx_one, 
			P->n, i, vec);
		}
}

#if TEXDOCU
INT perm_rep_vec2i(PERM_REP *P, INT *vec, INT *i)
#else
converts a vector in the 
affine / projective representation 
to the corresponding integer.
This function requires the vector to be normalized.
#endif
{
	if (P->f_affine)
		return z_affine_vec2i(P->q, P->n, vec, i);
	else
		return z_projective_vec2i(P->q, P->idx_zero, P->idx_one, 
			P->n, vec, i);
}

#if TEXDOCU
INT perm_rep_vec2i_with_normalization(PERM_REP *P, INT *vec, INT *i)
#else
Same as perm\_rep\_vec2i but if the representation is projective, 
the vector is normalized first.
#endif
{
	if (P->f_affine)
		return z_affine_vec2i(P->q, P->n, vec, i);
	else {
		if (P->f > 1) {
			z_vec_norm_projective_GFq_num(P->Z, vec, P->n);
			}
		else {
			z_vec_norm_projective_GFp(P->p, vec, P->n);
			}
		return z_projective_vec2i(P->q, P->idx_zero, P->idx_one, 
			P->n, vec, i);
		}
}

#if TEXDOCU
INT z_test_affine_rep(INT q, INT n)
#else
Tests the affine representation by converting 
vectors and integers back and forth.
#endif
{
	INT deg, i, ii;
	INT *vec;
	
	vec = (INT *) my_malloc(n * sizeof(INT), "z_test_affine_rep");
	deg = z_affine_degree(q, n);
	printf("z_test_affine_rep(): degree = %ld\n", deg);
	for (i = 0; i < deg; i++) {
		printf("%ld: ", i);
		fflush(stdout);
		z_affine_i2vec(q, n, i, vec);
		z_vec_i_print(vec, n);
		fflush(stdout);
		z_affine_vec2i(q, n, vec, &ii);
		printf(" %ld\n", ii);
		fflush(stdout);
		if (ii != i)
			return error("z_test_affine_rep() i != ii");
		}
	my_free(vec);
	return OK;
}

#if TEXDOCU
INT z_affine_degree(INT q, INT n)
#else
returns $q^n$.
#endif
{
	return i_power_j(q, n);
}

#if TEXDOCU
INT z_affine_i2vec(INT q, INT n, INT i, INT *vec)
#else
Converts an integer into the corresponding vector in 
the affine representation. Here, we use 
{\em positive} representatives mod p.
#endif
{
	INT j, r;
	
	j = 0;
	while (i != 0) {
		r = i % q;
		vec[j] = r;
		j++;
		i -= r;
		i /= q;
		if (j > n)
			return error("z_aff_i2vec() i too large !");
		}
	for ( ; j < n; j++)
		vec[j] = 0;
	return OK;
}

#if TEXDOCU
INT z_affine_vec2i(INT q, INT n, INT *vec, INT *i)
#else
Converst an vector in the affine representation to the 
corresponding integer.
The entries of the vector must be between $0$ and $q - 1$.
#endif
{
	INT a, b, l1;
	
	a = 0;
	for (l1 = n - 1; l1 >= 0; l1--) {
		b = vec[l1];
		if (b >= q || b < 0)
			return error("z_affine_vec2i(): b out of range!\n");
		a += b;
		if (l1 > 0)
			a *= q;
		}
	*i = a;
	return OK;
}

#if TEXDOCU
INT z_test_projective_rep(INT f, INT p, INT n)
#else
Tests the projective representation by converting back and forth between 
integers and vectors.
#endif
{
	ZECH_DATA *Z;
	INT idx_zero; /* the zero element: 0 * X^deg - 1 + ... 0 X + 0 */
	INT idx_one; /* the element 1 = X^0 */
	INT i, ii, deg, q;
	INT *vec;
	
	if (f > 1) {
		Z = zech_open(f, p, FALSE);
		q = Z->q;
		idx_zero = 0;
		idx_one = Z->Num_inv[0];
		}
	else {
		q = p;
		idx_zero = 0;
		idx_one = 1;
		}
	vec = (INT *) my_malloc(n * sizeof(INT), "z_test_projective_rep");
	if (vec == NULL)
		return error("z_test_projective_rep() no memory");
	
	deg = z_projective_degree(q, n);
	printf("z_test_projective_rep(): degree = %ld\n", deg);
	for (i = 0; i < deg; i++) {
		printf("%ld: ", i);
		fflush(stdout);
		z_projective_i2vec(q, idx_zero, idx_one, n, i, vec);
		z_vec_i_print(vec, n);
		fflush(stdout);
		z_projective_vec2i(q, idx_zero, idx_one, n, vec, &ii);
		printf(" %ld\n", ii);
		fflush(stdout);
		if (ii != i)
			return error("z_test_projective_rep() i != ii");
		}
	
	if (f > 1) {
		zech_free(Z);
		}
	my_free(vec);
	return OK;
}
#if TEXDOCU
INT z_projective_degree(INT q, INT n)
#else
Returns
\[
\frac{q^n - 1}{q-1} = \sum_{i=0}^{n-1} q^i 
\]
#endif
{
	INT qhl, l, deg;
	
	l = 0;
	qhl = 1;
	deg = 0;
	while (l < n) {
		deg += qhl;
		qhl *= q;
		l++;
		}	
	return deg;
}

#if TEXDOCU
INT z_projective_i2vec(INT q, INT idx_zero, INT idx_one, 
	INT n, INT i, INT *vec)
#else
Converts an integer into the corresponding vector 
in the projective representation.
#endif
{
	INT l, qhl, k, j, r;
	
	l = 0;
	qhl = 1;
	while (l < n) {
		if (i >= qhl) {
			i -= qhl;
			qhl *= q;
			l++;
			continue;
			}
		/* printf("l=%ld idx_zero = %ld idx_one = %ld\n", l, idx_zero, idx_one); */
		vec[l] = idx_one;
		for (k = l + 1; k < n; k++) {
			vec[k] = idx_zero;
			}
		j = 0;
		while (i != 0) {
			r = i % q;
			vec[j] = r;
			j++;
			i -= r;
			i /= q;
			}
		for ( ; j < l; j++)
			vec[j] = 0;
		return OK;
		}
	return error("z_projective_i2vec(): i too large");
}

#if TEXDOCU
INT z_projective_vec2i(INT q, INT idx_zero, INT idx_one, 
	INT n, INT *vec, INT *i)
#else
Converts a vector in the projective representation into 
the corresponding integer.
#endif
{
	INT l, l1, qhl, ii = 0, a, b;

	for (l = 0; l < n; l++) {
		b = vec[l];
		if (b >= q || b < 0)
			return error("z_projective_vec2i(): b out of range!\n");
		}
	l = n - 1;
	while (vec[l] == idx_zero && l >= 0)
		l--;
	if (l >= 0) {
		if (vec[l] != idx_one) {
			return error("z_projective_vec2i() vec not normalized !");
			}
		qhl = 1;
		for (l1 = 0; l1 < l; l1++) {
			ii += qhl;
			qhl *= q;
			}
		a = 0;
		for (l1 = l - 1; l1 >= 0; l1--) {
			a += vec[l1];
			if (l1 > 0)
				a *= q;
			}
		ii += a;
		*i = ii;
		return OK;
		}
	return error("z_projective_vec2i() zero vector !");
}

#if TEXDOCU
INT perm_rep_index_of_ei(PERM_REP *P, INT i)
#else
Determines the integer in the permutation representation 
PERM\_REP corresponding to the vector $e_i$ 
which has a 1 at its $i$-th position and zeros elsewhere.
#endif
{
	INT ii, j, qhii;
	
	if (i >= P->n)
		return error("perm_rep_index_of_ei(): i >= n");
	if (P->f_affine)
		j = i_power_j(P->q, i);
	else {
		j = 0;
		qhii = 1;
		for (ii = 0; ii < i; ii++) {
			j += qhii;
			qhii *= P->q;
			}
		}
	return j;
}

#if TEXDOCU
INT perm_rep_mtx(PERM_REP *P, INT *mtx, PERMUTATION_OP perm)
#else
Converts a matrix into a permutation describing the transformation 
as a permutation of the elements in the projective / affine 
representation.
#endif
{
	INT i, ii;
	
#ifdef DEBUG_REP_ZP
	printf("perm_rep_mtx: f_affine = %ld deg = %ld\n", 
		P->f_affine, P->deg);
	z_mtx_i_print(mtx, P->n);
#endif
	perm->m_il(P->deg);
	for (i = 0; i < P->deg; i++) {
		perm_rep_i2vec(P, i, P->vec1);

#ifdef DEBUG_PERM_REP_MTX
		printf("%ld: ", i);
		z_vec_i_print(P->vec1, P->n);
		fflush(stdout);
#endif
		if (P->f == 1)
			z_mtx_vec_mult_GFp(P->p, mtx, P->n, P->vec1, P->vec2);
		else
			z_mtx_vec_mult_GFq_num(P->Z, mtx, P->n, P->vec1, P->vec2);
		if (!P->f_affine) {
			if (P->f == 1)
				z_vec_norm_projective_GFp(P->p, P->vec2, P->n);
			else
				z_vec_norm_projective_GFq_num(P->Z, P->vec2, P->n);
			}
#ifdef DEBUG_PERM_REP_MTX
		z_vec_i_print(P->vec2, P->n);
		fflush(stdout);
#endif
		perm_rep_vec2i(P, P->vec2, &ii);
#ifdef DEBUG_PERM_REP_MTX
		printf(" %ld\n", ii);
		fflush(stdout);
#endif
		perm->m_ii(i, ii + 1);
		}
#ifdef DEBUG_PERM_REP_MTX
	perm->println();
#endif
	return OK;
}

#if TEXDOCU
INT perm_rep_to_mtx(PERM_REP *P, PERMUTATION_OP perm, INT *mtx)
#else
Computes the matrix determined by the transformation 
given by the permutation perm.
uses P->vec1 
#endif
{
	INT n = P->n;
	INT i, j, k, a;
	
	for (i = 0; i < n; i++) {
		/* compute the i-th column of mtx[]: */
		
#if 0
		if (P->f_affine) {
			j = i_power_j(q, i);
			}
		else {
			j = 0;
			qhii = 1;
			for (ii = 0; ii < i; ii++) {
				j += qhii;
				qhii *= q;
				}
			}
#endif
		j = perm_rep_index_of_ei(P, i);
		/* the vector e_i = (0 .. 0 1 0 .. 0) 
		 * (1 in i-th position) 
		 * is numbered by j in the permutation representation */
		k = perm->s_ii(j) - 1;
		perm_rep_i2vec(P, k, P->vec1);
		for (j = 0; j < n; j++) {
			a = P->vec1[j];
			mtx[j * n + i] = a;
			}
		}
#ifdef DEBUG_PERM_REP_TO_MATRIX
	printf("perm_rep_to_mtx:\n");
	perm->println();
	printf("is the matrix:\n");
	z_mtx_i_print(mtx, P->n);
#endif
	return OK;
}

#if TEXDOCU
INT perm_rep_determinante(PERM_REP *P, PERMUTATION_OP perm)
#else
derterminant of the matrix corresponding to the permutation perm 
in the representation.
Calls perm\_rep\_determinante\_v().
#endif
{
	return perm_rep_determinante_v(P, perm, FALSE);
}

#if TEXDOCU
INT perm_rep_determinante_v(PERM_REP *P, PERMUTATION_OP perm, INT f_verbose)
#else
derterminant of the matrix corresponding to the permutation perm 
in the representation.
#endif
{
	INT *mtx;
	INT det;
	INT n = P->n;
	
	mtx = (INT *) my_malloc(n * n * sizeof(INT), "perm_rep_determinante_v");
	perm_rep_to_mtx(P, perm, mtx);
	if (P->f > 1)
		det = z_mtx_determinante_GFq_num(P->Z, mtx, P->n);
	else
		det = z_mtx_determinante_GFp(P->p, mtx, P->n);
	if (f_verbose) {
		if (perm->s_li() < MAX_DEGREE_FOR_PRINT) {
			printf("the permutation");
			perm->println();
			printf("in matrix form:\n");
			z_mtx_i_print(mtx, P->n);
			printf("has determinant: %ld\n", det);
			}
		}
	my_free(mtx);
	return det;
}

#if TEXDOCU
INT perm_rep_translation_generators(PERM_REP *P, VECTOR_OP V, INT f_verbose)
#else
Computes generators for the translations. 
Generators are for each dimension and for each dimension 
of the additive group of the field over the prime field.
#endif
{
	INT V_len = 0;
	INT i, j;
	PERMUTATION_OB perm;
	
	V->m_il(0);
	for (i = 0; i < P->n; i++) {
		for (j = 0; j < P->f; j++) {
			if (P->f == 1)
				z_translation_rep_Zp(P->p, P->n, i, &perm);
			else
				z_translation_rep_Zq_num(P->Z, P->n, i, j, &perm);
			if (f_verbose) {
				printf("translation ei = %ld betaj = %ld\n", i, j);
				perm.println();
				fflush(stdout);
				}
			V->inc();
			perm.swap((PERMUTATION_OP) V->s_i(V_len));
			V_len++;
			}
		}
	return OK;
}

#if TEXDOCU
INT perm_rep_GL_generators(PERM_REP *P, 
	INT f_semi, INT f_special, INT f_translations, 
	VECTOR_OP V, INT f_verbose)
#else
affine or projective is already determined by PERM\_REP.
#endif
{
	PERM_REP *P1;
	INT i, j;
	INT **mtx, nb_gen = 0, det;
	INT *Mtx;
	PERMUTATION_OB perm;
	PERMUTATION_OP p_perm;
	INT V_len;
	
	// f_verbose = TRUE;
	if (f_verbose) {
		printf("perm_rep_GL_generators(\n"); fflush(stdout);
		printf("perm_rep_GL_generators(): \n"
			"f_affine = %ld f_semi = %ld f_special = %ld f_translations = %ld "
			"f = %ld p = %ld n = %ld\n", 
			P->f_affine, f_semi, f_special, f_translations, P->f, P->p, P->n);
		printf("as a group of permutations "
			"of degree %ld\n", P->deg);
		}
	if (f_special) {
		if (f_verbose) {
			printf("opening affine representation "
				"(only internally needed for GL2SL):\n");
			}
		P1 = open_perm_rep(P->q, 
			TRUE /* f_affine */ , P->n, f_verbose);
		}
	Mtx = (INT *) my_malloc(P->n * P->n * sizeof(INT), "perm_rep_GL_generators");
		
	/* get generators for GL: */
	if (f_verbose) {
		printf("generators as matrices:");
		fflush(stdout);
		}
	if (P->f > 1) {
		z_GL_n_q_gen_num(P->Z, &mtx, &nb_gen, P->n);
		}
	else {
		if (f_semi) {
			printf("perm_rep_GL_generators(): "
				"ignoring f_semi in case f = 1\n");
			f_semi = FALSE;
			}
		z_GL_n_p_gen(&mtx, &nb_gen, P->n, P->p);
		}
	if (f_verbose) {
		printf(" OK\n");
		fflush(stdout);
		}
	
	/* convert the generators into permutations
	 * (affine representation if f_special): */
	if (f_verbose) {
		printf("conversion into permutations (nb_gen=%ld):\n", nb_gen);
		fflush(stdout);
		}
	V->m_il(nb_gen);
	V_len = nb_gen;
	for (i = 0; i < nb_gen; i++) {
		if (f_verbose) {
			printf("%ld:\n", i);
			z_mtx_i_print(mtx[i], P->n);
			}
		if (f_special)
			perm_rep_mtx(P1, mtx[i], &perm);
		else
			perm_rep_mtx(P, mtx[i], &perm);
		if (f_verbose) {
			if (!f_special)
				perm.println();
			else
				printf("\n");
			fflush(stdout);
			}
		perm.copy((PERMUTATION_OP) V->s_i(i));
		
		if (f_special) {
			det = perm_rep_determinante(P1, (PERMUTATION_OP) V->s_i(i));
			if (f_verbose) {
				printf("det = %ld\n", det);
				}
			}
		}
	if (f_verbose) {
		printf(" OK\n");
		fflush(stdout);
		}
	
	if (f_special) {
		VECTOR_OB V1;
		
		if (f_verbose) {
			printf("calling GL2SL:");
			fflush(stdout);
			}
		perm_rep_GL2SL(P1, V /* GL_gen */, 
			&V1 /* SL_gen */, TRUE /* f_verbose */);
		V1.swap(V);
		V_len = V->s_li();
		printf("found %ld generators for SL\n", V_len);
		for (i = 0; i < V_len; i++) {
			p_perm = (PERMUTATION_OP) V->s_i(i);
			det = perm_rep_determinante_v(P1, p_perm, f_verbose);
			if (f_verbose) {
				printf("%ld: det = %ld\n", i, det);
				}
			}
		/* convert into projective representation again: */
		/* (or eventually affine) */
		if (f_verbose) {
			printf("GENERATORS:\n");
			}
		V1.m_il(V_len);
		for (i = 0; i < V_len; i++) {
			perm_rep_to_mtx(P1, (PERMUTATION_OP) V->s_i(i), Mtx);
			perm_rep_mtx(P, Mtx, (PERMUTATION_OP) V1.s_i(i));
			if (f_verbose) {
				V1.s_i(i)->println();
				}
			}
		V1.swap(V);
		if (f_verbose) {
			printf(" OK\n");
			fflush(stdout);
			}
		}
	
	if (f_semi) {
		z_Frobenius_rep_Zq_num(P->f_affine, P->Z, P->n, &perm);
		if (f_verbose) {
			printf("Frobenius map:\n");
			perm.println();
			fflush(stdout);
			}
		V->inc();
		perm.copy((PERMUTATION_OP) V->s_i(V_len));
		V_len++;
		}
	
	if (f_verbose) {
		print_generators(V);
		}
	
	/* generators for translations in the affine case: */
	if (P->f_affine && f_translations) {
		if (f_verbose) {
			printf("generators for translations:\n");
			fflush(stdout);
			}
		for (i = 0; i < P->n; i++) {
			for (j = 0; j < P->f; j++) {
				if (P->f == 1)
					z_translation_rep_Zp(P->p, P->n, i, &perm);
				else
					z_translation_rep_Zq_num(P->Z, P->n, i, j, &perm);
				if (f_verbose) {
					printf("translation ei = %ld betaj = %ld\n", i, j);
					perm.println();
					fflush(stdout);
					}
				V->inc();
				perm.swap((PERMUTATION_OP) V->s_i(V_len));
				V_len++;
				}
			}
		}
	if (f_verbose) {
		printf("finished with perm_rep_GL_generators(), freeing data\n"); fflush(stdout);
		}
	if (f_special) {
		perm_rep_free(P1);
		}
	my_free(Mtx);
	for (i = 0; i < nb_gen; i++) {
		my_free(mtx[i]);
		}
	my_free(mtx);
	if (f_verbose) {
		printf("exiting from perm_rep_GL_generators()\n"); fflush(stdout);
		}
	return OK;
}

static void print_generators(VECTOR_OP V)
{
	FILE *fp;
	INT V_len, i, deg;
	PERMUTATION_OP p_perm;
			
	V_len = V->s_li();
	printf("generators:\n");
	V->Print();
#if 0
	printf("GENERATORS, cycle notation (%ld):\n", V_len);
	for (i = 0; i < V_len; i++) {
		p_perm = (PERMUTATION_OP) V->s_i(i);
		p_perm->print_list();
		}
	printf("GENERATORS, list notation:\n");
	for (i = 0; i < V_len; i++) {
		p_perm = (PERMUTATION_OP) V->s_i(i);
		p_perm->print_list();
		deg = p_perm->s_li();
		}
#endif
	fp = fopen("generators", "w");
	fprintf(fp, "%ld %ld\n", V_len, deg);
	for (i = 0; i < V_len; i++) {
		p_perm = (PERMUTATION_OP) V->s_i(i);
		p_perm->fprint_list(fp);
		}
	fclose(fp);
	printf("written to file generators\n");
	fflush(stdout);
}

/* compare with
 * GL2SL
 * in ma.C !
 */

static INT trace(VECTOR_OP SVlast, VECTOR_OP SVgen, 
	VECTOR_OP G, INT i, SYM_OP p);

#if TEXDOCU
INT perm_rep_GL2SL(PERM_REP *P, 
	VECTOR_OP GL_gen, VECTOR_OP SL_gen, INT f_verbose)
#else
Computes generators for SL given generators for GL.
#endif
{
	PERMUTATION_OB M1, M2, M3, M4;
	PERMUTATION_OP G, VM1;
	INTEGER_OB det_M;
	VECTOR_OB V, det_V, SVlast, SVgen;
	INT nb_gen, i, j, k, next, det;
	
	nb_gen = GL_gen->s_li();
	if (nb_gen == 0)
		return error("perm_rep_GL2SL() no generators");
	if (GL_gen->s_i(0)->s_obj_k() != PERMUTATION)
		return error("perm_rep_GL2SL() not a permutation");
#ifdef DEBUG_GL2SL_GO
	printf("GL2SL: GL group order:\n");
	go(GL_gen);
#endif
	SL_gen->m_il(0);
	V.m_il(1);
	det_V.m_il(1);
	GL_gen->s_i(0)->copy(V.s_i(0));
	V.s_i(0)->one();
	det = perm_rep_determinante(P, (PERMUTATION_OP) V.s_i(0));
	det_M.m_i(det);
	det_M.copy((INTEGER_OP) det_V.s_i(0));
#if 0
	/* must use idx_zero if GFq ! */
	if (!det_V.s_i(0)->einsp())
		return error("perm_rep_GL2SL() det (I) != 1");
#endif
	SVlast.m_il(1);
	SVgen.m_il(1);
	SVlast.m_ii(0, -1);
	SVgen.m_ii(0, -1);
	i = 0;
#ifdef DEBUG_GL2SL_GO
	if (f_verbose) {
		printf("in perm_rep_GL2SL():\n");
		}
#endif
	while (i < V.s_li()) {
		for (j = 0; j < nb_gen; j++) {
			VM1 = (PERMUTATION_OP) V.s_i(i);
			G = (PERMUTATION_OP) GL_gen->s_i(j);
			VM1->mult(G, &M1);
			det = perm_rep_determinante(P, (PERMUTATION_OP) &M1);
			det_M.m_i(det);
			next = det_V.search_linear(&det_M);
			if (next < 0) {
				V.inc();
				det_V.inc();
				SVlast.inc();
				SVgen.inc();
				next = V.s_li() - 1;
#ifdef DEBUG_GL2SL_GO
				if (f_verbose)
					printf("%ld: det = %ld\n", V.s_li(), det);
#endif
				M1.swap(V.s_i(next));
				det_M.swap(det_V.s_i(next));
				SVlast.m_ii(next, i);
				SVgen.m_ii(next, j);
				}
			else {
				trace(&SVlast, &SVgen, GL_gen, next, &M2);
				M2.invers(&M3);
				M1.mult(&M3, &M4);
				k = SL_gen->search_linear(&M4);
				if (k == -1) {
					SL_gen->inc();
					M4.swap(SL_gen->s_i(
						SL_gen->s_li() - 1));
					}
				}
			}
		
		i++;
		}

#ifdef DEBUG_GL2SL_GO
	printf("GL2SL: SL group order:\n");
	go(SL_gen);
#endif

	return OK;
}

static INT trace(VECTOR_OP SVlast, VECTOR_OP SVgen, 
	VECTOR_OP G, INT i, SYM_OP p)
{
	SYM_OB p1;
	INT ii, prev, g;
	
	if (SVlast->s_li() < i)
		return error("trace: SVlast->s_li() < i");
	G->s_i(0)->copy(p);
	p->one();
	ii = i;
	while (TRUE) {
		prev = SVlast->s_ii(ii);
		g = SVgen->s_ii(ii);
		if (prev == -1)
			return OK;
		G->s_i(g)->mult(p, &p1);
		p1.swap(p);
		ii = prev;
		}
}

#if TEXDOCU
INT perm_rep_print_vectors(PERM_REP *P)
#else
Prints all vectors of the representation.
#endif
{
	INT i, ii, j, a, deg;
	INT vec[1000];
	
	deg = P->deg;
	printf("degree = %ld\n", P->deg);
	for (i = 0; i < deg; i++) {
		printf("%ld: ", i);
		fflush(stdout);
		perm_rep_i2vec(P, i, vec);
		perm_rep_vec2i(P, vec, &ii);
		z_vec_i_print(vec, P->n);
		fflush(stdout);
		printf("\n");
		fflush(stdout);
		if (ii != i)
			return error("perm_rep_print_vectors() i != ii");
		}
	for (i = 0; i < deg; i++) {
		printf("$%ld \\cong ", i);
		fflush(stdout);
		perm_rep_i2vec(P, i, vec);
		perm_rep_vec2i(P, vec, &ii);
		printf("(", i);
		for (j = 0; j < P->n; j++) {
			a = vec[j];
			printf("%ld", a);
			if (j < P->n - 1)
				printf(", ");
			}
		printf(")^t$\\\\\n");
		}
	return OK;
}

#if TEXDOCU
INT perm_rep_mult(PERM_REP *P, INT x, INT y)
#else
Multiplication of field elements.
#endif
{
	// printf("perm_rep_mult x=%ld y=%ld\n", x, y); fflush(stdout);
	if (P->f > 1)
		return z_mult_num(P->Z, x, y);
	else {
		return Asr(x * y, P->p);
		}
}

#if TEXDOCU
INT perm_rep_inverse(PERM_REP *P, INT x)
#else
Inversion of the field element x.
#endif
{
	INT y;
	
	if (P->f > 1)
		return z_inverse_num(P->Z, x);
	else {
		y = Inverse_mod(x, P->p);
		return Asr(y, P->p);
		}
}

#if TEXDOCU
INT perm_rep_negate(PERM_REP *P, INT x)
#else
Additive inversion of the field element.
#endif
{
	INT y;
	
	if (P->f > 1)
		return z_negate_num(P->Z, x);
	else {
		y = -x;
		return Asr(y, P->p);
		}
}

#if TEXDOCU
INT perm_rep_is_zero(PERM_REP *P, INT x)
#else
Tests if the field element x is zero.
#endif
{
	INT y;
	
	if (P->f > 1)
		return z_is_zero_num(P->Z, x);
	else {
		y = Asr(x, P->p);
		return (y == 0);
		}
}

#if TEXDOCU
INT perm_rep_zero(PERM_REP *P)
#else
Returns the zero field element.
#endif
{
	if (P->f > 1)
		return z_zero_num(P->Z);
	else {
		return 0;
		}
}

#if TEXDOCU
INT perm_rep_is_one(PERM_REP *P, INT x)
#else
Tests if the field element is one.
#endif
{
	INT y;
	
	if (P->f > 1)
		return z_is_one_num(P->Z, x);
	else {
		y = Asr(x, P->p);
		return (y == 1);
		}
}

#if TEXDOCU
INT perm_rep_one(PERM_REP *P)
#else
Returns the one field element.
#endif
{
	if (P->f > 1)
		return z_one_num(P->Z);
	else {
		return 1;
		}
}

#if TEXDOCU
INT perm_rep_mone(PERM_REP *P)
#else
Returns the field element minus one.
#endif
{
	if (P->f > 1)
		return z_mone_num(P->Z);
	else {
		if (P->p != 2)
			return Asr(P->p - 1, P->p);
		else
			return 1;
		}
}

#if TEXDOCU
INT perm_rep_add(PERM_REP *P, INT x, INT y)
#else
add field elements.
#endif
{
	INT z;
	
	if (P->f > 1)
		return z_add_num(P->Z, x, y);
	else {
		z = Asr(x + y, P->p);
		return z;
		}
}

#if TEXDOCU
INT perm_rep_power_int(PERM_REP *P, INT x, INT l)
#else
raises $x$ to the $l$-th power.
#endif
{
	INT a, b;
	
	// printf("perm_rep_power_int x=%ld l=%ld\n", x, l); fflush(stdout);
	if (l == 0) {
		return perm_rep_one(P);
		}
	a = perm_rep_one(P);
	b = x;
	if (l < 0) {
		b = perm_rep_inverse(P, b);
		l = -l;
		}
	while (l != 0) {
		if (ODD(l)) {
			/* Rule ii): 
			 * a := a * b; l--; */
			a = perm_rep_mult(P, a, b);
			l--;
			continue; /* l kann 0 geworden sein. */
			}
		/* now: EVEN(l) and l != 0 */
		/* Rule iii): 
		 * b := b * b; l := l / 2; */
		b = perm_rep_mult(P, b, b);
		l >>= 1;
		}
	return a;
}

#if TEXDOCU
INT perm_rep_gauss_n_m_rectangular(PERM_REP *P, 
	INT f_special, INT f_complete, MATRIX_OP M, VECTOR_OP base_cols, INT f_v)
#else
returns the rank.
f\_special TRUE: multiply (the row) with the inverse of the pivot to get the pivot to 1.
OUTPUT: coefficients absolutely $\le p/2$ (eventually negative) 
#endif
{
	INT n, m;
	INT i, j, r, k, jj;
	INT a, a1, a_inv, b, c, piv, piv_inv;
	
	base_cols->m_il(0);
	n = M->s_hi();
	m = M->s_li();
	i = 0;
	for (j = 0; j < m; j++) {
	
		/* search for pivot element: */
		for (k = i; k < n; k++) {
			if (!perm_rep_is_zero(P, M->s_iji(k, j))) {
				/* pivot element found: */
				if (k != i) {
					for (jj = j; jj < m; jj++) {
						M->s_ij(i, jj)->swap(M->s_ij(k, jj));
						}
					}
				break;
				} /* if != 0 */
			} /* next k */
		
		if (k == n) /* no pivot found */
			continue; /* increase j, leave i constant */
		
		if (f_v) {
			printf("i=%ld, pivot in row %ld col %ld\n", i, k, j);
			}
		base_cols->inc();
		base_cols->m_ii(i, j);
		if (!f_special) {
			/* make pivot to 1: */
			a = M->s_iji(i, j);
			a_inv = perm_rep_inverse(P, a);
			for (jj = j; jj < m; jj++) {
				b = perm_rep_mult(P, a_inv, M->s_iji(i, jj));
				M->m_iji(i, jj, b);
				}
			if (f_v) {
				printf("pivot=%ld pivot_inv=%ld made to one:%ld\n", 
					a, a_inv, M->s_iji(i, j));
				M->fprint_raw(stdout);
				}
			}
		else {
			piv = M->s_iji(i, j);
			piv_inv = perm_rep_inverse(P, piv);
			}
		
		/* do the gaussian elimination: */
		for (k = i + 1; k < n; k++) {
			a1 = M->s_iji(k, j);
			if (perm_rep_is_zero(P, a1))
				continue;
			M->m_iji(k, j, perm_rep_zero(P));
			for (jj = j + 1; jj < m; jj++) {
				a = M->s_iji(i, jj);
				b = M->s_iji(k, jj);
				if (f_special) {
					a = perm_rep_mult(P, a, piv_inv);
					}
				c = perm_rep_mult(P, a1, a);
				c = perm_rep_negate(P, c);
				c = perm_rep_add(P, b, c);
				M->m_iji(k, jj, c);
				}
			}
		if (f_v) {
			printf("after elimination\n");
			M->fprint_raw(stdout);
			}
		i++;
		} /* next j */
	r = i;
	if (f_complete) {
		for (i = r - 1; i >= 0; i--) {
			j = base_cols->s_ii(i);
			if (!f_special) {
				a = M->s_iji(i, j);
				}
			else {
				piv = M->s_iji(i, j);
				piv_inv = perm_rep_inverse(P, piv);
				}
			/* do the gaussian elimination in the upper part: */
			for (k = i - 1; k >= 0; k--) {
				a1 = M->s_iji(k, j);
				if (perm_rep_is_zero(P, a1))
					continue;
				M->m_iji(k, j, perm_rep_zero(P));
				for (jj = j + 1; jj < m; jj++) {
					a = M->s_iji(i, jj);
					b = M->s_iji(k, jj);
					if (f_special) {
						a = perm_rep_mult(P, a, piv_inv);
						}
					c = perm_rep_mult(P, a1, a);
					c = perm_rep_negate(P, c);
					c = perm_rep_add(P, b, c);
					M->m_iji(k, jj, c);
					}
				}
			}
		}
	return r;
}

#if TEXDOCU
INT perm_rep_get_kernel(PERM_REP *P, MATRIX_OP M, MATRIX_OP K, VECTOR_OP base_cols)
#else
#endif
{
	INT r, n, m, k, i, j, ii, iii, a, aa, b;
	VECTOR_OB kcol;
	
	n = M->s_hi();
	m = M->s_li();
	r = base_cols->s_li();
	k = m - r;
	K->m_ilih_n(k, m);
	kcol.m_il(k);
	ii = 0;
	j = 0;
	if (j < r)
		b = base_cols->s_ii(j);
	else
		b = -1;
	for (i = 0; i < m; i++) {
		if (i == b) {
			j++;
			if (j < r)
				b = base_cols->s_ii(j);
			else
				b = -1;
			}
		else {
			kcol.m_ii(ii, i);
			ii++;
			}
		}
	if (ii != k)
		return error("perm_rep_get_kernel ii != k");
	printf("kcol=");
	kcol.println();
	fflush(stdout);
	ii = 0;
	j = 0;
	if (j < r)
		b = base_cols->s_ii(j);
	else
		b = -1;
	for (i = 0; i < m; i++) {
		if (i == b) {
			for (iii = 0; iii < k; iii++) {
				a = kcol.s_ii(iii);
				aa = M->s_iji(j, a);
				K->m_iji(i, iii, aa);
				}
			j++;
			if (j < r)
				b = base_cols->s_ii(j);
			else
				b = -1;
			}
		else {
			for (iii = 0; iii < k; iii++) {
				if (iii == ii) {
					aa = perm_rep_mone(P);
					}
				else {
					aa = perm_rep_zero(P);
					}
				K->m_iji(i, iii, aa);
				}
			ii++;
			}
		}
	return OK;
}


#endif /* GFQ_TRUE  */




