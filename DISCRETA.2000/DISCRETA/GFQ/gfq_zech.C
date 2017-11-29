/* gfq_zech.C 
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

#undef DEBUG_REP_ZQ
#undef DEBUG_REP_ZP
#undef DEBUG_MTX_DETERMINANTE
#undef DEBUG_MTX_INVERSE

#if TEXDOCU
ZECH_DATA *zech_open(INT deg, INT chi, INT f_verbose)
#else
opens a ZECH\_DATA structure for a field extension of the finite field $GF(p)$.
#endif
{
	ZECH_DATA *Z;
	G_POLYNOM a, b, c;
	INT l, idx;
	BYTE str[1024];
	BYTE str1[1024];

	Z = (ZECH_DATA *) my_malloc(sizeof(ZECH_DATA), "zech_open");

	Z->chi = chi;
	Z->deg = deg;
	Z->q = i_power_j(chi, deg);
	Z->m_order = Z->q - 1;
	Z->m = g_singer(deg, chi);
	Z->idx_zero = Z->m_order /* + 1 */;
	if (EVEN(Z->m_order)) {
		/* ODD(chi) */
		Z->idx_m1 = Z->m_order >> 1;
		}
	else {
		Z->idx_m1 = 0;
		}

	/* alpha := the polynomial X */
	Z->alpha = g_alloc(deg);
	g_set_deg(Z->alpha, 1);
	Z->alpha[0] = 0;
	Z->alpha[1] = 1;

	Z->V = (G_POLYNOM *) my_malloc(sizeof(G_POLYNOM *) * (Z->m_order + 1), "zech_open V");
	Z->Z = (INT *) my_malloc(sizeof(INT) * (Z->m_order + 1), "zech_open Z");
	Z->Num = (INT *) my_malloc(sizeof(INT) * (Z->q + 1), "zech_open Num");
	Z->Num_inv = (INT *) my_malloc(sizeof(INT) * (Z->q + 1), "zech_open Num_inv");
	Z->Frob = (INT *) my_malloc(sizeof(INT) * (Z->q + 1), "zech_open Frob");
	Z->Frob_inv = (INT *) my_malloc(sizeof(INT) * (Z->q + 1), "zech_open Frob_inv");

	// a = g_copy(Z->alpha);
	// b = g_copy(Z->alpha);
	a = g_alloc(deg);
	b = g_alloc(deg);
	c = g_alloc(deg);
	g_move(Z->alpha, a);
	g_move(Z->alpha, b);
	if (f_verbose || TRUE) {
		printf("init ZECH for chi = %ld deg = %ld\n", chi, deg);
		printf("m = ");
		g_print(Z->m);
		}
	for (l = 1; l <= Z->m_order; l++) {
		if (f_verbose) {
			printf("X^%ld = ", l);
			if (l == Z->m_order)
				printf("X^0 = ");
			g_print(a);
			}
		if (l < Z->m_order) {
			if (g_is_one(a, chi)) {
				error("zech_open(): ERROR: m not primitive !!!");
				return NIL;
				}
			Z->V[l] = g_copy(a);
			}
		else {
			if (!g_is_one(a, chi)) {
				error("zech_open(): ERROR: Fermat violated, help !!!");
				return NIL;
				}
			Z->V[0] = g_copy(a);
			}
		g_move(Z->alpha, b);
		g_mult_mod(a, b, c, Z->m, chi);
		g_move(c, a);
		}

	for (l = 0; l < Z->m_order; l++) {
		g_move(Z->V[l], a);
		a[0] = g_asr(a[0] + 1, chi);
		idx = z_search(Z, a);
		if (idx < 0) {
			error("zech_open(): ERROR: u + 1 not found !");
			return NIL;
			}
		else {
			Z->Z[l] = idx;
			}
		}
	if (f_verbose) {
		printf("numeric elements (num, zlog, polynomial):\n");
		}
	for (l = 0; l < Z->q; l++) {
		g_numeric_pol(l, a, chi);
		str[0] = 0;
		g_sprint(a, str);
		idx = z_search(Z, a);
			/* can be idx_zero = m_order + 1 */
		if (idx < 0) {
			printf("l = %ld a = ", l);
			g_print(a);
			error("zech_open(): ERROR: numeric element not found !");
			return NIL;
			}
		else {
			Z->Num[l] = idx;
			Z->Num_inv[idx] = l;
			}
		if (f_verbose) {
			printf("(%3ld,%3ld,%s)\n", l, idx, str);
			}
		}
	if (f_verbose) {
		printf("Num[i] = zechlog_X(p_n X^n + \\ldots + p_1 X + p_0), \n");
		printf("where (p_n, \\ldots , p_0) is the "
			"p-adic representation of i\n");
		printf("Num_inv = Num^{-1}\n");
		printf("(");
		for (l = 0; l < Z->q; l++) {
			printf("%3ld ", Z->Num[l]);
			}
		printf(")\n");
		printf("(");
		for (l = 0; l < Z->q; l++) {
			printf("%3ld ", Z->Num_inv[l]);
			}
		printf(")\n\n");
		fflush(stdout);

		printf("frobenius map:\n");
		}
	for (l = 0; l < Z->q; l++) {
		str[0] = 0;
		str1[0] = 0;
		g_numeric_pol(l, a, chi);
		g_sprint(a, str);
		g_power_mod_apply(a, chi, Z->m, chi);
		g_sprint(a, str1);
		idx = z_search(Z, a);
			/* can be idx_zero = m_order (+ 1) */
		if (idx < 0) {
			printf("l = %ld a = ", l);
			g_print(a);
			error("zech_open(): ERROR: "
			"Frobenius image element not found !");
			return NIL;
			}
		else {
			idx = Z->Num_inv[idx];
			Z->Frob[l] = idx;
			Z->Frob_inv[idx] = l;
			}
		if (f_verbose) {
			printf("%ld -> %ld ; %s -> %s\n", l, idx, str, str1);
			}
		}
	if (f_verbose) {
		printf("Frob and Frob_inv (remember: "
			"deg = 2 => Frob = Frob^{-1}):\n");
		z_print_perm(Z->Frob, Z->q);
		printf("\n");
		z_print_perm(Z->Frob_inv, Z->q);
		printf("\n");
		fflush(stdout);
		}

	g_free(a);
	g_free(b);
	g_free(c);
	return Z;
}

#if TEXDOCU
void zech_free(ZECH_DATA *Z)
#else
Closes a ZECH\_DATA structure.
#endif
{
	INT l;
	
	g_free(Z->m);
	g_free(Z->alpha);
	for (l = 0; l < Z->m_order; l++) {
		g_free(Z->V[l]);
		}
	Z->m_order = 0;
	Z->q = 0;

	my_free(Z->V);
	my_free(Z->Z);
	my_free(Z->Num);
	my_free(Z->Num_inv);
	my_free(Z->Frob);
	my_free(Z->Frob_inv);
	my_free(Z);
}

#if TEXDOCU
void z_print_elt_num_verbose(ZECH_DATA *Z, INT i)
#else
#endif
{
	BYTE str[256];
	G_POLYNOM a;

	a = g_alloc(Z->deg);
	str[0] = 0;
	g_numeric_pol(i, a, Z->chi);
	g_sprint(a, str);
	printf("(%3ld,%3ld,%s)\n", i, Z->Num[i], str);
	g_free(a);
}

#if TEXDOCU
INT z_mult_zlog(ZECH_DATA *Z, INT i1, INT i2)
#else
#endif
{
	INT i3;
	
	if (i1 == Z->idx_zero || i2 == Z->idx_zero)
		return Z->idx_zero;
	i3 = i1 + i2;
	if (i3 >= Z->m_order)
		i3 -= Z->m_order;
	return i3;
}

#if TEXDOCU
INT z_mult_num(ZECH_DATA *Z, INT i1, INT i2)
#else
#endif
{
	INT ii1, ii2, ii3;
	
	ii1 = Z->Num[i1];
	ii2 = Z->Num[i2];
	if (ii1 == Z->idx_zero || ii2 == Z->idx_zero)
		return Z->Num_inv[Z->idx_zero];
	ii3 = ii1 + ii2;
	if (ii3 >= Z->m_order)
		ii3 -= Z->m_order;
	return Z->Num_inv[ii3];
}

#if TEXDOCU
INT z_inverse_zlog(ZECH_DATA *Z, INT i)
#else
#endif
{
	if (i == Z->idx_zero) {
		error("cannot invert 0");
		return 0;
		}
	if (i == 0)
		return 0;
	return Z->m_order - i;
}

#if TEXDOCU
INT z_inverse_num(ZECH_DATA *Z, INT i)
#else
#endif
{
	INT ii;
	
	ii = Z->Num[i];
#ifdef DEBUG_Z_INVERSE_NUM
	printf("z_inverse_num(): i = %ld ii = %ld\n", i, ii);
	fflush(stdout);
#endif
	ii = z_inverse_zlog(Z, ii);
#ifdef DEBUG_Z_INVERSE_NUM
	printf("inverse_zlog = %ld = ", ii);
	fflush(stdout);
#endif
	ii = Z->Num_inv[ii];
#ifdef DEBUG_Z_INVERSE_NUM
	printf("%ld\n", ii);
	fflush(stdout);
#endif
	return ii;
}

#if TEXDOCU
INT z_negate_zlog(ZECH_DATA *Z, INT i)
#else
#endif
{
	return z_mult_zlog(Z, i, Z->idx_m1);
}

#if TEXDOCU
INT z_negate_num(ZECH_DATA *Z, INT i)
#else
#endif
{
	INT ii;
	
	ii = Z->Num[i];
	ii = z_negate_zlog(Z, ii);
	return Z->Num_inv[ii];
}

#if TEXDOCU
INT z_is_zero_zlog(ZECH_DATA *Z, INT i)
#else
#endif
{
	if (i == Z->idx_zero)
		return TRUE;
	else
		return FALSE;
}

#if TEXDOCU
INT z_is_zero_num(ZECH_DATA *Z, INT i)
#else
#endif
{
	INT ii;
	
	ii = Z->Num[i];
	if (ii == Z->idx_zero)
		return TRUE;
	else
		return FALSE;
}

#if TEXDOCU
INT z_zero_zlog(ZECH_DATA *Z)
#else
#endif
{
	return Z->idx_zero;
}

#if TEXDOCU
INT z_zero_num(ZECH_DATA *Z)
#else
#endif
{
	return Z->Num_inv[ Z->idx_zero ];
}

#if TEXDOCU
INT z_is_one_zlog(ZECH_DATA *Z, INT i)
#else
#endif
{
	if (i == 0)
		return TRUE;
	else
		return FALSE;
}

#if TEXDOCU
INT z_is_one_num(ZECH_DATA *Z, INT i)
#else
#endif
{
	INT ii;
	
	ii = Z->Num[i];
	if (ii == 0)
		return TRUE;
	else
		return FALSE;
}

#if TEXDOCU
INT z_one_zlog(ZECH_DATA *Z)
#else
#endif
{
	return 0;
}

#if TEXDOCU
INT z_one_num(ZECH_DATA *Z)
#else
#endif
{
	return Z->Num_inv[ 0 ];
}

#if TEXDOCU
INT z_mone_zlog(ZECH_DATA *Z)
#else
#endif
{
	return Z->idx_m1;
}

#if TEXDOCU
INT z_mone_num(ZECH_DATA *Z)
#else
#endif
{
	return Z->Num_inv[ Z->idx_m1 ];
}

#if TEXDOCU
INT z_prim_elt_zlog(ZECH_DATA *Z)
#else
#endif
{
	return 1;
}

#if TEXDOCU
INT z_prim_elt_num(ZECH_DATA *Z)
#else
#endif
{
	return Z->Num_inv[1];
}

#if TEXDOCU
INT z_add_zlog(ZECH_DATA *Z, INT i1, INT i2)
#else
#endif
{
	INT idx_zero, tmp, w, i;
	
	idx_zero = Z->idx_zero;
	if (i1 == idx_zero) {
		return i2;
		}
	if (i2 == idx_zero) {
		return i1;
		}
	if (i1 > i2) {
		tmp = i1;
		i1 = i2;
		i2 = tmp;
		}
	/* now i1 <= i2 */
	w = i2 - i1;
	i = Z->Z[w];
	/* i can be the zero symbol: idx_zero */
	if (i == idx_zero) {
		return i;
		}
	i += i1;
	if (i >= Z->m_order)
		i -= Z->m_order;
	return i;
}

#if TEXDOCU
INT z_add_num(ZECH_DATA *Z, INT i1, INT i2)
#else
#endif
{
	INT ii1, ii2, i;
	
	ii1 = Z->Num[i1];
	ii2 = Z->Num[i2];
	i = z_add_zlog(Z, ii1, ii2);
	return Z->Num_inv[i];
}

#if TEXDOCU
INT z_apply_frob_zlog(ZECH_DATA *Z, INT i1)
#else
#endif
{
	INT i;

	if (i1 == Z->idx_zero)
		return Z->idx_zero;
	i = Z->Num_inv[i1];
	i = Z->Frob[i];
	i = Z->Num[i];
	return i;
}

#if TEXDOCU
INT z_apply_frob_num(ZECH_DATA *Z, INT i1)
#else
#endif
{
	INT i;

	i = Z->Frob[i1];
	return i;
}

#if TEXDOCU
INT z_apply_frob2_zlog(ZECH_DATA *Z, INT i1)
#else
#endif
{
	INT i, j, deg1;

	if (i1 == Z->idx_zero)
		return Z->idx_zero;
	deg1 = Z->deg >> 1;
	i = Z->Num_inv[i1];
	for (j = 0; j < deg1; j++)
		i = Z->Frob[i];
	i = Z->Num[i];
	return i;
}

#if TEXDOCU
INT z_apply_frob2_num(ZECH_DATA *Z, INT i1)
#else
#endif
{
	INT i, j, deg1;

	deg1 = Z->deg >> 1;
	i = i1;
	for (j = 0; j < deg1; j++)
		i = Z->Frob[i];
	return i;
}

#if TEXDOCU
INT z_apply_frob_sz_zlog(ZECH_DATA *Z, INT i1)
#else
$x \mapsto x^{2^{m+1}}$
Huppert and Blackburn, 
Finite Groups III~\cite{HuppertBlackburn82b}, p. 182.
#endif
{
	INT i, j, m;

	if (i1 == Z->idx_zero)
		return Z->idx_zero;
	m = (Z->deg - 1) >> 1;
	i = Z->Num_inv[i1];
	for (j = 0; j <= m; j++)
		i = Z->Frob[i];
	i = Z->Num[i];
	return i;
}

#if TEXDOCU
INT z_apply_frob_sz_num(ZECH_DATA *Z, INT i1)
#else
#endif
{
	INT i, j, m;

	m = (Z->deg - 1) >> 1;
	i = i1;
	for (j = 0; j <= m; j++)
		i = Z->Frob[i];
	return i;
}

#if TEXDOCU
INT z_apply_frob_sz0_zlog(ZECH_DATA *Z, INT i1)
#else
$x \mapsto x^{2^m}$
Huppert and Blackburn, 
Finite Groups III~\cite{HuppertBlackburn82b}, p. 183.
#endif
{
	INT i, j, m;

	if (i1 == Z->idx_zero)
		return Z->idx_zero;
	m = (Z->deg - 1) >> 1;
	i = Z->Num_inv[i1];
	for (j = 0; j < m; j++)
		i = Z->Frob[i];
	i = Z->Num[i];
	return i;
}

#if TEXDOCU
INT z_apply_frob_sz0_num(ZECH_DATA *Z, INT i1)
#else
#endif
{
	INT i, j, m;

	m = (Z->deg - 1) >> 1;
	i = i1;
	for (j = 0; j < m; j++)
		i = Z->Frob[i];
	return i;
}

#if TEXDOCU
INT z_trace_zlog(ZECH_DATA *Z, INT a)
#else
#endif
{
	INT i, r;

	r = a;
	for (i = 1; i < Z->deg; i++) {
		a = z_apply_frob_zlog(Z, a);
		r = z_add_zlog(Z, r, a);
		}
	return r;
}

#if TEXDOCU
INT z_trace_num(ZECH_DATA *Z, INT a)
#else
#endif
{
	INT i, r;

	r = a;
	for (i = 1; i < Z->deg; i++) {
		a = z_apply_frob_num(Z, a);
		r = z_add_num(Z, r, a);
		}
	return r;
}

#if TEXDOCU
INT z_trace2_zlog(ZECH_DATA *Z, INT a)
#else
trace over subfield of index 2 (deg needs to be even)
#endif
{
	INT i, r;

	r = a;
	for (i = 1; i < 2; i++) {
		a = z_apply_frob2_zlog(Z, a);
		r = z_add_zlog(Z, r, a);
		}
	return r;
}

#if TEXDOCU
INT z_trace2_num(ZECH_DATA *Z, INT a)
#else
#endif
{
	INT i, r;

	r = a;
	for (i = 1; i < 2; i++) {
		a = z_apply_frob2_num(Z, a);
		r = z_add_num(Z, r, a);
		}
	return r;
}

#if TEXDOCU
INT z_norm_zlog(ZECH_DATA *Z, INT a)
#else
#endif
{
	INT i, r;

	r = a;
	for (i = 1; i < Z->deg; i++) {
		a = z_apply_frob_zlog(Z, a);
		r = z_mult_zlog(Z, r, a);
		}
	return r;
}

#if TEXDOCU
INT z_norm_num(ZECH_DATA *Z, INT a)
#else
#endif
{
	INT i, r;

	r = a;
	for (i = 1; i < Z->deg; i++) {
		a = z_apply_frob_num(Z, a);
		r = z_mult_num(Z, r, a);
		}
	return r;
}

#if TEXDOCU
INT z_norm2_zlog(ZECH_DATA *Z, INT a)
#else
#endif
{
	INT i, r;

	r = a;
	for (i = 1; i < 2; i++) {
		a = z_apply_frob2_zlog(Z, a);
		r = z_mult_zlog(Z, r, a);
		}
	return r;
}

#if TEXDOCU
INT z_norm2_num(ZECH_DATA *Z, INT a)
#else
#endif
{
	INT i, r;

	r = a;
	for (i = 1; i < 2; i++) {
		a = z_apply_frob2_num(Z, a);
		r = z_mult_num(Z, r, a);
		}
	return r;
}

#if TEXDOCU
void z_print_perm(INT *p, INT len)
#else
#endif
{
	BYTE str[1024];

	str[0] = 0;
	z_sprint_perm(p, len, str);
	printf("%s", str);
}

#if TEXDOCU
INT z_sprint_perm(INT *p, INT length, BYTE *str)
#else
#endif
{
	BYTE have_seen[MAX_PERM_N + 1];
	INT l, l1, first, next, len, N;
	INT f_nothing_printed_at_all = TRUE;
	
	if (str == NIL) {
		return error("z_sprint_perm(): args NIL");
		}
	N = length;
	if (N > MAX_PERM_N) {
		return error("z_sprint_perm(): N > MAX_PERM_N");
		}
	for (l = 0; l < N; l++) {
		have_seen[l] = FALSE;
		}
	l = 0;
	while (TRUE) {
		if (l >= N) {
			if (f_nothing_printed_at_all) {
				sprintf(Eostr(str), "id");
				}
			return(OK);
			}
		if (have_seen[l]) {
			l++;
			continue;
			}
		/* Bearbeite Zyklus, beginnend mit l: */
		first = l;
		l1 = l;
		len = 1;
		while (TRUE) {
			have_seen[l1] = TRUE;
			next = p[l1];
			if (next > MAX_PERM_N) {
				printf("l1 = %ld next = %ld\n", 
				l1, next);
				return error("z_sprint_perm(): "
				"next > MAX_PERM_N");
				}
			if (next == first) {
				break;
				}
			if (have_seen[next]) {
				printf("first = %ld l1 = %ld "
					"next = %ld length = %ld\n", 
					first, l1, next, length);
				printf("%s\n", str);
				fflush(stdout);
				return error("z_sprint_perm() " 
				"have_seen[next]");
				}
			l1 = next;
			len++;
			}
		if (len == 1)
			continue;
		f_nothing_printed_at_all = FALSE;
		/* Drucke Zyklus, beginnend mit first: */
		l1 = first;
		sprintf(Eostr(str), "(");
		while (TRUE) {
			sprintf(Eostr(str), "%ld", l1);
			next = p[l1];
			if (next == first) {
				break;
				}
			sprintf(Eostr(str), " ");
			l1 = next;
			}
		sprintf(Eostr(str), ")");
		}
	return(OK);
}

#if TEXDOCU
INT z_search(ZECH_DATA *Z, G_POLYNOM p)
#else
works fine even for a = 0 
#endif
{
	INT l;

	if (g_is_zero(p, Z->chi))
		return Z->idx_zero;
	for (l = 0; l < Z->m_order; l++) {
		if (g_cmp(p, Z->V[l]) == 0) {
			return l;
			}
		}
	return -1;
}


#if TEXDOCU
void z_mtx_i_print(INT *mtx, INT dim)
#else
#endif
{
	INT i, j;
	
	for (i = 0; i < dim; i++) {
		printf("(");
		for (j = 0; j < dim; j++) {
			printf("%ld ", mtx[i * dim + j]);
			}
		printf(")\n");
		}
	printf("\n");
}

#if TEXDOCU
void z_mtx_print_nm(INT *mtx, INT n, INT m)
#else
#endif
{
	INT i, j;
	
	for (i = 0; i < n; i++) {
		printf("(");
		for (j = 0; j < m; j++) {
			printf("%ld ", mtx[i * m + j]);
			}
		printf(")\n");
		}
	printf("\n");
}

#if TEXDOCU
void z_vec_i_print(INT *vec, INT dim)
#else
#endif
{
	INT j, a;
	
	printf("(");
	for (j = 0; j < dim; j++) {
		a = vec[j];
		printf("%ld ", a);
		}
	printf(")");
}

#if TEXDOCU
void z_mtx3_print(ZECH_DATA *Z, INT *mtx)
#else
#endif
{
	z_mtx_i_print(mtx, 3);
}

#if TEXDOCU
void z_mtx4_print(ZECH_DATA *Z, INT *mtx)
#else
#endif
{
	z_mtx_i_print(mtx, 4);
}

#if TEXDOCU
void z_mtx_i_gen_diag_num(ZECH_DATA *Z, INT *mtx, INT *x, INT dim)
#else
#endif
{
	INT i, j, *p;

	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) {
			p = mtx + i * dim + j;
			if (i != j)
				*p = z_zero_num(Z);
			else
				*p = x[i];
			}
		}
}

#if TEXDOCU
void z_mtx3_gen_diag_num(ZECH_DATA *Z, INT *mtx, INT *x)
#else
#endif
{
	z_mtx_i_gen_diag_num(Z, mtx, x, 3);
}

#if TEXDOCU
void z_mtx4_gen_diag_num(ZECH_DATA *Z, INT *mtx, INT *x)
#else
#endif
{
	z_mtx_i_gen_diag_num(Z, mtx, x, 4);
}

/* 
 * affine and projective 
 * permutation representations
 */

#if TEXDOCU
INT z_mtx_vec_mult_GFq_num(ZECH_DATA *Z, INT *mtx, INT n, 
	INT *vec1, INT *vec2)
#else
#endif
{
	INT i, j, a, b, mij, vj;
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			mij = mtx[i * n + j];
			vj = vec1[j];
			a = z_mult_num(Z, mij, vj);
			if (j == 0)
				b = a;
			else
				b = z_add_num(Z, b, a);
			}
		vec2[i] = b;
		}
	return OK;
}

#if TEXDOCU
INT z_mtx_vec_mult_GFp(INT p, INT *mtx, INT n, 
	INT *vec1, INT *vec2)
#else
#endif
/* result is POSITIVE representative system mod p */
{
	INT i, j, a, b, mij, vj;
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			mij = mtx[i * n + j];
			vj = vec1[j];
			a = Asr(mij * vj, p);
			if (j == 0)
				b = a;
			else
				b = Asr(b + a, p);
			}
		if (b < 0)
			b += p;
		vec2[i] = b;
		}
	return OK;
}

#if TEXDOCU
INT z_vec_norm_projective_GFq_num(ZECH_DATA *Z, INT *vec, INT n)
#else
#endif
{
	INT idx_zero, idx_one, l, i, a, b;
	
	idx_zero = 0;
	idx_one = Z->Num_inv[0];
	for (l = n - 1; l >= 0; l--) {
		if (vec[l] == idx_zero)
			continue;
		a = z_inverse_num(Z, vec[l]);
		vec[l] = idx_one;
		for (i = 0; i < l; i++) {
			b = vec[i];
			if (b == idx_zero)
				continue;
			vec[i] = z_mult_num(Z, b, a);
			}
		return OK;
		}
	return error("z_vec_norm_projective_GFq_num() zero vector !");
}

#if TEXDOCU
INT z_vec_norm_projective_GFp(INT p, INT *vec, INT n)
#else
#endif
{
	INT idx_zero, idx_one, l, i, a, b, c;
	
	idx_zero = 0;
	idx_one = 1;
	for (l = n - 1; l >= 0; l--) {
		if (vec[l] == idx_zero)
			continue;
		a = Inverse_mod(vec[l], p);
		vec[l] = idx_one;
		for (i = 0; i < l; i++) {
			b = vec[i];
			if (b == idx_zero)
				continue;
			c = Asr(b * a, p);
			if (c < 0)
				c += p;
			vec[i] = c;
			}
		return OK;
		}
	return error("z_vec_norm_projective_GFp() zero vector !");
}

#if TEXDOCU
INT z_mtx_affine_rep_Zq_num(ZECH_DATA *Z, INT *mtx, INT n, 
	PERMUTATION_OP perm)
#else
#endif
{
	return z_mtx_rep_Zq_num(TRUE, Z, mtx, n, perm);
}

#if TEXDOCU
INT z_mtx_projective_rep_Zq_num(ZECH_DATA *Z, INT *mtx, INT n, 
	PERMUTATION_OP perm)
#else
#endif
{
	return z_mtx_rep_Zq_num(FALSE, Z, mtx, n, perm);
}

#if TEXDOCU
INT z_mtx_rep_Zq_num(INT f_affine, ZECH_DATA *Z, INT *mtx, INT n, 
	PERMUTATION_OP perm)
#else
#endif
{
	INT q, deg, i, ii;
	INT idx_zero, idx_one;
	INT *vec1, *vec2;
	
	q = Z->q;
	vec1 = (INT *) my_malloc(n * sizeof(INT), "z_mtx_rep_Zq_num");
	vec2 = (INT *) my_malloc(n * sizeof(INT), "z_mtx_rep_Zq_num");
	if (f_affine)
		deg = z_affine_degree(q, n);
	else {
		deg = z_projective_degree(q, n);
		idx_zero = 0;
		idx_one = Z->Num_inv[0];
		}
#ifdef DEBUG_REP_ZQ
	printf("z_mtx_rep_Zq_num: f_affine = %ld\n", f_affine);
	z_mtx_i_print(mtx, n);
#endif
	perm->m_il(deg);
	for (i = 0; i < deg; i++) {
		if (f_affine)
			z_affine_i2vec(q, n, i, vec1);
		else 
			z_projective_i2vec(q, idx_zero, idx_one, n, i, vec1);

#ifdef DEBUG_REP_ZQ
		z_vec_i_print(vec1, n);
		fflush(stdout);
#endif
		z_mtx_vec_mult_GFq_num(Z, mtx, n, vec1, vec2);
		if (!f_affine)
			z_vec_norm_projective_GFq_num(Z, vec2, n);
#ifdef DEBUG_REP_ZQ
		z_vec_i_print(vec2, n);
		fflush(stdout);
#endif
		if (f_affine)
			z_affine_vec2i(q, n, vec2, &ii);
		else 
			z_projective_vec2i(q, idx_zero, idx_one, n, vec2, &ii);
#ifdef DEBUG_REP_ZQ
		printf(" %ld\n", ii);
		fflush(stdout);
#endif
		perm->m_ii(i, ii + 1);
		}
	my_free(vec1);
	my_free(vec2);
#ifdef DEBUG_REP_ZQ
	perm->println();
#endif
	return OK;
}

#if TEXDOCU
INT z_mtx_affine_rep_Zp(INT p, INT *mtx, INT n, 
	PERMUTATION_OP perm)
#else
#endif
{
	return z_mtx_rep_Zp(TRUE, p, mtx, n, perm);
}

#if TEXDOCU
INT z_mtx_projective_rep_Zp(INT p, INT *mtx, INT n, 
	PERMUTATION_OP perm)
#else
#endif
{
	return z_mtx_rep_Zp(FALSE, p, mtx, n, perm);
}

#if TEXDOCU
INT z_mtx_rep_Zp(INT f_affine, INT p, INT *mtx, INT n, 
	PERMUTATION_OP perm)
#else
#endif
{
	INT deg, i, ii;
	INT idx_zero, idx_one;
	INT *vec1, *vec2;
	
	vec1 = (INT *) my_malloc(n * sizeof(INT), "z_mtx_rep_Zp");
	vec2 = (INT *) my_malloc(n * sizeof(INT), "z_mtx_rep_Zp");
	if (f_affine)
		deg = z_affine_degree(p, n);
	else {
		idx_zero = 0;
		idx_one = 1;
		deg = z_projective_degree(p, n);
		}
#ifdef DEBUG_REP_ZP
	printf("z_mtx_rep_Zp: f_affine = %ld\n", f_affine);
	z_mtx_i_print(mtx, n);
#endif
	perm->m_il(deg);
	for (i = 0; i < deg; i++) {
		if (f_affine)
			z_affine_i2vec(p, n, i, vec1);
		else
			z_projective_i2vec(p, idx_zero, idx_one, n, i, vec1);

#ifdef DEBUG_REP_ZP
		z_vec_i_print(vec1, n);
		fflush(stdout);
#endif
		z_mtx_vec_mult_GFp(p, mtx, n, vec1, vec2);
		if (!f_affine)
			z_vec_norm_projective_GFp(p, vec2, n);
#ifdef DEBUG_REP_ZP
		z_vec_i_print(vec2, n);
		fflush(stdout);
#endif
		if (f_affine)
			z_affine_vec2i(p, n, vec2, &ii);
		else
			z_projective_vec2i(p, idx_zero, idx_one, n, vec2, &ii);
#ifdef DEBUG_REP_ZP
		printf(" %ld\n", ii);
		fflush(stdout);
#endif
		perm->m_ii(i, ii + 1);
		}
	my_free(vec1);
	my_free(vec2);
#ifdef DEBUG_REP_ZP
	perm->println();
#endif
	return OK;
}

#if TEXDOCU
INT z_Frobenius_rep_Zq_num(INT f_affine, ZECH_DATA *Z, INT n, 
	PERMUTATION_OP perm)
#else
#endif
{
	INT q, deg, i, j, ii, a, b;
	INT idx_zero, idx_one;
	INT *vec1, *vec2;
	
	q = Z->q;
	vec1 = (INT *) my_malloc(n * sizeof(INT), "z_Frobenius_rep_Zq_num");
	vec2 = (INT *) my_malloc(n * sizeof(INT), "z_Frobenius_rep_Zq_num");
	if (f_affine)
		deg = z_affine_degree(q, n);
	else {
		deg = z_projective_degree(q, n);
		idx_zero = 0;
		idx_one = Z->Num_inv[0];
		}
#ifdef DEBUG_REP_ZQ
	printf("z_Frobenius_rep_Zq_num: f_affine = %ld\n", f_affine);
#endif
	perm->m_il(deg);
	for (i = 0; i < deg; i++) {
		if (f_affine)
			z_affine_i2vec(q, n, i, vec1);
		else 
			z_projective_i2vec(q, idx_zero, idx_one, n, i, vec1);

#ifdef DEBUG_REP_ZQ
		z_vec_i_print(vec1, n);
		fflush(stdout);
#endif
		for (j = 0; j < n; j++) {
			a = vec1[j];
			b = z_apply_frob_num(Z, a);
			vec2[j] = b;
			}
#ifdef DEBUG_REP_ZQ
		z_vec_i_print(vec2, n);
		fflush(stdout);
#endif
		if (f_affine)
			z_affine_vec2i(q, n, vec2, &ii);
		else 
			z_projective_vec2i(q, idx_zero, idx_one, n, vec2, &ii);
#ifdef DEBUG_REP_ZQ
		printf(" %ld\n", ii);
		fflush(stdout);
#endif
		perm->m_ii(i, ii + 1);
		}
	my_free(vec1);
	my_free(vec2);
#ifdef DEBUG_REP_ZQ
	perm->println();
#endif
	return OK;
}

#if TEXDOCU
INT z_translation_rep_Zq_num(ZECH_DATA *Z, INT n, INT ei, INT betaj, 
	PERMUTATION_OP perm)
#else
#endif
/* translation: 
 * x * ei \mapsto (x + betaj) * ei , 
 * ei the i-th basis vector (0 \le i < n), 
 * betaj the j-th basis vector in the basis for GFq over GFp, 
 * (0 \le betaj < deg). */
{
	INT q, p, p1, deg, i, j, ii, a;
	INT *vec1, *vec2;
	
	q = Z->q;
	p = Z->chi;
	if (betaj >= Z->deg)
		return error("z_translation_rep_Zq_num(): betaj >= Z->deg");
	p1 = i_power_j(p, betaj);
	if (ei >= n)
		return error("z_translation_rep_Zq_num(): ei >= n");
	vec1 = (INT *) my_malloc(n * sizeof(INT), "z_translation_rep_Zq_num");
	vec2 = (INT *) my_malloc(n * sizeof(INT), "z_translation_rep_Zq_num");
	deg = z_affine_degree(q, n);
#ifdef DEBUG_REP_ZQ
	printf("z_translation_rep_Zq_num: ei = %ld betaj = %ld\n", ei, betaj);
#endif
	perm->m_il(deg);
	for (i = 0; i < deg; i++) {
		z_affine_i2vec(q, n, i, vec1);

#ifdef DEBUG_REP_ZQ
		z_vec_i_print(vec1, n);
		fflush(stdout);
#endif
		for (j = 0; j < n; j++) {
			a = vec1[j];
			if (j == ei)
				a = z_add_num(Z, a, p1);
			vec2[j] = a;
			}
#ifdef DEBUG_REP_ZQ
		z_vec_i_print(vec2, n);
		fflush(stdout);
#endif
		z_affine_vec2i(q, n, vec2, &ii);
#ifdef DEBUG_REP_ZQ
		printf(" %ld\n", ii);
		fflush(stdout);
#endif
		perm->m_ii(i, ii + 1);
		}
	my_free(vec1);
	my_free(vec2);
#ifdef DEBUG_REP_ZQ
	perm->println();
#endif
	return OK;
}

#if TEXDOCU
INT z_translation_rep_Zp(INT p, INT n, INT ei, 
	PERMUTATION_OP perm)
#else
#endif
/* translation: 
 * x * ei \mapsto (x + 1) * ei , 
 * ei the i-th basis vector (0 \le i < n). */
{
	INT deg, i, j, ii, a;
	INT *vec1, *vec2;
	
	if (ei >= n)
		return error("z_translation_rep_Zp(): ei >= n");
	vec1 = (INT *) my_malloc(n * sizeof(INT), "z_translation_rep_Zp");
	vec2 = (INT *) my_malloc(n * sizeof(INT), "z_translation_rep_Zp");
	deg = z_affine_degree(p, n);
#ifdef DEBUG_REP_ZP
	printf("z_translation_rep_Zp: ei = %ld\n", ei);
#endif
	perm->m_il(deg);
	for (i = 0; i < deg; i++) {
		z_affine_i2vec(p, n, i, vec1);

#ifdef DEBUG_REP_ZP
		z_vec_i_print(vec1, n);
		fflush(stdout);
#endif
		for (j = 0; j < n; j++) {
			a = vec1[j];
			if (j == ei) {
				a++;
				if (a == p)
					a = 0;
				}
			vec2[j] = a;
			}
#ifdef DEBUG_REP_ZP
		z_vec_i_print(vec2, n);
		fflush(stdout);
#endif
		z_affine_vec2i(p, n, vec2, &ii);
#ifdef DEBUG_REP_ZP
		printf(" %ld\n", ii);
		fflush(stdout);
#endif
		perm->m_ii(i, ii + 1);
		}
	my_free(vec1);
	my_free(vec2);
#ifdef DEBUG_REP_ZP
	perm->println();
#endif
	return OK;
}

#if TEXDOCU
INT z_GL_n_q_gen_num(ZECH_DATA *Z, INT ***mtx, INT *len, INT n)
#else
#endif
{
	INT idx_zero = 0; /* the zero element: 0 * X^deg - 1 + ... 0 X + 0 */
	INT idx_one = Z->Num_inv[0]; /* the element 1 = X^0 */
	INT idx_alpha = Z->Num_inv[1]; /* the element X = X^1 */
	
	return z_GL_gen(mtx, len, n, idx_zero, idx_one, idx_alpha);
}

#if TEXDOCU
INT z_GL_n_p_gen(INT ***mtx, INT *len, INT n, INT p)
#else
#endif
{
	INT idx_zero = 0;
	INT idx_one = 1;
	INT idx_alpha = Asr(primitive_root(p), p);
	
	return z_GL_gen(mtx, len, n, idx_zero, idx_one, idx_alpha);	
}

#if TEXDOCU
INT z_GL_gen(INT ***mtx, INT *len, INT n, 
	INT idx_zero, INT idx_one, INT idx_alpha)
#else
#endif
/* alpha = primitive root of the field */
{
	INT i, j, k, k1, k2, l, a;
	INT **pp_mtx, *p_mtx;
	
	l = n + n * (n - 1);
	pp_mtx = (INT **) my_malloc(l * sizeof(INT *), "z_GL_gen");
	if (pp_mtx == NULL)
		return error("z_GL_gen() no memory");
	*len = l;
	*mtx = pp_mtx;
	for (i = 0; i < l; i++) {
		pp_mtx[i] = (INT *) my_malloc(n * n * sizeof(INT), "z_GL_gen pp_mtx[i]");
		if (pp_mtx[i] == NULL)
			return error("z_GL_gen() no memory");
		}
	/* n generators: diagonal matrix with 
	 * X at (i,i) position (X = generator of the 
	 * field (multiplicatively)) */
	for (k = 0; k < n; k++) {
		p_mtx = pp_mtx[k];
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				if (i == j) {
					if (i == k) {
						a = idx_alpha;
						}
					else {
						a = idx_one;
						}
					}
				else {
					a = idx_zero;
					}
				p_mtx[i * n + j] = a;
				}
			}
		}
		
	/* n * (n - 1) off diagonal positions with 1 */
	for (k1 = 0; k1 < n; k1++) {
		for (k2 = 0; k2 < n; k2++) {
			if (k1 == k2)
				continue;
			p_mtx = pp_mtx[k];
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					if (i == j) {
						a = idx_one;
						}
					else {
						if (i == k1 && j == k2) {
							a = idx_one;
							}
						else {
							a = idx_zero;
							}
						}
					p_mtx[i * n + j] = a;
					}
				}
			k++;
			}
		}
	return OK;
}

#if TEXDOCU
INT z_test_GL_perm_rep(INT f_affine, INT f_semi, INT f, INT p, INT n)
#else
#endif
{
	VECTOR_OB V;

	z_GL_n_q_perm_rep(f_affine, f_semi, f, p, n, &V, TRUE);
	return OK;
}

#if TEXDOCU
INT z_GL_n_q_perm_rep_info(INT f_affine, 
	INT f_semi, INT f_special, INT f, INT p, INT n, 
	INT *deg, BYTE *label, BYTE *label_tex)
#else
#endif
{
	INT q;
	BYTE str[64];

	q = i_power_j(p, f);
	if (f_affine)
		*deg = z_affine_degree(q, n);
	else
		*deg = z_projective_degree(q, n);
	label[0] = 0;
	label_tex[0] = 0;
	str[1] = 0;
	if (f_affine)
		str[0] = 'A';
	else
		str[0] = 'P';
	strcpy(label, str);
	strcpy(label_tex, str);

	/* special groups not yet implemented ! */
	if (f_special) {
		if (f_semi) {
			strcat(label, "SS");
			strcat(label_tex, "\\Sigma");
			}
		else {
			strcat(label, "S");
			strcat(label_tex, "S");
			}
		}
	else {
		if (f_semi) {
			strcat(label, "GG");
			strcat(label_tex, "\\Gamma");
			}
		else {
			strcat(label, "G");
			strcat(label_tex, "G");
			}
		}
	
	sprintf(str, "L_%ld_%ld", n, q);
	strcat(label, str);
	sprintf(str, "L_%ld(%ld)", n, q);
	strcat(label_tex, str);
	return OK;
}

#if TEXDOCU
INT z_GL_n_q_perm_rep(INT f_affine, INT f_semi, INT f, INT p, INT n, 
	VECTOR_OP V, INT f_verbose)
#else
#endif
{
	ZECH_DATA *Z;
	INT idx_zero; /* the zero element: 0 * X^deg - 1 + ... 0 X + 0 */
	INT idx_one; /* the element 1 = X^0 */
	INT i, j, deg, q;
	INT **mtx, nb_gen;
	PERMUTATION_OB perm;
	INT V_len;

	if (f_verbose) {
		printf("z_GL_n_q_perm_rep(): "
			"f_affine = %ld f_semi = %ld f = %ld p = %ld n = %ld\n", 
			f_affine, f_semi, f, p, n);
		}
	if (f > 1) {
		Z = zech_open(f, p, FALSE);
		q = Z->q;
		idx_zero = 0;
		idx_one = Z->Num_inv[0];
		z_GL_n_q_gen_num(Z, &mtx, &nb_gen, n);
		}
	else {
		if (f_semi) {
			printf("ignoring f_semi in case f = 1\n");
			f_semi = FALSE;
			}
		q = p;
		idx_zero = 0;
		idx_one = 1;
		z_GL_n_p_gen(&mtx, &nb_gen, n, p);
		}
	if (f_affine)
		deg = z_affine_degree(q, n);
	else
		deg = z_projective_degree(q, n);
	if (f_verbose) {
		printf("found %ld generators for GL(n,p^f) "
			"degree %ld\n", nb_gen, deg);
		fflush(stdout);
		}
	V->m_il(nb_gen);
	V_len = nb_gen;
	for (i = 0; i < nb_gen; i++) {
		if (f_verbose) {
			printf("%ld:\n", i);
			z_mtx_i_print(mtx[i], n);
			}
		if (f > 1) {
			z_mtx_rep_Zq_num(f_affine, Z, mtx[i], n, &perm);
			}
		else {
			z_mtx_rep_Zp(f_affine, p, mtx[i], n, &perm);
			}
		if (f_verbose) {
			perm.println();
			fflush(stdout);
			}
		perm.swap((PERMUTATION_OP) V->s_i(i));
		}
	if (f_semi) {
		z_Frobenius_rep_Zq_num(f_affine, Z, n, &perm);
		if (f_verbose) {
			printf("Frobenius map:\n");
			perm.println();
			fflush(stdout);
			}
		V->inc();
		perm.swap((PERMUTATION_OP) V->s_i(V_len));
		V_len++;
		}
	if (f_affine) {
		for (i = 0; i < n; i++) {
			for (j = 0; j < f; j++) {
				if (f == 1)
					z_translation_rep_Zp(p, n, i, &perm);
				else
					z_translation_rep_Zq_num(Z, n, i, j, &perm);
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
	
	if (f > 1) {
		zech_free(Z);
		}
	for (i = 0; i < nb_gen; i++) {
		my_free(mtx[i]);
		}
	my_free(mtx);
	return OK;
}

#if TEXDOCU
INT z_mtx_determinante_GFp(INT p, INT *mtx, INT n)
#else
#endif
{
	INT *M, i, j, a, det;
	
	M = (INT *) my_malloc(n * n * sizeof(INT), "z_mtx_determinante_GFp");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			M[i * n + j] = mtx[i * n + j];
			}
		}
#ifdef DEBUG_MTX_DETERMINANTE
	printf("z_mtx_determinante_GFp(): before gauss:\n");
	z_mtx_print_nm(M, n, n);
#endif
	z_gauss_n_m_GFp(p, TRUE /* f_special */, M, n, n);
#ifdef DEBUG_MTX_DETERMINANTE
	printf("z_mtx_determinante_GFp(): after gauss:\n");
	z_mtx_print_nm(M, n, n);
#endif
	for (i = 0; i < n; i++) {
		a = M[i * n + i];
		if (i == 0)
			det = a;
		else
			det = Asr(a * det, p);
		}
	
	my_free(M);
	return det;
}

#if TEXDOCU
INT z_mtx_determinante_GFq_num(ZECH_DATA *Z, INT *mtx, INT n)
#else
#endif
{
	INT *M, i, j, a, det;
	
	M = (INT *) my_malloc(n * n * sizeof(INT), "z_mtx_determinante_GFq_num");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			M[i * n + j] = mtx[i * n + j];
			}
		}
#ifdef DEBUG_MTX_DETERMINANTE
	printf("z_mtx_determinante_GFq_num(): before gauss:\n");
	z_mtx_print_nm(M, n, n);
#endif
	z_gauss_n_m_GFq_num(Z, TRUE /* f_special */, M, n, n);
#ifdef DEBUG_MTX_DETERMINANTE
	printf("z_mtx_determinante_GFq_num(): after gauss:\n");
	z_mtx_print_nm(M, n, n);
#endif
	for (i = 0; i < n; i++) {
		a = M[i * n + i];
		if (i == 0)
			det = a;
		else {
			det = z_mult_num(Z, a, det);
			}
		}
	
	my_free(M);
	return det;
}

#if TEXDOCU
INT z_mtx_inverse_GFp(INT p, INT *mtx, INT *mtx_inv, INT n)
#else
#endif
{
	INT *M, i, j, a, N;
	
	N = 2 * n;
	M = (INT *) my_malloc(n * N * sizeof(INT), "z_mtx_inverse_GFp");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			M[i * N + j] = mtx[i * n + j];
			}
		}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == j)
				M[i * N + N + j] = 1;
			else
				M[i * N + N + j] = 0;
			}
		}
#ifdef DEBUG_MTX_INVERSE
	printf("z_mtx_inverse_GFp(): before gauss:\n");
	z_mtx_print_nm(M, n, N);
#endif
	z_gauss_n_m_GFp(p, FALSE /* f_special */, M, n, N);
#ifdef DEBUG_MTX_INVERSE
	printf("z_mtx_inverse_GFp(): after gauss:\n");
	z_mtx_print_nm(M, n, N);
#endif
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a = M[i * N + N + j];
			if (a < 0)
				a += p;
			if (a < 0 || a >= p)
				return error("z_mtx_inverse_GFp() a < 0 || a >= p");
			mtx_inv[i * n + j] = a;
			}
		}
	my_free(M);
	return OK;
}

#if TEXDOCU
INT z_mtx_inverse_GFq_num(ZECH_DATA *Z, INT *mtx, INT *mtx_inv, INT n)
#else
#endif
{
	INT *M, i, j, a, N;
	INT idx_zero, idx_one;
	
	idx_zero = 0;
	idx_one = Z->Num_inv[0];
	N = 2 * n;
	M = (INT *) my_malloc(n * N * sizeof(INT), "z_mtx_inverse_GFq_num");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			M[i * N + j] = mtx[i * n + j];
			}
		}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == j)
				M[i * N + N + j] = idx_one;
			else
				M[i * N + N + j] = idx_zero;
			}
		}
#ifdef DEBUG_MTX_INVERSE
	printf("z_mtx_inverse_GFq_num(): before gauss:\n");
	z_mtx_print_nm(M, n, N);
#endif
	z_gauss_n_m_GFq_num(Z, FALSE /* f_special */, M, n, N);
#ifdef DEBUG_MTX_INVERSE
	printf("z_mtx_inverse_GFq_num(): after gauss:\n");
	z_mtx_print_nm(M, n, N);
#endif
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a = M[i * N + N + j];
			mtx_inv[i * n + j] = a;
			}
		}
	my_free(M);
	return OK;
}

#if TEXDOCU
INT z_gauss_n_m_GFp(INT p, INT f_special, INT *mtx, INT n, INT m)
#else
returns the rank.
f\_special TRUE: multiply (the row) with the inverse of the pivot to get the pivot to 1.
OUTPUT: coefficients absolutely $\le p/2$ (possibly negative) 
#endif
{
	INT i, j, k, jj, min_nm;
	INT a, a1, a_inv, b, c, piv, piv_inv;
	
	i = 0;
	min_nm = MIN(n, m);
	for (j = 0; j < min_nm; j++) {
	
		/* search for pivot element: */
		for (k = i; k < n; k++) {
			if (Asr(mtx[k * m + j], p) != 0) {
				/* pivot element found: */
				if (k != i) {
					for (jj = j; jj < m; jj++) {
						a = mtx[i * m + jj];
						mtx[i * m + jj] = mtx[k * m + jj];
						mtx[k * m + jj] = a;
						}
					}
				break;
				} /* if != 0 */
			} /* next k */
		
		if (k == n) /* no pivot found */
			continue; /* increase j, leave i constant */
		
		if (!f_special) {
			/* make pivot to 1: */
			a = mtx[i * m + j];
			a_inv = Inverse_mod(a, p);
			for (jj = j; jj < m; jj++) {
				mtx[i * m + jj] = Asr(a_inv * mtx[i * m + jj], p);
				}
			}
		else {
			piv = mtx[i * m + j];
			piv_inv = Inverse_mod(piv, p);
			}
		
		/* do the gaussian elimination: */
		for (k = i + 1; k < n; k++) {
			a1 = mtx[k * m + j];
			a1 = Asr(a1, p);
			if (a1 == 0)
				continue;
			mtx[k * m + j] = 0;
			for (jj = j + 1; jj < m; jj++) {
				a = mtx[i * m + jj];
				b = mtx[k * m + jj];
				if (f_special) {
					a = Asr(a * piv_inv, p);
					}
				c = Asr(b - a1 * a, p);
				mtx[k * m + jj] = c;
				}
			}
		i++;
		} /* next j */
	return i;
}

#if TEXDOCU
INT z_gauss_n_m_GFq_num(ZECH_DATA *Z, INT f_special, 
	INT *mtx, INT n, INT m)
#else
returns the rank.
f\_special TRUE: multiply (the row) with the inverse of the pivot to get the pivot to 1.
#endif
{
	INT i, j, k, jj, min_nm;
	INT a, a1, a_inv, b, c, piv, piv_inv;
	INT idx_zero, idx_one;
	
	idx_zero = 0;
	idx_one = Z->Num_inv[0];
	i = 0;
	min_nm = MIN(n, m);
	for (j = 0; j < min_nm; j++) {
	
		/* search for pivot element: */
		for (k = i; k < n; k++) {
			if (mtx[k * m + j] != idx_zero) {
				/* pivot element found: */
				if (k != i) {
					for (jj = j; jj < m; jj++) {
						a = mtx[i * m + jj];
						mtx[i * m + jj] = mtx[k * m + jj];
						mtx[k * m + jj] = a;
						}
					}
				break;
				} /* if != 0 */
			} /* next k */
		
		if (k == n) /* no pivot found */
			continue; /* increase j, leave i constant */
		
		if (!f_special) {
			/* make pivot to 1: */
			a = mtx[i * m + j];
			a_inv = z_inverse_num(Z, a);
			for (jj = j; jj < m; jj++) {
				mtx[i * m + jj] = z_mult_num(Z, a_inv, mtx[i * m + jj]);
				}
			}
		else {
			piv = mtx[i * m + j];
			piv_inv = z_inverse_num(Z, piv);
			}
		
		/* do the gaussian elimination: */
		for (k = i + 1; k < n; k++) {
			a1 = mtx[k * m + j];
			if (a1 == idx_zero)
				continue;
			mtx[k * m + j] = 0;
			for (jj = j + 1; jj < m; jj++) {
				a = mtx[i * m + jj];
				b = mtx[k * m + jj];
				if (f_special) {
					a = z_mult_num(Z, piv_inv, a);
					}
				c = z_mult_num(Z, a1, a);
				c = z_negate_num(Z, c);
				c = z_add_num(Z, b, c);
				/* c = b - a1 * a */
				mtx[k * m + jj] = c;
				}
			}
		i++;
		} /* next j */
	return i;
}

#if TEXDOCU
INT z_gauss_n_m_rectangular_GFp(INT p, INT f_special, INT *mtx, INT n, INT m, INT dim_m, 
	INT *base_cols)
#else
returns the rank.
f\_special TRUE: multiply (the row) with the inverse of the pivot to get the pivot to 1.
OUTPUT: coefficients absolutely $\le p/2$ (possibly negative) 
#endif
{
	INT i, j, k, jj;
	INT a, a1, a_inv, b, c, piv, piv_inv;
	
	i = 0;
	for (j = 0; j < m; j++) {
	
		/* search for pivot element: */
		for (k = i; k < n; k++) {
			if (Asr(mtx[k * dim_m + j], p) != 0) {
				/* pivot element found: */
				if (k != i) {
					for (jj = j; jj < m; jj++) {
						a = mtx[i * dim_m + jj];
						mtx[i * dim_m + jj] = mtx[k * dim_m + jj];
						mtx[k * dim_m + jj] = a;
						}
					}
				break;
				} /* if != 0 */
			} /* next k */
		
		if (k == n) /* no pivot found */
			continue; /* increase j, leave i constant */
		
		base_cols[i] = j;
		if (!f_special) {
			/* make pivot to 1: */
			a = mtx[i * dim_m + j];
			a_inv = Inverse_mod(a, p);
			for (jj = j; jj < m; jj++) {
				mtx[i * dim_m + jj] = Asr(a_inv * mtx[i * dim_m + jj], p);
				}
			}
		else {
			piv = mtx[i * dim_m + j];
			piv_inv = Inverse_mod(piv, p);
			}
		
		/* do the gaussian elimination: */
		for (k = i + 1; k < n; k++) {
			a1 = mtx[k * dim_m + j];
			a1 = Asr(a1, p);
			if (a1 == 0)
				continue;
			mtx[k * dim_m + j] = 0;
			for (jj = j + 1; jj < m; jj++) {
				a = mtx[i * dim_m + jj];
				b = mtx[k * dim_m + jj];
				if (f_special) {
					a = Asr(a * piv_inv, p);
					}
				c = Asr(b - a1 * a, p);
				mtx[k * dim_m + jj] = c;
				}
			}
		i++;
		} /* next j */
	return i;
}

#if TEXDOCU
INT z_gauss_n_m_rectangular_GFq_num(ZECH_DATA *Z, INT f_special, 
	INT *mtx, INT n, INT m, INT dim_m, INT *base_cols)
#else
returns the rank.
f\_special TRUE: multiply (the row) with the inverse of the pivot to get the pivot to 1.
#endif
{
	INT i, j, k, jj;
	INT a, a1, a_inv, b, c, piv, piv_inv;
	INT idx_zero, idx_one;
	
	idx_zero = 0;
	idx_one = Z->Num_inv[0];
	i = 0;
	for (j = 0; j < m; j++) {
	
		/* search for pivot element: */
		for (k = i; k < n; k++) {
			if (mtx[k * dim_m + j] != idx_zero) {
				/* pivot element found: */
				if (k != i) {
					for (jj = j; jj < dim_m; jj++) {
						a = mtx[i * dim_m + jj];
						mtx[i * dim_m + jj] = mtx[k * dim_m + jj];
						mtx[k * dim_m + jj] = a;
						}
					}
				break;
				} /* if != 0 */
			} /* next k */
		
		if (k == n) /* no pivot found */
			continue; /* increase j, leave i constant */
		
		base_cols[i] = j;
		if (!f_special) {
			/* make pivot to 1: */
			a = mtx[i * dim_m + j];
			a_inv = z_inverse_num(Z, a);
			for (jj = j; jj < m; jj++) {
				mtx[i * dim_m + jj] = z_mult_num(Z, a_inv, mtx[i * dim_m + jj]);
				}
			}
		else {
			piv = mtx[i * dim_m + j];
			piv_inv = z_inverse_num(Z, piv);
			}
		
		/* do the gaussian elimination: */
		for (k = i + 1; k < n; k++) {
			a1 = mtx[k * dim_m + j];
			if (a1 == idx_zero)
				continue;
			mtx[k * dim_m + j] = 0;
			for (jj = j + 1; jj < m; jj++) {
				a = mtx[i * dim_m + jj];
				b = mtx[k * dim_m + jj];
				if (f_special) {
					a = z_mult_num(Z, piv_inv, a);
					}
				c = z_mult_num(Z, a1, a);
				c = z_negate_num(Z, c);
				c = z_add_num(Z, b, c);
				/* c = b - a1 * a */
				mtx[k * dim_m + jj] = c;
				}
			}
		i++;
		} /* next j */
	return i;
}


#endif /* GFQ_TRUE */


