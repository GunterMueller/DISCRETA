/* tdo.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef GEO_TRUE

#include <DISCRETA/geo.h>
#include <DISCRETA/divs.h> // for STRING_OB

INT tdos_init_decomposition(TDO_SCHEME *tdos, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp)
{
	INT i, j, a, v, b, m, n;
	INT size;

	m = row_decomp->s_li();
	v = 0;
	for (i = 0; i < m; i++) {
		a = row_decomp->s_ii(i);
		v += a;
		}
	n = col_decomp->s_li();
	b = 0;
	for (i = 0; i < n; i++) {
		a = col_decomp->s_ii(i);
		b += a;
		}
	tdos->m = m;
	tdos->n = n;
	m++;
	n++;
	size = m * n * sizeof(INT);
	tdos->a = (INT *) my_malloc(size, "tdos_init_decomposition");
	if (tdos->a == NIL) {
		error("tdos_init_grid() no memory");
		return FALSE;
		}
	for (i = 0; i < n * m; i++) 
		tdos->a[i] = 0;
	
	for (i = 0; i < m - 1; i++) {
		a = row_decomp->s_ii(i);
		tdos->a[i * n + tdos->n] = a;
		}
	for (i = 0; i < n - 1; i++) {
		a = col_decomp->s_ii(i);
		tdos->a[tdos->m * n + i] = a;
		}
	return OK;
}

INT tdos_init_grid(TDO_SCHEME *tdos, 
	INT v, INT b, 
	INT nb_i_hbar, INT *i_hbar, 
	INT nb_i_vbar, INT *i_vbar)
{
	INT i, n, m, size, first, next, len;
	
	tdos->m = nb_i_hbar;
	tdos->n = nb_i_vbar;
	m = tdos->m + 1;
	n = tdos->n + 1;
	size = m * n * sizeof(INT);
	tdos->a = (INT *) my_malloc(size, "tdos_init_grid");
	if (tdos->a == NIL) {
		error("tdos_init_grid() no memory");
		return FALSE;
		}
	for (i = 0; i < n * m; i++) 
		tdos->a[i] = 0;
	
	for (i = 0; i < nb_i_hbar; i++) {
		first = i_hbar[i];
		if (i + 1 < nb_i_hbar)
			next = i_hbar[i + 1];
		else
			next = v;
		len = next - first;
		tdos->a[i * n + tdos->n] = len;
		}
	
	for (i = 0; i < nb_i_vbar; i++) {
		first = i_vbar[i];
		if (i + 1 < nb_i_vbar)
			next = i_vbar[i + 1];
		else
			next = b;
		len = next - first;
		tdos->a[tdos->m * n + i] = len;
		}
	
#if 0
	ii = 0;
	i_last = 0;
	for (i = 1; i <= v; i++) {
		if (i == v || i_hbar[i]) {
			len = i - i_last;
			i_last = i;
			tdos->a[ii * n + tdos->n] = len;
			ii++;
			}
		}
	jj = 0;
	j_last = 0;
	for (j = 1; j <= b; j++) {
		if (j == b || i_vbar[j]) {
			len = j - j_last;
			j_last = j;
			tdos->a[tdos->m * n + jj] = len;
			jj++;
			}
		}
#endif
	return TRUE;
}

INT tdos_free(TDO_SCHEME *tdos)
{
	if (tdos->a) {
		my_free(tdos->a);
		tdos->a = NIL;
		}
	return TRUE;
}

void tdos_nil(TDO_SCHEME *tdos)
{
	tdos->a = NIL;
}

INT tdos_m(TDO_SCHEME *tdos)
{
	return (tdos->m);
}

INT tdos_n(TDO_SCHEME *tdos)
{
	return (tdos->n);
}

INT tdos_Vi(TDO_SCHEME *tdos, INT i)
{
	INT n;
	
	n = tdos->n + 1;
	return (tdos->a[i * n + tdos->n]);
}

INT tdos_Bj(TDO_SCHEME *tdos, INT j)
{
	INT n;
	
	n = tdos->n + 1;
	return (tdos->a[tdos->m * n + j]);
}

INT tdos_ij(TDO_SCHEME *tdos, INT i, INT j)
{
	INT n;
	
	n = tdos->n + 1;
	return (tdos->a[i * n + j]);
}

void tdos_copy(TDO_SCHEME *t1, 
	TDO_SCHEME *t2)
{
	t2->m = t1->m;
	t2->n = t1->n;
	t2->a = t1->a;
}

void tdos_print(TDO_SCHEME *tdos)
{
	INT m, n, i, j;
	
	m = tdos->m + 1;
	n = tdos->n + 1;
	printf("      ");
	for (j = 0; j < tdos->n; j++)
		printf("%3ld ",  
			tdos->a[tdos->m * n + j]);
	printf("\n");
	for (i = 0; i < tdos->m; i++) {
		printf("%3ld | ",  
			tdos->a[i * n + n - 1]);
		for (j = 0; j < tdos->n; j++) {
			printf("%3ld",  
				tdos->a[i * n + j]);
			if (j < tdos->n - 1)
				printf(" ");
			}
		printf("\n");
		}
}

INT tdos_cmp(TDO_SCHEME *t1, 
	TDO_SCHEME *t2)
{
	INT m, n, i, l;
	
	if (t1->m < t2->m)
		return -1;
	if (t1->m > t2->m)
		return 1;
	if (t1->n < t2->n)
		return -1;
	if (t1->n > t2->n)
		return 1;
	m = t1->m + 1;
		/* inclusive Rand */
	n = t1->n + 1;
	l = m * n - 1;
		/* unten rechts steht eine 0, 
		 * weglassen */
	for (i = 0; i < l; i++) {
		if (t1->a[i] < t2->a[i])
			return -1;
		if (t1->a[i] > t2->a[i])
			return +1;
		}
	return 0;
}

TDOSS *tdoss_init_decomposition(VECTOR_OP row_decomp, VECTOR_OP col_decomp)
{
	TDO_SCHEME tdos;
	TDOSS *tdoss;

	tdos_init_decomposition(&tdos, row_decomp, col_decomp);
	tdos2tdoss(&tdos, NIL, NIL, &tdoss);
	tdos_free(&tdos);
	return tdoss;
}

void tdoss_free(TDOSS *t)
{
	SHORT *ddp_mult, *ddb_mult;
	
	ddp_mult = tdoss_get_ddp_mult(t);
	if (ddp_mult)
		my_free(ddp_mult);
	ddb_mult = tdoss_get_ddp_mult(t);
	if (ddb_mult)
		my_free(ddb_mult);
	my_free(t);
}

void tdoss_print(TDOSS *t)
{
	INT m, n, m1, n1, i, j;
	SHORT *ddp_mult, *ddb_mult;
	INT V, B, N;

	ddp_mult = tdoss_get_ddp_mult(t);
	ddb_mult = tdoss_get_ddb_mult(t);
	m = tdoss_m(t);
	n = tdoss_n(t);
	m1 = m + 1;
	n1 = n + 1;
	V = tdoss_V(t);
	B = tdoss_B(t);
	printf("V = %ld B = %ld\n", V, B);
	printf("      ");
	for (j = 0; j < n; j++)
		printf("%3ld ", tdoss_Bj(t, j));
	printf("\n");
	for (i = 0; i < m; i++) {
		printf("%3ld | ", tdoss_Vi(t, i));
		for (j = 0; j < n; j++) {
			printf("%3ld", tdoss_ij(t, i, j));
			if (j < n - 1)
				printf(" ");
			}
		printf("\n");
		}
	if (ddp_mult) {
		printf("mult (points):\n");
		N = (INT) ddp_mult[0];
		for (i = 0; i < N; i++) {
			printf("%ld ", 
			(INT) ddp_mult[1 + i]);
			}
		printf("\n");
		}
	if (ddb_mult) {
		printf("mult (blocks):\n");
		N = (INT) ddb_mult[0];
		for (i = 0; i < N; i++) {
			printf("%ld ", 
			(INT) ddb_mult[1 + i]);
			}
		printf("\n");
		}
	 fflush(stdout);
}

INT tdoss_cmp(TDOSS *t1, TDOSS *t2)
{
	INT m1, m2, n1, n2, i, l;
	SHORT *ddp1_mult, *ddp2_mult;
	SHORT *ddb1_mult, *ddb2_mult;
	INT N1, N2;
	
	m1 = tdoss_m(t1);
	m2 = tdoss_m(t2);
	if (m1 < m2)
		return -1;
	if (m1 > m2)
		return 1;
	n1 = tdoss_n(t1);
	n2 = tdoss_n(t2);
	if (n1 < n2)
		return -1;
	if (n1 > n2)
		return 1;
	l = (m1 + 1) * (n1 + 1) - 1;
	/* unten rechts steht eine 0, weglassen */
	t1 += TDOSS_OFF;
	t2 += TDOSS_OFF;
	t1 += 2;
	t2 += 2;
	for (i = 0; i < l; i++) {
		if (t1[i] < t2[i])
			return -1;
		if (t1[i] > t2[i])
			return +1;
		}
	t1 -= TDOSS_OFF;
	t2 -= TDOSS_OFF;
	t1 -= 2;
	t2 -= 2;
	ddp1_mult = tdoss_get_ddp_mult(t1);
	ddp2_mult = tdoss_get_ddp_mult(t2);
	/* Punktableitungen vorhanden ? */
	if (ddp1_mult && !ddp2_mult) {
		return error("tdoss_cmp() ddp1_mult && !ddp2_mult");
		}
	if (!ddp1_mult && ddp2_mult) {
		return error("tdoss_cmp() !ddp1_mult && ddp2_mult");
		}
	if (ddp1_mult) {
		N1 = (INT) ddp1_mult[0];
		N2 = (INT) ddp2_mult[0];
		if (N1 < N2)
			return -1;
		if (N1 > N2)
			return 1;
		for (i = 0; i < N1; i++) {
			if (ddp1_mult[1 + i] < ddp2_mult[1 + i])
				return -1;
			if (ddp1_mult[1 + i] > ddp2_mult[1 + i])
				return +1;
			}
		}
	ddb1_mult = tdoss_get_ddb_mult(t1);
	ddb2_mult = tdoss_get_ddb_mult(t2);
	/* Blockableitungen vorhanden ? */
	if (ddb1_mult && !ddb2_mult) {
		return error("tdoss_cmp() ddb1_mult && !ddb2_mult");
		}
	if (!ddb1_mult && ddb2_mult) {
		return error("tdoss_cmp() !ddb1_mult && ddb2_mult");
		}
	if (ddb1_mult) {
		N1 = (INT) ddb1_mult[0];
		N2 = (INT) ddb2_mult[0];
		if (N1 < N2)
			return -1;
		if (N1 > N2)
			return 1;
		for (i = 0; i < N1; i++) {
			if (ddb1_mult[1 + i] < ddb2_mult[1 + i])
				return -1;
			if (ddb1_mult[1 + i] > ddb2_mult[1 + i])
				return +1;
			}
		}
	return 0;
}

INT tdos2tdoss(TDO_SCHEME *tdos, 
	SHORT *ddp_mult, SHORT *ddb_mult, TDOSS **t)
{
	INT m, n, m1, n1, i, l, size;
	TDOSS *T;
	
	m = tdos->m;
	n = tdos->n;
	m1 = m + 1;
	n1 = n + 1;
	size = (TDOSS_OFF + 2 + m1 * n1) * sizeof(TDOSS);
	T = (TDOSS *) my_malloc(size, "tdos2tdoss");
	if (T == NIL) {
		return error("tdos2tdoss() no memory");
		}
	*t = T;
	tdoss_set_ddp_mult(T, ddp_mult);
	tdoss_set_ddb_mult(T, ddb_mult);
	T[TDOSS_OFF + 0] = m;
	T[TDOSS_OFF + 1] = n;
	l = m1 * n1;
	T += TDOSS_OFF;
	T += 2;
	for (i = 0; i < l; i++)
		T[i] = (TDOSS) tdos->a[i];
	return TRUE;
}

SHORT *tdoss_get_ddp_mult(TDOSS *t)
{
	SHORT **t1 = (SHORT **)t;

	return t1[0];
}

SHORT *tdoss_get_ddb_mult(TDOSS *t)
{
	SHORT **t1 = (SHORT **)t;

	return t1[1];
}

void tdoss_set_ddp_mult(TDOSS *t, SHORT *ddp_mult)
{
	SHORT **t1 = (SHORT **)t;

	t1[0] = ddp_mult;
}

void tdoss_set_ddb_mult(TDOSS *t, SHORT *ddb_mult)
{
	SHORT **t1 = (SHORT **)t;

	t1[1] = ddb_mult;
}

INT tdoss_m(TDOSS *t)
{
	return (INT) t[TDOSS_OFF + 0];
}

INT tdoss_n(TDOSS *t)
{
	return (INT) t[TDOSS_OFF + 1];
}

INT tdoss_Vi(TDOSS *t, INT i)
{
	INT m = tdoss_m(t);
	INT n = tdoss_n(t);
	INT n1 = n + 1;
	
	return ((INT) t[TDOSS_OFF + 2 + i * n1 + n]);
}

INT tdoss_Bj(TDOSS *t, INT j)
{
	INT m = tdoss_m(t);
	INT n = tdoss_n(t);
	INT n1 = n + 1;
	
	return ((INT) t[TDOSS_OFF + 2 + m * n1 + j]);
}

INT tdoss_V(TDOSS *t)
{
	INT i, m, V;
	
	V = 0;
	m = tdoss_m(t);
	for (i = 0; i < m; i++) {
		V += tdoss_Vi(t, i);
		}
	return V;
}

INT tdoss_B(TDOSS *t)
{
	INT j, n, B;
	
	B = 0;
	n = tdoss_n(t);
	for (j = 0; j < n; j++) {
		B += tdoss_Bj(t, j);
		}
	return B;
}

INT tdoss_ij(TDOSS *t, INT i, INT j)
{
	INT m = tdoss_m(t);
	INT n = tdoss_n(t);
	INT n1 = n + 1;
	
	return ((INT) t[TDOSS_OFF + 2 + i * n1 + j]);
}

INT tdoss_2_first(TDOSS *t, VECTOR_OP row_first, VECTOR_OP block_first)
{
	INT f, l, i, m, n;
	
	m = tdoss_m(t);
	row_first->m_il(m);
	f = 0;
	for (i = 0; i < m; i++) {
		row_first->m_ii(i, f);
		l = tdoss_Vi(t, i);
		f += l;
		}
	n = tdoss_n(t);
	block_first->m_il(n);
	f = 0;
	for (i = 0; i < n; i++) {
		block_first->m_ii(i, f);
		l = tdoss_Bj(t, i);
		f += l;
		}
	return OK;
}


// formerly tdo_scheme.C:

INT tdos_init_ntdo(TDO_SCHEME *tdos, NTDO *tdo, 
	NTDO_GRID *G0, NTDO_GRID *G1, INT f_derived)
{
	NTDO_GRID *Gpoints, *Gblocks;
	INT size;
	INT m, n, i, j, i1, first;
	
	tdos_free(tdos);
	if (G0->f_points) {
		Gpoints = G0;
		Gblocks = G1;
		}
	else {
		Gpoints = G1;
		Gblocks = G0;
		}
	tdos->m = Gpoints->G_max;
	tdos->n = Gblocks->G_max;
	m = tdos->m + 1;
	n = tdos->n + 1;
	size = m * n * sizeof(INT);
	tdos->a = (INT *) my_malloc(size, "tdos_init_ntdo");
	if (tdos->a == NIL) {
		return error("tdos_int(): no memory");
		}
	for (i = 0; i < tdos->m; i++) {
		first = Gpoints->first[i];
		i1 = Gpoints->type_idx[first];
		for (j = 0; j < tdos->n; j++) {
			tdos->a[i * n + j] = Gpoints->type[i1 * Gpoints->max_size + j];
			}
		}
	for (i = 0; i < tdos->m; i++)
		tdos->a[i * n + tdos->n] = Gpoints->len[i];
	for (j = 0; j < tdos->n; j++)
		tdos->a[tdos->m * n + j] = Gblocks->len[j];
	tdos->a[m * n - 1] = 0;
	return OK;
}

#undef DEBUG_STR_ALLOC

static void str_alloc(BYTE ***str, INT num, INT len)
{
	INT i;
	BYTE **p;
	
	p = (BYTE **) my_malloc(num * sizeof(BYTE *), "tdo.C str_alloc");
#ifdef DEBUG_STR_ALLOC
	printf("p allocated\n");
	fflush(stdout);
#endif
	for (i = 0; i < num; i++) {
		p[i] = (BYTE *) my_malloc(len, "tdo.C str_alloc p[i]");
#ifdef DEBUG_STR_ALLOC
		printf("p[%ld] allocated\n", i);
		fflush(stdout);
#endif
		}
	*str = p;
}

static void str_free(BYTE **str, INT num)
{
	INT i;
	
	for (i = 0; i < num; i++) {
#ifdef DEBUG_STR_ALLOC
		printf("my_free(%ld)\n", i);
		fflush(stdout);
#endif
		my_free(str[i]);
		}
	my_free(str);
}

INT tdos_print_theX_short(TDO_SCHEME *tdos, INT nb_X, SHORT *theX)
{
	INT i, j, first, k, a, v, b;
	INT vbar[NTDO_MAX_N + 1];
	INT hbar[NTDO_MAX_N + 1];
	INT M[NTDO_MAX_N][NTDO_MAX_N];
	BYTE **str;

	for (i = 0; i <= NTDO_MAX_N; i++) {
		vbar[i] = FALSE;
		hbar[i] = FALSE;
		}
	first = 0;
	for (i = 0; i < tdos_m(tdos); i++) {
		hbar[first] = TRUE;
		first += tdos_Vi(tdos, i);
		}
	v = first;
	if (v >= NTDO_MAX_N)
		return error("tdos_print_theX_short() v >= NTDO_MAX_N");
	hbar[first] = TRUE;
	first = 0;
	for (i = 0; i < tdos_n(tdos); i++) {
		vbar[first] = TRUE;
		first += tdos_Bj(tdos, i);
		}
	b = first;
	if (b >= NTDO_MAX_N)
		return error("tdos_print_theX_short() b >= NTDO_MAX_N");
	vbar[first] = TRUE;
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			M[i][j] = 0;
			}
		}
	for (k = 0; k < nb_X; k++) {
		a = theX[k];
		j = a % b;
		a -= j;
		i = a / b;
		M[i][j] = 1;
		}
	grid_alloc(&str, 3 * v, 3 * b);
	grid_prepare(str, v, b, NTDO_MAX_N, (INT *) M, vbar, hbar);
	grid_print(v, hbar, str);
	grid_free(str, 3 * v);
	return OK;
}

INT tdos_print(TDO_SCHEME *tdos, NTDO *tdo)
{
	NTDO_GRID *Gblocks, *Gpoints;
	INT i, j, ii, jj, first, k;
	INT vbar[NTDO_MAX_N + 1];
	INT hbar[NTDO_MAX_N + 1];
	INT M[NTDO_MAX_N][NTDO_MAX_N];
	BYTE **str;

	if (tdo->v >= NTDO_MAX_N)
		return error("tdos_print() tdo->v >= NTDO_MAX_N");
	if (tdo->b >= NTDO_MAX_N)
		return error("tdos_print() tdo->b >= NTDO_MAX_N");
	if (tdo->G_next->f_points) {
		Gpoints = tdo->G_next;
		Gblocks = tdo->G;
		}
	else {
		Gpoints = tdo->G;
		Gblocks = tdo->G_next;
		}
	for (i = 0; i <= NTDO_MAX_N; i++) {
		vbar[i] = FALSE;
		hbar[i] = FALSE;
		}
	for (i = 0; i < Gpoints->G_max; i++) {
		first = Gpoints->first[i];
		hbar[first] = TRUE;
		}
	hbar[tdo->v] = TRUE;
	for (i = 0; i < Gblocks->G_max; i++) {
		first = Gblocks->first[i];
		vbar[first] = TRUE;
		}
	vbar[tdo->b] = TRUE;

	for (i = 0; i < tdo->v; i++) {
		for (j = 0; j < tdo->b; j++) {
			M[i][j] = 0;
			}
		}
	for (k = 0; k < tdo->nb_X; k++) {
		i = tdo->theXi[k];
		j = tdo->theXj[k];
		ii = tdo->p.a[i];
		jj = tdo->q.a[j];
		M[ii][jj] = 1;
		}

	grid_alloc(&str, 3 * tdo->v, 3 * tdo->b);
	grid_prepare(str, tdo->v, tdo->b, NTDO_MAX_N, (INT *) M, vbar, hbar);
	grid_print(tdo->v, hbar, str);
	grid_free(str, 3 * tdo->v);
	// printf("tdos_print() finished()\n");
	return OK;
}

INT TD_char_inv_print(FILE *fp, INT f_tex, 
	INT nrow, INT ncol, INT nb_X, INT *X, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	VECTOR_OP row_char_first, VECTOR_OP col_char_first, 
	VECTOR_OP row_inv_first, VECTOR_OP col_inv_first)
{
	INT i, j, ii, jj, first, k, a;
	INT *vbar, *hbar, *M;
	// INT vbar[NTDO_MAX_N + 1];
	// INT hbar[NTDO_MAX_N + 1];
	// INT M[NTDO_MAX_N][NTDO_MAX_N];
	BYTE **str;
	
	vbar = (INT *) my_malloc((ncol + 1) * sizeof(INT), "TD_char_inv_print vbar");
	hbar = (INT *) my_malloc((nrow + 1) * sizeof(INT), "TD_char_inv_print hbar");
	M = (INT *) my_malloc(nrow * ncol * sizeof(INT), "TD_char_inv_print M");
	if (vbar == NIL || hbar == NIL || M == NIL)
		return error("TD_char_inv_print() no memory");
	for (i = 0; i <= ncol; i++)
		vbar[i] = 0;
	for (i = 0; i <= nrow; i++)
		hbar[i] = 0;
	for (i = 0; i < nrow; i++)
		for (j = 0; j < ncol; j++) 
			M[i * ncol + j] = 0;
#if 0
	if (nrow >= NTDO_MAX_N)
		return error("TD_char_inv_print() nrow >= NTDO_MAX_N");
	if (ncol >= NTDO_MAX_N)
		return error("TD_char_inv_print() ncol >= NTDO_MAX_N");
	for (i = 0; i <= NTDO_MAX_N; i++) {
		vbar[i] = 0;
		hbar[i] = 0;
		}
#endif
	for (i = 0; i < row_inv_first->s_li(); i++) {
		first = row_inv_first->s_ii(i);
		hbar[first] = 1;
		}
	hbar[nrow] = 1;
	for (i = 0; i < row_char_first->s_li(); i++) {
		first = row_char_first->s_ii(i);
		hbar[first] = 2;
		}
	hbar[nrow] = 2;
	for (i = 0; i < col_inv_first->s_li(); i++) {
		first = col_inv_first->s_ii(i);
		vbar[first] = 1;
		}
	vbar[ncol] = 1;
	for (i = 0; i < col_char_first->s_li(); i++) {
		first = col_char_first->s_ii(i);
		vbar[first] = 2;
		}
	vbar[ncol] = 2;
	hbar[0] = 3;
	vbar[0] = 3;
	hbar[nrow] = 3;
	vbar[ncol] = 3;
	
#if 0
	for (i = 0; i < nrow; i++) {
		for (j = 0; j < ncol; j++) {
			M[i * ncol + j] = 0;
			}
		}
#endif
	for (k = 0; k < nb_X; k++) {
		a = X[k];
		i = a / ncol;
		j = a % ncol;
		if (p)
			ii = p->s_ii(i) - 1;
		else
			ii = i;
		if (q)
			jj = q->s_ii(j) - 1;
		else
			jj = j;
		M[ii * ncol + jj] = 1;
		}

	grid_alloc(&str, 3 * nrow, 3 * ncol);
	grid_prepare2(str, nrow, ncol, ncol, M, vbar, hbar);
	grid_fprint(fp, f_tex, nrow, hbar, str);
	grid_free(str, 3 * nrow);
	// printf("TD_char_inv_print() finished()\n");
	my_free(M);
	my_free(vbar);
	my_free(hbar);
	return OK;
}
	
INT TD_print(FILE *fp, INT nrow, INT ncol, INT nb_X, INT *X, INT f_transposed, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	VECTOR_OP P_first, VECTOR_OP Q_first)
// p and q may be ommitted
{
	INT i, j, ii, jj, first, k, a;
	INT *vbar, *hbar, *M;
	// INT vbar[NTDO_MAX_N + 1];
	// INT hbar[NTDO_MAX_N + 1];
	// INT M[NTDO_MAX_N][NTDO_MAX_N];
	BYTE **str;
	PERMUTATION_OP pp, qq;
	VECTOR_OP PP_first, QQ_first;

	vbar = (INT *) my_malloc((ncol + 1) * sizeof(INT), "TD_print vbar");
	hbar = (INT *) my_malloc((nrow + 1) * sizeof(INT), "TD_print hbar");
	M = (INT *) my_malloc(nrow * ncol * sizeof(INT), "TD_print M");
	if (vbar == NIL || hbar == NIL || M == NIL)
		return error("TD_print() no memory");
	for (i = 0; i <= ncol; i++)
		vbar[i] = 0;
	for (i = 0; i <= nrow; i++)
		hbar[i] = 0;
	for (i = 0; i < nrow; i++)
		for (j = 0; j < ncol; j++) 
			M[i * ncol + j] = 0;
#if 0
	if (nrow >= NTDO_MAX_N)
		return error("TD_print() nrow >= NTDO_MAX_N");
	if (ncol >= NTDO_MAX_N)
		return error("TD_print() ncol >= NTDO_MAX_N");
	for (i = 0; i <= NTDO_MAX_N; i++) {
		vbar[i] = FALSE;
		hbar[i] = FALSE;
		}
	for (i = 0; i < nrow; i++) {
		for (j = 0; j < ncol; j++) {
			M[i][j] = 0;
			}
		}
#endif
	if (f_transposed) {
		PP_first = P_first;
		QQ_first = Q_first;
		pp = p;
		qq = q;
		}
	else {
		PP_first = Q_first;
		QQ_first = P_first;
		pp = q;
		qq = p;
		}
	for (i = 0; i < PP_first->s_li(); i++) {
		first = PP_first->s_ii(i);
		hbar[first] = TRUE;
		}
	hbar[nrow] = TRUE;
	
	for (i = 0; i < QQ_first->s_li(); i++) {
		first = QQ_first->s_ii(i);
		vbar[first] = TRUE;
		}
	vbar[ncol] = TRUE;

	for (k = 0; k < nb_X; k++) {
		a = X[k];
		i = a / ncol;
		j = a % ncol;
		if (pp)
			ii = pp->s_ii(i) - 1;
		else
			ii = i;
		if (qq)
			jj = qq->s_ii(j) - 1;
		else
			jj = j;
		M[ii * ncol + jj] = 1;
		}

	grid_alloc(&str, 3 * nrow, 3 * ncol);
	grid_prepare(str, nrow, ncol, ncol, M, vbar, hbar);
	grid_fprint(fp, FALSE, nrow, hbar, str);
	grid_free(str, 3 * nrow);
	// printf("TD_print() finished()\n");
	my_free(M);
	my_free(vbar);
	my_free(hbar);
	return OK;
}

void grid_print(INT v, INT *hbar, BYTE **str)
{
	INT i, ii;

	ii = 0;
	for (i = 0; i <= v; i++) {
		if (hbar[i]) {
			printf("%s\n", str[ii]);
			ii++;
			}
		printf("%s\n", str[ii]);
		ii++;
		}
}

void grid_fprint(FILE *fp, INT f_tex, INT v, INT *hbar, BYTE **str)
{
	INT i, ii;

	ii = 0;
	for (i = 0; i <= v; i++) {
		if (hbar[i]) {
			if (f_tex) {
				fprintf(fp, "%s\\\\[-4pt]\n", str[ii]);
				}
			else {
				fprintf(fp, "%s\n", str[ii]);
				}
			ii++;
			}
		if (f_tex) {
			fprintf(fp, "%s\\\\[-4pt]\n", str[ii]);
			}
		else {
			fprintf(fp, "%s\n", str[ii]);
			}
		ii++;
		}
}

void grid_alloc(BYTE ***str, INT num, INT len)
{
	INT i;

	str_alloc(str, num, len);

	for (i = 0; i < num; i++) {
		(*str)[i][0] = 0;
		}
}

void grid_free(BYTE **str, INT num)
{
	str_free(str, num);
}

void grid_prepare(BYTE **str, INT v, INT b, INT dim_M, INT *M, INT *vbar, INT *hbar)
{
	INT i, j, ii;
	
	for (j = 0; j <= b; j++) {
	
		if (vbar[j]) {

			/* vertical border: */
			ii = 0;
			for (i = 0; i <= v; i++) {
				if (hbar[i]) {
					if (i == 0 || i == v)
						strcat(str[ii], "#");
					else {
						if (j == 0 || j == b)
							strcat(str[ii], "#");
						else
							strcat(str[ii], "+");
						}
					ii++;
					}
				if (i < v) {
					if (j == 0 || j == b)
						strcat(str[ii], "#");
					else
						strcat(str[ii], "|");
					ii++;
					}
				} /* next i */

			} /* if vbar[j] */
		
		if (j < b) {
			ii = 0;
			for (i = 0; i <= v; i++) {

				if (hbar[i]) {
					if (i == 0 || i == v)
						strcat(str[ii], "#");
					else
						strcat(str[ii], "-");
					ii++;
					}
				
				if (i < v) {
					/* we are at incidence (i,j) : */
					if (M[i * dim_M + j])
						strcat(str[ii], "X");
					else
						strcat(str[ii], ".");
					ii++;
					}

				} /* next i */
			}

		} /* next j */

}

void grid_prepare2(BYTE **str, INT v, INT b, INT dim_M, INT *M, INT *vbar, INT *hbar)
{
	INT i, j, ii;
#if 0
	BYTE *v_sign[] = { "", "|", "#", "#" };
	BYTE *h_sign[] = { "", "-", "#", "#" };
	BYTE *x_sign[] = { "", "+", "#", "#" };
#endif
	BYTE *v_sign[] = { "", "|", "|", "|" };
	BYTE *h_sign[] = { "", "-", "-", "-" };
	BYTE *x_sign[][4] = {
		{ "", "", "", "" }, 
		{ "", "+", "+", "A" }, 
		{ "", "+", "+", "+" }, 
		{ "", "A", "+", "+" }, 
		};
	
	for (j = 0; j <= b; j++) {
	
		if (vbar[j]) {

			/* vertical border: */
			ii = 0;
			for (i = 0; i <= v; i++) {
				if (hbar[i]) {
					strcat(str[ii], x_sign[hbar[i]][vbar[j]] );
					ii++;
					}
				if (i < v) {
					strcat(str[ii], v_sign[vbar[j]]);
					ii++;
					}
				} /* next i */

			} /* if vbar[j] */
		
		if (j < b) {
			ii = 0;
			for (i = 0; i <= v; i++) {

				if (hbar[i]) {
					strcat(str[ii], h_sign[hbar[i]]);
					ii++;
					}
				
				if (i < v) {
					/* we are at incidence (i,j) : */
					if (M[i * dim_M + j])
						strcat(str[ii], "X");
					else
						strcat(str[ii], ".");
					ii++;
					}

				} /* next i */
			}

		} /* next j */

}


/*
 *
 */


INT tdoss_print_theX_short(TDOSS *tdoss, INT nb_X, SHORT *theX)
{
	INT i, j, first, k, a, v, b;
	INT vbar[NTDO_MAX_N + 1];
	INT hbar[NTDO_MAX_N + 1];
	INT M[NTDO_MAX_N][NTDO_MAX_N];
	BYTE **str;

	for (i = 0; i <= NTDO_MAX_N; i++) {
		vbar[i] = FALSE;
		hbar[i] = FALSE;
		}
	first = 0;
	for (i = 0; i < tdoss_m(tdoss); i++) {
		hbar[first] = TRUE;
		first += tdoss_Vi(tdoss, i);
		}
	v = first;
	if (v >= NTDO_MAX_N)
		return error("tdoss_print_theX_short() v >= NTDO_MAX_N");
	hbar[first] = TRUE;
	first = 0;
	for (i = 0; i < tdoss_n(tdoss); i++) {
		vbar[first] = TRUE;
		first += tdoss_Bj(tdoss, i);
		}
	b = first;
	if (b >= NTDO_MAX_N)
		return error("tdoss_print_theX_short() b >= NTDO_MAX_N");
	vbar[first] = TRUE;
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			M[i][j] = 0;
			}
		}
	for (k = 0; k < nb_X; k++) {
		a = theX[k];
		j = a % b;
		a -= j;
		i = a / b;
		M[i][j] = 1;
		}
	grid_alloc(&str, 3 * v, 3 * b);
	grid_prepare(str, v, b, NTDO_MAX_N, (INT *) M, vbar, hbar);
	grid_print(v, hbar, str);
	grid_free(str, 3 * v);
	return OK;
}

INT incma_latex_picture(FILE *fp, 
	INT width, INT width_10, 
	INT f_outline_thin, BYTE *unit_length, 
	BYTE *thick_lines, BYTE *thin_lines, BYTE *geo_line_width, 
	INT v, INT b, 
	INT V, INT B, INT *Vi, INT *Bj, 
	INT *R, INT *X, INT dim_X, 
	INT f_labelling_points, BYTE **point_labels, 
	INT f_labelling_blocks, BYTE **block_labels)
/* width for one box in 0.1mm 
 * width_10 is 1 10th of width
 * example: width = 40, width_10 = 4 */
{
	INT w, h, w1, h1;
	INT i, j, k, r, a;
	INT x0, y0, x1, y1;
	INT X0, Y0, X1, Y1;
	INT width_8, width_5;
	BYTE *tdo_line_width = thick_lines /* "0.7mm" */;
	BYTE *line_width = thin_lines /* "0.15mm" */;
	/* BYTE *geo_line_width = "0.25mm"; */
	
	width_8 = width - 2 * width_10;
	width_5 = width >> 1;
	fprintf(fp, "\\unitlength%s\n", unit_length);
	w = b * width;
	h = v * width;
	w1 = w;
	h1 = h;
	if (f_labelling_points)
		w1 += 2 * width;
	if (f_labelling_blocks)
		h1 += 2 * width;
	fprintf(fp, "\\begin{picture}(%ld,%ld)\n", w1, h1);

	/* the grid: */
	fprintf(fp, "\\linethickness{%s}\n", tdo_line_width);
	k = 0;
	for (i = -1; i < B; i++) {
		if (i >= 0) {
			a = Bj[i];
			k += a;
			}
		if (f_outline_thin) {
			if (i == -1 || i == B - 1)
				continue;
			}
		fprintf(fp, "\\put(%ld,0){\\line(0,1){%ld}}\n", 
			k * width, h);
		}
	if (k != b)
		return error("latex_picture(): k != b");
	k = 0;
	for (i = -1; i < V; i++) {
		if (i >= 0) {
			a = Vi[i];
			k += a;
			}
		if (f_outline_thin) {
			if (i == -1 || i == V - 1)
				continue;
			}
		fprintf(fp, "\\put(0,%ld){\\line(1,0){%ld}}\n", 
			h - k * width, w);
		}
	if (k != v)
		return error("latex_picture(): k != v");
	if (f_labelling_points) {
		for (i = 0; i < v; i++) {
			fprintf(fp, "\\put(0,%ld){\\makebox(0,0)[r]{%s$\\,$}}\n", 
				h - i * width - width_5, point_labels[i]);
			}
		}
	if (f_labelling_blocks) {
		for (i = 0; i < b; i++) {
			fprintf(fp, "\\put(%ld,%ld){\\makebox(0,0)[b]{%s}}\n", 
				i * width + width_5, h + width_5, block_labels[i]);
			}
		}

	fprintf(fp, "\\linethickness{%s}\n", line_width);
	fprintf(fp, "\\multiput(0,0)(%ld,0){%ld}{\\line(0,1){%ld}}\n", 
		width, b + 1, h);
	fprintf(fp, "\\multiput(0,0)(0,%ld){%ld}{\\line(1,0){%ld}}\n", 
		width, v + 1, w);

	/* the geometry: */
	fprintf(fp, "\\linethickness{%s}\n", geo_line_width);
	for (i = 0; i < v; i++) {
		y0 = h - i * width;
		y1 = h - (i + 1) * width;
		Y0 = y0 - width_10;
		Y1 = y1 + width_10;
		for (r = 0; r < R[i]; r++) {
			j = X[i * dim_X + r];
			// printf("%ld ", j);
			x0 = j * width;
			x1 = (j + 1) * width;
			X0 = x0 + width_10;
			X1 = x1 - width_10;
			/* hor. lines: */
			fprintf(fp, "\\put(%ld,%ld){\\line(1,0){%ld}}\n", 
				X0, Y0, width_8);
			fprintf(fp, "\\put(%ld,%ld){\\line(1,0){%ld}}\n", 
				X0, Y1, width_8);

			/* vert. lines: */
			fprintf(fp, "\\put(%ld,%ld){\\line(0,1){%ld}}\n", 
				X0, Y1, width_8);
			fprintf(fp, "\\put(%ld,%ld){\\line(0,1){%ld}}\n", 
				X1, Y1, width_8);

			}
		// printf("\n");
		}

	fprintf(fp, "\\end{picture}\n");
	return OK;
}

INT incma_latex(FILE *fp, 
	INT v, INT b, 
	INT V, INT B, INT *Vi, INT *Bj, 
	INT *R, INT *X, INT dim_X)
{
	incma_latex_picture(fp, 
		40 /* width */, 
		10 /* width_10 */,  
		FALSE /* f_outline_thin */, 
		"0.065mm" /* unit_length */, 
		"0.5mm" /* thick_lines */ , 
		"0.15mm" /* thin_lines */ , 
		"0.25mm" /* geo_line_width */ , 
		v, b, V, B, Vi, Bj, R, X, dim_X, 
		FALSE /* f_labelling_points */, NIL, 
		FALSE /* f_labelling_blocks */, NIL);
	
	return OK;
}

INT convert_number_string(INT a, BYTE *str)
{
	if (a == 10)
		a = 0;
	if (a < 10)
		sprintf(str, "%ld", a);
	else {
		a -= 11;
		if (a < 26) {
			str[0] = 'a' + a;
			str[1] = 0;
			}
		else {
			a -= 26;
			if (a < 26) {
				str[0] = 'A' + a;
				str[1] = 0;
				}
			else {
				a -= 26;
				if (a < 26) {
					sprintf(str, "${\\cal %c}$", (BYTE)('a' + a));
					}
				else {
					a -= 26;
					sprintf(str, "${\\cal %c}$", (BYTE)('A' + a));
					}
				}
			}
		}
	return OK;
}

INT incma_latex_integer_labels(FILE *fp, 
	INT v, INT b, 
	INT V, INT B, INT *Vi, INT *Bj, 
	INT *R, INT *X, INT dim_X, 
	VECTOR_OP point_labels, VECTOR_OP block_labels, INT offset)
{
	VECTOR_OB pl, bl;
	BYTE **pp, **qq;
	BYTE str[1024];
	STRING_OB s;
	STRING_OP ps;
	INT i, j, a;
	
	if (point_labels->s_li() != v) 
		return error("incma_latex_integer_labels() point_labels->s_li() != v");
	if (block_labels->s_li() != b) 
		return error("incma_latex_integer_labels() block_labels->s_li() != b");
	pp = (BYTE **) my_malloc(sizeof(BYTE *) * v, "incma_latex_integer_labels, pp");
	qq = (BYTE **) my_malloc(sizeof(BYTE *) * b, "incma_latex_integer_labels, qq");
	pl.m_il(v);
	for (i = 0; i < v; i++) {
		a = point_labels->s_ii(i) + offset;
		convert_number_string(a, str);
		s.init(str);
		s.swap((STRING_OP) pl.s_i(i));
		ps = (STRING_OP) pl.s_i(i);
		pp[i] = ps->s_str();
		}
	bl.m_il(b);
	for (j = 0; j < b; j++) {
		a = block_labels->s_ii(j) + offset;
		convert_number_string(a, str);
		s.init(str);
		s.swap((STRING_OP) bl.s_i(j));
		ps = (STRING_OP) bl.s_i(j);
		qq[j] = ps->s_str();
		}
	incma_latex_picture(fp, 
		40 /* width */, 
		10 /* width_10 */,  
		FALSE /* f_outline_thin */, 
		"0.065mm" /* unit_length */, 
		"0.5mm" /* thick_lines */ , 
		"0.15mm" /* thin_lines */ , 
		"0.25mm" /* geo_line_width */ , 
		v, b, V, B, Vi, Bj, R, X, dim_X, 
		TRUE /* f_labelling_points */, pp, 
		TRUE /* f_labelling_blocks */, qq);
	my_free(pp);
	my_free(qq);
	return OK;
}


#if 0
	/* lin_5_1: */
	INT theX[][1] = {
		{ 0 } , 
		{ 0 } , 
		{ 0 } , 
		{ 0 } , 
		{ 0 } , 
		};
	INT R[5] = { 1, 1, 1, 1, 1};
	INT v = 5, b = 1, r = 1;
	INT V = 1, B = 1;
	INT Vi[] = { 5 };
	INT Bj[] = { 1 };
	INT a[] = { 1 };
#endif
#endif /* GEO_TRUE */

