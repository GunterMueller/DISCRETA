/* fg_table.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SOLVABLE_TRUE

#include <DISCRETA/solvable.h>

#define MAX_NW 64

#define DEBUG_SHRINK

#if TEXDOCU
INT FG_OB::int_shrink(FG_OP G, VECTOR_OP embedding, INT i)
#endif
{
	SHORT a[MAX_NW];
	INT ii;

	int2nw(i, a);
	NW_shrink(a, s_nb_ze_i() - 1, embedding);
	G->nw2int(a, &ii);
	return ii;
}

#if TEXDOCU
INT FG_OB::i2gen_idx(INT i)
#endif
{
	ZE_OP ze;
	INT j, l, n0;

	l = s_nb_ze_i();
	for (j = 0; j < l; j++) {
		ze = s_ze_i(j);
		n0 = ze->s_n0_i();
		if (i == n0)
			return j;
		}
	return error("FG::i2gen_idx() not a generator !");
}

#if TEXDOCU
INT FG_OB::int2nw(INT ii, SHORT *nw)
#endif
{
	ZE_OP ze;
	INT i, j;

	for (i = s_nb_ze_i() - 1; i >= 0; i--) {
		ze = s_ze_i(i);
		if (ii >= ze->s_n_i())
			return error("int2nw(): ii >= ze->s_n_i()");
		j = ii / ze->s_n0_i();
		nw[i] = j;
		ii -= j * ze->s_n0_i();
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::nw2int(SHORT *nw, INT *ii)
#endif
{
	ZE_OP ze;
	INT i, j, iii;
	
	iii = 0;
	for (i = s_nb_ze_i() - 1; i >= 0; i--) {
		ze = s_ze_i(i);
		j = ze->s_n0_i() * nw[i];
		iii += j;
		}
	*ii = iii;
	return OK;
}

/* 
 *
 */

#if TEXDOCU
INT FG_OB::mult_op(INTEGER_OP a, INTEGER_OP b, INTEGER_OP c)
#endif
{
	INT i;
	
	i = gt_mult(a->s_i(), b->s_i());
	c->m_i(i);
	return OK;
}

#if TEXDOCU
INT FG_OB::inv_op(INTEGER_OP a, INTEGER_OP b)
#endif
{
	INT i;
	
	i = gt_inv(a->s_i());
	b->m_i(i);
	return OK;
}

#if TEXDOCU
INT FG_OB::conj_op(INTEGER_OP a, INTEGER_OP b, INTEGER_OP c)
#endif
{
	INT i;
	
	i = gt_conj(a->s_i(), b->s_i());
	c->m_i(i);
	return OK;
}

#if TEXDOCU
INT FG_OB::commutator_op(INTEGER_OP a, INTEGER_OP b, INTEGER_OP c)
#endif
{
	INT i;
	
	i = gt_commutator(a->s_i(), b->s_i());
	c->m_i(i);
	return OK;
}

#if TEXDOCU
INT FG_OB::power_op(INTEGER_OP a, INT exp, INTEGER_OP c)
#endif
{
	INT i;
	
	i = gt_power(a->s_i(), exp);
	c->m_i(i);
	return OK;
}

#if TEXDOCU
INT FG_OB::one_op(INTEGER_OP a)
#endif
{
	a->m_i(0);
	return OK;
}

#if TEXDOCU
INT FG_OB::onep_op(INTEGER_OP a)
#endif
{
	return (a->s_i() == 0);
}

#if TEXDOCU
INT FG_OB::gt_mult(INT a, INT b)
#endif
{
	INT c;

	c = s_theG_iji(a, b);
	return c;
}

#if TEXDOCU
INT FG_OB::gt_inv(INT a)
#endif
{
	INT a_inv;

	a_inv = s_theG_inv_i(a);
	return a_inv;
}

#if TEXDOCU
INT FG_OB::gt_conj(INT a, INT b)
#else
/* $b^{-1} a b$ */
#endif
{
	INT b_inv, c;
	
	b_inv = s_theG_inv_i(b);
	c = gt_mult(b_inv, a);
	c = gt_mult(c, b);
	return c;
}

#if TEXDOCU
INT FG_OB::gt_commutator(INT a, INT b)
#else
/* $a^{-1} b^{-1} a b$ */
#endif
{
	INT a_inv, b_inv, c;
	
	a_inv = s_theG_inv_i(a);
	b_inv = s_theG_inv_i(b);
	c = gt_mult(a_inv, b_inv);
	c = gt_mult(c, a);
	c = gt_mult(c, b);
	return c;
}

#if TEXDOCU
INT FG_OB::gt_power(INT a, INT exp)
#endif
{
	INT c = 0, b = a;
	
	while (exp > 0) {
		/* res = b^exp * c */
		if (ODD(exp)) {
			c = gt_mult(b, c);
			exp--;
			continue; /* exp == 0 possible */
			}
		if (EVEN(exp)) {
			b = gt_mult(b, b);
			exp >>= 1;
			}
		}
	return c;
}

#if TEXDOCU
INT FG_OB::gt_order(INT i, INT *ord)
#endif
{
	INT j, ord1;

	
	j = i;
	ord1 = 1;
	while (j != 0) {
#if 0
		if (j >= s_n_i() || j < 0) {
			return error("FG_OB::gt_order illegal j");
			}
#endif
		j = gt_mult(i, j);
		ord1++;
		}
	*ord = ord1;
	return(OK);
}

#if TEXDOCU
INT FG_OB::gt_order_if_prime(INT i, INT *ord)
#endif
{
	INT ord1;
	
	gt_order(i, &ord1);
	if (ord1 == 1) {
		*ord = 0;
		return(OK);
		}
	if (smallest_primedivisor(ord1) == ord1)
		*ord = ord1;
	else
		*ord = 0;
	return(OK);
}

#if TEXDOCU
INT FG_OB::gt_order_if_prime_power(INT i, INT *ord, INT *prime, INT *k)
#endif
{
	INT ord1, prime1, prime2, k1;

	gt_order(i, &ord1);
	if (ord1 == 1) {
		*ord = 0;
		return(OK);
		}
	*ord = ord1;
	prime1 = 0;
	k1 = 0;
	while (ord1 != 1) {
		prime2 = smallest_primedivisor(ord1);
		if (prime1) {
			if (prime2 != prime1) {
				/* different primes */
				*ord = 0;
				return(OK);
				}
			}
		prime1 = prime2;
		k1++;
		ord1 /= prime1;
		}
	*prime = prime1;
	*k = k1;
	return(OK);
}

#undef CALC_AUT_DEBUG
#undef APPLY_AUT_DEBUG

#if TEXDOCU
INT FG_OB::calc_aut_i(INT i, PERMUTATION_OP aut, INT f_use_table)
#endif
{
	SHORT bi[MAX_NW];
	ZE_OP ze;
	INT n, j, k;

#ifdef CALC_AUT_DEBUG
	printf("calc_aut_i() i = %ld\n", i);
	fflush(stdout);
#endif
	ze = s_ze_i(i);
	n = ze->s_n0_i();
#ifdef CALC_AUT_DEBUG
	printf("calc_aut_i() n = %ld\n", n);
	fflush(stdout);
#endif
	aut->m_il(n);
	for (j = 0; j < i; j++) {
		k = ze->s_A_ii(j);
		bi[j] = k;
		}
	for ( ; j < s_nb_ze_i(); j++)
		bi[j] = 0;
	for (j = 0; j < n; j++) {
		k = apply_aut(j, bi, f_use_table);
		aut->m_ii(j, k + 1);
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::apply_aut(INT i, SHORT *base_im, INT f_use_table)
#endif
{
	SHORT a[MAX_NW];
	INT ii, j, k, res;

	res = 0;
#ifdef APPLY_AUT_DEBUG
	printf("apply_aut i = %ld\n", i);
	fflush(stdout);
#endif
	int2nw(i, a);
#ifdef APPLY_AUT_DEBUG
	NW_print(a, s_nb_ze_i());
	NW_print(base_im, s_nb_ze_i());
	fflush(stdout);
#endif
	/* for (ii = s_nb_ze_i() - 1; ii >= 0; ii--) { */
	for (ii = 0; ii < s_nb_ze_i(); ii++) {
		j = a[ii];
		for (k = 0; k < j; k++) {
			if (f_use_table)
				res = gt_mult(res, base_im[ii]);
			else
				res = mult(res, base_im[ii]);
			}
		}
	return res;
}

#if TEXDOCU
INT FG_OB::bi2aut(SHORT *base_im, PERMUTATION_OP aut, INT f_use_table)
#else
/* does not use the automorphism table of the group ! */
#endif
{
	INT i, j, n;
	
	n = s_n_i();
	aut->m_il(n);
	for (i = 0; i < n; i++) {
		j = apply_aut(i, base_im, f_use_table);
		}
	aut->m_ii(i, j + 1);
	return OK;
}

#if TEXDOCU
INT FG_OB::inv(INT i)
#endif
{
	SHORT a[MAX_NW];
	SHORT b[MAX_NW];
	INT j;
	
	int2nw(i, a);
	nw_inv(s_nb_ze_i() - 1, a, b);
	nw2int(b, &j);
	return j;
}

#if TEXDOCU
INT FG_OB::mult(INT i, INT j)
#endif
{
	SHORT a[MAX_NW];
	SHORT b[MAX_NW];
	SHORT c[MAX_NW];
	INT k;
	
	int2nw(i, a);
	int2nw(j, b);
	nw_mult(s_nb_ze_i() - 1, a, b, c);
	nw2int(c, &k);
	return k;
}

#if TEXDOCU
INT FG_OB::conjugate(INT i, INT j)
#else
/* $j^{-1} * i * j$ */
#endif
{
	SHORT nw_i[MAX_NW];
	SHORT nw_j[MAX_NW];
	SHORT nw_jv[MAX_NW];
	SHORT a[MAX_NW];
	SHORT b[MAX_NW];
	INT k;
	
	int2nw(i, nw_i);
	int2nw(j, nw_j);
	nw_inv(s_nb_ze_i() - 1, nw_j, nw_jv);
	nw_mult(s_nb_ze_i() - 1, nw_jv, nw_i, a);
	nw_mult(s_nb_ze_i() - 1, a, nw_j, b);
	nw2int(b, &k);
	return k;
}

#if TEXDOCU
INT FG_OB::commutator(INT i, INT j)
#else
/* $i^{-1} * j^{-1} * i * j$ */
#endif
{
	SHORT nw_i[MAX_NW];
	SHORT nw_j[MAX_NW];
	SHORT nw_iv[MAX_NW];
	SHORT nw_jv[MAX_NW];
	SHORT a[MAX_NW];
	SHORT b[MAX_NW];
	INT k;
	
	int2nw(i, nw_i);
	int2nw(j, nw_j);
	nw_inv(s_nb_ze_i() - 1, nw_i, nw_iv);
	nw_inv(s_nb_ze_i() - 1, nw_j, nw_jv);
	nw_mult(s_nb_ze_i() - 1, nw_iv, nw_jv, a);
	nw_mult(s_nb_ze_i() - 1, a, nw_i, b);
	nw_mult(s_nb_ze_i() - 1, b, nw_j, a);
	nw2int(a, &k);
	return k;
}

#if TEXDOCU
INT FG_OB::power(INT i, INT exp)
#endif
{
	SHORT a[MAX_NW];
	SHORT b[MAX_NW];
	INT j;
	
	int2nw(i, a);
	nw_power(s_nb_ze_i() - 1, a, b, exp);
	nw2int(b, &j);
	return j;
}

#if TEXDOCU
INT FG_OB::nw_inv(INT i, SHORT *a, SHORT *b)
#endif
{
	SHORT c[MAX_NW];
	SHORT d[MAX_NW];
	SHORT e[MAX_NW];
	ZE_OP ze;
	INT j;
	
	if (i == 0) {
		ze = s_ze_i(0);
		b[0] = (ze->s_p_i() - a[0]) % ze->s_p_i();
		return OK;
		}
#ifdef DB_INV
	printf("FG::nw_inv()| a = "); NW_print(a, i);
#endif
	nw_reduce(i, a);
	NW_one(c, i);
	for (j = i; j >= 0; j--) {
		ze = s_ze_i(j);
		if (j > 0) {
			int2nw(ze->s_P_i(), e);
			nw_inv(j - 1, e, d);
			NW_pad(d, j - 1, i);
			}
		else {
			NW_one(d, i);
			}
		d[j] = (SHORT) (ze->s_p_i() - (INT) a[j]);
		nw_mult(i, c, d, e);
		NW_copy(e, c, i);
		}
	NW_copy(c, b, i);
#ifdef DB_INV
	printf("FG::nw_inv()| b = "); NW_print(b, i);
#endif
	return OK;
}

#if TEXDOCU
INT FG_OB::nw_mult(INT i, SHORT *a, SHORT *b, SHORT *c)
#endif
{
	INT k;
	SHORT d[MAX_NW];
	SHORT e[MAX_NW];
	ZE_OP ze;
	
	if (i == 0) {
		ze = s_ze_i(0);
		c[0] = (a[0] + b[0]) % ze->s_p_i();
		return OK;
		}
	nw_vorbei_v(b, d, i, (INT) a[i]);
	k = (INT) d[i];
	nw_mult(i - 1, a, d, e);
	e[i] = (SHORT) k;
	nw_reduce(i, e);
	NW_copy(e, c, i);
	return OK;
}

INT db_vorbei = FALSE;

#if TEXDOCU
INT FG_OB::nw_vorbei_v(SHORT *a, SHORT *b, INT i, INT vorne)
#endif
{
	INT j;
	SHORT c[MAX_NW];
	SHORT d[MAX_NW];

	if (vorne == 0) {
		NW_copy(a, b, i);
		return OK;
		}
#ifdef NS_DEBUG_NSW_VORBEI
	printf("FG::nw_vorbei_v()| vorne = %ld a = ", vorne);
	NW_print(a, i);
#endif
	NW_copy(a, c, i);
	for (j = 1; j <= vorne; j++) {
		nw_vorbei(c, d, i);
		NW_copy(d, c, i);
		}
	NW_copy(c, b, i);
	return OK;
}

#if TEXDOCU
INT FG_OB::nw_vorbei(SHORT *a, SHORT *b, INT i)
#endif
{
	INT j, k;
	SHORT c[MAX_NW];
	SHORT d[MAX_NW];
	SHORT e[MAX_NW];
	ZE_OP ze;
	
	if (db_vorbei) {
		printf("FG::nw_vorbei()| a = ");
		NW_print(a, i);
		}
	NW_one(c, i);
	ze = s_ze_i(i);
	for (j = 0; j < i; j++) {
		k = (INT) a[j];
		if (db_vorbei) {
			printf("FG::nw_vorbei()| j = %ld k = %ld\n", j, k);
			}
		int2nw(ze->s_Av_ii(j), e);
		nw_power(i - 1, e, d, k);
		if (db_vorbei) {
			NW_print(e, i - 1);
			printf("FG::vorbei()|nach nw_power(): d = ");
			NW_print(d, i - 1);
			}
		nw_mult(i - 1, c, d, e);
		if (db_vorbei) {
			printf("FG::nw_vorbei()|nach nw_mult(c, d, e): e = ");
			NW_print(e, i - 1);
			}
		NW_copy(e, c, i - 1);
		}
	c[i] = (SHORT) (a[i] + 1);
	nw_reduce(i, c);
	NW_copy(c, b, i);
	if (db_vorbei) {
		printf("FG::nw_vorbei()| b = ");
		NW_print(b, i);
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::nw_reduce(INT i, SHORT *a)
#endif
{
	SHORT b[MAX_NW];
	SHORT c[MAX_NW];
	SHORT d[MAX_NW];
	INT u_i, u_q, u_r;
	ZE_OP ze;
	
	ze = s_ze_i(i);
	if (i == 0) {
		a[0] = (SHORT) ((INT) a[0] % ze->s_p_i());
		return OK;
		}
	u_i = (INT) a[i];
	nw_reduce(i - 1, a);
	u_r = u_i % ze->s_p_i();
	u_q = (u_i - u_r) / ze->s_p_i();
	int2nw(ze->s_P_i(), d);
	nw_power(i - 1, d, b, u_q);
	nw_mult(i - 1, a, b, c);
	c[i] = (SHORT) u_r;
	NW_copy(c, a, i);
	return OK;
}

#if TEXDOCU
INT FG_OB::nw_power_ip(INT i, SHORT *a, INT exp)
#else
/* ip = in place */
#endif
{
	SHORT b[MAX_NW];
	
	nw_power(i, a, b, exp);
	NW_copy(b, a, i);
	return OK;
}

#if TEXDOCU
INT FG_OB::nw_power(INT i, SHORT *a, SHORT *res, INT exp)
#endif
{
	SHORT b[MAX_NW];
	SHORT c[MAX_NW];
	SHORT d[MAX_NW];
	
#ifdef DB_POWER
	printf("FG::nw_power()| exp = %ld a = ", exp);
	NW_print(a, i);
#endif
	NW_copy(a, b, i);
	NW_one(c, i);
	while (exp > 0) {
		/* res = c * b^exp */
		if (ODD(exp)) {
			nw_mult(i, b, c, d);
			NW_copy(d, c, i);
			exp--;
			continue; /* exp == 0 possible */
			}
		if (EVEN(exp)) {
			nw_mult(i, b, b, d);
			NW_copy(d, b, i);
			exp >>= 1;
			}
		}
	NW_copy(c, res, i);
#ifdef DB_POWER
	printf("FG::power()| res = ");
	NW_print(res, i);
#endif
	return OK;
}

/*
 * NW:
 */

#if TEXDOCU
INT NW_one(SHORT *a, INT i)
#endif
{
	INT j;
	
	for (j = 0; j <= i; j++)
		a[j] = 0;
	return OK;
}

#if TEXDOCU
INT NW_is_one(SHORT *a, INT i)
#endif
{
	INT j;
	
	for (j = 0; j <= i; j++) {
		if (a[j])
			return FALSE;
		}
	return TRUE;
}

#if TEXDOCU
INT NW_copy(SHORT *a, SHORT *b, INT i)
#endif
{
	INT j;
	
	for (j = 0; j <= i; j++)
		b[j] = a[j];
	return OK;
}

#if TEXDOCU
INT NW_2_str(SHORT *a, INT i, BYTE *str)
#else
/* haengt an str an */
#endif
{
	INT j;
	BYTE s1[64];
	
	if (NW_is_one(a, i)) {
		sprintf(str + strlen(str), "id");
		}
	else {
		for (j = 0; j <= i; j++) {
			if (a[j] == 0)
				continue;
			s1[0] = 'A' + j;
			s1[1] = 0;
			if (a[j] != 1) {
				sprintf(
				s1 + strlen(s1), "^%ld", (INT)a[j]);
				}
			strcat(str, s1);
			}
		}
	return OK;
}

#if TEXDOCU
INT NW_2_str_GAP(SHORT *a, INT i, INT nb_ze, BYTE *str)
#else
/* haengt an str an */
#endif
{
	INT j, f_first = TRUE;
	BYTE s1[1024];
	
	if (NW_is_one(a, i)) {
		sprintf(str + strlen(str), "IdWord");
		}
	else {
		for (j = 0; j <= i; j++) {
			if (a[j] == 0)
				continue;
			s1[0] = 0;
			if (!f_first)
				sprintf(s1 + strlen(s1), "*");
			sprintf(s1 + strlen(s1), "g%ld", nb_ze - j);
			if (a[j] != 1) {
				sprintf(
				s1 + strlen(s1), "^%ld", (INT)a[j]);
				}
			strcat(str, s1);
			f_first = FALSE;
			}
		}
	return OK;
}

#if TEXDOCU
INT NW_2_str_tex(SHORT *a, INT i, BYTE *str)
#else
/* haengt an str an */
#endif
{
	INT j;
	BYTE s1[64];
	
	if (NW_is_one(a, i)) {
		sprintf(str + strlen(str), "id");
		}
	else {
		for (j = 0; j <= i; j++) {
			if (a[j] == 0)
				continue;
			s1[0] = 'A' + j;
			s1[1] = 0;
			if (a[j] != 1) {
				sprintf(
				s1 + strlen(s1), "^{%ld}", (INT)a[j]);
				}
			strcat(str, s1);
			}
		}
	return OK;
}

#if TEXDOCU
INT NW_shrink(SHORT *a, INT i, VECTOR_OP embedding)
#endif
{
	INT k, k0, l;
	INT forgotten[MAX_NW];

	l = embedding->s_li();
#ifdef DEBUG_SHRINK
	for (k = 0; k <= i; k++) {
		forgotten[k] = TRUE;
		}
	for (k = 0; k < l; k++) {
		k0 = embedding->s_ii(k);
		forgotten[k0] = FALSE;
		}
	for (k = 0; k <= i; k++) {
		if (forgotten[k]) {
			if (a[k])
				return error("NW_shrink(): a[k] != 0 though k should be forgotten !");
			}
		}
#endif
	for (k = 0; k < l; k++) {
		k0 = embedding->s_ii(k);
		a[k] = a[k0];
		}
	for ( ; k <= i; k++)
		a[k] = 0;
	return OK;
}

#if TEXDOCU
INT NW_pad(SHORT *a, INT j, INT i)
#endif
{
	INT k;
	
	for (k = j + 1; k <= i; k++)
		a[k] = 0;
	return OK;
}

#if TEXDOCU
INT NW_print(SHORT *a, INT i)
#endif
{
	INT j;
	
	for (j = 0; j < i; j++) {
		printf("%ld ", (INT)a[j]);
		}
	printf("\n");
	return OK;
}

#if TEXDOCU
INT print_vec_of_fg(VECTOR_OP V)
#endif
{
	FG_OP G;
	INT i;
	
	for (i = 0; i < V->s_li(); i++) {
		G = (FG_OP) V->s_i(i);
		printf("%ld:\n", i);
		G->Print();
		}
	return OK;
}

/*
 *
 */

#if TEXDOCU
INT FG_OB::left_regular_representation(INT i, PERMUTATION_OP p)
#else
/* Multiplikationsdarstellung des Elementes i auf allen 
 * Gruppenelementen: 0 <= elemente < n.
 * Die Permutation hat Eintraege 1 <= eintraege <= n. */
#endif
{
	INT ii, j, jj, n;
	
	n = s_n_i();
	p->m_il(n);
	for (ii = 0; ii < n; ii++) {
		j = mult(i, ii);
		for (jj = 0; jj < ii; jj++) {
			if (p->s_ii(jj) - 1 == j)
				return error("left_regular_representation()|not a permutation");
			}
		p->m_ii(ii, j + 1);
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::right_regular_representation(INT i, PERMUTATION_OP p)
#endif
{
	INT ii, j, jj, n;
	
	n = s_n_i();
	p->m_il(n);
	for (ii = 0; ii < n; ii++) {
		j = mult(ii, i);
		for (jj = 0; jj < ii; jj++) {
			if (p->s_ii(jj) - 1 == j)
				return error("right_regular_representation()|not a permutation");
			}
		p->m_ii(ii, j + 1);
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::inn_generator(INT i, PERMUTATION_OP p, INT f_use_table)
#endif
{
	INT ii, j, jj, n;

	n = s_n_i();
	p->m_il(n);
	for (ii = 0; ii < n; ii++) {
		if (f_use_table) 
			j = gt_conj(ii, i);
		else
			j = conjugate(ii, i);
		p->m_ii(ii, j + 1);
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::theG(INT f_verbose)
#endif
{
	SHORT d[MAX_NW];
	INT n, i, j, k, dim_n1, idx_inv, idx_inv1;
	INT erg = OK;
	
	n = s_n_i();
	dim_n1 = n + 5;
	s_dim_n1()->m_i(dim_n1);
	erg += s_theG()->m_ilih_n(dim_n1, n);
	for (i = 0; i < n; i++) {
		idx_inv = -1;
		for (j = 0; j < n; j++) {
			k = mult(i, j);
			if (k == 0)
				idx_inv = j;
			s_theG_ij(i, j)->m_i(k);
			if (f_verbose)
				printf("%ld ", k);
			}
		if (idx_inv == -1)
			return error("FG::theG()|inverse element not found");
		idx_inv1 = inv(i);
		if (idx_inv1 != idx_inv) {
			printf("\ni = %ld ", i);
			print_int(i);
			printf("inv(i) = %ld ", idx_inv1);
			print_int(idx_inv1);
			printf(" idx_inv = %ld ", idx_inv);
			print_int(idx_inv);
			printf("\n");
			fflush(stdout);
			return error("FG::theG()|inv() wrong");
			}
		s_theG_inv(i)->m_i(idx_inv);
		if (f_verbose)
			printf(" %ld\n", idx_inv);
		}
	return erg;
}

#undef TEST_RECALC_TABLE

#if TEXDOCU
INT FG_OB::recalc_table(INT f_v)
#endif
{
	PERMUTATION_OB aut, autv;
	ZE_OP ze;
	INT n, i, j, ii, jj, i1, j1, idx, idx_inv, dim_n1;
	INT k, n0, p, P;
	MATRIX_OP G;

	n = s_n_i();
	dim_n1 = n + 5;
	s_dim_n1()->m_i(dim_n1);
	s_theG()->m_ilih_n(dim_n1, n);
	G = s_theG();
	G->m_iji(0, 0, 0);
	if (f_v) {
		printf("recalc_table\n");
		fflush(stdout);
		}
	for (k = 0; k < s_nb_ze_i(); k++) {
		if (f_v) {
			printf("k = %ld\n", k);
			fflush(stdout);
			}
		calc_aut_i(k, &aut, TRUE /* f_use_table */);
		if (f_v) {
			printf("aut = ");
			aut.println();
			fflush(stdout);
			}
		aut.invers(&autv);
		ze = s_ze_i(k);
		n0 = ze->s_n0_i();
		p = ze->s_p_i();
		P = ze->s_P_i();
		for (ii = 0; ii < p; ii++) {
			for (i = 0; i < n0; i++) {
				i1 = ii * n0 + i;
				idx_inv = -1;
				for (jj = 0; jj < p; jj++) {
					for (j = 0; j < n0; j++) {
						if (ii == 0 && jj == 0) {
							if (f_v)
								printf("    ");
							continue;
							}
						j1 = jj * n0 + j;
						idx = mult_Zp_extension2(n0, p, P, i, ii, j, jj, &aut, &autv);
						G->m_iji(i1, j1, idx);
						if (f_v)
							printf("%3ld ", idx);
						if (idx == 0)
							idx_inv = j1;
						}
					}
				if (ii != 0) {
					if (idx_inv == -1) {
						return error("recalc_table() inverse element not found");
						}
					s_theG_inv(i1)->m_i(idx_inv);
					if (f_v)
						printf(" %ld", idx_inv);
					}
				if (f_v)
					printf("\n");
				}
			}
		} /* next k */
#ifdef TEST_RECALC_TABLE
	{
		MATRIX_OB G0;
		INT a, b;
		
		s_theG()->swap(&G0);
		theG(FALSE);
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				a = G0.s_iji(i, j);
				b = s_theG_iji(i, j);
				if (a != b) {
					printf("recalc_table() tables differ at %ld %ld\n", i, j);
					G0.Print();
					s_theG()->Print();
					return error("recalc_table() tables differ");
					}
				}
			}
	}
#endif
	return OK;
}

#if TEXDOCU
INT FG_OB::theG_Zp_extension(FG_OP G0, PERMUTATION_OP A, INT f_verbose)
#endif
{
	PERMUTATION_OB Av;
	ZE_OP ze;
	INT n, i, j, ii, jj, i1, j1, idx, idx_inv, dim_n1;
	INT erg = OK;
	
	A->invers(&Av);
	n = s_n_i();
	dim_n1 = n + 5;
	s_dim_n1()->m_i(dim_n1);
	erg += s_theG()->m_ilih_n(dim_n1, n);
	ze = s_ze_i(s_nb_ze_i() - 1);
	for (ii = 0; ii < ze->s_p_i(); ii++) {
		for (i = 0; i < G0->s_n_i(); i++) {
			i1 = ii * G0->s_n_i() + i;
			idx_inv = -1;
			for (jj = 0; jj < ze->s_p_i(); jj++) {
				for (j = 0; j < G0->s_n_i(); j++) {
					j1 = jj * G0->s_n_i() + j;
					idx = mult_Zp_extension(G0, i, ii, j, jj, A, &Av);
					s_theG_ij(i1, j1)->m_i(idx);
					/* pi[i1 * dim_n1 + j1] = idx; */
					if (f_verbose)
						printf("%ld ", idx);
					if (idx == 0)
						idx_inv = j1;
					}
				}
			if (idx_inv == -1) {
				Srfs("fg_theG_Zp_extension", "inverse element not found");
				return ERROR;
				}
			s_theG_inv(i1)->m_i(idx_inv);
			if (f_verbose)
				printf(" %ld\n", idx_inv);
			}
		}
	return erg;
}

#if TEXDOCU
INT FG_OB::mult_Zp_extension(FG_OP G0, 
	INT i, INT ii, INT j, INT jj, 
	PERMUTATION_OP A, PERMUTATION_OP Av)
#else
/* berechnet $i * z^{ii} * j * z^{jj}$. */
#endif
{
	ZE_OP ze;
	INT i1, j1, jj1, l, P1, u_r, u_q;
	
	/* z^ii an j vorbeiziehen: 
	 * z^ii * j * z^jj = j1 * z^jj1 */
	j1 = j;
	for (l = 0; l < ii; l++) {
		j1 = Av->s_ii(j1) - 1;
		}
	i1 = G0->gt_mult(i, j1);
	jj1 = jj + ii;
	ze = s_ze_i(s_nb_ze_i() - 1);
	u_r = jj1 % ze->s_p_i();
	u_q = (jj1 - u_r) / ze->s_p_i();
	P1 = G0->gt_power(ze->s_P_i(), u_q);
	i1 = G0->gt_mult(i1, P1);
	i1 += u_r * G0->s_n_i();
	return i1;
}

#if TEXDOCU
INT FG_OB::mult_Zp_extension2(INT n0, INT p, INT P, 
	INT i, INT ii, INT j, INT jj, 
	PERMUTATION_OP A, PERMUTATION_OP Av)
#else
/* computes: $i * z^{ii} * j * z^{jj}$. */
#endif
{
	INT i1, j1, jj1, l, P1, u_r, u_q;
	
	/* z^ii an j vorbeiziehen: 
	 * z^ii * j * z^jj = j1 * z^jj1 */
	j1 = j;
	for (l = 0; l < ii; l++) {
		j1 = Av->s_ii(j1) - 1;
		}
	i1 = gt_mult(i, j1);
	jj1 = jj + ii;
	u_r = jj1 % p;
	u_q = (jj1 - u_r) / p;
	P1 = gt_power(P, u_q);
	i1 = gt_mult(i1, P1);
	i1 += u_r * n0;
	return i1;
}

#if TEXDOCU
INT FG_OB::is_associative(INT f_verbose)
#endif
{
	INT n, i, j, k, ij, jk, ijk1, ijk2;
	INT i_inv, id;
	
	n = s_n_i();
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			ij = gt_mult(i, j);
			for (k = 0; k < n; k++) {
				jk = gt_mult(j, k);
				ijk1 = gt_mult(ij, k);
				ijk2 = gt_mult(i, jk);
				if (ijk1 != ijk2) {
					if (f_verbose) {
						printf("not associative: "
							"i=%ld j=%ld k=%ld\n", 
							i, j, k);
						printf("(i*j)*k= %ld*k = %ld\n", 
							ij, ijk1);
						printf("i*(j*k)= i*%ld = %ld\n", 
							jk, ijk2);
						}
					return FALSE;
					}
				}
			}
		}
	for (i = 0; i < n; i++) {
		i_inv = gt_inv(i);
		id = gt_mult(i, i_inv);
		if (id != 0) {
			if (f_verbose)
				printf("not right-inverse: "
				"i=%ld i_inv=%ld i*i_inv=%ld\n", 
				i, i_inv, id);
				return FALSE;
			}
		id = gt_mult(i_inv, i);
		if (id != 0) {
			if (f_verbose)
				printf("not left-inverse: "
				"i=%ld i_inv=%ld i_inv*i=%ld\n", 
				i, i_inv, id);
				return FALSE;
			}
		}
	return TRUE;
}

#if TEXDOCU
INT FG_OB::Inn_generators(VECTOR_OP gen, INT f_use_table)
#endif
{
	PERMUTATION_OB p;
	INT i, l;
	
	l = s_nb_ze_i();
	gen->m_il(l);
	for (i = 0; i < l; i++) {
		inn_generator(g_i(i), &p, f_use_table);
		p.swap((PERMUTATION_OP) gen->s_i(i));
		}
	return OK;
}


#endif /* SOLVABLE_TRUE */


