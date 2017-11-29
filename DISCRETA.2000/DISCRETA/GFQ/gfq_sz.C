/* gfq_sz.C 
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
#include <DISCRETA/lb.h>

#if TEXDOCU
#else
generators for Sz(q)
according to Huppert and Blackburn, 
Finite Groups III~\cite{HuppertBlackburn82b}, p. 182.
#endif

void z_mtx4_gen_Mlambda_num(ZECH_DATA *Z, INT *mtx, INT lambda)
{
	INT x[4], lambda_inv, lambda_2m, lambda_inv_2m;

	lambda_inv = z_inverse_num(Z, lambda);
	lambda_2m = z_apply_frob_sz0_num(Z, lambda);
	lambda_inv_2m = z_apply_frob_sz0_num(Z, lambda_inv);
	x[0] = z_mult_num(Z, lambda, lambda_2m);
	x[1] = lambda_2m;
	x[2] = lambda_inv_2m;
	x[3] = z_mult_num(Z, lambda_inv, lambda_inv_2m);
	z_mtx_i_gen_diag_num(Z, mtx, x, 4);
#ifdef DEBUG_Z_MTX4_GEN_MLAMBDA_NUM
	z_mtx4_print(Z, mtx);
	fflush(stdout);
#endif
}

void z_mtx4_gen_Sab_num(ZECH_DATA *Z, INT *mtx, INT a, INT b)
{
	INT i, j, *p;
	INT a_pi, b_pi, aa, ab, d, c;

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			p = mtx + i * 4 + j;
			if (i < j)
				*p = z_zero_num(Z);
			else if (i == j)
				*p = z_one_num(Z);
			}
		}
	a_pi = z_apply_frob_sz_num(Z, a);
	b_pi = z_apply_frob_sz_num(Z, b);
	aa = z_apply_frob_num(Z, a);
	ab = z_mult_num(Z, a, b);
	d = z_mult_num(Z, aa, a_pi);
	d = z_add_num(Z, d, ab);
	d = z_add_num(Z, d, b_pi);
	c = z_mult_num(Z, a, a_pi);
	c = z_add_num(Z, c, b);
	mtx[1 * 4 + 0] = a;
	mtx[2 * 4 + 0] = b;
	mtx[2 * 4 + 1] = a_pi;
	mtx[3 * 4 + 0] = d;
	mtx[3 * 4 + 1] = c;
	mtx[3 * 4 + 2] = a;
#ifdef DEBUG_Z_MTX4_GEN_SAB_NUM
	z_mtx4_print(Z, mtx);
	fflush(stdout);
#endif
}

void z_mtx4_gen_T_num(ZECH_DATA *Z, INT *mtx)
{
	INT i, j, *p, elt_0, elt_1;

	elt_0 = z_zero_num(Z);
	elt_1 = z_one_num(Z);
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			p = mtx + i * 4 + j;
			*p = elt_0;
			}
		}
	mtx[0 * 4 + 3] = elt_1;
	mtx[1 * 4 + 2] = elt_1;
	mtx[2 * 4 + 1] = elt_1;
	mtx[3 * 4 + 0] = elt_1;
#ifdef DEBUG_Z_MTX4_GEN_T_NUM
	z_mtx4_print(Z, mtx);
	fflush(stdout);
#endif
}

TITS_OVOID_REP *tr_open(INT deg, INT chi)
{
	TITS_OVOID_REP *tr;
	INT k, len, j, idx;

	tr = (TITS_OVOID_REP *) my_malloc(sizeof(TITS_OVOID_REP), "gfq_sz.C: tr_open");
	tr->f_has_omega = FALSE;
	tr->Z = zech_open(deg, chi, TRUE);
	
	tr_gen_tits_ovoid(tr);
	printf("size of ovoid: %ld\n", tr->nb_o);

#if 0
	for (k = 0; k < tr->nb_o; k++) {
		printf("%3ld: %2ld %2ld %2ld %2ld\n", k, 
			tr->O[k * 4 + 0], 
			tr->O[k * 4 + 1], 
			tr->O[k * 4 + 2], 
			tr->O[k * 4 + 3]);
		}
#endif
	tr_gen_omega(tr);
#if 1
	for (k = 0; k < tr->omega_h; k++) {
		printf("---orbit %ld:\n", k);
		len = tr->omega_len[k];
		for (j = 0; j < len; j++) {
			idx = tr->omega[k * tr->omega_l + j];
			printf("%3ld: %2ld %2ld %2ld %2ld\n", idx, 
				tr->O[idx * 4 + 0], 
				tr->O[idx * 4 + 1], 
				tr->O[idx * 4 + 2], 
				tr->O[idx * 4 + 3]);
			}
		}
#endif
	tr->f_has_omega = FALSE;
	
	return tr;
}

void tr_close(TITS_OVOID_REP *tr)
{
	if (tr->O)
		my_free(tr->O);
	tr->nb_o = 0;
	if (tr->f_has_omega) {
		my_free(tr->omega_idx);
		my_free(tr->omega);
		my_free(tr->omega_len);
		}
	zech_free(tr->Z);
	my_free(tr);
}

INT tr_gen_tits_ovoid(TITS_OVOID_REP *tr)
{
	ZECH_DATA *Z = tr->Z;
	INT *O, nb_o;
	INT x, y, x_pi, y_pi, xx, xy, d;

	nb_o = 0;
	O = (INT *) my_malloc(4 * (Z->q * Z->q + 1) * sizeof(INT), "gfq_sz.C: tr_gen_tits_ovoid");
	O[nb_o * 4 + 0] = z_one_num(Z);
	O[nb_o * 4 + 1] = z_zero_num(Z);
	O[nb_o * 4 + 2] = z_zero_num(Z);
	O[nb_o * 4 + 3] = z_zero_num(Z);
	nb_o++;
	for (x = 0; x < Z->q; x++) {
		xx = z_apply_frob_num(Z, x);
		x_pi = z_apply_frob_sz_num(Z, x);
		for (y = 0; y < Z->q; y++) {
			y_pi = z_apply_frob_sz_num(Z, y);
			xy = z_mult_num(Z, x, y);
			d = z_mult_num(Z, x_pi, xx);
			d = z_add_num(Z, d, xy);
			d = z_add_num(Z, d, y_pi);
			O[nb_o * 4 + 0] = d;
			O[nb_o * 4 + 1] = y;
			O[nb_o * 4 + 2] = x;
			O[nb_o * 4 + 3] = z_one_num(Z);
			nb_o++;
			}
		}
	tr->O = O;
	tr->nb_o = nb_o;
	return OK;
}

INT tr_gen_omega(TITS_OVOID_REP *tr)
{
	ZECH_DATA *Z = tr->Z;
	INT *omega_idx, *omega, *omega_len;
	INT omega_l, omega_h, deg, i, j, k, idx;
	INT x[4];
	
	deg = Z->deg;
	omega_idx = (INT *) my_malloc(tr->nb_o * sizeof(INT), "tr_gen_omega");
	omega = (INT *) my_malloc(tr->nb_o * deg * sizeof(INT), "tr_gen_omega");
	omega_len = (INT *) my_malloc(tr->nb_o * sizeof(INT), "tr_gen_omega");
	for (i = 0; i < tr->nb_o; i++) {
		omega_idx[i] = -1;
		omega_len[i] = 0;
		}
	omega_h = 0;
	omega_l = deg;
	for (i = 0; i < tr->nb_o; i++) {
		if (omega_idx[i] != -1) /* already processed ? */
			continue; /* next i */
		x[0] = tr->O[i * 4 + 0];
		x[1] = tr->O[i * 4 + 1];
		x[2] = tr->O[i * 4 + 2];
		x[3] = tr->O[i * 4 + 3];
		omega_idx[i] = omega_h;
		omega[omega_h * deg + 0] = i;
		omega_len[omega_h] = 1;
		for (j = 0; j < deg; j++) {
			for (k = 0; k < 4; k++)
				x[k] = z_apply_frob_num(Z, x[k]);
			/* no need to normalize x[] ! */
			idx = tr_search(tr, x);
			if (idx < 0)
				return error("tr_gen_omega(): "
					"vector not found");
			if (idx == i) /* we are through */
				break;
			if (omega_len[omega_h] >= deg)
				return error("tr_gen_omega(): "
					"orbit too large !");
			if (omega_idx[idx] != -1)
				return error("tr_gen_omega(): "
					"omega_idx[idx] != -1");
			omega[omega_h * deg + omega_len[omega_h]] = idx;
			omega_len[omega_h]++;
			omega_idx[idx] = omega_h;
			}

		omega_h++;
		}
	printf("tr_gen_omega(): found %ld orbits of "
		"x \\mapsto x^2 on the ovoid !\n", omega_h);
	tr->omega_idx = omega_idx;
	tr->omega = omega;
	tr->omega_len = omega_len;
	tr->omega_l = omega_l;
	tr->omega_h = omega_h;
	tr->f_has_omega = TRUE;
	return OK;
}

INT tr_cmp_vec(INT *v1, INT *v2)
/* the projective elements have normalized 
 * coordinates;
 * we assume that we have numeric form of the elements 
 * so that the lexicographic sorting reflects the 
 * sorting according to the p-adic representation */
{
	INT i, a, b;

	for (i = 3; i >= 0; i--) {
		a = v1[i];
		b = v2[i];
		if (a < b)
			return -1;
		if (a > b)
			return 1;
		}
	return 0;
}

INT tr_search(TITS_OVOID_REP *tr, INT *x)
{
	INT i, j, mid, res;
	
	i = 0;
	j = tr->nb_o - 1;
	while (i <= j) {
		mid = i + (j - i) / 2;
		res = tr_cmp_vec(tr->O + mid * 4, x);
		if (res == 0)
			return mid;
		else if (res < 0)
			i = mid + 1;
		else
			j = mid - 1;
		}
	return -1;
}

void tr_vec4_normalize_num(TITS_OVOID_REP *tr, INT *v)
{
	ZECH_DATA *Z = tr->Z;
	INT a;

	if (!z_is_zero_num(Z, v[3])) {
		a = z_inverse_num(Z, v[3]);
		v[3] = z_mult_num(Z, a, v[3]);
		v[2] = z_mult_num(Z, a, v[2]);
		v[1] = z_mult_num(Z, a, v[1]);
		v[0] = z_mult_num(Z, a, v[0]);
		return;
		}
	if (!z_is_zero_num(Z, v[2])) {
		a = z_inverse_num(Z, v[2]);
		v[2] = z_mult_num(Z, a, v[2]);
		v[1] = z_mult_num(Z, a, v[1]);
		v[0] = z_mult_num(Z, a, v[0]);
		return;
		}
	if (!z_is_zero_num(Z, v[1])) {
		a = z_inverse_num(Z, v[1]);
		v[1] = z_mult_num(Z, a, v[1]);
		v[0] = z_mult_num(Z, a, v[0]);
		return;
		}
	if (!z_is_zero_num(Z, v[0])) {
		a = z_inverse_num(Z, v[0]);
		v[0] = z_mult_num(Z, a, v[0]);
		return;
		}
}

INT tr_mtx4_perm(TITS_OVOID_REP *tr, 
	INT *mtx, PERMUTATION_OP perm)
{
	ZECH_DATA *Z = tr->Z;
	INT i, j, k, k1;
	INT v4[4], w4[4], x0, a;

	perm->m_il(tr->nb_o);
	x0 = z_zero_num(Z);
	for (k = 0; k < tr->nb_o; k++) {
		v4[0] = tr->O[k * 4 + 0];
		v4[1] = tr->O[k * 4 + 1];
		v4[2] = tr->O[k * 4 + 2];
		v4[3] = tr->O[k * 4 + 3];
		for (j = 0; j < 4; j++) {
			w4[j] = x0;
			for (i = 0; i < 4; i++) {
				a = z_mult_num(Z, v4[i], mtx[i * 4 + j]);
				w4[j] = z_add_num(Z, w4[j], a);
				}
			}
#ifdef DEBUG_MTX4_PERM_VERBOSE
		printf("w4 = ");
		printf("%2ld %2ld %2ld %2ld\n", w4[0], w4[1], w4[2], w4[3]);
#endif
		tr_vec4_normalize_num(tr, w4);
#ifdef DEBUG_MTX4_PERM_VERBOSE
		printf("%2ld %2ld %2ld %2ld\n", w4[0], w4[1], w4[2], w4[3]);
#endif
		k1 = tr_search(tr, w4);
		if (k1 < 0)
			return error("tr_mtx4_perm(): "
				"vector not found");
		perm->m_ii(k, k1 + 1);
		}
#ifdef DEBUG_MTX4_PERM
	perm->print();
#endif
	return OK;
}

INT tr_mtx4_perm_omega(TITS_OVOID_REP *tr, 
	INT *mtx, PERMUTATION_OP perm)
{
	ZECH_DATA *Z = tr->Z;
	INT i, j, k, k1, idx1, idx2;
	INT v4[4], w4[4], x0, a;

	perm->m_il(tr->omega_h);
	x0 = z_zero_num(Z);
	for (k = 0; k < tr->omega_h; k++) {
		((INTEGER_OP) perm->s_i(k))->m_i(-1);
		}
	for (k = 0; k < tr->omega_h; k++) {
		idx1 = tr->omega[k * tr->omega_l + 0];
		v4[0] = tr->O[idx1 * 4 + 0];
		v4[1] = tr->O[idx1 * 4 + 1];
		v4[2] = tr->O[idx1 * 4 + 2];
		v4[3] = tr->O[idx1 * 4 + 3];
		for (j = 0; j < 4; j++) {
			w4[j] = x0;
			for (i = 0; i < 4; i++) {
				a = z_mult_num(Z, v4[i], mtx[i * 4 + j]);
				w4[j] = z_add_num(Z, w4[j], a);
				}
			}
#ifdef DEBUG_MTX4_PERM_VERBOSE
		printf("w4 = ");
		printf("%2ld %2ld %2ld %2ld\n", w4[0], w4[1], w4[2], w4[3]);
#endif
		tr_vec4_normalize_num(tr, w4);
#ifdef DEBUG_MTX4_PERM_VERBOSE
		printf("%2ld %2ld %2ld %2ld\n", w4[0], w4[1], w4[2], w4[3]);
#endif
		idx2 = tr_search(tr, w4);
		if (idx2 < 0)
			return error("tr_mtx4_perm(): "
				"vector not found");
		k1 = tr->omega_idx[idx2];
		perm->m_ii(k, k1 + 1);
		}
	for (k = 0; k < tr->omega_h; k++) {
		if (perm->s_ii(k) == -1)
			return error("tr_mtx4_perm_omega(): "
				"not a permutation");
		}
#ifdef DEBUG_MTX4_PERM
	perm->println();
#endif
	return OK;
}

INT tr_generator_T(TITS_OVOID_REP *tr, 
	PERMUTATION_OP T_perm, INT f_verbose)
{
	INT T[64];

	z_mtx4_gen_T_num(tr->Z, T);
	
	if (f_verbose) {
		printf("T = \n");
		z_mtx4_print(tr->Z, T);
		}
	
	if (tr->f_has_omega)
		tr_mtx4_perm_omega(tr, T, T_perm);
	else
		tr_mtx4_perm(tr, T, T_perm);
	/* T_perm.print(); */
	
	return OK;
}

INT tr_generator_Mlambda(TITS_OVOID_REP *tr, 
	PERMUTATION_OP M_perm, INT f_verbose)
{
	INT lambda;
	INT M[64];

	lambda = z_prim_elt_num(tr->Z);
	z_mtx4_gen_Mlambda_num(tr->Z, M, lambda);
	
	if (f_verbose) {
		printf("M(%ld) = \n", lambda);
		z_mtx4_print(tr->Z, M);
		}
	
	if (tr->f_has_omega)
		tr_mtx4_perm_omega(tr, M, M_perm);
	else
		tr_mtx4_perm(tr, M, M_perm);
	/* M_perm.print(); */
	
	return OK;
}

INT tr_generators_S(TITS_OVOID_REP *tr, 
	VECTOR_OP V, INT f_verbose)
{
	INT a, b, i;
	INT S[64];
	PERMUTATION_OP S_perm;

	V->m_il(tr->nb_o - 1);
	for (i = 1; i < tr->nb_o; i++) {
		a = tr->O[i * 4 + 2];
		b = tr->O[i * 4 + 1];
		z_mtx4_gen_Sab_num(tr->Z, S, a, b);
		
		if (f_verbose) {
			printf("S(%ld,%ld) = \n", a, b);
			z_mtx4_print(tr->Z, S);
			}
	
		S_perm = (PERMUTATION_OP) V->s_i(i - 1);
		if (tr->f_has_omega)
			tr_mtx4_perm_omega(tr, S, S_perm);
		else
			tr_mtx4_perm(tr, S, S_perm);
		/* S_perm.print(); */
		}

	return OK;
}

INT Sz_q_generators(INT q, VECTOR_OP V)
{
	VECTOR_OB vp, ve;
	INT chi, deg, m;
	TITS_OVOID_REP *tr;
	INT nb_gen_S;
	PERMUTATION_OP p;
	SYM_OB go;

	factor_integer(q, &vp, &ve);
	if (vp.s_li() != 1)
		return error("Sz_q_generators(): "
			"q not a prime power");
	chi = vp.s_ii(0);
	if (chi != 2)
		return error("Sz_q_generators(): "
			"characteristic must be 2");
	deg = ve.s_ii(0);
	if (!ODD(deg))
		return error("Sz_q_generators(): "
			"degree over GF(2) must be odd");
	m = (deg - 1) >> 1;
	printf("Sz_q_generators(): chi = %ld deg = %ld m = %ld\n", chi, deg, m);
	fflush(stdout);

	tr = tr_open(deg, chi);
#if 0
	if (tr->f_has_omega)
		p_deg = tr->omega_h;
	else
		p_deg = tr->nb_o;
	p_kind = ik_permutation(p_deg);
	*pp_deg = p_deg;
	*pp_kind = p_kind;
#endif
	
	tr_generators_S(tr, V, TRUE);
	
	reduce_generators(V, &go, TRUE /* f_verbose */);
	printf("group order S(a,b):\n");
	go.println();
	nb_gen_S = V->s_li();
	printf("number of generators vor S: "
		"%ld\n", nb_gen_S);
	fflush(stdout);

	V->inc();
	p = (PERMUTATION_OP) V->s_i(nb_gen_S);
	tr_generator_Mlambda(tr, p, TRUE /* f_verbose */);
	
	V->inc();
	p = (PERMUTATION_OP) V->s_i(nb_gen_S + 1);
	tr_generator_T(tr, p, TRUE /* f_verbose */);
	
	printf("total number of generators for Sz(%ld): "
		"%ld\n", q, V->s_li());
	fflush(stdout);

	tr_close(tr);
	return OK;
}

INT gen_Sz(INT q)
{
	VECTOR_OB V;
	SYM_OB go;
	INT go_int;
	VECTOR_OB vp, ve;
	BYTE str[1024];
	LABRA_OB L;
	
	Sz_q_generators(q, &V);
	reduce_generators_labra(&V, &go, 
		TRUE /* f_verbose */, &L);
	printf("group order Sz(%ld) = ", q);
	go.print();
	fflush(stdout);
	
	go_int = ((INTEGER_OP) &go)->s_i();
	factor_integer(go_int, &vp, &ve);
	print_factorization(&vp, &ve, str);
	printf(" = %s\n", str);
	fflush(stdout);
	return OK;
}

#if 0
int main(int argc, char **argv)
{
	anfang();
	{

	/* gen_PSU(2); */
	/* gen_PSU(3); */
	/* gen_PSU(4); */
	/* gen_PSU(5); */
	/* gen_PSU(7); */
	/* gen_PSU(8); */
	/* gen_PSU(9); */
	gen_Sz(8);
	/* gen_Sz(32); */
	
	}
	ende();
	return 0;
}
#endif

#endif /* GFQ_TRUE */


