/* gfq_psu.C 
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
generators for $\PSU_3(q^2)$
according to Huppert, Endliche Gruppen I~\cite{Huppert67}, p. 244.
#endif

void z_mtx3_gen_Qab_num(ZECH_DATA *Z, INT *mtx, INT a, INT b)
{
	INT i, j, *p, a1;

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			p = mtx + i * 3 + j;
			if (i > j)
				*p = z_zero_num(Z);
			else if (i == j)
				*p = z_one_num(Z);
			}
		}
#ifdef DEBUG_Z_MTX3_GEN_QAB_NUM
	printf("a = ");
	z_print_elt_num_verbose(Z, a);
	fflush(stdout);
#endif
	a1 = z_apply_frob2_num(Z, a);
#ifdef DEBUG_Z_MTX3_GEN_QAB_NUM
	printf("a1 = F(a) = ");
	z_print_elt_num_verbose(Z, a1);
	fflush(stdout);
#endif
	a1 = z_negate_num(Z, a1);
#ifdef DEBUG_Z_MTX3_GEN_QAB_NUM
	printf("a1' = -F(a) = ");
	z_print_elt_num_verbose(Z, a1);
	fflush(stdout);
	printf("b = ");
	z_print_elt_num_verbose(Z, b);
	fflush(stdout);
#endif
	mtx[0 * 3 + 1] = a;
	mtx[1 * 3 + 2] = a1;
	mtx[0 * 3 + 2] = b;
}

void z_mtx3_gen_Hk_num(ZECH_DATA *Z, INT *mtx, INT k)
{
	INT i, j, *p, k_inv, k_tau, a, b;

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			p = mtx + i * 3 + j;
			if (i != j)
				*p = z_zero_num(Z);
			}
		}
#ifdef DEBUG_Z_MTX3_GEN_HK_NUM
	printf("k = ");
	z_print_elt_num_verbose(Z, k);
	fflush(stdout);
#endif
	k_inv = z_inverse_num(Z, k);
#ifdef DEBUG_Z_MTX3_GEN_HK_NUM
	printf("k_inv = ");
	z_print_elt_num_verbose(Z, k_inv);
	fflush(stdout);
#endif
	k_tau = z_apply_frob2_num(Z, k);
#ifdef DEBUG_Z_MTX3_GEN_HK_NUM
	printf("k_tau = ");
	z_print_elt_num_verbose(Z, k_tau);
	fflush(stdout);
#endif
	a = z_apply_frob2_num(Z, k_inv);
#ifdef DEBUG_Z_MTX3_GEN_HK_NUM
	printf("a = ");
	z_print_elt_num_verbose(Z, a);
	fflush(stdout);
#endif
	b = z_mult_num(Z, k_tau, k_inv);
#ifdef DEBUG_Z_MTX3_GEN_HK_NUM
	printf("b = ");
	z_print_elt_num_verbose(Z, b);
	fflush(stdout);
#endif
	mtx[0 * 3 + 0] = a;
	mtx[1 * 3 + 1] = b;
	mtx[2 * 3 + 2] = k;
#ifdef DEBUG_Z_MTX3_GEN_HK_NUM
	z_mtx3_print(Z, mtx);
	fflush(stdout);
#endif
}

void z_mtx3_gen_T_num(ZECH_DATA *Z, INT *mtx)
{
	INT i, j;

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			mtx[i * 3 + j] = z_zero_num(Z);
			}
		}
	mtx[0 * 3 + 2] = z_one_num(Z);
	mtx[1 * 3 + 1] = z_mone_num(Z);
	mtx[2 * 3 + 0] = z_one_num(Z);
}

ISOTROPIC_REP *ir_open(INT deg, INT chi)
{
	ISOTROPIC_REP *ir;

	ir = (ISOTROPIC_REP *) my_malloc(sizeof(ISOTROPIC_REP), "gfq_psu.C: ir_open");
	ir->Z = zech_open(deg, chi, TRUE);
	
	ir->q1 = i_power_j(chi, (deg >> 1));
	ir->IV = (INT *) my_malloc(ir->Z->q * ir->Z->q * 3 * sizeof(INT), "gfq_psu.C: ir_open()");
	ir->nb_iv = 0;
	ir->IV[ir->nb_iv * 3 + 0] = 0;
	ir->IV[ir->nb_iv * 3 + 1] = 0;
	ir->IV[ir->nb_iv * 3 + 2] = 1;
	ir->nb_iv++;
	printf("added isotropic vector (0,0,1)\n");
	z_collect_tr_plus_nr_zero(ir->Z, ir->IV, &ir->nb_iv);
	printf("isotropic vectors: %ld\n", ir->nb_iv);

#if 0
	INT k;
	for (k = 0; k < ir->nb_iv; k++) {
		printf("%3ld: %2ld %2ld %2ld\n", k, 
			ir->IV[k * 3 + 0], 
			ir->IV[k * 3 + 1], 
			ir->IV[k * 3 + 2]);
		}
#endif
	
	return ir;
}

void ir_close(ISOTROPIC_REP *ir)
{
	if (ir->IV)
		my_free(ir->IV);
	ir->nb_iv = 0;

	zech_free(ir->Z);
	my_free(ir);
}

INT ir_cmp_vec(INT *v1, INT *v2)
/* the projective elements have normalized 
 * coordinates;
 * we assume that we have numeric form of the elements 
 * so that the lexicographic sorting reflects the 
 * sorting according to the p-adic representation */
{
	INT i, a, b;

	for (i = 0; i < 3; i++) {
		a = v1[i];
		b = v2[i];
		if (a < b)
			return -1;
		if (a > b)
			return 1;
		}
	return 0;
}

INT ir_search(ISOTROPIC_REP *ir, INT *x)
{
	INT i, j, mid, res;
	
	i = 0;
	j = ir->nb_iv - 1;
	while (i <= j) {
		mid = i + (j - i) / 2;
		res = ir_cmp_vec(ir->IV + mid * 3, x);
		if (res == 0)
			return mid;
		else if (res < 0)
			i = mid + 1;
		else
			j = mid - 1;
		}
	return -1;
}


INT z_collect_tr_plus_nr_zero(ZECH_DATA *Z, 
	INT *T, INT *nb_sol)
/* computes a list of all isotropic vectors of P_2(K) 
 * where K = GF(q^2). 
 * The ground-field is that which z_trace2_num and z_norm2_num use:
 * the unique subfield of index 2 (= GF(q) );
 * The K-elements are held in numeric (not zech log) form. 
 * The vector is sorted. */
{
	INT a, b, a1, b1, trace_b, norm_a, c;

	for (a = 0; a < Z->q; a++) {
		a1 = Z->Num[a];
		norm_a = z_norm2_zlog(Z, a1);
		for (b = 0; b < Z->q; b++) {
			b1 = Z->Num[b];
			trace_b = z_trace2_zlog(Z, b1);
			c = z_add_zlog(Z, norm_a, trace_b);
			if (c == Z->idx_zero) {
				T[*nb_sol * 3 + 0] = 1;
				T[*nb_sol * 3 + 1] = a;
				T[*nb_sol * 3 + 2] = b;
				(*nb_sol)++;
				}
			}
		}
	return OK;
}

INT ir_vec3_check_if_normalized_num(ISOTROPIC_REP *ir, INT *v)
{
	if (!ir_vec3_is_normalized_num(ir, v)) {
		return error("ir_check_if_normalized3(): "
			"vector not normalized");
		}
	return OK;
}

INT ir_vec3_is_normalized_num(ISOTROPIC_REP *ir, INT *v)
{
	ZECH_DATA *Z = ir->Z;
	INT a;

	a = v[0];
	if (!z_is_zero_num(Z, a)) {
		if (!z_is_one_num(Z, a)) {
			return FALSE;
			}
		return TRUE;
		}
	a = v[1];
	if (!z_is_zero_num(Z, v[1])) {
		if (!z_is_one_num(Z, a)) {
			return FALSE;
			}
		return TRUE;
		}
	a = v[2];
	if (!z_is_zero_num(Z, v[2])) {
		if (!z_is_one_num(Z, a)) {
			return FALSE;
			}
		return TRUE;
		}
	return error("ir_vec3_is_normalized_num(): zero vector");
}

void ir_vec3_normalize_num(ISOTROPIC_REP *ir, INT *v)
{
	ZECH_DATA *Z = ir->Z;
	INT a;

	if (!z_is_zero_num(Z, v[0])) {
		a = z_inverse_num(Z, v[0]);
		v[0] = z_mult_num(Z, a, v[0]);
		v[1] = z_mult_num(Z, a, v[1]);
		v[2] = z_mult_num(Z, a, v[2]);
		return;
		}
	if (!z_is_zero_num(Z, v[1])) {
		a = z_inverse_num(Z, v[1]);
		v[1] = z_mult_num(Z, a, v[1]);
		v[2] = z_mult_num(Z, a, v[2]);
		return;
		}
	if (!z_is_zero_num(Z, v[2])) {
		a = z_inverse_num(Z, v[2]);
		v[2] = z_mult_num(Z, a, v[2]);
		return;
		}
	error("ir_vec3_normalize_num(): zero vector");
}

void ir_mtx3_normalize_num(ISOTROPIC_REP *ir, INT *v)
{
	ZECH_DATA *Z = ir->Z;
	INT a, i, j;

	for (i = 0; i < 3; i++) {
		a = v[i];
		if (!z_is_zero_num(Z, a)) {
			a = z_inverse_num(Z, a);
			for (j = i; j < 9; j++) {
				v[j] = z_mult_num(Z, a, v[j]);
				}
			return;
			}
		}
}

INT ir_mtx3_perm(ISOTROPIC_REP *ir, 
	INT *mtx, PERMUTATION_OP perm)
{
	ZECH_DATA *Z = ir->Z;
	INT i, j, k, k1;
	INT v3[3], w3[3], x0, a;

	perm->m_il(ir->nb_iv);
	x0 = z_zero_num(Z);
	for (k = 0; k < ir->nb_iv; k++) {
		v3[0] = ir->IV[k * 3 + 0];
		v3[1] = ir->IV[k * 3 + 1];
		v3[2] = ir->IV[k * 3 + 2];
		for (j = 0; j < 3; j++) {
			w3[j] = x0;
			for (i = 0; i < 3; i++) {
				a = z_mult_num(Z, v3[i], mtx[i * 3 + j]);
				w3[j] = z_add_num(Z, w3[j], a);
				}
			}
#ifdef DEBUG_MTX3_PERM_VERBOSE
		printf("w3 = ");
		printf("%2ld %2ld %2ld\n", w3[0], w3[1], w3[2]);
#endif
		ir_vec3_normalize_num(ir, w3);
#ifdef DEBUG_MTX3_PERM_VERBOSE
		printf("%2ld %2ld %2ld\n", w3[0], w3[1], w3[2]);
#endif
		k1 = ir_search(ir, w3);
		if (k1 < 0)
			return error("ir_mtx3_perm(): "
				"vector not found");
		perm->m_ii(k, k1 + 1);
		}
#ifdef DEBUG_MTX3_PERM
	perm->print();
#endif
	return OK;
}

INT ir_generator_T(ISOTROPIC_REP *ir, 
	PERMUTATION_OP T_perm, INT f_verbose)
{
	INT T[64];

	z_mtx3_gen_T_num(ir->Z, T);
	
	if (f_verbose) {
		printf("T = \n");
		z_mtx3_print(ir->Z, T);
		}
	
	ir_mtx3_perm(ir, T, T_perm);
	/* H_perm.print(); */
	
	return OK;
}

INT ir_generator_Hk(ISOTROPIC_REP *ir, 
	PERMUTATION_OP H_perm, INT f_verbose)
{
	INT k;
	INT H[64];

	k = z_prim_elt_num(ir->Z);
	z_mtx3_gen_Hk_num(ir->Z, H, k);
	
	if (f_verbose) {
		printf("H(%ld) = \n", k);
		z_mtx3_print(ir->Z, H);
		}
	
	ir_mtx3_perm(ir, H, H_perm);
	/* H_perm.print(); */
	
	return OK;
}

INT ir_generators_Q(ISOTROPIC_REP *ir, 
	VECTOR_OP V, INT f_verbose)
{
	INT a, b, i;
	INT Q[64];
	PERMUTATION_OP Q_perm;

	V->m_il(ir->nb_iv - 1);
	for (i = 1; i < ir->nb_iv; i++) {
		a = ir->IV[i * 3 + 1];
		b = ir->IV[i * 3 + 2];
		z_mtx3_gen_Qab_num(ir->Z, Q, a, b);
		
		if (f_verbose) {
			printf("Q(%ld,%ld) = \n", a, b);
			z_mtx3_print(ir->Z, Q);
			}
	
		Q_perm = (PERMUTATION_OP) V->s_i(i - 1);
		ir_mtx3_perm(ir, Q, Q_perm);
		/* Q_perm.print(); */
		}

	return OK;
}

INT PSU_3_q2_generators(INT q, VECTOR_OP V)
{
	VECTOR_OB vp, ve;
	INT chi, deg1, deg2;
	ISOTROPIC_REP *ir;
	INT nb_gen_Q, nb_gen_PSU;
	PERMUTATION_OP p;
	SYM_OB go;

	factor_integer(q, &vp, &ve);
	if (vp.s_li() != 1)
		return error("PSU_3_q2_generators(): "
			"q not a prime power");
	chi = vp.s_ii(0);
	deg1 = ve.s_ii(0);
	deg2 = deg1 << 1;
	printf("PSU_3_q2_generators(): "
		"chi = %ld deg1 = %ld\n", chi, deg1);
	fflush(stdout);

	ir = ir_open(deg2, chi);
	/* p_deg = ir->nb_iv;
	p_kind = ik_permutation(p_deg); */
	
	ir_generators_Q(ir, V, FALSE);
	
	reduce_generators(V, &go, TRUE /* f_verbose */);
	printf("group order Q(a,b):\n");
	go.println();
	nb_gen_Q = V->s_li();
	printf("number of generators vor Q: %ld\n", nb_gen_Q);
	fflush(stdout);

	V->inc();
	p = (PERMUTATION_OP) V->s_i(nb_gen_Q);
	ir_generator_Hk(ir, p, TRUE /* f_verbose */);
	
	V->inc();
	p = (PERMUTATION_OP) V->s_i(nb_gen_Q + 1);
	ir_generator_T(ir, p, TRUE /* f_verbose */);
	
	printf("total number of generators for PSU(3,%ld^2): "
		"%ld\n", ir->q1, V->s_li());
	fflush(stdout);

	nb_gen_PSU = V->s_li();
	ir_close(ir);
	return OK;
}


INT gen_PSU(INT q)
{
	VECTOR_OB V;
	SYM_OB go_PSU;
	INT go_PSU_int;
	VECTOR_OB vp, ve;
	BYTE str[1024];
	LABRA_OB L;
	
	PSU_3_q2_generators(q, &V);
	reduce_generators_labra(&V, &go_PSU, TRUE /* f_verbose */, &L);
	printf("group order PSU_3(%ld^2) = ", q);
	go_PSU.print();
	fflush(stdout);
	
	go_PSU_int = ((INTEGER_OP) &go_PSU)->s_i();
	factor_integer(go_PSU_int, &vp, &ve);
	print_factorization(&vp, &ve, str);
	printf(" = %s\n", str);
	fflush(stdout);
	return OK;
}

#endif /* GFQ_TRUE */


