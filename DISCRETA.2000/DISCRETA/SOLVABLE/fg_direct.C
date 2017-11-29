/* fg_direct.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SOLVABLE_TRUE

#include <DISCRETA/solvable.h>

#define MAX_NW 64

#if TEXDOCU
INT FG_OB::gen_is_maximal(INT i)
#endif
{
	SHORT nw[MAX_NW];
	INT j, l, P;
	ZE_OP ze;
	
	l = s_nb_ze_i();
#if 0
	gi = g_i(i);
	for (j = 0; j < l; j++) {
		if (j == i)
			continue;
		ze = s_ze_i(j);
		P = ze->s_P_i();
		if (P == gi)
			return FALSE;
		}
#endif
	for (j = 0; j < l; j++) {
		ze = s_ze_i(j);
		P = ze->s_P_i();
		int2nw(P, nw);
		if (nw[i] != 0)
			return FALSE;
		}
	return TRUE;
}

#if TEXDOCU
INT FG_OB::gen_is_central(INT i)
#endif
{
	INT j, l, gi, gj, a;
	ZE_OP ze;
	
	gi = g_i(i);
	l = s_nb_ze_i();
	ze = s_ze_i(i);
	for (j = 0; j < i; j++) {
		gj = g_i(j);
		a = ze->s_A_ii(j);
		if (a != gj)
			return FALSE;
		}
	for (j = i + 1; j < l; j++) {
		ze = s_ze_i(j);
		a = ze->s_A_ii(i);
		if (a != gi)
			return FALSE;
		}
	return TRUE;
}

#if TEXDOCU
INT FG_OB::gen_is_direct(INT i)
#endif
{
	SHORT nw[MAX_NW];
	INT j, k, l, P, a;
	ZE_OP ze;
	
	l = s_nb_ze_i();
	for (j = 0; j < l; j++) {
		ze = s_ze_i(j);
		for (k = 0; k < j; k++) {
			if (k == i)
				continue;
			a = ze->s_A_ii(k);
			int2nw(a, nw);
			if (nw[i] != 0)
				return FALSE;
			}
		}
	return TRUE;
}

#if TEXDOCU
INT FG_OB::adf_decomposition(FILE *fp_txt, INT f_v, INT f_vv)
#endif
{
	VECTOR_OB exp, Ord, Gidx, embedding, primes;
	INT l, r, n, i, e;
	FG_OB R, G1;

	n = s_n_i();
	abelian_direct_factors(&exp, &Ord, &Gidx, 
		&embedding, &primes, f_vv);
	l = exp.s_li();
	r = n;
	for (i = 0; i < l; i++) {
		e = exp.s_ii(i);
		r /= e;
		}
	if (f_v) {
		fprintf(fp_txt, "%ld x ", r);
		exp.fprint(fp_txt);
		}
	if (r == n) {
		if (f_v) {
			fprintf(fp_txt, " irreducible !\n");
			fflush(fp_txt);
			}
		}
	else {
		if (f_v) {
			fprintf(fp_txt, " not irreducible\n");
			fflush(fp_txt);
			}
		}
	residual_group(&R, &embedding, FALSE /* f_v */ );
	if (f_v) {
		R.fprint_gen(fp_txt);
		fflush(fp_txt);
		}
	R.complete_residual_group(fp_txt, &G1, &Ord, f_v);
	G1.swap(this);
	return OK;
}

#if TEXDOCU
INT FG_OB::complete_residual_group(FILE *fp_txt, FG_OP G, 
	VECTOR_OP Ord, INT f_v)
#endif
{
	INT Primes[MAX_NW];
	VECTOR_OP Ord_;
	ZE_OP ze0, ze;
	INT nb_ze, nb_ze0, p, P, l, ll, i, j, jj, a, b, k, gi;
	
	nb_ze0 = s_nb_ze_i();
	for (i = 0; i < nb_ze0; i++) {
		ze0 = s_ze_i(i);
		p = ze0->s_p_i();
		Primes[i] = p;
		}
	k = nb_ze0;
	l = Ord->s_li();
	for (i = 0; i < l; i++) {
		Ord_ = (VECTOR_OP) Ord->s_i(i);
		ll = Ord_->s_li();
		for (j = ll - 1; j >= 0; j--) {
			a = Ord_->s_ii(j);
			Primes[k++] = a;
			}
		}
	nb_ze = k;
#if 0
	if (f_v) {
		printf("FG::complete_residual_group(): ");
		for (i = 0; i < nb_ze; i++) {
			printf("%ld ", Primes[i]);
			}
		printf("\n");
		fflush(stdout);
		}
#endif
	G->init(nb_ze, Primes);
	for (i = 0; i < nb_ze0; i++) {
		ze0 = s_ze_i(i);
		ze = G->s_ze_i(i);
		P = ze0->s_P_i();
		ze->s_P()->m_i(P);
		for (j = 0; j < i; j++) {	
			a = ze0->s_A_ii(j);
			b = ze0->s_Av_ii(j);
			ze->s_A_i(j)->m_i(a);
			ze->s_Av_i(j)->m_i(b);
			}
		}
	k = nb_ze0;
	for (i = 0; i < l; i++) {
		Ord_ = (VECTOR_OP) Ord->s_i(i);
		ll = Ord_->s_li();
		for (j = 0; j < ll; j++) {
#if 0
			if (f_v) {
				printf("i = %ld j = %ld k = %ld\n", i, j, k);
				fflush(stdout);
				}
#endif
			ze = G->s_ze_i(k);
			if (j == 0)
				P = 0;
			else
				P = G->g_i(k - 1);
			ze->s_P()->m_i(P);
			for (jj = 0; jj < k; jj++) {
				gi = G->g_i(jj);
				ze->s_A_i(jj)->m_i(gi);
				ze->s_Av_i(jj)->m_i(gi);
				}
			k++;
			}
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::residual_group(FG_OP G, VECTOR_OP embedding, INT f_v)
#endif
{
	INT primes[MAX_NW];
	INT i, j, j0, k, k0, l, p, P0, P, a, a0;
	ZE_OP ze, ze0;

	l = embedding->s_li();
	for (i = 0; i < l; i++) {
		k0 = embedding->s_ii(i);
		ze0 = s_ze_i(k0);
		p = ze0->s_p_i();
		primes[i] = p;
		}
	G->init(l, primes);
	if (f_v) {
		printf("residual_group() residual group of order %ld\n", G->s_n_i());
		fflush(stdout);
		}
	for (i = 0; i < l; i++) {
		k0 = embedding->s_ii(i);
		ze0 = s_ze_i(k0);
		ze = G->s_ze_i(i);
		P0 = ze0->s_P_i();
		P = int_shrink(G, embedding, P0);
		ze->s_P()->m_i(P);
		if (f_v) {
			printf("residual_group() %c (originally %c): P0 = %ld becomes P = %ld\n", 
				(BYTE)('A' + i), (BYTE)('A' + k0), P0, P);
			fflush(stdout);
			}
		for (j = 0; j < i; j++) {
			j0 = embedding->s_ii(j);
			a0 = ze0->s_A_ii(j0);
			a = int_shrink(G, embedding, a0);
			ze->s_A_i(j)->m_i(a);
			if (f_v) {
				printf("residual_group() a0 = %ld becomes a = %ld\n", a0, a);
				fflush(stdout);
				}
			a0 = ze0->s_Av_ii(j0);
			a = int_shrink(G, embedding, a0);
			ze->s_Av_i(j)->m_i(a);
			}
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::complete_residue(VECTOR_OP embedding, 
	VECTOR_OP primes, VECTOR_OP Gidx, VECTOR_OP Ord)
#endif
{
	VECTOR_OP Gidx_, Ord_;
	INT i, j, l, ll, a, b, k;

	k = embedding->s_li();
	l = Gidx->s_li();
	for (i = 0; i < l; i++) {
		Gidx_ = (VECTOR_OP) Gidx->s_i(i);
		Ord_ = (VECTOR_OP) Ord->s_i(i);
		ll = Gidx_->s_li();
		for (j = ll - 1; j >= 0; j--) {
			a = Gidx_->s_ii(j);
			b = Ord_->s_ii(j);
			embedding->inc();
			embedding->m_ii(k, a);
			primes->inc();
			primes->m_ii(k, b);
			k++;
			}
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::gen_power_ord_and_index(INT i, VECTOR_OP Ord, VECTOR_OP Gidx)
#endif
{
	INT k, l, p, P, P_idx;
	ZE_OP ze;

	k = 0;
	Ord->m_il(0);
	Gidx->m_il(0);
	while (i >= 0) {
		ze = s_ze_i(i);
		p = ze->s_p_i();
		Ord->inc();
		Gidx->inc();
		Ord->m_ii(k, p);
		Gidx->m_ii(k, i);
		k++;
		P = ze->s_P_i();
		if (P == 0) {
			i = -1;
			continue;
			}
		P_idx = i2gen_idx(P);
		if (P_idx >= i)
			return error("gen_power_ord_and_index() P_idx >= i");
		i = P_idx;
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::is_adf_irreducible(FILE *fp_txt, INT f_v, INT f_vv)
#endif
{
	VECTOR_OB exp, Ord, Gidx, embedding, primes;
	INT l, r, n, i, e;

	n = s_n_i();
	abelian_direct_factors(&exp, &Ord, &Gidx, &embedding, &primes, f_vv);
	l = exp.s_li();
	r = n;
	for (i = 0; i < l; i++) {
		e = exp.s_ii(i);
		r /= e;
		}
	if (f_v) {
		fprintf(fp_txt, "%ld x ", r);
		exp.fprint(fp_txt);
		}
	if (r == n) {
		if (f_v) {
			fprintf(fp_txt, " irreducible !\n");
			fflush(fp_txt);
			}
		return TRUE;
		}
	if (f_v) {
		fprintf(fp_txt, " not irreducible\n");
		fflush(fp_txt);
		}
	return FALSE;
}

#if TEXDOCU
INT FG_OB::abelian_direct_factors(VECTOR_OP exp, VECTOR_OP Ord, VECTOR_OP Gidx, 
	VECTOR_OP embedding, VECTOR_OP primes, INT f_v)
#endif
{
	INT n, p, P, P_idx, i, l, r, gi, ord, idx, f_found;
	INT ii, jj, ll;
	INTEGER_OB k_ob, ord_ob;
	INT f_central[MAX_NW];
	INT f_direct[MAX_NW];
	INT f_residue[MAX_NW];
	ZE_OP ze;
	VECTOR_OB Ord_, Gidx_;
	
	n = s_n_i();
	l = s_nb_ze_i();
	for (i = 0; i < l; i++) {
		f_residue[i] = TRUE;
		}
	r = 0;
	k_ob.m_i(0);
	exp->m_il(0);
	Ord->m_il(0);
	Gidx->m_il(0);
	for (i = 0; i < l; i++) {
		if (f_v) {
			printf("gen %c: ", (BYTE)('A' + i));
			fflush(stdout);
			}
		if (!gen_is_central(i)) {
			if (f_v) {
				printf("is not central !\n");
				fflush(stdout);
				}
			f_central[i] = FALSE;
			continue;
			}
		f_central[i] = TRUE;
		if (!gen_is_direct(i)) {
			if (f_v) {
				printf("is not direct !\n");
				fflush(stdout);
				}
			f_direct[i] = FALSE;
			continue;
			}
		f_direct[i] = TRUE;
		gi = g_i(i);
		ze = s_ze_i(i);
		p = ze->s_p_i();
		P = ze->s_P_i();
		if (P) {
			P_idx = i2gen_idx(P);
			if (P_idx >= i)
				return error("abelian direct factors() P_idx >= i");
			if (!f_central[P_idx]) {
				if (f_v) {
					printf("a power is not central !\n");
					fflush(stdout);
					}
				f_central[i] = FALSE;
				continue;
				}
			if (!f_direct[P_idx]) {
				if (f_v) {
					printf("a power is not direct !\n");
					fflush(stdout);
					}
				f_direct[i] = FALSE;
				continue;
				}
			}
		if (!gen_is_maximal(i)) {
			if (f_v) {
				printf("is not maximal !\n");
				fflush(stdout);
				}
			continue;
			}
		gt_order(gi, &ord);
		ord_ob.m_i(ord);
		exp->search(k_ob.s_i(), FALSE, &ord_ob, &idx, &f_found);
		exp->insert_at(&k_ob, idx, &ord_ob);
		gen_power_ord_and_index(i, &Ord_, &Gidx_);
		k_ob.dec();
		Ord->insert_at(&k_ob, idx, &Ord_);
		k_ob.dec();
		Gidx->insert_at(&k_ob, idx, &Gidx_);
		ll = Gidx_.s_li();
		for (ii = 0; ii < ll; ii++) {
			jj = Gidx_.s_ii(ii);
			f_residue[jj] = FALSE;
			}
#if 0
		exp->inc();
		g_idx->inc();
		exp->m_ii(k, ord);
		g_idx->m_ii(k, i);
#endif
		if (f_v) {
			printf("order = %ld\n", ord);
			fflush(stdout);
			}
		}
	exp->realloc_z(k_ob.s_i());
	Ord->realloc_z(k_ob.s_i());
	Gidx->realloc_z(k_ob.s_i());
	r = 0;
	embedding->m_il(l);
	primes->m_il(l);
	for (i = 0; i < l; i++) {
		if (!f_residue[i])
			continue;
		embedding->m_ii(r, i);
		primes->m_ii(r, p);
		r++;
		}
	embedding->realloc_z(r);
	primes->realloc_z(r);
	return OK;
}

#endif /* SOLVABLE_TRUE */


