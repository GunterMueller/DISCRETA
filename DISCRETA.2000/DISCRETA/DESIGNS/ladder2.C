/* dc2.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef LADDER_TRUE

#include <DISCRETA/perm.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/lb.h>
#include <DISCRETA/cp.h>
#include <DISCRETA/divs.h>
#include <DISCRETA/ladder.h>

#define PRINT_TRANS_REPS
#define PRINT_STABILIZER

#define MAX_STEP 64

#if TEXDOCU
INT dc_get_transversals(DCY_OP dc, SYM_OP go, INT n, INT step, INT dc_no, 
	VECTOR_OP dc_rep, VECTOR_OP dc_idx, VECTOR_OP dc_trans_idx, 
	VECTOR_OP dc_omega, VECTOR_OP oiti, VECTOR_OP trans, VECTOR_OP trans_idx, INT f_v)
#endif
{
	// VECTOR_OB dc_rep, dc_idx, dc_trans_idx;
	// VECTOR_OB dc_omega;
	// VECTOR_OB oiti; // omega - image - transversal - index
	// VECTOR_OB trans, trans_idx;
	VECTOR_OB tr, ti, oi;
	VECTOR_OB TD;
	PERMUTATION_OP d, dm1, p;
	DCY_OP D, Dm1;
	INT i, j, ii, i0, s, a, prev, index, old_nb_dc;
	INT int_dc_idx;
	SYM_OB fusel;
	INT k, k1, omega0, omega, omega1;
	INT type = 0;
	void *data = NIL;

	dc_find_path(dc, go, step, dc_no, dc_rep, dc_idx, dc_trans_idx, f_v);
	
	if (f_v)
		printf("\ndc_get_transversals()\n");
	trans->m_il(step + 1);
	trans_idx->m_il(step + 1);
	dc_omega->m_il(step + 1);
	for (s = 0; s <= step; s++)
		dc_omega->m_ii(s, -1);
	oiti->m_il(step + 1);
	
	for (s = 1; s <= step; s++) {
		if (f_v) {
			printf("\n");
			printf("step %2ld:\n", s);
			printf("========\n");
			}
		a = dc_idx->s_ii(s);
		prev = dc_idx->s_ii(s - 1);
		dc_print_dc(dc, go, s, a);
		D = dc + s;
		Dm1 = dc + s - 1;
		index = D->s_T()->s_li();
		old_nb_dc = Dm1->s_D()->s_li();
		if (f_v) {
			printf("step %ld ", s);
			if (D->s_fDown_i())
				printf("(down) ");
			else
				printf("(up) ");
			printf("index = %ld old_nb_dc = %ld\n", index, old_nb_dc);
			printf("dc_rep step %ld = ", s);
			
				{
				INT sf = f_perm_print_start_with_zero;
				f_perm_print_start_with_zero = TRUE;
				
				dc_rep->s_i(s)->println();
				f_perm_print_start_with_zero = sf;
				}
			}
		if (D->s_fDown_i()) {
			if (f_v) {
				printf("on step %ld (down): "
					"previous of dc %ld is %ld\n", s, a, prev);
				}
			k = s / 2;
			k++;
			omega0 = n - k;
			if (f_v) {
				printf("omega0 = %ld\n", omega0); fflush(stdout);
				}
			d = (PERMUTATION_OP) dc_rep->s_i(s);
			omega = d->s_ii(omega0) - 1;
			if (f_v) {
				printf("omega = omega0^d = %ld\n", omega); fflush(stdout);
				}
			k = dc_trans_idx->s_ii(s - 1);
			p = (PERMUTATION_OP) D->s_T_i(k);
			dm1 = (PERMUTATION_OP) dc_rep->s_i(s - 1);
			omega1 = p->s_ii(omega0) - 1;
			omega1 = dm1->s_ii(omega1) - 1;
			if (f_v) {
				printf("omega1 = omega0^td_{s-1} = %ld\n", omega1); fflush(stdout);
				}
			if (omega1 != omega)
				return error("dc_get_transversals() omega1 != omega");
			dc_omega->m_ii(s, omega);
		
			tr.m_il(index);
			ti.m_il(index);
			oi.m_il(n);
			i0 = 0;
			for (ii = 0; ii < n; ii++)
				oi.m_ii(ii, -1);
			for (ii = 0; ii < index; ii++) {
				if (D->s_TDidx_iji(ii, prev) == a) {
					D->s_TDfusel_ij(ii, prev)->copy(tr.s_i(i0));
					ti.m_ii(i0, ii);
					p = (PERMUTATION_OP) tr.s_i(i0);
					k = p->s_ii(omega) - 1;
					p = (PERMUTATION_OP) D->s_T_i(ii);
					k1 = p->s_ii(omega0) - 1;
					k1 = dm1->s_ii(k1) - 1;
					if (k1 != k)
						return error("dc_get_transversals(): k != k1");
					oi.m_ii(k, i0);
					i0++;
					}
				}
			tr.realloc_z(i0);
			ti.realloc_z(i0);
			tr.copy((VECTOR_OP) trans->s_i(s));
			ti.copy((VECTOR_OP) trans_idx->s_i(s));
			oi.copy((VECTOR_OP) oiti->s_i(s));
			if (f_v) {
				printf("a subgroup of index %ld\n", i0);
				}
			}
		else { // upstep
			if (f_v) {
				printf("on step %ld (up): "
					"previous of dc %ld is %ld\n", s, a, prev);
				fflush(stdout);
				}
					
			tr.m_il(index);
			ti.m_il(index);
			i0 = 0;
			TD.m_il(index);
			vec_translate(D->s_T(), Dm1->s_D_i(prev), &TD, type, data);
			for (j = 0; j < index; j++) {
				if (f_v) {
					printf("j = %ld calling dc_trace_coset()\n", j);
					}
				if (dc_trace_coset(dc, TD.s_i(j), s - 1, 
					&int_dc_idx, &fusel, type, data) != OK)
					return error("dc_get_transversals(): error in dc_trace_coset");
				if (f_v) {
					printf("int_dc_idx = %ld\n", int_dc_idx);
					}
				if (int_dc_idx == prev) {
					fusel.copy(tr.s_i(i0));
					ti.m_ii(i0, j);
					i0++;
					}
				} // next j
					
			tr.realloc_z(i0);
			ti.realloc_z(i0);
			tr.copy((VECTOR_OP) trans->s_i(s));
			ti.copy((VECTOR_OP) trans_idx->s_i(s));
			if (f_v) {
				printf("of index %ld in an overgroup\n", i0);
				}
			} // else (upstep)
		
		if (f_v) {
			;
#ifdef PRINT_TRANS_REPS
		{
		VECTOR_OP p_tr, p_ti, p_oi;
		PERMUTATION_OP p;
		INT omega, image;
		INT sf = f_perm_print_start_with_zero;
		f_perm_print_start_with_zero = TRUE;
		
		p_tr = (VECTOR_OP) trans->s_i(s);
		p_ti = (VECTOR_OP) trans_idx->s_i(s);
		p_oi = (VECTOR_OP) oiti->s_i(s);
		omega = dc_omega->s_ii(s);
		if (omega >= 0)
			printf("omega = %ld\n", omega);
		printf("the transversal:\n");
		for (i = 0; i < p_tr->s_li(); i++) {
			printf("%ld (in coset %ld): ", i, p_ti->s_ii(i));
			if (omega >= 0) {
				p = (PERMUTATION_OP) p_tr->s_i(i);
				image = p->s_ii(omega) - 1;
				printf("omega image = %ld ", image);
				if (p_oi->s_ii(image) != i)
					return error("p_oi->s_ii(image) != i");
				}
			p_tr->s_i(i)->println();
			}
		f_perm_print_start_with_zero = sf;
		}
#endif
#ifdef PRINT_STABILIZER
		{
		INT dc_no = dc_idx->s_ii(s);
		VECTOR_OP A = D->s_Ad_i(dc_no);
		INT l;
		INT sf = f_perm_print_start_with_zero;
		f_perm_print_start_with_zero = TRUE;

		printf("stabilizer of dc %ld:\n", dc_no);
		l = A->s_li();
		for (ii = 0; ii < l; ii++) {
			printf("stab_gen %ld: ", ii);
			A->s_i(ii)->println();
			}
		f_perm_print_start_with_zero = sf;
		}
#endif
		} // if (f_v)
		}
	
	return OK;
}

#if TEXDOCU
INT dc_find_path(DCY_OP dc, SYM_OP go, INT step, INT dc_no, 
	VECTOR_OP dc_rep, VECTOR_OP dc_idx, VECTOR_OP dc_trans_idx, INT f_v)
#endif
{
	DCY_OP D, Dm1;
	INT i, j, s, a, b, index, old_nb_dc;
	SYM_OB fusel;

	if (f_v)
		printf("dc_find_path()\n");
	dc_idx->m_il(step + 1);
	dc_idx->m_ii(step, dc_no);
	
	dc_rep->m_il(step + 1);
	D = dc + step;
	D->s_D_i(dc_no)->copy(dc_rep->s_i(step));
	
	dc_trans_idx->m_il(step + 1);

	for (s = step; s > 0; s--) {
		if (f_v) {
			printf("\n");
			printf("step %2ld:\n", s);
			printf("========\n");
			}
		a = dc_idx->s_ii(s);
		if (f_v) {
			dc_print_dc(dc, go, s, a);
			}
		D = dc + s;
		Dm1 = dc + s - 1;
		index = D->s_T()->s_li();
		old_nb_dc = Dm1->s_D()->s_li();
		if (f_v) {
			printf("step %ld ", s);
			if (D->s_fDown_i())
				printf("(down) ");
			else
				printf("(up) ");
			printf("index = %ld old_nb_dc = %ld\n", index, old_nb_dc);
			}
		if (D->s_fDown_i()) {
			b = -1;
			for (j = 0; j < old_nb_dc; j++) {
				for (i = 0; i < index; i++) {
					b = D->s_TDidx_iji(i, j);
					if (b == a) {
						if (f_v) {
							printf("on step %ld (down): "
								"previous of dc %ld is %ld\n", s, a, j);
							}
						dc_idx->m_ii(s - 1, j);
						dc_trans_idx->m_ii(s - 1, i);
						Dm1->s_D_i(j)->copy(dc_rep->s_i(s - 1));
						break;
						}
					} // next i
				if (b == a)
					break;
				} // next j
			}
		else { // upstep
			for (i = 0; i < old_nb_dc; i++) {
				b = D->s_TDidx_iji(0, i);
				if (b == a) {
					if (f_v) {
						printf("on step %ld (up): "
							"previous of dc %ld is %ld\n", s, a, i);
						fflush(stdout);
						}
					dc_idx->m_ii(s - 1, i);
					dc_trans_idx->m_ii(s - 1, -1);
					Dm1->s_D_i(i)->copy(dc_rep->s_i(s - 1));
					break;
					}
				} // next i
			}
		}
	return OK;
	
}

#if TEXDOCU
INT dc2_print_k_set(VECTOR_OP R, SYM_OP stab_go, SYM_OP ol, INT alpha)
#endif
{
	SYM_OB go1;
	INT l, j;
	BYTE s1[256];
	BYTE s2[256];
	
	l = R->s_li();
	printf(" { ");
	for (j = 0; j < l; j++) {
		printf("%ld", R->s_ii(j));
		if (j < l - 1)
			printf(", ");
		}
	printf(" }");
	s1[0] = 0;
	s2[0] = 0;
	stab_go->sprint(s1);
	ol->sprint(s2);
	printf("_%s_%s alpha=%ld\n", s1, s2, alpha);
	return OK;
}

#if TEXDOCU
INT dc2_calc_intersection_indicator(DCY_OP dc, SYM_OP go, INT step, 
	PERMUTATION_OP d0, PERMUTATION_OP d1, INT k2, UNIPOLY_OP ii)
#endif
{
	DCY_OP Dc;
	INT i, l;
	VECTOR_OB R;
	SYM_OB ago, ol;
	INT alpha;
	VECTOR_OB V;
	SYM_OB tmp;

	V.m_il_n(k2 + 1);
	Dc = dc + step;
	l = Dc->s_D()->s_li();
	for (i = 0; i < l; i++) {
		dc2_get_orbit_length_and_intersection(dc, go, step, i, 
			d0, d1, k2, &R, &ago, &ol, &alpha);
		if (alpha > k2)
			return error("dc2_calc_intersection_indicator(): alpha > k2");
		// dc2_print_k_set(&R, &ago, &ol, alpha); fflush(stdout);
		V.s_i(alpha)->add(&ol, &tmp);
		tmp.swap(V.s_i(alpha));
		}
	ii->m_v(&V);
	ii->degree();
	return OK;
}

#if TEXDOCU
INT dc2_print_orbits(DCY_OP dc, SYM_OP go, INT step, 
	PERMUTATION_OP d0, PERMUTATION_OP d1, INT k2)
#endif
{
	DCY_OP Dc;
	INT i, l;
	VECTOR_OB R;
	SYM_OB ago, ol;
	INT alpha;

	Dc = dc + step;
	l = Dc->s_D()->s_li();
	for (i = 0; i < l; i++) {
		dc2_get_orbit_length_and_intersection(dc, go, step, i, 
			d0, d1, k2, &R, &ago, &ol, &alpha);
		dc2_print_k_set(&R, &ago, &ol, alpha);
		}
	return OK;
}

#if TEXDOCU
INT dc2_get_orbit_length_and_intersection(DCY_OP dc, SYM_OP go, INT step, INT dc_no, 
	PERMUTATION_OP d0, PERMUTATION_OP d1, INT k2, 
	VECTOR_OP R, SYM_OP ago, SYM_OP ol, INT *alpha)
#endif
// compare with kramer_mesner.C dc_print_dc()
{
	DCY_OP Dc;
	PERMUTATION_OP d;
	SYM_OB aut;
	LABRA_OB labra_A1;
	INT type;
	void *data;
	PERMUTATION_OB d1v, q;
	INT n, i, l, alpha_, a, b, k;

	type = 0;
	data = NIL;
	Dc = dc + step;
	// printf("step %ld dc no %ld: ", step, dc_no);
	d = (PERMUTATION_OP) Dc->s_D()->s_i(dc_no);
	reduce_generators_labra(Dc->s_Ad_i(dc_no), ago, 
		FALSE /* f_verbose */, &labra_A1);
	if (step == 0)
		k = 0;
	else if (step == 1)
		k = 1;
	else
		k = step / 2 + 1;
	d0->mult(d, &q);
	dc_get_k_set(&q, R, k, FALSE);
	go->ganzdiv(ago, ol);
	n = d1->s_li();
	d1->invers(&d1v);
	alpha_ = 0;
	l = R->s_li();
	for (i = 0; i < l; i++) {
		a = R->s_ii(i);
		b = d1v.s_ii(a) - 1;
		if (b >= n - k2)
			alpha_++;
		}
	*alpha = alpha_;
	return OK;
}

#if TEXDOCU
INT dc2_intersection(PERMUTATION_OP d, 
	PERMUTATION_OP d1, PERMUTATION_OP d2, INT k1, INT k2)
#endif
{
	PERMUTATION_OB q, d2v;
	VECTOR_OB R;
	INT i, l, n, alpha, a, b;
	
	d1->mult(d, &q);
	dc_get_k_set(&q, &R, k1, FALSE);
	n = d1->s_li();
	d2->invers(&d2v);
	alpha = 0;
	l = R.s_li();
	for (i = 0; i < l; i++) {
		a = R.s_ii(i);
		b = d2v.s_ii(a) - 1;
		if (b >= n - k2)
			alpha++;
		}
	return alpha;
}

#if TEXDOCU
INT dc2_calc_intersection_M(DCY_OP dc, SYM_OP go, INT k_max, MATRIX_OP M)
#endif
{
	INT up_to_step, nb_d;
	INT n, i, j, s1, s2, dcno1, dcno2;
	PERMUTATION_OP perm;
	UNIPOLY_OB ii;
	DCY_OB dc2[MAX_STEP];

	perm = (PERMUTATION_OP) dc->s_D_i(0);
	n = perm->s_li();
	if (k_max == 0)
		up_to_step = 0;
	else
		up_to_step = 2 * k_max - 1;
	nb_d = dc_calc_nb_d(dc, up_to_step);
	M->m_ilih_n(nb_d, nb_d);
	for (i = 0; i < nb_d; i++) {
		dc_dc_no_to_dc_idx(dc, i, &s1, &dcno1);
		for (j = 0; j < nb_d; j++) {
			dc_dc_no_to_dc_idx(dc, j, &s2, &dcno2);
			double_cosets2(dc, dc2, go, n, s1, dcno1, s2, dcno2, &ii);
			ii.swap(M->s_ij(i, j));
			}
		}
	return OK;
}

#if TEXDOCU
INT dc2_calc_intersection_M_tk(DCY_OP dc, SYM_OP go, INT t, INT k, MATRIX_OP M)
#endif
{
	DCY_OP Dc1, Dc2;
	INT s1, s2;
	INT nb_d1, nb_d2;
	DCY_OB dc2[MAX_STEP];
	INT i, j, n;
	PERMUTATION_OP perm;
	UNIPOLY_OB ii;

	if (k < t)
		return error("dc2_calc_intersection_M_tk() k < t");
	s1 = dc_k_to_step(t);
	s2 = dc_k_to_step(k);
	perm = (PERMUTATION_OP) dc->s_D_i(0);
	n = perm->s_li();
	Dc1 = dc + s1;
	Dc2 = dc + s2;
	nb_d1 = Dc1->s_D()->s_li();
	nb_d2 = Dc2->s_D()->s_li();
	M->m_ilih(nb_d2, nb_d1);
	for (i = 0; i < nb_d1; i++) {
		for (j = 0; j < nb_d2; j++) {
			double_cosets2(dc, dc2, go, n, s1, i, s2, j, &ii);
			// double_cosets_tg(dc, go, n, s1, i, s2, j, &ii);
			// intersections_by_enumeration(dc, go, n, s1, i, s2, j, &ii);
			printf("intersection indicator between %ld/%ld and %ld/%ld: \n", t, i, k, j);
			ii.println();
			ii.swap(M->s_ij(i, j));
			}
		}
	return OK;
}

#if TEXDOCU
INT intersections_by_enumeration(DCY_OP dc, SYM_OP go, INT n, 
	INT s1, INT dcno1, INT s2, INT dcno2, UNIPOLY_OP ii)
#endif
{
	VECTOR_OP A, U;
	LABRA_OB labra_A, labra_U;
	SYM_OB goA, goU;
	DCY_OP Dc2;
	INT k1, k2;
	SINGLE_COSET_WORK *scw;
	VECTOR_OB ii_vec;
	PERMUTATION_OB d1, d2, rep;
	INT alpha, no;
	
	// printf("intersection between %ld/%ld and %ld/%ld\n", 
	// 	s1, dcno1, s2, dcno2);
	ii_vec.m_il_n(n + 1);
	k1 = dc_step_to_k(s1);
	k2 = dc_step_to_k(s2);
	A = dc->s_Ad_i(0);
	// printf("got group A with %ld generators\n", A->s_li());
	reduce_generators_labra(A, &goA, FALSE /* f_verbose */, &labra_A);
	// printf("group order: ");
	// goA.println();

	Dc2 = dc + s2;
	U = Dc2->s_Ad_i(dcno2);
	// printf("got group U with %ld generators\n", U->s_li());
	reduce_generators_labra(U, &goU, FALSE /* f_verbose */, &labra_U);
	// printf("group order: ");
	// goU.println();
	
	Dc2 = dc + s1;
	Dc2->s_D_i(dcno1)->copy(&d1);
	Dc2 = dc + s2;
	Dc2->s_D_i(dcno2)->copy(&d2);

	scw = single_coset_open(&labra_A, &labra_U, FALSE);

	no = 0;
	single_coset_first(scw, &rep);
	do {
		no++;
		if (no % 100 == 0)
			printf("%ld ", no);
		if (no % 1000 == 0)
			printf("\n", no);
		fflush(stdout);
		alpha = dc2_intersection(&rep, &d1, &d2, k1, k2);
		// printf("(%ld) ", alpha);
		ii_vec.s_i(alpha)->inc();
	} while (single_coset_next(scw, &rep));
	printf("\n");
	
	single_coset_free(scw);
	ii->m_v(&ii_vec);
	return OK;
}

#if TEXDOCU
INT double_cosets2(DCY_OP dc, DCY_OP dc2, SYM_OP go, INT n, 
	INT s1, INT dcno1, INT s2, INT dcno2, UNIPOLY_OP ii)
#endif
{
	VECTOR_OP A;
	LABRA_OB labra_A;
	SYM_OB goA;
	SYM_OB id;
	DCY_OP Dc2;
	VECTOR_OB dc_rep, dc_idx, dc_trans_idx;
	VECTOR_OB dc_omega, oiti, trans, trans_idx;
	INT type;
	void *data;
	INT i, l, up_to_step, o, k1, k2;
	VECTOR_OP T, oiti_vec;
	PERMUTATION_OB d1, d2;
	
#if 0
	if (s2 < s1)
		return error("s2 < s1");
#endif
	type = DO_TYPE_SYM;
	data = NIL;
	printf("double cosets between %ld/%ld and %ld/%ld\n", 
		s1, dcno1, s2, dcno2);
	k1 = dc_step_to_k(s1);
	k2 = dc_step_to_k(s2);
	Dc2 = dc + s2;
	A = Dc2->s_Ad_i(dcno2);
	printf("got group A with %ld generators\n", A->s_li());
	reduce_generators_labra(A, &goA, FALSE /* f_verbose */, &labra_A);
	printf("group order: ");
	goA.println();

	dc_get_transversals(dc, go, n, s1, dcno1, 
		&dc_rep, &dc_idx, &dc_trans_idx, 
		&dc_omega, &oiti, &trans, &trans_idx, FALSE /* f_v */);

#if 0
	printf("dc_idx = ");
	dc_idx.println();
	printf("dc_omega = ");
	dc_omega.println();
	printf("oiti = ");
	oiti.Print();
#endif
	
	up_to_step = s1;
	do_copy(A->s_i(0), &id, type, data);
	do_one(&id, type, data);
	for (i = 0; i <= up_to_step; i++) {
		o = dc_omega.s_ii(i);
		oiti_vec = (VECTOR_OP) oiti.s_i(i);
		// printf("calling initialize_arbitrary step = %ld\n", i); fflush(stdout);
		T = (VECTOR_OP) trans.s_i(i);
		dc2[i].initialize_arbitrary(T, n /* deg */, i /* step */, o, oiti_vec, 
			0 /* type */, NIL /* data */);
#if 0
		dc2[i].initialize_Young(&T, n /* deg */, i /* step */, 
			0 /* type */, NIL /* data */);
#endif
		T->copy(dc2[i].s_T());
		// printf("index = %ld\n", T->s_li());
		if (i <= up_to_step)
			dc2[i].s_reduce_generators_mode()->m_i(2);
				/* 0: no reduction 
				 * 1: dimino 
				 * 2: labra */
		}
		
	// printf("initializing D\n"); fflush(stdout);
	dc2[0].s_D()->m_il(1);
	do_copy(&id, dc2[0].s_D_i(0), type, data);
	
	// printf("initializing Ad\n"); fflush(stdout);
	dc2[0].s_Ad()->m_il(1);
	((SYM_OP) A)->copy(dc2[0].s_Ad()->s_i(0));
	dc2[0].s_D()->println();
	dc2[0].s_Ad()->println();

	for (i = 1; i <= up_to_step; i++) {
		// printf("calling dc_do_step() step = %ld\n", i); fflush(stdout);
		dc_do_step(&dc2[0], i, &id, n /* deg */, 
			FALSE /* f_verbose */, type, data);
		}
	Dc2 = dc + s1;
	Dc2->s_D_i(dcno1)->copy(&d1);
	Dc2 = dc + s2;
	Dc2->s_D_i(dcno2)->copy(&d2);

	Dc2 = dc2 + s1;
	l = Dc2->s_D()->s_li();
	printf("double cosets2(): found %ld double cosets "
		"N_A(d_%ld,%ld) \\ A / N_A(d_%ld,%ld)\n", l, s1, dcno1, s2, dcno2);

	dc2_print_orbits(dc2, &goA, s1, &d1, &d2, k2);
#if 1
	{
	INT sf = f_perm_print_start_with_zero;
	f_perm_print_start_with_zero = TRUE;
	printf("d1 = ");
	d1.println();
	printf("d2 = ");
	d2.println();
	for (i = 0; i < l; i++) {
		printf("dc rep %ld: ", i);
		Dc2->s_D_i(i)->println();
		}
	f_perm_print_start_with_zero = sf;
	}
#endif
	dc2_calc_intersection_indicator(dc2, &goA, s1, &d1, &d2, k2, ii);
	printf("intersection indicator: ");
	ii->println();
	printf("\n");
#if 0
	for (i = 0; i < l; i++) {
		dc_print_dc_(dc2, &goA, s1, i, TRUE, &d1, &d2, s2);
		}
#endif
	
	
	return OK;
}

#endif /* LADDER_TRUE */

