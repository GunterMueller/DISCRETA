/* intersection_aijk.C */

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


#if TEXDOCU
INT Intersection_M_via_Aijk(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	MATRIX_OP M, MATRIX_OP Aijk)
#endif
{
	INT nb_d, i, j;
	UNIPOLY_OB ii;

	nb_d = Aijk->s_li();
	M->m_ilih_n(nb_d, nb_d);
	for (i = 0; i < nb_d; i++) {
		for (j = 0; j < nb_d; j++) {
			Intersection_via_Aijk(go, stab_go, K, Aijk, i, j, &ii);
			ii.swap(M->s_ij(i, j));
			}
		}
	return OK;
}

#if TEXDOCU
INT Intersection_via_Aijk(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, MATRIX_OP Aijk, 
	INT dcno1, INT dcno2, UNIPOLY_OP ii)
#endif
{
	INT nb_dc, row;
	INT a, b, c, k, max_layer, size, oli, ol1;
	VECTOR_OB ii_vec;
	SYM_OB ol;
	SYM_OP ago;

	nb_dc = Aijk->s_li();
	row = dcno1 * nb_dc + dcno2;
	
	max_layer = K->s_ii(nb_dc - 1);
	ii_vec.m_il_n(max_layer + 1);
	
	ago = stab_go->s_i(dcno1);
	go->ganzdiv(ago, &ol);
	if (ol.s_obj_k() != INTEGER)
		return error("Intersection_via_Aijk() orbit length not an integer");
	ol1 = ((INTEGER_OP) &ol)->s_i();
	
	for (k = 0; k < nb_dc; k++) {
		a = Aijk->s_iji(row, k);
		size = K->s_ii(k);
		ago = stab_go->s_i(k);
		go->ganzdiv(ago, &ol);
		if (ol.s_obj_k() != INTEGER)
			return error("intersection_via_Aijk() orbit length not an integer");
		oli = ((INTEGER_OP) &ol)->s_i();
		a *= oli;
		b = a / ol1;
		if (b * ol1 != a)
			return error("Intersection_via_Aijk() b * ol1 != a");
		c = ii_vec.s_ii(size);
		c += b;
		ii_vec.m_ii(size, c);
		}
	ii->m_v(&ii_vec);
	ii->degree();
	return OK;
}

#if TEXDOCU
INT intersection_M_via_Aijk(DCY_OP dc, SYM_OP go, INT k_max, MATRIX_OP M, 
	MATRIX_OP Aijk, INT nb_d)
#endif
{
	INT i, j;
	UNIPOLY_OB ii;

	M->m_ilih_n(nb_d, nb_d);
	for (i = 0; i < nb_d; i++) {
		for (j = 0; j < nb_d; j++) {
			intersection_via_Aijk(dc, go, Aijk, nb_d, i, j, &ii);
			ii.swap(M->s_ij(i, j));
			}
		}
	return OK;
}

#if TEXDOCU
INT intersection_via_Aijk(DCY_OP dc, SYM_OP go, MATRIX_OP Aijk, INT nb_dc, 
	INT dcno1, INT dcno2, UNIPOLY_OP ii)
#endif
{
	DCY_OP Dc;
	VECTOR_OP Ad;
	INT row = dcno1 * nb_dc + dcno2;
	INT a, b, c, k, s, idx, max_layer, size, oli, ol1;
	VECTOR_OB ii_vec;
	SYM_OB ago, ol;
	LABRA_OB labra_A1;

	dc_dc_no_to_dc_idx(dc, nb_dc - 1, &s, &idx);
	max_layer = dc_step_to_k(s);
	ii_vec.m_il_n(max_layer + 1);
	
	dc_dc_no_to_dc_idx(dc, dcno1, &s, &idx);
	Dc = dc + s;
	Ad = Dc->s_Ad_i(idx);
	reduce_generators_labra(Ad, &ago, FALSE /* f_verbose */, &labra_A1);
	go->ganzdiv(&ago, &ol);
	if (ol.s_obj_k() != INTEGER)
		return error("intersection_via_Aijk() orbit length not an integer");
	ol1 = ((INTEGER_OP) &ol)->s_i();
	
	for (k = 0; k < nb_dc; k++) {
		a = Aijk->s_iji(row, k);
		dc_dc_no_to_dc_idx(dc, k, &s, &idx);
		size = dc_step_to_k(s);
		Dc = dc + s;
		Ad = Dc->s_Ad_i(idx);
		reduce_generators_labra(Ad, &ago, FALSE /* f_verbose */, &labra_A1);
		go->ganzdiv(&ago, &ol);
		if (ol.s_obj_k() != INTEGER)
			return error("intersection_via_Aijk() orbit length not an integer");
		oli = ((INTEGER_OP) &ol)->s_i();
		a *= oli;
		b = a / ol1;
		if (b * ol1 != a)
			return error("intersection_via_Aijk() b * ol1 != a");
		c = ii_vec.s_ii(size);
		c += b;
		ii_vec.m_ii(size, c);
		}
	ii->m_v(&ii_vec);
	ii->degree();
	return OK;
}

#if TEXDOCU
INT intersection_M_tk_via_Aijk(DCY_OP dc, SYM_OP go, INT t, INT k, 
	MATRIX_OP M, MATRIX_OP Aijk, INT nb_d)
#endif
{
	DCY_OP Dc1, Dc2;
	INT s1, s2;
	INT nb_d1, nb_d2;
	INT n, i, j, dcno1, dcno2;
	UNIPOLY_OB ii;
	PERMUTATION_OP perm;
	
	if (k < t)
		return error("intersection_M_tk_via_Aijk() k < t");
	if (nb_d * nb_d != Aijk->s_hi() || nb_d != Aijk->s_li())
		return error("intersection_M_tk_via_Aijk() Aijk has wrong size");
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
		dc_dc_idx_to_dc_no(dc, s1, i, &dcno1);
		for (j = 0; j < nb_d2; j++) {
			dc_dc_idx_to_dc_no(dc, s2, j, &dcno2);
			intersection_via_Aijk(dc, go, Aijk, nb_d, dcno1, dcno2, &ii);
			
			printf("intersection indicator between %ld/%ld(%ld) and %ld/%ld(%ld): \n", t, i, dcno1, k, j, dcno2);
			ii.println();
			ii.swap(M->s_ij(i, j));
			}
		}
	return OK;
}

#if TEXDOCU
INT calc_Aijk(MATRIX_OP Ainf, MATRIX_OP Aijk)
#else
// Ainf is an upper triangular matrix
#endif
{
	MATRIX_OB Ainf_inv;
	MATRIX_OB A, B;
	INT m, n, i, j, k, a, b, c;
	
	m = Ainf->s_hi();
	n = Ainf->s_li();
	if (m != n)
		return error("calc_Aijk() Ainf not quadratic !");
	Ainf->invers(&Ainf_inv);
	A.m_ilih(1, m);
	Aijk->m_ilih(n, n * n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < m; k++) {
				a = Ainf->s_iji(k, i);
				b = Ainf->s_iji(k, j);
				c = a * b;
				A.m_iji(k, 0, c);
				}
			Ainf_inv.mult(&A, &B);
			for (k = 0; k < m; k++) {
				a = B.s_iji(k, 0);
				Aijk->m_iji(i * n + j, k, a);
				}
			}
		}
	return OK;
}

#endif /* LADDER_TRUE */

