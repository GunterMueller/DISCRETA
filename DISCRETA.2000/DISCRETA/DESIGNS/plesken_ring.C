/* plesken_ring.C */

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



/*
 * Plesken ring structure:
 */

#if TEXDOCU
INT plesken_product_block_wise(MATRIX_OP Ainf_block_wise, MATRIX_OP coeff, 
	VECTOR_OP orbit_length, VECTOR_OP args, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, INT k_min, 
	VECTOR_OP type, INT f_v, INT f_vv)
#else
Computes the plesken product of $G$-orbits with index $i$, $i \in args$.
Result is type. 
The result is collected layer-wise. 
First comes highest\_layer, altogether num\_layers downwards. 
#endif
{
	INT i, j, k;
	VECTOR_OB aiaj, aijk;
	SYM_OB c, tmp;

	if (f_vv) {
		printf("plesken_product() args");
		args->println();
		}
	k = K_first->s_li() - 1;
	if (highest_layer - (num_layers - 1) != k_min)
		return error("plesken_product_block_wise() highest_layer - (num_layers - 1) != k_min");
	if (highest_layer > k)
		return error("plesken_product_block_wise() highest_layer > k");
	plesken_multiply_columns(Ainf_block_wise, 
		orbit_length, args, 
		K_first, K_len, 
		k_min, highest_layer, k, 
		&aiaj, f_v, f_vv);
	if (f_vv) {
		printf("aiaj = ");
		aiaj.Print();
		printf("\n");
		}
	type->m_il(num_layers);
	for (j = 0; j < num_layers; j++) {
		i = highest_layer - j;
		plesken_decompose_layer(Ainf_block_wise, coeff, 
			orbit_length, 
			K_first, K_len, 
			k_min, highest_layer, k, i, type->s_i(j), 
			&aiaj, f_v, f_vv);
		}
	if (f_v) {
		printf("type of ");
		args->print();
		printf(" = ");
		type->print();
		printf("\n");
		}

	return OK;
}

#if TEXDOCU
INT plesken_decompose_layer(MATRIX_OP Ainf_block_wise, MATRIX_OP coeff, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT k_min, INT k_max, INT k, INT i, SYM_OP ni, 
	VECTOR_OP aiaj, INT f_v, INT f_vv)
#endif
{
	INT first, last, len, first_k_orbit;
	INT fi, li, ii, j, jj, fj, lj;
	SYM_OB tmp1, tmp2, tmp3, sum;
	SYM_OP pa, pb, pc;
	MATRIX_OP pM;
	
	if (i < k_min || i > k_max)
		return error("plesken_decompose_layer() i < k_min || i > k_max");
	first = K_first->s_ii(k_min);
	last = K_first->s_ii(k_max) + K_len->s_ii(k_max);
	len = last - first;
	first_k_orbit = K_first->s_ii(k);
	fi = K_first->s_ii(i);
	li = K_len->s_ii(i);
	ni->freeself();
	ni->m_i_i(0);
	for (ii = 0; ii < li; ii++) {
		// multiply with the identity matrix M_{i,i}:
		aiaj->s_i(fi + ii - first)->copy(&sum);
		// multiply with M_{i,j}:
		for (j = i + 1; j <= k_max; j++) {
			fj = K_first->s_ii(j);
			lj = K_len->s_ii(j);
			pM = (MATRIX_OP) Ainf_block_wise->s_ij(i - k_min, j - k_min);
			if (pM->s_hi() != li)
				return error("plesken_decompose_layer() pM->s_hi() != li");
			if (pM->s_li() != lj)
				return error("plesken_decompose_layer() pM->s_li() != lj");
			pc = coeff->s_ij(i, j);
			for (jj = 0; jj < lj; jj++) {
				pa = pM->s_ij(ii, jj);
				pb = aiaj->s_i(fj + jj - first);
				tmp1.freeself();
				tmp2.freeself();
				tmp3.freeself();
				pa->mult(pb, &tmp1);
				tmp1.mult(pc, &tmp2);
				tmp2.add(&sum, &tmp3);
				tmp3.swap(&sum);
				}
			}
		pa = orbit_length->s_i(fi + ii);
		tmp1.freeself();
		tmp2.freeself();
		sum.mult(pa, &tmp1);
		tmp1.add(ni, &tmp2);
		tmp2.swap(ni);
		}
	
	return OK;
}

#if TEXDOCU
INT plesken_multiply_columns(MATRIX_OP Ainf_block_wise, 
	VECTOR_OP orbit_length, VECTOR_OP args, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT k_min, INT k_max, INT k, 
	VECTOR_OP aiaj, INT f_v, INT f_vv)
#endif
{
	INT first, last, len, i, ii, j, l, m, idx, num_args, first_k_orbit;
	INT f_is_constant;
	MATRIX_OP pM;
	SYM_OB tmp1, tmp2, tmp3;
	SYM_OP pa;
	
	first = K_first->s_ii(k_min);
	last = K_first->s_ii(k_max) + K_len->s_ii(k_max);
	len = last - first;
	first_k_orbit = K_first->s_ii(k);
	if (aiaj->s_obj_k() == VECTOR && aiaj->s_li() == len)
		;
	else
		aiaj->m_il(len);
	num_args = args->s_li();
	if (num_args < 1)
		return error("plesken_multiply_columns() num_args < 1");
	f_is_constant = TRUE;
	for (i = 1; i < num_args; i++) {
		if (args->s_ii(i) != args->s_ii(i - 1)) {
			f_is_constant = FALSE;
			break;
			}
		}
	ii = 0;
	for (i = k_min; i <= k_max; i++) {
		l = K_len->s_ii(i);
		pM = (MATRIX_OP) Ainf_block_wise->s_ij(i - k_min, k - k_min);
		if (pM->s_hi() != l)
			return error("plesken_multiply_columns() pM->s_hi() != l");
		for (j = 0; j < l; j++) {
			tmp1.freeself();
			tmp2.freeself();
			for (m = 0; m < num_args; m++) {
				idx = args->s_ii(m) - first_k_orbit;
				pa = pM->s_ij(j, idx);
				if (m == 0) 
					pa->copy(&tmp1);
				else {
					tmp2.freeself();
					pa->mult(&tmp1, &tmp2);
					tmp2.swap(&tmp1);
					}
				}
			if (f_is_constant) {
				tmp2.freeself();
				tmp3.freeself();
				idx = args->s_ii(0) - first_k_orbit;
				pa = pM->s_ij(j, idx);
				pa->addinvers(&tmp2);
				tmp1.add(&tmp2, &tmp3);
				tmp3.swap(&tmp1);
				}
			tmp1.swap(aiaj->s_i(ii));
			ii++;
			}
		}
	if (ii != len)
		return error("plesken_multiply_columns() ii != len");
	return OK;
}
	
#if TEXDOCU
INT plesken_product(MATRIX_OP Ainf, MATRIX_OP coeff, 
	VECTOR_OP orbit_length, VECTOR_OP args, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, INT k_min, 
	VECTOR_OP type, INT f_v, INT f_vv)
#else
Computes the plesken product of $G$-orbits with index $i$, $i \in args$.
Result is type. 
The result is collected layer-wise. 
First comes highest\_layer, altogether num\_layers downwards. 
#endif
{
	INT nb_dc, i, j, f, l, ll, ii, idx, d0;
	VECTOR_OB aiaj, aijk;
	SYM_OP pa;
	SYM_OB c, tmp;

	if (f_v) {
		printf("plesken_product() args");
		args->println();
		}
	d0 = K_first->s_ii(k_min);
	nb_dc = Ainf->s_li() + d0;
	if (highest_layer - (num_layers - 1) < k_min)
		return error("plesken_product() highest_layer - (num_layers - 1) < k_min");
	f = K_first->s_ii(highest_layer - num_layers + 1);
	i = K_first->s_ii(highest_layer) + K_len->s_ii(highest_layer);
	l = i - f;
	ll = nb_dc - f;
	aiaj.m_il_n(nb_dc);
	for (i = 0; i < ll; i++) {
		ii = f + i;
		c.freeself();
		for (j = 0; j < args->s_li(); j++) {
			idx = args->s_ii(j) - d0;
			pa = Ainf->s_ij(ii - d0, idx);
			if (j == 0) {
				pa->copy(&c);
				}
			else {
				tmp.freeself();
				pa->mult(&c, &tmp);
				tmp.swap(&c);
				}
			} // next j
		c.swap(aiaj.s_i(ii));
		}
	if (f_vv) {
		printf("aiaj[%ld,+%ld] = ", f, ll);
		for (i = 0; i < ll; i++) {
			ii = f + i;
			aiaj.s_i(ii)->print();
			if (i < ll - 1)
				printf(", ");
			}
		printf("\n");
		}
	aijk.m_il_n(nb_dc);
	
	plesken_decompose(Ainf, coeff, K_first, K_len, 
		&aiaj, &aijk, f, l, d0, f_v);
	
	if (f_v) {
		printf("aijk = ");
		aijk.println();
		}

	plesken_decomposition_collect_layer_wise(orbit_length, 
		&aijk, K_first, K_len, 
		highest_layer, num_layers, type, f_v, f_vv);

	return OK;
}

#if TEXDOCU
INT plesken_decompose(MATRIX_OP Ainf, MATRIX_OP coeff, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	VECTOR_OP aiaj, VECTOR_OP aijk, INT first, INT len, INT d0, INT f_v)
#else
aijk is already allocated.
#endif
{
	INT i, ii, j, nb_dc, layer_i, layer_j;
	SYM_OP a, b, c;
	SYM_OB d, tmp, tmp1;

	nb_dc = Ainf->s_hi() + d0;
	for (i = 0; i < len; i++) {
		ii = first + i;
		d.freeself();
		for (j = ii; j < nb_dc; j++) {
			a = Ainf->s_ij(ii - d0, j - d0);
			b = aiaj->s_i(j);
			layer_i = dc_plesken_layer_idx(K_first, K_len, ii);
			layer_j = dc_plesken_layer_idx(K_first, K_len, j);
			c = coeff->s_ij(layer_i, layer_j);
			if (j == ii) {
				a->mult(b, &tmp1);
				tmp1.mult(c, &d);
				}
			else {
				tmp.freeself();
				tmp1.freeself();
				a->mult(b, &tmp1);
				tmp1.mult(c, &tmp);
				d.add(&tmp, &tmp1);
				tmp1.swap(&d);
				}
			} // next j
		d.swap(aijk->s_i(ii));
		} // next i
	if (f_v) {
		printf("aijk[%ld,+%ld] = ", first, len);
		for (i = 0; i < len; i++) {
			ii = first + i;
			aijk->s_i(ii)->print();
			if (i < len)
				printf(", ");
			}
		printf("\n");
		}
	return OK;
}

#if TEXDOCU
INT plesken_decompose_with_inverse(MATRIX_OP Ainf_inv, 
	VECTOR_OP aiaj, VECTOR_OP aijk, INT first, INT len, INT d0, INT f_v)
#else
aijk is already allocated.
#endif
{
	INT i, ii, j, nb_dc;
	SYM_OP a, b;
	SYM_OB c, tmp, tmp1;

	nb_dc = Ainf_inv->s_hi() + d0;
	for (i = 0; i < len; i++) {
		ii = first + i;
		c.freeself();
		for (j = ii; j < nb_dc; j++) {
			a = Ainf_inv->s_ij(ii - d0, j - d0);
			b = aiaj->s_i(j);
			if (j == ii) {
				a->mult(b, &c);
				}
			else {
				tmp.freeself();
				tmp1.freeself();
				a->mult(b, &tmp);
				c.add(&tmp, &tmp1);
				tmp1.swap(&c);
				}
			} // next j
		c.swap(aijk->s_i(ii));
		} // next i
	if (f_v) {
		printf("aijk[%ld,+%ld] = ", first, len);
		for (i = 0; i < len; i++) {
			ii = first + i;
			aijk->s_i(ii)->print();
			if (i < len)
				printf(", ");
			}
		printf("\n");
		}
	return OK;
}

#if TEXDOCU
INT plesken_decomposition_collect_layer_wise(VECTOR_OP orbit_length, 
	VECTOR_OP aijk, VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, VECTOR_OP type, INT f_v, INT f_vv)
#endif
{
	INT i, layer, f, l;
	SYM_OP pa;
	
	type->m_il_n(num_layers);
	for (i = 0; i < num_layers; i++) {
		layer = highest_layer - i;
		f = K_first->s_ii(layer);
		l = K_len->s_ii(layer);
		pa = type->s_i(i);
		plesken_decomposition_collect_layer(orbit_length, 
			aijk, f, l, pa, f_vv);
		}
	if (f_v) {
		printf("plesken_decomposition_collect_layer_wise() "
			"highest_layer=%ld num_layers=%ld, type=", highest_layer, num_layers);
		type->println();
		}
	return OK;
}

#if TEXDOCU
INT plesken_decomposition_collect_layer(VECTOR_OP orbit_length, 
	VECTOR_OP aijk, INT first, INT len, SYM_OP sum, INT f_v)
#else
Collects the intersections according to the layers of the lattice, 
e.g. according to the values in K.
#endif
{
	INT i, ii;
	SYM_OP pa, pb;
	SYM_OB tmp, tmp1;
	
	sum->m_i_i(0);
	for (i = 0; i < len; i++) {
		ii = first + i;
		pa = aijk->s_i(ii);
		if (pa->nullp())
			continue;
		pb = orbit_length->s_i(ii);
		tmp.freeself();
		pa->mult(pb, &tmp);
		sum->add(&tmp, &tmp1);
		tmp1.swap(sum);
		tmp1.freeself();
		}
	if (f_v) {
		printf("plesken_decomposition_collect_layer() [%ld,%ld] = ", first, len);
		sum->println();
		}
	return OK;
}



#endif /* LADDER_TRUE */



