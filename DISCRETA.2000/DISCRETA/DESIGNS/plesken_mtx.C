/* plesken_mtx.C */

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


#undef TEST_AINF_INV

#if TEXDOCU
INT dc_plesken_layer_idx(VECTOR_OP K_first, VECTOR_OP K_len, INT i)
#endif
{
	INT i1, len, f, l;
	
	len = K_len->s_li();
	for (i1 = 0; i1 < len; i1++) {
		f = K_first->s_ii(i1);
		l = K_len->s_ii(i1);
		if (i >= f && i < f + l) {
			break;
			}
		}
	if (i1 >= len) {
		return error("dc_plesken_layer_idx() i out of range");
		}
	return i1;
}

#if TEXDOCU
INT dc_plesken_prepare(INT k_min, INT k_max, VECTOR_OP MM, 
	MATRIX_OP Ainf_block_wise, MATRIX_OP coeff, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT f_v, INT f_vv)
#endif
{
	if (f_v) {
		printf("dc_plesken_prepare (k_min = %ld k_max = %ld)\n", k_min, k_max);
		fflush(stdout);
		}
	dc_calc_K_first_len_from_MM(MM, k_max, K_first, K_len);
	if (f_v) {
		printf("K_first: ");
		K_first->println();
		printf("K_len: ");
		K_len->println();
		}
	
	dc_Ainf_inv_coeff_ma(k_max, coeff, f_vv);
	printf("coeff=\n");
	coeff->Print();
	fflush(stdout);
	dc_Ainf_block_wise_via_MM(k_min, k_max, Ainf_block_wise, MM, f_vv);
	printf("Ainf_block_wise computed\n");
	fflush(stdout);
	return OK;
}

#if TEXDOCU
INT dc_plesken_matrix_prepare(INT k_min, INT k_max, VECTOR_OP MM, 
	MATRIX_OP Ainf_t, MATRIX_OP coeff, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT f_v, INT f_vv)
#endif
{
	MATRIX_OB Ainf_block_wise;
	MATRIX_OP Mij;
	SYM_OP c1;
	INT i, j, fi, li, fj, lj, nb_d, d0;
	INT ii, jj;
	
	if (f_v) {
		printf("preparing plesken matrix (k_max = %ld) \n", k_max);
		fflush(stdout);
		}
	dc_calc_K_first_len_from_MM(MM, k_max, K_first, K_len);
	d0 = K_first->s_ii(k_min);
	nb_d = K_first->s_ii(k_max) + K_len->s_ii(k_max)  - d0;
	if (f_v) {
		printf("number of orbits from k_min = %ld to k_max = %ld: %ld\n", 
			k_min, k_max, d0, nb_d);
		printf("K_first: ");
		K_first->println();
		printf("K_len: ");
		K_len->println();
		}
	
	dc_Ainf_inv_coeff_ma(k_max, coeff, f_vv);
	printf("coeff=\n");
	coeff->Print();
	fflush(stdout);
	dc_Ainf_block_wise_via_MM(k_min, k_max, &Ainf_block_wise, MM, f_vv);
	printf("Ainf_block_wise computed\n");
	fflush(stdout);
	
	Ainf_t->m_ilih_n(nb_d, nb_d);
	for (i = 0; i < nb_d; i++) {
		Ainf_t->m_iji(i, i, 1);
		}
	for (i = k_min; i <= k_max; i++) {
		fi = K_first->s_ii(i);
		li = K_len->s_ii(i);
		for (j = i + 1; j <= k_max; j++) {
			fj = K_first->s_ii(j);
			lj = K_len->s_ii(j);
			printf("i=%ld j=%ld\n", i, j);
			fflush(stdout);
			Mij = (MATRIX_OP) Ainf_block_wise.s_ij(i - k_min, j - k_min);
			if (Mij->s_hi() != li)
				return error("dc_plesken_matrix_prepare() Mij->s_hi() != li");
			if (Mij->s_li() != lj)
				return error("dc_plesken_matrix_prepare() Mij->s_li() != lj");
			if (f_vv) {
				printf("i=%ld,j=%ld\n", i, j);
				fflush(stdout);
				}
			for (ii = 0; ii < li; ii++) {
				for (jj = 0; jj < lj; jj++) {
					c1 = Ainf_t->s_ij(fi - d0 + ii, fj - d0 + jj);
					Mij->s_ij(ii, jj)->copy(c1);
					}
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT dc_plesken_matrices_prepare(INT k_min, INT k_max, VECTOR_OP MM, 
	MATRIX_OP Ainf_block_wise, MATRIX_OP Ainf_t, MATRIX_OP Ainf_t_inv, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT f_v, INT f_vv)
#else
Computes Ainf\_block\_wise, Ainf\_t, Ainf\_t\_inv, K\_first, K\_len.
The input for this data is MM, the vector of neighbouring 
Kramer Mesner matrices $M_{i,i+1}$ for $0 \le i < k$ ($k = k\_max$).
Ainf\_block\_wise has $M_{ij}$ in its $(i,j)$-th position.
$M_{ij}$ is the Kramer Mesner matrix of all $i$-orbits and $j$-orbits.
This matrix contains the data used for building up 
$(A^\wedge)^t$ and $((A^\wedge)^t)^{-1}$.
K\_first and K\_len contain the index of the first orbit of G on $i$-sets 
and the number of $i$-orbits respectively.
#endif
{
	MATRIX_OB coeff;
	MATRIX_OP Mij;
	SYM_OP c1, c2;
	SYM_OB b;
	INT i, j, fi, li, fj, lj, nb_d, d0;
	INT ii, jj;
	
	if (f_v) {
		printf("preparing plesken matrices (k_max = %ld) \n", k_max);
		fflush(stdout);
		}
	dc_calc_K_first_len_from_MM(MM, k_max, K_first, K_len);
	d0 = K_first->s_ii(k_min);
	nb_d = K_first->s_ii(k_max) + K_len->s_ii(k_max)  - d0;
	if (f_v) {
		printf("number of orbits from k_min = %ld to k_max = %ld: %ld\n", 
			k_min, k_max, d0, nb_d);
		printf("K_first: ");
		K_first->println();
		printf("K_len: ");
		K_len->println();
		}
	
	dc_Ainf_inv_coeff_ma(k_max, &coeff, f_vv);
	printf("coeff=\n");
	coeff.Print();
	fflush(stdout);
	dc_Ainf_block_wise_via_MM(k_min, k_max, Ainf_block_wise, MM, f_vv);
	printf("Ainf_block_wise\n");
	fflush(stdout);
	
	Ainf_t->m_ilih_n(nb_d, nb_d);
	Ainf_t_inv->m_ilih_n(nb_d, nb_d);
	for (i = 0; i < nb_d; i++) {
		Ainf_t->m_iji(i, i, 1);
		Ainf_t_inv->m_iji(i, i, 1);
		}
	for (i = k_min; i <= k_max; i++) {
		fi = K_first->s_ii(i);
		li = K_len->s_ii(i);
		for (j = i + 1; j <= k_max; j++) {
			fj = K_first->s_ii(j);
			lj = K_len->s_ii(j);
			printf("i=%ld j=%ld\n", i, j);
			fflush(stdout);
			Mij = (MATRIX_OP) Ainf_block_wise->s_ij(i - k_min, j - k_min);
			if (Mij->s_hi() != li)
				return error("dc_plesken_matrices_prepare() Mij->s_hi() != li");
			if (Mij->s_li() != lj)
				return error("dc_plesken_matrices_prepare() Mij->s_li() != lj");
			b.m_i_i(coeff.s_iji(i, j));
			if (f_vv) {
				printf("i=%ld,j=%ld, coeff= %ld\n", i, j, b.s_i_i());
				fflush(stdout);
				}
			for (ii = 0; ii < li; ii++) {
				for (jj = 0; jj < lj; jj++) {
					c1 = Ainf_t->s_ij(fi - d0 + ii, fj - d0 + jj);
					c2 = Ainf_t_inv->s_ij(fi - d0 + ii, fj - d0 + jj);
					Mij->s_ij(ii, jj)->copy(c1);
					Mij->s_ij(ii, jj)->mult(&b, c2);
					}
				}
			}
		}
#ifdef TEST_AINF_INV
	{
	MATRIX_OB I;
	
	Ainf_t->mult(Ainf_t_inv, &I);
	if (!I.einsp())
		return error("dc_plesken_matrices_prepare() Ainf_t has wrong inverse!");
	printf("dc_plesken_matrices_prepare() Ainf_t inverse checked and OK !\n");
	fflush(stdout);
	}
#endif
	return OK;
}

#if TEXDOCU
INT dc_Ainf_inv_coeff_ma(INT k_max, MATRIX_OP coeff, INT f_v)
#else
Computes a utility matrix which is used for determination of $(A^\wedge)^{-1}$.
In fact, this matrix is an upper triangular matrix with 
$(-1)^{i+j}$ in its entries (upper triangle).
#endif
{
	INT i, k, t, r, bb;
	SYM_OP n;
	SYM_OB c, b, tmp, tmp1;
	
	coeff->m_ilih_n(k_max + 1, k_max + 1);
	for (i = 0; i <= k_max; i++) {
		coeff->m_iji(i, i, 1);
		}
	for (k = k_max; k >= 1; k--) {
		for (t = k - 1; t >= 0; t--) {
			c.m_i_i(-1);
			if (f_v) {
				printf("t = %ld k = %ld\n", t, k);
				}
			for (r = k - 1; r > t; r--) {
				// Binomial(k - t, k - r, &b);
				bb = n_choose_k(k - t, k - r);
				b.m_i_i(bb);
				n = coeff->s_ij(r, k);
				b.mult(n, &tmp);
				tmp.addinvers_apply();
				c.add(&tmp, &tmp1);
				if (f_v) {
					printf("r = %ld\n", r);
					printf("n:\n");
					n->println();
					printf("k-t atop k-r = %ld\n", bb);
					printf("c:\n");
					c.println();
					}
				tmp1.swap(&c);
				}
			c.swap(coeff->s_ij(t, k));
			}
		}
	if (f_v) {
		printf("dc_Ainf_inv_coeff_ma() k_max = %ld\n", k_max);
		coeff->fprint_raw(stdout);
		fflush(stdout);
		}
	return OK;
}

#if TEXDOCU
INT dc_Ainf_block_wise_via_MM(INT k_min, INT k_max, MATRIX_OP Ainf_block_wise, VECTOR_OP MM, INT f_v)
#else
Fast computation of Ainf\_block\_wise out of MM. 
Used by dc\_plesken\_matrices\_prepare.
#endif
{
	INT l, d, d_max, i, i0, t, k, r, m, n;
	MATRIX_OP pMtk, pMtr, pMrk;
	
	l = k_max - k_min + 1;
	Ainf_block_wise->m_ilih(l, l);
	for (i = k_min; i < k_max; i++) {
		MM->s_i(i)->copy(Ainf_block_wise->s_ij(i - k_min, i - k_min + 1));
		}
	d_max = k_max;
	for (d = 2; d <= d_max; d++) {
		printf("d=%ld\n", d);
		fflush(stdout);
		for (i = k_min; i <= k_max - d; i++) {
			i0 = i - k_min;
			t = i0;
			k = i0 + d;
			r = i0 + d - 1;
			pMtk = (MATRIX_OP) Ainf_block_wise->s_ij(t, k);
			pMtr = (MATRIX_OP) Ainf_block_wise->s_ij(t, r);
			pMrk = (MATRIX_OP) Ainf_block_wise->s_ij(r, k);
			m = pMtr->s_hi();
			n = pMrk->s_li();
			printf(" M_{%ld,%ld} via M_{%ld,%ld} and M_{%ld,%ld} "
				"of size %ld x %ld\n", t, k, t, r, r, k, m, n);
			fflush(stdout);
			dc_Mtk_via_Mtr_Mrk(t, r, k, pMtr, pMrk, pMtk, f_v);
			}
		}
	printf("finished!\n");
	return OK;
}

#if TEXDOCU
INT dc_Mtk_via_Mtr_Mrk(INT t, INT r, INT k, MATRIX_OP Mtr, MATRIX_OP Mrk, MATRIX_OP Mtk, INT f_v)
#else
Computes $M_{tk}$ via a recursion formula:
$M_{tk} = {{k - t} \choose {k - r}} \cdot M_{t,r} \cdot M_{r,k}$.
#endif
{
	MATRIX_OB M;
	SYM_OB s, h1, h2, h3;
	SYM_OP h;
	INT i, j;
	
	if (f_v) {
		printf("Mtk_via_Mtr_Mrk(): t = %ld, r = %ld, k = %ld\n", t, r, k);
		fflush(stdout);
		}
	if (Mtr->s_obj_k() != MATRIX)
		return error("dc_Mtk_via_Mtr_Mrk() Mtr->s_obj_k() != MATRIX");
	if (Mrk->s_obj_k() != MATRIX)
		return error("dc_Mtk_via_Mtr_Mrk() Mrk->s_obj_k() != MATRIX");
	Mtr->mult(Mrk, &M);
		/* M := (k - t) atop (k - r) * M_t,k */
	Binomial(k - t, k - r, &s);
	for (i = 0; i < M.s_hi(); i++) {
		for (j = 0; j < M.s_li(); j++) {
			h = M.s_ij(i, j);
			h->ganzdiv(&s, &h1);
				
			h1.mult(&s, &h2);
			h2.addinvers_apply();
			h2.add(h, &h3);
			if (!h3.nullp()) {
			// if (h1 * s != h) {
				printf("t = %ld\n", t);
				printf("r = %ld\n", r);
				printf("k = %ld\n", k);
				printf("h = ");
					h->println();
				printf("s = ");
					s.println();
				printf("h1 = ");
					h1.println();
				printf("i = %ld\n", i);
				printf("j = %ld\n", j);
				printf("M:\n");
				M.Print();
				return error("dc_Mtk_via_Mtr_Mrk()|h1 * s != h");
				}
			h1.swap(M.s_ij(i, j));
			}
		}
	M.swap(Mtk);
	return OK;
}



#if TEXDOCU
INT dc_calc_K_first_len(VECTOR_OP K, INT k_max, VECTOR_OP K_first, VECTOR_OP K_len)
#endif
{
	INT i, f, l;
	
	K_first->m_il(k_max + 1);
	K_len->m_il(k_max + 1);
	for (i = 0; i <= k_max; i++) {
		dc_find_first_len(K, i, &f, &l);
		K_first->m_ii(i, f);
		K_len->m_ii(i, l);
		}
	return OK;
}

#if TEXDOCU
INT dc_calc_K_first_len_from_MM(VECTOR_OP MM, INT k_max, VECTOR_OP K_first, VECTOR_OP K_len)
#endif
{
	MATRIX_OP KM;
	INT i, f = 0, l;
	
	K_first->m_il(k_max + 1);
	K_len->m_il(k_max + 1);
	for (i = 0; i < k_max; i++) {
		KM = (MATRIX_OP) MM->s_i(i);
		l = KM->s_hi();
		K_first->m_ii(i, f);
		K_len->m_ii(i, l);
		f += l;
		if (i == k_max - 1) {
			l = KM->s_li();
			K_first->m_ii(i + 1, f);
			K_len->m_ii(i + 1, l);
			}
		}
	return OK;
}

#if TEXDOCU
INT dc_find_first_len(VECTOR_OP K, INT t, INT *first, INT *len)
#endif
{
	INT i, j, nb_dc;
	
	nb_dc = K->s_li();
	for (i = 0; i < nb_dc; i++) {
		if (K->s_ii(i) >= t) {
			*first = i;
			for (j = i + 1; j < nb_dc; j++)
				if (K->s_ii(j) > t)
					break;
			*len = j - i;
			// printf("%ld sets: first = %ld len = %ld\n", t, *first, *len);
			return OK;
			}
		}
	return error("find_first_len(): not found!");
}


#endif /* LADDER_TRUE */

