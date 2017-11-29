/* intersection.C */

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

static INT compare_vector_equality(MATRIX_OP A, VECTOR_OP B);
static INT unipoly_koefficient(UNIPOLY_OP p, INT i);

#if TEXDOCU
INT manage_block_alpha_i_data(VECTOR_OP block_alpha_i, 
	INTEGER_OP bai_len, 
	VECTOR_OP bai_ref, VECTOR_OP bai_refv, 
	VECTOR_OP bai_data, INT *ref, INT f_v)
#else
bai\_ref is a permutation, bai\_refv its inverse !
#endif
{
	INT i, a, l, idx, f_found, ref1;
	INTEGER_OB int_ob;
	
	l = bai_len->s_i();
	if (l == 0) {
		bai_ref->m_il(VECTOR_OVERSIZE);
		bai_refv->m_il(VECTOR_OVERSIZE);
		bai_data->m_il(VECTOR_OVERSIZE);
		}
	bai_data->search(l, TRUE, block_alpha_i, &idx, &f_found);
	if (f_found) {
		idx--;
		ref1 = bai_refv->s_ii(idx);
		if (f_v)
			printf("block_alpha_i ref %ld found !\n", ref1);
		*ref = ref1;
		}
	else {
		bai_data->insert_at(bai_len, idx, block_alpha_i);
		bai_len->dec();
		int_ob.m_i(l);
		bai_refv->insert_at(bai_len, idx, &int_ob);
		bai_len->dec();
		int_ob.m_i(idx);
		bai_ref->insert_at(bai_len, l, &int_ob);
		for (i = 0; i < l; i++) {
			a = bai_ref->s_ii(i);
			if (a >= idx) {
				a++;
				bai_ref->m_ii(i, a);
				}
			}
		if (f_v)
			printf("block_alpha_i ref %ld newly defined !\n", l);
		*ref = l;
		}
	return OK;
}

#if TEXDOCU
INT total_intersection_3_tuples(MATRIX_OP X, INT k,
	SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	MATRIX_OP Ainf, 
	VECTOR_OP intersection_3_tuples)
#endif
{
	INT a, i, j, kk, first, len, n, deg;
	SYM_OP pa, pb;
	SYM_OB tmp;
	UNIPOLY_OB ii;
	
	dc_find_first_len(K, k, &first, &len);
	intersection_3_tuples->m_il_n(k + 1);
	n = X->s_hi();
	for (i = 0; i < n; i++) {
		if (X->s_iji(i, 0) == 0)
			continue;
		for (j = 0; j < n; j++) {
			if (X->s_iji(j, 0) == 0)
				continue;
			for (kk = 0; kk < n; kk++) {
				if (X->s_iji(kk, 0) == 0)
					continue;
				Intersection_of3_via_Ainf(go, stab_go, K, Ainf, 
					first + i, first + j, first + kk, 
					&ii, FALSE /* f_v */, FALSE /* f_vv */);
				deg = ii.degree();
				for (a = 0; a <= deg; a++) {
					pa = ii.s_i(a);
					pb = intersection_3_tuples->s_i(a);
					pa->add(pb, &tmp);
					tmp.swap(pb);
					}
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT total_intersection_3_sets(MATRIX_OP X, INT k,
	SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	MATRIX_OP Ainf, 
	VECTOR_OP intersection_3_sets)
#endif
{
	INT i, ii, j, kk, first, len, n;
	// SYM_OP pa, pb;
	// SYM_OB tmp;
	VECTOR_OB inter_tmp, inter_tmp2; 
	SYM_OB two_ob, six_ob, minus_three_ob;
	
	dc_find_first_len(K, k, &first, &len);
	intersection_3_sets->m_il_n(k + 1);
	n = X->s_hi();
	two_ob.m_i_i(2);
	six_ob.m_i_i(6);
	minus_three_ob.m_i_i(-3);

	// first case: 
	// all blocks belong to {\em different} block orbits:
	printf("%% different block orbits: ");
	ii = 0;
	for (i = 0; i < n; i++) {
		if (X->s_iji(i, 0) == 0)
			continue;
		ii++;
		}
	printf(" (%ld orbits) \n", ii);
	printf("%% ");
	ii = 0;
	for (i = 0; i < n; i++) {
		if (X->s_iji(i, 0) == 0)
			continue;
		if ((ii + 1) % 10 == 0) {
			if ((ii + 1) % 50 == 0)
				printf(", %ld\n%%", ii);
			else
				printf(",");
			}
		else
			printf(".");
		ii++;
		fflush(stdout);
		for (j = i + 1; j < n; j++) {
			if (X->s_iji(j, 0) == 0)
				continue;
			for (kk =  j + 1; kk < n; kk++) {
				if (X->s_iji(kk, 0) == 0)
					continue;
				intersection_of3_and_add(
					first + i, first + j, first + kk, 
					go, stab_go, K, Ainf, 
					intersection_3_sets, 1 /* sign */);
#if 0
				Intersection_of3_via_Ainf(
					go, stab_go, K, Ainf, 
					first + i, first + j, first + kk, 
					&ii, TRUE /* f_v */, FALSE /* f_vv */);
				deg = ii.degree();
				for (a = 0; a <= deg; a++) {
					pa = ii.s_i(a);
					pb = intersection_3_sets->s_i(a);
					pa->add(pb, &tmp);
					tmp.swap(pb);
					}
#endif
				}
			}
		}
	printf("\n");

	// second case:
	// the case when 2 (of the 3) blocks belong to the same block-orbit (i):
	// let x and y be the elements of block orbit i, 
	// z the element of orbit j.
	// the triples $(x,y,z)$ with $x \neq y$ are counted twice
	// the triples $(x,x,z)$ are counted once.
	printf("%% two intersections (different orbits): ");
	fflush(stdout);
	inter_tmp.m_il_n(k + 1);
	for (i = 0; i < n; i++) {
		if (X->s_iji(i, 0) == 0)
			continue;
		for (j = 0; j < n; j++) {
			if (j == i)
				continue;
			if (X->s_iji(j, 0) == 0)
				continue;
			intersection_of3_and_add(
				first + i, first + i, first + j, 
				go, stab_go, K, Ainf, 
				&inter_tmp, 1);
			}
		}
	// subtract the intersections for those cases when 
	// the first two elements are the same element in the i-th block orbit:
	// the case $(x,x,z)$.
	for (i = 0; i < n; i++) {
		if (X->s_iji(i, 0) == 0)
			continue;
		for (j = 0; j < n; j++) {
			if (j == i)
				continue;
			if (X->s_iji(j, 0) == 0)
				continue;
			intersection_of3_and_add(
				first + i, first + j, -1, 
				go, stab_go, K, Ainf, 
				&inter_tmp, -1);
			}
		}
	// and divide by two because the elements in the i-th block orbit 
	// are counted twice (reordering !)
	inter_tmp.divide_out(&two_ob);

	// add inter_tmp to intersection_3_sets:
	intersection_3_sets->add_apply_elementwise(&inter_tmp, 1);
	printf("\n");


	// third case:
	// now the third case, namely when all elements belong to the same block orbit.
	printf("%% three elements of the same orbit: ");
	fflush(stdout);
	inter_tmp.m_il_n(k + 1);
	for (i = 0; i < n; i++) {
		if (X->s_iji(i, 0) == 0)
			continue;
		intersection_of3_and_add(
			first + i, first + i, first + i, 
			go, stab_go, K, Ainf, 
			&inter_tmp, 1);
		// counts all cases of elements $x,y,z$: 
		// all three elements distinct (6 times because of the ordering), 
		// two elements coincide ($x=y$)  (3 times for each fixed z) and
		// all three elements coincide (once).

		// we have to subtract the cases with coincidences and 
		// divide out the factor 6.

		// we dont want the $x=y=z$ cases:
		inter_tmp.s_i(k)->m_i_i(0);
		
		inter_tmp2.m_il_n(k + 1);
		// all cases of intersections of 2 elements:
		// two different elements (2 times) and
		// two equal elements (once).
		intersection_of3_and_add(
			first + i, first + i, -1, 
			go, stab_go, K, Ainf, 
			&inter_tmp2, 1);
		
		// the highest coefficient gives the orbit length !
		
		// we do not want the $x=y=z$ cases:
		inter_tmp2.s_i(k)->m_i_i(0);
		
		// // divide out the factor 2:
		// divide_out(&inter_tmp2, &two_ob);
		
		inter_tmp2.multiply_elementwise(&minus_three_ob);
		
		inter_tmp.add_apply_elementwise(&inter_tmp2, 1);
		
		}
	printf("dividing out 6");
	fflush(stdout);
	// divide out the factor 6:
	inter_tmp.divide_out(&six_ob);
	// add inter_tmp to intersection_3_sets:
	intersection_3_sets->add_apply_elementwise(&inter_tmp, 1);
	printf("\n");
	return OK;
}

#if TEXDOCU
INT intersection_of3_and_add(INT i, INT j, INT k, 
	SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	MATRIX_OP Ainf, 
	VECTOR_OP intersection_3_sets, INT sign)
#else
#endif
{
	SYM_OP pa, pb;
	SYM_OB tmp;
	UNIPOLY_OB ii;
	INT deg, a;
	
	Intersection_of3_via_Ainf(go, stab_go, K, Ainf, 
		i, j, k, 
		&ii, FALSE /* f_v */, FALSE /* f_vv */);
	deg = ii.degree();
	for (a = 0; a <= deg; a++) {
		pa = ii.s_i(a);
		if (sign < 0)
			pa->addinvers_apply();
		pb = intersection_3_sets->s_i(a);
		pa->add(pb, &tmp);
		tmp.swap(pb);
		tmp.freeself();
		}
	return OK;
}

#if TEXDOCU
INT total_intersection_2_tuples_to_2_sets(VECTOR_OP intersection_2_tuples, 
	VECTOR_OP intersection_2_sets)
#endif
{
	SYM_OB two_ob;
	INT i, l;
	
	two_ob.m_i_i(2);
	intersection_2_tuples->copy(intersection_2_sets);
	l = intersection_2_sets->s_li();
	intersection_2_sets->s_i(l - 1)->zero();
		// no pairs of the form $(x,x)$
	for (i = 0; i < l - 1; i++) {
		(*intersection_2_sets)[i] /= two_ob;
		}
	return OK;
}

#if TEXDOCU
INT total_intersections_2_tuples(MATRIX_OP X, INT k, 
	SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	VECTOR_OP block_types, VECTOR_OP intersection_2_tuples, 
	INTEGER_OP bai_len, VECTOR_OP bai_ref, VECTOR_OP bai_refv,
	VECTOR_OP bai_data, INT f_v)
#endif
{
	INT i, j, idx, n, n1, block_type, first, len;
	VECTOR_OP T;
	SYM_OP ago, pa, pb;
	SYM_OB ol, tmp1, tmp2;
	
	dc_find_first_len(K, k, &first, &len);
	intersection_2_tuples->m_il_n(k + 1);
	n = X->s_hi();
	n1 = 0;
	for (i = 0; i < n; i++) {
		if (X->s_iji(i, 0) == 0)
			continue;
		block_type = block_types->s_ii(n1);
		if (block_type > bai_len->s_i())
			return error("total_intersections_2_tuples() block_type > bai_len->s_i()");
		// printf("block-type = %ld\n", block_type);
		// fflush(stdout);
		idx = bai_ref->s_ii(block_type);
		if (idx > bai_len->s_i())
			return error("total_intersections_2_tuples() idx > bai_len->s_i()");
		// printf("idx = %ld\n", idx);
		// fflush(stdout);

		ago = stab_go->s_i(first + i);
		go->ganzdiv(ago, &ol);
		// printf("ol = ");
		// ol.print();
		// printf("\\\\\n");
		
		T = (VECTOR_OP) bai_data->s_i(idx);
		if (T->s_li() != k + 1) 
			return error("total_intersections_2_tuples() T->s_li() != k + 1");
		// T->println();
		// printf("\\\\\n");
		// fflush(stdout);
		for (j = 0; j <= k; j++) {
			pa = intersection_2_tuples->s_i(j);
			pb = T->s_i(j);
			pb->mult(&ol, &tmp1);
			tmp1.add(pa, &tmp2);
			tmp2.swap(pa);
			}

		n1++;
		}
	if (n1 != block_types->s_li())
		return error("total_intersections_2_tuples() n1 != block_types->s_li()");
	return OK;
}

#if TEXDOCU
INT vec_sum_of_all_elements(VECTOR_OP V, SYM_OP S)
#endif
{
	INT i, l;
	SYM_OP pa;
	SYM_OB tmp;

	l = V->s_li();
	S->m_i_i(0);
	for (i = 0; i < l; i++) {
		pa = V->s_i(i);
		S->add(pa, &tmp);
		tmp.swap(S);
		}
	return OK;
}

#if TEXDOCU
INT check_intersection_eqns(MATRIX_OP X, MATRIX_OP Mendelsohn, 
	VECTOR_OP Eqns, VECTOR_OP RHS, VECTOR_OP block_types, 
	INTEGER_OP bai_len, VECTOR_OP bai_ref, VECTOR_OP bai_refv, 
	VECTOR_OP bai_data, INT f_v)
#else
This function checks the intersection equations for all orbits 
of the design (this means it checks them for all entries '1' in 
the solution vector X). 
Eqns is the vector of intersection equations, that is a vector of matrices 
of dimension $(k+1) \times l$ where $l$ is the number of $k$-orbits.
RHS is the right hand side of the intersection equations, e.g. the vector 
with ${k \choose i} \cdot \lambda_i$ in its $i$-th entry. 
For each block orbit, the function check\_intersection\_eqn is called. 

If bai\_len is given (the pointer is not NIL), a statistic of intersection 
types is computed. bai\_data is a sorted list of all intersection types 
and bai\_ref and bai\_refv are references to it (bai\_ref and bai\_refv 
are a pair of inverse permutations). 

The block\_intersection types are managed 
via the function manage\_block\_alpha\_i\_data.
Note that bai\_data need not be empty: it may contain block intersection vectors 
from previous solution vectors of the Kramer-Mesner system. 
BUT: bai\_data needs to be cleared if the value $\lambda$ of the design changes.
This is because the intersection equations (RHS) depend on $\lambda$ !

#endif
{
	INT a, i, n, n1, ref;
	MATRIX_OP eqn;
	VECTOR_OB block_alpha_i;

	n = X->s_hi();
	if (f_v) {
		printf("check_intersection_eqns()\n");
		printf("X=");
		for (i = 0; i < n; i++)  {
			X->s_ij(i, 0)->print();
			printf(" ");
			}
		printf("\n");
		}
	
	if (bai_len) {
		block_types->m_il(n);
		n1 = 0;
		}
	for (a = 0; a < n; a++) {
		if (X->s_iji(a, 0) == 0)
			continue;
		eqn = (MATRIX_OP) Eqns->s_i(a);
		if (f_v) {
			printf("checking for orbit %ld\n", a);
			printf("eqn is a %ld x %ld matrix\n", eqn->s_hi(), eqn->s_li());
			// eqn->Print();
			}
		check_intersection_eqn(a /* block_orbit */, X, Mendelsohn, 
			eqn, RHS, &block_alpha_i, f_v, FALSE /* f_vv */);
		if (bai_len) {
			manage_block_alpha_i_data(&block_alpha_i, bai_len, 
				bai_ref, bai_refv, bai_data, &ref, FALSE /* f_v */);
			block_types->m_ii(n1, ref);
			n1++;
			}
		}
	if (bai_len) {
		block_types->realloc_z(n1);
		}
	return OK;
}

#if TEXDOCU
INT check_intersection_eqn(INT block_orbit, MATRIX_OP X, MATRIX_OP Mendelsohn, 
	MATRIX_OP eqn, VECTOR_OP RHS, VECTOR_OP block_alpha_i, INT f_v, INT f_vv)
#else
eqn is the part of the Intersection system which describes 
the Intersection between a representative of the orbit block\_orbit 
and the set of all k-orbits (all k-orbits). 
eqn is a $(k + 1) \times l$ matrix  where $l$ is the numer of $k$-orbits of the group.
X is the solution vector (as a matrix of a single column)


This function has two modi whether or not 
Mendelsohn is given (the pointer is not NIL).
Assume first that Mendelsohn is given.
Then, $eqn \cdot X =$ block\_alpha\_i $= (\alpha_0, \ldots , \alpha_k)^t$. 
and block\_alpha\_i is returned. 
This function checks if $Mendelsohn \cdot$ block\_alpha\_i $=$ RHS.

In the case that Mendelsohn is NIL one assumes that mendelsohn 
already has been applied to eqn. In this case, 
the equality $eqn \cdot X = $RHS is checked. 
In this case, block\_alpha\_i is unused.
#endif
{
	INT i, l, n, ret;
	MATRIX_OB Y, Z;

	n = X->s_hi();
	if (f_v) {
		printf("check_intersection_eqn() block_orbit = %ld\n", block_orbit);
		printf("eqn is a %ld x %ld matrix\n", eqn->s_hi(), eqn->s_li());
		}
	
	eqn->mult(X, &Y);
	if (f_vv) {
		printf("Y = ");
		Y.Print();
		}
	if (Mendelsohn) {
		if (block_alpha_i) {
			l = Y.s_hi();
			block_alpha_i->m_il(l);
			for (i = 0; i < l; i++) 
				Y.s_ij(i, 0)->copy(block_alpha_i->s_i(i));
			// Y.copy(block_alpha_i);
			}
		if (f_v) {
			printf("applying Mendelsohn system\n");
			}
		Mendelsohn->mult(&Y, &Z);
		Z.swap(&Y);
		if (f_vv) {
			printf("Y = ");
			Y.Print();
			}
		}
	ret = compare_vector_equality(&Y, RHS);
	if (!ret) {
		printf("WARNING !: intersection eqn for block orbit %a not satisfied\n", block_orbit);
		Y.println();
		RHS->println();
		}
	return OK;
}

#if TEXDOCU
static INT compare_vector_equality(MATRIX_OP A, VECTOR_OP B)
#endif
{
	SYM_OP pa, pb;
	SYM_OB c, d;
	INT i, l;

	l = A->s_hi();
	if (B->s_li() != l)
		return error("compare_vector() different size");
	for (i = 0; i < l; i++) {
		pa = A->s_ij(i, 0);
		pb = B->s_i(i);
		pb->copy(&c);
		c.addinvers_apply();
		pa->add(&c, &d);
		if (!d.nullp())
			return FALSE;
		}
	return TRUE;
}

#if TEXDOCU
INT eqns_apply(VECTOR_OP Eqns, MATRIX_OP M)
#else
Eqns is a vector of matrices.
applies M from the left to all matrices in Eqns.
#endif
{
	MATRIX_OP peqn;
	MATRIX_OB eqn;
	INT i, l;

	l = Eqns->s_li();
	for (i = 0; i < l; i++) {
		peqn = (MATRIX_OP) Eqns->s_i(i);
		M->mult(peqn, &eqn);
		eqn.swap((MATRIX_OP) Eqns->s_i(i));
		}
	return OK;
}

#if TEXDOCU
INT intersection_eqns(VECTOR_OP Eqns, MATRIX_OP Ikk, INT k, INT f_v)
#else
Ikk is the Matrix over unipoly describing the intersections between 
all pairs of orbits of $k$-sets. (the matrix is symmetric up to scalar factors). 
The unipolys are used just as carriers for the intersection coefficients, 
so the unipoly $1X^k + \ldots + 2X^1 + 3X^0$ 
stands for $1 \alpha_k + \ldots + 2 \alpha_2 + 3 \alpha_0$, e.g. 
for 1 intersection of cardinality $k$, ... 
2 intersections of size 1 and 3 disjoint blocks. Note that 
$\alpha_k \neq 0$ is only possible on the diagonal of Ikk, i.e. 
for intersections of a block orbits with itself !
(we are dealing with simple designs only).

This routine transforms the Ikk data into Eqns, a vector of length $l$ 
of $(k+1) \times l$ -matrices ($l$ the number of $k$-orbits under $G$).
#endif
{
	MATRIX_OB eqn;
	INT i, j, a, l, b;
	UNIPOLY_OP ii;
	
	l = Ikk->s_hi();
	if (l != Ikk->s_li())
		return error("intersection_eqns() Ikk not quadratic !");
	Eqns->m_il(l);
	for (i = 0; i < l; i++) {
		eqn.m_ilih_n(l, k + 1);
		for (j = 0; j < l; j++) {
			ii = (UNIPOLY_OP) Ikk->s_ij(i, j);
			for (a = 0; a <= k; a++) {
				b = unipoly_koefficient(ii, a);
				eqn.m_iji(a, j, b);
				}
			}
		if (f_v) {
			printf("eqn %ld (of %ld):\n", i + 1, l);
			eqn.Print();
			}
		eqn.swap((MATRIX_OP) Eqns->s_i(i));
		}
	return OK;
}

#if TEXDOCU
static INT unipoly_koefficient(UNIPOLY_OP p, INT i)
#else
returns the coefficient of $x^i$ in the unipoly. 
#endif
{
	if (i >= p->s_li())
		return 0;
	return p->s_ii(i);
}

#if TEXDOCU
INT Intersection_Mtk_via_Ainf(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	INT t, INT k, MATRIX_OP Mtk, MATRIX_OP Ainf, INT f_v, INT f_vv)
#else
Computes $I_{t,k}$, e.g. the intersection matrix between $t$ and $k$-orbits. 
$t = k$ gives Ikk. 
The intersections are computed via Plesken theory, e.g. via Ainf.
#endif
{
	INT i, j, t_first, t_len, k_first, k_len, nb_dc;
	UNIPOLY_OB ii;

	printf("Intersection matrix %ld / %ld sets:\n", t, k); fflush(stdout);
	nb_dc = K->s_li();
	dc_find_first_len(K, t, &t_first, &t_len);
	dc_find_first_len(K, k, &k_first, &k_len);
	Mtk->m_ilih(k_len, t_len);
	for (i = 0; i < t_len; i++) {
		for (j = 0; j < k_len; j++) {
			Intersection_via_Ainf(go, stab_go, K, Ainf, 
				t_first + i, k_first + j, &ii, f_vv, f_vv);
			ii.swap(Mtk->s_ij(i, j));
			}
		if ((i + 1) % 10 == 0)
			printf(",");
		else
			printf(".");
		if ((i + 1) % 50 == 0)
			printf("%ld\n", i + 1);
		fflush(stdout);
		}
	return OK;
}

#if TEXDOCU
INT Intersection_M_via_Ainf(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	MATRIX_OP M, MATRIX_OP Ainf, INT f_v, INT f_vv)
#else
Computes intersection matrix labelled by all $G$-orbits 
(on sets of size $\le k$). 
#endif
{
	INT nb_d, i, j;
	UNIPOLY_OB ii;

	nb_d = Ainf->s_li();
	M->m_ilih_n(nb_d, nb_d);
	for (i = 0; i < nb_d; i++) {
		for (j = 0; j < nb_d; j++) {
			Intersection_via_Ainf(go, stab_go, K, Ainf, i, j, &ii, f_v, f_vv);
			ii.swap(M->s_ij(i, j));
			}
		}
	return OK;
}

#if TEXDOCU
INT Intersection_via_Ainf(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, MATRIX_OP Ainf, 
	INT dcno1, INT dcno2, UNIPOLY_OP ii, INT f_v, INT f_vv)
#else
Cmputes the intersections between 
a representative of a dcno1 orbit and (all blocks of) dcno2 orbit of $G$. 
Result is a unipoly ii. 
Decomposes the pointwise product of the columns dcno1 and dcno2 of Ainf 
as a sum of columns of Ainf.
Collects the intersections according to the layers of the lattice, 
e.g. according to the values in K.
#endif
{
	INT nb_dc, i;
	VECTOR_OB aiaj, aijk;
	SYM_OP pa, pb, pc;

	if (f_v) {
		printf("Intersection_via_Ainf() dcno1 = %ld dcno2 = %ld\n", dcno1, dcno2);
		}
	nb_dc = Ainf->s_li();
	aiaj.m_il(nb_dc);
	for (i = 0; i < nb_dc; i++) {
		pa = Ainf->s_ij(i, dcno1);
		pb = Ainf->s_ij(i, dcno2);
		pc = aiaj.s_i(i);
		pa->mult(pb, pc);
		}
	if (f_vv) {
		printf("aiaj = ");
		aiaj.println();
		}
	aijk.m_il_n(nb_dc);
	Plesken_decompose(Ainf, &aiaj, &aijk, f_vv);

	collect_intersections_by_layer(go, stab_go, K, 
		Ainf, &aijk, ii, f_vv);

	unipoly_divide_out_oli(ii, go, stab_go, dcno1);
	
	return OK;
}

#if TEXDOCU
INT Intersection_of3_via_Ainf(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, MATRIX_OP Ainf, 
	INT i1, INT i2, INT i3, UNIPOLY_OP ii, INT f_v, INT f_vv)
#else
Cmputes the intersections between triples of blocks from $G$-orbits 
i1, i2 and i3 respectively. 
Result is a unipoly ii. 
Decomposes the pointwise product of the columns (i1,i2,i3) of Ainf 
as a sum of columns of Ainf.
Collects the intersections according to the layers of the lattice, 
e.g. according to the values in K.
New: if $i3 < 0$ this routine decomposes only the product of columns i1 and i2 !
#endif
{
	INT nb_dc, i;
	VECTOR_OB aiaj, aijk;
	SYM_OP pa, pb, pc;
	SYM_OB tmp;

	if (f_v) {
		printf("Intersection_of3_via_Ainf() i1 = %ld i2 = %ld i3 = %ld\n", i1, i2, i3);
		}
	nb_dc = Ainf->s_li();
	aiaj.m_il(nb_dc);
	for (i = 0; i < nb_dc; i++) {
		pa = Ainf->s_ij(i, i1);
		pb = Ainf->s_ij(i, i2);
		pc = aiaj.s_i(i);
		pa->mult(pb, pc);
		if (i3 >= 0) {
			pc->mult(Ainf->s_ij(i, i3), &tmp);
			tmp.swap(pc);
			tmp.freeself();
			// aiaj[i] *= *(Ainf->s_ij(i, i3));
			}
		}
	if (f_vv) {
		printf("aiaj = ");
		aiaj.println();
		}
	aijk.m_il_n(nb_dc);
	Plesken_decompose(Ainf, &aiaj, &aijk, f_vv);
	if (f_v) {
		printf("aijk = ");
		aijk.println();
		}

	collect_intersections_by_layer(go, stab_go, K, 
		Ainf, &aijk, ii, f_vv);

	return OK;
}

#if TEXDOCU
INT Plesken_decompose(MATRIX_OP Ainf, VECTOR_OP aiaj, VECTOR_OP aijk, INT f_v)
#else
Decomposes the hadamard product of vectors $a_ia_j$ in aiaj 
into $\sum_k a_{ijk} a_k$. 
$a_k$ denotes the $k$-th column of Ainf. 
This routine knows (assumes) that Ainf is upper triangular.
#endif
{
	SYM_OP pa, pb;
	SYM_OB a0, tmp1, tmp2;
	INT nb_dc, i, j;
	
	nb_dc = Ainf->s_li();
	for (i = nb_dc - 1; i >= 0; i--) {
		aiaj->s_i(i)->copy(&a0);
		if (a0.nullp())
			continue;
		if (f_v) {
			printf("i=%ld a0 = ", i);
			a0.println();
			}
		a0.copy(aijk->s_i(i));
		a0.addinvers_apply();
		
		// aijk.m_ii(i, -a);
		if (!Ainf->s_ij(i, i)->einsp()) 
			return error("Plesken_decompose() Ainf->s_iji(i, i) != 1");
		for (j = i; j >= 0; j--) {
			pb = Ainf->s_ij(j, i);
			a0.mult(pb, &tmp1);
#if 1
			pa = aiaj->s_i(j);
			pa->add(&tmp1, &tmp2);
			tmp2.swap(pa);
			tmp2.freeself();
#else
			(*aiaj)[j] += tmp1;
#endif
			// c = aiaj.s_ii(j) - a * b;
			// aiaj.m_ii(j, c);
			}
		}
	if (f_v) {
		printf("aijk = ");
		aijk->println();
		}
	return OK;
}

#if TEXDOCU
INT unipoly_divide_out_oli(UNIPOLY_OP ii, SYM_OP go, VECTOR_OP stab_go, INT i)
#endif
{
	INT l, k;
	SYM_OB oli;
	SYM_OP ago;
	
	ago = stab_go->s_i(i);
	go->ganzdiv(ago, &oli);
	
	l = ii->degree();
	for (k = l; k >= 0; k--) {
		*(ii->s_i(k)) /= oli;
		}
	return OK;
}

#if TEXDOCU
INT collect_intersections_by_layer(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	MATRIX_OP Ainf, VECTOR_OP aijk, UNIPOLY_OP ii, INT f_v)
#else
Collects the intersections according to the layers of the lattice, 
e.g. according to the values in K.
#endif
{
	INT nb_dc, max_layer, k, size;
	VECTOR_OB ii_vec;
	SYM_OB tmp1, tmp2, olk;
	SYM_OP ago, pa, pc;
	
	nb_dc = Ainf->s_li();
	max_layer = K->s_ii(nb_dc - 1);
	ii_vec.m_il_n(max_layer + 1);
	
	for (k = 0; k < nb_dc; k++) {
		pa = aijk->s_i(k);
		size = K->s_ii(k);
	
		ago = stab_go->s_i(k);
		go->ganzdiv(ago, &olk);
		
		pa->mult(&olk, &tmp1);
		pc = ii_vec.s_i(size);
		pc->add(&tmp1, &tmp2);
		tmp2.swap(pc);
		}
	ii->m_v(&ii_vec);
	ii->degree();
	if (f_v) {
		printf("ii = ");
		ii->println();
		}
	return OK;
}

#endif /* LADDER_TRUE */

