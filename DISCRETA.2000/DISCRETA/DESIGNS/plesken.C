/* plesken.C */

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

#if 0
#ifdef SYM_GEO
#include <DISCRETA/geo.h>
#include <DISCRETA/geo_data.h>
#include <DISCRETA/geo_canon.h>
#endif
#endif

#if TEXDOCU
INT dc_calc_Ainf_t_via_MM(INT k_max, MATRIX_OP Ainf_t, VECTOR_OP MM, INT f_verbose)
#else
outdated, because dc\_plesken\_matrices\_prepare() is much faster. 
$MM$ is a vector of KM-matrices $M^{t,t+1}$ for $0 \le t < kmax$.
This routine computes the Plesken matrix $A^\wedge$ (transposed). 
This is an upper triangular matrix with 1s on the diagonal and with 
the KM-matrices $M_{t,k}$ as block submatrices. 
#endif
{
	MATRIX_OB A, Ainf, Ainf_block_wise;
	INT nb_d0, nb_d;
	INT nb_D, nb_D1, t_sets, k_sets;
	INT i, k, t, ii, jj;
	
	printf("dc_calc_Ainf_t_via_MM()\n");
	fflush(stdout);
	nb_d0 = dc_calc_nb_d_via_MM(MM, k_max);
	if (f_verbose) {
		printf("nb_d = %ld\n", nb_d0);
		fflush(stdout);
		}
	Ainf_t->m_ilih_n(nb_d0, nb_d0);
	
	dc_Ainf_block_wise_via_MM(0 /* k_min */, k_max, &Ainf_block_wise, MM, f_verbose);
	
	for (i = 0; i < nb_d0; i++)
		Ainf_t->m_iji(i, i, 1);
	nb_D = 0; 
		/* number of dc's before the t-sets */
	for (t = 0; t <= k_max; t++) {
		nb_d = nb_D;
			/* number of dc's before the k sets */
		for (k = t + 1; k <= k_max; k++) {
			// dc_Mtk_via_MM(t, k, &A, MM, f_verbose);
			A.freeself();
			Ainf_block_wise.s_ij(t, k)->swap(&A);
			if (f_verbose) {
				printf("M_%ld,%ld is a %ld \\times %ld matrix\n", t, k, A.s_hi(), A.s_li());
				fflush(stdout);
				// A.Print();
				}
			if (k == t + 1) {
				t_sets = A.s_hi();
				nb_D1 = nb_D + A.s_hi();
				nb_d += t_sets; /* important ! */
				}
			else {
				if (A.s_hi() != t_sets)
					return error("dc_calc_Ainf_t_via_MM() A.s_hi() != t_sets");
				}
			k_sets = A.s_li();
			for (ii = 0; ii < t_sets; ii++) {
				for (jj = 0; jj < k_sets; jj++) {
					Ainf_t->m_iji(nb_D + ii, nb_d + jj, A.s_iji(ii, jj));
					}
				}
			nb_d += k_sets;
			}
		nb_D = nb_D1;
		}
	if (f_verbose) {
		Ainf_t->latex_upper_tri(stdout);
		// Ainf_t->transpose(&Ainf);
		// Ainf.latex_lower_tri(stdout);
		fflush(stdout);
		}
	return OK;
}

#if TEXDOCU
INT dc_Mtk_via_MM(INT t, INT k, MATRIX_OP M, VECTOR_OP MM, INT f_v)
#else
outdated !
#endif
{
	MATRIX_OP pM;
	MATRIX_OB M1, M2;
	SYM_OP h;
	SYM_OB h1, h2, h3, s_ob;
	INT t1, i, j, s;
	
	printf("dc_Mtk_via_MM() t = %ld k = %ld\n", t, k);
	fflush(stdout);
#if 0
	if (t == 0)
		return error("dc_Mtk_via_MM() t == 0");
#endif
	if (k < t) 
		return error("dc_Mtk_via_MM(): k < t !");
	if (k == t) {
		INT nb_d;

		pM = (MATRIX_OP) MM->s_i(t);
		if (pM->s_obj_k() != MATRIX) {
			pM = (MATRIX_OP) MM->s_i(t - 1);
			nb_d = pM->s_li();
			}
		else
			nb_d = pM->s_hi();
		M->m_ilih_n(nb_d, nb_d);
		M->one();
		return OK;
		}
	
	MM->s_i(t)->copy((SYM_OP) M);
	printf("dc_Mtk_via_MM() M^{%ld,%ld} is a %ld x %ld matrix\n", 
		t, t + 1, M->s_hi(), M->s_li());
	fflush(stdout);
		/* M := M_t,t+1 */
	for (t1 = t + 1; t1 < k; t1++) {
		/* now M = M_t,t1 */
		MM->s_i(t1)->copy(&M1);
			/* M1 := M_t1,t1+1 */
		printf("dc_Mtk_via_MM() M is a %ld x %ld matrix\n", M->s_hi(), M->s_li());
		printf("dc_Mtk_via_MM() M^{%ld,%ld} is a %ld x %ld matrix\n", 
			t1, t1 + 1, M1.s_hi(), M1.s_li());
		fflush(stdout);
		M->mult(&M1, &M2);
		
		M2.copy(M);
			/* M := (k - t) atop 
			        (k - t1) * 
			                M_t,k   (k = t1 + 1)

			 *    = (t1 + 1 - t) atop 
			        (t1 + 1 - t1) * 
			                M_t,t1+1 

			 *    = (t1 + 1 - t) * M_t,t1+1  */
		s = t1 + 1 - t;
		for (i = 0; i < M->s_hi(); i++) {
			for (j = 0; j < M->s_li(); j++) {
				h = M->s_ij(i, j);
				s_ob.m_i_i(s);
				h->ganzdiv(&s_ob, &h1);
				// h1 = h / s;
				
				h1.mult(&s_ob, &h2);
				h2.addinvers_apply();
				h2.add(h, &h3);
				if (!h3.nullp()) {
				// if (h1 * s != h) {
					INT ii;
					MATRIX_OP A;
					
					printf("t = %ld\n", t);
					printf("k = %ld\n", k);
					printf("h = ");
						h->println();
					printf("s = ");
						s_ob.println();
					printf("h1 = ");
						h1.println();
					printf("i = %ld\n", i);
					printf("j = %ld\n", j);
					printf("M:\n");
					M->Print();
					for (ii = 0; ii < MM->s_li(); ii++) {
						A = (MATRIX_OP) MM->s_i(ii);
						printf("MM(%ld):\n", ii);
						A->Print();
						}
					return error("dc_Mtk_via_MM()|h1 * s != h");
					}
				h1.copy(M->s_ij(i, j));
				// M->m_iji(i, j, h1);
				}
			}
		}
	return OK;
}


/*
 *  working with plesken rings:
 */

#if TEXDOCU
INT design_2intersections(INT k, MATRIX_OP Ik2, MATRIX_OP Ainf, MATRIX_OP Ainf_inv, 
	VECTOR_OP orbit_length, VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, VECTOR_OP orbits, 
	INT f_multiplicities, 
	INT f_use_complement, VECTOR_OP orbits_c, 
	VECTOR_OP type2, INT f_print_block_types)
#else
Computes the essential global intersection numbers of 2 sets of the design 
and possibly of its complement. 
This is the vector $\alpha_i^{(2)}({\cal D})$ with $t < i < k$. 
The matrix Ik2 stores data which speeds up the computation a lot. 
If Ik2 is NIL, the computation is slow and uses Ainf\_inv instead.
If f\_use\_complement is TRUE, 
both vectors $\alpha_i^{(2)}({\cal D})$ and $\alpha_i^{(2)}({\cal D}^c)$ 
are collected in a vector of length 2 and returned in type.
Otherwise, only the vector $\alpha_i^{(2)}({\cal D})$ is returned in type.
type can be used as an isomorphism invariant of the design.
#endif
{
	VECTOR_OB block_types_sorted, orbit_idx;
	VECTOR_OB gl_type_2_sets, gl_type_2_sets_c;
	VECTOR_OB T;
	
	all_block_intersection_types(k, Ik2, Ainf, Ainf_inv, 
			orbit_length, K_first, K_len, 
			highest_layer, num_layers, 
			orbits, 
			&block_types_sorted, &orbit_idx, FALSE /* f_v */, FALSE /* f_vv */);

	if (f_print_block_types) {
		print_block_intersection_types(&block_types_sorted, 
			&orbit_idx, highest_layer, num_layers);
		}
	
	global_intersection_2_sets(k, 
			orbit_length, K_first, K_len, 
			highest_layer, num_layers, 
			&block_types_sorted, &orbit_idx, 
			f_multiplicities, &gl_type_2_sets, TRUE /* f_v */);

#if 0
	if (f_multiplicities) {
		VECTOR_OP gl_type, type_vec, mult_vec;
		
		gl_type = (VECTOR_OP) gl_type_2_sets.s_i(0);
		type_vec = (VECTOR_OP) gl_type_2_sets.s_i(1);
		mult_vec = (VECTOR_OP) gl_type_2_sets.s_i(2);
		printf("%% ");
		print_type_with_multiplicities(gl_type, type_vec, mult_vec);
		}
#endif
	
	if (!f_use_complement) {
		gl_type_2_sets.swap(type2);
		return OK;
		}
	
	all_block_intersection_types(k, Ik2, Ainf, Ainf_inv, 
			orbit_length, K_first, K_len, 
			highest_layer, num_layers, 
			orbits_c, 
			&block_types_sorted, &orbit_idx, FALSE /* f_v */, FALSE /* f_vv */);

	if (f_print_block_types) {
		printf("%% block intersection types of complementary design:\n");
		print_block_intersection_types(&block_types_sorted, 
			&orbit_idx, highest_layer, num_layers);
		}
	
	global_intersection_2_sets(k, 
			orbit_length, K_first, K_len, 
			highest_layer, num_layers, 
			&block_types_sorted, &orbit_idx, 
			f_multiplicities, &gl_type_2_sets_c, TRUE /* f_v */);
	
	T.m_il(2);
	gl_type_2_sets.swap((VECTOR_OP) T.s_i(0));
	gl_type_2_sets_c.swap((VECTOR_OP) T.s_i(1));
	T.swap(type2);

	return OK;
}
	
#if TEXDOCU
INT print_block_intersection_types(VECTOR_OP block_types_sorted, 
	VECTOR_OP orbit_idx, INT highest_layer, INT num_layers)
#else
Prints the block intersection types and 
the block orbit indices where they occure.
#endif
{
	VECTOR_OP type, oi;
	INT i, l, j, ll, a, tl;

	// printf("%% ");
	printf("block intersection types:");
	printf(" $(");
	for (i = 0; i < num_layers; i++) {
		printf("\\alpha_{%ld}(B_h)", highest_layer - i);
		if (i < num_layers - 1)
			printf(", ");
		}
	printf(")$ ");
	printf("\\\\\n");
	l = block_types_sorted->s_li();
	for (i = 0; i < l; i++) {
		type = (VECTOR_OP) block_types_sorted->s_i(i);
		oi = (VECTOR_OP) orbit_idx->s_i(i);
		ll = oi->s_li();
		// printf("%% ");
		tl = type->s_li();
		if (tl == 1) {
			printf("$\\alpha_{%ld}(B_h) = ", highest_layer);
			type->s_i(0)->print();
			printf("$ ");
			}
		else
			type->print();
		printf("for $h \\in \\{$");
		for (j = 0; j < ll; j++) {
			a = oi->s_ii(j);
			printf("%ld", a + 1);
			if (j < ll - 1)
				printf(", ");
			}
		printf("$\\}$");
		if (i < l - 1)
			printf(", \n");
		}
	printf(". \\\\\n");
	printf("%% end of %ld types\n", l);
	return OK;
}

#if TEXDOCU
INT print_gl_type2_intersections(VECTOR_OP type, INT f_multiplicities, INT f_complements)
#endif
{
	VECTOR_OP type1;
	
	if (f_complements) {
		type1 = (VECTOR_OP) type->s_i(0);
		}
	else {
		type1 = type;
		}
	if (f_multiplicities) {
		VECTOR_OP gl_type, type_vec, mult_vec;
		
		gl_type = (VECTOR_OP) type1->s_i(0);
		type_vec = (VECTOR_OP) type1->s_i(1);
		mult_vec = (VECTOR_OP) type1->s_i(2);
		// printf("%% ");
		printf("$\\alpha^{(2)}({\\cal D}):$ ");
		print_type_with_multiplicities(gl_type, type_vec, mult_vec);
		}
	else {
		// printf("%% ");
		printf("$\\alpha^{(2)}({\\cal D}):$ ");
		type1->latex(stdout);
		printf("\n");
		}
	return OK;
}

#if TEXDOCU
INT print_type_with_multiplicities(VECTOR_OP gl_type, 
	VECTOR_OP type_vec, VECTOR_OP mult_vec)
#endif
{
	INT type_len, i, l;
	SYM_OP y;
	VECTOR_OP x;
	VECTOR_OB mult_vec1;
	SYM_OB g;

	type_len = gl_type->s_li();
	printf("$");
	if (type_len == 1) {
		gl_type->s_i(0)->latex(stdout);
		}
	else {
		gl_type->latex(stdout);
		}
	printf(" = \\frac{1}{2} ");
	
	mult_vec->gcd_all_elements_and_divide_out(&g, &mult_vec1);
	if (!g.einsp()) {
		g.latex(stdout);
		printf(" \\times (");
		}
	
	l = type_vec->s_li();
	for (i = 0; i < l; i++) {
		x = (VECTOR_OP) type_vec->s_i(i);
		if (x->s_li() != type_len)
			return error("print_type_with_multiplicities() x->s_li() != type_len");
		y = mult_vec1.s_i(i);
		y->latex(stdout);
		printf(" \\times ");
		
		if (type_len == 1) {
			x->s_i(0)->latex(stdout);
			}
		else {
			x->latex(stdout);
			}
			
		if (i < l - 1)
			printf(" + ");
		}
	if (!g.einsp()) {
		printf(")");
		}
	printf("$\n");
	return OK;
}

#if TEXDOCU
INT print_invariants_and_classes(VECTOR_OP types_sorted, VECTOR_OP classes, 
	INT f_multiplicities, INT f_complements)
#else
prints the types (of types\_sorted) together with the class 
(the elements where they occure). classes is a vector of vectors.
#endif
{
	VECTOR_OB class_size;
	VECTOR_OP type, type1, cl, gl_type, type_vec, mult_vec;
	INT i, l, j, ll, a, max_class_size, n;
	double x_quer, sigma, varianz, x;

	printf("%% types:\\\\\n");
	l = types_sorted->s_li();
	n = 0;
	max_class_size = 0;
	for (i = 0; i < l; i++) {
		cl = (VECTOR_OP) classes->s_i(i);
		a = cl->s_li();
		max_class_size = MAXIMUM(a, max_class_size);
		n += a;
		}
	class_size.m_il_n(max_class_size + 1);
	for (i = 0; i < l; i++) {
		cl = (VECTOR_OP) classes->s_i(i);
		a = cl->s_li();
		class_size.s_i(a)->inc();
		}
	printf("n=%ld objects in %ld classes\\\\\n", n, l);
	for (i = 0; i <= max_class_size; i++) {
		a = class_size.s_ii(i);
		printf("size %ld: %ld classes\\\\\n", i, a);
		}
	x_quer = (double) n / (double) l;
	printf("$\\overline{x} = %lf$\\\\\n", x_quer);
	sigma = 0.;
	for (i = 0; i <= max_class_size; i++) {
		a = class_size.s_ii(i);
		x = (double)i - x_quer;
		x = x * x;
		x = x * (double) a;
		sigma += x;
		}
	varianz = sigma / (double) (l - 1);
	sigma = sqrt(varianz);
	printf("$v = %lf, \\sigma=%lf$\\\\\n", varianz, sigma);
	
	for (i = 0; i < l; i++) {
		type = (VECTOR_OP) types_sorted->s_i(i);
		cl = (VECTOR_OP) classes->s_i(i);
		ll = cl->s_li();
		printf("%% ");
		if (f_multiplicities) {
			if (f_complements) {
				type1 = (VECTOR_OP) type->s_i(0);
				gl_type = (VECTOR_OP) type1->s_i(0);
				type_vec = (VECTOR_OP) type1->s_i(1);
				mult_vec = (VECTOR_OP) type1->s_i(2);
				}
			else {
				gl_type = (VECTOR_OP) type->s_i(0);
				type_vec = (VECTOR_OP) type->s_i(1);
				mult_vec = (VECTOR_OP) type->s_i(2);
				}
			print_type_with_multiplicities(gl_type, 
				type_vec, mult_vec);

#if 0
			gl_type->latex(stdout);
			printf(" = ");
			lll = type_vec->s_li();
			for (ii = 0; ii < lll; ii++) {
				x = type_vec->s_i(ii);
				y = mult_vec->s_i(ii);
				x->latex(stdout);
				printf(" \\times ");
				y->latex(stdout);
				if (ii < lll - 1)
					printf(" + ");
				}
#endif
			}
		else {
			printf("$ ");
			type->latex(stdout);
			printf("$ \n");
			}
		printf("%% for $\\{");
		for (j = 0; j < ll; j++) {
			a = cl->s_ii(j);
			printf("%ld", a);
			if (j < ll - 1)
				printf(", ");
			}
		printf("\\}$\n");
		}
	return OK;
}

#define F_CHECK_MULTIPLICITIES 

#if TEXDOCU
INT global_intersection_2_sets(INT k_max, 
	VECTOR_OP orbit_length, VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	VECTOR_OP block_types_sorted, VECTOR_OP orbit_idx, 
	INT f_multiplicities, VECTOR_OP gl_type_2_sets, INT f_v)
#else
Computes $\alpha_i^{(2)}({\cal D})$ via the well known formula
$\alpha_i^{(2)}({\cal D}) = 1/2 \cdot \sum_{j=1}^l |\tilde{K_j}| \cdot \alpha_i^{(2)}(K_j)$
where $K_1, \ldots, K_l$ are representing sets for the orbits of the design, 
$\tilde{K_j}$ is the corresponding $G$-orbit.
type contains the $\alpha_i^{(2)}(K_j)$ and orbit\_idx the orbital indices 
$j$ for which this value occurs.
#endif
{
	VECTOR_OP type, oi;
	VECTOR_OB t;
	SYM_OP ol;
	SYM_OB two_ob;
	VECTOR_OB block_types_multiplicities;
	SYM_OP x;
	SYM_OB y;
	INT i, l, j, ll, a, first_k_orbit;

	printf("%% global_intersection_2_sets()\n");
	two_ob.m_i_i(2);
	first_k_orbit = K_first->s_ii(k_max);
	l = block_types_sorted->s_li();
	gl_type_2_sets->m_il_n(num_layers);
	if (f_multiplicities) {
		block_types_multiplicities.m_il_n(l);
		}
	for (i = 0; i < l; i++) {
		type = (VECTOR_OP) block_types_sorted->s_i(i);
		oi = (VECTOR_OP) orbit_idx->s_i(i);
		ll = oi->s_li();
		for (j = 0; j < ll; j++) {
			a = oi->s_ii(j);
			ol = orbit_length->s_i(first_k_orbit + a);
			type->copy(&t);
			t.multiply_elementwise(ol);
			gl_type_2_sets->add_apply_elementwise(&t, +1 /* sign */);
			if (f_multiplicities) {
				x = block_types_multiplicities.s_i(i);
				x->add(ol, &y);
				y.swap(x);
				}
			}
		}
	gl_type_2_sets->divide_out(&two_ob);
	if (f_v) {
		printf("%% gl_type_2_sets: ");
		gl_type_2_sets->println();
		}
	if (f_multiplicities) {
		if (f_v && 0) {
			printf("%% multiplicities: \n");
			printf("%% $");
			for (i = 0; i < l; i++) {
				type = (VECTOR_OP) block_types_sorted->s_i(i);
				x = block_types_multiplicities.s_i(i);
				type->latex(stdout);
				printf(" \\times ");
				x->latex(stdout);
				if (i < l - 1)
					printf(", ");
				}
			printf("$\n");
			}
#ifdef F_CHECK_MULTIPLICITIES
#endif
		y.freeself();
		gl_type_2_sets->swap(&y);
		gl_type_2_sets->m_il(3);
		y.swap(gl_type_2_sets->s_i(0));
		block_types_sorted->copy((VECTOR_OP) gl_type_2_sets->s_i(1));
		block_types_multiplicities.copy((VECTOR_OP) gl_type_2_sets->s_i(2));
		}
	return OK;
}

#if TEXDOCU
INT all_block_intersection_types(INT k_max, 
	MATRIX_OP Ik2, MATRIX_OP Ainf, MATRIX_OP Ainf_inv, 
	VECTOR_OP orbit_length, VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	VECTOR_OP design_orbits, 
	VECTOR_OP block_types_sorted, VECTOR_OP orbit_idx, INT f_v, INT f_vv)
#else
Computes all $\alpha_i^{(2)}(K_j)$ of orbits of the design. 
The numbers of orbits in the design are listed in design\_orbits.
The output is block\_types\_sorted and orbit\_idx.
#endif
{
	INTEGER_OB len;
	VECTOR_OB type, empty_vec;
	VECTOR_OP oi;
	INT do_idx, l, idx, f_found;
	
	l = design_orbits->s_li();
	block_types_sorted->m_il(l);
	orbit_idx->m_il(l);
	empty_vec.m_il(0);
	len.m_i(0);
	printf("%%");
	printf("computing block types of %ld orbits", l);
	if (Ik2) 
		printf(" (using Ik2 matrix)");
	printf(":\n");
	printf("%%");
	fflush(stdout);
	for (do_idx = 0; do_idx < l; do_idx++) {
		if (Ik2) {
			get_block_intersection_type(Ik2, k_max, 
				design_orbits, do_idx, orbit_length, K_first, K_len, 
				highest_layer, num_layers, &type);
			}
		else {
			block_intersection_type(k_max, Ainf, Ainf_inv, 
				orbit_length, K_first, K_len, 
				highest_layer, num_layers, 
				do_idx, design_orbits, 
				&type, FALSE /* f_v */, FALSE /* f_vv */);
			}
		block_types_sorted->search(len.s_i(), TRUE, &type, &idx, &f_found);
		if (f_found) {
			idx--;
			}
		else {
			block_types_sorted->insert_at(&len, idx, &type);
			len.dec();
			orbit_idx->insert_at(&len, idx, &empty_vec);
			}
		oi = (VECTOR_OP) orbit_idx->s_i(idx);
		oi->inc();
		oi->m_ii(oi->s_li() - 1, design_orbits->s_ii(do_idx));
		if ((do_idx + 1) % 10 == 0)
			printf(",");
		else
			printf(".");
		fflush(stdout);
		}
	block_types_sorted->realloc_z(len.s_i());
	orbit_idx->realloc_z(len.s_i());
	printf("\n");
	return OK;
}

#if TEXDOCU
INT block_intersection_type(INT k_max, MATRIX_OP Ainf, MATRIX_OP Ainf_inv, 
	VECTOR_OP orbit_length, VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	INT do_idx, VECTOR_OP design_orbits, 
	VECTOR_OP type, INT f_v, INT f_vv)
#else
Computes $\alpha_i^{(2)}(K_j)$ where $j$ is design\_orbits$->$s\_ii(do\_idx).
That means, this routine computes the essential block intersection type 
of the design orbit do\_idx.
#endif
{
	INT i, j, oi, oj, l, pos_layer_k, first_k_orbit;
	VECTOR_OB args, t, t0;
	SYM_OP ol;

	l = design_orbits->s_li();
	i = do_idx;
	args.m_il(2);
	oi = design_orbits->s_ii(do_idx);
	first_k_orbit = K_first->s_ii(k_max);
	args.m_ii(0, first_k_orbit + oi);
	for (j = 0; j < l; j++) {
		oj = design_orbits->s_ii(j);
		args.m_ii(1, first_k_orbit + oj);
		if (f_vv) {
			printf("i=%ld j=%ld\n", i, j);
			}
		plesken_product(Ainf, Ainf_inv, orbit_length, &args, 
			K_first, K_len, highest_layer, num_layers, 0 /* k_min */, 
			&t0, f_vv, f_vv);
		if (i == j) {
			pos_layer_k = highest_layer - k_max;
			if (pos_layer_k >= 0 && pos_layer_k < num_layers) {
				t0.m_ii(pos_layer_k, 0);
				}
			}
		if (j == 0) {
			t0.swap(&t);
			}
		else {
			t.add_apply_elementwise(&t0, +1 /* sign */);
			}
		}
	ol = orbit_length->s_i(first_k_orbit + oi);
	t.divide_out(ol);
	if (f_v) {
		printf("block type orbit %ld: ", oi);
		t.println();
		}
	t.swap(type);
	return OK;
}




#endif /* LADDER_TRUE */


