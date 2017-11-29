/* plesken_Ik2.C */

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
INT calc_intersections_Iknm(MATRIX_OP Ainf_block_wise, MATRIX_OP coeff, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	MATRIX_OP I, INT k_min, INT k, 
	VECTOR_OP Sel1, VECTOR_OP Sel2, INT f_v, INT f_vv)
#else
Computes the Iknm matrix useful for determining intersections of pairs of orbits 
(the orbit numbers are listed in Sel1 and Sel2).
#endif
{
	INT t0, t1, user_time;
	BYTE str[256];
	INT i, ii, j, jj, N, k_first, k_len, r, l;
	VECTOR_OB args, type, I_vec;
	SYM_OP ol;
	INT nb_sel1 = 0;
	INT nb_sel2 = 0;
	
	t0 = os_ticks();
	k_first = K_first->s_ii(k);
	k_len = K_len->s_ii(k);
	if (f_v) {
		printf("calc_intersections_Iknm() k_min = %ld k = %ld k_first = %ld k_len = %ld\n", k_min, k, k_first, k_len);
		printf("highest_layer = %ld num_layers = %ld\n", highest_layer, num_layers);
		fflush(stdout);
		}
	if (highest_layer - (num_layers  - 1) < k_min)
		return error("calc_intersections_Iknm() highest_layer - (num_layers - 1) < k_min");
	nb_sel1 = Sel1->s_li();
	nb_sel2 = Sel2->s_li();
	printf("Sel1 contains %ld elements\n", nb_sel1);
	printf("Sel2 contains %ld elements\n", nb_sel2);
	N = nb_sel1 * nb_sel2;
	I->m_ilih_n(num_layers, N);
	args.m_il(2);
	l = 0;
	for (i = 0; i < nb_sel1; i++) {
		ii = k_first + Sel1->s_ii(i);
		printf("i=%ld,block_orbit=%ld:\n", i, ii);
		fflush(stdout);
		for (j = 0; j < nb_sel2; j++) {
			jj = k_first + Sel2->s_ii(j);
			args.m_ii(0, ii);
			args.m_ii(1, jj);
			plesken_product_block_wise(Ainf_block_wise, coeff, 
				orbit_length, &args, K_first, K_len, 
				highest_layer, num_layers, k_min, 
				&type, TRUE /* f_v */, FALSE /* f_vv */);
			if (type.s_li() != num_layers)
				return error("calc_intersections_Iknm() type.s_li() != num_layers");
			ol = orbit_length->s_i(ii);
			type.divide_out(ol);
			for (r = 0; r < num_layers; r++) {
				type.s_i(r)->copy(I->s_ij(l, r));
				}
			l++;
			if ((j + 1) % 10 == 0) {
				printf(",");
				if ((j + 1) % 50 == 0)
					printf(" %ld\n", j + 1);
				}
			else
				printf(".");
			fflush(stdout);
			}
		printf("\n");
		}
	printf(";\n");
	fflush(stdout);
	if (f_v) {
		printf("calc_intersections_Iknm() \n");

		intersections_Iknm_blow_up(I, Sel1, Sel2, 
			k, &I_vec, orbit_length, K_first, K_len, 
			highest_layer, num_layers, TRUE /* f_v */);
		fflush(stdout);
		}
	t1 = os_ticks();
	
	user_time = t1 - t0;
	strcpy(str, "Running time for calc_intersections_Iknm(): ");
	print_delta_time(user_time, str);
	printf("%s\n", str);

	return OK;
}

#if TEXDOCU
INT calc_intersections_Iknn(MATRIX_OP Ainf, MATRIX_OP coeff, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	MATRIX_OP I, INT k_min, INT k, 
	VECTOR_OP Sel, INT f_v, INT f_vv)
#else
Computes the Iknn matrix useful for determining intersections of pairs of orbits 
(the orbit numbers are listed in Sel).
#endif
{
	INT t0, t1, user_time;
	BYTE str[256];
	INT i, ii, j, N, k_first, k_len, r;
	INT i1, i2, ii1, ii2;
	VECTOR_OB args, type, I_vec;
	SYM_OP ol;
	INT nb_sel = 0;
	
	t0 = os_ticks();
	k_first = K_first->s_ii(k);
	k_len = K_len->s_ii(k);
	if (f_v) {
		printf("calc_intersections_Iknn() k_min = %ld k = %ld k_first = %ld k_len = %ld\n", k_min, k, k_first, k_len);
		printf("highest_layer = %ld num_layers = %ld\n", highest_layer, num_layers);
		fflush(stdout);
		}
	if (highest_layer - (num_layers  - 1) < k_min)
		return error("calc_intersections_Iknn() highest_layer - (num_layers - 1) < k_min");
	nb_sel = Sel->s_li();
	printf("Sel contains %ld elements\n", nb_sel);
	N = nb_sel + ((nb_sel * (nb_sel - 1)) >> 1);
	I->m_ilih_n(num_layers, N);
	args.m_il(2);
	for (i = 0; i < nb_sel; i++) {
		ii = k_first + Sel->s_ii(i);
		args.m_ii(0, ii);
		args.m_ii(1, ii);
		plesken_product(Ainf, coeff, 
			orbit_length, &args, K_first, K_len, 
			highest_layer, num_layers, k_min, 
			&type, FALSE /* f_v */, FALSE /* f_vv */);
		if (type.s_li() != num_layers)
			return error("calc_intersections_Iknn() type.s_li() != num_layers");
		ol = orbit_length->s_i(ii);
		type.divide_out(ol);
		for (r = 0; r < num_layers; r++) {
			type.s_i(r)->copy(I->s_ij(i, r));
			}
		if ((i + 1) % 10 == 0) {
			printf(",");
			if ((i + 1) % 50 == 0)
				printf(" %ld\n", i + 1);
			}
		else
			printf(".");
		fflush(stdout);
		}
	printf(";\n");
	fflush(stdout);
	j = 0;
	for (i1 = 0; i1 < nb_sel; i1++) {
		ii1 = k_first + Sel->s_ii(i1);
		args.m_ii(0, ii1);
		for (i2 = i1 + 1; i2 < nb_sel; i2++) {
			if (ij2k(i1, i2, nb_sel) != j)
				return error("calc_intersections_Iknn() ij2k(i1, i2, nb_sel) != j");
			ii2 = k_first + Sel->s_ii(i2);
			args.m_ii(1, ii2);
			plesken_product(Ainf, coeff, 
				orbit_length, &args, K_first, K_len, 
				highest_layer, num_layers, k_min, 
				&type, FALSE /* f_v */, FALSE /* f_vv */);
			if (type.s_li() != num_layers)
				return error("calc_intersections_Iknn() type.s_li() != num_layers");
			ol = orbit_length->s_i(ii1);
			type.divide_out(ol);
			for (i = 0; i < num_layers; i++) {
				type.s_i(i)->copy(I->s_ij(nb_sel + j, i));
				}
			j++;
			}
		if ((i1 + 1) % 10 == 0) {
			printf(",");
			if ((i1 + 1) % 50 == 0)
				printf(" %ld\n", i1 + 1);
			}
		else
			printf(".");
		fflush(stdout);
		}
	printf(";\n");
	fflush(stdout);
	if (f_v) {
		printf("calc_intersections_Iknn() \n");

		intersections_Iknn_blow_up(I, Sel, 
			k, &I_vec, orbit_length, K_first, K_len, 
			highest_layer, num_layers, TRUE /* f_v */);
		fflush(stdout);
		}
	t1 = os_ticks();
	
	user_time = t1 - t0;
	strcpy(str, "Running time for calc_intersections_Iknn(): ");
	print_delta_time(user_time, str);
	printf("%s\n", str);

	return OK;
}

/* 
 * the matrix Ik2 stores the plesken ring information
 * for products of two elements
 */

#if TEXDOCU
INT calc_intersections_Ik2(MATRIX_OP Ainf, MATRIX_OP coeff, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, INT k_min, 
	MATRIX_OP Ik2, INT k, INT f_v, INT f_vv)
#else
Computes the Ik2 matrix useful for determining intersections of pairs of orbits.
#endif
{
	INT t0, t1, user_time;
	BYTE str[256];
	INT i, ii, j, N, k_first, k_len, r;
	INT i1, i2, ii1, ii2;
	VECTOR_OB args, type, Ik2_vec;
	SYM_OP ol;
	
	t0 = os_ticks();
	k_first = K_first->s_ii(k);
	k_len = K_len->s_ii(k);
	if (f_v) {
		printf("calc_intersections_Ik2() k = %ld k_first = %ld k_len = %ld\n", k, k_first, k_len);
		printf("highest_layer = %ld num_layers = %ld\n", highest_layer, num_layers);
		fflush(stdout);
		}
	N = k_len + ((k_len * (k_len - 1)) >> 1);
	Ik2->m_ilih_n(num_layers, N);
	args.m_il(2);
	for (i = 0; i < k_len; i++) {
		ii = k_first + i;
		args.m_ii(0, ii);
		args.m_ii(1, ii);
		plesken_product(Ainf, coeff, 
			orbit_length, &args, K_first, K_len, 
			highest_layer, num_layers, 0 /* k_min */, 
			&type, FALSE /* f_v */, FALSE /* f_vv */);
		if (type.s_li() != num_layers)
			return error("calc_intersections_Ik2() type.s_li() != num_layers");
		ol = orbit_length->s_i(ii);
		type.divide_out(ol);
		for (r = 0; r < num_layers; r++) {
			type.s_i(r)->copy(Ik2->s_ij(i, r));
			}
		if ((i + 1) % 10 == 0) {
			printf(",");
			if ((i + 1) % 50 == 0)
				printf(" %ld\n", i + 1);
			}
		else
			printf(".");
		fflush(stdout);
		}
	printf(";\n");
	fflush(stdout);
	j = 0;
	for (i1 = 0; i1 < k_len; i1++) {
		ii1 = k_first + i1;
		args.m_ii(0, ii1);
		for (i2 = i1 + 1; i2 < k_len; i2++) {
			if (ij2k(i1, i2, k_len) != j)
				return error("calc_intersections_Ik2() ij2k(i1, i2, k_len) != j");
			ii2 = k_first + i2;
			args.m_ii(1, ii2);
			plesken_product(Ainf, coeff, 
				orbit_length, &args, K_first, K_len, 
				highest_layer, num_layers, 0 /* k_min */, 
				&type, FALSE /* f_v */, FALSE /* f_vv */);
			if (type.s_li() != num_layers)
				return error("calc_intersections_Ik2() type.s_li() != num_layers");
			ol = orbit_length->s_i(ii1);
			type.divide_out(ol);
			for (i = 0; i < num_layers; i++) {
				type.s_i(i)->copy(Ik2->s_ij(k_len + j, i));
				}
			j++;
			}
		if ((i1 + 1) % 10 == 0) {
			printf(",");
			if ((i1 + 1) % 50 == 0)
				printf(" %ld\n", i1 + 1);
			}
		else
			printf(".");
		fflush(stdout);
		}
	printf(";\n");
	fflush(stdout);
	if (f_v) {
		printf("calc_intersections_Ik2() \n");
		intersections_Ik2_blow_up(Ik2, k, &Ik2_vec, 
			orbit_length, K_first, K_len, 
			highest_layer, num_layers, TRUE /* f_v */);
		// Ik2->Print();
		fflush(stdout);
		}
	t1 = os_ticks();
	
	user_time = t1 - t0;
	strcpy(str, "Running time for calc_intersections_Ik2(): ");
	print_delta_time(user_time, str);
	printf("%s\n", str);

	return OK;
}

#if TEXDOCU
INT get_block_intersection_type(MATRIX_OP Ik2, INT k, 
	VECTOR_OP design_orbits, INT do_idx, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	VECTOR_OP type)
#endif
{
	INT i, l, do1, do2;
	VECTOR_OB t0;
	
	type->m_il_n(num_layers);
	l = design_orbits->s_li();
	do1 = design_orbits->s_ii(do_idx);
	for (i = 0; i < l; i++) {
		do2 = design_orbits->s_ii(i);
		get_intersections_of_2(Ik2, k, do1, do2, 
			orbit_length, K_first, K_len, 
			highest_layer, num_layers, 
			&t0);
		type->add_apply_elementwise(&t0, +1);
		}
	return OK;
}


#if TEXDOCU
INT intersections_Iknm_blow_up(MATRIX_OP I, VECTOR_OP Sel1, VECTOR_OP Sel2, 
	INT k, VECTOR_OP I_vec, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, INT f_v)
#else
Decodes the Iknm matrix into a vector I\_vec of matrices.
#endif
{
	INT i, i1, i2, k_first, nb_sel1, nb_sel2;
	MATRIX_OB I1;
	MATRIX_OP M;
	VECTOR_OB type;
	
	if (f_v) {
		printf("blowing up Iknm\n");
		fflush(stdout);
		}
	nb_sel1 = Sel1->s_li();
	nb_sel2 = Sel2->s_li();
	k_first = K_first->s_ii(k);
	I_vec->m_il(num_layers);
	for (i = 0; i < num_layers; i++) {
		I1.m_ilih_n(nb_sel2, nb_sel1);
		I1.swap((MATRIX_OP) I_vec->s_i(i));
		}
	for (i1 = 0; i1 < nb_sel1; i1++) {
		for (i2 = 0; i2 < nb_sel2; i2++) {
			Iknm_get_intersections_of_2(I, Sel1, Sel2, i1, i2, 
				num_layers, &type);
			for (i = 0; i < num_layers; i++) {
				M = (MATRIX_OP) I_vec->s_i(i);
				type.s_i(i)->swap(M->s_ij(i1, i2));
				}
			}
		}
	if (f_v) {
		printf("%%Iknm matrix, k=%ld highest_layer=%ld num_layers=%ld\n", 
			k, highest_layer, num_layers);
		for (i = 0; i < num_layers; i++) {
			printf("layer=%ld\n", highest_layer - i);
			M = (MATRIX_OP) I_vec->s_i(i);
			M->Print();
			}
		}
	return OK;
}

#if TEXDOCU
INT intersections_Iknn_blow_up(MATRIX_OP I, VECTOR_OP Sel, 
	INT k, VECTOR_OP I_vec, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, INT f_v)
#else
Decodes the Iknn matrix into a vector I\_vec of matrices.
#endif
{
	INT i, i1, i2, k_first, nb_sel;
	MATRIX_OB I1;
	MATRIX_OP M;
	VECTOR_OB type;
	
	if (f_v) {
		printf("blowing up Iknn\n");
		fflush(stdout);
		}
	nb_sel = Sel->s_li();
	k_first = K_first->s_ii(k);
	I_vec->m_il(num_layers);
	for (i = 0; i < num_layers; i++) {
		I1.m_ilih_n(nb_sel, nb_sel);
		I1.swap((MATRIX_OP) I_vec->s_i(i));
		}
	for (i1 = 0; i1 < nb_sel; i1++) {
		for (i2 = 0; i2 < nb_sel; i2++) {
			Iknn_get_intersections_of_2(I, Sel, k, i1, i2, 
				orbit_length, K_first, K_len, 
				highest_layer, num_layers, 
				&type);
			for (i = 0; i < num_layers; i++) {
				M = (MATRIX_OP) I_vec->s_i(i);
				type.s_i(i)->swap(M->s_ij(i1, i2));
				}
			}
		}
	if (f_v) {
		printf("%%Iknn matrix, k=%ld highest_layer=%ld num_layers=%ld\n", 
			k, highest_layer, num_layers);
		for (i = 0; i < num_layers; i++) {
			printf("layer=%ld\n", highest_layer - i);
			M = (MATRIX_OP) I_vec->s_i(i);
			M->Print();
			}
		}
	return OK;
}

#if TEXDOCU
INT Iknm_get_intersections_of_2(MATRIX_OP I, VECTOR_OP Sel1, VECTOR_OP Sel2, 
	INT i1, INT i2, 
	INT num_layers, 
	VECTOR_OP type)
#else
reads out the intersection type between orbits Sel[i1] and Sel[i2] from I (which is a Iknn). 
#endif
{
	INT j, r;
	INT nb_sel1, nb_sel2;
	
	nb_sel1 = Sel1->s_li();
	nb_sel2 = Sel2->s_li();
	type->m_il(num_layers);
	
	j = i1 * nb_sel2 + i2;
	for (r = 0; r < num_layers; r++) 
		I->s_ij(j, r)->copy(type->s_i(r));
	return OK;
}

#if TEXDOCU
INT Iknn_get_intersections_of_2(MATRIX_OP I, VECTOR_OP Sel, 
	INT k, INT i1, INT i2, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	VECTOR_OP type)
#else
reads out the intersection type between orbits Sel[i1] and Sel[i2] from I (which is a Iknn). 
Only the essential intersections are stored, used for computing 
$\alpha_i$ with $t < i < k$.
#endif
{
	SYM_OP ol;
	INT k_first, j, r;
	INT nb_sel, ii1, ii2;
	
	k_first = K_first->s_ii(k);
	nb_sel = Sel->s_li();
	type->m_il(num_layers);
	if (i1 == i2) {
		j = i1;
		for (r = 0; r < num_layers; r++) 
			I->s_ij(j, r)->copy(type->s_i(r));
		return OK;
		}
	
	j = nb_sel + ij2k(i1, i2, nb_sel);
	for (r = 0; r < num_layers; r++) 
		I->s_ij(j, r)->copy(type->s_i(r));
	
	if (i1 < i2) {
		return OK;
		}
	
	ii1 = Sel->s_ii(i1);
	ii2 = Sel->s_ii(i2);
	ol = orbit_length->s_i(k_first + ii2);
	type->multiply_elementwise(ol);
	ol = orbit_length->s_i(k_first + ii1);
	type->divide_out(ol);
	return OK;
}

#if TEXDOCU
INT intersections_Ik2_blow_up(MATRIX_OP Ik2, INT k, VECTOR_OP Ik2_vec, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, INT f_v)
#else
Decodes the Ik2 matrix into a vector Ik2\_vec of matrices.
#endif
{
	INT i, i1, i2, k_first, k_len;
	MATRIX_OB Ik2_i;
	MATRIX_OP M;
	VECTOR_OB type;
	
	if (f_v) {
		printf("blowing up Ik2\n");
		fflush(stdout);
		}
	k_first = K_first->s_ii(k);
	k_len = K_len->s_ii(k);
	Ik2_vec->m_il(num_layers);
	for (i = 0; i < num_layers; i++) {
		Ik2_i.m_ilih_n(k_len, k_len);
		Ik2_i.swap((MATRIX_OP) Ik2_vec->s_i(i));
		}
	for (i1 = 0; i1 < k_len; i1++) {
		for (i2 = 0; i2 < k_len; i2++) {
			get_intersections_of_2(Ik2, k, i1, i2, 
				orbit_length, K_first, K_len, 
				highest_layer, num_layers, 
				&type);
			for (i = 0; i < num_layers; i++) {
				M = (MATRIX_OP) Ik2_vec->s_i(i);
				type.s_i(i)->swap(M->s_ij(i1, i2));
				}
			}
		}
	if (f_v) {
		printf("%%Ik2 matrix, k=%ld highest_layer=%ld num_layers=%ld\n", 
			k, highest_layer, num_layers);
		for (i = 0; i < num_layers; i++) {
			printf("layer=%ld\n", highest_layer - i);
			M = (MATRIX_OP) Ik2_vec->s_i(i);
			M->Print();
			}
		}
	return OK;
}

#if TEXDOCU
INT get_intersections_of_2(MATRIX_OP Ik2, INT k, INT i1, INT i2, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	VECTOR_OP type)
#else
reads out the intersection type between orbits i1 and i2 from Ik2. 
Only the essential intersections are stored, used for computing 
$\alpha_i$ with $t < i < k$.
#endif
{
	SYM_OP ol;
	INT k_first, k_len, j, r;
	
	k_first = K_first->s_ii(k);
	k_len = K_len->s_ii(k);
	type->m_il(num_layers);
	if (i1 == i2) {
		j = i1;
		for (r = 0; r < num_layers; r++) 
			Ik2->s_ij(j, r)->copy(type->s_i(r));
		return OK;
		}
	
	j = k_len + ij2k(i1, i2, k_len);
	for (r = 0; r < num_layers; r++) 
		Ik2->s_ij(j, r)->copy(type->s_i(r));
	
	if (i1 < i2) {
		return OK;
		}
	
	ol = orbit_length->s_i(k_first + i2);
	type->multiply_elementwise(ol);
	ol = orbit_length->s_i(k_first + i1);
	type->divide_out(ol);
	return OK;
}

#if TEXDOCU
INT Ik2_print(MATRIX_OP Ik2, INT k, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, INT f_v)
#endif
{
	INT i, j, r, i1, i2, k_first, k_len;
	
	k_first = K_first->s_ii(k);
	k_len = K_len->s_ii(k);
	printf("Ik2 print, k =%ld k_first = %ld k_len = %ld\n", k, k_first, k_len);
	for (i = 0; i < k_len; i++) {
		printf("i=%ld (with itself): ", i);
		for (r = 0; r < num_layers; r++) {
			Ik2->s_ij(i, r)->print();
			if (r < num_layers - 1)
				printf(", ");
			}
		printf("\n");
		}
	for (i1 = 0; i1 < k_len; i1++) {
		for (i2 = 0; i2 < k_len; i2++) {
			j = k_len + ij2k(i1, i2, k_len);
			printf("i1=%ld i2=%ld: ", i1, i2);
			for (r = 0; r < num_layers; r++) {
				Ik2->s_ij(j, r)->print();
				if (r < num_layers - 1)
					printf(", ");
				}
			printf("\n");
			}
		}
	return OK;
}

/* 
 * the matrix Ik3 stores the plesken ring information
 * for products of three elements
 */

/* 
 * Ik3
 * contains data for intersections of all 3-sets of orbits 
 */

#if TEXDOCU
INT calc_intersections_Ik3(MATRIX_OP Ainf, MATRIX_OP coeff, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	MATRIX_OP Ik3, INT k, INT f_v, INT f_vv)
#else
Computes the Ik3 matrix useful for determining intersections of 3-sets of orbits.
#endif
{
	INT t0, t1, user_time;
	BYTE str[256];
	INT i, ii, j, N, k_first, k_len, r;
	INT i1, i2, i3, ii1, ii2, ii3;
	INT k_len2, k_len_atop_3;
	VECTOR_OB args, type, Ik3_vec;
	SYM_OP ol;
	
	t0 = os_ticks();
	k_first = K_first->s_ii(k);
	k_len = K_len->s_ii(k);
	if (f_v) {
		printf("calc_intersections_Ik3() k = %ld k_first = %ld k_len = %ld\n", k, k_first, k_len);
		printf("highest_layer = %ld num_layers = %ld\n", highest_layer, num_layers);
		fflush(stdout);
		}
	k_len2 = k_len * k_len;
	k_len_atop_3 = (((k_len * (k_len - 1)) >> 1) * (k_len - 2) ) / 3;
	N = k_len + k_len2 + k_len_atop_3;
	Ik3->m_ilih_n(num_layers, N);
	args.m_il(3);
	for (i = 0; i < k_len; i++) {
		ii = k_first + i;
		args.m_ii(0, ii);
		args.m_ii(1, ii);
		args.m_ii(2, ii);
		plesken_product(Ainf, coeff, 
			orbit_length, &args, K_first, K_len, 
			highest_layer, num_layers, 0 /* k_min */, 
			&type, FALSE /* f_v */, FALSE /* f_vv */);
		if (type.s_li() != num_layers)
			return error("calc_intersections_Ik3() type.s_li() != num_layers");
		ol = orbit_length->s_i(ii);
		type.divide_out(ol);
		for (r = 0; r < num_layers; r++) {
			type.s_i(r)->copy(Ik3->s_ij(i, r));
			}
		if ((i + 1) % 10 == 0) {
			printf(",");
			if ((i + 1) % 50 == 0)
				printf(" %ld\n", i + 1);
			}
		else
			printf(".");
		fflush(stdout);
		}
	printf(";\n");
	fflush(stdout);
	
	j = 0;
	for (i1 = 0; i1 < k_len; i1++) {
		ii1 = k_first + i1;
		args.m_ii(0, ii1);
		args.m_ii(1, ii1);
		for (i2 = 0; i2 < k_len; i2++) {
			if ((i1 * k_len + i2) != j)
				return error("calc_intersections_Ik3() (i1 * k_len + i2) != j");
			ii2 = k_first + i2;
			args.m_ii(2, ii2);
			plesken_product(Ainf, coeff, 
				orbit_length, &args, K_first, K_len, 
				highest_layer, num_layers, 0 /* k_min */, 
				&type, FALSE /* f_v */, FALSE /* f_vv */);
			if (type.s_li() != num_layers)
				return error("calc_intersections_Ik3() type.s_li() != num_layers");
			ol = orbit_length->s_i(ii1);
			type.divide_out(ol);
			for (i = 0; i < num_layers; i++) {
				type.s_i(i)->copy(Ik3->s_ij(k_len + j, i));
				}
			j++;
			}
		if ((i1 + 1) % 10 == 0) {
			printf(",");
			if ((i1 + 1) % 50 == 0)
				printf(" %ld\n", i1 + 1);
			}
		else
			printf(".");
		fflush(stdout);
		}
	printf(";\n");
	fflush(stdout);
	
	j = 0;
	for (i1 = 0; i1 < k_len; i1++) {
		ii1 = k_first + i1;
		args.m_ii(0, ii1);

		for (i2 = i1 + 1; i2 < k_len; i2++) {
			ii2 = k_first + i2;
			args.m_ii(1, ii2);
			
			for (i3 = i2 + 1; i3 < k_len; i3++) {
#if 0
				if (i1 != j)
					return error("calc_intersections_Ik3() (i1 * k_len + i2) != j");
#endif
				ii3 = k_first + i3;
				args.m_ii(2, ii3);
				plesken_product(Ainf, coeff, 
					orbit_length, &args, K_first, K_len, 
					highest_layer, num_layers, 0 /* k_min */, 
					&type, FALSE /* f_v */, FALSE /* f_vv */);
				if (type.s_li() != num_layers)
					return error("calc_intersections_Ik3() type.s_li() != num_layers");
				ol = orbit_length->s_i(ii1);
				type.divide_out(ol);
				for (i = 0; i < num_layers; i++) {
					type.s_i(i)->copy(Ik3->s_ij(k_len + k_len2 + j, i));
					}
				j++;
				} // next i3
			} // next i2
		if ((i1 + 1) % 10 == 0) {
			printf(",");
			if ((i1 + 1) % 50 == 0)
				printf(" %ld\n", i1 + 1);
			}
		else
			printf(".");
		fflush(stdout);
		} // next i1
	printf(";\n");
	fflush(stdout);
	
	if (f_v) {
		printf("calc_intersections_Ik3() \n");
#if 0
		intersections_Ik3_blow_up(Ik2, k, &Ik2_vec, 
			orbit_length, K_first, K_len, 
			highest_layer, num_layers, TRUE /* f_v */);
#endif
		// Ik3->Print();
		fflush(stdout);
		}
	t1 = os_ticks();
	
	user_time = t1 - t0;
	strcpy(str, "Running time for calc_intersections_Ik3(): ");
	print_delta_time(user_time, str);
	printf("%s\n", str);

	return OK;
}




#endif /* LADDER_TRUE */


