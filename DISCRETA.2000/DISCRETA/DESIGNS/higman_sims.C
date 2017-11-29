/* higman_sims.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1997
 */


#include <DISCRETA/discreta.h>

#ifdef LADDER_TRUE

#include <DISCRETA/perm.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/lb.h>
#include <DISCRETA/cp.h>
#include <DISCRETA/divs.h>
#include <DISCRETA/ladder.h>

#ifdef SYM_GEO
#include <DISCRETA/geo.h>
#endif

static INT get_blocks(VECTOR_OP blocks, VECTOR_OP block_selection, 
	VECTOR_OP selected_blocks);

#if TEXDOCU
INT construct_Higman_Sims_176(VECTOR_OP HS_gen)
#endif
{
	INT t0, t1, user_time;
	BYTE s[256];
	VECTOR_OB W24_blocks, M24_gen, stab_gen, stab1_gen, stab2_gen;
	VECTOR_OB V_on_points, V_on_blocks;
	VECTOR_OB H_point_idx;
	VECTOR_OB H_block_idx;
	MATRIX_OB G, H;
	INT i, j, i0, j0, k, n, a, l;
	INT nrow, ncol, nb_X, *theX = NIL, back_to, block_idx;
	INT point_a, point_b;
	VECTOR_OP B0;
	LABRA_OB aut;
	SYM_OB stab_go, tmp1;
	PERMUTATION_OB p, q;
	VECTOR_OB aut_gen;
	LABRA_OB L;
	VECTOR_OB point_blocks, block_blocks;
	
	t0 = os_ticks();
	
	printf("construct_Higman_Sims_176()\n");
	construct_W24(&W24_blocks, &M24_gen, &stab_gen, &block_idx, TRUE);
	B0 = (VECTOR_OP) W24_blocks.s_i(block_idx);
	printf("B0 = ");
	B0->println();
		{
		INT sf = f_perm_print_start_with_zero;
		f_perm_print_start_with_zero = TRUE;
		stab_gen.Print();
		f_perm_print_start_with_zero = sf;
		}
	point_a = B0->s_ii(0);
	point_b = B0->s_ii(1);
	printf("point a = %ld point b = %ld\n", point_a, point_b);
	printf("stabilizing a and b\n");
	fflush(stdout);
	vec_generators_stabilize_a_point(&stab_gen, point_a, &stab1_gen);
	vec_generators_stabilize_a_point(&stab1_gen, point_b, &stab2_gen);
	reduce_generators_labra(&stab2_gen, &stab_go, 
		FALSE /* f_verbose */, &L);
	printf("group order ");
	stab_go.println();
	fflush(stdout);
	L.reduced_generating_set(&stab_gen, FALSE /* f_bottom_up */, TRUE /* f_v */);
		{
		INT sf = f_perm_print_start_with_zero;
		f_perm_print_start_with_zero = TRUE;
		stab_gen.Print();
		f_perm_print_start_with_zero = sf;
		}
	fflush(stdout);
	
	construct_Higman_design(&W24_blocks, point_a, point_b, 
		&H_point_idx, &H_block_idx, &G);
	nrow = G.s_hi();
	ncol = G.s_li();
	
	get_blocks(&W24_blocks, &H_point_idx, &point_blocks);
	get_blocks(&W24_blocks, &H_block_idx, &block_blocks);
	// important: point_blocks / block_blocks is sorted !
	printf("get_blocks finished()!\n");
	stab_gen.copy(&V_on_points);
	vec_induce_action_on_blocks(&V_on_points, &point_blocks);
	printf("V on points:\n");
	V_on_points.Print();
	fflush(stdout);
	write_file_of_generators(&V_on_points, "V_on_points.txt");
	
	stab_gen.copy(&V_on_blocks);
	vec_induce_action_on_blocks(&V_on_blocks, &block_blocks);
	printf("V on blocks:\n");
	V_on_blocks.Print();
	fflush(stdout);
	write_file_of_generators(&V_on_blocks, "V_on_blocks.txt");

	{
		VECTOR_OB P_tda_ordering, P_orbit_first, P_orbit_length;
		VECTOR_OB Q_tda_ordering, Q_orbit_first, Q_orbit_length;
		LABRA_OB P_lab, Q_lab;
		PERMUTATION_OB p, pv, q, qv;
	
		printf("computing P_lab\n"); fflush(stdout);
		reduce_generators_labra(&V_on_points, &tmp1, 
			FALSE /* f_verbose */, &P_lab);
		printf("computing Q_lab\n"); fflush(stdout);
		reduce_generators_labra(&V_on_blocks, &tmp1, 
			FALSE /* f_verbose */, &Q_lab);
		{
			VECTOR_OB P_O_labra;
			P_lab.orbits(&P_tda_ordering, &P_orbit_first, &P_orbit_length, 
				FALSE /* f_calc_labras */, &P_O_labra, TRUE /* f_v */);
			printf("finished with P-orbits\n");
		}
		{
			VECTOR_OB Q_O_labra;
			Q_lab.orbits(&Q_tda_ordering, &Q_orbit_first, &Q_orbit_length, 
				FALSE /* f_calc_labras */, &Q_O_labra, TRUE /* f_v */);
			printf("finished with Q-orbits\n");
		}
		p.m_il(nrow);
		for (i = 0; i < nrow; i++) {
			j = P_tda_ordering.s_ii(i);
			p.m_ii(j, i + 1);
			}
		p.invers(&pv);
		q.m_il(ncol);
		for (i = 0; i < ncol; i++) {
			j = Q_tda_ordering.s_ii(i);
			q.m_ii(j, i + 1);
			}
		q.invers(&qv);
		
		for (i = 0; i < nrow; i++) {
			i0 = pv.s_ii(i) - 1;
			for (j = 0; j < ncol; j++) {
				j0 = qv.s_ii(j) - 1;
				a = G.s_iji(i0, j0);
				if (a)
					printf("X");
				else
					printf(".");
				}
			printf("\n");
			}
		printf("\n");
		printf("P_orbit_first = \n");
		P_orbit_first.println();
		printf("Q_orbit_first = \n");
		Q_orbit_first.println();
		fflush(stdout);

		if (P_orbit_first.s_li() != Q_orbit_first.s_li())
			return error("P_orbit_first.s_li() != Q_orbit_first.s_li()");
		for (k = 0; k < P_orbit_first.s_li(); k++) {
			i0 = P_orbit_first.s_ii(k);
			l = P_orbit_length.s_ii(k);
			H.m_ilih(l, l);
			for (i = 0; i < l; i++) {
				for (j = 0; j < l; j++) {
					a = G.s_iji(i0 + i, i0 + j);
					H.m_iji(i, j, a);
					}
				}
			printf("orbit %ld:\n", k);
			for (i = 0; i < l; i++) {
				for (j = 0; j < l; j++) {
					a = H.s_iji(i, j);
					if (a)
						printf("X");
					else
						printf(".");
					}
				printf("\n");
				}
			printf("\n");
			fflush(stdout);
			for (i = 0; i < l; i++) {
				n = 0;
				for (j = 0; j < l; j++) {
					a = H.s_iji(i, j);
					if (a)
						n++;
					}
				printf("sum row %ld: %ld\n", i, n);
				}
			
			}
		
	}

#if 0
	printf("\n");
	printf("%ld %ld\n", G.s_hi(), G.s_li());
	for (i = 0; i < G.s_hi(); i++) {
		for (j = 0; j < G.s_li(); j++) {
			printf("%ld ", G.s_iji(i, j));
			}
		printf("\n");
		}
	printf("\n");
#endif
	nb_X = 0;
	for (i = 0; i < nrow; i++) {
		n = 0;
		for (j = 0; j < ncol; j++) {
			a = G.s_iji(i, j);
			if (a)
				n++;
			}
		// printf("row %ld has sum %ld\n", i, n);
		nb_X += n;
		}
	theX = (INT *) my_malloc(nb_X * sizeof(INT), "theX");
	
	n = 0;
	for (i = 0; i < nrow; i++) {
		for (j = 0; j < ncol; j++) {
			a = G.s_iji(i, j);
			if (a == 0)
				continue;
			theX[n++] = i * ncol + j;
			}	
		}
#ifdef SYM_GEO
	geo_canon_simple(FALSE, &back_to, nrow, ncol, nb_X, TRUE /* f_print_dots */, 
		theX, &p, &q, FALSE /* f_transposed */, 
		TRUE /* f_get_aut_group */, &aut, 
		FALSE /* f_canon_v */,  FALSE /* f_canon_vv */);
#else
	return error("construct_Higman_Sims_176() GEOLIB not available !");
#endif
	// changes theX !
	printf("\n");
	// printf("automorphism group order ");
	// ago.println();
	// aut.save("aut_higman_design.dsc");
	// aut.generators(&aut_gen, &ago);
	aut.reduced_generating_set(&aut_gen, FALSE /* f_bottom_up */, TRUE /* f_v */);
	aut_gen.save("aut_higman_design_generators.dsc");
	write_file_of_generators(&aut_gen, "aut_higman_design_generators.txt");
	fflush(stdout);
	aut_gen.swap(HS_gen);

	
	my_free(theX);

	t1 = os_ticks();
	user_time = t1 - t0;
	s[0] = 0;
	print_delta_time(user_time, s);
	printf("total computing time: %s\n", s);
	fflush(stdout);
	return OK;
}

#if TEXDOCU
INT construct_Higman_design(VECTOR_OP W24_blocks, INT point_a, INT point_b, 
	VECTOR_OP H_point_idx, VECTOR_OP H_block_idx, MATRIX_OP I)
#endif
{
	VECTOR_OP B, B1;
	INT i, j, idx1, idx2, n, l, ll, llp = 0, llb = 0;
	INT idx, f_found_a, f_found_b;
	INTEGER_OB inta, intb;

	inta.m_i(point_a);
	intb.m_i(point_b);
	l = W24_blocks->s_li();
	H_point_idx->m_il(l);
	H_block_idx->m_il(l);
	for (i = 0; i < l; i++) {
		B = (VECTOR_OP) W24_blocks->s_i(i);
		if (B->s_li() != 8)
			return error("construct_Higman_design() B was wrong length");

		// The elements are numbered 0..23 !
		B->search(B->s_li(), TRUE, &inta, &idx, &f_found_a);
		B->search(B->s_li(), TRUE, &intb, &idx, &f_found_b);
	
		if (f_found_a && !f_found_b) {
			H_block_idx->m_ii(llb, i);
			llb++;
			}
		if (!f_found_a && f_found_b) {
			H_point_idx->m_ii(llp, i);
			llp++;
			}
		}
	printf("construct_Higman_design() found %ld points and %ld blocks\n", llp, llb); fflush(stdout);
	if (llp != 176)
		return error("construct_Higman_design() llp!= 176");
	if (llb != 176)
		return error("construct_Higman_design() llb!= 176");
	H_point_idx->realloc_z(llp);
	H_block_idx->realloc_z(llb);
	ll = llp;
	I->m_ilih_n(ll, ll);
	for (i = 0; i < ll; i++) {
		idx1 = H_point_idx->s_ii(i);
		B = (VECTOR_OP) W24_blocks->s_i(idx1);
		for (j = 0; j < ll; j++) {
			idx2 = H_block_idx->s_ii(j);
			B1 = (VECTOR_OP) W24_blocks->s_i(idx2);
			n = block_intersection(B, B1);
			if (n == 0 || n == 4)
				I->m_iji(i, j, 1);
			}
		}
	return OK;
}

#if TEXDOCU
INT construct_Higman_Sims_graph(MATRIX_OP G)
#endif
{
	VECTOR_OB W22_blocks;
	// INT s, i, j, k;
	
	construct_W22(&W22_blocks);
	printf("found a 3-(22,6,1) Steiner system with b=%ld blocks !\n"
		"this is W22 !\n", W22_blocks.s_li());
	construct_Higman_Sims_graph2(&W22_blocks, G);
#if 0
	printf("\n");
	printf("%ld %ld\n", G->s_hi(), G->s_li());
	for (i = 0; i < 100; i++) {
		for (j = 0; j < 100; j++) {
			printf("%ld ", G->s_iji(i, j));
			}
		printf("\n");
		}
	printf("\n");
#endif
	
#if 0
	G->Print();
	for (i = 0; i < 100; i++) {
		s = 0;
		for (j = 0; j < 100; j++) {
			k = get_graph_value(G, i, j);
			if (k)
				s++;
			}
		printf("row sum %ld: %ld\n", i, s);
		}
	for (j = 0; j < 100; j++) {
		s = 0;
		for (i = 0; i < 100; i++) {
			k = get_graph_value(G, i, j);
			if (k)
				s++;
			}
		printf("col sum %ld: %ld\n", j, s);
		}
#endif
	return OK;
}

#if TEXDOCU
INT construct_Higman_Sims_graph2(VECTOR_OP W22_blocks, MATRIX_OP G)
#endif
{
	VECTOR_OP B, B1, B2;
	INT i, j, k, kk, ii, a;
	
	if (W22_blocks->s_li() != 77)
		return error("construct_Higman_Sims_graph2(): W22_blocks->s_li() != 77");
	G->m_ilih_n(100, 100);

	// the additional point is joined with the points of W22:
	for (j = 1; j <= 22; j++)
		G->m_iji(0, j, 1);
	for (k = 0; k < 77; k++) {
		B = (VECTOR_OP) W22_blocks->s_i(k);
		if (B->s_li() != 6)
			return error("construct_Higman_Sims_graph2(): B->s_li() != 6");
		for (ii = 0; ii < 6; ii++) {
			a = B->s_ii(ii);
			// block k contains point a
			i = a + 1;
			j = 23 + k;
			G->m_iji(i, j, 1);
			}
		
		}
	for (k = 0; k < 77; k++) {
		for (kk = k + 1; kk < 77; kk++) {
			B1 = (VECTOR_OP) W22_blocks->s_i(k);
			B2 = (VECTOR_OP) W22_blocks->s_i(kk);
			if (blocks_disjoint(B1, B2)) {
				i = 23 + k;
				j = 23 + kk;
				G->m_iji(i, j, 1);
				}
			}
		}
	return OK;
}

#if TEXDOCU
static INT get_blocks(VECTOR_OP blocks, VECTOR_OP block_selection, 
	VECTOR_OP selected_blocks)
#endif
{
	// INT b = blocks->s_li();
	INT b1 = block_selection->s_li();
	INT i, idx;

	selected_blocks->m_il(b1);
	for (i = 0; i < b1; i++) {
		idx = block_selection->s_ii(i);
		blocks->s_i(idx)->copy(selected_blocks->s_i(i));
		}
	return OK;
}

#endif /* LADDER_TRUE */


