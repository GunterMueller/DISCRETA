/* conway.C */

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

#if TEXDOCU
INT construct_Conway3_276(VECTOR_OP Co3_gen)
#endif
{
	INT t0, t1, user_time;
	BYTE s[256];
	VECTOR_OB W23_blocks, M23_gen, stab_gen;
	VECTOR_OP B0;
	MATRIX_OB G, H;
	INT i, j, n, a;
	INT nrow, ncol, nb_X, nb_X1, *theX = NIL, back_to, block_idx;
	LABRA_OB aut;
	SYM_OB stab_go, tmp1;
	PERMUTATION_OB p, q;
	VECTOR_OB aut_gen;
	LABRA_OB L;
	VECTOR_OB B;
	
	t0 = os_ticks();
	
	printf("construct_Conway3_276()\n");
	construct_W23(&W23_blocks, &M23_gen, &stab_gen, &block_idx, TRUE);
	B0 = (VECTOR_OP) W23_blocks.s_i(block_idx);
	printf("B0 = ");
	B0->println();
		{
		INT sf = f_perm_print_start_with_zero;
		f_perm_print_start_with_zero = TRUE;
		stab_gen.Print();
		f_perm_print_start_with_zero = sf;
		}
	
	construct_graph_from_design(23 /* v */, &W23_blocks, &G, 3);
	printf("calling construct_two_graph()\n");
	fflush(stdout);
	construct_two_graph(&G, &B);
	printf("the two-graph (l=%ld):\n", B.s_li());
	B.Print();
	fflush(stdout);

#if 1
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
	nrow = G.s_hi();
	ncol = G.s_li();
	
	for (i = 0; i < nrow; i++) {
		n = 0;
		for (j = 0; j < ncol; j++) {
			if (j == i)
				continue;
			a = G.s_iji(i, j);
			if (a)
				n++;
			}
		printf("row %ld has sum %ld\n", i, n);
		}
	
	nb_X = 0;
	for (i = 0; i < nrow; i++) {
		n = 0;
		for (j = i + 1; j < ncol; j++) {
			a = G.s_iji(i, j);
			if (a)
				n++;
			}
		// printf("row %ld has sum %ld\n", i, n);
		nb_X += n;
		}
	printf("nb_X = %ld\n", nb_X);
	fflush(stdout);
	theX = (INT *) my_malloc(2 * nb_X * sizeof(INT), "theX");
	
	n = 0;
	nb_X1 = 0;
	for (i = 0; i < nrow; i++) {
		for (j = i + 1; j < ncol; j++) {
			a = G.s_iji(i, j);
			if (a == 0)
				continue;
			theX[2 * nb_X1] = i;
			theX[2 * nb_X1 + 1] = j;
			nb_X1++;
			}	
		}
#ifdef SYM_GEO
	geo_canon_simple(FALSE, &back_to, nb_X, ncol, nb_X, TRUE /* f_print_dots */, 
		theX, &p, &q, FALSE /* f_transposed */, 
		TRUE /* f_get_aut_group */, &aut, 
		FALSE /* f_canon_v */,  FALSE /* f_canon_vv */);
#else
	return error("construct_Conway3_276() GEOLIB not available !");
#endif
	// changes theX !
	printf("\n");
	// printf("automorphism group order ");
	// ago.println();
	aut.save("conway3_276.dsc");
	// aut.generators(&aut_gen, &ago);
	aut.reduced_generating_set(&aut_gen, FALSE /* f_bottom_up */, TRUE /* f_v */);
	aut_gen.save("conway3_276.dsc");
	write_file_of_generators(&aut_gen, "conway3_276.txt");
	fflush(stdout);
	aut_gen.swap(Co3_gen);

	
	my_free(theX);

	t1 = os_ticks();
	user_time = t1 - t0;
	s[0] = 0;
	print_delta_time(user_time, s);
	printf("total computing time: %s\n", s);
	fflush(stdout);
	return OK;
}

#endif /* LADDER_TRUE */

