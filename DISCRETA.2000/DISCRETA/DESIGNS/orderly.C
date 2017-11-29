/* orderly.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#include <stdlib.h> // for system

#ifdef LADDER_TRUE

#include <DISCRETA/perm.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/lb.h>
#include <DISCRETA/cp.h>
#include <DISCRETA/divs.h>
#include <DISCRETA/ladder.h>
#include <DISCRETA/geo.h>


ORDERLY_SETS *orderly_sets_open(INT k, LABRA_OP G, INT f_stabilizer_action_verbose)
{
	ORDERLY_SETS *O = NIL;
	INT n;
	
	O = (ORDERLY_SETS *) my_malloc(sizeof(ORDERLY_SETS), "orderly_sets_open O");
	O->n = n = G->s_degree_i();
	O->k = k;
	O->l = 0;
	O->G = G;
	O->TG = (MATRIX_OP) callocobject("orderly_sets_open TG");
	O->go = (SYM_OP) callocobject("orderly_sets_open go");
	O->set = (INT *) my_malloc((k + 1) * sizeof(INT), "orderly_sets_open: set");
	O->nb_digits = (INT *) my_malloc((k + 1) * sizeof(INT), "orderly_sets_open: nb_digits");
	O->cur_digit = (INT *) my_malloc((k + 1) * sizeof(INT), "orderly_sets_open: cur_digit");
	O->digits = (INT *) my_malloc((k + 1) * n * sizeof(INT), "orderly_sets_open: digits");
	G->calc_transversal_matrix(O->TG, FALSE /* f_v */);
	G->group_order(O->go);
	O->f_stabilizer_action_verbose = f_stabilizer_action_verbose;
	return O;
}

INT orderly_sets_close(ORDERLY_SETS * O)
{
	freeall(O->TG);
	freeall(O->go);
	my_free(O->set);
	my_free(O->nb_digits);
	my_free(O->cur_digit);
	my_free(O->digits);
	my_free(O);
	return OK;
}

INT orderly_sets_is_canonical(ORDERLY_SETS *O, LABRA_OP aut, SYM_OP ago)
{
	INT r;
	
	r = geo_maxtest_set(O->n, O->l, O->set, 
		O->G, O->TG, 
		aut, ago, FALSE /* f_v */, FALSE /* f_vv */);
	// printf("orderly_sets_is_canonical = %ld\n", r);
	if (r == -1)
		return TRUE;
	return FALSE;
}

INT orderly_sets_first(ORDERLY_SETS *O, LABRA_OP stab, SYM_OP stab_order)
{
	INT is_c;
	
	orderly_sets_first_set(O);
	is_c = orderly_sets_is_canonical(O, stab, stab_order);
	if (!is_c)
		return error("orderly_sets_first: empty set is not canonical!");
	orderly_sets_compute_good_digits(O, stab);
	return OK;
}

INT orderly_sets_next(ORDERLY_SETS *O, LABRA_OP stab, SYM_OP stab_order)
{
	INT is_c, f_next;
	
	f_next = orderly_sets_next_set(O);
	if (!f_next)
		return FALSE;
	
	while (TRUE) {
		// print_map(X, l);
		
		is_c = orderly_sets_is_canonical(O, stab, stab_order);
		if (is_c) {
			// printf("is canonical\n");
			if (O->l < O->k) {
				orderly_sets_compute_good_digits(O, stab);
				}
			return TRUE;
			}
		// printf("not canonical\n");
		if (is_c)
			f_next = orderly_sets_next_set(O);
		else
			f_next = orderly_sets_next_set_no_descend(O);
		if (!f_next)
			break;
		}
	return FALSE;
}

INT orderly_sets_first_set(ORDERLY_SETS *O)
{
	O->l = 0;
	return TRUE;
}

INT orderly_sets_next_set(ORDERLY_SETS *O)
{
	if (O->l < O->k) {
		if (orderly_sets_first_digit(O, O->l)) {
			O->l++;
			return TRUE;
			}
		}
	while (O->l > 0) {
		if (orderly_sets_next_digit(O, O->l - 1)) {
			return TRUE;
			}
		O->l--;
		}
	return FALSE;
}

INT orderly_sets_next_set_no_descend(ORDERLY_SETS *O)
{
	while (O->l > 0) {
		if (orderly_sets_next_digit(O, O->l - 1)) {
			return TRUE;
			}
		O->l--;
		}
	return FALSE;
}

INT orderly_sets_first_digit(ORDERLY_SETS *O, INT i)
{
	INT nb, d;

	nb = O->nb_digits[i];
	if (nb == 0)
		return FALSE;
	O->cur_digit[i] = 0;
	d = O->digits[i * O->n + 0];
	O->set[i] = d;
#if 0
	INT f;
	
	if (i == 0) {
		f = 0;
		}
	else {
		f = O->set[i - 1] + 1;
		}
	if (f == O->n)
		return FALSE;
	O->set[i] = f;
#endif
	return TRUE;
}

INT orderly_sets_next_digit(ORDERLY_SETS *O, INT i)
{
	INT d;
	
	if (O->cur_digit[i] < O->nb_digits[i] - 1) {
		d = O->digits[i * O->n + ++O->cur_digit[i]];
		O->set[i] = d;
		return TRUE;
		}
	else
		return FALSE;
#if 0
	INT f;
	
	O->set[i]++;
	if (O->set[i] == O->n)
		return FALSE;
	return TRUE;
#endif
}

#define LET_STABILIZER_ACT
#define STABILIZER_ACTION_VERBOSE

INT orderly_sets_compute_good_digits(ORDERLY_SETS *O, LABRA_OP aut)
{
	VECTOR_OB gen;
	SYM_OB ago;
	INT i, l, f;
	INT *orbit = NIL;
	INT *orbit_min = NIL;
	INT *q = NIL, nb_q;
	INT ff, cur, next, r, nb_orbits;
	PERMUTATION_OP p;
	
	aut->generators(&gen);
	aut->group_order(&ago);
	r = gen.s_li();
	l = O->l;
	if (l == 0) {
		f = 0;
		}
	else {
		f = O->set[l - 1] + 1;
		}
#ifdef LET_STABILIZER_ACT
	if (ago.einsp()) {
		O->nb_digits[l] = 0;
		for (i = f; i < O->n; i++) {
			O->digits[l * O->n + O->nb_digits[l]++] = i;
			}
		if (O->f_stabilizer_action_verbose) {
			printf("found %ld good digits (trivial stabilizer action)\n", 
				O->nb_digits[l]);
			}
		return TRUE;
		}
	
	orbit = (INT *) my_malloc(O->n * sizeof(INT), "orderly_sets_compute_good_digits: orbit");
	orbit_min = (INT *) my_malloc(O->n * sizeof(INT), "orderly_sets_compute_good_digits: orbit_min");
	q = (INT *) my_malloc(O->n * sizeof(INT), "orderly_sets_compute_good_digits: q");
	for (i = 0; i < O->n; i++) {
		orbit[i] = -1;
		}
	nb_orbits = 0;
	if (O->f_stabilizer_action_verbose) {
		printf("computing stabilizer action on [%ld..%ld]\n", f, O->n - 1);
		}
	for (ff = f; ff < O->n; ff++) {
		if (orbit[ff] != -1)
			continue;
		nb_q = 1;
		q[0] = ff;
		orbit_min[nb_orbits] = ff;
		orbit[ff] = nb_orbits;
		
		while (nb_q) {
			cur = q[--nb_q];
			// printf("cur=%ld orbit_min=%ld\n", cur, orbit_min[nb_orbits]);
		
			for (i = 0; i < r; i++) {
				p = (PERMUTATION_OP) gen.s_i(i);
				next = p->s_ii(cur) - 1;
				if (orbit[next] == -1) {
					// printf("next=%ld\n", next);
					orbit[next] = nb_orbits;
					q[nb_q++] = next;
					if (next < orbit_min[nb_orbits]) {
						orbit_min[nb_orbits] = next;
						// printf("orbit_min set to %ld\n", next);
						}
					}
				} // next i
			} // while (nb_q)
		if (O->f_stabilizer_action_verbose) {
			printf("found orbit %ld with minimal element %ld\n", 
				nb_orbits, orbit_min[nb_orbits]);
			}
		nb_orbits++;
		} // next ff
	
	// now we extract the "good" digits:
	// the "good" digits are the minimal elements of orbits which are 
	// greater than or equal to f
	O->nb_digits[l] = 0;
	for (i = 0; i < nb_orbits; i++) {
		if (orbit_min[i] >= f) {
			O->digits[l * O->n + O->nb_digits[l]++] = orbit_min[i];
			}
		}
	if (O->f_stabilizer_action_verbose) {
		printf("found %ld good digits\n", O->nb_digits[l]);
		}
#else /* LET_STABILIZER_ACT */
	O->nb_digits[l] = 0;
	for (i = f; i < O->n; i++) {
		O->digits[l * O->n + O->nb_digits[l]++] = i;
		}
	printf("found %ld good digits (no stabilizer action)\n", O->nb_digits[l]);
#endif
	my_free(orbit);
	my_free(orbit_min);
	my_free(q);
	return OK;
}



#endif /* LADDER_TRUE */

