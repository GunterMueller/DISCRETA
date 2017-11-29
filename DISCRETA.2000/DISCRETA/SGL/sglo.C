/* sglo.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SGL_TRUE

#include <DISCRETA/ma.h>
#include <DISCRETA/perm.h>
#include <DISCRETA/dimino.h>
#include <DISCRETA/sgl.h>

#undef SGO_EXTRA_VERBOSE
#undef SGO_EXTRA_VERY_VERBOSE
#undef DEBUG_TRACE_ON
#undef DEBUG_ORBIT_ALGORITHM_SGO
#undef DEBUG_SPLIT_ORBIT_ALGORITHM_SGO

#undef CALC_TABLE

INT SGO_OB::field_name(INT i, INT j, BYTE *str)
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "f_has_aut"; break;
	case 1: s = "o_len"; break;
	case 2: s = "so_len"; break;
	case 3: s = "go"; break;
	case 4: s = "gen_Zidx"; break;
	case 5: s = "first"; break;
	case 6: s = "Zidx"; break;
	case 7: s = "SVlast"; break;
	case 8: s = "SVgen"; break;
	case 9: s = "Aut_SVlast"; break;
	case 10: s = "Aut_SVgen"; break;
	case 11: s = "p"; break;
	case 12: s = "pv"; break;
	case 13: s = "G_gen_on_orbit"; break;
	case 14: s = "Aut_gen_on_orbit"; break;
	case 15: s = "Nlayer"; break;
	case 16: s = "Norbit"; break;
	case 17: s = "Nrep"; break;
	default:
		return error("SGO::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

INT SGO_OB::Print(SGL_OP L, INT f_zuppos_expanded, INT type, void *data)
{
	VECTOR_OP Zidx;
	SYM_OP z;
	INT i, j, l, nb_groups, nb_zuppos;
	INT z_idx, f_has_Nreps;
	
	printf("ORBIT of (Zuppo generators):\n");
	L->print_zuppo_vector(s_gen_Zidx(), FALSE /* f_vertically */, type, data);
	printf("go = %ld f_has_aut = %ld "
		"o_len = %ld so_len = %ld first: %ld "
		"Nlayer = %ld Norbit = %ld\n", 
		s_go_i(), s_f_has_aut_i(), 
		s_o_len_i(), s_so_len_i(), s_first_i(), 
		s_Nlayer_i(), s_Norbit_i());
	f_has_Nreps = !s_Nrep()->emptyp();
	nb_groups = s_o_len_i();
	printf("table of (conjugacy class of) groups:\n");
	if (s_f_has_aut_i()) {
		printf("(Aut_SVlast, Aut_SVgen) ");
		}
	else {
		printf("(SVlast, SVgen) ");
		}
	if (f_has_Nreps)
		printf("nrep ");
	printf(" - list of zuppos contained in the group\n");
	printf("(the list of zuppos describes the subgroup uniquely)\n");
	for (i = 0; i < nb_groups; i++) {
		if (s_f_has_aut_i()) {
			printf("(%ld, %ld) ", s_Aut_SVlast_ii(i), s_Aut_SVgen_ii(i));
			}
		else {
			printf("(%ld, %ld) ", s_SVlast_ii(i), s_SVgen_ii(i));
			}
		fflush(stdout);
		if (f_has_Nreps) {
			printf("nrep %ld ", s_Nrep_ii(i));
			}
		printf(" - ");
		fflush(stdout);
		if (f_zuppos_expanded) {
			Zidx = s_Zidx_i(i);
			nb_zuppos = Zidx->s_li();
			for (j = 0; j < nb_zuppos; j++) {
				z_idx = Zidx->s_ii(j);
				if (f_zuppos_expanded) {
					z = L->s_Z_i(z_idx);
					do_print(z, type, data);
					fflush(stdout);
					if (j < nb_zuppos - 1)
						printf(", ");
					}
				else {
					printf("%ld ", z_idx);
					}
				}
			}
		printf("\n");
		fflush(stdout);
		}
	if (s_f_has_aut_i() || s_p()->s_obj_k() != EMPTY) {
		// printf("p = ");
		// s_p()->println();
		// printf("pv = ");
		// s_pv()->println();
		fflush(stdout);
		}
#if 0
	printf("generators for G on this orbit:\n");
	fflush(stdout);
	l = s_G_gen_on_orbit()->s_li();
	for (i = 0; i < l; i++)
		s_G_gen_on_orbit_i(i)->println();
	if (s_f_has_aut_i()) {
		printf("generators for Aut(G) on this orbit:\n");
		fflush(stdout);
		l = s_Aut_gen_on_orbit()->s_li();
		for (i = 0; i < l; i++)
			s_Aut_gen_on_orbit_i(i)->println();
		}
#endif
	printf("\n");

#ifdef CALC_TABLE
	{
		MATRIX_OB T;
		VECTOR_OB generators, embedding, embedding_inv;
		
		calc_group_table(L, s_first_i(), &T, &generators, 
			&embedding, &embedding_inv, type, data);
		/* T.Print();
		fflush(stdout); */
	}
#endif
	return OK;
}

INT SGO_OB::Init(
	VECTOR_OP H_gen, VECTOR_OP H, 
	VECTOR_OP N_gen, VECTOR_OP N, 
	SGL_OP L, INT f_verbose, 
	INT type, void *data)
{
	INT f_very_verbose = FALSE;
	INT o_len, f_has_aut;
	
	m_il(18);
	c_obj_k(SGO);
	
#ifdef SGO_EXTRA_VERBOSE
	f_verbose = TRUE;
#endif
#ifdef SGO_EXTRA_VERY_VERBOSE
	f_very_verbose = TRUE;
#endif
	f_has_aut = L->s_f_has_aut_group_i();
	s_f_has_aut()->m_i(f_has_aut);
	s_o_len()->m_i(0);
	s_so_len()->m_i(0);
	s_G_gen_on_orbit()->m_il(0);
	s_Aut_gen_on_orbit()->m_il(0);
	/* s_Nrep()->m_il(0); */
	s_go()->m_i(H->s_li());
	
	s_Nlayer()->m_i(-1);
	s_Norbit()->m_i(-1);
	if (L->Z2Zidx(H_gen, s_gen_Zidx(), type, data) != OK) {
		printf("SGO_OB::Init() error in L->Z2Zidx()\n");
		fflush(stdout);
		return ERROR;
		}
#if 0
	H_gen->copy(s_gen());
#endif
	if (f_very_verbose || f_verbose) {
		printf("SGO::Init() --- New orbit: H_gen = ");
		v_do_println(H_gen, 
			FALSE /* f_numerated */, 
			2 /* n_on_a_row */, 
			type, data);
		printf("group order = %ld f_has_aut = %ld\n", H->s_li(), f_has_aut);
		fflush(stdout);
		}
	H_gen->copy(N_gen);
	H->copy(N);
	if (f_has_aut) {
		calc_orbit_Aut_SV(L, H, f_verbose, f_very_verbose);
		o_len = s_o_len_i();
		s_Zidx()->v_shorten(o_len);
		s_Aut_SVlast()->v_shorten(o_len);
		s_Aut_SVgen()->v_shorten(o_len);
		
		gen_on_orbit(L, L->s_Aut_gen_on_zuppos(), 
			s_Aut_gen_on_orbit(), f_very_verbose);
		
		calc_splitting_orbits(L, N, N_gen, 
			f_verbose, f_very_verbose, type, data);
		}
	else {
		return error("without automorphisms not yet implemented");
		}
	
	if (f_verbose) {
		printf("found orbit of %ld conjugated groups\n", o_len);
		}
	if (f_very_verbose || f_verbose)
		Print(L, TRUE /* f_zuppos_expanded */, type, data);
	if (f_very_verbose || f_verbose) {
#if 0
		printf("N_gen = ");
		v_do_println(N_gen, 
			FALSE /* f_numerated */, 
			2 /* n_on_a_row */, 
			type, data);
		/* N_gen->println(); */
		printf("normalizer group order = %ld\n", N->s_li());
#endif
		}

	return OK;
}

INT SGO_OB::calc_splitting_orbits(SGL_OP L, 
	VECTOR_OP N, VECTOR_OP N_gen, 
	INT f_verbose, INT f_very_verbose, 
	INT type, void *data)
{
	INT *Q = NIL, nb_Q = 0, nb_gen, o_len;
	INT split, so_len = 0, os, k, cur, j, idx, f_found;
	VECTOR_OB H_zidx1, gpp;
	SYM_OB g1, g2, g3, gv;

	nb_gen = L->s_G_gen()->s_li();
	o_len = s_o_len_i();
	Q = (INT *) my_malloc(o_len * sizeof(INT), "sglo.C: calc_splitting_orbits");
	if (Q == NIL)
		return error("SGO::calc_splitting_orbits() no memory for Q");
	s_SVlast()->m_il(o_len);
	s_SVgen()->m_il(o_len);
	s_p()->m_il(o_len);
	s_pv()->m_il(o_len);
	for (j = 0; j < o_len; j++) {
		s_SVlast()->m_ii(j, -2);
		s_SVgen()->m_ii(j, -2);
		}
	for (split = 0; ;split++) {
		if (split == 0) {
			j = s_first_i();
			}
		else {
			for (j = 0; j < o_len; j++)
				if (s_SVlast_ii(j) == -2)
						break;
			if (j == o_len) {
				if (f_very_verbose) {
					printf("finished with orbit split -- "
						"o_len = %ld so_len = %ld split = %ld\n", 
						o_len, s_so_len_i(), split);
					fflush(stdout);
					}
				break;
				}
			}
		s_SVlast()->m_ii(j, -1);
		s_SVgen()->m_ii(j, -1);
		s_p()->m_ii(split * so_len + 0, j + 1);
		s_pv()->m_ii(j, split * so_len + 0 + 1);
		Q[0] = j;
		nb_Q = 1;
		os = 1;
		while (nb_Q) {
			/* apply all generators 
			 * to all groups of the orbit: */
			cur = Q[0];
			for (k = 1; k < nb_Q; k++)
				Q[k - 1] = Q[k];
			nb_Q--;
			trace_G_gen(L, cur, &g1, type, data);
			for (j = 0; j < nb_gen; j++) {
#ifdef DEBUG_SPLIT_ORBIT_ALGORITHM_SGO
				printf("split_orbit: cur = %ld j = %ld\n", cur, j);
#endif
				L->conjugate_zuppos_by_perm(s_Zidx_i(cur), 
					L->s_G_gen_on_zuppos_i(j), &H_zidx1);
				if (s_Zidx()->search(o_len, TRUE, 
					&H_zidx1, &idx, &f_found) != OK)
					return error("SGO::calc_splitting_orbits(): "
						"error in s_Zidx->search(H_zidx1)");
				if (!f_found)
					return error("SGO::calc_splitting_orbits(): not found");
				idx--;
				if (s_SVlast_ii(idx) == -2) {
					/* new orbit element */
#ifdef DEBUG_SPLIT_ORBIT_ALGORITHM_SGO
					printf("split_orbit: found new orbit element %ld\n", idx);
#endif
					s_SVlast()->m_ii(idx, cur);
					s_SVgen()->m_ii(idx, j);
					s_p()->m_ii(split * so_len + os, idx + 1);
					s_pv()->m_ii(idx, split * so_len + os + 1);
					os++;
					Q[nb_Q] = idx;
					nb_Q++;
					}
				else {
					if (split == 0) {
						/* Schreier generators for the 
						 * stabilizer of the first group only: */
						do_mult(&g1, L->s_G_gen_i(j), &g2, type, data);
						trace_G_gen(L, idx, &g3, type, data);
						do_invers(&g3, &gv, type, data);
						do_mult(&g2, &gv, &g3, type, data);
						gpp.m_il(0);
						prime_power_parts(&g3, &gpp, type, data);
						for (k = 0; k < gpp.s_li(); k++)
							dimino_extend(N, N_gen, gpp.s_i(k), NIL, 0, NIL, type, data);
						}
					}
				}
			}
		if (split == 0) {
			so_len = os;
			s_so_len()->m_i(os);
			}
		else {
			if (os != so_len)
				return error("os != so_len");
			}
#ifdef DEBUG_SPLIT_ORBIT_ALGORITHM_SGO
		printf("split_orbit: split = %ld so_len = %ld\n", split, so_len);
		printf("SVlast = ");
		s_SVlast()->println();
		printf("SVgen = ");
		s_SVgen()->println();
		fflush(stdout);
#endif
		} /* next split */

#ifdef DEBUG_SPLIT_ORBIT_ALGORITHM_SGO
	// printf("p = ");
	// s_p()->println();
	// printf("pv = ");
	// s_pv()->println();
	fflush(stdout);
#endif
	
	if (Q) {
		my_free(Q);
		Q = NIL;
		}
	return OK;
}

INT SGO_OB::calc_orbit_Aut_SV(SGL_OP L, VECTOR_OP H, 
	INT f_verbose, INT f_very_verbose)
{
	INT *Q = NIL, *Q_new, nb_Q = 0, nb_gen, ii, os, k, cur, j, idx, f_found, v_len;
	VECTOR_OB H_zidx, H_zidx1;
	PERMUTATION_OB g1;

	v_len = VECTOR_OVERSIZE;
	Q = (INT *) my_malloc(v_len * sizeof(INT), "sglo.C: calc_orbit_Aut_SV");
	if (Q == NIL)
		return error("SGO::calc_orbit_Aut_SV() no memory for Q");

	s_Zidx()->m_il(v_len);
	s_Aut_SVlast()->m_il(v_len);
	s_Aut_SVgen()->m_il(v_len);
	nb_gen = L->s_Aut_gen()->s_li();
	L->calc_zidx(H, &H_zidx);	
	/* H_zidx.println(); */
	H_zidx.copy(s_Zidx_i(0));
	s_Aut_SVlast()->m_ii(0, -1);
	s_Aut_SVgen()->m_ii(0, -1);
	s_first()->m_i(0);
	Q[0] = 0;
	nb_Q = 1;
	os = 1;
	while (nb_Q) {
		/* apply all generators 
		 * to all groups of the orbit: */
		cur = Q[0];
		for (k = 1; k < nb_Q; k++)
			Q[k - 1] = Q[k];
		nb_Q--;
		trace_Aut_gen_on_zuppos(L, cur, &g1);
		for (j = 0; j < nb_gen; j++) {
#ifdef DEBUG_ORBIT_ALGORITHM_SGO
			printf("cur = %ld j = %ld\n", cur, j);
			fflush(stdout);
#endif
			L->conjugate_zuppos_by_perm(s_Zidx_i(cur), 
				L->s_Aut_gen_on_zuppos_i(j), &H_zidx1);
			if (s_Zidx()->search(os, TRUE, &H_zidx1, &idx, &f_found) != OK)
				return error("SGO::Init(): error in s_Zidx->search(H_zidx1)");
			if (!f_found) {
				if (f_very_verbose) {
					printf("found new group: ");
					H_zidx1.println();
					fflush(stdout);
					}
				if (os >= v_len) {
					v_len += VECTOR_OVERSIZE;
					s_Zidx()->realloc_z(v_len);
					s_Aut_SVlast()->realloc_z(v_len);
					s_Aut_SVgen()->realloc_z(v_len);
					Q_new = (INT *) my_malloc(v_len * sizeof(INT), "sglo.C: Init");
					if (Q_new == NIL)
						return error("SGO::calc_orbit_Aut_SV() no memory for Q_new");
					for (ii = 0; ii < nb_Q; ii++)
						Q_new[ii] = Q[ii];
					my_free(Q);
					Q = Q_new;
					Q_new = NIL;
#ifdef DEBUG_ORBIT_ALGORITHM_SGO
					printf("realloc: new length = %ld\n", v_len);
					fflush(stdout);
#endif

					}
				/* neue Untergruppe an der Stelle 
				 * idx einfuegen, 
				 * first, SVlast und Q[] und 
				 * cur updaten: */
				if (s_first_i() >= idx)
					s_first()->inc();
				for (k = 0; k < os; k++) {
					if (s_Aut_SVlast_ii(k) >= idx)
						s_Aut_SVlast_i(k)->inc();
					}
				for (k = 0; k < nb_Q; k++)
					if (Q[k] >= idx)
						Q[k]++;
				if (cur >= idx)
					cur++;

				for (k = os; k > idx; k--) {
					s_Zidx_i(k)->swap(s_Zidx_i(k - 1));
					s_Aut_SVlast_i(k)->swap(s_Aut_SVlast_i(k - 1));
					s_Aut_SVgen_i(k)->swap(s_Aut_SVgen_i(k - 1));
					}
				H_zidx1.copy(s_Zidx_i(idx));
				s_Aut_SVlast()->m_ii(idx, cur);
					/* cur bereits aktualisiert */
				s_Aut_SVgen()->m_ii(idx, j);
				os++;
				/* printf("os = %ld\n", os); */
				Q[nb_Q++] = idx;
				/* keine Schreier Erzeuger 
				 * in diesem Falle */
				}
#if 0
	/* we do not want schreier generators here, 
	 * they would generate the stabilizer in the holomorph */
			else {
				idx--;
				do_mult(&g1, L->s_G_gen_i(j), &g2, 
					type, data);
				trace(L, idx, &g3, type, data);
				do_invers(&g3, &gv, type, data);
				do_mult(&g2, &gv, &g3, type, data);
				gpp.m_il(0);
				prime_power_parts(&g3, &gpp, type, data);
				for (k = 0; k < gpp.s_li(); k++)
					dimino_extend(N, N_gen, 
					gpp.s_i(k), NIL, 0, NIL, type, data);
				}
#endif
			}
		}
	if (Q) {
		my_free(Q);
		Q = NIL;
		}
	s_o_len()->m_i(os);
	return OK;
}

INT SGO_OB::calc_nreps(SGO_OP N_orbit, INT nrep0, INT f_verbose)
{
	PERMUTATION_OB g;
	INT i, len, nrep;
	
	len = s_Zidx()->s_li();
	s_Nrep()->m_il(len);
	for (i = 0; i < len; i++) {
		if (s_f_has_aut_i()) {
			trace_on(TRUE, i, &g, N_orbit->s_Aut_gen_on_orbit(), DO_TYPE_PERM, NIL);
			}
		else {
			trace_on(FALSE, i, &g, N_orbit->s_G_gen_on_orbit(), DO_TYPE_PERM, NIL);
			}
		nrep = g.s_ii(nrep0) - 1;
		s_Nrep()->m_ii(i, nrep);
		}
	return OK;
}

INT SGO_OB::recalc_Zidx(SGL_OP L, INT type, void *data, INT f_v)
{
	VECTOR_OP V;
	INT o_len, i;

	if (s_Zidx()->s_obj_k() != EMPTY)
		return OK;
	o_len = s_o_len_i();
	if (f_v) {
		printf("recalc_Zidx() recalculating Zidx for orbit of length %ld\n", o_len);
		fflush(stdout);
		}
	s_Zidx()->m_il(o_len);
	for (i = 0; i < o_len; i++) {
		if (f_v) {
			printf("%ld ", i);
			fflush(stdout);
			}
		V = s_Zidx_i(i);
		calc_Zidx(L, i, V, type, data);
		}
	if (f_v) {
		printf("\n");
		fflush(stdout);
		}
	return OK;
}

#define DEBUG_CALC_GROUP_TABLE
#undef DEBUG_CALC_GROUP_TABLE_VERBOSE

INT SGO_OB::calc_group_table(SGL_OP L, INT rep, 
	MATRIX_OP T, VECTOR_OP generators, 
	VECTOR_OP embedding, VECTOR_OP embedding_inv, 
	INT type, void *data)
{
	PERMUTATION_OB p, pv;
	VECTOR_OB ge;
	SYM_OB id, a;
	INT i, j, k, ii, jj, kk, l, idx, f_found, ll, go;
	
#ifdef DEBUG_CALC_GROUP_TABLE
	printf("SGO_OB::calc_group_table():");
	printf("calling calc_group_elements() rep = %ld\n", rep);
	fflush(stdout);
#endif
	calc_group_elements(L, rep, &ge, generators, type, data);
	l = ge.s_li();
	go = L->s_go_i();
#ifdef DEBUG_CALC_GROUP_TABLE
	printf("finished ! (subgroup order / "
		"group order = %ld / %ld)\n", l, go);
	fflush(stdout);
#endif
	L->s_G_gen_i(0)->copy(&id);
	do_one(&id, type, data);
	ge.search(l, TRUE /* f_ascending */, &id, &idx, &f_found);
	if (!f_found)
		return error("SGO_OB::calc_group_table() id not found");
	idx--;
#ifdef DEBUG_CALC_GROUP_TABLE
	printf("id element found at %ld (in a group of order %ld)\n", idx, l);
	fflush(stdout);
#endif
	p.m_il(l);
	pv.m_il(l);
	p.one();
	if (idx != 0) {
		p.m_ii(0, idx + 1);
		p.m_ii(idx, 0 + 1);
		ll = generators->s_li();
		for (i = 0; i < ll; i++) {
			k = generators->s_ii(i);
			k = p.s_ii(k) - 1;
			generators->m_ii(i, k);
			}
		printf("id element permuted at front position\n");
		}
	p.invers(&pv);
#ifdef DEBUG_CALC_GROUP_TABLE
	// printf("p = ");
	// p.println();
	// printf("pv = ");
	// pv.println();
	fflush(stdout);
#endif
	embedding->m_il_n(l);
	embedding_inv->m_il(go);
	for (i = 0; i < go; i++) {
		embedding_inv->m_ii(i, -1);
		}
	T->m_ilih(l, l);
	for (i = 0; i < l; i++) {
		ii = pv.s_ii(i) - 1;
		
		if (ge.s_i(ii)->s_obj_k() == INTEGER) {
			j = ge.s_ii(ii);
		
			embedding->m_ii(i, j);
			embedding_inv->m_ii(j, i);
			}
		}
	for (i = 0; i < l; i++) {
		ii = pv.s_ii(i) - 1;
		for (j = 0; j < l; j++) {
			jj = pv.s_ii(j) - 1;
			do_mult(ge.s_i(ii), ge.s_i(jj), &a, type, data);
			ge.search(l, TRUE /* f_ascending */, &a, &idx, &f_found);
			if (!f_found)
				return error("SGO_OB::calc_group_table() ii * jj not found");
			idx--;
			kk = idx;
			k = p.s_ii(kk) - 1;
#ifdef DEBUG_CALC_GROUP_TABLE_VERBOSE
			printf("%ld * %ld = %ld\n", i, j, k);
			fflush(stdout);
#endif
			T->m_iji(i, j, k);
			}
		}
#ifdef DEBUG_CALC_GROUP_TABLE_VERBOSE
	T->Print();
	fflush(stdout);
#endif
	return OK;
}

#undef DEBUG_CALC_GROUP_ELEMENTS

INT SGO_OB::calc_group_elements(SGL_OP L, INT rep, 
	VECTOR_OP ge, VECTOR_OP generators, INT type, void *data)
{
	VECTOR_OB gen_Zidx;

#ifdef DEBUG_CALC_GROUP_ELEMENTS
	printf("SGO_OB::calc_group_elements() calling calc_gen_Zidx()\n");
	fflush(stdout);
#endif
	calc_gen_Zidx(L, rep, &gen_Zidx);
#ifdef DEBUG_CALC_GROUP_ELEMENTS
	printf("SGO_OB::calc_group_elements() calling L->Zidx2Z()\n");
	fflush(stdout);
#endif
	L->Zidx2Z(&gen_Zidx, ge);
	ge->copy(generators);
	if (ge->s_li() == 0) {
		ge->inc();
		L->s_G_gen_i(0)->copy(ge->s_i(0));
		do_one(ge->s_i(0), type, data);
		}
#ifdef DEBUG_CALC_GROUP_ELEMENTS
	printf("SGO_OB::calc_group_elements() calling gruppen_elemente1()\n");
	fflush(stdout);
#endif
	gruppen_elemente1(ge, type, data);
#ifdef DEBUG_CALC_GROUP_ELEMENTS
	printf("SGO_OB::calc_group_elements() finished with gruppen_elemente1()\n");
	fflush(stdout);
#endif
	return OK;
}

INT SGO_OB::calc_Zidx(SGL_OP L, INT rep, 
	VECTOR_OP Zidx, INT type, void *data)
/* computes Zidx for the group 'rep' of the orbit 
 * WITHOUT using SGO's Zidx. */
{
	VECTOR_OB gen_Zidx;

	calc_gen_Zidx(L, rep, &gen_Zidx);
	L->Zidx2Z(&gen_Zidx, Zidx);
	gruppen_elemente1(Zidx, type, data);
	return OK;
}

INT SGO_OB::calc_gen_Zidx(SGL_OP L, INT rep, VECTOR_OP gen_Zidx)
/* computes generators for group number rep 
 * in this orbit.
 * the generators (computed into gen_Zidx) 
 * are images of the generators of the 'first' group
 * under conjugation by a fusing element. */
{
	PERMUTATION_OB g;
	INT i, j, k, len, idx, f_found;
	
	if (s_f_has_aut_i())
		trace_Aut_gen_on_zuppos(L, rep, &g);
	else
		trace_G_gen_on_zuppos(L, rep, &g);
	len = s_gen_Zidx()->s_li();
	gen_Zidx->m_il(len);
	for (i = 0; i < len; i++) {
		j = g.s_ii(s_gen_Zidx_ii(i)) - 1;
		if (gen_Zidx->search_and_insert_int(i, j) != OK)
			return error("SGO::calc_gen_Zidx(): duplicate generator");
		}
	return OK;
}

INT SGO_OB::gen_on_orbit(SGL_OP L, 
	VECTOR_OP gen_on_zuppos, VECTOR_OP gen_on_orbit, INT f_verbose)
{
	PERMUTATION_OP g_j, p;
	VECTOR_OB H_zidx;
	INT i, j, nb_gen, nb_groups, idx, f_found;
	
	nb_gen = gen_on_zuppos->s_li();
	nb_groups = s_Zidx()->s_li();
	gen_on_orbit->m_il(nb_gen);
	for (j = 0; j < nb_gen; j++) {
		g_j = (PERMUTATION_OP) gen_on_zuppos->s_i(j);
		p = (PERMUTATION_OP) gen_on_orbit->s_i(j);
		p->m_il(nb_groups);
		for (i = 0; i < nb_groups; i++) {
			L->conjugate_zuppos_by_perm(s_Zidx_i(i), g_j, &H_zidx);
			if (s_Zidx()->search(nb_groups, TRUE, &H_zidx, &idx, &f_found) != OK)
				return error("SGO::gen_on_orbit(): error in "
					"s_Zidx->search(H_zidx)");
			if (!f_found)
				return error("SGO::gen_on_orbit(): !f_found");
			idx--;
			p->m_ii(i, idx + 1);
			}
		if (f_verbose) {
			printf("SGO::gen_on_orbit(): generator %ld: \n", j);
			g_j->println();
			printf("on this orbit by conjugation: ");
			p->println();
			}
		}
	return OK;
}

INT SGO_OB::find_group(VECTOR_OP Zidx, INT *idx, INT *f_found)
{
	if (s_Zidx()->s_li() <= 0)
		return error("SGO::find_group(): s_Zidx()->s_li() <= 0");

	/* the number of zuppos must be equal: */
	if (s_Zidx_i(0)->s_li() != Zidx->s_li()) {
		*f_found = FALSE;
		return OK;
		}
	if (s_Zidx()->search(s_Zidx()->s_li(), TRUE, Zidx, idx, f_found) != OK)
		return error("SGO::find_group(): error in s_Zidx->search(Zidx)");
	if (*f_found)
		(*idx)--;
	return OK;
}

INT SGO_OB::trace_on(INT f_Aut_SV, INT i, SYM_OP g, VECTOR_OP Gen, INT type, void *data)
{
	SYM_OB g1;
	INT last, gen;
	
	if (Gen->s_li() <= 0)
		return error("SGO::trace_on(): no generators");
	Gen->s_i(0)->copy(g);
	do_one(g, type, data);
	while (TRUE) {
		if (f_Aut_SV) {
			last = s_Aut_SVlast_ii(i);
			gen = s_Aut_SVgen_ii(i);
			}
		else {
			last = s_SVlast_ii(i);
			gen = s_SVgen_ii(i);
			}
#ifdef DEBUG_TRACE_ON
		printf("trace_on: i = %ld last = %ld gen = %ld\n", i, last, gen);
		fflush(stdout);
#endif
		if (last == -1)
			break;
		do_mult(Gen->s_i(gen), g, &g1, type, data);
		g1.swap(g);
		i = last;
		}
	return OK;
}

INT SGO_OB::trace_G_gen(SGL_OP L, INT i, SYM_OP g, INT type, void *data)
{
	return trace_on(FALSE /* f_Aut_SV */, i, g, L->s_G_gen(), type, data);
}

INT SGO_OB::trace_G_gen_on_zuppos(SGL_OP L, INT i, PERMUTATION_OP g)
{
	return trace_on(FALSE /* f_Aut_SV */, i, g, L->s_G_gen_on_zuppos(), DO_TYPE_PERM, NIL);
}

INT SGO_OB::trace_Aut_gen_on_zuppos(SGL_OP L, INT i, PERMUTATION_OP g)
{
	return trace_on(TRUE /* f_Aut_SV */, i, g, L->s_Aut_gen_on_zuppos(), DO_TYPE_PERM, NIL);
}

INT SGO_OB::generators(SGL_OP L, VECTOR_OP gen, INT type, void *data)
{
	L->Zidx2Z(s_gen_Zidx(), gen);
	return OK;
}

#endif /* SGL_TRUE */

