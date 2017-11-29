/* sgls.C */

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

INT SGL_OB::field_name(INT i, INT j, BYTE *str)
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "Z"; break;
	case 1: s = "Zidx"; break;
	case 2: s = "Zorder"; break;
	case 3: s = "Zprime"; break;
	case 4: s = "Zk"; break;
	case 5: s = "PP"; break;
	case 6: s = "PPidx"; break;
	case 7: s = "zuppos_on_zuppos"; break;
	case 8: s = "nb_Layers"; break;
	case 9: s = "theOrbits"; break;
	case 10: s = "nb_subgroups"; break;
	case 11: s = "total_nb_orbits"; break;
	case 12: s = "total_nb_subgroups"; break;
	case 13: s = "G_gen"; break;
	case 14: s = "G_gen_on_zuppos"; break;
	case 15: s = "f_has_aut_group"; break;
	case 16: s = "Aut_gen"; break;
	case 17: s = "Aut_gen_on_zuppos"; break;
	case 18: s = "go"; break;
	case 19: s = "f_verbose"; break;
	case 20: s = "f_very_verbose"; break;
	default:
		return error("SGL::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

INT SGL_OB::Init_with_Aut_generators(VECTOR_OP G_gen, VECTOR_OP G, VECTOR_OP Aut_gen, 
	INT f_verbose, INT f_very_verbose, INT type, void *data)
{
	VECTOR_OB G_gen1, Zidx;

	if (Init1(G, f_verbose, f_very_verbose, type, data) != OK)
		return error("SGL::Init_with_Aut_generators(): error in Init1()");

	Z2Zidx_repetitions_allowed(G_gen, &Zidx, type, data);
	Zidx2Z(&Zidx, &G_gen1);

	Init2(&G_gen1, type, data);

	Aut_gen->copy(s_Aut_gen());
	printf("SGL_OB::Init_with_Aut_generators() calling calc_aut_on_zuppos()\n");
	fflush(stdout);
	if (calc_aut_on_zuppos(G) != OK)
		return error("SGL::Init_with_Aut_generators(): "
			"error in calc_aut_on_zuppos()");
	s_f_has_aut_group()->m_i(TRUE);
	add_trivial_subgroups(&G_gen1, G, f_verbose, type, data);
	return OK;
}

INT SGL_OB::Init(VECTOR_OP G_gen, VECTOR_OP G, 
	INT f_verbose, INT f_very_verbose, INT type, void *data)
{
	VECTOR_OB G_gen1, Zidx;

	Init1(G, f_verbose, f_very_verbose, type, data);

	Z2Zidx_repetitions_allowed(G_gen, &Zidx, type, data);
	Zidx2Z(&Zidx, &G_gen1);
	
	Init2(&G_gen1, type, data);

	add_trivial_subgroups(&G_gen1, G, f_verbose, type, data);
	/* PrintOrbit(0, 0);
	PrintOrbit(nb_primes_G, 0); */

	return OK;
}

INT SGL_OB::Init1(VECTOR_OP G, INT f_verbose, INT f_very_verbose, 
	INT type, void *data)
{
	INT nb_primes_G, k, go;
	INT erg = OK;
	
	m_il(21);
	c_obj_k(SGL);
	s_f_has_aut_group()->m_i(FALSE);
	s_Aut_gen()->m_il(0);
	s_Aut_gen_on_zuppos()->m_il(0);
	s_f_verbose()->m_i(f_verbose);
	s_f_very_verbose()->m_i(f_very_verbose);
		
	go = G->s_li();
	s_go()->m_i(go);
	nb_primes_G = nb_primes(go);
	s_total_nb_orbits()->m_i(0);
	s_total_nb_subgroups()->m_i(0);
	s_nb_Layers()->m_i(nb_primes_G);
	s_theOrbits()->m_il(nb_primes_G + 1);
	s_nb_subgroups()->m_il_n(nb_primes_G + 1);
	for (k = 0; k <= nb_primes_G; k++)
		s_theOrbits_i(k)->m_il(0);

	if (zuppos(G, s_f_very_verbose_i(), type, data) != OK)
		return error("SGL::Init1(): error in zuppos()");
	
	return OK;
}

INT SGL_OB::Init2(VECTOR_OP G_gen, INT type, void *data)
{
	if (G_gen->s_li() <= 0)
		return error("SGL::Init2(): no generators");
	
	G_gen->copy(s_G_gen());
	calc_G_gen_on_zuppos(type, data);
	
	/* PrintOrbit(0, 0);
	PrintOrbit(nb_primes_G, 0); */

	return OK;
}

INT SGL_OB::add_trivial_subgroups(VECTOR_OP G_gen, VECTOR_OP G, 
	INT f_verbose, INT type, void *data)
{
	VECTOR_OB E_gen, E;
	
	/* the one - group: */
	E_gen.m_il(0);
	E.m_il(1);
	G_gen->s_i(0)->copy(E.s_i(0));
	do_one(E.s_i(0), type, data);
	
	printf("SGL_OB::add_trivial_subgroups() adding G\n");
	fflush(stdout);
	if (add_group(G_gen, G, f_verbose, type, data) != OK)
		return error("SGL_OB::add_trivial_subgroups error in add_group(G_gen, G)");
	printf("SGL_OB::add_trivial_subgroups() finished with adding G\n");
	fflush(stdout);

	printf("SGL_OB::add_trivial_subgroups() adding 1 (trivial subgroup)\n");
	fflush(stdout);
	if (add_group(&E_gen, &E, f_verbose, type, data) != OK)
		return error("SGL_OB::add_trivial_subgroups error in add_group(E_gen, E)");
	printf("SGL_OB::add_trivial_subgroups() finished with adding 1\n");
	fflush(stdout);
	
	return OK;
}

#if 0
static INT ex2(VECTOR_OP V, VECTOR_OP G);

INT SGL_OB::exchange_enum_elements(VECTOR_OP G)
{
	ex2(s_Z(), G);
	ex2(s_PP(), G);
	ex2(s_G_gen(), G);
	return OK;
}

static INT ex2(VECTOR_OP V, VECTOR_OP G)
{
	INTEGER_OP int_op;
	INT i, j, len;

	len = V->s_li();
	for (i = 0; i < len; i++) {
		int_op = (INTEGER_OP) V->s_i(i);
		if (int_op->s_obj_k() != INTEGER)
			return error("SGL::ex2(): not an INTEGER object");
		j = int_op->s_i();
		G->s_i(j)->copy(V->s_i(i));
		}
	return OK;
}
#endif

INT SGL_OB::PrintOrbits(INT f_zuppos_expanded, INT type, void *data)
{
	INT nb_layers, nb_orbits, i, j;
	
	nb_layers = s_theOrbits()->s_li();
	for (i = 0; i < nb_layers; i++) {
		nb_orbits = s_theOrbits_i(i)->s_li();
		for (j = 0; j < nb_orbits; j++) {
			PrintOrbit(i, j, f_zuppos_expanded, type, data);
			}
		}
	return OK;
}

INT SGL_OB::PrintOrbit(INT layer, INT orbit, 
	INT f_zuppos_expanded, INT type, void *data)
{
	SGO_OP O;
	
	printf("layer %ld orbit %ld:\n", layer, orbit);
	O = s_theOrbits_ij(layer, orbit);
	O->Print(this, f_zuppos_expanded, type, data);
	return OK;
}

INT SGL_OB::PrintStatistics()
{
	INT i, len;
	
	len = s_nb_Layers_i();
	printf("{");
	for (i = 0; i <= len; i++) {
		printf("%ld ", s_theOrbits_i(i)->s_li());
		}
	printf("} total_nb_orbits = %ld\n", 
		s_total_nb_orbits_i());
	printf("{");
	for (i = 0; i <= len; i++) {
		printf("%ld ", s_nb_subgroups_ii(i));
		}
	printf("} total_nb_subgroups = %ld\n", 
		s_total_nb_subgroups_i());
	return OK;
}

INT SGL_OB::sprint_statistics(BYTE *str1, BYTE *str2, BYTE *str3)
{
	INT i, j, l, len, o_len, so_len, split;
	INT total_split_orbits, split_orbits;
	SGO_OP sgo;
	VECTOR_OP O;
	
	len = s_nb_Layers_i();
	sprintf(str1 + strlen(str1), "{ ");
	for (i = 0; i <= len; i++) {
		sprintf(str1 + strlen(str1), "%ld ", s_theOrbits_i(i)->s_li());
		}
	sprintf(str1 + strlen(str1), "} %ld", s_total_nb_orbits_i());
	
	total_split_orbits = 0;
	sprintf(str2 + strlen(str2), "{ ");
	for (i = 0; i <= len; i++) {
		O = s_theOrbits_i(i);
		l = O->s_li();
		split_orbits = 0;
		for (j = 0; j < l; j++) {
			sgo = s_theOrbits_ij(i, j);
			o_len = sgo->s_o_len_i();
			so_len = sgo->s_so_len_i();
			split = o_len / so_len;
			split_orbits += split;
			total_split_orbits += split;
			}
		sprintf(str2 + strlen(str2), "%ld ", split_orbits);
		}
	sprintf(str2 + strlen(str2), "} %ld", total_split_orbits);

	sprintf(str3 + strlen(str3), "{ ");
	for (i = 0; i <= len; i++) {
		sprintf(str3 + strlen(str3), "%ld ", s_nb_subgroups_ii(i));
		}
	sprintf(str3 + strlen(str3), "} %ld", s_total_nb_subgroups_i());
	return OK;
}

INT SGL_OB::sprint_statistics_tex(BYTE *str1, BYTE *str2, BYTE *str3)
{
	INT i, j, l, len, o_len, so_len, split;
	INT total_split_orbits, split_orbits;
	SGO_OP sgo;
	VECTOR_OP O;
	
	len = s_nb_Layers_i();
	sprintf(str1 + strlen(str1), "\\{ ");
	for (i = 0; i <= len; i++) {
		sprintf(str1 + strlen(str1), "%ld ", s_theOrbits_i(i)->s_li());
		}
	sprintf(str1 + strlen(str1), "\\} %ld", s_total_nb_orbits_i());
	
	total_split_orbits = 0;
	sprintf(str2 + strlen(str2), "\\{ ");
	for (i = 0; i <= len; i++) {
		O = s_theOrbits_i(i);
		l = O->s_li();
		split_orbits = 0;
		for (j = 0; j < l; j++) {
			sgo = s_theOrbits_ij(i, j);
			o_len = sgo->s_o_len_i();
			so_len = sgo->s_so_len_i();
			split = o_len / so_len;
			split_orbits += split;
			total_split_orbits += split;
			}
		sprintf(str2 + strlen(str2), "%ld ", split_orbits);
		}
	sprintf(str2 + strlen(str2), "\\} %ld", total_split_orbits);

	sprintf(str3 + strlen(str3), "\\{ ");
	for (i = 0; i <= len; i++) {
		sprintf(str3 + strlen(str3), "%ld ", s_nb_subgroups_ii(i));
		}
	sprintf(str3 + strlen(str3), "\\} %ld", s_total_nb_subgroups_i());
	return OK;
}

INT SGL_OB::get_orbit_generators(VECTOR_OP Gen, INT type, void *data)
{
	INT nb_layers, nb_orbits, i, j, l, k;
	SGO_OP sgo;
	
	l = s_total_nb_orbits_i();
	Gen->m_il(l);
	k = 0;
	nb_layers = s_theOrbits()->s_li();
	for (i = 0; i < nb_layers; i++) {
		nb_orbits = s_theOrbits_i(i)->s_li();
		for (j = 0; j < nb_orbits; j++) {
			sgo = s_theOrbits_ij(i, j);
			sgo->generators(this, (VECTOR_OP) Gen->s_i(k), type, data);
			k++;
			}
		}
	return OK;
}

INT SGL_OB::add_group(VECTOR_OP U_gen, VECTOR_OP U, INT f_verbose, INT type, void *data)
{
	SGO_OB U_orbit;
	SGO_OP pU_orbit;
	VECTOR_OB N_gen, N, U_zidx, N_zidx;
	INT erg = OK, f_found, f_selfnormalizing;
	INT k, nb_groups;
	INT Ulayer, Uorbit, Urep;
	INT Nlayer, Norbit, Nrep;
	INT f_very_verbose = FALSE;
	
	if (f_very_verbose) {
		printf("SGL::add_group():\n");
		}
	calc_zidx(U, &U_zidx);
	Ulayer = nb_primes(U->s_li());
	find_group(&U_zidx, Ulayer, &Uorbit, &Urep, &f_found);
	if (f_found)
		return error("SGL::add_group(): U found in lattice");
	if (U_orbit.Init(U_gen, U, &N_gen, &N, this, f_verbose, type, data) != OK) {
		printf("SGL_OB::add_group() error in U_orbit.Init(U_gen)\n");
		printf("U_gen = ");
		U_gen->println();
		fflush(stdout);
		return ERROR;
		}
	if (N.s_li() == U->s_li())
		f_selfnormalizing = TRUE;
	else
		f_selfnormalizing = FALSE;
	if (f_selfnormalizing) {
		Nlayer = Ulayer;
		Nrep = U_orbit.s_first_i();
		}
	else {
		calc_zidx(&N, &N_zidx);
		Nlayer = nb_primes(N.s_li());
		while (TRUE) {
			find_group(&N_zidx, Nlayer, &Norbit, &Nrep, &f_found);
			if (f_found)
				break;
			if (add_group(&N_gen, &N, f_verbose, type, data) != OK) {
				printf("SGL_OB::add_group() error in add_group(N_gen, N)\n");
				printf("N_gen = ");
				N_gen.println();
				fflush(stdout);
				return ERROR;
				}
			/* we must search the normalizer again, 
			 * it needs not to be an orbit leader ! */
			}
		if (f_very_verbose) {
			printf("SGL::add_group(): found normalizer at %ld / %ld / %ld\n", 
				Nlayer, Norbit, Nrep);
			}
		}
	/* U_orbit anfuegen: */
	/* nb_groups = U_orbit.s_Zidx()->s_li(); */
	nb_groups = U_orbit.s_o_len_i();
	s_theOrbits_i(Ulayer)->inc();
	k = s_theOrbits_i(Ulayer)->s_li() - 1;
	U_orbit.swap(s_theOrbits_ij(Ulayer, k));
	pU_orbit = s_theOrbits_ij(Ulayer, k);
	s_nb_subgroups_i(Ulayer)->m_i(s_nb_subgroups_ii(Ulayer) + nb_groups);
	s_total_nb_subgroups()->m_i(s_total_nb_subgroups_i() + nb_groups);
	s_total_nb_orbits()->inc();
	if (f_selfnormalizing) {
		Uorbit = s_theOrbits_i(Ulayer)->s_li() - 1;
		Norbit = Uorbit;
		}
	
	/* Nreps berechnen: */
	pU_orbit->s_Nlayer()->m_i(Nlayer);
	pU_orbit->s_Norbit()->m_i(Norbit);
	pU_orbit->calc_nreps(s_theOrbits_ij(Nlayer, Norbit), Nrep, f_very_verbose);
	if (f_verbose) {
		PrintStatistics();
		fflush(stdout);
		}
	return erg;
}

INT SGL_OB::all_layers(INT f_verbose, INT type, void *data)
{
	INT i, len;
	
	len = s_nb_Layers_i();
	for (i = 0; i < len; i++) {
		if (extend_layer(i, f_verbose, 
			type, data) != OK)
			return error("SGL::all_layers(): error in extend_layer()");
		}
	return OK;
}

INT SGL_OB::extend_layer(INT layer, INT f_verbose, INT type, void *data)
{
	VECTOR_OB Gamma;
	INT i, len;
	
	len = s_theOrbits_i(layer)->s_li();
	for (i = 0; i < len; i++) {
		if (f_verbose) {
			printf("SGL::extend_layer() layer = %ld orbit = %ld of %ld\n", 
				layer, i + 1, len);
			fflush(stdout);
			}
		extend_init_gamma(layer, i, &Gamma);
		if (extend_group(layer, i, &Gamma, f_verbose, type, data) != OK)
			return error("SGL::extend_layer(): error in extend_group()");
		}
	return OK;
}

INT SGL_OB::extend_group(INT layer, INT orbit, VECTOR_OP Gamma, INT f_verbose, 
	INT type, void *data)
{
	SGO_OP U_orbit, V_orbit;
	VECTOR_OP U_Zidx, V1_Zidx;
	VECTOR_OB U_gen, U;
	VECTOR_OB V_gen, V, V_zidx, bad_list;
	INTEGER_OB bad_list_len;
	SYM_OB g1;
	INT i, k, len, trial;
	INT Vlayer, Vorbit, Vrep;
	INT f_found, rel_order;
	INT l, len1, f_subseteq;
	INT f_found_bad_element;
	INT f_very_verbose = FALSE;
	
	U_orbit = s_theOrbits_ij(layer, orbit);
	U_Zidx = U_orbit->s_Zidx_i(U_orbit->s_first_i());
	Zidx2Z(U_orbit->s_gen_Zidx(), &U_gen);
	if (U_gen.s_li() == 0) {
		U.m_il(1);
		s_G_gen_i(0)->copy(U.s_i(0));
		do_one(U.s_i(0), type, data);
		}
	else {
		U_gen.copy(&U);
		gruppen_elemente1(&U, type, data);
		}
	if (sylow_normalizing_zuppos(layer, orbit, U.s_li(), Gamma, 
		s_f_very_verbose_i(), type, data) != OK) {
		printf("SGL::extend_group(): error in sylow_normalizing_zuppos()\n");
		}
	bad_list.m_il(VECTOR_OVERSIZE);
	bad_list_len.m_i(0);

	trial = 0;
	if (f_verbose) 
		printf("SGL::externd_group(): Gamma len = %ld\n", Gamma->s_li());
	while (Gamma->s_li()) {
		s_Z_i(Gamma->s_ii(0))->copy(&g1);
		for (k = 1; k < Gamma->s_li(); k++)
			Gamma->s_i(k - 1)->swap(Gamma->s_i(k));
		Gamma->s_l()->dec();
		U_gen.copy(&V_gen);
		U.copy(&V);
		dimino_extend(&V, &V_gen, &g1, &bad_list, bad_list_len.s_i(), 
			&f_found_bad_element, type, data);
		if (!f_found_bad_element) {
			/* printf("|V| = %ld\n", V.s_li()); */
			trial++;
			Vlayer = nb_primes(V.s_li());
			calc_zidx(&V, &V_zidx);
			find_group(&V_zidx, Vlayer, &Vorbit, &Vrep, &f_found);
			if (!f_found) {
				rel_order = V.s_li() / U.s_li();
				if (s_f_very_verbose_i()) {
					printf("trial = %ld\n", trial);
					/* printf("Gamma len = %ld\n", Gamma->s_li()); */
					fflush(stdout);
					}
				if (add_group(&V_gen, &V, f_verbose, type, data) != OK) {
					printf("SGL::extend_group(): error in add_group() "
						"layer = %ld orbit = %ld\n", layer, orbit);
					fflush(stdout);
					return error("stop");
					}
				if (nb_primes(rel_order) == 1) {
					V_orbit = s_theOrbits_ij(layer + 1, 
						s_theOrbits_i(layer + 1)->s_li() - 1);
					/* len1 = V_orbit->s_Zidx()->s_li(); */
					len1 = V_orbit->s_o_len_i();
					for (l = 0; l < len1; l++) {
						V1_Zidx = V_orbit->s_Zidx_i(l);
						U_Zidx->subseteq(V1_Zidx, 
							U_Zidx->s_li() /* lenA */, 
							V1_Zidx->s_li() /* lenB */, 
							TRUE /* f_ascending */, 
							&f_subseteq);
						if (f_subseteq) {
							Gamma->v_minus_v_ip(V1_Zidx);
							}
						}
					}
				trial = 0;
				if (f_verbose) {
					PrintStatistics();
					fflush(stdout);
					}
				}
			} /* if (!f_found_bad_element) */
		/* add g1 to bad_list: */
		if (bad_list.append_element(&bad_list_len, &g1) != OK)
			return error("SGL::extend_group(): error in bad_list.append_element()");
		}
	if (s_f_very_verbose_i()) {
		printf("trial = %ld\n", trial);
		fflush(stdout);
		}
	return OK;
}

INT SGL_OB::extend_init_gamma(INT layer, INT orbit, VECTOR_OP Gamma)
{
	SGO_OP U_orbit, V_orbit;
	VECTOR_OP U_Zidx, V_Zidx;
	INT nb_zuppos, k, len, l, len1, f_subseteq;
	
	nb_zuppos = s_Z()->s_li();
	Gamma->m_il(nb_zuppos);
	for (k = 0; k < nb_zuppos; k++)
		Gamma->m_ii(k, k);
	U_orbit = s_theOrbits_ij(layer, orbit);
	U_Zidx = U_orbit->s_Zidx_i(U_orbit->s_first_i());
	Gamma->v_minus_v_ip(U_Zidx);
	
	/* avoid subgroups already in the lattice: */
	len = s_theOrbits_i(layer + 1)->s_li();
	for (k = 0; k < len; k++) {
		V_orbit = s_theOrbits_ij(layer + 1, k);
		/* len1 = V_orbit->s_Zidx()->s_li(); */
		len1 = V_orbit->s_o_len_i();
		for (l = 0; l < len1; l++) {
			V_Zidx = V_orbit->s_Zidx_i(l);
			U_Zidx->subseteq(V_Zidx, 
				U_Zidx->s_li() /* lenA */, 
				V_Zidx->s_li() /* lenB */, 
				TRUE /* f_ascending */, &f_subseteq);
			if (f_subseteq) {
				Gamma->v_minus_v_ip(V_Zidx);
				}
			}
		}
	return OK;
}

INT SGL_OB::sylow_normalizing_zuppos(INT layer, INT orbit, 
	INT U_go, VECTOR_OP Gamma, INT f_verbose, 
	INT type, void *data)
{
	INTEGER_OB int_ob;
	SGO_OP U_orbit;
	VECTOR_OP U_Zidx;
	VECTOR_OB gp, ge, up, ue, Gamma1;
	VECTOR_OB P_gen, P_Zidx, P;
	INT i, i0, nb_primes_G, p;
	INT j, k, z, idx, f_found;
	INT max_p_idx, ny_p0, ny_p1;
	PERMUTATION_OP zperm;
	INT pz, pz1, ro;
	BYTE str[256];
	BYTE str1[256];
	BYTE str2[256];

	U_orbit = s_theOrbits_ij(layer, orbit);
	U_Zidx = U_orbit->s_Zidx_i(U_orbit->s_first_i());
	factor_integer(s_go_i(), &gp, &ge);
	factor_integer(U_go, &up, &ue);
	if (f_verbose) {
		sprintf(str, "U order = %ld = ", U_go);
		print_factorization(&up, &ue, str + strlen(str));
		printf("%s\n", str);
		}
	Gamma1.m_il(0);
	nb_primes_G = gp.s_li();
	i0 = 0; /* naechste Primzahl in up */
	for (i = 0; i < nb_primes_G; i++) {
		p = gp.s_ii(i);
		if (i0 < up.s_li() && up.s_ii(i0) == p) {
			/* U hat eine p - Sylowgruppe: */
			ny_p0 = ue.s_ii(i0);
			if (find_max_p_element(U_Zidx, p, 
				ny_p0, &max_p_idx, &ny_p1, 
				type, data) != OK) {
				str1[0] = 0;
				str2[0] = 0;
				print_factorization(&gp, &ge, str1);
				print_factorization(&up, &ue, str2);
				printf("s_go_i() = %ld = %s U_go = %ld = %s\n", 
					s_go_i(), str1, U_go, str2);
				
				printf("SGL::sylow_normalizing_zuppos(): "
					"error in find_max_p_element()\n");
				return ERROR;
				}
			z = U_Zidx->s_ii(max_p_idx);
			if (f_verbose) {
				printf("found max p element: "
					"p = %ld ny_p = %ld ny_p1 = %ld "
					"element = ", p, ny_p0, ny_p1);
				s_Z_i(z)->println();
				}
			P_gen.m_il(1);
			s_Z_i(z)->copy(P_gen.s_i(0));
			P_gen.copy(&P);
			gruppen_elemente1(&P, type, data);
			P_Zidx.m_il(1);
			P_Zidx.m_ii(0, z);
			while (ny_p1 < ny_p0) {
				find_normalizing_p_zuppo(&P_Zidx, &P, U_Zidx, p, 
					&idx, f_verbose, type, data);
				z = U_Zidx->s_ii(idx);
				if (f_verbose) {
					printf("found normalizing zuppo:");
					s_Z_i(z)->println();
					}
				dimino_extend_normal(&P, &P_gen, s_Z_i(z), type, data);
				ny_p1 = ny_p(P.s_li(), p);
				if (f_verbose) {
					printf("new P group order: %ld "
						"ny_p = %ld\n", P.s_li(), ny_p1);
					}
				if (!is_prime_power(P.s_li(), p))
					return error("SGL::sylow_normalizing_zuppos():"
						"order of P not a prime power of p");
				calc_zidx(&P, &P_Zidx);
				}
			if (f_verbose) {
				printf("found p - Sylow subgroup (p = %ld):\n", p);
				/* P.println(); */
				print_zuppo_vector(&P_Zidx, 
					FALSE /* f_vertically */, 
					type, data);
				}
			/* uebernehme jetzt alle p - zuppos 
			 * aus Gamma nach Gamma1, 
			 * welche P normalisieren und 
			 * Relativordnung p - Potenz haben: */
			for (j = 0; j < Gamma->s_li(); j++) {
				z = Gamma->s_ii(j);
				if (s_Zprime_ii(z) != p)
					continue;
				zperm = s_zuppos_on_zuppos_i(z);
				for (k = 0; k < P_Zidx.s_li(); k++) {
					pz = P_Zidx.s_ii(k);
					pz1 = zperm->s_ii(pz) - 1;
					int_ob.m_i(pz1);
					if (P_Zidx.search(P_Zidx.s_li(), TRUE, 
						&int_ob, &idx, &f_found) != OK)
						return error("SGL::sylow_normalizing_zuppos(): error in search()");
					if (!f_found)
						break;
					}
				if (k < P_Zidx.s_li())
					continue; /* next j */
		
				/* Test, ob Relativordnung p Potenz ist: */
				RelativeOrder(&P, s_Z_i(z), &ro, type, data);
				if (!is_prime_power(ro, p))
					continue;
				
				/* Zuppo z in Gamma1 einfuegen: */
				Gamma1.inc();
				if (Gamma1.search_and_insert_int(
					Gamma1.s_li() - 1, z) != OK)
					return error("SGL::sylow_normalizing_zuppos(): z found in Gamma1");

				} /* next j */
			i0++;
			}
		else {
			/* keine p - Sylowgruppe in U;
			 * alle p - Zuppos aus Gamma muessen 
			 * uebernommen werden: */
			for (j = 0; j < Gamma->s_li(); j++) {
				z = Gamma->s_ii(j);
				if (s_Zprime_ii(z) != p)
					continue;
				Gamma1.inc();
				if (Gamma1.search_and_insert_int(
					Gamma1.s_li() - 1, z) != OK)
					return error("SGL::sylow_normalizing_zuppos(): z found in Gamma1");
				}
			}
		}
	if (s_f_verbose_i()) {
		printf("|Gamma1| = %ld |Gamma| = %ld\n", Gamma1.s_li(), Gamma->s_li());
#if 0
		printf("Gamma1:\n");
		print_zuppo_vector(&Gamma1, 
			FALSE /* f_vertically */, 
			type, data);
#endif
		}
	Gamma1.swap(Gamma);
	return OK;
}

INT SGL_OB::find_normalizing_p_zuppo(VECTOR_OP U_Zidx, VECTOR_OP U, 
	VECTOR_OP G_Zidx, INT p, INT *idx, INT f_verbose, INT type, void *data)
/* liefert p - Zuppo Element aus G_Zidx 
 * ausserhalb von U_Zidx, 
 * welches U_Zidx normalisiert 
 * und Relativordnung p bzgl. U hat;
 * in idx steht dessen Position in G_Zidx. */
{
	PERMUTATION_OP gperm;
	INTEGER_OB int_ob;
	INT i, j, Ulen, Glen, gz;
	INT idx1, uz, uz1, f_found, ro;

	Glen = G_Zidx->s_li();
	Ulen = U_Zidx->s_li();
	for (i = 0; i < Glen; i++) {
		gz = G_Zidx->s_ii(i);
		if (s_Zprime_ii(gz) != p)
			continue;
		int_ob.m_i(gz);
		if (U_Zidx->search(Ulen, TRUE, 
			&int_ob, &idx1, &f_found) != OK)
			return error("SGL::find_normalizing_p_zuppo(): "
				"error in search()");
		if (f_found)
			continue;
		/* Test, ob gz U_Zidx normalisiert: */
		gperm = s_zuppos_on_zuppos_i(gz);
		for (j = 0; j < Ulen; j++) {
			uz = U_Zidx->s_ii(j);
			uz1 = gperm->s_ii(uz) - 1;
			int_ob.m_i(uz1);
			if (U_Zidx->search(Ulen, TRUE, 
				&int_ob, &idx1, &f_found) != OK)
				return error("SGL::find_normalizing_p_zuppo(): "
					"error in search()");
			if (!f_found)
				break;
			}
		if (j < Ulen)
			continue;
		
		/* Test, ob Relativordnung 
		 * p Potenz ist: */
		RelativeOrder(U, s_Z_i(gz), &ro, type, data);
		if (!is_prime_power(ro, p))
			continue;
		*idx = i;
		return OK;
		}
	return error("SGL::find_normalizing_p_zuppo(): no element found");
}

INT SGL_OB::find_max_p_element(VECTOR_OP Z_idx, INT p, INT ny_p0, 
	INT *idx, INT *ny_p1, INT type, void *data)
/* sucht maximalen p - Zuppo in Z_idx; 
 * Dessen Position in Z_idx 
 * wird in idx zurueckgegeben, 
 * sein p - Exponent in ny_p1. */
{
	INT i, len, z, ny_p_max = 0;
	
	*idx = -1;
	len = Z_idx->s_li();
	for (i = 0; i < len; i++) {
		z = Z_idx->s_ii(i);
		if (s_Zprime_ii(z) != p)
			continue;
		if (s_Zk_ii(z) > ny_p_max) {
			*idx = i;
			ny_p_max = s_Zk_ii(z);
			if (ny_p_max == ny_p0) {
				*ny_p1 = ny_p_max;
				return OK;
				}
			}
		}
	if (*idx == -1) {
		printf("Z_idx = \n");
		print_zuppo_vector(Z_idx, 
			FALSE /* f_vertically */, 
			type, data);
		printf("p = %ld, ny_p = %ld\n", p, ny_p0);
		printf("SGL::find_max_p_element(): no max p element found\n");
		return ERROR;
		}
	*ny_p1 = ny_p_max;
	return OK;
}

INT SGL_OB::find_group(VECTOR_OP Zidx, INT layer, INT *orbit, INT *rep, INT *f_found)
{
	SGO_OP Orbit;
	INT l, len;
	
	len = s_theOrbits_i(layer)->s_li();
	for (l = 0; l < len; l++) {
		Orbit = s_theOrbits_ij(layer, l);
		Orbit->find_group(Zidx, rep, f_found);
		if (*f_found) {
			*orbit = l;
			return OK;
			}
		}
	*f_found = FALSE;
	return OK;
}

INT SGL_OB::representation_on_zuppos(SYM_OP g, PERMUTATION_OP p, INT type, void *data)
{
	PERMUTATION_OB pv, q;
	INT i, j, len, idx, f_found;
	SYM_OP z;
	SYM_OB z1, z2, g_inv;
	
	len = s_Z()->s_li();
	p->m_il(len);
	pv.m_il(len);
	for (i = 0; i < len; i++) {
		p->m_ii(i, -1);
		pv.m_ii(i, -1);
		}
	do_invers(g, &g_inv, type, data);
	for (i = 0; i < len; i++) {
		z = s_Z_i(i);
		/* z2 := g_inv * z * g */
		do_mult(&g_inv, z, &z1, type, data);
		do_mult(&z1, g, &z2, type, data);
		if (s_PP()->search(s_PP()->s_li(), TRUE, &z2, &idx, &f_found) != OK)
			return error("SGL::representation_on_zuppos(): "
				"error in search(PP, z2)");
		if (!f_found)
			return error("SGL::representation_on_zuppos(): "
				"z2 not found in PP");
		idx--;
		j = s_PPidx_ii(idx);
		/* now: i -> j under conjugation by g. */
		if (pv.s_ii(j) != -1)
			return error("SGL::representation_on_zuppos(): "
				"not a permutation");
		p->m_ii(i, j + 1);
		pv.m_ii(j, i + 1);
		}
	p->mult(&pv, &q);
	if (!q.einsp())
		return error("SGL::representation_on_zuppos(): !q->einsp()");
	return OK;
}

INT SGL_OB::Aut_representation_on_zuppos(VECTOR_OP G, 
	PERMUTATION_OP aut, PERMUTATION_OP p)
{
	PERMUTATION_OB pv, q;
	INT i, j, len, idx, f_found;
	INTEGER_OP z;
	SYM_OB z2;
	
	len = s_Z()->s_li();
	p->m_il(len);
	pv.m_il(len);
	for (i = 0; i < len; i++) {
		p->m_ii(i, -1);
		pv.m_ii(i, -1);
		}
	for (i = 0; i < len; i++) {
		z = (INTEGER_OP) s_Z_i(i);
		if (z->s_obj_k() != INTEGER) {
			printf("zuppo %ld z = ", i);
			z->println();
			fflush(stdout);
			// return error("SGL::Aut_representation_on_zuppos(): z->s_obj_k() != INTEGER");
			G->search(G->s_li(), TRUE /* f_ascending */, z, &idx, &f_found);
			if (!f_found)
				return error("SGL_OB::Aut_representation_on_zuppos(): z not found in G");
			idx--;
			printf("is element %ld in G\n", idx);
			fflush(stdout);
			j = idx;
			j = aut->s_ii(j) - 1;
			G->s_i(j)->copy(&z2);
			}
		else {
			j = z->s_i();
			((INTEGER_OP) &z2)->m_i(aut->s_ii(j) - 1);
			}
		/* z2 := z^aut */
		if (s_PP()->search(s_PP()->s_li(), TRUE, &z2, &idx, &f_found) != OK)
			return error("SGL::Aut_representation_on_zuppos(): "
				"error in search(PP, z2)");
		if (!f_found) {
			printf("PP=");
			s_PP()->println();
			printf("z2=");
			z2.println();
			printf("G=");
			G->println();
			printf("aut=");
			aut->println();
			printf("Aut_gen=");
			s_Aut_gen()->println();
			return error("SGL::Aut_representation_on_zuppos(): "
				"z2 not found in PP");
			}
		idx--;
		j = s_PPidx_ii(idx);
		/* now: i -> j under conjugation by aut. */
		if (pv.s_ii(j) != -1)
			return error("SGL::Aut_representation_on_zuppos(): "
				"not a permutation");
		p->m_ii(i, j + 1);
		pv.m_ii(j, i + 1);
		}
	p->mult(&pv, &q);
	if (!q.einsp())
		return error("SGL::Aut_representation_on_zuppos(): !q->einsp()");
	return OK;
}

INT SGL_OB::conjugate_zuppos_by_perm(VECTOR_OP z_idx, PERMUTATION_OP g, VECTOR_OP z_idx1)
{
	INT l, len, z, z1;
	
	len = z_idx->s_li();
	z_idx1->m_il(len);
	for (l = 0; l < len; l++) {
		z = z_idx->s_ii(l);
		z1 = g->s_ii(z) - 1;
		z_idx1->m_ii(l, z1);
		}
	z_idx1->quicksort(len, TRUE /* f_ascending */);
	return OK;
}

INT SGL_OB::conjugate_z_idx(VECTOR_OP z_idx, SYM_OP g, VECTOR_OP z_idx1, 
	INT type, void *data)
{
	INT k, l, len, zidx, zidx1, idx, f_found;
	SYM_OP z;
	SYM_OB z1, z2, g_inv;
	
	len = z_idx->s_li();
	z_idx1->m_il(len);
	do_invers(g, &g_inv, type, data);
	for (l = 0; l < len; l++) {
		zidx = z_idx->s_ii(l);
		z = s_Z_i(zidx);
		/* z2 := g_inv * z * g 
		 * Beachte: vormals: z2 := g * z * g_inv !!! */
		do_mult(&g_inv, z, &z1, type, data);
		do_mult(&z1, g, &z2, type, data);
		if (s_PP()->search(s_PP()->s_li(), TRUE, 
			&z2, &idx, &f_found) != OK)
			return error("SGL::conjugate_z_idx(): "
				"error in search(PP, z2)");
		if (!f_found)
			return error("SGL::conjugate_z_idx(): "
				"z2 not found in PP");
		idx--;
		zidx1 = s_PPidx_ii(idx);
		if (z_idx1->search_and_insert_int(l, zidx1) != OK)
			return error("SGL::conjugate_z_idx(): "
				"found zidx1 in z_idx1");
		}
	return OK;
}

INT SGL_OB::Aut_conjugate_z_idx(VECTOR_OP z_idx, PERMUTATION_OP aut, VECTOR_OP z_idx1)
{
	INT j, k, l, len, zidx, zidx1, idx, f_found;
	INTEGER_OP z;
	INTEGER_OB z2;
	
	len = z_idx->s_li();
	z_idx1->m_il(len);
	for (l = 0; l < len; l++) {
		zidx = z_idx->s_ii(l);
		z = (INTEGER_OP) s_Z_i(zidx);
		if (z->s_obj_k() != INTEGER)
			return error("SGL::Aut_conjugate_z_idx(): "
				"z->s_obj_k() != INTEGER");
		j = z->s_i();
		/* z2 := z^aut */
		z2.m_i(aut->s_ii(j) - 1);
		if (s_PP()->search(s_PP()->s_li(), TRUE, &z2, &idx, &f_found) != OK)
			return error("SGL::Aut_conjugate_z_idx(): "
				"error in search(PP, z2)");
		if (!f_found)
			return error("SGL::Aut_conjugate_z_idx(): "
			"z2 not found in PP");
		idx--;
		zidx1 = s_PPidx_ii(idx);
		if (z_idx1->search_and_insert_int(l, zidx1) != OK)
			return error("SGL::Aut_conjugate_z_idx(): "
				"found zidx1 in z_idx1");
		}
	return OK;
}

#endif /* SGL_TRUE */

