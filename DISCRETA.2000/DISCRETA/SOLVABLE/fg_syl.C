/* fg_syl.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SOLVABLE_TRUE

#include <DISCRETA/solvable.h>
#include <DISCRETA/sgl.h>

#define MAX_NW 64

#undef SGL_DEBUG

#if TEXDOCU
INT FG_OB::sgl(INT f_v, INT f_vv, INT f_vvv)
#endif
{
	VECTOR_OB G_gen, G;
	ZE_OP ze;
	SGL_OB L;
	INT i, g_i;
	INT sz1, sz2;
	
	
	if (f_v || 1) {
		printf("in FG::sgl(group = %s", s_label_s());
		printf("\n");
		Print_gen();
		fflush(stdout);
		}
	G_gen.m_il(s_nb_ze_i());
	for (i = 0; i < s_nb_ze_i(); i++) {
		ze = s_ze_i(i);
		g_i = ze->s_n0_i();
		G_gen.m_ii(i, g_i);
		}
	{
		VECTOR_OB G1;
		
		if (f_vvv || 1) {
			printf("G_gen = \n");
			G_gen.println();
			printf("calling vec_to_vec_pp()\n");
			fflush(stdout);
			}
		vec_to_vec_pp(&G_gen, &G1, 
			DO_TYPE_FG /* type */, this /* data */);
		G1.swap(&G_gen);
		if (f_vvv || 1) {
			printf("finished !\nG_gen =\n");
			G_gen.println();
			fflush(stdout);
			}
	}
	

	G.m_il(s_n_i());
	for (i = 0; i < s_n_i(); i++) {
		G.m_ii(i, i);
		}
	
	if (f_vvv) {
		printf("L.init: theG=\n");
		fflush(stdout);
		s_theG()->Print();
		}
	if (!s_f_has_aut_i()) {
		printf("FG_OB::sgl() aut group not calculated, we use inner automorphisms\n");
		fflush(stdout);
		L.Init(&G_gen, &G, 
			f_vv /* f_verbose */, 
			f_vvv /* f_very_verbose */, 
			DO_TYPE_FG, this);
		}
	else {
		VECTOR_OB Aut_gen;
		PERMUTATION_OP p;
		BYTE str[1024];
		INT i, l;
		
		get_aut_generators(&Aut_gen);
		l = Aut_gen.s_li();
		for (i = 0; i < l; i++) {
			p = (PERMUTATION_OP) Aut_gen.s_i(i);
			str[0] = 0;
			sprint_aut_on_base(p, str);
			printf("%ld-th automorphism: %s\n", i, str);
			}
		if (f_v) {
			printf("found %ld generators of the automorphism group\n", 
				Aut_gen.s_li());
			printf("calling SGL::Init_with_Aut_generators()\n");
			fflush(stdout);
			}
		L.Init_with_Aut_generators(
			&G_gen, &G, &Aut_gen, 
			f_vv /* f_verbose */, 
			f_vvv /* f_very_verbose */, DO_TYPE_FG, this);
		}
	if (f_vvv) {
		printf("Z_info:\n");
		fflush(stdout);
		L.latex_Z_info(stdout, 
			"FG-group" /* tex_group_name */, DO_TYPE_FG, this);
		printf("calling L.all_layers:\n");
		fflush(stdout);
		}
	if (f_v) {
		printf("calling SGL::all_layers()...\n");
		fflush(stdout);
		}
	L.all_layers(f_vv /* f_verbose */, DO_TYPE_FG, this);
	sz1 = L.calc_size_on_file();
	printf("SGL size on file = %ld\n", sz1);
	fflush(stdout);
	if (f_v || TRUE) {
		printf("calling SGL::shrink()...\n");
		fflush(stdout);
		}
	L.shrink();
	sz2 = L.calc_size_on_file();
	printf("SGL size on file = %ld\n", sz2);
	fflush(stdout);
	L.swap(s_sgl());
	s_f_has_sgl()->m_i(TRUE);
	/* s_C()->add_object("SGL", &L); */
	if (f_v) {
		printf(")\n");
		fflush(stdout);
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::sylow_type(INT f_verbose)
#endif
{
	VECTOR_OB G, P, P_gen;
	VECTOR_OB N, N_gen, primes, exponents;
	BYTE str[256];
	INT erg = OK, i, p;
	INT G_order, N_order, p_type;
	
	erg += factor_integer(
	s_n_i(), &primes, &exponents);
	str[0] = 0;
	print_factorization(&primes, &exponents, str);
	if (f_verbose)
		printf("%ld = %s\n", s_n_i(), str);
	G_order = s_n_i();
	G.m_il(G_order);
	for (i = 0; i < G_order; i++)
		G.m_ii(i, i);
	s_sylow_type()->m_il(primes.s_li());
	for (i = 0; i < primes.s_li(); i++) {
		p = primes.s_ii(i);
		erg += p_Sylow(p, &G, &P, &P_gen, f_verbose);
		erg += p_Normalizer(0, &P, &P_gen, &G, &N, &N_gen);
		N_order = N.s_li();
		p_type = G_order / N_order;
		if (f_verbose)
			printf("G_order = %ld N_order = %ld p_type = %ld\n", 
				G_order, N_order, p_type);
		s_sylow_type()->m_ii(i, p_type);
		}
	s_f_has_sylow_type()->m_i(TRUE);
	if (f_verbose) {
		printf("sylow_type: ");
		s_sylow_type()->println();
		}
	return erg;
}

#if TEXDOCU
INT FG_OB::max_p_Element(INT p, INT *idx, INT *the_k, INT *f_found)
#endif
{
	INT i, len, nyp, order, prime, k_max, k;
	INT erg = OK;
	
	len = s_n_i();
	nyp = ny_p(len, p);
	k_max = 0;
	*f_found = FALSE;
	for (i = 0; i < len; i++) {
		erg += gt_order_if_prime_power(
			i, &order, &prime, &k);
		if (order == 0 || order == 1)
			/* Zusammengesetze Ordnung oder id */
			continue;
		if (prime != p)
			continue;
		if (k != 1) /* !!! min_p_Element instead */
			continue;
		if (k > k_max) {
			*idx = i;
			*the_k = k;
			k_max = k;
			*f_found = TRUE;
			if (k_max == nyp)
				return erg;
			}
		}
	return erg;
}

#if TEXDOCU
INT FG_OB::p_DiminoExtend(VECTOR_OP p, VECTOR_OP gen, INT g, 
	INT f_p, INT prime, INT *f_not_added)
#else
/* g wird an gen angefuegt (alle Vektoren ueber INT). */
#endif
{
	INTEGER_OB int_ob;
	VECTOR_OB p_tmp;
	VECTOR_OB reps;
	INT the_one = 0, g1, h, h1;
	INT erg = OK, i1, j, k, l, l1;
	INT idx, f_found, order;
	INT prime1, k1, f_p_element;
	
	gen->inc();
	gen->m_ii(gen->s_li() - 1, g);
	erg += p->copy(&p_tmp);
	erg += reps.m_il(1);
	reps.m_ii(0, the_one);
	/* the dimino induction step: 
	 * apply all generators to all coset reps. */
	for (j = 0; j < reps.s_li(); j++) {
		for (i1 = 0; i1 < gen->s_li(); i1++) {
			if (j == 0 && i1 < gen->s_li() - 1)
				continue;
			g1 = gen->s_ii(i1);
			h = gt_mult(reps.s_ii(j), g1);
			int_ob.m_i(h);
			erg += p->search(p->s_li(), 
				TRUE /* f_ascending */, 
				&int_ob, &idx, &f_found);
			if (f_found)
				continue;
			reps.inc();
			reps.m_ii(reps.s_li() - 1, h);
			l1 = p->s_li();
			p->v_realloc(l1 + p_tmp.s_li());
			for (k = 0; k < p_tmp.s_li(); k++) {
				h1 = gt_mult(p_tmp.s_ii(k), h);
				if (f_p) {
					erg += gt_order_if_prime_power(
						h1, &order, &prime1, &k1);
					f_p_element = TRUE;
					if (order == 0 || order == 1)
						f_p_element = FALSE;
					else if (prime1 != prime)
						f_p_element = FALSE;
					if (!f_p_element) {
						gen->dec();
						*f_not_added = TRUE;
						p_tmp.swap(p);
						return erg;
						}
					}
				/* p hat jetzt l1 + k benuetzte Elemente */
				int_ob.m_i(h1);
				erg += p->search(l1 + k, 
					TRUE /* f_ascending */, 
					&int_ob, &idx, &f_found);
				if (f_found) {
					printf("error in FG::p_DiminoExtend(): h1 found\n");
					return ERROR;
					}
				for (l = l1 + k; l > idx; l--) {
					p->s_i(l - 1)->swap(p->s_i(l));
					}
				p->m_ii(idx, h1);
				} /* next k */
			} /* next i1 */
		} /* next j */
	if (f_p) {
		*f_not_added = FALSE;
		}
	return erg;
}

#if TEXDOCU
INT FG_OB::DiminoExtend(VECTOR_OP p, VECTOR_OP gen, INT g)
#else
/* g wird an gen angefuegt. */
#endif
{
	INT erg = OK, i;
	
	erg += p_DiminoExtend(p, gen, g, 
		FALSE /* f_p */, 
		0 /* prime */, &i /* f_not_added */);
	return erg;
}

#if TEXDOCU
INT FG_OB::Dimino(VECTOR_OP p)
#else
/* p Vektor von Generatoren 
 * einer Permutationsgruppe;
 * liefert im Vektor p in sortierter 
 * Reihenfolge alle Gruppenelemente. */
#endif
{
	VECTOR_OB G, G1;
	INT the_one = 0;
	INT erg = OK, i, nb_G, g;
	
	erg += p->copy(&G);
	nb_G = G.s_li();
	if (nb_G == 0)
		return erg;
	G1.m_il(0);
	p->m_il(1);
	p->m_ii(0, the_one);
	for (i = 0; i < nb_G; i++) {
		g = G.s_ii(i);
		erg += DiminoExtend(p, &G1, g);
		} /* next i */
	return erg;
}

#if TEXDOCU
INT FG_OB::p_Normalizer(INT p, 
	VECTOR_OP H, VECTOR_OP H_gen, 
	VECTOR_OP G, VECTOR_OP N, VECTOR_OP N_gen)
#else
/* p = 0: ordinary normalizer */
#endif
{
	INTEGER_OB int_ob;
	VECTOR_OB Gamma;
	INT s1, s2, g, gv;
	INT erg = OK, i, j, l;
	INT idx, f_found, order, prime, k;
	INT f_del_coset, f_not_added, f_p;
	
	if (p)
		f_p = TRUE;
	else
		f_p = FALSE;
	G->copy(&Gamma);
	H->copy(N);
	H_gen->copy(N_gen);
	Gamma.v_minus_v_ip(N);
	while (Gamma.s_li()) {
		g = Gamma.s_ii(0);
		f_del_coset = FALSE;
		if (p) {
			erg += gt_order_if_prime_power(
				g, &order, &prime, &k);
			if (order == 0 || order == 1)
				f_del_coset = TRUE;
			else if (prime != p)
				f_del_coset = TRUE;
			}
		if (!f_del_coset) {
			gv = gt_inv(g);
			for (i = 0; i < H_gen->s_li(); i++) {
				/* s2 := gv * h_i * g */
				s1 = gt_mult(gv, H_gen->s_ii(i));
				s2 = gt_mult(s1, g);
				int_ob.m_i(s2);
				erg += H->search(H->s_li(), 
					TRUE /* f_ascending */, 
					&int_ob, &idx, &f_found);
				if (!f_found)
					break;
				}
			if (i < H_gen->s_li())
				f_del_coset = TRUE;
			}
		if (!f_del_coset) { /* g normalizes H */
			p_DiminoExtend(N, N_gen, g, f_p, 
				p /* prime */, &f_not_added);
			if (f_p && f_not_added) {
				f_del_coset = TRUE;
				}
			else {
				/* Gamma -= N: */
				Gamma.v_minus_v_ip(N);
				}
			}
		if (f_del_coset) {
			/* g does not normalize H or 
			 * <N,g> not p group. */
			/* Gamma -= N g: */
			l = Gamma.s_li();
			for (j = 0; j < N->s_li(); j++) {
				s1 = gt_mult(N->s_ii(j), g);
				int_ob.m_i(s1);
				erg += Gamma.search(l, 
					TRUE /* f_ascending */, 
					&int_ob, &idx, &f_found);
				if (!f_found)
					continue;
				idx--;
				Gamma.v_del_ith2(l, idx);
				l--;
				}
			Gamma.v_shorten(l);
			}
		}
	return erg;
}

#if TEXDOCU
INT FG_OB::p_Sylow(INT p, VECTOR_OP G, 
	VECTOR_OP P, VECTOR_OP P_gen, INT f_verbose)
#endif
{
	INT i, len, nyp, idx, the_k, f_found, g;
	INT erg = OK;
	VECTOR_OB H, H_gen, N, N_gen;
	BYTE str[256];
	
	len = G->s_li();
	nyp = ny_p(len, p);
	erg += max_p_Element(p, &idx, &the_k, &f_found);
	if (!f_found) {
		return error("FG::p_Sylow()|no p element found");
		}
	if (f_verbose)
		printf("FG::p_Sylow()|len = %ld nyp = %ld; "
		"found a max_p element: p = %ld k = %ld\n", 
			len, nyp, p, the_k);
	str[0] = 0;
	g = G->s_ii(idx);
	if (f_verbose) {
		printf("%ld\n", g);
		fflush(stdout);
		}
	N_gen.m_il(1);
	N_gen.m_ii(0, idx);
	N_gen.copy(&N);
	Dimino(&N);
	if (f_verbose) {
		printf("nach Dimino:\n");
		N.println();
		}
	while (the_k != nyp) {
		N.swap(&H);
		N_gen.swap(&H_gen);
		N.freeself();
		N_gen.freeself();
		erg += p_Normalizer(p, &H, &H_gen, G, &N, &N_gen);
		if (f_verbose) {
			printf("nach p_Normalizer:\n");
			N.println();
			}
		the_k = ny_p(N.s_li(), p);
		if (f_verbose) {
			printf("normalizer: %ld hat nyp %ld\n", N.s_li(), the_k);
			fflush(stdout);
			}
		/* printf("N_gen = \n");
		orbit_p_vec_print(N_gen);
		fflush(stdout); */
		}
	N.swap(P);
	N_gen.swap(P_gen);
	return erg;
}

#endif /* SOLVABLE_TRUE */


