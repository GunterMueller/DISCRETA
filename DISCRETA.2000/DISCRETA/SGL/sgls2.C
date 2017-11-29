/* sgls2.C */

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

INT SGL_OB::zuppos(VECTOR_OP G, INT f_verbose, INT type, void *data)
{
	VECTOR_OP Z = s_Z(); 
	VECTOR_OP Zidx = s_Zidx();
	VECTOR_OP Zorder = s_Zorder();
	VECTOR_OP Zprime = s_Zprime();
	VECTOR_OP Zk = s_Zk();
	VECTOR_OP PP = s_PP();
	VECTOR_OP PPidx = s_PPidx();
	SYM_OB g1, g2;
	SYM_OP g;
	INT len, i, j;
	INT prime_power, prime, k, idx, f_found;
	INT nb_zuppos = 0, nb_PP = 0;
	BYTE *f_visited = NIL;
	
	len = G->s_li();
	f_visited = (BYTE *) my_malloc(len * sizeof(BYTE), "sgls2.C: zuppos");
	if (f_visited == NIL)
		return error("zuppos(): no memory");
	for (i = 0; i < len; i++)
		f_visited[i] = (BYTE) FALSE;
	Z->m_il(len);
	Zidx->m_il(len);
	Zorder->m_il(len);
	Zprime->m_il(len);
	Zk->m_il(len);
	PP->m_il(len);
	PPidx->m_il(len);
	for (i = 0; i < len; i++) {
		if (f_visited[i])
			continue;
		g = G->s_i(i);
		if (do_order_if_prime_power(g, &prime_power, &prime, &k, type, data) != OK)
			return error("zuppos(): error in do_order_if_prime_power()");
		if (prime_power == 1) { /* id */
			continue;
			}
		if (prime_power == 0) {
			/* not a prime power */
			continue;
			}
		
		/* Eintragen in Zuppos: 
		 * g / i / prime_power / prime / k: */
		g->copy(Z->s_i(nb_zuppos));
		Zidx->m_ii(nb_zuppos, i);
		Zorder->m_ii(nb_zuppos, prime_power);
		Zprime->m_ii(nb_zuppos, prime);
		Zk->m_ii(nb_zuppos, k);
		if (f_verbose) {
			printf("zuppo no %ld: ", nb_zuppos + 1);
			g->println();
			}
		
		g->copy(&g1);
		for (j = 1; j < prime_power; j++) {
			if (j > 1) {
				do_mult(g, &g1, &g2, type, data);
				g2.swap(&g1);
				/* jetzt g1 = g^j */
				}
			if ((j % prime) == 0)
				continue;
			if (G->search(len, TRUE, &g1, &idx, &f_found) != OK)
				return error("zuppos(): error in search()");
			if (!f_found)
				return error("zuppos(): g1 not f_found in G");
			idx--;
			f_visited[idx] = (BYTE) TRUE;

			if (PP->search(nb_PP, TRUE, &g1, &idx, &f_found) != OK)
				return error("zuppos(): error in PP->search(&g1)");
			if (f_found)
				return error("zuppos(): found element in PP");
			for (k = nb_PP; k > idx; k--) {
				PP->s_i(k)->swap(PP->s_i(k - 1));
				PPidx->s_i(k)->swap(PPidx->s_i(k - 1));
				}
			g1.copy(PP->s_i(idx));
			PPidx->m_ii(idx, nb_zuppos);
			nb_PP++;
			}
		nb_zuppos++;
		}
	Z->v_shorten(nb_zuppos);
	Zidx->v_shorten(nb_zuppos);
	Zorder->v_shorten(nb_zuppos);
	Zprime->v_shorten(nb_zuppos);
	Zk->v_shorten(nb_zuppos);
	PP->v_shorten(nb_PP);
	PPidx->v_shorten(nb_PP);
	if (s_f_verbose_i()) {
		printf("SGL::zuppos(): nb_zuppos = %ld nb_PP = %ld\n", nb_zuppos, nb_PP);
		printf("SGL::zuppos(): calculating zuppos on zuppos");
		fflush(stdout);
		}

	zuppos_on_zuppos(type, data);
	
	if (s_f_verbose_i()) {
		printf(".\n");
		fflush(stdout);
		}
	
l_exit:
	if (f_visited) {
		my_free(f_visited);
		f_visited = NIL;
		}
	return OK;
}

INT SGL_OB::zuppos_on_zuppos(INT type, void *data)
{
	INT nb_zuppos, k;
	
	nb_zuppos = s_Z()->s_li();
	s_zuppos_on_zuppos()->m_il(nb_zuppos);
	for (k = 0; k < nb_zuppos; k++) {
		representation_on_zuppos(s_Z_i(k), s_zuppos_on_zuppos_i(k), type, data);
		/* i -> g^{-1} i g for g := Z[k] 
		 * and i an arbitrary zuppo element */
		if (s_f_very_verbose_i()) {
			printf("zuppo %ld:\n", k);
			s_Z_i(k)->println();
			printf("as permutation on the zuppos:\n");
			s_zuppos_on_zuppos_i(k)->println();
			}
		}
	return OK;
}

INT SGL_OB::calc_G_gen_on_zuppos(INT type, void *data)
{
	INT k;
	
	s_G_gen_on_zuppos()->m_il(s_G_gen()->s_li());
	for (k = 0; k < s_G_gen()->s_li(); k++) {
		representation_on_zuppos(s_G_gen_i(k), s_G_gen_on_zuppos_i(k), type, data);
		if (s_f_very_verbose_i()) {
			printf("generator %ld:\n", k);
			s_G_gen_i(k)->println();
			printf("as permutation on the zuppos:\n");
			s_G_gen_on_zuppos_i(k)->println();
			}
		}
	return OK;
}

INT SGL_OB::calc_aut_on_zuppos(VECTOR_OP G)
{
	INT k, l;
	
	l = s_Aut_gen()->s_li();
	s_Aut_gen_on_zuppos()->m_il(l);
	for (k = 0; k < l; k++) {
		Aut_representation_on_zuppos(G, s_Aut_gen_i(k), s_Aut_gen_on_zuppos_i(k));
		if (s_f_very_verbose_i()) {
			printf("Aut generator %ld:\n", k);
			s_Aut_gen_i(k)->println();
			printf("as permutation on the zuppos:\n");
			s_Aut_gen_on_zuppos_i(k)->println();
			}
		}
	return OK;
}

INT SGL_OB::latex_Z_info(FILE *fp, BYTE *tex_group_name, INT type, void *data)
{
	INT i, j, k, Z_len, G_len;
	
	Z_len = s_Z()->s_li();
	G_len = s_G_gen()->s_li();
	fprintf(fp, "\\begin{center} {\\Large $%s$ } "
		"\\end{center}\\par\n", tex_group_name);
	fprintf(fp, "\\begin{center} {\\bf Zuppo info}: "
		"\\end{center}\\par\n");
	fprintf(fp, "\\[\n");
	fprintf(fp, "\\begin{array}{cccccc}\n");
	fprintf(fp, "\\mbox{no.} & \\mbox{Zuppo} & "
		"\\mbox{idx in G} & "
		"\\mbox{order}\\; = & \\mbox{p} & "
		"\\mbox{k} \\\\ \n");
	for (i = 0; i < Z_len; i++) {
		fprintf(fp, "%ld & ", i);
		do_latex(s_Z_i(i), fp, type, data);
		fprintf(fp, " & %ld & %ld & %ld & %ld\\\\\n", 
			s_Zidx_ii(i), s_Zorder_ii(i), 
			s_Zprime_ii(i), s_Zk_ii(i));
		}
	fprintf(fp, "\\end{array}\n");
	fprintf(fp, "\\]\n");
	fprintf(fp, "\\begin{center} "
		"{\\bf Zuppo conjugation table} \\\\\n");
	fprintf(fp, "$(i, j)$-entry is $z_i^{z_j}$ "
		"(first part of the table)\\\\\n");
	fprintf(fp, "$(i, j)$-entry is $z_i^{g_j}$ "
		"(second part):\\\\\\end{center}\n");
	fprintf(fp, "\\[\n");
	fprintf(fp, "\\begin{array}{r|");
	for (i = 0; i < Z_len; i++)
		fprintf(fp, "c");
	fprintf(fp, "||");
	for (i = 0; i < G_len; i++)
		fprintf(fp, "c");
	fprintf(fp, "}\n");
	fprintf(fp, " & ");
	for (i = 0; i < Z_len; i++) {
		fprintf(fp, " %ld ", i);
		if (i < Z_len - 1)
			fprintf(fp, " & ");
		}
	fprintf(fp, "\\\\ \\hline \n");
	for (i = 0; i < Z_len; i++) {
		fprintf(fp, " %ld & ", i);
		for (j = 0; j < Z_len; j++) {
			k = s_zuppos_on_zuppos_i(j)->s_ii(i) - 1;
			fprintf(fp, " %ld ", k);
			if (j < Z_len - 1)
				fprintf(fp, " & ");
			}
		fprintf(fp, " & ");
		for (j = 0; j < G_len; j++) {
			k = s_G_gen_on_zuppos_i(j)->s_ii(i) - 1;
			fprintf(fp, " %ld ", k);
			if (j < G_len - 1)
				fprintf(fp, " & ");
			}
		fprintf(fp, "\\\\\n");
		}
	fprintf(fp, "\\end{array}\n");
	fprintf(fp, "\\]\n");
	fflush(fp);
	return OK;
}

INT SGL_OB::print_zuppo_vector(VECTOR_OP V, INT f_vertically, INT type, void *data)
{
	INT i, z, len;

	len = V->s_li();
	printf("(");
	for (i = 0; i < len; i++) {
		z = V->s_ii(i);
		do_print(s_Z_i(z), type, data);
		if (i < len - 1) {
			printf(", ");
			if (f_vertically)
				printf("\n");
			}
		}
	printf(")\n");
	return OK;
}

INT SGL_OB::calc_zidx(VECTOR_OP H, VECTOR_OP H_zidx)
/* Berechnet Vektor H_zidx, 
 * welcher die Indizes (nach Z) 
 * der Zuppos enthaelt, 
 * welche in H auftreten. */
{
	SYM_OP z;
	INT k, l, len, idx, f_found, nb_H_zidx;
	
	len = s_Z()->s_li();
	H_zidx->m_il(len);
	nb_H_zidx = 0;
	for (l = 0; l < len; l++) {
		z = s_Z_i(l);
		if (H->search(H->s_li(), TRUE, z, &idx, &f_found) != OK)
			return error("SGL::calc_zidx(): error in search(H)");
		if (f_found) {
			if (H_zidx->search_and_insert_int(nb_H_zidx, l) != OK)
				return error("SGL::calc_zidx(): found zidx");
			nb_H_zidx++;
			}
		}
	H_zidx->v_shorten(nb_H_zidx);
	return OK;
}

INT SGL_OB::Z2Zidx_repetitions_allowed(VECTOR_OP V, VECTOR_OP V_Zidx, INT type, void *data)
{
	INT i, k, len, idx, idx1, idx2, f_found, len1;

	len = V->s_li();
	V_Zidx->m_il(len);
	len1 = 0;
	for (i = 0; i < len; i++) {
		s_PP()->search(s_PP()->s_li(), TRUE, V->s_i(i), &idx, &f_found);
		if (!f_found) {
			do_print(V->s_i(i), type, data);
			return error("SGL::Z2Zidx(): element not found in PP");
			}
		idx--;
		idx1 = s_PPidx_ii(idx);
		V_Zidx->search(V_Zidx->s_li(), TRUE, s_PPidx_i(idx), &idx2, &f_found);
		if (f_found) {
			continue;
			}
		if (V_Zidx->search_and_insert_int(len1, idx1) != OK) {
			printf("V = ");
			V->println();
			fflush(stdout);
			printf("PP = ");
			s_PP()->println();
			fflush(stdout);
			printf("PPidx = ");
			s_PPidx()->println();
			fflush(stdout);
			printf("V_Zidx = ");
			V_Zidx->println();
			printf("SGL::Z2Zidx(): duplicate zuppo index\n");
			fflush(stdout);
			return ERROR;
			/* return error("SGL::Z2Zidx(): duplicate zuppo index"); */
			}
		len1++;
		}
	V_Zidx->realloc_z(len1);
	return OK;
}

INT SGL_OB::Z2Zidx(VECTOR_OP V, VECTOR_OP V_Zidx, INT type, void *data)
/* before: vec_to_Zidx_vec() */
/* Vektor von Gruppenelementen (Zuppos) 
 * wird umgerechnet 
 * in Vektor ueber integer, bestehend aus 
 * den zugehoerigen Zuppo - Indices 
 * (aufsteigend sortiert); 
 * zugehoerig heisst, dass 
 * das jeweilige Element 
 * in der zyklischen Untergruppe 
 * des angegebenen Zuppos liegt, 
 * es wird also ueber 
 * PP und PPidx gegangen. */
{
	INT i, k, len, idx, idx1, f_found;

	len = V->s_li();
	V_Zidx->m_il(len);
	for (i = 0; i < len; i++) {
		s_PP()->search(s_PP()->s_li(), TRUE, V->s_i(i), &idx, &f_found);
		if (!f_found) {
			do_print(V->s_i(i), type, data);
			return error("SGL::Z2Zidx(): element not found in PP");
			}
		idx--;
		if (V_Zidx->search_and_insert_int(i, s_PPidx_ii(idx)) != OK) {
			printf("V = ");
			V->println();
			fflush(stdout);
			printf("PP = ");
			s_PP()->println();
			fflush(stdout);
			printf("PPidx = ");
			s_PPidx()->println();
			fflush(stdout);
			printf("V_Zidx = ");
			V_Zidx->println();
			printf("SGL::Z2Zidx(): duplicate zuppo index\n");
			fflush(stdout);
			return ERROR;
			/* return error("SGL::Z2Zidx(): duplicate zuppo index"); */
			}
		}
	return OK;
}


INT SGL_OB::Zidx2Z(VECTOR_OP V_Zidx, VECTOR_OP V)
/* before Zidx_vec_to_vec() */
{
	SYM_OP z;
	INT i, j, len, idx, f_found;

	len = V_Zidx->s_li();
	V->m_il(len);
	for (i = 0; i < len; i++) {
		z = s_Z_i(V_Zidx->s_ii(i));
		if (V->search_and_insert(i, z) != OK)
			return error("SGL::Zidx2Z(): duplicate zuppo indices");
		}
	return OK;
}

INT SGL_OB::IndexInNormalizer(INT layer, INT orbit, INT *idx_in_normalizer)
{
	SGO_OP U_orbit, N_orbit;
	INT l, Nlayer, Norbit;
	
	U_orbit = s_theOrbits_ij(layer, orbit);
	Nlayer = U_orbit->s_Nlayer_i();
	Norbit = U_orbit->s_Norbit_i();
	N_orbit = s_theOrbits_ij(Nlayer, Norbit);
	if (U_orbit->s_go_i() == 0)
		return error("SGL::IndexInNormalizer()|go == 0");
	l = N_orbit->s_go_i() / U_orbit->s_go_i();
	if (l * U_orbit->s_go_i() != N_orbit->s_go_i())
		return error("SGL::IndexInNormalizer()|no integer division !");
	*idx_in_normalizer = l;
	return OK;
}

INT SGL_OB::orbit2lo(INT orbit, INT *l, INT *o)
{
	INT l1, total_orbits;
	
	total_orbits = 0;
	for (l1 = 0; l1 <= s_nb_Layers_i(); l1++) {
		if (total_orbits + s_theOrbits_i(l1)->s_li() > orbit) {
			*l = l1;
			*o = orbit - total_orbits;
			return OK;
			}
		total_orbits += s_theOrbits_i(l1)->s_li();
		}
	return error("SGL::orbit2lo(): orbit number illegal");
}

INT SGL_OB::lo2orbit(INT l, INT o, INT *orbit)
{
	INT l1, orbit1;
	
	orbit1 = 0;
	for (l1 = 0; l1 < l; l1++) {
		orbit1 += s_theOrbits_i(l1)->s_li();
		}
	orbit1 += o;
	*orbit = orbit1;
	return OK;
}

INT SGL_OB::IsSubgroup(
	INT layer1, INT orbit1, INT rep1, 
	INT layer2, INT orbit2, INT rep2, 
	INT *f_is_subgroup)
{
	SGO_OP H_orbit, K_orbit;
	VECTOR_OP H_Zidx, K_Zidx;
	
	H_orbit = s_theOrbits_ij(layer1, orbit1);
	K_orbit = s_theOrbits_ij(layer2, orbit2);
	if (H_orbit->s_Zidx()->s_obj_k() == EMPTY)
		return error("SGL::IsSubgroup() H_orbit->s_Zidx() is EMPTY");
	H_Zidx = H_orbit->s_Zidx_i(rep1);
	if (K_orbit->s_Zidx()->s_obj_k() == EMPTY)
		return error("SGL::IsSubgroup() K_orbit->s_Zidx() is EMPTY");
	K_Zidx = K_orbit->s_Zidx_i(rep2);
	H_Zidx->subseteq(K_Zidx, H_Zidx->s_li(), K_Zidx->s_li(), TRUE, f_is_subgroup);
	return(OK);
}

INT SGL_OB::IsSubgroup_with_recalc(
	INT layer1, INT orbit1, INT rep1, 
	INT layer2, INT orbit2, INT rep2, 
	INT *f_is_subgroup, INT type, void *data)
{
	SGO_OP H_orbit, K_orbit;
	VECTOR_OP H_Zidx, K_Zidx;
	VECTOR_OB H_Zidx_, K_Zidx_;
	
	H_orbit = s_theOrbits_ij(layer1, orbit1);
	K_orbit = s_theOrbits_ij(layer2, orbit2);
	if (H_orbit->s_Zidx()->s_obj_k() != EMPTY) {
		H_Zidx = H_orbit->s_Zidx_i(rep1);
		}
	else {
		H_orbit->calc_Zidx(this, rep1, &H_Zidx_, type, data);
		H_Zidx = &H_Zidx_;
		}
	if (K_orbit->s_Zidx()->s_obj_k() != EMPTY) {
		K_Zidx = K_orbit->s_Zidx_i(rep2);
		}
	else {
		K_orbit->calc_Zidx(this, rep2, &K_Zidx_, type, data);
		K_Zidx = &K_Zidx_;
		}
	H_Zidx->subseteq(K_Zidx, H_Zidx->s_li(), K_Zidx->s_li(), TRUE, f_is_subgroup);
	return(OK);
}

INT SGL_OB::Asup(MATRIX_OP A)
{
	INT nb_orbits;
	SGO_OP H_orbit, K_orbit;
	INT i, j, l1, o1, l2, o2, i_len;
	INT k, a_ij, f_Leq;
	
	nb_orbits = s_total_nb_orbits_i();
	if (A->m_ilih(nb_orbits, nb_orbits) != OK)
		return error("SGL::Asup(): error in m_ilih_m()");
	for (i = 0; i < nb_orbits; i++) {
		orbit2lo(i, &l1, &o1);
		H_orbit = s_theOrbits_ij(l1, o1);
		/* i_len = H_orbit->s_Zidx()->s_li(); */
		i_len = H_orbit->s_o_len_i();
		for (j = 0; j < nb_orbits; j++) {
			if (j == i) {
				A->m_iji(i, j, 1);
				continue;
				}
			orbit2lo(j, &l2, &o2);
			K_orbit = s_theOrbits_ij(l2, o2);
			/* Berechne a_ij sup: */
			a_ij = 0;
			for (k = 0; k < i_len; k++) {
				if (IsSubgroup(l1, o1, k, l2, o2, 0, &f_Leq) != OK)
					return error("SGL::Asup(): error in IsSubgroup()");
				if (f_Leq) {
					a_ij++;
					}
				}
			A->m_iji(i, j, a_ij);
			}
		}
	return OK;
}

INT SGL_OB::Ni(MATRIX_OP D)
{
	INT nb_orbits, i, j, l1, o1, d_ij;
	
	nb_orbits = s_total_nb_orbits_i();
	if (D->m_ilih_n(nb_orbits, nb_orbits) != OK)
		return error("SGL::Ni(): error in m_ilih_m()");
	for (i = 0; i < nb_orbits; i++) {
		orbit2lo(i, &l1, &o1);
		for (j = 0; j < nb_orbits; j++) {
			if (j != i) {
				continue;
				}
			if (IndexInNormalizer(l1, o1, &d_ij) != OK)
				return error("SGL::Ni(): error in sgl_IndexInNormalizer()");
			D->m_iji(i, j, d_ij);
			}
		}
	return OK;
}

INT SGL_OB::OrbitSizes(VECTOR_OP v)
{
	INT nb_orbits, i, l1, o1, s;
	SGO_OP Orbit;
	
	nb_orbits = s_total_nb_orbits_i();
	if (v->m_il_n(nb_orbits) != OK)
		return error("SGL::OrbitSizes(): error in m_il_nv()");
	for (i = 0; i < nb_orbits; i++) {
		orbit2lo(i, &l1, &o1);
		Orbit = s_theOrbits_ij(l1, o1);
		/* s = Orbit->s_Zidx()->s_li(); */
		s = Orbit->s_o_len_i();
		v->m_ii(i, s);
		}
	return OK;
}

INT SGL_OB::shrink()
{
	SGO_OP Orbit;
	INT nb_orbits, i, l, o;

	nb_orbits = s_total_nb_orbits_i();
	for (i = 0; i < nb_orbits; i++) {
		orbit2lo(i, &l, &o);
		Orbit = s_theOrbits_ij(l, o);
		Orbit->s_Zidx()->freeself();
		}
	return OK;
}

INT SGL_OB::grow(INT type, void *data, INT f_v)
{
	SGO_OP Orbit;
	INT nb_orbits, i, l, o;

	nb_orbits = s_total_nb_orbits_i();
	for (i = 0; i < nb_orbits; i++) {
		orbit2lo(i, &l, &o);
		Orbit = s_theOrbits_ij(l, o);
		Orbit->recalc_Zidx(this, type, data, f_v);
		}
	return OK;
}

INT SGL_OB::Burnside_info(VECTOR_OP orbit_size, 
	MATRIX_OP As, MATRIX_OP Ai, MATRIX_OP D, 
	MATRIX_OP M, MATRIX_OP B, SYM_OP d, INT f_v)
{
	OrbitSizes(orbit_size);
	Asup(As);
	if (f_v) {
		printf("Asup =\n");
		As->Print();
		}
	As->Asup2Ainf(Ai);
	if (f_v) {
		printf("Ainf =\n");
		Ai->Print();
		}
	Ni(D);
	if (f_v) {
		printf("D =\n");
		D->Print();
		}
	Ai->mult(D, M);
	if (f_v) {
		printf("M =\n");
		M->Print();
		}
	M->invers(B);
	D->s_ij(0, 0)->copy(d);
	// M->determinante(d);
	B->mult_apply_scalar_matrix(d);
	if (f_v) {
		printf("B = 1/");
		d->println();
		B->Print();
		}
	return OK;
}

#endif /* SGL_TRUE */

