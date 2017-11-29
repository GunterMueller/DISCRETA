/* sgl.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SGL_TRUE

#include <DISCRETA/divs.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/fga.h>
#include <DISCRETA/perm.h>
#include <DISCRETA/dimino.h>
#include <DISCRETA/sgl.h>

INT calc_sgl(VECTOR_OP G_gen, VECTOR_OP G, VECTOR_OP Inn_gen, SGL_OP L)
{
	VECTOR_OB H_gen, H;
	VECTOR_OB N_gen, N;
	
	VECTOR_OB G1_gen, G1;
	INT t0, t1, user_time;
	BYTE str[256];
	VECTOR_OB Aut_gen;
		
	t0 = os_ticks();


	gruppen_elemente(G_gen, G, 0, NIL);
	
		
	printf("calling get_Inn_G()\n");
	fflush(stdout);
	get_Inn_G(G_gen, G, &Aut_gen);
	Aut_gen.copy(Inn_gen);
	printf("found %ld generators of the (inner) automorphism group\n", 
		Aut_gen.s_li());
	printf("calling SGL::Init_with_Aut_generators()\n");
	fflush(stdout);
	L->Init_with_Aut_generators(
		G_gen, G, &Aut_gen, 
		TRUE /* f_verbose */, 
		FALSE /* f_very_verbose */, 0, NIL);
	printf("finished SGL::Init_with_Aut_generators()\n");
	fflush(stdout);

#if 0
	if (f_Z_info) {
		printf("calling L->latex_Z_info()\n");
		fflush(stdout);
		L->latex_Z_info(z_tex, tex_group_name, 0, NIL);
		}
#endif
	printf("calling L->all_layers()\n");
	fflush(stdout);
	L->all_layers(TRUE /* f_verbose */, 0, NIL);
#if 0
	L->PrintOrbits(TRUE /* f_zuppos_expanded */, 0, NIL);
#endif
	L->PrintStatistics();

	t1 = os_ticks();
	user_time = t1 - t0;
	strcpy(str, "Running time: ");
	print_delta_time(user_time, str);
	printf("%s\n", str);
	

	return OK;
}

INT get_Inn_G(VECTOR_OP G_gen, VECTOR_OP G, VECTOR_OP Aut_gen)
{
	PERMUTATION_OB p;
	SYM_OB gen, gen_v, a, b;
	INT i, j, l, n, idx, f_found;
	
	l = G_gen->s_li();
	n = G->s_li();
	
	Aut_gen->m_il(l);
	for (i = 0; i < l; i++) {
		p.m_il(n);
		G_gen->s_i(i)->copy(&gen);
		gen.invers(&gen_v);
		for (j = 0; j < n; j++) {
			gen_v.mult(G->s_i(j), &a);
			a.mult(&gen, &b);
			G->search(n, TRUE /* f_ascending */, &b, &idx, &f_found);
			if (!f_found)
				return error("get_Inn_G(): not found");
			idx--;
			p.m_ii(j, idx + 1);
			}
		// p.println();
		// fflush(stdout);
		p.swap((PERMUTATION_OP) Aut_gen->s_i(i));
		}
	return OK;
}

SG_LATTICE_LOCAL *open_sgll(SGL_OP L, 
	INT f_L_allocated, 
	INT f_show_generators, 
	INT f_with_perm, INT draw_lines_type, 
	INT f_draw_sgl, INT f_burnside_tex, 
	BYTE *path_name, BYTE *group_name, 
	BYTE *tex_group_name, INT f_verbose)
/* draw_lines_type: 0 = Asup, 1 = Ainf, 2 = all, 3 = poset of orbits */
{
	MATRIX_OP Asup = (MATRIX_OP) callocobject("open_sgll");
	MATRIX_OP Ainf = (MATRIX_OP) callocobject("open_sgll");
	MATRIX_OP D = (MATRIX_OP) callocobject("open_sgll");
	MATRIX_OP M = (MATRIX_OP) callocobject("open_sgll");
	MATRIX_OP B = (MATRIX_OP) callocobject("open_sgll");
	MATRIX_OP Acover = (MATRIX_OP) callocobject("open_sgll");
	VECTOR_OP nl = (VECTOR_OP) callocobject("open_sgll");
	VECTOR_OP orbit_size = (VECTOR_OP) callocobject("open_sgll");
	SYM_OP d = callocobject("open_sgll");
	MEM_OP vbp_plaz = (MEM_OP) callocobject("open_sgll");
	MEM_OP vbp_o_dx = (MEM_OP) callocobject("open_sgll");
	SG_LATTICE_LOCAL *p = NIL;
	BYTE tex_fname[1024];
	FILE *lp;
	
	
	sprintf(tex_fname, "%s%s.tex", path_name, group_name);
	
	p = (SG_LATTICE_LOCAL *) my_malloc(sizeof(SG_LATTICE_LOCAL), "open_sgll");
	if (p == NIL)
		return NIL;
	L->Burnside_info(orbit_size, Asup, Ainf, D, M, B, d, f_verbose);

	p->plaz = vbp_plaz;
	p->o_dx = vbp_o_dx;
	p->L = L;
	p->f_L_allocated = f_L_allocated;
	p->Asup = Asup;
	p->nl = nl;
	p->orbit_size = orbit_size;
	p->d = d;
	p->subgroups = NIL;
	p->generators = NIL;
	sgll_init_generators(p);
	Asup->Asup2Acover(Acover);
	if (f_verbose) {
		printf("Acover =\n");
		Acover->Print();
		}
	Acover->Acover2nl(nl);
	if (f_verbose) {
		printf("nl =\n");
		nl->println();
		}
	vbp(nl, orbit_size, 
		SGL_VBP_X_PIX, SGL_VBP_Y_PIX, 
		vbp_plaz, vbp_o_dx, TRUE /* f_upside_down */);
	
	p->extrema[0] = 0.;
	p->extrema[1] = (double) SGL_VBP_X_PIX /* * 2. */;
	p->extrema[2] = 0.;
	p->extrema[3] = (double) SGL_VBP_Y_PIX /* * 2. */;
	p->extrema[4] = -1.;
	p->extrema[5] = 1.;

	/* enlarge the CO-system to shorten the graphics output ! */

	p->f_show_generators = f_show_generators;
	p->f_with_perm = f_with_perm;
	p->draw_lines_type = draw_lines_type;
		/* 0 = Asup, 1 = Ainf, 2 = all, 3 = Poset of orbits */
	strcpy(p->path_name, path_name);
	strcpy(p->group_name, group_name);
	strcpy(p->tex_group_name, tex_group_name);
	if (f_draw_sgl) {
		BYTE caption_sgl[1000];
		BYTE fname[1000];

		sprintf(fname, "%s%s.mp", path_name, group_name);
		sprintf(caption_sgl, "$%s$", p->tex_group_name);
		draw_mp_file(p, &p->vdev, fname, 1000 /* factor 1000 */, sgll_draw_func);
		// draw_PS_file(p, &p->vdev, SGL_VBP_X_PIX, SGL_VBP_Y_PIX, ps_name, sgll_draw_func);
		// draw_epic_file(p, &p->vdev, epic_name, 666 /* factor_1000 */, caption_sgl /* caption */, sgll_draw_func);
		}
	if (f_burnside_tex) {
		lp = fopen(tex_fname, "w");
		fprintf(lp, "\\documentclass{article}\n");
		fprintf(lp, "\\usepackage{epsfig}\n");
		fprintf(lp, "\\begin{document}\n");

		fprintf(lp, "\\begin{center}\n");
		fprintf(lp, "\\epsfxsize=150mm\n");
		fprintf(lp, "$\\epsffile{%s%s.1}$\n", path_name, group_name);
		fprintf(lp, "\\end{center}\n");
		
		fprintf(lp, "\\clearpage\n\n");

		sgll_output_layers(p, lp);
		
		fprintf(lp, "\\[\nA^\\vee \\; = \\; \n");
		Asup->latex_upper_tri(lp);
		fprintf(lp, "\\]\n\\par\n");
#ifdef BRUCHTRUE
		fprintf(lp, "\\[\nA^\\wedge \\; = \\; \n");
		Ainf->latex_upper_tri(lp);
		fprintf(lp, "\\]\n\\par\n");
		fprintf(lp, "\\[M \\; = \\; \n");
		M->latex_upper_tri(lp);
		fprintf(lp, "\\]\n\\par\n");
		fprintf(lp, "\\[\nB(%s) \\; = \\; \n", tex_group_name);
		fprintf(lp, "\\frac{1}{");
		d->latex(lp);
		fprintf(lp, "}\n");
		B->latex_upper_tri(lp);
		fprintf(lp, "\\]\n\\par\n");
#endif


		fprintf(lp, "\\end{document}\n\n");
		fclose(lp);
		}
	freeall(Ainf);
	freeall(Acover);
	freeall(D);
	freeall(M);
	freeall(B);
	return p;
}

INT free_sgll(SG_LATTICE_LOCAL *p)
{
	if (p == NIL)
		return OK;
	if (p->plaz) {
		freeall(p->plaz);
		p->plaz = NIL;
		}
	if (p->o_dx) {
		freeall(p->o_dx);
		p->o_dx = NIL;
		}
	if (p->f_L_allocated) {
		if (p->L) {
			freeall(p->L);
			p->L = NIL;
			}
		}
	if (p->Asup) {
		freeall(p->Asup);
		p->Asup = NIL;
		}
	if (p->nl) {
		freeall(p->nl);
		p->nl = NIL;
		}
	if (p->orbit_size) {
		freeall(p->orbit_size);
		p->orbit_size = NIL;
		}
	if (p->subgroups) {
		freeall(p->subgroups);
		p->subgroups = NIL;
		}
	if (p->generators) {
		freeall(p->generators);
		p->generators = NIL;
		}
	my_free(p);
	return OK;
}

INT sgll_output_layers(SG_LATTICE_LOCAL *sgll, FILE *fp)
{
	VECTOR_OP Gen;
	VECTOR_OP gen;
	SGL_OP L = sgll->L;
	INT i, j, jj, l, ll;
	PERMUTATION_OP p;
	INT nb_layers, nb_orbits, k;
	SGO_OP sgo;
	
	L = sgll->L;
	Gen = sgll->generators;
	l = L->s_total_nb_orbits_i();
	k = 0;
	nb_layers = L->s_theOrbits()->s_li();
	fprintf(fp, "\\begin{center} \n");
	fprintf(fp, "\\begin{tabular}{cccl} \n");
	fprintf(fp, "l/o & go & n & generators \\\\\n");
	fprintf(fp, "\\hline\n");
	fprintf(fp, "\\hline\n");
	
	for (i = 0; i < nb_layers; i++) {
		nb_orbits = L->s_theOrbits_i(i)->s_li();
		for (j = 0; j < nb_orbits; j++) {
			sgo = L->s_theOrbits_ij(i, j);
			gen = (VECTOR_OP) Gen->s_i(k);
			ll = gen->s_li();
			fprintf(fp, "%ld/%ld: & ", i, j);
			fprintf(fp, "%ld & %ld/%ld & ", 
				sgo->s_go_i(), sgo->s_Nlayer_i(), sgo->s_Norbit_i());
			fprintf(fp, "$\\langle ");
			for (jj = 0; jj < ll; jj++) {
				p = (PERMUTATION_OP) gen->s_i(jj);
				p->latex(fp);
				if (jj < ll - 1)
					fprintf(fp, "$, \\\\\n & & & $");
				}
			fprintf(fp, "\\rangle$ ");
			fprintf(fp, "\\\\ \n");
			fprintf(fp, "\\hline\n");
			k++;
			}
		fprintf(fp, "\\hline\n");
		}
	fprintf(fp, "\\end{tabular}\\\\\n");
	fprintf(fp, "l=layer, o=orbit, go=group order, n = normalizer layer/orbit\\\\\n");
	fprintf(fp, "\\end{center} \n");
	fflush(fp);
	return OK;
}

INT sgll_init_generators(SG_LATTICE_LOCAL *sgll)
{
	VECTOR_OP Gen = (VECTOR_OP) callocobject("sgll_init_generators");
	INT i, j, l;
	SGL_OP L = sgll->L;
	INT nb_layers, nb_orbits, k;
	SGO_OP sgo;

	sgll->generators = Gen;
	l = L->s_total_nb_orbits_i();
	Gen->m_il(l);
	k = 0;
	nb_layers = L->s_theOrbits()->s_li();
	for (i = 0; i < nb_layers; i++) {
		nb_orbits = L->s_theOrbits_i(i)->s_li();
		for (j = 0; j < nb_orbits; j++) {
			sgo = L->s_theOrbits_ij(i, j);
			sgo->generators(L, (VECTOR_OP) Gen->s_i(k), 0, NIL);
			k++;
			}
		}
	return OK;
}

#endif /* SGL_TRUE */

