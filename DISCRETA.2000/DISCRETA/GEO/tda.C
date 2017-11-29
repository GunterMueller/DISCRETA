/* tda.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */



#include <DISCRETA/discreta.h>
#include <DISCRETA/ma.h> // for multivalued_geo
#include <DISCRETA/lb.h> // for TDA 
#include <DISCRETA/fga.h> // for  report_group_latex_stdout()

#include <DISCRETA/geo.h>




#define F_GENERATORS_BOTTOM_UP TRUE

void tda_geo_tex_head(FILE *fp)
{
	fprintf(fp, "\\documentclass[11pt,a4paper]{report}\n");
	fprintf(fp, "\\usepackage{rotating}\n");
	fprintf(fp, "\\newfont{\\mytt}{cmtt12 scaled 1000}\n");
	fprintf(fp, "\\textheight=640pt\n");
	fprintf(fp, "\\textwidth=440pt\n");
	fprintf(fp, "\\topmargin=0pt\n");
	fprintf(fp, "\\headsep=0pt\n");
	fprintf(fp, "\\footskip=10pt\n");
	fprintf(fp, "\\mathsurround=1pt\n");
	fprintf(fp, "\\evensidemargin=15pt\n");
	fprintf(fp, "\\oddsidemargin=15pt\n");

	fprintf(fp, "\n\\pagestyle{empty}\n\n");

	fprintf(fp, "\\begin{document}\n\n");
}

void tda_geo_tex_foot(FILE *fp)
{
	fprintf(fp, "\\end{document}\n\n");
}

INT calc_TDO(FILE *fp, INT nrow, INT ncol, INT nb_X, INT *X, 
	INT f_calc_second_tdo, TDOSS **tdoss)
{
	printf("tda.C : calc_TDO() disactivated!\n");
	exit(1);
#if 0
	SHORT Y[GEO_MAX_N << 1];
	SHORT *theY;
	TDOSS *tdo;
	INT i;
	
	tdo = calc_ntdo(nrow, ncol, nb_X, X, f_calc_second_tdo, 
		TRUE /* f_calc_theY */, Y, 
		FALSE /* gd->f_calc_ddp */, FALSE /* gd->f_calc_ddb */, 
		FALSE /* gd->f_tdo_v */, FALSE /* gd->f_tdo_vv */,
		FALSE /* gd->f_dd_v */, FALSE /* gd->f_dd_vv */);
	theY = the_Y_get_the_X(Y);
	for (i = 0; i < nb_X; i++) {
		X[i] = theY[i];
		}
	*tdoss = tdo;
	// tdoss_free(tdo);
	ntdo_the_Y_free(Y); /* eventually free ddp, ddb ! */
#endif
}

INT calc_TDA(FILE *TDA_fp, FILE *TEX_fp, INT nb_geo, 
	INT nrow, INT ncol, INT nb_X, INT *X, 
	INT f_transposed, LABRA_OP A, 
	INT *nb_row_blocks_TDO, INT *nb_col_blocks_TDO, 
	INT *nb_row_blocks_TDA, INT *nb_col_blocks_TDA, 
	INT f_v, INT f_vv, INT f_incidences)
{
	VECTOR_OB P_tda_ordering;
	VECTOR_OB Q_tda_ordering;
	VECTOR_OB P_orbit_first, P_orbit_length;
	VECTOR_OB Q_orbit_first, Q_orbit_length;
	VECTOR_OB P_O_labra;
	VECTOR_OB Q_O_labra;
	VECTOR_OB row_char_first, col_char_first;
	PERMUTATION_OB p, pv, q, qv;
		// p takes the points into tda-ordering, pv = p^{-1}
		// q the blocks
	PERMUTATION_OB qlex, qlexv;
	PERMUTATION_OB tmp1, tmp2;
	INT i, j, l, first, length, v, b;
	VECTOR_OB B, gen, gen0, gen1, gen2;
	LABRA_OB A1;
	SYM_OB ago, go;
	TDOSS *tdoss;
	VECTOR_OP row_inv_first, col_inv_first;
	PERMUTATION_OP pp, qq;
	PERMUTATION_OB ppp, qqq;
	
	calc_TDO(TDA_fp, nrow, ncol, nb_X, X, 
		TRUE /* f_calc_second_tdo */, &tdoss);
	tdoss_2_first(tdoss, &row_char_first, &col_char_first);
	
	geo_Canonicize(FALSE /* f_maxtest */, NIL /* back_to */, 
		nrow, ncol, nb_X, 
		FALSE /* f_print_backtrack_points */, 
		X, tdoss, f_transposed, 
		FALSE /* gd->f_calc_ddp */, NIL /* ddp */, 
		FALSE /* gd->f_calc_ddb */, NIL /* ddb */, 
		&ppp, &qqq, 
		TRUE /* f_get_aut_group */, A, &ago, 
		FALSE /* gd->f_canon_v */, FALSE /* gd->f_canon_vv */);
	tdoss_free(tdoss);
	
	// A->generators(&gen, &go);
	A->strong_generating_set(&gen, TRUE /* f_v */);
	// A->reduced_generating_set(&gen, F_GENERATORS_BOTTOM_UP, TRUE /* f_v */);
	if (f_v) {
		// report_group_latex_stdout(&gen);
		}
	A->group_order(&go);
	if (f_v) {
		printf("group order = ");
		go.println();
		}
	
	A->orbits(&P_tda_ordering, &P_orbit_first, &P_orbit_length, 
		FALSE /* f_calc_labras */, &P_O_labra, FALSE /* f_v */);
	if (f_transposed) {
		if (nrow != P_tda_ordering.s_li())
			return error("calc_TDA() f_transposed && nrow != P_tda_ordering.s_li()");
		v = nrow;
		b = ncol;
		}
	else {
		if (ncol != P_tda_ordering.s_li())
			return error("calc_TDA() !f_transposed && ncol != P_tda_ordering.s_li()");
		v = ncol;
		b = nrow;
		}
	
	p.m_il(v);
	for (i = 0; i < v; i++) {
		j = P_tda_ordering.s_ii(i);
		p.m_ii(j, i + 1);
		}
	p.invers(&pv);
	if (f_vv) {
		printf("tda-ordering of the points:\n");
		P_tda_ordering.println();
		}
	l = P_orbit_first.s_li();
	for (i = 0; i < l; i++) {
		first = P_orbit_first.s_ii(i);
		length = P_orbit_length.s_ii(i);
		if (f_vv) {
			printf("%ld/%ld\n", first, length);
			}
		}

	Q_orbit_first.m_il(1);
	Q_orbit_length.m_il(1);
	Q_orbit_first.m_ii(0, 0);
	// Q_orbit_len.m_ii(0, b);

#if 0
	if (TDA_fp != stdout) {
		TD_print(TDA_fp, nrow, ncol, nb_X, X, f_transposed, 
			&p, NIL, &P_orbit_first, &Q_orbit_first);
		}
#endif
	if (f_vv) {
		printf("point orbits: \n");
		TD_print(stdout, nrow, ncol, nb_X, X, f_transposed, 
			&p, NIL, &P_orbit_first, &Q_orbit_first);
		}
	
	if (f_vv) {
		printf("calculating block orbits\n");
		}
#if 0
	A->stab_generators(0, &gen, FALSE /* f_v */);
	A->generators(&gen, &go);
#endif
	
	if (f_vv) {
		printf("permutation of the points; tda_ordering =\n");
		P_tda_ordering.println();
		printf("p=");
		p.println();
		}
	calc_blocks(nrow, ncol, nb_X, X, f_transposed, 
		&p, &pv, &qlex, &qlexv, &B, f_vv);

	l = gen.s_li();
	gen0.m_il(l);
	for (i = 0; i < l; i++) {
		pv.mult((PERMUTATION_OP) gen.s_i(i), &tmp1);
		tmp1.mult(&p, &tmp2);
		tmp2.copy((PERMUTATION_OP) gen0.s_i(i));
		}
	gen0.copy(&gen1);
	vec_induce_action_on_blocks(&gen1, &B);
	l = gen1.s_li();
	gen2.m_il(l);
	for (i = 0; i < l; i++) {
		qlex.mult((PERMUTATION_OP) gen1.s_i(i), &tmp1);
		tmp1.mult(&qlexv, &tmp2);
		tmp2.copy((PERMUTATION_OP) gen2.s_i(i));
		}
	if (f_vv) {
		l = gen.s_li();
		for (i = 0; i < l; i++) {
			printf("generator %ld on points: ", i);
			gen.s_i(i)->println();
			printf("on TDA ordering of points: ");
			gen0.s_i(i)->println();
			printf("on lexicographic blocks: ");
			gen1.s_i(i)->println();
			printf("on blocks: ");
			gen2.s_i(i)->println();
			}
		}
	A1.init_by_generators(&gen2, b, FALSE);
	if (f_v) {
		printf("the group induced on the blocks ");
		A1.print_group_order();
		}
	A1.orbits(&Q_tda_ordering, &Q_orbit_first, &Q_orbit_length, 
		FALSE /* f_calc_labras */, &Q_O_labra, FALSE /* f_v */);
	
	if (b != Q_tda_ordering.s_li())
		return error("calc_TDA() b != Q_tda_ordering.s_li()");
	q.m_il(b);
	for (i = 0; i < b; i++) {
		j = Q_tda_ordering.s_ii(i);
		q.m_ii(j, i + 1);
		}
	q.invers(&qv);
	if (f_vv) {
		printf("tda-ordering of the blocks:\n");
		Q_tda_ordering.println();
		}
	l = Q_orbit_first.s_li();
	for (i = 0; i < l; i++) {
		first = Q_orbit_first.s_ii(i);
		length = Q_orbit_length.s_ii(i);
		if (f_vv) {
			printf("%ld/%ld\n", first, length);
			}
		}

#if 0
	if (TDA_fp != stdout) {
		TD_print(TDA_fp, nrow, ncol, nb_X, X, f_transposed, 
			&p, &q, &P_orbit_first, &Q_orbit_first);
		}
	if (f_v) {
		TD_print(stdout, nrow, ncol, nb_X, X, f_transposed, 
			&p, &q, &P_orbit_first, &Q_orbit_first);
		}
#endif

	if (f_transposed) {
		row_inv_first = &P_orbit_first;
		col_inv_first = &Q_orbit_first;
		pp = &p;
		qq = &q;
		}
	else {
		row_inv_first = &Q_orbit_first;
		col_inv_first = &P_orbit_first;
		pp = &q;
		qq = &p;
		}
	*nb_row_blocks_TDO = row_char_first.s_li();
	*nb_col_blocks_TDO = col_char_first.s_li();
	*nb_row_blocks_TDA = row_inv_first->s_li();
	*nb_col_blocks_TDA = col_inv_first->s_li();
	
	if (TDA_fp != stdout) {
		geo_output_char_inv_decomposition(TDA_fp, FALSE /* f_tex */, nb_geo, 
			nrow, ncol, nb_X, X, 
			pp, qq, 
			&row_char_first, &col_char_first, 
			row_inv_first, col_inv_first, 
			&ago, &gen0, f_incidences);
		}
	if (TEX_fp) {
		geo_output_char_inv_decomposition(TEX_fp, TRUE /* f_tex */, nb_geo, 
			nrow, ncol, nb_X, X, 
			pp, qq, 
			&row_char_first, &col_char_first, 
			row_inv_first, col_inv_first, 
			&ago, &gen0, f_incidences);
		}
	geo_output_char_inv_decomposition(stdout, FALSE /* f_tex */, nb_geo, 
		nrow, ncol, nb_X, X, 
		pp, qq, 
		&row_char_first, &col_char_first, 
		row_inv_first, col_inv_first, 
		&ago, &gen0, f_incidences);
	
	return OK;
}

INT geo_output_char_inv_decomposition(FILE *fp, INT f_tex, INT nb_geo, 
	INT nrow, INT ncol, INT nb_X, INT *X, 
	PERMUTATION_OP pp, PERMUTATION_OP qq, 
	VECTOR_OP row_char_first, 
	VECTOR_OP col_char_first, 
	VECTOR_OP row_inv_first, 
	VECTOR_OP col_inv_first, 
	SYM_OP ago, VECTOR_OP generators, INT f_incidences)
{
	INT i, l;
	BYTE str[10000];
	
	if (f_tex) {
		fprintf(fp, "\\noindent\n");
		fprintf(fp, "\\parbox[t]{16cm}{\n");
		fprintf(fp, "geo no %ld:\\\\\n", nb_geo);
		fprintf(fp, "{\\mytt\n", nb_geo);
		
		TD_char_inv_print(fp, TRUE /* f_tex */, nrow, ncol, nb_X, X, pp, qq, 
			row_char_first, col_char_first, 
			row_inv_first, col_inv_first);
		fprintf(fp, "}\n%% end mytt\n");
		fprintf(fp, "}\\\\\n");
		fprintf(fp, "${\\rm Aut} = \\langle \n");
		l = generators->s_li();
		for (i = 0; i < l; i++) {
			str[0] = 0;
			((PERMUTATION_OP) generators->s_i(i))->sprint_latex(str);
			fprintf(fp, "$%s$", str);
			if (i < l - 1) {
				fprintf(fp, ", ");
				fprintf(fp, "\\newline");
				}
			}
		fprintf(fp, "\\rangle$\\newline");
		fprintf(fp, "$|{\\rm Aut}| = ");
		str[0] = 0;
		ago->sprint(str);
		fprintf(fp, "%s$\\newline\n", str);
		if (f_incidences) {
			for (i = 0; i < nb_X; i++) {
				fprintf(fp, "%ld ", X[i]);
				}
			fprintf(fp, "\\newline\n");
			}
		}
	else {
		fprintf(fp, "geo no %ld:\n", nb_geo);
		
		TD_char_inv_print(fp, FALSE /* f_tex */, nrow, ncol, nb_X, X, pp, qq, 
			row_char_first, col_char_first, 
			row_inv_first, col_inv_first);
		fprintf(fp, "automorphism group order = ");
		ago->fprintln(fp);
		fprintf(fp, "generators:\n");
		l = generators->s_li();
		for (i = 0; i < l; i++) {
			generators->s_i(i)->fprintln(fp);
			}
		if (f_incidences) {
			for (i = 0; i < nb_X; i++) {
				fprintf(fp, "%ld ", X[i]);
				}
			fprintf(fp, " %%GEODATA %ld\n", nb_geo);
			}
	
		l = row_char_first->s_li();
		fprintf(fp, "%ld ", l);
		for (i = 0; i < l; i++) {
			fprintf(fp, "%ld ", row_char_first->s_ii(i));
			}
		fprintf(fp, " %%GEOrow_char\n");
		l = col_char_first->s_li();
		fprintf(fp, "%ld ", l);
		for (i = 0; i < l; i++) {
			fprintf(fp, "%ld ", col_char_first->s_ii(i));
			}
		fprintf(fp, " %%GEOcol_char\n");
		l = row_inv_first->s_li();
		fprintf(fp, "%ld ", l);
		for (i = 0; i < l; i++) {
			fprintf(fp, "%ld ", row_inv_first->s_ii(i));
			}
		fprintf(fp, " %%GEOrow_inv\n");
		l = col_inv_first->s_li();
		fprintf(fp, "%ld ", l);
		for (i = 0; i < l; i++) {
			fprintf(fp, "%ld ", col_inv_first->s_ii(i));
			}
		fprintf(fp, " %%GEOcol_inv\n");
		fprintf(fp, "\n");
		}
		
	fflush(fp);
	return OK;
}

INT calc_blocks(INT nrow, INT ncol, INT nb_X, INT *X, INT f_transposed, 
	PERMUTATION_OP p, PERMUTATION_OP pv, 
	PERMUTATION_OP q, PERMUTATION_OP qv, 
	VECTOR_OP B, INT f_v)
{
	INT i, j, k, a, v, b;
	INT nb_empty, idx, f_found;
	VECTOR_OP block;
	VECTOR_OB B1;
	INT *Xi = NIL;
	INT *Xj = NIL;

	Xi = (INT *) my_malloc(nb_X * sizeof(INT), "calc_blocks(): Xi");
	Xj = (INT *) my_malloc(nb_X * sizeof(INT), "calc_blocks(): Xi");
	if (Xi == NIL || Xj == NIL)
		return error("calc_blocks() no memory");
	for (k = 0; k < nb_X; k++) {
		a = X[k];
		i = a / ncol;
		j = a % ncol;
		if (f_transposed) {
			Xi[k] = p->s_ii(i) - 1;
			Xj[k] = j;
			}
		else {
			Xi[k] = p->s_ii(j) - 1;
			Xj[k] = i;
			}
		}
	if (f_transposed) {
		v = nrow;
		b = ncol;
		}
	else {
		v = ncol;
		b = nrow;
		}
	B->m_il(b);
	for (j = 0; j < b; j++) {
		block = (VECTOR_OP) B->s_i(j);
		block->m_il(0);
		}
	for (k = 0; k < nb_X; k++) {
		i = Xi[k];
		j = Xj[k];
		block = (VECTOR_OP) B->s_i(j);
		block->inc();
		block->m_ii(block->s_li() - 1, i);
		}
	for (j = 0; j < b; j++) {
		block = (VECTOR_OP) B->s_i(j);
		block->quicksort(block->s_li(), TRUE /* f_ascending */);
		}
	B->copy(&B1);
	B->quicksort(B->s_li(), TRUE /* f_ascending */);
	if (f_v) {
		printf("the original blocks:\n");
		for (i = 0; i < b; i++) {
			block = (VECTOR_OP) B1.s_i(i);
			printf("block %ld: ", i);
			block->println();
			}
		printf("the lexicographic blocks:\n");
		for (i = 0; i < b; i++) {
			block = (VECTOR_OP) B->s_i(i);
			printf("block %ld: ", i);
			block->println();
			}
		}
	q->m_il(b);
	nb_empty = 0;
	for (i = 0; i < b; i++) {
		block = (VECTOR_OP) B1.s_i(i);
		if (block->s_li() == 0) {
			idx = nb_empty;
			nb_empty++;
			}
		else {
			B->search(b, TRUE, block, &idx, &f_found);
			if (!f_found)
				return error("calc_blocks(): not found");
			idx--;
			}
		q->m_ii(i, idx + 1);
		if (f_v) {
			printf("original block %ld is lexicographic block %ld\n", i, idx);
			}
		}
	q->invers(qv);
	if (f_v) {
		printf("found %ld blocks:\n", b);
		for (j = 0; j < b; j++) {
			block = (VECTOR_OP) B->s_i(j);
			printf("%ld (original block no %ld): ", j, qv->s_ii(j) - 1);
			block->println();
			}
		}
	my_free(Xi);
	my_free(Xj);
	return OK;
}


