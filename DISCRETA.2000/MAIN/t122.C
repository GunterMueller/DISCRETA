// t122.C
//
// 
// computes global intersection numbers of designs.
// 
//
// 
//
// Anton Betten
// Bayreuth, 16.05.1999

#include <DISCRETA/discreta.h>
#include <DISCRETA/ladder.h>

#define COLUMN_SPLIT  6000

INT do_it(BYTE *KM_fname, INT lambda, INT s, INT f_tex, INT f_v);

FILE *fp_tex;

int main(int argc, char **argv)
{
	INT t0, t1, user_time;
	BYTE s[256];

	if( argc < 3 ) {
		printf("usage: t122.out [options] km-file-name lambda s\n");
		printf("where options can be:\n");
		printf("-v         : verbose\n");
		printf("-tex       : produce tex output in the file\n");
		printf("           : report_<km-file-name>_<lambda>\n");
		return 1;
		}
	discreta_init();
	{
	// INT t0, t1, user_time;
	BYTE str[256];
	BYTE *km;
	INT lambda;
	INT i, s;
	INT f_v = FALSE;
	INT f_tex = FALSE;
	
	t0 = os_ticks();
	
	for (i = 1; i < argc - 3; i++) {
		if (strcmp(argv[i], "-v") == 0) {
			f_v = TRUE;
			}
		else if (strcmp(argv[i], "-tex") == 0) {
			f_tex = TRUE;
			}
		}
	km = argv[argc - 3];
	sscanf(argv[argc - 2], "%ld", &lambda);
	sscanf(argv[argc - 1], "%ld", &s);


	if (f_tex) {
		BYTE fname[10000];
		
		sprintf(fname, "report_%s_%ld.tex", km, lambda);
		fp_tex = fopen(fname, "w");
		fprintf(fp_tex, "%% tex report for %s with lambda = %ld\n", km, lambda);
		fprintf(fp_tex, "%% created by DISCRETA t122.out\n");
		fflush(fp_tex);
		}
	do_it(km, lambda, s, f_tex, f_v);
	if (f_tex) {
		fclose(fp_tex);
		}

	t1 = os_ticks();
	user_time = t1 - t0;
	str[0] = 0;
	print_delta_time_100(user_time, str);
	printf("total computing time: %s\n", str);
	fflush(stdout);
	
	}
	discreta_exit();
	return 0;
}

INT do_it(BYTE *KM_fname, INT lambda, INT s, INT f_tex, INT f_v)
{
	VECTOR_OB stab_go;
	INT i, j, ii, l, d, v, t, k, nb_sol, no;
	VECTOR_OB S;
	MATRIX_OB P;
	MATRIX_OP pP, P1, P2, P3;
	VECTOR_OP pS, pz;
	VECTOR_OB zz, z, z1, z2, vi, orbits;
	MATRIX_OB Bv, S1, S1t, L, Y, D, As1, As2;
	SYM_OB go, b, c;
	
	
	km_read_ascii_vtk(KM_fname, &v, &t, &k);
	nb_sol = km_nb_of_solutions(KM_fname, lambda);
	printf("found %ld solutions for lambda = %ld\n", nb_sol, lambda);

	if (f_tex) {
		fprintf(fp_tex, "\\section{%ld Designs with $\\lambda=%ld$}\n", nb_sol, lambda);
		}
	km_get_solutions(KM_fname, lambda, 0 /* from */, nb_sol /* len */, &S);
	P.m_ilih(k + 1, k + 1);
	{
		MATRIX_OB M;
		VECTOR_OB G_gen, RR, MM, Orbits_below1, Orbits_below2;
		
		km_read_ascii(KM_fname, &M, &v, &t, &k, 
			&G_gen, &RR, &MM, &stab_go, &Orbits_below1, &Orbits_below2);
		for (i = 0; i < k; i++) {
			MM.s_i(i)->swap(P.s_ij(i, i + 1));
			}
	}
	((VECTOR_OP) stab_go.s_i(0))->s_i(0)->copy(&go);
	printf("the group order is ");
	go.println();
	fflush(stdout);

	for (i = 0; i < k; i++) {
		pP = (MATRIX_OP) P.s_ij(i, i + 1);
		printf("matrix %ld,%ld of size %ld,%ld\n", i, i + 1, 
			pP->s_hi(), pP->s_li());
		}
	fflush(stdout);
#if 0
	for (d = 2; d <= k - t; d++) {
		for (i = t; i < k; i++) {
			j = i + d;
			if (j > k)
				continue;
			printf("d=%ld: i=%ld j=%ld\n", d, i, j); fflush(stdout);
			P1 = (MATRIX_OP) P.s_ij(i, j - 1);
			P2 = (MATRIX_OP) P.s_ij(j - 1, j);
			P3 = (MATRIX_OP) P.s_ij(i, j);
			dc_Mtk_via_Mtr_Mrk(i, j - 1, j, P1, P2, P3, f_v);
			}
		}
	for (i = t; i < k; i++) {
		for (j = i + 1; j <= k; j++) {
			P1 = (MATRIX_OP) P.s_ij(i, j);
			printf("matrix %ld,%ld of size %ld,%ld\n", i, j, 
				P1->s_hi(), P1->s_li());
			if (f_tex) {
				fprintf(fp_tex, "$P_{%ld,%ld}^\\cap$ of size $%ld \\times %ld$\\\\\n", 
					i, j, P1->s_hi(), P1->s_li());
				}
#if 1
			if (P1->s_li() > COLUMN_SPLIT) {
				INT ii;
				VECTOR_OB S;
				P1->split_column_wise(&S, COLUMN_SPLIT /* col_width */);
				for (ii = 0; ii < S.s_li(); ii++) {
					((MATRIX_OP) S.s_i(ii))->integer_print_dense();
					}
				}
			else {
				P1->integer_print_dense();
				}
#endif
			}
		}
#endif



	// compute orbit length matrices:
	for (i = 0; i <= k; i++) {
		P1 = (MATRIX_OP) P.s_ij(0, i);
		pS = (VECTOR_OP) stab_go.s_i(i);
		l = pS->s_li();
		P1->m_ilih(l, 1);
		for (j = 0; j < l; j++) {
			go.ganzdiv_integral(pS->s_i(j), P1->s_ij(0, j));
			}
		}
	printf("orbit length matrices computed\n");
	fflush(stdout);
	
#if 0
	{
	MATRIX_OB tmp;
	
	tmp.m_ilih(1, 1);
	tmp.m_iji(0, 0, 1);
	tmp.swap((MATRIX_OP) P.s_ij(0, 0));
	for (j = 2; j <= k; j++) {
		P1 = (MATRIX_OP) P.s_ij(0, j - 1);
		P2 = (MATRIX_OP) P.s_ij(j - 1, j);
		P3 = (MATRIX_OP) P.s_ij(0, j);
		dc_Mtk_via_Mtr_Mrk(0, j - 1, j, P1, P2, P3, f_v);
		}	
	}
#endif
	for (j = 0; j <= k; j++) {
		P.s_ij(0, j)->println();
		}
	if (f_tex) {
		fprintf(fp_tex, "orbit lengths:\n");
		for (j = 0; j <= k; j++) {
			fprintf(fp_tex, "%ld-orbits: ", j);
			pP = (MATRIX_OP) P.s_ij(0, j);
			l = pP->s_li();
			for (ii = 0; ii < l; ii++) {
				pP->s_ij(0, ii)->latex(fp_tex);
				if (ii < l - 1)
					fprintf(fp_tex, ", ");
				}
			fprintf(fp_tex, "\\\\\n");
			// fprintf(fp_tex, "\\[\n");
			// P.s_ij(0, j)->latex(fp_tex);
			// fprintf(fp_tex, "\\]\n");
			}
		fflush(fp_tex);
		}
	
	calc_lambda_ijs_matrix(t, v, k, lambda, 1 /* s */, &L);
	printf("lambda_{ij}:\n");
	L.fprint_raw(stdout);
	if (f_tex) {
		fprintf(fp_tex, "the parameters $\\lambda_{i,j}$ for $i+j \\le t=%ld$:\n", t);
		fprintf(fp_tex, "\\[\n");
		L.latex(fp_tex);
		fprintf(fp_tex, "\\]\n");
		}
	
	L.s_ij(0, 0)->copy(&b);
	vi.m_il(t + 1);
	for (i = 0; i <= t; i++) {
		Binomial(v, i, vi.s_i(i));
		}
	printf("{v \\choose i}:\n");
	vi.println();
	if (f_tex) {
		fprintf(fp_tex, "${v \\choose i}$ for $i \\le t=%ld$:\n", t);
		fprintf(fp_tex, "\\[\n");
		vi.latex(fp_tex);
		fprintf(fp_tex, "\\]\n");
		}
	
	Y.m_ilih(s, k + 1);
	for (i = 0; i <= t; i++) {
		for (j = 0; j < s; j++) {
			L.s_ij(i, 0)->copy(&c);
			c.power_int_apply(j + 1);
			vi.s_i(i)->mult(&c, Y.s_ij(i, j));
			}
		}
	for (j = 0; j < s; j++) {
		b.copy(Y.s_ij(k, j));
		}
	printf("Y:\n");
	Y.fprint_raw(stdout);
	if (f_tex) {
		fprintf(fp_tex, "$Y^{[%ld]}$:\n", s);
		fprintf(fp_tex, "\\[\n");
		Y.latex(fp_tex);
		fprintf(fp_tex, "\\]\n");
		}
	
// INT calc_lambda_max(INT t, INT v, INT k, INT lambda, SYM_OP l_max);
// INT calc_lambda_ijs_matrix(INT t, INT v, INT k, INT lambda, INT s, MATRIX_OP M);
// INT calc_lambda_ijs(INT t, INT v, INT k, INT lambda, INT s, INT i, INT j, SYM_OP lijs);
// INT calc_lambda_ij(INT t, INT v, INT k, INT lambda, INT i, INT j, SYM_OP lij);

	Bv.binomial(k, TRUE /* f_extended */, 
		TRUE /* f_inverse */, FALSE /* f_v */);
	printf("Bv:\n");
	Bv.fprint_raw(stdout);
	if (f_tex) {
		fprintf(fp_tex, "$B^{-1}$:\n");
		fprintf(fp_tex, "\\[\n");
		Bv.latex(fp_tex);
		fprintf(fp_tex, "\\]\n");
		}
	
	S1.stirling_first(s, FALSE /* f_extended */, 
		FALSE /* f_signless */, FALSE /* f_v */);
	S1.transpose(&S1t);
	printf("S1 transposed:\n");
	S1t.fprint_raw(stdout);
	if (f_tex) {
		fprintf(fp_tex, "${S^{(1)}}^\\top$:\n");
		fprintf(fp_tex, "\\[\n");
		S1t.latex(fp_tex);
		fprintf(fp_tex, "\\]\n");
		}

	D.m_ilih_n(s, s);
	for (i = 0; i < s; i++) {
		D.s_ij(i, i)->fakul_int(i + 1);
		// D.s_ij(i, i)->invers_apply();
		}
	printf("D:\n");
	D.fprint_raw(stdout);
	fflush(stdout);
	
	{
	VECTOR_OB Inv, inv, inv1, inv_val, inv_classes;
	VECTOR_OB inv_mult, inv_mult_val, inv_mult_mult;
	VECTOR_OP p_class;
	BYTE str[10000];
	
	Inv.m_il(nb_sol);
	
	for (no = 0; no < nb_sol; no++) {
		pS = (VECTOR_OP) S.s_i(no);
		
		printf("solution no %ld:\n", no + 1);
		if (f_tex) {
			fprintf(fp_tex, "solution no %ld:\\\\\n", no + 1);
			}
		design_orbits_vector(pS, &orbits, 
			FALSE /* f_complement */);
		// printf("len = %ld\n", pS->s_li());
		printf("\\frakD_{%ld}:\n", no + 1);
		for (i = 0; i < orbits.s_li(); i++) {
			printf("%ld", orbits.s_ii(i) + 1);
			if (i < orbits.s_li() - 1)
				printf(", ");
			if ((i + 1) % 10 == 0)
				printf("\n");
			}
		printf(".\n");
		if (f_tex) {
			fprintf(fp_tex, "${\\mathfrak D}_{%ld}:$\n", no + 1);
			for (i = 0; i < orbits.s_li(); i++) {
				fprintf(fp_tex, "%ld", orbits.s_ii(i) + 1);
				if (i < orbits.s_li() - 1)
					fprintf(fp_tex, ", ");
				if ((i + 1) % 10 == 0)
					fprintf(fp_tex, "\n");
				}
			fprintf(fp_tex, ".\n\n");
			}
		
		zz.m_il(k + 1);
		pz = (VECTOR_OP) zz.s_i(k);
		pP = (MATRIX_OP) P.s_ij(0, k);
		l = pP->s_li();
		pz->m_il_n(l);
		for (i = 0; i < orbits.s_li(); i++) {
			pz->m_ii(orbits.s_ii(i), 1);
			}
		printf("z_%ld:\n", k);
		for (ii = 0; ii < pz->s_li(); ii++) {
			printf("%ld", pz->s_ii(ii));
			if (ii < pz->s_li() - 1)
				printf(", ");
			}
		printf("\n");
		fflush(stdout);
		for (i = k - 1; i > t; i--) {
			pz = (VECTOR_OP) zz.s_i(i);
			((MATRIX_OP) P.s_ij(i, i + 1))->mult_vector(
				(VECTOR_OP) zz.s_i(i + 1), pz);
			b.m_i_i(k - i);
			for (ii = 0; ii < pz->s_li(); ii++) {
				pz->s_i(ii)->ganzdiv_integral(&b, &c);
				c.swap(pz->s_i(ii));
				}
			
			printf("z_%ld:\n", i);
			for (ii = 0; ii < pz->s_li(); ii++) {
				printf("%ld", pz->s_ii(ii));
				if (ii < pz->s_li() - 1)
					printf(", ");
				}
			printf("\n");
			fflush(stdout);
			if (f_tex) {
				fprintf(fp_tex, "$P_{%ld,%ld}^\\cap \\cdot {\\mathfrak x}^\\top=$ (\n", i, k);
				for (ii = 0; ii < pz->s_li(); ii++) {
					fprintf(fp_tex, "%ld", pz->s_ii(ii));
					if (ii < pz->s_li() - 1)
						fprintf(fp_tex, ", ");
					}
				fprintf(fp_tex, ")\n\n");
				}
			// pz->latex(stdout);
			for (j = 0; j < s; j++) {
				l = pz->s_li();
				pz->copy(&z1);
				for (ii = 0; ii < l; ii++) {
					z1.s_i(ii)->power_int_apply(j + 1);
					}
				((MATRIX_OP) P.s_ij(0, i))->mult_vector(&z1, &z2);
				if (z2.s_li() != 1)
					return error("z2.s_li() != 1");
				z2.s_i(0)->copy(Y.s_ij(i, j));
				printf("Y_{%ld,%ld}=", i, j);
				Y.s_ij(i, j)->println();
				fflush(stdout);
				}
			}
		printf("Y:\n");
		Y.fprint_raw(stdout);
		if (f_tex) {
			MATRIX_OB YY;
			INT i, j, m;
			
			m = k - (t + 1);
			YY.m_ilih(s, m);
			for (i = 0; i < m; i++) {
				for (j = 0; j < s; j++) {
					Y.s_ij(t + 1 + i, j)->copy(YY.s_ij(i, j));
					}
				}
			fprintf(fp_tex, "$Y^{[%ld]}_{[%ld\\ldots %ld]}$:\n", s, t + 1, k - 1);
			fprintf(fp_tex, "\\[\n");
			YY.latex(fp_tex);
			fprintf(fp_tex, "\\]\n");
			}
		
		Bv.mult(&Y, &As1);

		printf("As1:\n");
		As1.fprint_raw(stdout);

		As1.mult(&S1t, &As2);
			{
			SYM_OP fact;
			SYM_OB tmp;
			
			for (j = 1; j < s; j++) {
				fact = D.s_ij(j, j);
				for (i = 0; i <= k; i++) {
					As2.s_ij(i, j)->ganzdiv_integral(fact, &tmp);
					tmp.swap(As2.s_ij(i, j));
					}
				}
			}
		printf("As2:\n");
		As2.fprint_raw(stdout);
		
		inv.m_il(k - 1 - t);
		for (i = 0; i < k - 1 - t; i++) {
#if 1
			inv1.m_il(s - 1);
			for (j = 1; j < s; j++) 
				As2.s_ij(t + 1 + i, j)->copy(inv1.s_i(j - 1));
			inv1.swap(inv.s_i(i));
#else
			((MATRIX_OP) P.s_ij(t + 1 + i, k))->mult_vector(pS, &z);
			//printf("z_%ld:\n", t + 1 + i);
			z.copy((VECTOR_OP) inv.s_i(i));
#endif
			}
		if (f_tex) {
			fprintf(fp_tex, "invariants:\n");
			fprintf(fp_tex, "\\[\n");
			inv.latex(fp_tex);
			fprintf(fp_tex, "\\]\n");
			}
		inv.swap(Inv.s_i(no));

		if (f_tex) {
			INT u;
			
			fprintf(fp_tex, "$A^{(%ld)}$:\n", s);
			fprintf(fp_tex, "\\begin{align*}\n");
			fprintf(fp_tex, "\\begin{array}{r|*{%ld}{r}}\n", s);
			fprintf(fp_tex, "i ");
			for (u = 1; u <= s; u++) {
				fprintf(fp_tex, "& \\alpha_i^{(%ld)}({\\cal D}) ", u);
				} 
			fprintf(fp_tex, "\\\\\n");
			fprintf(fp_tex, "\\hline\n");
			for (i = 0; i <= k; i++) {
				fprintf(fp_tex, "%ld ", i);
				for (u = 0; u < s; u++) {
					fprintf(fp_tex, "& ");
					As2.s_ij(i, u)->latex(fp_tex);
					} 
				fprintf(fp_tex, "\\\\\n");
				}
			fprintf(fp_tex, "\\end{array}\n", s);
			fprintf(fp_tex, "\\end{align*}\n", s);
			// fprintf(fp_tex, "\\[\n");
			// As2.latex(fp_tex);
			// fprintf(fp_tex, "\\]\n");
			}
		
		} // next no
	
	printf("design invariants:\n");
	for (no = 0; no < nb_sol; no++) {
		printf("%ld: ", no + 1);
		Inv.s_i(no)->println();
		}
	Inv.classify(&inv_val, &inv_classes);
	for (i = 0; i < inv_val.s_li(); i++) {
		inv_val.s_i(i)->print();
		printf("for $\\{");
		p_class = (VECTOR_OP) inv_classes.s_i(i);
		for (j = 0; j < p_class->s_li(); j++) {
			printf("%ld", p_class->s_ii(j) + 1);
			if (j < p_class->s_li() - 1)
				printf(", ");
			}
		printf("\\}$\n");
		}
	if (f_tex) {
		fprintf(fp_tex, "classification by $A^{(%ld)}$:\\\\\n", s);
		for (i = 0; i < inv_val.s_li(); i++) {
			fprintf(fp_tex, "$");
			inv_val.s_i(i)->latex(fp_tex);
			fprintf(fp_tex, "$");
			fprintf(fp_tex, "for $\\{");
			p_class = (VECTOR_OP) inv_classes.s_i(i);
			for (j = 0; j < p_class->s_li(); j++) {
				fprintf(fp_tex, "%ld", p_class->s_ii(j) + 1);
				if (j < p_class->s_li() - 1)
					fprintf(fp_tex, ", ");
				}
			fprintf(fp_tex, "\\}$\\\\\n");
			}
		}
	inv_classes.norm(&inv_mult);
	inv_mult.multiplicities(&inv_mult_val, &inv_mult_mult);
	str[0] = 0;
	inv_mult_val.sprint_multiplicities(&inv_mult_mult, str);
	printf("multiplicities: %s\n", str);
	if (f_tex) {
		fprintf(fp_tex, "multiplicities: $%s$\\\\\n", str);
		}
	}

	return OK;
}




