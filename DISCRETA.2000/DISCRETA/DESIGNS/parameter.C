/* parameter.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef LADDER_TRUE

#include <DISCRETA/perm.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/lb.h>
#include <DISCRETA/cp.h>
#include <DISCRETA/divs.h>
#include <DISCRETA/ladder.h>

#include <stdlib.h>
#include <ctype.h> // for isspace()



#define DP_MIN_T 4
#define DP_ALLOC_LEN 5
#define MAX_DP 200
#define DP_MAX_LAMBDA 150000000
#define BUFSIZE 20000


static INT db_dp_HTML_d1(DATABASE_OP db, BAYERTREE_OP btree);
static INT db_dp_HTML_d2(DATABASE_OP db, BAYERTREE_OP btree);
static INT db_dp_HTML_d3(DATABASE_OP db, BAYERTREE_OP btree);
static INT db_dp_HTML_d4(DATABASE_OP db, BAYERTREE_OP btree);
void calc_prev_string(BYTE *s, INT prev);
static INT write_head_HTML(FILE *fp, DESIGN_PARAMETER_OP dp);
static INT write_source_HTML(FILE *fp, DESIGN_PARAMETER_OP dp);
void print_page_head(FILE *fp);
void copyright(FILE *fp);



#if TEXDOCU
INT calc_lambda_max(INT t, INT v, INT k, INT lambda, SYM_OP l_max)
#endif
{
	Binomial(v - t, k - t, l_max);
	return OK;
}

#if TEXDOCU
INT calc_lambda_ijs_matrix(INT t, INT v, INT k, INT lambda, INT s, MATRIX_OP M)
#endif
{
	INT i, j;
	
	M->m_ilih(t + 1, t + 1);
	for (i = 0; i <= t; i++) {
		for (j = 0; j <= t - i; j++) {
			calc_lambda_ijs(t, v, k, lambda, s, i, j, M->s_ij(i, j));
			}
		}
	return OK;
}

#if TEXDOCU
INT calc_lambda_ijs(INT t, INT v, INT k, INT lambda, INT s, INT i, INT j, SYM_OP lijs)
#else
\[
\lambda_{i,j}^{(s)} = \sum_{h=0}^j (-1)^h {j \choose h} {\lambda_{i+h} \choose s}
\]
cf. Wilson, Van Lint~\cite{VanLintWilson92}.
#endif
{
	SYM_OB res, a, b, c, d, e;
	INT h;
	
	res.m_i_i(0);
	for (h = 0; h <= j; h++) {
		Binomial(j, h, &a);
		if (ODD(h))
			a.addinvers_apply();
		calc_lambda_ij(t, v, k, lambda, i + h, 0, &b);
		N_choose_K(&b, s, &c);
		a.mult(&c, &d);
		res.add(&d, &e);
		e.swap(&res);
		}
	res.copy(lijs);
	return OK;
}

#if TEXDOCU
INT calc_lambda_ij(INT t, INT v, INT k, INT lambda, INT i, INT j, SYM_OP lij)
#else
\[
\lambda_{i,j} = \lambda \frac{{v-i-j \choose k-i}}{{v-t \choose k-t}}
\]
cf. Wilson, Van Lint~\cite{VanLintWilson92}.
#endif
{
	SYM_OB a, b, c, d;
	
	Binomial(v - i - j, k - i, &a);
	b.m_i_i(lambda);
	a.mult(&b, &c);
	Binomial(v - t, k - t, &d);
	c.ganzdiv(&d, lij);
	// printf("calc_lambda_ij: i=%ld j=%ld \\lambda_{i,j}=", i, j);
	// lij->println();
	return OK;
}

#if TEXDOCU
INT calc_mendelsohn_coefficient_matrix(INT t, INT m, MATRIX_OP M)
#else
The Mendelsohn equations for any $t$-$(v,k,\lambda)$ design $\cD = (\cV, \cB)$  
and any $m$-subset $M \subseteq \cV$ are for $s \ge 1$:
\[
\sum_{j=i}^m {m \choose j} \alpha_j^{(s)}(M) = 
{\lambda_i \choose s} {m \choose i} \quad \text{for} i=0,\ldots,t 
\]
cf. Mendelsohn~\cite{Mendelsohn71}.
#endif
{
	INT i, j;
	
	M->m_ilih_n(m + 1, t + 1);
	for (i = 0; i <= t; i++) {
		for (j = i; j <= m; j++) {
			Binomial(j, i, M->s_ij(i, j));
			}
		}
	return OK;
}

#if TEXDOCU
INT calc_mendelsohn_rhs(INT v, INT t, INT k, INT lambda, INT m, INT s, VECTOR_OP rhs)
#endif
{
	INT i;
	SYM_OB a, b, c;
	
	rhs->m_il_n(t + 1);
	for (i = 0; i <= t; i++) {
		Binomial(m, i, &a);
		calc_lambda_ij(t, v, k, lambda, i, 0, &b);
		N_choose_K(&b, s, &c);
		a.mult(&c, rhs->s_i(i));
		}
	return OK;
}



#if TEXDOCU
INT calc_and_print_design_parameter(INT t, INT v, INT k, INT lambda)
#endif
{
	SYM_OB b_ob, r_ob, bb_ob, bbb_ob, two_ob, three_ob, b2_ob, b3_ob;
	
	design_calc_b(t, v, k, lambda, &b_ob);
	printf("b = ");
	b_ob.print();
	printf("\\\\\n");
	fflush(stdout);
	
	design_calc_r(v, k, &b_ob, &r_ob);
	printf("r = ");
	r_ob.print();
	printf("\\\\\n");
	fflush(stdout);

	b_ob.mult(&b_ob, &bb_ob);
	printf("$b^2$ = ");
	bb_ob.print();
	printf("\\\\\n");
	fflush(stdout);
	
	two_ob.m_i_i(2);
	n_choose_k_ob(&b_ob, &two_ob, &b2_ob);
	printf("${b \\choose 2}$ = ");
	b2_ob.print();
	printf("\\\\\n");
	fflush(stdout);

	bb_ob.mult(&b_ob, &bbb_ob);
	printf("$b^3$ = ");
	bbb_ob.print();
	printf("\\\\\n");
	fflush(stdout);
	
	three_ob.m_i_i(3);
	n_choose_k_ob(&b_ob, &three_ob, &b3_ob);
	printf("${b \\choose 3}$ = ");
	b3_ob.print();
	printf("\\\\\n");
	fflush(stdout);
	return OK;
}
		
#if TEXDOCU
INT design_calc_b(INT t, INT v, INT k, INT lambda, SYM_OP b)
#endif
{
	SYM_OB x, y, v_ob, t_ob, k_ob, lambda_ob, tmp;

	v_ob.m_i_i(v);
	t_ob.m_i_i(t);
	k_ob.m_i_i(k);
	n_choose_k_ob(&v_ob, &t_ob, &x);
	lambda_ob.m_i_i(lambda);
	x.mult(&lambda_ob, &tmp);
	n_choose_k_ob(&k_ob, &t_ob, &y);
	tmp.ganzdiv(&y, b);
	return OK;
}

#if TEXDOCU
INT design_calc_r(INT v, INT k, SYM_OP b, SYM_OP r)
#endif
{
	SYM_OB v_ob, k_ob, tmp;

	v_ob.m_i_i(v);
	k_ob.m_i_i(k);
	k_ob.mult(b, &tmp);
	tmp.ganzdiv(&v_ob, r);
	return OK;
}

#if TEXDOCU
INT design_print_Mendelsohn_and_Koehler(INT v, INT t, INT k, MATRIX_OP Mendelsohn_mtx, 
	VECTOR_OP Mendel_RHS, MATRIX_OP LAmbda, VECTOR_OP alpha_i, 
	VECTOR_OP Coeff, VECTOR_OP Constant_term, 
	VECTOR_OP Koehler_Constant_term_v, VECTOR_OP Koehler_Coeff_v)
#else
Koehlers equations are due to Koehler~\cite{Koehler89}. 
They have been generalized by Tran van Trung et al.~\cite{TranVanTrungetal96}.
#endif
{
	VECTOR_OP RHS_1, RHS_s;
	VECTOR_OP constant_term, coeff;
	INT i, j, m1 = k, s, s_max;
	
	s_max = Mendel_RHS->s_li() - 1;
	RHS_1 = (VECTOR_OP) Mendel_RHS->s_i(1);
	printf("The {\\sc Mendelsohn} system:\\\\\n");
	printf("\\begin{equation}\n");
	Mendelsohn_mtx->latex_upper_tri(stdout);
	printf("\\left(\n");
	printf("{\\arraycolsep=2pt\n");
	printf("\\begin{array}{r}\n");
	for (i = 0; i <= m1; i++) {
		if (i < 10)
			printf("\\alpha_%ld\\\\\n", i);
		else
			printf("\\alpha_{%ld}\\\\\n", i);
		}
	printf("\\end{array}\n");
	printf("}\n");
	printf("\\right)\n");
	printf("\\; = \\; \n");
	printf("\\left(\n");
	printf("\\begin{array}{c}\n");
	for (i = 0; i <= t; i++) {
		RHS_1->s_i(i)->latex(stdout);
		printf("\\\\\n");
		}
	printf("\\end{array}\n");
	printf("\\right)\n");
	printf("\\; = \\; \n");
	printf("\\left(\n");
	printf("\\begin{array}{rr}\n");
	for (i = 0; i <= t; i++) {
		printf("{ %ld \\choose %ld } & ", m1, i);
		LAmbda->s_ij(i, 0)->latex(stdout);
		printf("\\\\\n");
		}
	printf("\\end{array}\n");
	printf("\\right)\n");
	printf("\\end{equation}\n");

	if (k == t + 1) {
		solve_Mendelsohn(v, t, Mendelsohn_mtx, RHS_1, alpha_i);
		printf("The (unique) solution:\n");
		printf("\\begin{equation}\n");
		printf("\\left(\n");
		printf("\\begin{array}{r}\n");
		for (i = 0; i <= m1; i++) {
			alpha_i->s_i(i)->latex(stdout);
			printf("\\\\\n");
			}
		printf("\\end{array}\n");
		printf("\\right)\n");
		printf("\\end{equation}\n");
		}
	
	printf("The triangle\\\\\n");
	printf("{\\small\n");
	printf("{\\arraycolsep=2pt\n");
	printf("\\begin{equation}\n");
	// LAmbda->latex(stdout);
	printf("\\begin{array}{*{%ld}r}\n", t + 1);
	for (i = 0; i <= t; i++) {
		for (j = 0; j <= t; j++) {
			if (j <= t - i)
				LAmbda->s_ij(i, j)->latex(stdout);
			if (j < t)
				printf(" & ");
			}
		printf("\\\\\n");
		}
	printf("\\end{array}\n");
	printf("\\end{equation}\n");
	printf("}%%\n");
	printf("}%%\n");
	
	for (s = 1; s <= 1 /* s_max */; s++) {
		printf("The {\\sc K\\\"ohler} equations (for $B \\in {\\cal B}$ and $s=%ld$):\\\\\n", s);
		printf("\\begin{equation}\n");
		printf("\\begin{array}{*{%ld}r}\n", m1 - t + 3);
		for (j = 0; j <= t; j++) {
			constant_term = (VECTOR_OP) Constant_term->s_i(j);
			coeff = (VECTOR_OP) Coeff->s_i(j);
			Koehler_eqn_print(v, t, m1, s, j, constant_term, coeff, 
				"\\alpha", "(B)");
			printf("\\\\\n");
			}
		printf("\\end{array}\n");
		printf("\\end{equation}\n");
		} // next s
	
	if (s_max > 1) {
		printf("The generalized Mendelsohn systems have "
			"the following right side ($s \\le %ld$)\\\\\n", s_max);
		printf("\\begin{equation}\n");
		for (s = 2; s <= s_max; s++) {
			RHS_s = (VECTOR_OP) Mendel_RHS->s_i(s);
			printf("\\left(\n");
			printf("\\begin{array}{r}\n");
			for (i = 0; i <= t; i++) {
				RHS_s->s_i(i)->latex(stdout);
				printf("\\\\\n");
				}
			printf("\\end{array}\n");
			printf("\\right)");
			if (s < s_max)
				printf(", \n");
			}
		printf("\\end{equation}\n");
		}
	
	for (s = 1; s <= s_max; s++) {
		printf("The {\\sc K\\\"ohler} equations applied to "
			"$M=V$ (with $m=v$) and $s=%ld$ are ", s);
		printf("(the $\\alpha_j({\\cal B})$-terms with "
			"$j=k+1,\\ldots,v$ left out):\\\\\n");
		printf("\\begin{equation}\n");
		printf("\\begin{array}{*{%ld}r}\n", k /* v */ - t + 3);
		for (j = 0; j <= t; j++) {
			constant_term = (VECTOR_OP) Koehler_Constant_term_v->s_i(j);
			coeff = (VECTOR_OP) Koehler_Coeff_v->s_i(j);
			Koehler_eqn_print(v, t, k /* v */, s, j, 
				constant_term, coeff, "\\alpha", "({\\cal B})");
			printf("\\\\\n");
			}
		printf("\\end{array}\n");
		printf("\\end{equation}\n");
		} // next s
	
	return OK;
}

#if TEXDOCU
INT solve_Mendelsohn(INT v, INT t, MATRIX_OP M, 
	VECTOR_OP RHS, VECTOR_OP alpha_i)
#endif
{
	INT i, j;
	SYM_OP pb;
	SYM_OB tmp1, tmp2, coeff;
	VECTOR_OB bb;
	
	if (M->s_li() != t + 2)
		return error("solve_Mendelsohn() M->s_li() != t + 2");
	if (RHS->s_li() != t + 1)
		return error("solve_Mendelsohn() RHS->s_li() != t + 1");
	RHS->copy(&bb);
	bb.inc();
	alpha_i->m_il(t + 2);
	alpha_i->m_ii(t + 1, 1);
	bb.m_ii(t + 1, 1);
	for (j = t + 1; j >= 0; j--) {
		bb.s_i(j)->copy(alpha_i->s_i(j));
		alpha_i->s_i(j)->copy(&coeff);
		coeff.addinvers_apply();
		for (i = j; i >= 0; i--) {
			if (i < t + 1)
				coeff.mult(M->s_ij(i, j), &tmp1);
			else
				coeff.copy(&tmp1);
			pb = bb.s_i(i);
			pb->add(&tmp1, &tmp2);
			tmp2.swap(pb);
			}
		}
	return OK;
}

#if TEXDOCU
INT Mendelsohn_generalized_RHS(INT v, INT t, INT k, INT m, INT s, INT lambda, 
	VECTOR_OP RHS, INT f_v)
#endif
{
	INT i;
	SYM_OB a;
	
	if (s < 1)
		return error("Mendelsohn_generalized_RHS() s < 1");
	calc_Lambda(v, t, k, lambda, RHS);
	if (s == 1)
		return OK;
	for (i = 0; i <= t; i++) {
		N_choose_K(RHS->s_i(i), s, &a);
		a.swap(RHS->s_i(i));
		}
	return OK;
}

#if TEXDOCU
INT Mendelsohn(INT v, INT t, INT k, INT m, INT lambda, 
	MATRIX_OP M, VECTOR_OP RHS, MATRIX_OP LAmbda, INT f_v)
#endif
{
	INT i, j, a;
	SYM_OB c, d;

	if (f_v) {
		printf("Mendelsohn():\n");
		printf("v = %ld, t = %ld, k = %ld, m = %ld, lambda = %ld\n", 
			v, t, k, m, lambda);
		}
	calc_Lambda(v, t, k, lambda, RHS);
	calc_LAmbda(v, t, k, lambda, LAmbda);
	if (f_v) {
		printf("Lambda = ");
		RHS->println();
		}
	for (j = 0; j <= t; j++) {
		a = n_choose_k(m, j);
		c.m_i_i(a);
		RHS->s_i(j)->mult(&c, &d);
		d.swap(RHS->s_i(j));
		// b = RHS->s_ii(j) * a;
		// RHS->m_ii(j, b);
		}
	if (f_v) {
		printf("RHS = ");
		RHS->println();
		}
	M->m_ilih_n(m + 1, t + 1);
	for (j = 0; j <= t; j++) {
		for (i = j; i <= m; i++) {
			a = n_choose_k(i, j);
			M->m_iji(j, i, a);
			}
		}
	if (f_v) {
		printf("M = \n");
		// M->Print();
		M->fprint_raw(stdout);
		printf("LAmbda = \n");
		LAmbda->fprint_raw(stdout);
		}
	return OK;
}

#if TEXDOCU
INT calc_LAmbda(INT v, INT t, INT k, INT lambda, MATRIX_OP LAmbda)
#endif
{
	VECTOR_OB Lambda;
	INT i, j;
	SYM_OB a, b, c;

	calc_Lambda(v, t, k, lambda, &Lambda);
	LAmbda->m_ilih_n(t + 1, t + 1);
	for (i = 0; i <= t; i++)
		Lambda.s_i(i)->copy(LAmbda->s_ij(i, 0));
	for (j = 1; j <= t; j++) {
		for (i = 0; i <= t - j; i++) {
			LAmbda->s_ij(i, j - 1)->copy(&a);
			LAmbda->s_ij(i + 1, j - 1)->copy(&b);
			b.addinvers_apply();
			a.add(&b, &c);
			c.copy(LAmbda->s_ij(i, j));
			// c = a - b;
			// LAmbda->m_iji(i, j, c);
			}
		}
	return OK;
}

#if TEXDOCU
INT calc_Lambda(INT v, INT t, INT k, INT lambda, VECTOR_OP Lambda)
#endif
{
	INT i;
	SYM_OB a, g, x, y, c, d;
	
	Lambda->m_il(t + 1);
	Lambda->m_ii(t, lambda);
	for (i = t; i >= 0; i--) {
		if (i == t) {
			c.m_i_i(lambda);
			}
		else {
			x.m_i_i(v - i);
			y.m_i_i(k - i);
			c.mult(&x, &a);
			a.quores(&y, &c, &d);
			// ggt_iipi(a, b, &g);
			// a /= g;
			// b /= g;
			if (!d.nullp()) {
				printf("non-integer division !\nv=%ld t=%ld k=%ld lambda=%ld\n", v, t, k, lambda);
				fflush(stdout);
				return error("calc_Lambda() d != 1, parameter set not admissible");
				}
			c.copy(Lambda->s_i(i));
			}
		}
	return OK;
}

#if TEXDOCU
INT calc_delta_lambda(INT v, INT t, INT k)
#endif
{
	INT lambda;
	INT i, a, b, g, rhs_a = 1, rhs_b = 1, delta_lambda = 1, dl;

	lambda = 1;
	printf("calc_delta_lambda(): v=%ld t=%ld k=%ld lambda=%ld\n", 
		v, t, k, lambda);
	fflush(stdout);
	for (i = t; i >= 0; i--) {
		if (i == t) {
			rhs_a = lambda;
			rhs_b = 1;
			delta_lambda = 1;
			}
		else {
			a = rhs_a * (v - i);
			b = rhs_b * (k - i);
			ggt_iipi(a, b, &g);
			a /= g;
			b /= g;
			kgv_iipi(delta_lambda, b, &dl);
			delta_lambda = dl;
			printf("t'=%ld lambda'=%ld/%ld delta_lambda=%ld\n", 
				i, a, b, delta_lambda);
			fflush(stdout);
			rhs_a = a;
			rhs_b = b;
			}
		}
	return delta_lambda;
}

#if TEXDOCU
INT Koehler_eqns_for_blocks(INT v, INT t, INT k, INT lambda, INT s_max, 
	INT max_intersection, 
	VECTOR_OP Coeff, VECTOR_OP Constant_term, BYTE *variable_name, INT f_v)
#endif
{
	INT j;
	VECTOR_OB constant_term, coeff, Lambda;
	SYM_OP alpha_k, c;
	SYM_OB a;

	calc_Lambda(v, t, k, lambda, &Lambda);
	Coeff->m_il(t + 1);
	Constant_term->m_il(t + 1);
	for (j = 0; j <= t; j++) {
		Koehler_j(v, t, k /* m */, s_max, j, &Lambda, &constant_term, &coeff);
		alpha_k = coeff.s_i(k);
		c = constant_term.s_i(1); /* for s = 1 */
		alpha_k->add(c, &a);
		a.swap(c);
		coeff.realloc_z(max_intersection + 1);
		if (f_v) {
			Koehler_eqn_print(v, t, max_intersection /* m */, s_max, j, 
				&constant_term, &coeff, variable_name, "" /* argument */);
			printf("\n");
			}
		constant_term.swap((VECTOR_OP) Constant_term->s_i(j));
		coeff.swap((VECTOR_OP) Coeff->s_i(j));
		}
	return OK;
}

#if TEXDOCU
INT Koehler(INT v, INT t, INT k, INT lambda, INT m, INT s_max, 
	VECTOR_OP Coeff, VECTOR_OP Constant_term, INT f_v)
#endif
{
	INT j;
	VECTOR_OB constant_term, coeff, Lambda;

	calc_Lambda(v, t, k, lambda, &Lambda);
	Coeff->m_il(t + 1);
	Constant_term->m_il(t + 1);
	for (j = 0; j <= t; j++) {
		Koehler_j(v, t, m, s_max, j, &Lambda, &constant_term, &coeff);
		if (f_v) {
			Koehler_eqn_print(v, t, m, s_max, j, 
				&constant_term, &coeff, "\\alpha", "" /* argument */);
			printf("\n");
			}
		constant_term.swap((VECTOR_OP) Constant_term->s_i(j));
		coeff.swap((VECTOR_OP) Coeff->s_i(j));
		}
	return OK;
}

#if TEXDOCU
INT Koehler_j(INT v, INT t, INT m, INT s_max, INT j, VECTOR_OP Lambda, 
	VECTOR_OP constant_term, VECTOR_OP coeff)
#endif
{
	INT s, r, sign;
	SYM_OB tmp, tmp1, tmp2, c0, lambda_atop_s;
	SYM_OB m_ob, r_ob, a_ob, b_ob, c_ob;
	
	if (j < 0)
		return error("Koehler(): j < 0");
	if (j > t)
		return error("Koehler(): j > t");
	m_ob.m_i_i(m);
	constant_term->m_il_n(s_max + 1);
	coeff->m_il_n(m + 1);
	for (s = 1; s <= s_max; s++) {
		c0.m_i_i(0);
		for (r = j; r <= t; r++) {
			r_ob.m_i_i(r);
			if (EVEN(r + j))
				sign = 1;
			else
				sign = -1;
			// a = n_choose_k(r, j);
			// b = n_choose_k(m, r);
			N_choose_K(&r_ob, j, &a_ob);
			N_choose_K(&m_ob, r, &b_ob);
			a_ob.mult(&b_ob, &c_ob);
			if (sign == -1) {
				c_ob.addinvers_apply();
				}
			// c = sign * a * b;
			// tmp.m_i_i(c);
			N_choose_K(Lambda->s_i(r), s, &lambda_atop_s);
			c_ob.mult(&lambda_atop_s, &tmp1);
			c0.add(&tmp1, &tmp2);
			tmp2.swap(&c0);
			}
		c0.swap(constant_term->s_i(s));
		}
	for (r = 0; r < m - t; r++) {
		if (EVEN(t + j + 1))
			sign = 1;
		else
			sign = -1;
		tmp1.m_i_i(t + r - j);
		N_choose_K(&tmp1, r, &a_ob);
		// a = n_choose_k(t + r - j, r);
		tmp1.m_i_i(t + r + 1);
		N_choose_K(&tmp1, j, &b_ob);
		// b = n_choose_k(t + r + 1, j);
		a_ob.mult(&b_ob, &c_ob);
		if (sign == -1) {
			c_ob.addinvers_apply();
			}
		// c = sign * a * b;
		c_ob.copy(coeff->s_i(t + r + 1));
		// coeff->m_ii(t + r + 1, c);
		}
	return OK;
}

#if TEXDOCU
INT Koehler_eqn_print(INT v, INT t, INT m, INT s_max, INT j, 
	VECTOR_OP constant_term, VECTOR_OP coeff, 
	BYTE *variable_name, BYTE *argument)
#endif
{
	INT s, r, a;
	
#if 0
	if (s == 1) {
		printf("\\alpha_{%ld}%s  & =  & ", j, argument);
		}
	else {
		printf("\\alpha_{%ld}^{(%ld)}%s  & =  & ", j, s, argument);
		}
#endif
	printf("%s_{%ld}^{(s)}%s  & =  & ", variable_name, j, argument);
	if (s_max > 1)
		printf("(");
	for (s = 1; s <= s_max; s++) {
		constant_term->s_i(s)->print();
		if (s < s_max)
			printf(", ");
		}
	if (s_max > 1)
		printf(")");
	for (r = 0; r < m - t; r++) {
		a = coeff->s_ii(t + r + 1);
		printf(" & ");
		if (a >= 0)
			printf(" + ");
		printf("%ld \\, %s_", a, variable_name);
		if (t + r + 1 >= 10)
			printf("{%ld}", t + r + 1);
		else
			printf("%ld", t + r + 1);
		// if (s > 1)
		//	printf("^{(%ld)}", s);
		printf("^{(s)}");
		printf("%s ", argument);
		}
	// printf("\n");
	return OK;
}

#if TEXDOCU
void design_print_solution_vector(BYTE *buf, INT no, INT width, BYTE *offset)
#endif
{
	INT i, l, l1;
	BYTE c;
	
	l = strlen(buf);
	for (i = 0; i < l; i+=width) {
		l1 = l - i;
		if (l1 > width)
			l1 = width;
		c = buf[i + l1];
		buf[i + l1] = 0;
		printf("%%{\\tiny\n");
		if (i == 0)
			printf("%%%ld: ", no);
		if (i == 0 && l > width)
			printf("\\\\[%s]\n", offset);
		printf("%%%s}\\\\[%s]\n", buf + i, offset);
		buf[i + l1] = c;
		}
}


#if TEXDOCU
void design_print_solution_vector_numerical(BYTE *buf, INT no, BYTE *offset)
#endif
{
	INT i, l, f_first = TRUE;
	BYTE c;
	
	// printf("%% ");
	printf("\n\n");
	printf("${\\mathfrak D}_{%ld}$: ", no);
	l = strlen(buf);
	for (i = 0; i < l; i++) {
		c = buf[i];
		if (buf[i] == '0')
			continue;
		if (f_first)
			f_first = FALSE;
		else
			printf(", ");
		printf("%ld", i + 1);
		}
	printf("\\\\\n");
}

/*
 * DESIGN_PARAMETER 
 */

#if TEXDOCU
#else
The DESIGN\_PARAMETER class implements a class for design parameter sets.
This means that we store the values t, v, k, lambda of a design 
together with an id (a unique number for the specific parameter set) 
and the construction which led to the design. A string comment may also be given. 
In the discreta application, we have a file 
design.txt which contains the initial design parameter sets. 
Via various constructions, this basis is used to produce a lot of 
designs and parameter sets. The constructions are, for example, 
derived design, residual design, design with respect to smaller t and so. 
Complementary designs are treated in aspecial way: 
Not that we have two ways of forming the complementary design: 
We can choose the block set ${ V \choose k} \backslash {\cal B}$ 
with $\lambda = \sigma_Z - \lambda$ 
$\sigma_Z := {{v - t} \choose {k - t}}$ 
(construction CI). 
It is also possible to choose $\{ V - B | B \in {\cal B} \}$ as the set of 
blocks of a design (construction CII). 
In the database, only CI is considered. Because the complementary design 
gives nothing new, only one of a pair of complementary designs is stored. 
If a design is selfcomplementary ($\lambda = \sigma_Z / 2$), 
this is indicated in the comment.

Not all constructions are implemented or used at the moment: 
for instance, Alltops constructions are available but not used at the moment.
Alltops construction is described in Alltop~\cite{Alltop75}.

Sometimes, large integers seem to be a difficulty for handling with 
design parameter sets. Note carefully the implementation of design 
parameter calculation using DISCRETA objects. Thus it is possible 
to compute with long integers, too, because of the facility of 
switching to long integers in DISCRETA. 

There should be more information about special properties of the designs 
in the comment line. One could imagine \lq tight design\rq or so.

The main construction routine for derived / residual / etc. parameter sets 
is implemented twice. Once, here, without a database. 
The second time, using a database to store the data. 
This is the version which is used for discreta.
#endif

#if TEXDOCU
INT DESIGN_PARAMETER_SOURCE_OB::field_name(INT i, INT j, BYTE *str)
#endif
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "prev"; break;
	case 1: s = "rule"; break;
	case 2: s = "comment"; break;
	case 3: s = "references"; break;
	default:
		return error("DESIGN_PARAMETER_SOURCE::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

#if TEXDOCU
INT DESIGN_PARAMETER_SOURCE_OB::init()
#endif
{
	INT erg = OK;
	
	erg += m_il(4);
	c_obj_k(DESIGN_PARAMETER_SOURCE_KIND);
	s_prev()->m_i(-1);
	s_rule()->m_i(-1);
	s_comment()->init("");
	s_references()->m_il(0);
	return OK;
}

#if TEXDOCU
INT DESIGN_PARAMETER_SOURCE_OB::sprint(BYTE *s)
#endif
{
	char str[10000];
	char str0[10000];
	char str1[10000];
	char str2[10000];
	INT prev;
	// INT rule, i, l;


	sprint_text012(str0, str1, str2);
	prev = s_prev_i();
	sprintf(str, "%s", str1);
	if (prev != -1) {
		sprintf(str + strlen(str), " %ld", prev);
		}
	sprintf(str + strlen(str), " %s", str2);

#if 0
	rule = s_rule_i();
	prev = s_prev_i();
	str[0] = 0;
	if (rule == DP_RULE_COMPLEMENT) {
		sprintf(str + strlen(str), "complementary design of number %ld", prev);
		}
	else if (rule == DP_RULE_REDUCED_T) {
		sprintf(str + strlen(str), "design of number %ld with respect to smaller t", prev);
		}
	else if (rule == DP_RULE_DERIVED) { 
		sprintf(str + strlen(str), "derived from number %ld", prev);
		}
	else if (rule == DP_RULE_RESIDUAL) {
		sprintf(str + strlen(str), "residual design of number %ld", prev);
		}
	else if (rule == DP_RULE_ALLTOP) {
		sprintf(str + strlen(str), "Alltop construction for design number %ld", prev);
		}
	else if (rule == DP_RULE_TRUNG_SUPPLEMENTARY) {
		sprintf(str + strlen(str), "Tran van Trung construction with supplementary design for design number %ld", prev);
		}
	else if (rule == DP_RULE_COMPL_REDUCED_T) {
		sprintf(str + strlen(str), "complementary design of number %ld with respect to smaller t", prev);
		}
	else if (rule == DP_RULE_COMPL_DERIVED) { 
		sprintf(str + strlen(str), "derived from complement of number %ld", prev);
		}
	else if (rule == DP_RULE_COMPL_RESIDUAL) {
		sprintf(str + strlen(str), "residual design of complement of number %ld", prev);
		}
	else if (rule == DP_RULE_COMPL_ALLTOP) {
		sprintf(str + strlen(str), "Alltop construction for complementary design of number %ld", prev);
		}
	else {
		sprintf(str + strlen(str), "no rule ");
		}
	l = s_references()->s_li();
	if (l > 0) {
		sprintf(str + strlen(str), "( ");
		for (i = 0; i < l; i++) {
			sprintf(str + strlen(str), "%s", s_references_is(i));
			if (i < l - 1)
				sprintf(str + strlen(str), ", ");
			}
		sprintf(str + strlen(str), " )");
		}
#endif
	strcat(s, str);
	return OK;
}

#if TEXDOCU
INT DESIGN_PARAMETER_SOURCE_OB::sprint_text012(BYTE *s0, BYTE *s1, BYTE *s2)
#else
special print function for the design construction strings: 
three strings are generated which give parts of an english sentence. 
The wholes between the strings may be filled with numbers (ids) of 
parameter sets. In the html page writing routines, these numbers 
are printed using htmls capability of including links directly to 
the definition of the other parameter set.
#endif
{
	INT rule;

	s1[0] = 0;
	if (strlen(s_comment_s()) > 0) {
		sprintf(s1 + strlen(s1), "%s ", s_comment_s());
		}
	rule = s_rule_i();
	if (rule == DP_RULE_COMPLEMENT) {
		sprintf(s1 + strlen(s1), "complementary design of number");
		sprintf(s2, "");
		}
	else if (rule == DP_RULE_REDUCED_T) {
		sprintf(s1 + strlen(s1), "design of number");
		sprintf(s2, "with respect to smaller t");
		}
	else if (rule == DP_RULE_DERIVED) { 
		sprintf(s1 + strlen(s1), "derived from number");
		sprintf(s2, "");
		}
	else if (rule == DP_RULE_RESIDUAL) {
		sprintf(s1 + strlen(s1), "residual design of number");
		sprintf(s2, "");
		}
	else if (rule == DP_RULE_ALLTOP) {
		sprintf(s1 + strlen(s1), "Alltop construction for design number");
		sprintf(s2, "");
		}
	else if (rule == DP_RULE_SUPPLEMENTARY) {
		sprintf(s1 + strlen(s1), "supplementary design of number");
		sprintf(s2, "");
		}
	else if (rule == DP_RULE_TRUNG_SUPPLEMENTARY) {
		sprintf(s1 + strlen(s1), "Tran van Trung construction with supplementary design for design number");
		sprintf(s2, "");
		}
	else if (rule == DP_RULE_COMPL_REDUCED_T) {
		sprintf(s1 + strlen(s1), "complementary design of number");
		sprintf(s2, "with respect to smaller t");
		}
	else if (rule == DP_RULE_COMPL_DERIVED) { 
		sprintf(s1 + strlen(s1), "derived from complement of number");
		sprintf(s2, "");
		}
	else if (rule == DP_RULE_COMPL_RESIDUAL) {
		sprintf(s1 + strlen(s1), "residual design of complement of number");
		sprintf(s2, "");
		}
	else if (rule == DP_RULE_COMPL_ALLTOP) {
		sprintf(s1 + strlen(s1), "Alltop construction for complementary design of number");
		sprintf(s2, "");
		}
	else {
		/* sprintf(s1 + strlen(s1), "no rule"); */
		sprintf(s2, "");
		}
	return OK;
}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::field_name(INT i, INT j, BYTE *str)
#endif
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "id"; break;
	case 1: s = "v"; break;
	case 2: s = "t"; break;
	case 3: s = "k"; break;
	case 4: s = "lambda"; break;
	case 5: s = "source"; break;
	default:
		return error("DESIGN_PARAMETER::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::init()
#endif
{
	INT erg = OK;
	
	erg += m_il(8);
	c_obj_k(DESIGN_PARAMETER_KIND);
	s_id()->m_i(-1);
	s_v()->m_i(-1);
	s_t()->m_i(-1);
	s_k()->m_i(-1);
	s_lambda()->m_i(-1);
	s_source()->m_il(0);
	return erg;
}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::sprint(BYTE *s)
#endif
{
	BYTE str[10000];
	BYTE str1[10000];
	INT i, l, lc;

	str[0] = 0;
	sprintf(str + strlen(str), "%ld: %ld-(%ld,%ld,%ld) ", 
		s_id_i(), s_t_i(), s_v_i(), s_k_i(), s_lambda_i());
	if (selfcomplementary()) {
		sprintf(str + strlen(str), "(selfcomplementary) ");
		}
	else {
		lc = lambda_compl();
		sprintf(str + strlen(str), "(complement has lambda %ld) ", lc);
		}
	l = s_source()->s_li();
	for (i = 0; i < l; i++) {
		str1[0] = 0;
		s_source_i(i)->sprint(str1);
		sprintf(str + strlen(str), "%s", str1);
		if (i < l - 1) {
			sprintf(str + strlen(str), ", ");
			}
		}

	//if (strlen(s) + strlen(str) < 200)
		strcat(s, str);
	return OK;
}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::selfcomplementary()
#endif
{
	DESIGN_PARAMETER_OB c;
	
	complement(&c);
	if (c.s_lambda_i() == s_lambda_i())
		return TRUE;
	else
		return FALSE;
}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::lambda_compl()
#endif
{
	DESIGN_PARAMETER_OB c;
	
	complement(&c);
	return c.s_lambda_i();
}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::complement(DESIGN_PARAMETER_OP p)
#endif
{
	INT i, id, v, t, k, lambda;
	INT nom, denom, n, d, g;
	INT lambda_new;
	DESIGN_PARAMETER_SOURCE_OB S;

	id = s_id_i();
	v = s_v_i();
	t = s_t_i();
	k = s_k_i();
	lambda = s_lambda_i();
	nom = 1;
	denom = 1;
	n = v - t;
	d = k - t;
	for (i = 0; i < k - t; i++) {
		nom *= n;
		denom *= d;
		n--;
		d--;
		ggt_iipi(nom, denom, &g);
		if (g != 1 && g != -1) {
			nom /= g;
			denom /= g;
			}
		}
	if (denom != 1)
		return error("DP::complement() error: denom != 1");
	lambda_new = nom - lambda;
	if (lambda_new <= 0)
		return error("DP::complement() error: lambda_new <= 0");
	p->init();
	p->s_v()->m_i(v);
	p->s_t()->m_i(t);
	p->s_k()->m_i(k);
	p->s_lambda()->m_i(lambda_new);
	S.init();
	S.s_prev()->m_i(id);
	S.s_rule()->m_i(DP_RULE_COMPLEMENT);
	p->s_source()->append_in_place(&S);
	return OK;
}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::reduce_t(DESIGN_PARAMETER_OP p)
#endif
{
	INT id, v, t, k, lambda;
	INT lambda_new;
	DESIGN_PARAMETER_SOURCE_OB S;

	id = s_id_i();
	v = s_v_i();
	t = s_t_i();
	k = s_k_i();
	lambda = s_lambda_i();
	lambda_new = (lambda * (v - t + 1)) / (k - t + 1);
	p->init();
	p->s_v()->m_i(v);
	p->s_t()->m_i(t - 1);
	p->s_k()->m_i(k);
	p->s_lambda()->m_i(lambda_new);
	S.init();
	S.s_prev()->m_i(id);
	S.s_rule()->m_i(DP_RULE_REDUCED_T);
	p->s_source()->append_in_place(&S);
	return OK;
}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::compl_reduce_t(DESIGN_PARAMETER_OP p)
#endif
{
	INT id, v, t, k, lambda;
	INT lambda_new;
	DESIGN_PARAMETER_SOURCE_OB S;
	DESIGN_PARAMETER_OB c;

	complement(&c);
	id = s_id_i();
	v = c.s_v_i();
	t = c.s_t_i();
	k = c.s_k_i();
	lambda = c.s_lambda_i();
	lambda_new = (lambda * (v - t + 1)) / (k - t + 1);
	p->init();
	p->s_v()->m_i(v);
	p->s_t()->m_i(t - 1);
	p->s_k()->m_i(k);
	p->s_lambda()->m_i(lambda_new);
	S.init();
	S.s_prev()->m_i(id);
	S.s_rule()->m_i(DP_RULE_COMPL_REDUCED_T);
	p->s_source()->append_in_place(&S);
	return OK;
}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::derive(DESIGN_PARAMETER_OP p)
#endif
{
	INT id, v, t, k, lambda;
	DESIGN_PARAMETER_SOURCE_OB S;

	id = s_id_i();
	v = s_v_i();
	t = s_t_i();
	k = s_k_i();
	lambda = s_lambda_i();
	p->init();
	p->s_v()->m_i(v - 1);
	p->s_t()->m_i(t - 1);
	p->s_k()->m_i(k - 1);
	p->s_lambda()->m_i(lambda);
	S.init();
	S.s_prev()->m_i(id);
	S.s_rule()->m_i(DP_RULE_DERIVED);
	p->s_source()->append_in_place(&S);
	return OK;
}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::compl_derive(DESIGN_PARAMETER_OP p)
#endif
{
	INT id, v, t, k, lambda;
	DESIGN_PARAMETER_SOURCE_OB S;
	DESIGN_PARAMETER_OB c;

	complement(&c);
	id = s_id_i();
	v = c.s_v_i();
	t = c.s_t_i();
	k = c.s_k_i();
	lambda = c.s_lambda_i();
	p->init();
	p->s_v()->m_i(v - 1);
	p->s_t()->m_i(t - 1);
	p->s_k()->m_i(k - 1);
	p->s_lambda()->m_i(lambda);
	S.init();
	S.s_prev()->m_i(id);
	S.s_rule()->m_i(DP_RULE_COMPL_DERIVED);
	p->s_source()->append_in_place(&S);
	return OK;
}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::residual(DESIGN_PARAMETER_OP p)
#endif
{
	INT id, v, t, k, lambda;
	INT lambda_new;
	DESIGN_PARAMETER_SOURCE_OB S;

	id = s_id_i();
	v = s_v_i();
	t = s_t_i();
	k = s_k_i();
	lambda = s_lambda_i();
	lambda_new = (lambda * (v - t + 1)) / (k - t + 1) - lambda;
	p->init();
	p->s_v()->m_i(v - 1);
	p->s_t()->m_i(t - 1);
	p->s_k()->m_i(k);
	p->s_lambda()->m_i(lambda_new);
	S.init();
	S.s_prev()->m_i(id);
	S.s_rule()->m_i(DP_RULE_RESIDUAL);
	p->s_source()->append_in_place(&S);
	return OK;
}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::compl_residual(DESIGN_PARAMETER_OP p)
#endif
{
	INT id, v, t, k, lambda;
	INT lambda_new;
	DESIGN_PARAMETER_SOURCE_OB S;
	DESIGN_PARAMETER_OB c;

	complement(&c);
	id = s_id_i();
	v = c.s_v_i();
	t = c.s_t_i();
	k = c.s_k_i();
	lambda = c.s_lambda_i();
	lambda_new = (lambda * (v - t + 1)) / (k - t + 1) - lambda;
	p->init();
	p->s_v()->m_i(v - 1);
	p->s_t()->m_i(t - 1);
	p->s_k()->m_i(k);
	p->s_lambda()->m_i(lambda_new);
	S.init();
	S.s_prev()->m_i(id);
	S.s_rule()->m_i(DP_RULE_COMPL_RESIDUAL);
	p->s_source()->append_in_place(&S);
	return OK;
}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::trung_supplementary(DESIGN_PARAMETER_OP p)
#endif
{
	INT id, v, t, k, lambda;
	INT lambda_new;
	DESIGN_PARAMETER_SOURCE_OB S;

	id = s_id_i();
	v = s_v_i();
	t = s_t_i();
	k = s_k_i();
	lambda = s_lambda_i();
	if (v == 2 * k + 1) {
		// lambda_new = (lambda * (v - t)) / (k - t);
		lambda_new = (lambda * (2 * k + 2 - t)) / (k + 1 - t);
		p->init();
		p->s_v()->m_i(v + 1);
		p->s_t()->m_i(t);
		p->s_k()->m_i(k + 1);
		p->s_lambda()->m_i(lambda_new);
		S.init();
		S.s_prev()->m_i(id);
		S.s_rule()->m_i(DP_RULE_TRUNG_SUPPLEMENTARY);
		p->s_source()->append_in_place(&S);
		return TRUE;
		}
	return FALSE;
}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::supplementary(DESIGN_PARAMETER_OP p)
#endif
{
	INT id, v, t, k, lambda;
	DESIGN_PARAMETER_SOURCE_OB S;
	SYM_OB lambda_new;

	id = s_id_i();
	v = s_v_i();
	t = s_t_i();
	k = s_k_i();
	lambda = s_lambda_i();
	calc_lambda_ijs(t, v, k, lambda, 1 /* s */, 0 /* i */, t /* j */, &lambda_new);
	p->init();
	p->s_v()->m_i(v);
	p->s_t()->m_i(t);
	p->s_k()->m_i(v - k);
	lambda_new.copy(p->s_lambda());
	S.init();
	S.s_prev()->m_i(id);
	S.s_rule()->m_i(DP_RULE_SUPPLEMENTARY);
	p->s_source()->append_in_place(&S);
	return TRUE;

}

#if TEXDOCU
INT DESIGN_PARAMETER_OB::alltop(DESIGN_PARAMETER_OP p)
#else
Alltop~\cite{Alltop75}.
#endif
{
	INT id, v, t, k, lambda;
	DESIGN_PARAMETER_SOURCE_OB S;
	SYM_OB l, lmax, l2, two;

	id = s_id_i();
	v = s_v_i();
	t = s_t_i();
	k = s_k_i();
	lambda = s_lambda_i();
	if (v == 2 * k + 1 && EVEN(t)) {
		p->init();
		p->s_v()->m_i(v + 1);
		p->s_t()->m_i(t + 1);
		p->s_k()->m_i(k + 1);
		p->s_lambda()->m_i(lambda);
		S.init();
		S.s_prev()->m_i(id);
		S.s_rule()->m_i(DP_RULE_ALLTOP);
		p->s_source()->append_in_place(&S);
		return TRUE;
		}
	if (v == 2 * k + 1 && ODD(t)) {
		Binomial(v - t, k - t, &lmax);
		two.m_i_i(2);
		lmax.ganzdiv(&two, &l2);
		l.m_i_i(lambda);
		if (l.eq(&l2)) {
			p->init();
			p->s_v()->m_i(v + 1);
			p->s_t()->m_i(t + 1);
			p->s_k()->m_i(k + 1);
			p->s_lambda()->m_i(lambda);
			S.init();
			S.s_prev()->m_i(id);
			S.s_rule()->m_i(DP_RULE_ALLTOP);
			p->s_source()->append_in_place(&S);
			return TRUE;
			}
		}
	return FALSE;
}

#if TEXDOCU
INT design_parameter_produce(VECTOR_OP V, 
	INT t, INT v, INT k, INT lambda, BYTE *comment, INT f_v)
#else
Produces design parameter sets out of the given $t$-$(v,k,\lambda)$ one. 
The parameter sets are stored in a vector V.
The vector V is not sorted, so a linear search must be applied.
f\_v TRUE gives text output during processing.
#endif
{
	DESIGN_PARAMETER_OB p, q;
	DESIGN_PARAMETER_OP pp;
	INT max_len = DP_ALLOC_LEN;
	INT i, j, old_id, id = 0;

	p.init();
	p.s_id()->m_i(id);
	p.s_v()->m_i(v);
	p.s_t()->m_i(t);
	p.s_k()->m_i(k);
	p.s_lambda()->m_i(lambda);
#if 0
	// p.s_prev()->m_i(-1);
	// p.s_rule()->m_i(-1);
	if (comment)
		p.s_comment()->init(comment);
#endif
	V->m_il(max_len);
	p.swap(V->s_i(0));
	id++;
	i = 0;
	j = 1;
	while (i < j) {
		if (j > MAX_DP) {
			V->realloc_z(j);
			return OK;
			}
		if (j + 5 > max_len) {
			max_len += DP_ALLOC_LEN;
			V->realloc_z(max_len);
			}
		pp = (DESIGN_PARAMETER_OP) V->s_i(i);
		old_id = pp->s_id_i();

#if 0
		pp->complement(&q);
		q.s_id()->m_i(id);
		q.s_prev()->m_i(old_id);
		if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
			if (dp_search(V, j, &q) == -1) {
				q.swap(V->s_i(j));
				j++;
				id++;
				}
			}
#endif

		pp->reduce_t(&q);
		q.s_id()->m_i(id);
		// q.s_prev()->m_i(old_id);
		if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
			if (dp_search(V, j, &q) == -1) {
				q.swap(V->s_i(j));
				j++;
				id++;
				}
			}

		pp->derive(&q);
		q.s_id()->m_i(id);
		// q.s_prev()->m_i(old_id);
		if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
			if (dp_search(V, j, &q) == -1) {
				q.swap(V->s_i(j));
				j++;
				id++;
				}
			}

		pp->residual(&q);
		q.s_id()->m_i(id);
		// q.s_prev()->m_i(old_id);
		if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
			if (dp_search(V, j, &q) == -1) {
				q.swap(V->s_i(j));
				j++;
				id++;
				}
			}

		// trung_supplementary missing
		
		i++;
		}
	V->realloc_z(j);
	return OK;
}

#if TEXDOCU
INT dp_search(VECTOR_OP V, INT len, DESIGN_PARAMETER_OP p)
#else
Searches the vector V for a specific parameter set.
Searches for v, t, k, lambda, 
but not for the id or so.
This is only a linear search. 
#endif
{
	DESIGN_PARAMETER_OP pp;
	INT i, v, t, k, lambda;

	v = p->s_v_i();
	t = p->s_t_i();
	k = p->s_k_i();
	lambda = p->s_lambda_i();
	for (i = 0; i < len; i++) {
		pp = (DESIGN_PARAMETER_OP) V->s_i(i);
		if (pp->s_v_i() != v)
			continue; 
		if (pp->s_t_i() != t)
			continue; 
		if (pp->s_k_i() != k)
			continue; 
		if (pp->s_lambda_i() != lambda)
			continue; 
		return i;
		}
	return -1;
}



#if TEXDOCU
#else
This file implements a database of design parameter sets. 
A database of solution vectors is also implemented.
#endif

#undef DEBUG_BT

#undef DP_DEBUG_LOAD_ID
#undef DP_DEBUG_LOAD_ID_VERBOSE
#undef DP_DEBUG_SEARCH
#undef DP_DEBUG_SEARCH_VERBOSE


#if TEXDOCU
INT i_bitvec_db(DATABASE_OP *db, BYTE *db_prefix)
#else
bitvec\_db is a database of design 0,1 solution vectors. 
Here, we create a database object for the database. 
database objects are used for accessing the elements in the 
database or for storing new elements in the database.
#endif
{
	DATABASE_OP db1;
	BAYERTREE_OP bt;
	BYTE fn_db[1024];
	BYTE fn_idx[1024];

	sprintf(fn_db, "%sbitvec.db", db_prefix);
	sprintf(fn_idx, "%sbitvec0.idx", db_prefix);

	db1 = (DATABASE_OP) callocobject("dp1.C: i_bitvec_db");
	db1->init(fn_db, BITVEC_KIND);
	db1->add_btree(fn_idx, TRUE /* f_duplicatekeys */, 0 /* btree_idx */ );
	bt = db1->s_btree_i(0);
	bt->add_key_INT4_VEC(0, 0, 
		0 /* vec_fst */, 12 /* 8 */ /*  vec_len */);

	*db = db1;
	return OK;
}

#if TEXDOCU
void e_bitvec_db(DATABASE_OP db)
#else
Free the object.
#endif
{
	freeall(db);
}

#if TEXDOCU
INT do_db_bitvec_create_and_close(BYTE *db_prefix)
#else
create the database, close is (an empty database). 
#endif
{
	DATABASE_OP db;
	
	if (i_bitvec_db(&db, db_prefix) != OK) {
		return error("do_db_bitvec_create_and_close() i_bitvec_db");
		}
	if (db->create() != OK) {
		return error("do_db_bitvec_create_and_close() DB::create");
		}
	if (db->close() != OK) {
		return error("do_db_bitvec_create_and_close() DB::close");
		}
	e_bitvec_db(db);
	return OK;
}

#if TEXDOCU
INT read_01_file(BYTE *db_prefix, BYTE *fname)
#else
Reads the solution file \lq fname\rq and creates and fills a database. 
The database has \lq db\_prefix\rq as ists prefix in all its filenames.
#endif
{
	DATABASE_OP db;
	BAYERTREE_OP bt;
	BYTE str[BUFSIZE];
	FILE *fp;
	INT l;
	BITVEC_OB B;
		
	sprintf(str, "%sbitvec.db", db_prefix);
	if (file_size(str) <= 0) {
		do_db_bitvec_create_and_close(db_prefix);
		}
	if (i_bitvec_db(&db, db_prefix) != OK) {
		return error("read_01_file() i_bitvec_db");
		}
	bt = db->s_btree_i(0);
	if (db->open() != OK) {
		return error("read_01_file() DB::open");
		}
	fp = fopen(fname, "r");
	while (TRUE) {
		if (fgets(str, BUFSIZE, fp) == NULL)
			break;
		l = strlen(str);
		str[l - 1] = '0';
		printf("%s\n", str);
		B.init_ascii(str);
		B.println();
		if (db->add_op(&B) != OK) {
			return error("read_01_file() DB::add_op");
			}
		
		}
	fclose(fp);
	bt->len(&l);
	printf("len = %ld\n", l);
	fflush(stdout);
	if (db->close() != OK) {
		return error("read_01_file() DB::close");
		}
	
	e_bitvec_db(db);
	return OK;
}

/* 
 * btree #0: id, v, t, k, lambda
 * btree #1: v, t, k, lambda, id
 * btree #2: t, v, k, lambda, id
 * btree #3: lambda, v, t, k, id
 */

#if TEXDOCU
INT i_db_dp(DATABASE_OP *p_db, BYTE *path)
#else
Creates a database object for the database of design parameter sets.
The design parameter set database has 4 btree access paths: \\
btree 0 is access via id, v, t, k, lambda \\
btree 1 is access via v, t, k, lambda, id \\
btree 2 is access via t, v, k, lambda, id \\
btree 3 is access via lambda, v, t, k, id \\
path without trailing slash.
#endif
{
	DATABASE_OP db;
	BAYERTREE_OP bt;
	BYTE file_db[1024];
	BYTE file_idx0[1024];
	BYTE file_idx1[1024];
	BYTE file_idx2[1024];
	BYTE file_idx3[1024];

	strcpy(file_db, path);
	strcat(file_db, "/dp.db");
	strcpy(file_idx0, path);
	strcat(file_idx0, "/dp0.idx");
	strcpy(file_idx1, path);
	strcat(file_idx1, "/dp1.idx");
	strcpy(file_idx2, path);
	strcat(file_idx2, "/dp2.idx");
	strcpy(file_idx3, path);
	strcat(file_idx3, "/dp3.idx");
	db = (DATABASE_OP) callocobject("DATABASE_OB i_db_dp()");
	db->init(file_db, DESIGN_PARAMETER_KIND);

	db->add_btree(file_idx0, 
		TRUE /* f_duplicatekeys */, 
		0 /* btree_idx */ );
	bt = db->s_btree_i(0);
	bt->add_key_INT4(0 /* id */, 0);
	bt->add_key_INT4(1 /* v */, 0);
	bt->add_key_INT4(2 /* t */, 0);
	bt->add_key_INT4(3 /* k */, 0);
	bt->add_key_INT4(4 /* lambda */, 0);

	db->add_btree(file_idx1, 
		TRUE /* f_duplicatekeys */, 
		1 /* btree_idx */ );
	bt = db->s_btree_i(1);
	bt->add_key_INT4(1 /* v */, 0);
	bt->add_key_INT4(2 /* t */, 0);
	bt->add_key_INT4(3 /* k */, 0);
	bt->add_key_INT4(4 /* lambda */, 0);
	bt->add_key_INT4(0 /* id */, 0);

	db->add_btree(file_idx2, 
		TRUE /* f_duplicatekeys */, 
		2 /* btree_idx */ );
	bt = db->s_btree_i(2);
	bt->add_key_INT4(2 /* t */, 0);
	bt->add_key_INT4(1 /* v */, 0);
	bt->add_key_INT4(3 /* k */, 0);
	bt->add_key_INT4(4 /* lambda */, 0);
	bt->add_key_INT4(0 /* id */, 0);

	db->add_btree(file_idx3, 
		TRUE /* f_duplicatekeys */, 
		3 /* btree_idx */ );
	bt = db->s_btree_i(3);
	bt->add_key_INT4(4 /* lambda */, 0);
	bt->add_key_INT4(1 /* v */, 0);
	bt->add_key_INT4(2 /* t */, 0);
	bt->add_key_INT4(3 /* k */, 0);
	bt->add_key_INT4(0 /* id */, 0);

	*p_db = db;
	return OK;
}

#if TEXDOCU
void e_db_dp(DATABASE_OP db)
#else
Frees the object.
#endif
{
	freeall(db);
}

#if TEXDOCU
INT do_db_dp_create_and_close(BYTE *path)
#else
Creates an empty database and closes it.
#endif
{
	DATABASE_OP db;
	
	if (i_db_dp(&db, path) != OK) {
		return error("do_db_dp_create_and_close() i_db_dp");
		}
	if (db->create() != OK) {
		return error("do_db_dp_create() DB::create");
		}
	if (db->close() != OK) {
		return error("do_db_dp_create() DB::close");
		}
	e_db_dp(db);
	return OK;
}

#if TEXDOCU
INT do_db_dp_info(BYTE *path, INT btree_idx)
#else
Calls the btree function \lq print\_pages\rq for the btree index \lq btree\_idx\rq.
#endif
{
	BAYERTREE_OP btree;
	DATABASE_OP db;
	
	if (i_db_dp(&db, path) != OK) {
		return error("do_db_dp_info() i_db_dp");
		}
	/* da_print(db, 0); */
	btree = db->s_btree_i(btree_idx);
	btree->print_pages();
	e_db_dp(db);
	return OK;
}

#if TEXDOCU
INT db_dp_highest_id(DATABASE_OP db, INT *highest_id)
#else
returns the highest id which is used in the database.
#endif
{
	BAYERTREE_OP btree;
	INT btree_idx;
	KEYTYPE key;
	INT4 *pi;
	INT id, len;
	DATATYPE data;

	btree_idx = 0;
	btree = db->s_btree_i(btree_idx);
	if (btree->open() != OK) {
		return error("db_dp_highest_id() BT::open");
		}
	if (btree->len(&len) != OK) {
		return error("db_dp_highest_id() BT::len");
		}
	if (len == 0) {
		btree->close();
		*highest_id = -1;
		return OK;
		}
	if (btree->ith(len - 1, &key, &data) != OK) {
		return error("db_dp_highest_id() BT::ith");
		}
	pi = (INT4 *) &key;
	id = pi[0];
	*highest_id = id;
	if (btree->close() != OK) {
		return error("db_dp_highest_id() BT::close");
		}
	return OK;
}

#if TEXDOCU
INT db_dp_load_id_data(DATABASE_OP db, INT id, DESIGN_PARAMETER_OP p, DATATYPE *data)
#else
Loads the database parameter set by id.
#endif
{
	BAYERTREE_OP btree;
	INT btree_idx;
	KEYTYPE key, key1;
	INT4 *pi, *pi1;
	INT idx, f_found;

#ifdef DP_DEBUG_LOAD_ID
	printf("db_dp_load_id_data() loading id = %ld\n", id);
#endif
	btree_idx = 0; /* id, v, t, k, lambda */
	btree = db->s_btree_i(btree_idx);
	if (btree->open() != OK) {
		return error("db_dp_load_id_data() BT::open");
		}
	pi = (INT4 *) &key;
	pi1 = (INT4 *) &key1;
	pi[0] = id + 1;
	pi[1] = 0;
	pi[2] = 0;
	pi[3] = 0;
	pi[4] = 0;
	if (btree->search(pi /* pSearchKey */, 
		data, &idx, &f_found) != OK) {
		return error("db_dp_load_id_data() BT::search");
		}
#ifdef DP_DEBUG_LOAD_ID_VERBOSE
	printf("db_dp_load_id_data(): f_found = %ld %ld "
		"idx = %ld\n", 
		f_found, (INT) pi[0], 
		idx);
#endif
	if (f_found)
		return error("db_dp_load_id_data(): found id + 1, 0, 0, 0, 0");

	if (db->open_DB() != OK) {
		return error("db_dp_load_id_data() DB::open_DB");
		}
	if (btree->ith(idx, &key1, data) != OK) {
		return error("db_dp_load_id_data() BT::ith");
		}
	if (pi1[0] != id) {
		return error("db_dp_load_id_data() pi1[0] != id (id not found)");
		}
	if (db->get_op(data, p) != OK) {
		return error("db_dp_load_id_data() DB::get_op");
		}
	if (db->close_DB() != OK) {
		return error("db_dp_load_id_data() DB::close_DB");
		}
	if (btree->close() != OK) {
		return error("db_dp_load_id_data() BT::close");
		}
	return OK;
}

#if TEXDOCU
INT db_dp_load_id(DATABASE_OP db, INT id, DESIGN_PARAMETER_OP p)
#else
Loads the database parameter set by id.
#endif
{
	DATATYPE data;

	return db_dp_load_id_data(db, id, p, &data);
}

#if TEXDOCU
INT db_dp_search(DATABASE_OP db, 
	INT t, INT v, INT k, INT lambda, 
	INT *idx, INT *len, INT *f_found)
#else
Searches for a database parameter set $t$-$(v,k,lambda)$.
Returns idx / len /f\_found.
#endif
{
	BAYERTREE_OP btree;
	INT btree_idx;
	KEYTYPE key;
	INT4 *pi;
	INT idx_first, idx_last, len1, f_found1;
	DATATYPE data;

#ifdef DP_DEBUG_SEARCH
	printf("in db_dp_search() %ld-(%ld,%ld,%ld)\n", t, v, k, lambda);
	fflush(stdout);
#endif
	btree_idx = 2; /* t, v, k, lambda, id */
	btree = db->s_btree_i(btree_idx);
	if (btree->open() != OK) {
		return error("db_dp_search() BT::open");
		}
	pi = (INT4 *) &key;
	pi[0] = t;
	pi[1] = v;
	pi[2] = k;
	pi[3] = lambda + 1;
	pi[4] = -1; /* id */
	if (btree->search(pi /* pSearchKey */, 
		&data, &idx_last, &f_found1) != OK) {
		return error("db_dp_search() BT::search");
		}
#ifdef DP_DEBUG_SEARCH_VERBOSE
	printf("  f_found = %ld %ld-(%ld,%ld,%ld) id = %ld idx_last = %ld\n", 
		f_found1, (INT) pi[0], (INT) pi[1], (INT) pi[2], (INT) pi[3], (INT) pi[4], 
		idx_last);
#endif
	if (f_found1)
		return error("db_dp_search(): id = 0 found");

	pi[3] = lambda;
	if (btree->search(pi /* pSearchKey */, 
		&data, &idx_first, &f_found1) != OK) {
		return error("db_dp_search() BT::search");
		}
#ifdef DP_DEBUG_SEARCH_VERBOSE
	printf("  f_found = %ld %ld-(%ld,%ld,%ld) id = %ld idx_first = %ld\n", 
		f_found1, (INT) pi[0], (INT) pi[1], (INT) pi[2], (INT) pi[3], (INT) pi[4], 
		idx_first);
#endif
	if (f_found1)
		return error("db_dp_search(): id = 0 found");
	idx_first++;

	len1 = idx_last + 1 - idx_first;
	if (len1 > 0)
		*f_found = TRUE;
	else
		*f_found = FALSE;
	*idx = idx_first;
	*len = len1;
#ifdef DP_DEBUG_SEARCH
	printf("in db_dp_search() f_found = %ld\n", *f_found);
	fflush(stdout);
#endif
	if (btree->close() != OK) {
		return error("db_dp_search() BT::close");
		}
	return OK;
}

#if TEXDOCU
INT db_dp_parameter_add(DATABASE_OP db, 
	DESIGN_PARAMETER_OP p, INT *id, INT f_verbose, 
	INT *f_added, FILE *fp_tex)
#else
Adds the design parameter set object p if the parameter set is new, i.e. 
not already contained in the database (nor the complement). 
f\_found and f\_found\_complement are set to TRUE / FALSE whether or not 
the design or its complement have been found in the database.
#endif
{
	DESIGN_PARAMETER_OB p_compl;
	SYM_OB lmax, two, l2;
	INT idx, len;
	INT f_found;
	BYTE s2[10000];

	*f_added = FALSE;
	
	if (p->s_lambda()->s_obj_k() == LONGINT) {
		return OK;
		}
	Binomial(p->s_v_i() - p->s_t_i(), p->s_k_i() - p->s_t_i(), &lmax);
	two.m_i_i(2);
	lmax.ganzdiv(&two, &l2);
	if (p->s_lambda()->gt(&l2)) {
		*f_added = FALSE;
		return OK;
		}
	
	if (db_dp_search(db, p->s_t_i(), p->s_v_i(), p->s_k_i(), p->s_lambda_i(), 
		&idx, &len, &f_found) != OK)
		return error("db_dp_parameter_add_if_new() error in db_dp_search()");
	if (!f_found) {

#if 0
		if (p->s_v_i() == 2 * p->s_k_i() + 1 || p->s_v_i() == 2 * p->s_k_i() - 1) {
			;
			}
		else {
			p->complement(&p_compl);
			if (db_dp_search(db, p_compl.s_t_i(), p_compl.s_v_i(), 
				p_compl.s_k_i(), p_compl.s_lambda_i(), 
				&idx, &len, f_found_complement) != OK)
				return error("db_dp_parameter_add_if_new() error in db_dp_search(compl)");
			if (*f_found_complement) {
				if (f_verbose) {
					printf("found the complement %ld-(%ld,%ld,%ld), not added !\n", 
						p_compl.s_t_i(), p_compl.s_v_i(), 
						p_compl.s_k_i(), p_compl.s_lambda_i());
					fflush(stdout);
					}
				return OK;
				}
			}
#endif


		if (db->open() != OK) {
			return error("db_dp_parameter_add_if_new() error in DB::open()");
			}
		p->s_id()->m_i(*id);
		if (fp_tex) {
			fprintf(fp_tex, "%ld: & \\mbox{$%ld$-$(%ld,%ld,%ld)$} & \\mbox{", 
				p->s_id_i(), p->s_t_i(), p->s_v_i(), p->s_k_i(), p->s_lambda_i());
#if 0
			p->sprint_text012(s0, s1, s2);
			fprintf(fp_tex, "%s ", s1);
			if (p->s_prev_i() >= 0)
				fprintf(fp_tex, "%ld", p->s_prev_i());
#endif
			fprintf(fp_tex, " %s} \\\\\n", s2);
			}
		(*id)++;
		if (f_verbose) {
				printf("vor db->add_op() ");
				p->println();
				fflush(stdout);
				}
		if (db->add_op(p) != OK) {
			db->close();
			return error("db_dp_parameter_add_if_new() error in DB::add_op");
			}
		if (db->close() != OK) {
			return error("db_dp_parameter_add_if_new() error in DB::close()");
			}
		*f_added = TRUE;
		}
	else {
		INT btree_idx;
		KEYTYPE key;
		DATATYPE data;
		DESIGN_PARAMETER_OB p1;
		VECTOR_OB V;
		
		btree_idx = 2; // t, v, k, lambda
		if (db->open() != OK) {
			return error("db_dp_parameter_add_if_new() error in DB::open()");
			}
		db->ith(btree_idx, idx, &key, &data);
		db->get_op(&data, &p1);
		p1.s_source()->append(p->s_source(), &V);
		V.swap(p1.s_source());
		db->del_op(&p1, data.datref, data.data_size);
		if (db->add_op(&p1) != OK) {
			db->close();
			return error("db_dp_parameter_add_if_new() error in DB::add_op");
			}
		if (db->close() != OK) {
			return error("db_dp_parameter_add_if_new() error in DB::close()");
			}
		}
	return OK;
}

#if TEXDOCU
static INT dp_add(VECTOR_OP dp_data, INTEGER_OP dp_data_length, 
	INT t, INT v, INT k, INT lambda, INT id)
#endif
{
	VECTOR_OB dp;
	
	dp.m_il(5);
	dp.m_ii(0, t);
	dp.m_ii(1, v);
	dp.m_ii(2, k);
	dp.m_ii(3, lambda);
	dp.m_ii(4, id);
	dp_data->append_element(dp_data_length, &dp);
	return OK;
}

#if TEXDOCU
INT db_dp_design_parameter_produce(DATABASE_OP db, 
	INT t, INT v, INT k, INT lambda, BYTE *comment, INT f_v, FILE *fp_tex, 
	MATRIX_OP data)
#endif
{
	DESIGN_PARAMETER_OB p, q, pp, pp_compl;
	INT i, old_id, id = 0;
	INT f_added;
	VECTOR_OB dp_data;
	VECTOR_OP dp;
	INTEGER_OB dp_data_length;
	INT l;
	INT t_, v_, k_, lambda_, id_;

	data->m_ilih(0, 0);
	dp_data.m_il(VECTOR_OVERSIZE);
	dp_data_length.m_i(0);
	
	db_dp_highest_id(db, &id);
	id++;
	p.init();
	p.s_id()->m_i(id);
	p.s_v()->m_i(v);
	p.s_t()->m_i(t);
	p.s_k()->m_i(k);
	p.s_lambda()->m_i(lambda);
	// p.s_prev()->m_i(-1);
	// p.s_rule()->m_i(-1);
	if (comment) {
		DESIGN_PARAMETER_SOURCE_OB S;
		
		S.init();
		S.s_comment()->init(comment);
		p.s_source()->append_in_place(&S);
		}
	
	i = id;
	if (db_dp_parameter_add(db, 
		&p, &id, f_v /* f_verbose */ , 
		&f_added, fp_tex) != OK)
		return error("db_dp_design_parameter_produce() "
			"error in db_dp_parameter_add_if_new()");
	if (id == i)
		return OK;
	// now id = i + 1 !
	dp_add(&dp_data, &dp_data_length, t, v, k, lambda, id - 1 /* the old id */);
	
	while (i < id) {
		printf("db_dp_design_parameter_produce() i = %ld id = %ld\n", i, id);
		if (db_dp_load_id(db, i, &pp) != OK)
			return error("db_dp_design_parameter_produce() "
				"error in db_dp_load_id()");
		old_id = pp.s_id_i();
		if (old_id != i)
			return error("db_dp_design_parameter_produce() old_id != i");

		if (pp.s_v_i() == 2 * pp.s_k_i() + 1) {
			pp.supplementary(&q);
			q.s_id()->m_i(id);
			// q.s_prev()->m_i(old_id);
			if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
				if (db_dp_parameter_add(db, 
					&q, &id, f_v /* f_verbose */ , 
					&f_added, fp_tex) != OK)
					return error("db_dp_design_parameter_produce() "
						"error in db_dp_parameter_add_if_new()");
				if (f_added)
					dp_add(&dp_data, &dp_data_length, 
						q.s_t_i(), q.s_v_i(), q.s_k_i(), q.s_lambda_i(), q.s_id_i());
				}
			}
		
		pp.reduce_t(&q);
		q.s_id()->m_i(id);
		// q.s_prev()->m_i(old_id);
		if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
			if (db_dp_parameter_add(db, 
				&q, &id, f_v /* f_verbose */ , 
				&f_added, fp_tex) != OK)
				return error("db_dp_design_parameter_produce() "
					"error in db_dp_parameter_add_if_new()");
			if (f_added)
				dp_add(&dp_data, &dp_data_length, 
					q.s_t_i(), q.s_v_i(), q.s_k_i(), q.s_lambda_i(), q.s_id_i());
			}

		pp.derive(&q);
		q.s_id()->m_i(id);
		// q.s_prev()->m_i(old_id);
		if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
			if (db_dp_parameter_add(db, 
				&q, &id, f_v /* f_verbose */ , 
				&f_added, fp_tex) != OK)
				return error("db_dp_design_parameter_produce() "
					"error in db_dp_parameter_add_if_new()");
			if (f_added)
				dp_add(&dp_data, &dp_data_length, 
					q.s_t_i(), q.s_v_i(), q.s_k_i(), q.s_lambda_i(), q.s_id_i());
			}

		pp.residual(&q);
		q.s_id()->m_i(id);
		// q.s_prev()->m_i(old_id);
		if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
			if (db_dp_parameter_add(db, 
				&q, &id, f_v /* f_verbose */ , 
				&f_added, fp_tex) != OK)
				return error("db_dp_design_parameter_produce() "
					"error in db_dp_parameter_add_if_new()");
			if (f_added)
				dp_add(&dp_data, &dp_data_length, 
					q.s_t_i(), q.s_v_i(), q.s_k_i(), q.s_lambda_i(), q.s_id_i());
			}

		if (pp.trung_supplementary(&q)) {
			q.s_id()->m_i(id);
			// q.s_prev()->m_i(old_id);
			if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
				if (db_dp_parameter_add(db, 
					&q, &id, f_v /* f_verbose */ , 
					&f_added, fp_tex) != OK)
					return error("db_dp_design_parameter_produce() "
						"error in db_dp_parameter_add_if_new()");
				if (f_added)
					dp_add(&dp_data, &dp_data_length, 
						q.s_t_i(), q.s_v_i(), q.s_k_i(), q.s_lambda_i(), q.s_id_i());
				}
			}

		if (pp.alltop(&q)) {
			q.s_id()->m_i(id);
			// q.s_prev()->m_i(old_id);
			if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
				if (db_dp_parameter_add(db, 
					&q, &id, f_v /* f_verbose */ , 
					&f_added, fp_tex) != OK)
					return error("db_dp_design_parameter_produce() "
						"error in db_dp_parameter_add_if_new()");
				if (f_added)
					dp_add(&dp_data, &dp_data_length, 
						q.s_t_i(), q.s_v_i(), q.s_k_i(), q.s_lambda_i(), q.s_id_i());
				}
			}


		pp.compl_reduce_t(&q);
		q.s_id()->m_i(id);
		// q.s_prev()->m_i(old_id);
		if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
			if (db_dp_parameter_add(db, 
				&q, &id, f_v /* f_verbose */ , 
				&f_added, fp_tex) != OK)
				return error("db_dp_design_parameter_produce() "
					"error in db_dp_parameter_add_if_new()");
			if (f_added)
				dp_add(&dp_data, &dp_data_length, 
					q.s_t_i(), q.s_v_i(), q.s_k_i(), q.s_lambda_i(), q.s_id_i());
			}

		pp.compl_derive(&q);
		q.s_id()->m_i(id);
		// q.s_prev()->m_i(old_id);
		if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
			if (db_dp_parameter_add(db, 
				&q, &id, f_v /* f_verbose */ , 
				&f_added, fp_tex) != OK)
				return error("db_dp_design_parameter_produce() "
					"error in db_dp_parameter_add_if_new()");
			if (f_added)
				dp_add(&dp_data, &dp_data_length, 
					q.s_t_i(), q.s_v_i(), q.s_k_i(), q.s_lambda_i(), q.s_id_i());
			}

		pp.compl_residual(&q);
		q.s_id()->m_i(id);
		// q.s_prev()->m_i(old_id);
		if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
			if (db_dp_parameter_add(db, 
				&q, &id, f_v /* f_verbose */ , 
				&f_added, fp_tex) != OK)
				return error("db_dp_design_parameter_produce() "
					"error in db_dp_parameter_add_if_new()");
			if (f_added)
				dp_add(&dp_data, &dp_data_length, 
					q.s_t_i(), q.s_v_i(), q.s_k_i(), q.s_lambda_i(), q.s_id_i());
			}




#if 0
		/* apply all constructions also for the complementary design: */
		pp.complement(&pp_compl);

		pp_compl.reduce_t(&q);
		// q.s_rule()->m_i(DP_RULE_COMPL_REDUCED_T);
		q.s_id()->m_i(id);
		// q.s_prev()->m_i(old_id);
		if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
			if (db_dp_parameter_add_if_new(db, 
				&q, &id, f_v /* f_verbose */ , 
				&f_found, &f_found_complement, fp_tex) != OK)
				return error("db_dp_design_parameter_produce() "
					"error in db_dp_parameter_add_if_new()");
			if (!f_found && !f_found_complement)
				dp_add(&dp_data, &dp_data_length, 
					q.s_t_i(), q.s_v_i(), q.s_k_i(), q.s_lambda_i(), q.s_id_i());
			}

		pp_compl.derive(&q);
		// q.s_rule()->m_i(DP_RULE_COMPL_DERIVED);
		q.s_id()->m_i(id);
		// q.s_prev()->m_i(old_id);
		if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
			if (db_dp_parameter_add_if_new(db, 
				&q, &id, f_v /* f_verbose */ , 
				&f_found, &f_found_complement, fp_tex) != OK)
				return error("db_dp_design_parameter_produce() "
					"error in db_dp_parameter_add_if_new()");
			if (!f_found && !f_found_complement)
				dp_add(&dp_data, &dp_data_length, 
					q.s_t_i(), q.s_v_i(), q.s_k_i(), q.s_lambda_i(), q.s_id_i());
			}

		pp_compl.residual(&q);
		// q.s_rule()->m_i(DP_RULE_COMPL_RESIDUAL);
		q.s_id()->m_i(id);
		// q.s_prev()->m_i(old_id);
		if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
			if (db_dp_parameter_add_if_new(db, 
				&q, &id, f_v /* f_verbose */ , 
				&f_found, &f_found_complement, fp_tex) != OK)
				return error("db_dp_design_parameter_produce() "
					"error in db_dp_parameter_add_if_new()");
			if (!f_found && !f_found_complement)
				dp_add(&dp_data, &dp_data_length, 
					q.s_t_i(), q.s_v_i(), q.s_k_i(), q.s_lambda_i(), q.s_id_i());
			}

		pp_compl.alltop(&q);
		// q.s_rule()->m_i(DP_RULE_COMPL_RESIDUAL);
		q.s_id()->m_i(id);
		// q.s_prev()->m_i(old_id);
		if (q.s_t_i() >= DP_MIN_T && q.s_lambda_i() <= DP_MAX_LAMBDA) {
			if (db_dp_parameter_add_if_new(db, 
				&q, &id, f_v /* f_verbose */ , 
				&f_found, &f_found_complement, fp_tex) != OK)
				return error("db_dp_design_parameter_produce() "
					"error in db_dp_parameter_add_if_new()");
			if (!f_found && !f_found_complement)
				dp_add(&dp_data, &dp_data_length, 
					q.s_t_i(), q.s_v_i(), q.s_k_i(), q.s_lambda_i(), q.s_id_i());
			}
#endif

		
		i++;
		}

	l = dp_data_length.s_i();
	data->m_ilih(5, l);
	for (i = 0; i < l; i++) {
		dp = (VECTOR_OP) dp_data.s_i(i);
		t_ = dp->s_ii(0);
		v_ = dp->s_ii(1);
		k_ = dp->s_ii(2);
		lambda_ = dp->s_ii(3);
		id_ = dp->s_ii(4);
		data->m_iji(i, 0, t_);
		data->m_iji(i, 1, v_);
		data->m_iji(i, 2, k_);
		data->m_iji(i, 3, lambda_);
		data->m_iji(i, 4, id_);
		}
	return OK;
}


#if TEXDOCU
INT db_dp_read_design_txt(BYTE *fname_txt, BYTE *fname_tex, 
	BYTE *db_path, INT f_clear_always, VECTOR_OP vec_data)
#endif
{
	FILE *fp, *fp_tex;
	INT i, l, t, v, k, lambda, line = 0, f_has_comment;
	DATABASE_OP db;
	BYTE buf[BUFSIZE], *p, *comment;
	INT f_v = TRUE;
	MATRIX_OB data;
	INT vec_data_len;

	vec_data->m_il(0);
	vec_data_len = 0;

	fp = fopen(fname_txt, "r");
	fp_tex = fopen(fname_tex, "w");
	i_db_dp(&db, db_path);

	call_system("rm a");
	date_as_string(buf);
	fprintf(fp_tex, "%% created by db_dp_read_design_txt() at %s\n", buf);
	
	while (TRUE) {
		buf[0] = 0;
		line++;
		if (fgets(buf, BUFSIZE, fp) == NULL) {
			break;
			}
		l = strlen(buf);
		if (l) {
			if (buf[l - 1] == '\n')
				buf[l - 1] = 0;
			}
		printf("%ld: %s\n", line, buf);
		fflush(stdout);
		p = buf;
		
		if (!s_scan_int(&p, &t)) {
			printf("can't read t !\n");
			break;
			}
		if (t <= 0)
			break;
		if (!s_scan_int(&p, &v)) {
			printf("can't read t !\n");
			break;
			}
		if (!s_scan_int(&p, &k)) {
			printf("can't read t !\n");
			break;
			}
		if (!s_scan_int(&p, &lambda)) {
			printf("can't read t !\n");
			break;
			}
		f_has_comment = FALSE;
		for (i = 0; i < (INT) strlen(p); i++) {
			if (!isspace(p[i])) {
				f_has_comment = TRUE;
				break;
				}
			}
		printf("*** %ld-(%ld,%ld,%ld)", t, v, k, lambda);

		if (f_clear_always) {
			do_db_dp_create_and_close(db_path);
			}
	
		if (f_has_comment) {
			comment = p;
			printf(" %s:\n", comment);
			}
		else {
			comment = NIL;
			printf(":\n");
			}
		fprintf(fp_tex, "\\noindent\n{\\bf Parameter sets obtained from %ld-(%ld,%ld,%ld)}\n", 
			t, v, k, lambda);
		fprintf(fp_tex, "\\[\n\\begin{array}{rll}\n");
		db_dp_design_parameter_produce(db, t, v, k, lambda, comment, f_v, fp_tex, &data);
		fprintf(fp_tex, "\\end{array}\n\\]\n");
		if (data.s_hi() > 0) {
			vec_data->inc();
			((SYM_OP) &data)->copy(vec_data->s_i(vec_data_len));	
			vec_data_len++;
			printf("vec_data->s_li() = %ld vec_data_len = %ld\n", 
				vec_data->s_li(), vec_data_len);
			data.freeself();
			}
		}
	e_db_dp(db);
	fclose(fp);
	fclose(fp_tex);
	return OK;
}




INT db_dp_HTML(BYTE *path)
{
	BAYERTREE_OP btree0, btree1, btree2, btree3;
	DATABASE_OP db;
	
	if (i_db_dp(&db, path) != OK)
		return error("db_dp_HTML(): error in i_db_dp()");
	btree0 = db->s_btree_i(0);
	btree1 = db->s_btree_i(1);
	btree2 = db->s_btree_i(2);
	btree3 = db->s_btree_i(3);
	/* btree->print_pages(); */


	if (db->open() != OK)
		return error("db_dp_HTML(): error in db->open()");

	db_dp_HTML_d1(db, btree2);
	db_dp_HTML_d2(db, btree1);
	db_dp_HTML_d3(db, btree3);
	db_dp_HTML_d4(db, btree0);
	
	
	db->close();
	e_db_dp(db);
	return OK;
}

#define DB_DP_HTML_FNAME1 "d1.html"
#define DB_DP_HTML_FNAME2 "d2.html"
#define DB_DP_HTML_FNAME3 "d3.html"
#define DB_DP_HTML_FNAME4 "d4.html"

#define T_MIN_FOR_HTML 5

static INT db_dp_HTML_d1(DATABASE_OP db, BAYERTREE_OP btree)
{
	BYTE *fname = DB_DP_HTML_FNAME1;
	FILE *fp;
	KEYTYPE key;
	DATATYPE data;
	DESIGN_PARAMETER_OB dp;
	INT i, len1, old_t, t, old_v, v, k, lambda, id;

	fp = fopen(fname, "w");
	
	fputs("<html>\n", fp);
	fputs("<head>\n", fp);
	fputs("<title>\n", fp);
	fputs("Large Designs: t, v, k, lambda\n", fp);
	fputs("</title>\n", fp);
	fputs("</head>\n", fp);
	
	fputs("<body>\n", fp);
	fputs("<h1>\n", fp);
	fputs("Large Designs: t, v, k, lambda\n", fp);
	fputs("</h1>\n", fp);
	fputs("<hr>\n", fp);
	print_page_head(fp);
	
	btree->len(&len1);
	old_t = -1;
	old_v = -1;
	for (i = 0; i < len1; i++) {
		btree->ith(i, &key, &data);
		db->get_op(&data, &dp);
		t = dp.s_t_i();
		if (t < T_MIN_FOR_HTML)
			continue;
		v = dp.s_v_i();
		k = dp.s_k_i();
		lambda = dp.s_lambda_i();
		id = dp.s_id_i();
		if (t != old_t) {
			fprintf(fp, "<p>\n");
			fprintf(fp, "<hr>\n");
			fprintf(fp, "<h2> %ld-Designs:</h2>\n", t);
			old_v = -1;
			}
		old_t = t;
			
		if (old_v != -1 && v != old_v)
			fprintf(fp, "<p>\n");
		else
			fprintf(fp, "<br>\n");

		write_head_HTML(fp, &dp);
		write_source_HTML(fp, &dp);
		

		old_v = v;
		}
	fprintf(fp, "<p>\n");
	copyright(fp);
	fputs("</body>\n", fp);
	fputs("</html>\n\n", fp);
	fflush(fp);
	fclose(fp);
	
	printf("wrote file %s of size %ld\n", fname, file_size(fname));
	fflush(stdout);
	return OK;
}

static INT db_dp_HTML_d2(DATABASE_OP db, BAYERTREE_OP btree)
{
	BYTE *fname = DB_DP_HTML_FNAME2;
	FILE *fp;
	KEYTYPE key;
	DATATYPE data;
	DESIGN_PARAMETER_OB dp;
	INT i, len1, old_v, v, old_t, t;

	fp = fopen(fname, "w");
	
	fputs("<html>\n", fp);
	fputs("<head>\n", fp);
	fputs("<title>\n", fp);
	fputs("Large Designs: v, t, k, lambda\n", fp);
	fputs("</title>\n", fp);
	fputs("</head>\n", fp);
	
	fputs("<body>\n", fp);
	fputs("<h1>\n", fp);
	fputs("Large Designs: v, t, k, lambda\n", fp);
	fputs("</h1>\n", fp);
	fputs("<hr>\n", fp);
	print_page_head(fp);
	
	btree->len(&len1);
	old_v = -1;
	old_t = -1;
	for (i = 0; i < len1; i++) {
		btree->ith(i, &key, &data);
		db->get_op(&data, &dp);
		v = dp.s_v_i();
		t = dp.s_t_i();
		if (t < T_MIN_FOR_HTML)
			continue;
		if (v != old_v) {
			fprintf(fp, "<p>\n");
			fprintf(fp, "<hr>\n");
			fprintf(fp, "<h2> on %ld points:</h2>\n", v);
			old_t = -1;
			}
		old_v = v;

		if (old_t != -1 && t != old_t)
			fprintf(fp, "<p>\n");
		else
			fprintf(fp, "<br>\n");

		write_head_HTML(fp, &dp);
		write_source_HTML(fp, &dp);


		
		old_t = t;
		}
	fprintf(fp, "<p>\n");
	copyright(fp);
	fputs("</body>\n", fp);
	fputs("</html>\n\n", fp);
	fflush(fp);
	fclose(fp);
	printf("wrote file %s of size %ld\n", fname, file_size(fname));
	fflush(stdout);
	return OK;
}

static INT db_dp_HTML_d3(DATABASE_OP db, BAYERTREE_OP btree)
{
	BYTE *fname = DB_DP_HTML_FNAME3;
	FILE *fp;
	KEYTYPE key;
	DATATYPE data;
	DESIGN_PARAMETER_OB dp;
	INT i, len1, old_v, v, old_t, t;

	fp = fopen(fname, "w");
	
	fputs("<html>\n", fp);
	fputs("<head>\n", fp);
	fputs("<title>\n", fp);
	fputs("Large Designs: lambda, v, t, k\n", fp);
	fputs("</title>\n", fp);
	fputs("</head>\n", fp);
	
	fputs("<body>\n", fp);
	fputs("<h1>\n", fp);
	fputs("Large Designs: lambda, v, t, k\n", fp);
	fputs("</h1>\n", fp);
	fputs("<hr>\n", fp);
	print_page_head(fp);
	
	btree->len(&len1);
	old_v = -1;
	old_t = -1;
	for (i = 0; i < len1; i++) {
		btree->ith(i, &key, &data);
		db->get_op(&data, &dp);
		v = dp.s_v_i();
		t = dp.s_t_i();
		if (t < T_MIN_FOR_HTML)
			continue;
#if 0
		if (v != old_v) {
			fprintf(fp, "<p>\n");
			fprintf(fp, "<hr>\n");
			fprintf(fp, "<h2> on %ld points:</h2>\n", v);
			old_t = -1;
			}
		old_v = v;

		if (old_t != -1 && t != old_t)
			fprintf(fp, "<p>\n");
		else
			fprintf(fp, "<br>\n");
#endif
		fprintf(fp, "<br>\n");
		write_head_HTML(fp, &dp);
		write_source_HTML(fp, &dp);


		
		// old_t = t;
		}
	fprintf(fp, "<p>\n");
	copyright(fp);
	fputs("</body>\n", fp);
	fputs("</html>\n\n", fp);
	fflush(fp);
	fclose(fp);
	printf("wrote file %s of size %ld\n", fname, file_size(fname));
	fflush(stdout);
	return OK;
}

static INT db_dp_HTML_d4(DATABASE_OP db, BAYERTREE_OP btree)
{
	BYTE *fname = DB_DP_HTML_FNAME4;
	FILE *fp;
	KEYTYPE key;
	DATATYPE data;
	DESIGN_PARAMETER_OB dp;
	INT i, len1, id;

	fp = fopen(fname, "w");
	
	fputs("<html>\n", fp);
	fputs("<head>\n", fp);
	fputs("<title>\n", fp);
	fputs("Large Designs: by id\n", fp);
	fputs("</title>\n", fp);
	fputs("</head>\n", fp);
	
	fputs("<body>\n", fp);
	fputs("<h1>\n", fp);
	fputs("Large Designs: by id\n", fp);
	fputs("</h1>\n", fp);
	fputs("<hr>\n", fp);
	print_page_head(fp);
	
	btree->len(&len1);
	for (i = 0; i < len1; i++) {
		btree->ith(i, &key, &data);
		db->get_op(&data, &dp);
		id = dp.s_id_i();

		fprintf(fp, "<a name=\"%ld\"> <code>%3ld:</code> </a> ", id, id);
		write_head_HTML(fp, &dp);
		write_source_HTML(fp, &dp);


		if ((i + 1) % 10 == 0)
			fprintf(fp, "<p>\n");
		else
			fprintf(fp, "<br>\n");
		}
	fprintf(fp, "<p>\n");
	copyright(fp);
	fputs("</body>\n", fp);
	fputs("</html>\n\n", fp);
	fflush(fp);
	fclose(fp);
	printf("wrote file %s of size %ld\n", fname, file_size(fname));
	fflush(stdout);
	return OK;
}

void calc_prev_string(BYTE *s, INT prev)
{
		if (prev >= 0)
			sprintf(s, "<a href=\"%s#%ld\"> %ld </a>", DB_DP_HTML_FNAME4, prev, prev);
		else
			s[0] = 0;
}

static INT write_head_HTML(FILE *fp, DESIGN_PARAMETER_OP dp)
{
	BYTE s1[10000];
	DESIGN_PARAMETER_OB p_compl;
	INT lc;
	
	calc_prev_string(s1, dp->s_id_i());
	fprintf(fp, "%s %ld-(%ld,%ld,%ld) ", s1, 
		dp->s_t_i(), dp->s_v_i(), dp->s_k_i(), dp->s_lambda_i());
	dp->complement(&p_compl);
	lc = p_compl.s_lambda_i();
	if (lc != dp->s_lambda_i()) {
		fprintf(fp, "(complement has $\\lambda=%ld$) ", lc);
		}
	else {
		fprintf(fp, "(selfcomplementary) ");
		}
	return OK;
}

static INT write_source_HTML(FILE *fp, DESIGN_PARAMETER_OP dp)
{
	INT l, ii, prev;
	DESIGN_PARAMETER_SOURCE_OP S;
	BYTE s[1000];
	BYTE str0[10000];
	BYTE str1[10000];
	BYTE str2[10000];
	
	l = dp->s_source()->s_li();
	fprintf(fp, "<table>");
	for (ii = 0; ii < l; ii++) {
		S = dp->s_source_i(ii);
		S->sprint_text012(str0, str1, str2);
		prev = S->s_prev_i();
		calc_prev_string(s, prev);
		fprintf(fp, "<tr><td>%s %s %s</td></tr>\n", str1, s, str2);
		fprintf(fp, "\n");
		}
	fprintf(fp, "</table>");
	return OK;
}

void print_page_head(FILE *fp)
{
	BYTE str[BUFSIZE];
	FILE *fp1;
	
	fputs("<h2>\n", fp);
	fputs("Lehrstuhl II f&uuml;r Mathematik</h2>\n", fp);
	fputs("<strong>Universit&auml;t Bayreuth</strong>\n", fp);
	fputs("<br>\n", fp);
	fputs("<strong> D-95440 Bayreuth </strong>\n", fp);
	fputs("<p>\n", fp);
#if 0
	fputs("Dipl. Math. Anton Betten<br>\n", fp);
	fputs("Prof. Dr. Reinhard Laue<br>\n", fp);
	fputs("Dr. Alfred Wassermann<br>\n", fp);
#else
	fputs("<a href=\"http://www.mathe2.uni-bayreuth.de/betten/anton.html\">Anton Betten</a><br>\n", fp);
	fputs("<a href=\"http://www.mathe2.uni-bayreuth.de/people/laue.html\">Reinhard Laue</a><br>\n", fp);
	fputs("<a href=\"http://did.mat.uni-bayreuth.de/wassermann/wassermann.html\">Alfred Wassermann</a><br>\n", fp);
#endif
	fputs("<br>\n", fp);
	fputs("our <A href=\"designs.bib\">references in bibtex format</a><br>\n", fp);
	fputs("as these pages are a bit large you may want to download the compressed version and browse locally:<br>\n", fp);
	fputs("the file <A href=\"DATA/d1.html.gz\"> d1.html.gz (107 K) </a><br>\n", fp);
	fputs("the file <A href=\"DATA/d2.html.gz\"> d2.html.gz (107 K) </a><br>\n", fp);
	fputs("the file <A href=\"DATA/d3.html.gz\"> d3.html.gz (107 K) </a><br>\n", fp);
	fputs("the file <A href=\"DATA/d4.html.gz\"> d4.html.gz (304 K) </a><br>\n", fp);
	fputs("something about <A href=\"http://www.mathe2.uni-bayreuth.de/betten/DISCRETA/designs.html\">DISCRETA and designs</a><br>\n", fp);
	fputs("go to the <A href=\"http://www.mathe2.uni-bayreuth.de/betten/DISCRETA/Index.html\">DISCRETA home page</a><br>\n", fp);
	fputs("<p>\n", fp);
	system("rm a");
	system("date >a");
	fp1 = fopen("a", "r");
	fgets(str, BUFSIZE, fp1);
	fclose(fp1);
	fprintf(fp, "creation date: %s<p>\n", str);
	
	fprintf(fp, "WARNING: Completeness of know results "
		"is attempted only for t >= 6 and v <= 40 in the following list.<p>\n");
}

void copyright(FILE *fp)
{
#if 0
	fputs("<address>responsable for this page: \n", fp);
	fputs("<a href=\"mailto:Anton.Betten@uni-bayreuth.de\">\n", fp);
	fputs("Anton Betten\n", fp);
	fputs("</address></a><p>\n", fp);
	fputs("<a href=\"http://www.mathe2.uni-bayreuth.de/\">\n", fp);
	fputs("Back to the home page of <strong> Lehrstuhl II</strong></a>\n", fp);
#endif
}


#endif /* LADDER_TRUE */


