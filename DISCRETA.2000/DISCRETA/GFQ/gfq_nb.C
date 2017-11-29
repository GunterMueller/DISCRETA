/* gfq_nb.C 
 * 
 * normal bases for GF_q
 *
 * Anton Betten 
 * May 17, 1996 
 */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */

#include <DISCRETA/discreta.h>

#ifdef GFQ_TRUE

#include <DISCRETA/perm.h>
#include <DISCRETA/gfq.h>
#include <DISCRETA/part.h>
#include <DISCRETA/ma.h>

#undef DEBUG_MINPOL
#undef DEBUG_DEP
#undef DEBUG_ORD_IDEAL
#undef DEBUG_GENERATOR

static INT kung_mue_i(PARTITION_OP lambda_exp, INT i, INT m);
static INT degree_of(VECTOR_OP V, INT i);

#if TEXDOCU
INT gfq_GL_classes_bi(VECTOR_OP R, INT chi, INT k, INT f_v, INT f_vv)
#else
bi stands for \lq base image\rq i.~e. the 
matrices are represented by integers describing the 
image of the basis elements.
#endif
{
	VECTOR_OB bi;
	VECTOR_OB V1, V, V_mult, V_part;
	MATRIX_OB M;
	SYM_OB gl_order, tmp;
	SYM_OB co, cl, cl_sum;
	BYTE str[256];
	INT deg_min = 2;
	INT deg_max = k;
	INT l, i, a0, a;

	R->m_il(0);
	l = 0;
	((INTEGER_OP) &gl_order)->m_i(1);
	((INTEGER_OP) &cl_sum)->m_i(0);
	a0 = i_power_j(chi, k);
	for (i = 0; i < k; i++) {
		a = a0 - i_power_j(chi, i);
		((INTEGER_OP) &tmp)->m_i(a);
		tmp.mult_apply(&gl_order);
		}
	str[0] = 0;
	gl_order.sprint(str);
	printf("order GL_%ld(%ld) = %s\n", k, chi, str);
	fflush(stdout);

	/* the polynomial 1: */
	V1.m_il(1);
	V1.m_ii(0, 1);
	V.m_il(1);
	V1.swap((VECTOR_OP) V.s_i(0));

	// the polynomials X+a, a \in GF(q)^* (e.g. a != 0):
	for (i = 1; i < chi; i++) {
		V1.m_il(2);
		V1.m_ii(0, i);
		V1.m_ii(1, 1);
		V.inc();
		V1.swap((VECTOR_OP) V.s_i(i + 1 - 1));
		}
#if 0
	V.m_il(0);
#endif
	printf("calling gfq_irred_polynomials chi=%ld deg_min=%ld deg_max=%ld\n", chi, deg_min, deg_max);
	fflush(stdout);
	gfq_irred_polynomials(&V, chi, deg_min, deg_max, f_v, f_vv);
	printf("found %ld irreducible polynomials of degree %ld - %ld "
		"(plus the polynomial 1, plus X+a, a != 0)\n", 
		V.s_li(), deg_min, deg_max);
	fflush(stdout);

	if (gfq_gl_class_first(&V, &V_mult, &V_part, k)) {
		V_mult.println();
		V_part.println();
		fflush(stdout);
		gfq_gl_class2matrix(chi, &V, &V_mult, &V_part, k, &M);
		M.Print();
		fflush(stdout);
		M.latex(stdout);
		fflush(stdout);

		gfq_gl_class_centralizer_order_kung(chi, &V, &V_mult, &V_part, &co);
		str[0] = 0;
		co.sprint(str);
		printf("centralizer order %s\n", str);
		fflush(stdout);
		gl_order.ganzdiv(&co, &cl);
		str[0] = 0;
		cl.sprint(str);
		printf("class length %s\n", str);
		fflush(stdout);
		cl_sum.add(&cl, &tmp);
		tmp.swap(&cl_sum);
		str[0] = 0;
		cl_sum.sprint(str);
		printf("class length sum %s\n", str);
		fflush(stdout);

		gfq_mtx2bi(&M, &bi, chi);
		bi.println();
		R->inc();
		bi.swap((VECTOR_OP) R->s_i(l));
		l++;
		fflush(stdout);

		while (gfq_gl_class_next(&V, &V_mult, &V_part, k)) {
			V_mult.println();
			V_part.println();
			fflush(stdout);
			gfq_gl_class2matrix(chi, &V, &V_mult, &V_part, k, &M);
			M.Print();
			fflush(stdout);
			M.latex(stdout);
			fflush(stdout);

			gfq_gl_class_centralizer_order_kung(chi, &V, &V_mult, &V_part, &co);
			str[0] = 0;
			co.sprint(str);
			printf("centralizer order %s\n", str);
			fflush(stdout);
			gl_order.ganzdiv(&co, &cl);
			str[0] = 0;
			cl.sprint(str);
			printf("class length %s\n", str);
			fflush(stdout);
			cl_sum.add(&cl, &tmp);
			tmp.swap(&cl_sum);
			str[0] = 0;
			cl_sum.sprint(str);
			printf("class length sum %s\n", str);
			fflush(stdout);

			gfq_mtx2bi(&M, &bi, chi);
			bi.println();
			R->inc();
			bi.swap((VECTOR_OP) R->s_i(l));
			l++;
			fflush(stdout);

			}
		}
	return OK;
}

#if TEXDOCU
INT gfq_bi2mtx(MATRIX_OP M, VECTOR_OP bi, INT q, INT k)
#else
#endif
{
	INT i, j, n, a;
	
	M->m_ilih(k, k);
	for (j = 0; j < k; j++) {
		n = bi->s_ii(j);
		for (i = 0; i < k; i++) {
			a = n % q;
			n /= q;
			M->m_iji(i, j, a);
			}
		}
	return OK;
}

#if TEXDOCU
INT gfq_mtx2bi(MATRIX_OP M, VECTOR_OP bi, INT q)
#else
#endif
{
	INT i, j, k, n0, n, a;

	k = M->s_li();
	bi->m_il(k);
	for (j = 0; j < k; j++) {
		n0 = 1;
		n = 0;
		for (i = 0; i < k; i++) {
			a = M->s_iji(i, j);
			n += a * n0;
			n0 *= q;
			}
		bi->m_ii(j, n);
		}
	return OK;
}

#if TEXDOCU
INT gfq_gl_classes(VECTOR_OP V, INT q, INT k)
#else
#endif
{
	VECTOR_OB V_mult, V_part;
	MATRIX_OB M;

	if (gfq_gl_class_first(V, &V_mult, &V_part, k)) {
		V_mult.println();
		V_part.println();
		fflush(stdout);
		gfq_gl_class2matrix(q, V, &V_mult, &V_part, k, &M);
		M.Print();
		fflush(stdout);
		while (gfq_gl_class_next(V, &V_mult, &V_part, k)) {
			V_mult.println();
			V_part.println();
			fflush(stdout);
			gfq_gl_class2matrix(q, V, &V_mult, &V_part, k, &M);
			M.Print();
			fflush(stdout);
			}
		}
	return OK;
}

#if TEXDOCU
INT gfq_gl_class_centralizer_order_kung(INT q, 
	VECTOR_OP V, VECTOR_OP V_mult, VECTOR_OP V_part, SYM_OP co)
#else
Computes the centraizer order according to Kungs formula~\cite{Kung81}.
#endif
{
	/* VECTOR_OP pol; */
	PARTITION_OB lambda_exp, lambda_;
	PARTITION_OP lambda;
	SYM_OB int_ob, co1;
	INT a, l, m, d, i, j, b, mue_i, aa, bb, cc;
	
	l = V->s_li();
	((INTEGER_OP) co)->m_i(1);
	for (a = l - 1; a >= 0; a--) { /* for all polynomials: */
		m = V_mult->s_ii(a); /* this is called c by Fripertinger */
		if (m) {
			d = degree_of(V, a);
			if (m > 1) {
				lambda = (PARTITION_OP) V_part->s_i(a);
				lambda->t_VECTOR_EXPONENT(&lambda_exp);
				}
			else {
				lambda_.first(1);
				lambda_.t_VECTOR_EXPONENT(&lambda_exp);
				}
			
			/* here comes Kung's formula: */
			((INTEGER_OP) &co1)->m_i(1);
			for (i = 1; i <= m; i++) {
				b = lambda_exp.exp_nb_i_parts(i);
				if (b == 0)
					continue;
				for (j = 1; j <= b; j++) {
					mue_i = kung_mue_i(&lambda_exp, i, m);
					aa = i_power_j(q, d * mue_i);
					bb = i_power_j(q, d * (mue_i - j));
					cc = aa - bb;
					((INTEGER_OP) &int_ob)->m_i(cc);
					int_ob.mult_apply(&co1);
					}
				}
			co1.mult_apply(co);
			
			} /* if m */
		}
	return OK;
}

static INT kung_mue_i(PARTITION_OP lambda_exp, INT i, INT m)
{
	INT k, mue;
	
	mue = 0;
	for (k = 1; k <= i; k++) {
		mue += lambda_exp->exp_nb_i_parts(k) * k;
		}
	for (k = i + 1; k <= m; k++) {
		mue += lambda_exp->exp_nb_i_parts(k) * i;
		}
	return mue;
}

#if TEXDOCU
INT gfq_gl_class2matrix(INT q, VECTOR_OP V, VECTOR_OP V_mult, VECTOR_OP V_part, INT k, MATRIX_OP M)
#else
#endif
{
	VECTOR_OP pol;
	PARTITION_OP lambda;
	INT i, j, i0, l, a, aa, m, d, coef;
	INT sign;

	if (q == 2)
		sign = 1;
	else
		sign = -1;
	l = V->s_li();
	i0 = 0;
	M->m_ilih_n(k, k);
	for (a = l - 1; a >= 0; a--) {
		m = V_mult->s_ii(a);
		if (m) {
			pol = (VECTOR_OP) V->s_i(a);
			/* d = pol->s_li() - 1; */
			d = degree_of(V, a);
			for (aa = 0; aa < m; aa++) {
				/* fill in a companion matrix of type pol: */

				/* right hand side column: */
				for (i = 0; i < d; i++) {
					coef = sign * pol->s_ii(i);
					M->m_iji(i0 + i, i0 + d - 1, coef);
					}
				/* lower diagonal: */
				for (j = 0; j < d - 1; j++)
					M->m_iji(i0 + j + 1, i0 + j, 1);
				i0 += d;
				}
			}
		}
	if (i0 != k)
		return error("gfq_gl_class2matrix() i0 != k (first)");
	i0 = 0;
	for (a = l - 1; a >= 0; a--) {
		m = V_mult->s_ii(a);
		if (m) {
			pol = (VECTOR_OP) V->s_i(a);
			/* d = pol->s_li() - 1; */
			d = degree_of(V, a);
			if (m > 1) {
				INT ll, li, la;

				lambda = (PARTITION_OP) V_part->s_i(a);
				ll = lambda->s_li();
				for (li = ll - 1; li >= 0; li--) {
					la = lambda->s_ii(li);
					/* we have a block of la times the same 
					 * polynomial, join them by '1': */
					for (i = 0; i < la; i++) {
						if (i < la - 1) {
							M->m_iji(i0 + d, i0 + d - 1, 1);
							}
						i0 += d;
						}
					}
				}
			else { /* m == 1 */
				i0 += d;
				}
			}
		}
	if (i0 != k)
		return error("gfq_gl_class2matrix() i0 != k (second)");
	return OK;
}

#if TEXDOCU
INT gfq_gl_class_first(VECTOR_OP V, VECTOR_OP V_mult, VECTOR_OP V_part, INT k)
#else
#endif
{
	if (!gfq_choose_pol_first(V, V_mult, k))
		return FALSE;
	while (TRUE) {
		if (gfq_pol_part_first(V_mult, V_part))
			return TRUE;
		if (!gfq_choose_pol_next(V, V_mult))
			return FALSE;
		}
}

#if TEXDOCU
INT gfq_gl_class_next(VECTOR_OP V, VECTOR_OP V_mult, VECTOR_OP V_part, INT k)
#else
#endif
{
	if (gfq_pol_part_next(V_mult, V_part))
		return TRUE;
	while (TRUE) {
		if (!gfq_choose_pol_next(V, V_mult))
			return FALSE;
		if (gfq_pol_part_first(V_mult, V_part))
			return TRUE;
		}
}

#if TEXDOCU
INT gfq_pol_part_first(VECTOR_OP V_mult, VECTOR_OP V_part)
#else
#endif
{
	INT i, l, m;
	PARTITION_OP p;

	l = V_mult->s_li();
	V_part->m_il_n(l);
	for (i = l - 1; i >= 0; i--) {
		m = V_mult->s_ii(i);
		if (m > 1) {
			p = (PARTITION_OP) V_part->s_i(i);
			/* printf("calling p->first(%ld)\n", m);
			fflush(stdout); */
			p->first(m);
			/* printf("finished\n");
			fflush(stdout); */
			}
		}
	return TRUE;
}

#if TEXDOCU
INT gfq_pol_part_next(VECTOR_OP V_mult, VECTOR_OP V_part)
#else
#endif
{
	INT i, l, m, ret;
	PARTITION_OP p;

	l = V_mult->s_li();
	for (i = l - 1; i >= 0; i--) {
		m = V_mult->s_ii(i);
		if (m > 1) {
			p = (PARTITION_OP) V_part->s_i(i);
			ret = p->next_VECTOR(p);
			if (ret)
				return TRUE;
			p->first(m);
			}
		}
	return FALSE;
}

#if TEXDOCU
INT gfq_choose_pol(VECTOR_OP V, INT k)
#else
#endif
{
	VECTOR_OB V_mult;

	if (gfq_choose_pol_first(V, &V_mult, k)) {
		V_mult.println();
		while (gfq_choose_pol_next(V, &V_mult)) {
			V_mult.println();
			}
		}
	return OK;
}

#if TEXDOCU
INT gfq_choose_pol_first(VECTOR_OP V, VECTOR_OP V_mult, INT k)
#else
#endif
{
	INT i, l, k1 = k, d, m;

	l = V->s_li();
	V_mult->m_il_n(l);
	for (i = l - 1; i >= 0; i--) {
		d = degree_of(V, i);
		m = k1 / d;
		V_mult->m_ii(i, m);
		k1 -= m * d;
		if (k1 == 0)
			return TRUE;
		}
	if (k1 == 0)
		return TRUE;
	else
		return FALSE;
}

#if TEXDOCU
INT gfq_choose_pol_next(VECTOR_OP V, VECTOR_OP V_mult)
#else
#endif
{
	INT i, ii, l, k1 = 0, d, m;

	l = V->s_li();
	k1 = V_mult->s_ii(0) * degree_of(V, 0);
	V_mult->m_ii(0, 0);
	do {
		for (i = 1; i < l; i++) {
			m = V_mult->s_ii(i);
			if (m) {
				k1 += degree_of(V, i);
				m--;
				V_mult->m_ii(i, m);
				break;
				}
			}
		if (i == l)
			return FALSE;
		for (ii = i - 1; ii >= 0; ii--) {
			d = degree_of(V, ii);
			m = k1 / d;
			V_mult->m_ii(ii, m);
			k1 -= m * d;
			if (k1 == 0)
				return TRUE;
			}
		k1 += V_mult->s_ii(0) * degree_of(V, 0);
		V_mult->m_ii(0, 0);
	} while (k1);
	return FALSE;
}

static INT degree_of(VECTOR_OP V, INT i)
{
	INT d;

	d = ((VECTOR_OP) V->s_i(i))->s_li() - 1;
	if (d == 0)
		d++;
	/* printf("degree of %ld is %ld\n", i, d); */
	return d;
}

#if TEXDOCU
INT gfq_irred_polynomials(VECTOR_OP V, INT chi, INT deg_min, INT deg_max, INT f_v, INT f_vv)
#else
Computes all irreducible polynomials over $GF(chi)$ 
of degree deg with 
deg\_min $\le$ deg $\le$ deg\_max. 
The polynomials are stored in the vector $V$.
#endif
{
	INT i;

	/* V->m_il(0); */
	for (i = deg_min; i <= deg_max; i++) {
		gfq_irred_pol(V, chi, i, f_v, f_vv);
		}
	return OK;
}

#if TEXDOCU
INT gfq_irred_pol(VECTOR_OP V, INT chi, INT deg, INT f_v, INT f_vv)
#else
Computes all irreducible polynomials of degree deg over $GF(chi)$. 
The polynomials are stored in the vector $V$.
#endif
{
	VECTOR_OB V1;
	INT *nb, dim_nb;
	INT *v, *g;
	G_POLYNOM g_pol, mue;
	INT q, i, j, n;
	GALOIS GFq;
	G_POLYNOM m;
	G_INT_TYPE *M;
	INT f_vvv = f_vv;

	m = g_singer(deg, chi);
	if (f_vv) {
		printf("m:\n");
		g_print(m);
		fflush(stdout);
		}
	mtx_frobenius(&M, m, chi, f_vv);
	GFq.chi = chi;
	GFq.deg = deg;
	GFq.m = m;
	GFq.Frob = M;
	dim_nb = deg + 1;
	nb = (INT *) my_malloc(dim_nb * dim_nb * sizeof(INT), "test_gfq_nb() nb");
	for (i = 0; i < dim_nb * dim_nb; i++)
		nb[i] = 0;
	
	gfq_generator(&GFq, nb, dim_nb, f_v, f_vv);


	v = (INT *) my_malloc(deg * sizeof(INT), "test_gfq_nb() v");
	g = (INT *) my_malloc(deg * sizeof(INT), "test_gfq_nb() g");
	g_pol = g_alloc(deg + 1);
	mue = g_alloc(deg + 1);
	q = i_power_j(chi, deg);
	
	n = 0;
	if (gfq_first_irred_pol(q, chi, deg, v)) {

		do {
			n++;

			if (f_vvv) {
				printf("gfq_irred_pol: ");
				for (i = 0; i < deg; i++)
					printf("%ld ", v[i]);
				printf("\n");
				}
			
			for (i = 0; i < deg; i++) {
				g[i] = 0;
				for (j = 0; j < deg; j++) {
					g[i] = g_asr(g[i] + v[j] * nb[i * dim_nb + j], chi);
					}
				}
			/* now g is the element in the standard base */
			g_pol[-1] = deg;
			for (i = 0; i < deg; i++) {
				g_pol[i] = g[i];
				}
			g_pol[deg] = 0;
			
			gfq_minpol(&GFq, g_pol, mue);
			if (f_vv) {
				printf("gfq_irred_pol minpol= ");
				g_print(mue);
				fflush(stdout);
				}
			g_print_latex(mue);
			V1.m_il(deg + 1);
			for (i = 0; i <= deg; i++) {
				V1.m_ii(i, mue[i]);
				}
			V->inc();
			V1.swap((VECTOR_OP) V->s_i(V->s_li() - 1));
			

			} while (gfq_next_irred_pol(q, chi, deg, v));

		} /* ! first */
	

	if (f_v) {
		printf("found %ld irred. polynomials of degree %ld over GF_%ld\n", n, deg, chi);
		fflush(stdout);
		}
	
	my_free(M);

	g_free(g_pol);
	g_free(mue);
	my_free(v);
	my_free(g);
	my_free(nb);
	return OK;
}

#if TEXDOCU
void print_mtx(INT *M, INT dim_M, INT deg)
#else
#endif
{
	INT i, j, a;
	
	for (i = 0; i < deg; i++) {
		for (j = 0; j < deg; j++) {
			a = M[i * dim_M + j];
			printf("%ld ", a);
			}
		printf("\n");
		}
}

#if TEXDOCU
void test_gfq_nb(INT chi, INT deg, INT f_v)
#else
#endif
{
	INT *nb, dim_nb;
	INT *v, *g;
	G_POLYNOM g_pol, mue;
	INT q, i, j, n;
	GALOIS GFq;
	G_POLYNOM m;
	G_INT_TYPE *M;
	INT f_vv = FALSE;
	
	m = g_singer(deg, chi);
	if (m) {
		printf("m:\n");
		g_print(m);
		fflush(stdout);
		}
	mtx_frobenius(&M, m, chi, TRUE /* f_v */);
	GFq.chi = chi;
	GFq.deg = deg;
	GFq.m = m;
	GFq.Frob = M;
	dim_nb = deg + 1;
	nb = (INT *) my_malloc(dim_nb * dim_nb * sizeof(INT), "test_gfq_nb() nb");
	for (i = 0; i < dim_nb * dim_nb; i++)
		nb[i] = 0;
	
	gfq_generator(&GFq, nb, dim_nb, f_v, f_vv);


	v = (INT *) my_malloc(deg * sizeof(INT), "test_gfq_nb() v");
	g = (INT *) my_malloc(deg * sizeof(INT), "test_gfq_nb() g");
	g_pol = g_alloc(deg + 1);
	mue = g_alloc(deg + 1);
	q = i_power_j(chi, deg);
	
	n = 0;
	if (gfq_first_irred_pol(q, chi, deg, v)) {

		do {
			n++;

			if (f_vv) {
				for (i = 0; i < deg; i++)
					printf("%ld ", v[i]);
				printf("\n");
				}
			
			for (i = 0; i < deg; i++) {
				g[i] = 0;
				for (j = 0; j < deg; j++) {
					g[i] = g_asr(g[i] + v[j] * nb[i * dim_nb + j], chi);
					}
				}
			/* now g is the element in the standard base */
			g_pol[-1] = deg;
			for (i = 0; i < deg; i++) {
				g_pol[i] = g[i];
				}
			g_pol[deg] = 0;
			
			gfq_minpol(&GFq, g_pol, mue);
			if (f_v || 1) {
				g_print(mue);
				}

			} while (gfq_next_irred_pol(q, chi, deg, v));

		} /* ! first */
	
	
	printf("found %ld irred. polynomials of degree %ld over GF_%ld\n", n, deg, chi);
	fflush(stdout);
	
	my_free(M);

	g_free(g_pol);
	g_free(mue);
	my_free(v);
	my_free(g);
	my_free(nb);
}

#define MAX_DEG 10000

#if TEXDOCU
INT gfq_minpol(GALOIS *GFq, G_POLYNOM g, G_POLYNOM mue)
#else
Lueneburg~\cite{Lueneburg87a}, p. 112.
#endif
{
	INT chi, deg, i, j, k, ii, sign, notyet;
	G_POLYNOM h, h2, sigma[MAX_DEG];
	
	deg = GFq->deg;
	chi = GFq->chi;
	g_calc_deg(g, chi);
#ifdef DEBUG_MINPOL
	printf("minpol of ");
	g_print(g);
	fflush(stdout);
#endif
	for (i = 0; i <= deg + 1; i++) {
		sigma[i] = g_alloc(deg + 1);
		}
	h = g_alloc(deg + 1);
	h2 = g_alloc(deg + 1);
	sigma[0][-1] = 0;
	sigma[0][0] = 1;
	g_move(g, sigma[1]);
	i = 1;
	notyet = (sigma[1][-1] > 0);
	while (notyet) {
		i++;
		gfq_frob_apply(GFq, g);
#ifdef DEBUG_MINPOL
		printf("%ld: ", i);
		g_print(g);
		fflush(stdout);
#endif
		g_mult_mod(g, sigma[i - 1], sigma[i], GFq->m, chi);
		for (j = i - 1; j >= 1; j--) {
			g_mult_mod(g, sigma[j - 1], h, GFq->m, chi);
			g_add(sigma[j], h, h2, chi);
			g_move(h2, sigma[j]);
			}
		k = i;
		while (k > 0 && sigma[k][-1] <= 0)
			k--;
		notyet = (sigma[k][-1] > 0);
		}
	mue[-1] = i;
	sign = -1;
	for (ii = i; ii >= 0; ii--) {
		if (sigma[i - ii][-1] == 0 && sigma[i - ii][0] == 0) {
			mue[ii] = 0;
			}
		else {
			mue[ii] = g_asr(sign * sigma[i - ii][0], chi);
			}
		}
#ifdef DEBUG_MINPOL
		printf("minpol mue = ");
		g_print(mue);
		fflush(stdout);
#endif
	
	for (i = 0; i <= deg + 1; i++) {
		g_free(sigma[i]);
		}
	g_free(h);
	g_free(h2);
	return OK;
}

#if TEXDOCU
INT gfq_dep(GALOIS *GFq, INT *v, INT *a, INT dim_a, INT m, INT *rho)
#else
Lueneburg~\cite{Lueneburg87a} p. 104.
#endif
{
	INT i, j, k, deg, chi, f_null, pp;
	G_INT_TYPE a00;
	
	deg = GFq->deg;
	chi = GFq->chi;
#ifdef DEBUG_DEP
	printf("gfq_dep() m = %ld\n", m);
	print_mtx(a, dim_a, deg);
	fflush(stdout);
#endif
	/* a[m] = v^rho: */
	for (j = 0; j < deg; j++) {
		a[m * dim_a + j] = v[rho[j]];
		}
#ifdef DEBUG_DEP
	printf("a[%ld] = ", m);
	for (j = 0; j < deg; j++) {
		printf("%ld ", a[m * dim_a + j]);
		}
	printf("\n");
	fflush(stdout);
#endif
	for (k = 0; k < m; k++) {
	
#ifdef DEBUG_DEP
		printf("k = %ld\n", k);
		fflush(stdout);
#endif
		
		for (j = k + 1; j < deg; j++) {
		
#ifdef DEBUG_DEP
			printf("j = %ld\n", j);
			fflush(stdout);
#endif
			a[m * dim_a + j] = g_asr(a[k * dim_a + k] * a[m * dim_a + j], chi);
			a[m * dim_a + j] = g_asr(a[m * dim_a + j] - a[m * dim_a + k] * a[k * dim_a + j], chi);
			if (k > 0) {
				a00 = g_inv_mod(a[(k - 1) * dim_a + k - 1], chi);
				a[m * dim_a + j] = g_asr(a[m * dim_a + j] * a00, chi);
				}
			} /* next j */
		} /* next k */
#ifdef DEBUG_DEP
	printf("gfq_dep() now:\n");
	print_mtx(a, dim_a, deg);
	fflush(stdout);
#endif
	
	f_null = (m == deg);
	if (!f_null) {
	
		/* search for an non-zero entry in row m starting in column m.
		 * permute that column into column m, change the col-permutation rho */
		j = m;
		while ((a[m * dim_a + j] == 0) && (j < deg - 1))
			j++;
		f_null = (a[m * dim_a + j] == 0);
		if (!f_null && j > m) {
#ifdef DEBUG_DEP
			printf("choosing column %ld\n", j);
			fflush(stdout);
#endif
			for (i = 0; i <= m; i++) {
				pp = a[i * dim_a + m];
				a[i * dim_a + m] = a[i * dim_a + j];
				a[i * dim_a + j] = pp;
				} /* next i */
			pp = rho[m];
			rho[m] = rho[j];
			rho[j] = pp;
			}
		}
	return f_null;
}

#if TEXDOCU
INT gfq_ord_ideal(GALOIS *GFq, INT i, G_POLYNOM mue, INT *nb, INT dim_nb)
#else
Lueneburg~\cite{Lueneburg87a} p. 105.
#endif
{
	INT deg, chi;
	INT *v, *v1, *rho;
	INT ii, j, m, f_null;
	
	deg = GFq->deg;
	chi = GFq->chi;
#ifdef DEBUG_ORD_IDEAL
	printf("gfq_ord_ideal() i = %ld\n", i);
	fflush(stdout);
	/* for (i = 0; i < dim_nb * deg; i++)
		nb[i] = 0; */
#endif
	v = (INT *) my_malloc(deg * sizeof(INT), "gfq_ord_ideal() v");
	v1 = (INT *) my_malloc(deg * sizeof(INT), "gfq_ord_ideal() v1");
	rho = (INT *) my_malloc(deg * sizeof(INT), "gfq_ord_ideal() rho");
	for (ii = 0; ii < deg; ii++) {
		if (ii == i)
			v[ii] = 1;
		else
			v[ii] = 0;
		}
	for (ii = 0; ii < deg; ii++)
		rho[ii] = ii;
	
	m = 0;
	f_null = gfq_dep(GFq, v, nb, dim_nb, m, rho);
	while (!f_null) {
	
		/* apply frobenius (the images are written into the columns): */
		for (ii = 0; ii < deg; ii++) {
			v1[ii] = 0;
			for (j = 0; j < deg; j++) {
				v1[ii] = g_asr(v1[ii] + GFq->Frob[ii * deg + j] * v[j], chi);
				}
			}
		for (ii = 0; ii < deg; ii++)
			v[ii] = v1[ii];
		
		m++;
		f_null = gfq_dep(GFq, v, nb, dim_nb, m, rho);
		if (m == deg && !f_null)
			return error("gfq_ord_ideal() m == deg && ! f_null");
		}
	
	mue[-1] = m;
	mue[m] = 1;
	for (j = m - 1; j >= 0; j--) {
		mue[j] = nb[m * dim_nb + j];
		for (ii = m - 1; ii >= j + 1; ii--)
			mue[j] = g_asr(mue[j] + mue[ii] * nb[ii * dim_nb + j], chi);
		mue[j] = g_asr(- mue[j] * - g_inv_mod(nb[j * dim_nb + j], chi), chi);
		}
	
	my_free(v);
	my_free(v1);
	my_free(rho);
	return OK;
}

#if TEXDOCU
INT gfq_generator(GALOIS *GFq, INT *nb, INT dim_nb, INT f_v, INT f_vv)
#else
Lueneburg~\cite{Lueneburg87a} p. 106.
#endif
{
	INT i, j, deg, chi, d;
	G_POLYNOM mue, mue1, b, b1, v, v1, v2, ggt, q, r, r1, r2, tmp;
	
	deg = GFq->deg;
	chi = GFq->chi;
	d = deg + 1;
	mue = g_alloc(d);
	mue1 = g_alloc(d);
	b = g_alloc(d);
	b1 = g_alloc(d);
	v = g_alloc(d);
	v1 = g_alloc(d);
	v2 = g_alloc(d);
	ggt = g_alloc(d);
	q = g_alloc(d);
	r = g_alloc(d);
	r1 = g_alloc(d);
	r2 = g_alloc(d);
	tmp = g_alloc(d);
	
	i = 0;
	gfq_ord_ideal(GFq, i, mue, nb, dim_nb);
	if (f_vv) {
		printf("gfq_generator:\n");
		printf("%ld-th order ideal is generated by ", i);
		g_print(mue);
		fflush(stdout);
		}
	
	v[-1] = 0;
	v[0] = 1;
	while (mue[-1] < deg) {
		i++;
		gfq_ord_ideal(GFq, i, mue1, nb, dim_nb);
		if (f_vv) {
			printf("%ld-th order ideal is generated by mue1 = ", i);
			g_print(mue1);
			printf("mue = ");
			g_print(mue);
			printf("generator v = ");
			g_print(v);
			fflush(stdout);
			}
	
		g_gcd(mue, mue1, ggt, chi);
		if (f_vv) {
			printf("ggt = ");
			g_print(ggt);
			fflush(stdout);
			}
		
		if (ggt[-1] < mue1[-1]) {
		
			/* b = (0 0 \ldots 0 1) */
			b[-1] = i;
			for (j = 0; j < i; j++)
				b[j] = 0;
			b[i] = 1;
			
			g_div_rem(mue1, ggt, q, tmp, chi);
			if (f_vv) {
				printf("q = mue1 / ggt = ");
				g_print(q);
				fflush(stdout);
				}
			g_pol_r(mue, q, r, chi);
			if (f_vv) {
				printf("r = pol_r(mue, q) = ");
				g_print(r);
				fflush(stdout);
				}
			g_div_rem(mue, r, q, tmp, chi);
			if (f_vv) {
				printf("q = mue / r = ");
				g_print(q);
				fflush(stdout);
				}
			g_move(v, v1);
			
			/* Frob Modul structure: apply q to v (q is monic): */
			for (j = q[-1] - 1; j >= 0; j--) {
				gfq_frob_apply(GFq, v1);
				g_mult_scalar(v, v2, q[j], chi);
				g_add(v1, v2, tmp, chi);
				g_move(tmp, v1);
				} /* next j */
			/* now: ord(v1) = Ideal(r) 
			 * v1 = v *(mue/r)(Frob) = v * q (Frob)
			 */
			
			g_div_rem(mue, ggt, q, tmp, chi);
			g_pol_r(mue1, q, r1, chi);
			g_gcd(r, r1, ggt, chi);
			g_div_rem(r1, ggt, r2, tmp, chi);
			/* now: gcd(r, r2) = 1
			 * r * r2 = lcm(mue,mue1) 
			 */
			g_div_rem(mue1, r2, q, tmp, chi);
			g_move(b, b1);
			
			/* Frob Modul structure: apply q to b (q is monic): */
			for (j = q[-1] - 1; j >= 0; j--) {
				gfq_frob_apply(GFq, b1);
				g_mult_scalar(b, v2, q[j], chi);
				g_add(b1, v2, tmp, chi);
				g_move(tmp, b1);
				} /* next j */
			/* now: ord(b1) = Ideal(r2) 
			 * b1 = b *(mue1/r2)(Frob) = v * q (Frob)
			 */
			g_add(v1, b1, v, chi);
			/* ord(v) = Ideal(r * r2), 
			 * gcd(r, r2) = 1
			 */
			g_mult(r, r2, mue, chi);
			} /* if */
		} /* while */
	
	gfq_calc_nb(GFq, v, nb, dim_nb, f_v);
	if (f_vv) {
		printf("gfq_generator finished!\n");
		fflush(stdout);
		}
	
	g_free(mue);
	g_free(mue1);
	g_free(b);
	g_free(b1);
	g_free(v);
	g_free(v1);
	g_free(v2);
	g_free(ggt);
	g_free(q);
	g_free(r);
	g_free(r1);
	g_free(r2);
	g_free(tmp);
	return OK;
}

#if TEXDOCU
INT gfq_calc_nb(GALOIS *GFq, G_POLYNOM p, INT *nb, INT dim_nb, INT f_v)
#else
Computes a normal basis
#endif
{
	INT deg, chi, i, j, jj;
	INT *v1, *v2;
	
	deg = GFq->deg;
	chi = GFq->chi;
	v1 = (INT *) my_malloc(deg * sizeof(INT), "gfq_calc_nb() v1");
	v2 = (INT *) my_malloc(deg * sizeof(INT), "gfq_calc_nb() v2");
	for (i = 0; i < deg; i++) {
		if (i <= p[-1])
			v1[i] = p[i];
		else
			v1[i] = 0;
		}
	for (j = 0; j <= deg; j++) {
		if (f_v) {
			printf("nb[%ld]: ", j);
			for (i = 0; i < deg; i++) {
				printf("%ld ", v1[i]);
				}
			printf("\n");
			}
		if (j == deg)
			break;
		for (i = 0; i < deg; i++) {
			nb[i * dim_nb + j] = v1[i];
			}
		/* apply frobenius to v1: */
		for (i = 0; i < deg; i++) {
			v2[i] = 0;
			for (jj = 0; jj < deg; jj++) {
				v2[i] = g_asr(v2[i] + GFq->Frob[i * deg + jj] * v1[jj], chi);
				}
			}
		for (i = 0; i < deg; i++)
			v1[i] = v2[i];
		}
	my_free(v1);
	my_free(v2);
	return OK;
}

#if TEXDOCU
INT gfq_frob_apply(GALOIS *GFq, G_POLYNOM p)
#else
Applies the frobenius endomorphism.
compare FROBAUT, 
Lueneburg~\cite{Lueneburg87a} p. 102.
#endif
{
	INT deg, chi;
	INT *v1, *v2;
	INT i, j, d;

	deg = GFq->deg;
	chi = GFq->chi;
	v1 = (INT *) my_malloc(deg * sizeof(INT), "gfq_frob_apply() v1");
	v2 = (INT *) my_malloc(deg * sizeof(INT), "gfq_frob_apply() v2");
	d = p[-1];
	for (i = 0; i < deg; i++) {
		if (i <= d)
			v1[i] = p[i];
		else
			v1[i] = 0;
		}
	for (i = 0; i < deg; i++) {
		v2[i] = 0;
		for (j = 0; j < deg; j++) {
			v2[i] = g_asr(v2[i] + GFq->Frob[i * deg + j] * v1[j], chi);
			}
		}
	p[-1] = deg - 1;
	for (i = 0; i < deg; i++)
		p[i] = v2[i];
	g_calc_deg(p, chi);
	my_free(v1);
	my_free(v2);
	return OK;
}

#if TEXDOCU
INT gfq_is_regular_word(INT *v, INT n, INT q)
#else
Returns TRUE if the word v of length n is regular, i.~e. 
lies in an orbit of length $n$ under the action of the cyclic group 
$C_n$ acting on the coordinates. 
Lueneburg~\cite{Lueneburg87a} p. 118.
v is a vector over $\{0, 1, \ldots , q-1\}$
#endif
{
	INT i, k, ipk, f_rg;
	
	if (n == 1)
		return TRUE;
	k = 1;
	do {
		i = 0;
		ipk = i + k;
		while (v[ipk] == v[i] && i < n - 1) {
			i++;
			if (ipk == n - 1)
				ipk = 0;
			else
				ipk++;
			}
		f_rg = (v[ipk] < v[i]);
		k++;
	} while (f_rg && k <= n - 1);
	return f_rg;
}

#if TEXDOCU
INT gfq_first_irred_pol(INT q, INT p, INT deg, INT *v)
#else
Assumes $q = p^{deg}$. Computes the first irreducible polynomial 
of degree deg over $GF(p)$ described by a regular word.
#endif
{
	G_POLYNOM pol;
	INT i, ii;
	
	pol = g_alloc(deg + 1);
	for (i = 0; i < q; i++) {
		g_numeric_pol(i, pol, p);
		for (ii = 0; ii < deg; ii++) {
			if (ii <= pol[-1])
				v[ii] = pol[ii];
			else
				v[ii] = 0;
			}
		if (gfq_is_regular_word(v, deg, p)) {
			g_free(pol);
			return TRUE;
			}
		}
	g_free(pol);
	return FALSE;
}

#if TEXDOCU
INT gfq_next_irred_pol(INT q, INT p, INT deg, INT *v)
#else
Assumes $q = p^{deg}$. 
This function computes the next irreducible polynomial of degree deg 
over $GF(p)$
described by a regular word of length deg over $GF(p)$.
#endif
{
	G_POLYNOM pol;
	INT i, ii;
	
	pol = g_alloc(deg + 1);
	for (ii = 0; ii < deg; ii++) 
		pol[ii] = v[ii];
	pol[-1] = deg - 1;
	i = g_pol_numeric(pol, p);
	for (i++; i < q; i++) {
		g_numeric_pol(i, pol, p);
		for (ii = 0; ii < deg; ii++) {
			if (ii <= pol[-1])
				v[ii] = pol[ii];
			else
				v[ii] = 0;
			}
		if (gfq_is_regular_word(v, deg, p)) {
			g_free(pol);
			return TRUE;
			}
		}
	g_free(pol);
	return FALSE;
}


#endif /* GFQ_TRUE */


