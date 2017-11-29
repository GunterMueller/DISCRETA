/* singer.C 
 * Anton Betten 
 * Jul 22, 1995 
 * 
 * this program replaces an earlier 
 * version written by Matthias Scholz;
 * his program was not free from bugs:
 * p=2, deg=8 gave a core dump. 
 */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */



#include <ctype.h>

#include <DISCRETA/discreta.h>

#ifdef SINGER_TRUE

#include <DISCRETA/gfq.h>

#if TEXDOCU
G_POLYNOM g_alloc(INT deg)
#endif
{
	G_POLYNOM p;

	p = (G_POLYNOM) my_malloc(sizeof(G_INT_TYPE) * (deg + 2), "singer.C: g_alloc");
	p[0] = deg;
	return p + 1;
}

#if TEXDOCU
void g_free(G_POLYNOM p)
#endif
{
	if (p == NIL)
		error("g_free(): NIL pointer");
	my_free(p - 1);
}

#if TEXDOCU
void g_zero(G_POLYNOM p)
#endif
{
	p[-1] = 0;
	p[0] = 0;
}

#if TEXDOCU
void g_one(G_POLYNOM p)
#endif
{
	p[-1] = 0;
	p[0] = 1;
}

#if TEXDOCU
void g_zero_coeffs(G_POLYNOM p)
#endif
{
	INT d, i;
	
	d = g_deg(p);
	for (i = 0; i <= d; i++)
		p[i] = 0;
}

#if TEXDOCU
void g_sprint(G_POLYNOM p, BYTE *str)
#endif
{
	G_INT_TYPE coef;
	INT i, deg, f_prev = FALSE;
	BYTE str1[256];
	
	strcat(str, "(");
	deg = g_deg(p);
	for (i = deg; i >= 0; i--) {
		coef = p[i];
		str1[0] = 0;
		if (coef == 0)
			continue;
		if (coef < 0) {
			strcat(str1, " -");
			coef = - coef;
			}
		else {
			if (f_prev)
				strcat(str1, " +");
			}
		if (!((coef == 1) && (i != 0)))
			sprintf(str1 + strlen(str1), 
				" %ld", (INT) coef);
		if (i == 1)
			strcat(str1, " X");
		else if (i > 1) {
			sprintf(str1 + strlen(str1), " X^%ld", i);
			}
		if (strlen(str) + strlen(str1) < 200)
			strcat(str, str1);
		else
			return;
		f_prev = TRUE;
		}
	strcat(str, ")");
}

#if TEXDOCU
void g_sprint_latex(G_POLYNOM p, BYTE *str)
#endif
{
	G_INT_TYPE coef;
	INT i, deg, f_prev = FALSE;
	BYTE str1[256];
	
	deg = g_deg(p);
	for (i = deg; i >= 0; i--) {
		coef = p[i];
		str1[0] = 0;
		if (coef == 0)
			continue;
		if (coef < 0) {
			strcat(str1, " -");
			coef = - coef;
			}
		else {
			if (f_prev)
				strcat(str1, " +");
			}
		if (!((coef == 1) && (i != 0)))
			sprintf(str1 + strlen(str1), 
				" %ld", (INT) coef);
		if (i == 1)
			strcat(str1, " X");
		else if (i > 1) {
			sprintf(str1 + strlen(str1), " X^{%ld}", i);
			}
		if (strlen(str) + strlen(str1) < 200)
			strcat(str, str1);
		else
			return;
		f_prev = TRUE;
		}
	strcat(str, "");
}

#if TEXDOCU
void g_print(G_POLYNOM p)
#endif
{
	BYTE str[256];

	str[0] = 0;
	g_sprint(p, str);
	printf("%s\n", str);
}

#if TEXDOCU
void g_print_latex(G_POLYNOM p)
#endif
{
	BYTE str[256];

	str[0] = 0;
	g_sprint_latex(p, str);
	printf("%s\n", str);
}

#if TEXDOCU
void g_derive(G_POLYNOM p, G_POLYNOM q)
#endif
{
	INT deg, i;
	G_INT_TYPE a;

	deg = g_deg(p);
#ifdef TEST_DEG
	if (g_deg(q) < deg - 1)
		error("g_derive(): g_deg(q) < deg - 1");
#endif
	for (i = deg; i > 0; i--) {
		a = p[i];
		a *= i;
		q[i - 1] = a;
		}
	g_set_deg(q, deg - 1);
}

#if TEXDOCU
void g_calc_deg(G_POLYNOM p, G_INT_TYPE chi)
#endif
{
	INT deg, i;
	G_INT_TYPE a;

	deg = g_deg(p);
	for (i = deg; i >= 0; i--) {
		a = g_asr(p[i], chi);
		p[i] = a;
		if (a != 0) {
			g_set_deg(p, i);
			return;
			}
		if (i == 0) {
			g_set_deg(p, 0);
			return;
			}			
		}
}

#if TEXDOCU
INT g_is_const(G_POLYNOM p, G_INT_TYPE chi)
#endif
{
	g_calc_deg(p, chi);
	return (g_deg(p) == 0);
}

#if TEXDOCU
INT g_is_zero(G_POLYNOM p, INT chi)
#endif
{
	g_calc_deg(p, chi);
	if (g_deg(p) > 0)
		return FALSE;
	if (g_asr(p[0], chi) == 0)
		return TRUE;
	return FALSE;
}

#if TEXDOCU
INT g_is_one(G_POLYNOM p, INT chi)
#endif
{
	g_calc_deg(p, chi);
	if (g_deg(p) > 0)
		return FALSE;
	if (g_asr(p[0], chi) == 1)
		return TRUE;
	return FALSE;
}

#if TEXDOCU
G_INT_TYPE g_as_const(G_POLYNOM p, G_INT_TYPE chi)
#endif
{
	INT deg;
	
	g_calc_deg(p, chi);
	deg = g_deg(p);
	if (deg > 1)
		return error("g_as_const(): not constant");
	return p[0];
}

#if TEXDOCU
INT g_cmp(G_POLYNOM p, G_POLYNOM q)
#endif
{
	INT deg_p, deg_q, l, a, b;

	deg_p = g_deg(p);
	deg_q = g_deg(q);
	if (deg_p < deg_q)
		return -1;
	if (deg_p > deg_q)
		return 1;
	for (l = deg_p; l >= 0; l--) {
		a = p[l];
		b = q[l];
		if (a < b)
			return -1;
		if (a > b)
			return 1;
		}
	return 0;
}

#if TEXDOCU
G_POLYNOM g_copy(G_POLYNOM p)
#else
returns a copy of $p$.
#endif
{
	G_POLYNOM q;
	INT d, i;
	
	d = g_deg(p);
	q = g_alloc(d);
	for (i = 0; i <= d; i++) {
		q[i] = p[i];
		}
	return q;
}

#if TEXDOCU
void g_move(G_POLYNOM p, G_POLYNOM q)
#else
$q := p$.
#endif
{
	INT d, i;

	d = g_deg(p);
#ifdef TEST_DEG
	if (g_deg(q) < d) {
		error("g_move(): g_deg(q) < d");
		return NIL;
		}
#endif
	for (i = d; i >= 0; i--) {
		q[i] = p[i];
		}
	g_set_deg(q, d);
}

#if TEXDOCU
INT g_add(G_POLYNOM a, G_POLYNOM b, G_POLYNOM c, G_INT_TYPE chi)
#else
$c := a + b$ in $GF(chi)$.
#endif
{
	INT a_deg, b_deg, c_deg;
	INT i;
	G_INT_TYPE c1;

	a_deg = g_deg(a);
	b_deg = g_deg(b);
	c_deg = MAXIMUM(a_deg, b_deg);
#ifdef TEST_DEG
	if (g_deg(c) < c_deg) {
		error("g_add(): deg(c) too small");
		return NIL;
		}
#endif
	g_set_deg(c, c_deg);
	g_zero_coeffs(c);
	for (i = 0; i <= c_deg; i++) {
		if (i <= a_deg && i <= b_deg) {
			c1 = g_asr(a[i] + b[i], chi);
			c[i] = c1;
			}
		else if (i <= a_deg)
			c[i] = a[i];
		else
			c[i] = b[i];
		}
	g_calc_deg(c, chi);
	return OK;
}

#if TEXDOCU
INT g_rank(G_INT_TYPE *M, INT m, INT n, G_INT_TYPE chi)
#else
returns the rank of the $m \times n$ matrix $M$ over $GF(chi)$.
#endif
{
	INT i, j, k, jj;
	G_INT_TYPE a, a1, a_inv, b, c;

	i = 0;
	for (j = 0; j < n; j++) {
		
		/* search pivot element: */
		for (k = i; k < m; k++) {
			if (g_asr(M[k * n + j], chi) != 0) {
				/* pivot element found: */
				if (k != i) {
					/* swap pivot row: */
					for (jj = j; jj < n; jj++) {
						a = M[i * n + jj];
						M[i * n + jj] = M[k * n + jj];
						M[k * n + jj] = a;
						}
					}
				break;
				}
			}

		if (k == m) {
			/* no pivot element found */
			/* leave i untouched */
			continue; /* next j */
			}
		/* make pivot to 1: */
		a = M[i * n + j];
		a_inv = g_inv_mod(a, chi);
		for (jj = j; jj < n; jj++) {
			M[i * n + jj] = 
				g_asr(a_inv * M[i * n + jj], chi);
			}

		/* do the gaussian elimination: */
		for (k = i + 1; k < m; k++) {
			a1 = M[k * n + j];
			a1 = g_asr(a1, chi);
			if (a1 == 0)
				continue;
			M[k * n + j] = 0;
			for (jj = j + 1; jj < n; jj++) {
				a = M[i * n + jj];
				b = M[k * n + jj];
				c = g_asr(b - a1 * a, chi);
				M[k * n + jj] = c;
				}
			}
		i++;
		}
	return i;
}

#if TEXDOCU
G_INT_TYPE *g_berlekamp_matrix(G_POLYNOM m, G_INT_TYPE chi)
#else
computes the berlekamp matrix. 
This matrix has $X^{chi^j}$ in its $j$-th column.
as the coefficients of the 
polynomial representation modulo $m$.
The matrix is stored in a $d \times d$ array. The pointer to 
this array is returned (with $d = deg(m)$).
#endif
{
	INT d, i, j, db;
	G_INT_TYPE *M, mii, mij;
	G_POLYNOM a, b, c;
	
	d = g_deg(m);
	if (d < 2) {
		error("g_berlekamp_matrix(): d < 2");
		return NIL;
		}

#ifdef DEBUG_BERLEKAMP
	printf("m = \n");
	g_print(m);
	fflush(stdout);
#endif
	a = g_alloc(d - 1);
	b = g_alloc(d - 1);
	c = g_alloc(d - 1);
	M = (G_INT_TYPE *) my_malloc(sizeof(G_INT_TYPE) * d * d, "singer.C: g_berlekamp_matrix");
	g_zero_coeffs(a);
	a[1] = 1; /* a := X */
	g_set_deg(a, 1);

	/* i = 0: (1 0 0 0 ...) */
	M[0] = 1;
	for (j = 1; j < d; j++)
		M[0 * d + j] = 0;
	
	/* i = 1: X^chi */
	g_power_mod_apply(a, chi, m, chi);
	for (j = 0; j < d; j++)
		M[1 * d + j] = a[j];

	g_move(a, b);

	for (i = 2; i < d; i++) {
		g_mult_mod(a, b, c, m, chi);
#ifdef DEBUG_BERLEKAMP
		/* printf("a,b,c=\n");
		g_print(a);
		g_print(b);
		g_print(c);
		fflush(stdout); */
#endif
		g_move(c, b);
		db = g_deg(b);
		for (j = 0; j < d; j++) {
			if (j <= db)
				mij = b[j];
			else
				mij = 0;
			M[i * d + j] = mij;
			}
		}
	
	/* M := M - 1 I , I identity matrix: */
	for (i = 0; i < d; i++) {
		mii = g_asr(M[i * d + i] - 1, chi);
		M[i * d + i] = mii;
		}
	
#ifdef DEBUG_BERLEKAMP
	for (i = 0; i < d; i++) {
		for (j = 0; j < d; j++) {
			printf("%3ld ", (INT) M[i * d + j]);
			}
		printf("\n");
		}
	printf("\n");
	fflush(stdout);
#endif
	g_free(a);
	g_free(b);
	g_free(c);
	return M;
}

#undef DEBUG_RUSS_POWER_MOD

#if TEXDOCU
INT g_power_mod_apply(G_POLYNOM p, INT l, G_POLYNOM m, G_INT_TYPE chi)
#else
Lueneburg~\cite{Lueneburg87a}:
$RussPower(a, b, l) := a * b^l; \mbox{ mod } m$. 
where $m$ is a polynomial over $GF(chi)$.
\begin{enumerate}
\item
$RussPower(a, b, - l) = RussPower(a, 1 / b, l); (PolyInverseMod(b, m))$
\item
$RussPower(a, b, l + 1) = a * b^{l + 1} = a * b * b^l = RussPower(a * b, b, l);$
\item
$RussPower(a, b, 2 * l) = a * b^{2 * l} = a * {(b^2)}^l = RussPower(a, b^2, l);$
\item
$RussPower(a, b, 0) = a * b^0 = a;$
\end{enumerate}
this algorithm is sometimes also called
\lq repeated squaring and multiplying\rq 
#endif
{
	INT d;
	G_POLYNOM a, b;

	d = g_deg(m) - 1;
	a = g_alloc(d);
	b = g_alloc(d);
	g_one(a);
	g_move(p, b);
	g_russ_power_mod(a, b, l, m, chi);
	g_move(a, p);
#ifdef DEBUG_RUSS_POWER_MOD
	printf("finished with g_russ_power_mod(): a=\n");
	g_print(a);
	printf("p=\n");
	g_print(p);
	fflush(stdout);
#endif
	return OK;
}

#if TEXDOCU
INT g_russ_power_mod(G_POLYNOM a, G_POLYNOM b, INT l, G_POLYNOM m, G_INT_TYPE chi)
#endif
{
	INT d;
	G_POLYNOM q;
	
	d = g_deg(m);
	g_reduce(a, m, chi);
	g_reduce(b, m, chi);
	q = g_alloc(d - 1);
	if (l <= 0) {
		return error("g_russ_power_mod(): l <= 0");
		}
	while (l != 0) {
#ifdef DEBUG_RUSS_POWER_MOD
		printf("l = %ld\n", l);
		printf("a=\n");
		g_print(a);
		printf("b=\n");
		g_print(b);
		fflush(stdout);
#endif
		if (ODD(l)) {
			/* Rule ii): 
			 * a := a * b; l--; */
			g_mult_mod(a, b, q, m, chi);
			g_move(q, a);
			l--;
			continue; /* l kann 0 geworden sein. */
			}
		/* now: EVEN(l) and l != 0 */
		/* Rule iii): 
		 * b := b * b; l := l / 2; */
		g_mult_mod(b, b, q, m, chi);
		g_move(q, b);
		l >>= 1;
		}
	g_free(q);
	return OK;
}

#if TEXDOCU
INT g_mult_mod(G_POLYNOM a, G_POLYNOM b, G_POLYNOM c, G_POLYNOM m, G_INT_TYPE chi)
#else
$c := a * b \mbox{ mod } m$. where $m$ is a polynomial over $GF(chi)$.
Lueneburg~\cite{Lueneburg87a}, p. 36.
#endif
{
	G_POLYNOM q;
	INT dm, da, db;
	INT i, j, l;
	G_INT_TYPE al, leit;
	
#ifdef DEBUG_MULT_MOD
	printf("in g_mult_mod(): a, b, m:\n");
	g_print(a);
	g_print(b);
	g_print(m);
#endif
	dm = g_deg(m);
	g_reduce(a, m, chi);
	g_reduce(b, m, chi);
#ifdef DEBUG_MULT_MOD
	printf("nach reduce: a,b=\n");
	g_print(a);
	g_print(b);
#endif
	q = g_alloc(dm);
	g_zero_coeffs(q);
	da = g_deg(a);
	db = g_deg(b);

	/* 
	 * q(x) := a[da] * b(x) 
	 */
	al = a[da];
	for (j = db; j >= 0; j--) {
		q[j] = g_asr(al * b[j], chi);
		}

	for (l = da - 1; l >= 0; l--) {
		/* 
		 * q(x) := X * q(x) mod m(x) 
		 */
		leit = q[dm - 1];
		if (leit == 0) {
			for (i = dm - 1; i > 0; i--)
				q[i] = q[i - 1];
			q[0] = 0;
			}
		else {
			for (i = dm - 1; i > 0; i--)
				q[i] = g_asr(q[i - 1] - 
					m[i] * leit, chi);
			q[0] = g_asr( - m[0] * leit, chi);
			}

		/*
		 * q(x) := q(x) + a[l] * b(x)
		 */
		al = a[l];
		for (i = 0; i <= db; i++) {
			q[i] = g_asr(q[i] + al * b[i], chi);
			}
		}
	g_set_deg(q, dm - 1);
	g_calc_deg(q, chi);
	g_move(q, c);
#ifdef DEBUG_MULT_MOD
	printf("c=\n");
	g_print(c);
#endif
	return OK;
}

#if TEXDOCU
INT g_mult(G_POLYNOM a, G_POLYNOM b, G_POLYNOM c, G_INT_TYPE chi)
#else
$c := a * b$ over $GF(chi)$.
#endif
{
	INT a_deg, b_deg, c_deg;
	INT i, j, ij;
	G_INT_TYPE c1;

	a_deg = g_deg(a);
	b_deg = g_deg(b);
	c_deg = a_deg + b_deg;
#ifdef TEST_DEG
	if (g_deg(c) < c_deg)
		return error("g_mult(): deg(c) too small");
#endif
	g_set_deg(c, c_deg);
	g_zero_coeffs(c);
	for (i = 0; i <= a_deg; i++) {
		for (j = 0, ij = i; j <= b_deg; j++, ij++) {
			c1 = g_asr(a[i] * b[j] + c[ij], chi);
			c[ij] = c1;
			}
		}
	return OK;
}

#if TEXDOCU
INT g_mult_scalar(G_POLYNOM a, G_POLYNOM b, G_INT_TYPE s, INT chi)
#else
$b := a * s$ with $s \in GF(chi)$.
#endif
{
	INT a_deg, i;

	a_deg = g_deg(a);
	b[-1] = a_deg;
	for (i = 0; i <= a_deg; i++) {
		b[i] = g_asr(a[i] * s, chi);
		}
	return OK;
}

#if TEXDOCU
INT g_mult_apply_scalar(G_POLYNOM a, G_INT_TYPE s, INT chi)
#else
$a := a * s$ with $s \in GF(chi)$.
#endif
{
	INT a_deg, i;

	a_deg = g_deg(a);
	for (i = 0; i <= a_deg; i++) {
		a[i] = g_asr(a[i] * s, chi);
		}
	return OK;
}

#if TEXDOCU
INT g_reduce(G_POLYNOM a, G_POLYNOM m, G_INT_TYPE chi)
#else
$a := a \mbox{ mod } m$ (all over $GF(chi)$).
#endif
{
	INT d;
	G_POLYNOM q, r;

	d = g_deg(a);
	if (d < g_deg(m))
		return OK;
	q = g_alloc(d);
	r = g_alloc(d);
	if (g_div_rem(a, m, q, r, chi) != OK)
		return ERROR;
	g_move(r, a);
	g_free(q);
	g_free(r);
	return OK;
}

#if TEXDOCU
INT g_div_rem(G_POLYNOM m, G_POLYNOM n, 
	G_POLYNOM q, G_POLYNOM r, G_INT_TYPE chi)
#else
q and r need to be already allocated for the correct size.
#endif
{
	INT m_deg, n_deg, q_deg, i, j, ii, jj;
	G_INT_TYPE a, a_inv, b, ba_inv;
	
#ifdef DEBUG_DIV_REM
	printf("div_rem: m, n:\n");
	g_print(m);
	g_print(n);
	fflush(stdout);
#endif
	g_calc_deg(m, chi);
	g_calc_deg(n, chi);
	m_deg = g_deg(m);
	n_deg = g_deg(n);
	if (n_deg == 0) {
		if (n[0] == 0)
			return error("g_div_rem(): division by zero");
		}
	if (n_deg > m_deg) {
		g_zero(q);
		g_move(m, r);
		return OK;
		}
	q_deg = m_deg - n_deg;
#ifdef TEST_DEG
	if (g_deg(q) < q_deg)
		return error("g_as_const(): g_deg(q) < q_deg");
#endif
	g_move(m, r);
	a = n[n_deg];
	a_inv = g_inv_mod(a, chi);
	for (i = m_deg, j = q_deg; i >= n_deg; i--, j--) {
		b = r[i];
		ba_inv = g_asr(b * a_inv, chi);
		q[j] = ba_inv;
		for (ii = i, jj = n_deg; jj >= 0; ii--, jj--) {
			r[ii] = g_asr(r[ii] - ba_inv * n[jj], chi);
			if (ii == i && r[ii] != 0)
				return error("g_div_rem(): ii == i && r[ii] != 0");
			}
		}
	g_set_deg(q, q_deg);
	g_set_deg(r, MAXIMUM(n_deg - 1, 0));
#ifdef DEBUG_DIV_REM
	printf("div_rem: q, r:\n");
	g_print(q);
	g_print(r);
	fflush(stdout);
#endif
	return OK;
}

#if TEXDOCU
INT g_is_squarefree(G_POLYNOM p, INT chi)
#else
checks if $deg(gcd(p, pd)) = 0$
where $pd$ is the formal derivative of $p$.
The trick is that repeated factors of $p$ 
must occur also in the formal derivative and so also in the gcd.
#endif
{
	G_POLYNOM pd, g;
	INT d, f_is_const;

	d = g_deg(p);
	pd = g_alloc(d);
	g = g_alloc(d);
	g_derive(p, pd);
	g_gcd(p, pd, g, chi);
	f_is_const = g_is_const(g, chi);
#ifdef DEBUG_IS_SQUAREFREE
	printf("p, pd, gcd =\n");
	g_print(p);
	g_print(pd);
	g_print(g);
	if (!f_is_const) {
		G_POLYNOM q, r;
		
		q = g_alloc(d);
		r = g_alloc(d);

		g_div_rem(p, g, q, r, chi);
		printf("p / g = , r = \n");
		g_print(q);
		g_print(r);
		g_free(q);
		g_free(r);
		}
	fflush(stdout);
#endif
	g_free(pd);
	g_free(g);
	return f_is_const;
}

#if TEXDOCU
INT g_is_irreducible(G_POLYNOM p, INT chi)
#else
checks if 
\begin{enumerate}
\item
$p$ is squarefree and
\item
the rank of berlekamp-matrix is maximal, e.g. $d - 1$
where $d$ is the degree of $p$.
\end{enumerate}
#endif
{
	G_INT_TYPE *M, d, r;
	
	if (!g_is_squarefree(p, chi)) {
#ifdef DEBUG_IS_IRREDUCIBLE
		printf("g_is_irreducible(): not squarefree\n");
		fflush(stdout);
#endif
		return FALSE;
		}

	d = g_deg(p);
	M = g_berlekamp_matrix(p, chi);
	r = g_rank(M, d, d, chi);
#ifdef DEBUG_IS_IRREDUCIBLE
	printf("g_is_irreducible(): deg = %ld rank = %ld\n", d, r);
	fflush(stdout);
#endif
	my_free(M);
	if (r == d - 1)
		return TRUE;
	return FALSE;
}

#if TEXDOCU
INT g_is_primitive(G_POLYNOM p, G_INT_TYPE chi, INT m, VECTOR_OP vp)
#else
returns TRUE iff the polynomial $X$ 
has order m modulo $p(x)$ over GF(chi); 
the prime factorization of $m$ is given 
(only the primes, not the exponents); 
we use the fact:
$X$ is primitive (has order $m$) iff 
$X^{m/p} != 1$ for all $p | m$.
#endif
{
	INT d, i, l, m1;
	G_POLYNOM x;

	d = g_deg(p);
	x = g_alloc(d);
	l = vp->s_li();
	for (i = 0; i < l; i++) {
		
		/* x := the polynomial X: */
		g_set_deg(x, 1);
		x[0] = 0;
		x[1] = 1;
		
		m1 = m / vp->s_ii(i);
		g_power_mod_apply(x, m1, p, chi);
		if (g_is_one(x, chi)) {
			g_free(x);
			return FALSE;
			}
		}
	g_free(x);
	return TRUE;
}

#if TEXDOCU
INT g_gcd(G_POLYNOM m, G_POLYNOM n, G_POLYNOM g, INT chi)
#else
$g := gcd(m,n)$ over $GF(chi)$.
#endif
{
	INT max_deg;
	G_POLYNOM m1, n1, q, r, tmp;
	
	if (g == NIL) {
		return error("g_gcd(): args NIL");
		}
	g_calc_deg(m, chi);
	g_calc_deg(n, chi);
	if (g_deg(m) < g_deg(n)) {
		return g_gcd(n, m, g, chi);
		}
	/* ab jetzt deg(n) <= deg(m) */
	if (g_is_zero(n, chi)) {
		/* printf("n is zero\n");
		fflush(stdout); */
		g_move(m, g);
		return OK;
		}
	max_deg = g_deg(m);
	m1 = g_copy(m);
	n1 = g_copy(n);
	q = g_alloc(max_deg);
	r = g_alloc(max_deg);
	tmp = g_alloc(max_deg);
	/* m1 = m; 
	 * n1 = n;  */
	while (TRUE) {
		if (g_div_rem(m1, n1, q, r, chi) != OK) {
			return error("g_gcd() div_rem");
			}
#ifdef DEBUG_GCD
		/* g_calc_deg(r, chi); */
		printf("deg(r) = %ld\n", (INT) g_deg(r));
		printf("m1, n1, q, r=\n");
		g_print(m1);
		g_print(n1);
		g_print(q);
		g_print(r);
		fflush(stdout);
#endif
		if (g_is_zero(r, chi)) {
			g_move(n1, g);
			goto l_free;
			}
		g_move(n1, m1);
		g_move(r, n1);
		}
l_free:
	g_free(m1);
	g_free(n1);
	g_free(q);
	g_free(r);
	g_free(tmp);
	return OK;
}

#if TEXDOCU
INT g_bezout(
	G_POLYNOM m, G_POLYNOM n, 
	G_POLYNOM u, G_POLYNOM v, 
	G_POLYNOM g, INT chi)
#else
$g := gcd(m, n) = u * m + v * n$ over $GF(chi)$.
#endif
{
	INT max_deg;
	G_POLYNOM m1, n1, q, r, tmp;
	G_POLYNOM u1, u2, u3, v1, v2, v3;
	
	if (u == NIL || v == NIL || g == NIL) {
		return error("g_bezout(): args NIL");
		}
	g_calc_deg(m, chi);
	g_calc_deg(n, chi);
	if (g_deg(m) < g_deg(n)) {
		return g_bezout(n, m, v, u, g, chi);
		}
	/* ab jetzt deg(n) <= deg(m) */
	if (g_is_zero(n, chi)) {
		/* printf("n is zero\n");
		fflush(stdout); */
		g_one(u);
		g_zero(v);
		g_move(m, g);
		return OK;
		}
	max_deg = g_deg(m);
	m1 = g_copy(m);
	n1 = g_copy(n);
	u1 = g_alloc(max_deg);
	u2 = g_alloc(max_deg);
	u3 = g_alloc(max_deg);
	v1 = g_alloc(max_deg);
	v2 = g_alloc(max_deg);
	v3 = g_alloc(max_deg);
	q = g_alloc(max_deg);
	r = g_alloc(max_deg);
	tmp = g_alloc(max_deg);
	g_one(u1);
	g_zero(u2);
	g_zero(v1);
	g_one(v2);
	/* m1 = m; u1 = 1L; v1 = 0L;
	n1 = n; u2 = 0L; v2 = 1L; */
	while (TRUE) {
		if (g_div_rem(m1, n1, q, r, chi) != OK) {
			return error("g_bezout() div_rem");
			}
		/* g_calc_deg(r, chi);
		printf("deg(r) = %ld\n", (INT) g_deg(r));
		fflush(stdout); */
		if (g_is_zero(r, chi)) {
			g_move(u2, u);
			g_move(v2, v);
			g_move(n1, g);
			goto l_free;
			}
		/* u3 := u1 - q * u2;
		 * u1 := u2; u2 := u3; */
		g_mult(q, u2, tmp, chi);
		g_mult_apply_scalar(tmp, (G_INT_TYPE) -1, chi);
		g_add(u1, tmp, u3, chi);
		g_move(u2, u1);
		g_move(u3, u2);
		
		/* v3 := v1 - q * v2;
		 * v1 := v2; v2 := v3; */
		g_mult(q, v2, tmp, chi);
		g_mult_apply_scalar(tmp, (G_INT_TYPE) -1, chi);
		g_add(v1, tmp, v3, chi);
		g_move(v2, v1);
		g_move(v3, v2);
		
		g_move(n1, m1);
		g_move(r, n1);
		
		/* u3 = u1 - q * u2;
		v3 = v1 - q * v2;
		m1 = n1; n1 = r;
		u1 = u2; u2 = u3;
		v1 = v2; v2 = v3; */
		}
l_free:
	g_free(m1);
	g_free(n1);
	g_free(u1);
	g_free(u2);
	g_free(u3);
	g_free(v1);
	g_free(v2);
	g_free(v3);
	g_free(q);
	g_free(r);
	g_free(tmp);
	return OK;
}

#if TEXDOCU
void g_numeric_pol(INT i, G_POLYNOM p, G_INT_TYPE chi)
#else
$p := \sum_j a_j X^j$ where 
$i = \sum_j a_j chi^j$ is the representation of $i$ to the base $chi$.
Useful to enumerate all elements of $GF(q)$ over $GF(p)$.
#endif
{
	INT j, r;
	
	if (i < 0) {
		error("g_numeric_pol(): i < 0");
		return;
		}
	if (i == 0) {
		g_zero(p);
		return;
		}
	j = 0;
	while (i != 0) {
		r = i % chi;
		p[j] = g_asr(r, chi);
		j++;
		i -= r;
		i /= chi;
		}
	j--;
	g_set_deg(p, j);
}

#if TEXDOCU
INT g_pol_numeric(G_POLYNOM p, G_INT_TYPE chi)
#else
returns $i := \sum_j a_j chi^j$ where $p = \sum_j a_j X^j$.
#endif
{
	INT i, n, n0, d, a;
	
	d = p[-1];
	n = 0;
	n0 = 1;
	for (i = 0; i <= d; i++) {
		a = p[i];
		if (a < 0)
			a += chi;
		n += a * n0;
		n0 *= chi;
		}
	return n;
}

#if TEXDOCU
G_POLYNOM g_singer(INT deg, G_INT_TYPE chi)
#else
computes an irreducible polynomial of degree deg over $GF(chi)$ 
which has the property that one of its roots (and therefore all its roots) 
are primitive elements of $GF(chi^{deg})$ over $GF(chi)$.
Such polynomials exist.
#endif
{
	G_POLYNOM p;
	INT a, b, m, i, low, high;
	VECTOR_OB vp, ve;

	if (chi <= 1) {
		error("g_singer(): chi <= 1");
		return NIL;
		}
	if (!is_prime(chi)) {
		error("g_singer(): chi not prime !");
		return NIL;
		}
	m = i_power_j(chi, deg) - 1;
	factor_integer(m, &vp, &ve);
#ifdef DEBUG_SINGER
	print_factorization(&vp, &ve, str);
	printf("m = %ld^%ld - 1 = %ld = %s\n", 
		chi, deg, m, str);
#endif
	a = primitive_root(chi);
#ifdef DEBUG_SINGER
	printf("g_singer(): a primitive root "
		"of %ld is %ld.\n", chi, a);
#endif

	p = g_alloc(deg);

	for (b = 0; b < chi; b++) {
		g_singer_kandidat(p, deg, b, a);
		if (g_singer_test(p, chi, m, &vp)) {
			return p;
			}
		}

	low = m + 1;
	high = low << 1; /* = 2 * low */
	for (i = low; i <= high; i++) {
		g_numeric_pol(i, p, chi);
		if (g_singer_test(p, chi, m, &vp)) {
			return p;
			}
		}

	g_free(p);
	printf("g_singer(): no primitive irreducible "
		"polynomial found for p = %ld deg = %ld\n", chi, deg);
	fflush(stdout);
	return NIL;
}

#if TEXDOCU
INT g_singer_test(G_POLYNOM p, G_INT_TYPE chi, INT m, VECTOR_OP vp)
#endif
{
	INT i;

#ifdef DEBUG_SINGER
	printf("singer testing:\n");
	g_print(p);
	fflush(stdout);
#endif	
	i = g_is_irreducible(p, chi);

#ifdef DEBUG_SINGER
	printf("is_irreducible = %ld\n", i);
	fflush(stdout);
#endif
	if (i == 0)
		return FALSE;

	i = g_is_primitive(p, chi, m, vp);

#ifdef DEBUG_SINGER
	printf("is_primitive = %ld\n", i);
	fflush(stdout);
#endif
	if (i == 0)
		return FALSE;
	return TRUE;
}

#if TEXDOCU
void g_singer_kandidat(G_POLYNOM p, INT k, G_INT_TYPE b, G_INT_TYPE a)
#else
generates the polynomial $X^k + X^{k-1} + b X + a$.
#endif
{
	INT i;

	if (g_deg(p) < k)
		error("g_singer_kandidat(): g_deg(p) < k");
	p[k] = (G_INT_TYPE) 1;
	if (k > 0)
		p[k - 1] = (G_INT_TYPE) 1;
	
	for (i = k - 2; i > 1; i--) {
		p[i] = (G_INT_TYPE) 0;
		}
	p[1] = b;
	p[0] = a;
}

#if TEXDOCU
G_INT_TYPE g_inv_mod(G_INT_TYPE a, G_INT_TYPE chi)
#else
computes $a^{-1} \mbox{ mod } chi$.
#endif
{
	INT a1 = (INT) a;
	INT p = (INT) chi;
	INT a_inv;

	if (a == 0)
		return error("g_inv_mod(): a == 0");
	inverse_mod_integer(a1, p, &a_inv);
	return (G_INT_TYPE) a_inv;
}

#if TEXDOCU
G_INT_TYPE g_asr(G_INT_TYPE a, G_INT_TYPE chi)
#else
absolutely smallest remainder.
#endif
{
	INT a1 = (INT) a;
	INT p = (INT) chi;
	INT a_rem;

	if (a1 == 0)
		return 0;
	asr(a1, p, &a_rem);
	return (G_INT_TYPE) a_rem;
}

#if TEXDOCU
int singer(long chi, long deg, long *pp)
#endif
{
	G_POLYNOM p;
	INT i;
	
	p = g_singer(deg, chi);
	if (p == NIL) {
		return error("no irreducible primitive polynomial found, help !");
		}
	/* printf("singer(): found polynomial\n");
	g_print(p);
	fflush(stdout); */
	for (i = 0; i <= deg; i++) {
		pp[i] = (long) p[i];
		if (pp[i] < 0)
			pp[i] += chi;
		}
	g_free(p);
	return OK;
}


#if 0
int main(int argc, char **argv)
{
	anfang();

	{
	G_POLYNOM p, q, M, n, r, g, pp;
	INT a, b, chi = 5, deg = 5, m, i;
	BYTE str[256];
	G_INT_TYPE *B;


	p = g_alloc(deg);
	q = g_alloc(deg);
	M = g_alloc(deg);
	n = g_alloc(deg);
	r = g_alloc(deg);
	g = g_alloc(deg);

#if 0
	for (b = 0; b < chi; b++) {
		g_singer_kandidat(p, deg, b, a);
		g_derive(p, q);
		g_print(p);
		g_print(q);
		printf("\n");
		fflush(stdout);
		}
#endif

#if 0
	g_zero_coeffs(M);
	g_zero_coeffs(n);
	M[deg] = 1;
	M[0] = -1;
	n[1] = 1;
	n[0] = -1;
	g_calc_deg(n, chi);
	g_div_rem(M, n, q, r, chi);
	printf("M, n, q, r:\n");
	g_print(M);
	g_print(n);
	g_print(q);
	g_print(r);
	printf("\n");
	fflush(stdout);
#endif

#if 0
	g_zero_coeffs(M);
	g_zero_coeffs(n);
	M[deg] = 1;
	M[0] = -1;
	n[3] = 1;
	n[0] = -1;
	g_calc_deg(n, chi);
	g_bezout(M, n, q, r, g, chi);
	printf("M, n, q, r, g:\n");
	g_print(M);
	g_print(n);
	g_print(q);
	g_print(r);
	g_print(g);
	printf("\n");
	fflush(stdout);
#endif
	
#if 0
	g_zero_coeffs(M);
	g_zero_coeffs(p);
	g_zero_coeffs(q);
	M[deg] = 1;
	g_one(p);
	q[1] = 1;
	q[0] = 1;
	g_russ_power_mod(p, q, 5, M, chi);
	printf("M, p, q:\n");
	g_print(M);
	g_print(p);
	g_print(q);
	printf("\n");
	fflush(stdout);
#endif

#if 0
	g_zero_coeffs(M);
	M[deg] = 1;
	M[2] = 3;
	M[1] = 2;
	M[0] = 2;
	B = g_berlekamp_matrix(M, chi);
#endif

#if 0
	g_zero_coeffs(M);
	M[deg] = 1;
	M[4] = 1;
	M[3] = -1;
	M[1] = -1;
	M[0] = 2;
	i = g_is_irreducible(M, chi);
	printf("is_irreducible = %ld\n", i);
	fflush(stdout);
	i = g_is_primitive(M, chi, m, &vp);
	printf("is_primitive = %ld\n", i);
	fflush(stdout);
#endif

	chi = 2;
	deg = 5;

	pp = g_singer(deg, chi);
	if (pp) {
		printf("pp:\n");
		g_print(pp);
		fflush(stdout);
		}

	g_free(p);
	g_free(q);
	g_free(M);
	g_free(n);
	g_free(r);
	g_free(g);
	}

	ende();
	return 0;
}
#endif

#if TEXDOCU
INT g_pol_r(G_POLYNOM a, G_POLYNOM b, G_POLYNOM r, G_INT_TYPE chi)
#else
computes the polynomial $r$ with
\begin{enumerate}
\item
$r$ divides $a$
\item
$gcd(r,b) = 1$ and
\item
each irreducible polynomial dividing $a/r$ divides $b$.
Lueneburg~\cite{Lueneburg87a}, p. 37.
\end{enumerate}
#endif
{
	G_POLYNOM rr, ggt, tmp;
	G_INT_TYPE linv;
	INT d;
	
	d = g_deg(a);
	rr = g_alloc(d);
	ggt = g_alloc(d);
	tmp = g_alloc(d);
	g_move(a, r);
	g_gcd(r, b, ggt, chi);
	while (g_deg(ggt)) {
		if (g_div_rem(r, ggt, rr, tmp, chi) != OK) {
			return error("g_pol_r() error in div_rem()");
			}
		g_move(rr, r);
		g_gcd(r, b, ggt, chi);
		}
	d = g_deg(r);
	linv = g_inv_mod(r[d], chi);
	g_mult_apply_scalar(r, linv, chi);
	g_free(rr);
	g_free(ggt);
	g_free(tmp);
	return OK;
}

#undef DEBUG_MTX_FROB

#if TEXDOCU
INT mtx_frobenius(G_INT_TYPE **M, G_POLYNOM m, G_INT_TYPE chi, INT f_v)
#else
cf. LOADFROB, Lueneburg~\cite{Lueneburg87a}, p. 101.
#endif
{
	G_INT_TYPE *M1;
	G_POLYNOM p, p1, p2;
	INT d = g_deg(m), i, j, d1;
	
	p = g_alloc(d);
	p1 = g_alloc(d);
	p2 = g_alloc(d);
#ifdef DEBUG_MTX_FROB
	printf("mtx_frobenius() deg = %ld p = %ld\n", d, (INT) chi);
	fflush(stdout);
#endif
	M1 = (G_INT_TYPE *) my_malloc(d * d * sizeof(G_INT_TYPE), "mtx_frobenius");
	/* first column: 1 0 0 \ldots 0 */
	M1[0] = (G_INT_TYPE) 1;
	for (i = 1; i < d; i++)
		M1[i * d + 0] = (G_INT_TYPE) 0;
	/* X^chi nach p: */
	p[-1] = 1;
	p[0] = 0;
	p[1] = 1;
	g_power_mod_apply(p, (INT) chi, m, chi);
#ifdef DEBUG_MTX_FROB
	printf("mtx_frobenius() X^%ld (mod m) =\n", (INT) chi);
	g_print(p);
	fflush(stdout);
#endif
	g_one(p1);
	for (j = 1; j < d; j++) {
		g_mult_mod(p, p1, p2, m, chi);
#ifdef DEBUG_MTX_FROB
		printf("mtx_frobenius() p * p1 =\n");
		g_print(p2);
		fflush(stdout);
#endif
		g_move(p2, p1);
		d1 = g_deg(p1);
		for (i = 0; i < d; i++) {
			if (i <= d1)
				M1[i * d + j] = p1[i];
			else
				M1[i * d + j] = (G_INT_TYPE) 0;
			}
		}
	if (f_v) {
		printf("mtx_frobenius:\n");
		for (i = 0; i < d; i++) {
			for (j = 0; j < d; j++) {
				printf("%ld ", (INT) M1[i * d + j]);
				}
			printf("\n");
			}
		}

	*M = M1;
	g_free(p);
	g_free(p1);
	g_free(p2);
	return OK;
}


#endif /* SINGER_TRUE */


