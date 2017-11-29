/* in1.C: */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>
#include <DISCRETA/part.h>
#include <DISCRETA/ma.h> // for Binomial()

static INT Binomial_using_table(INT n, INT k, MATRIX_OP T, SYM_OP res);

#if TEXDOCU
INT ord_to_order_if_prime(INT ord1, INT *ord)
#endif
{
	if (ord1 == 1) {
		*ord = 0;
		return(OK);
		}
	if (smallest_primedivisor(ord1) == ord1)
		*ord = ord1;
	else
		*ord = 0;
	return OK;
}

#if TEXDOCU
INT ord_to_order_if_prime_power(INT ord1, INT *ord, INT *prime, INT *k)
#endif
{
	INT prime1, prime2, k1;
	
	if (ord1 == 1) {
		*ord = 0;
		return(OK);
		}
	*ord = ord1;
	prime1 = 0;
	k1 = 0;
	while (ord1 != 1) {
		prime2 = smallest_primedivisor(ord1);
		if (prime1) {
			if (prime2 != prime1) {
				/* different primes */
				*ord = 0;
				return(OK);
				}
			}
		prime1 = prime2;
		k1++;
		ord1 /= prime1;
		}
	*prime = prime1;
	*k = k1;
	return OK;
}

#if TEXDOCU
INT kgv_iipi(INT m, INT n, INT *t)
#else
Computes the lcm (kgV) of $m$ and $n$ into $t$.
#endif
{
	INT g, m1;

	ggt_iipi(m, n, &g);
	m1 = m / g;
	*t = m1 * n;
	return OK;
}

#if TEXDOCU
INT ggt_iipi(INT m, INT n, INT *t)
#else
Computes the gcd (ggT) of $m$ and $n$ into $t$.
#endif
{
	INT r, s;
	
	if (n > m) {
		r = m;
		m = n;
		n = r;
		}
	if (n == 0) {
		*t = m;
		return OK;
		}
	while (TRUE) {
		s = m / n;
		r = m - (s * n);
		if (r == 0) {
			*t = n;
			return(OK);
			}
		m = n;
		n = r;
		}
}

#if TEXDOCU
INT smallest_primedivisor(INT n)
#else
Computes the smallest prime dividing $n$. 
The algorithm is based on Lueneburg.
#endif
/* Heinz Lueneburg: 
 * On The Rational Normal Form 
 * of Endomorphisms. */
{
	INT flag, i, q;
	
	if (EVEN(n))
		return(2);
	if ((n % 3) == 0)
		return(3);
	i = 5;
	flag = 0;
	while (TRUE) {
		q = n / i;
		if (n == q * i)
			return(i);
		if (q < i)
			return(n);
		if (flag)
			i += 4;
		else
			i += 2;
		flag = !flag;
		}
}

#if TEXDOCU
INT sp_ge(INT n, INT p_min)
#else
Computes the smalles prime dividing $n$ 
which is greater than or equal to p\_min.
#endif
/* smallest prime greater or equal to */
/* AB 230594 */
{
	INT flag, i, q;
	
	if (p_min == 0)
		p_min = 2;
	if (p_min < 0)
		p_min = - p_min;
	if (p_min <= 2) {
		if (EVEN(n))
			return 2;
		p_min = 3;
		}
	if (p_min <= 3) {
		if ((n % 3) == 0)
			return 3;
		p_min = 5;
		}
	if (EVEN(p_min))
		p_min--;
	i = p_min;
	if (EVEN((p_min - 1) >> 1))
		/* p_min cong 1 mod 4 ? */
		flag = 0;
	else
		flag = 1;
	while (TRUE) {
		q = n / i;
		if (n == q * i)
			return(i);
		if (q < i)
			return(n);
		if (flag)
			i += 4;
		else
			i += 2;
		flag = !flag;
		}
}

#if TEXDOCU
INT ny_p(INT n, INT p)
#else
Returns the integer $k$ with $n=p^k n'$ and $p \nmid n'$.
#endif
{
	INT ny_p;
	
	if (n == 0) {
		error("ny_p()|n == 0");
		return 0;
		}
	if (n < 0)
		n = -n;
	ny_p = 0;
	while (n != 1) {
		if ((n % p) != 0)
			break;
		n /= p;
		ny_p++;
		}
	return ny_p;
}

#if TEXDOCU
INT nb_primes(INT n)
#else
Returns the number of primes in the prime factorization 
of $n$ (including multiplicities).
#endif
{
	INT i = 0;
	INT d;
	
	if (n < 0)
		n = -n;
	while (n != 1) {
		d = smallest_primedivisor(n);
		i++;
		n /= d;
		}
	return i;
}

#if TEXDOCU
INT is_prime(INT n)
#else
TRUE if and only if $n$ is prime.
#endif
{
	if (smallest_primedivisor(n) == n)
		return TRUE;
	else
		return FALSE;
}

#if TEXDOCU
INT is_prime_power(INT n, INT p)
#else
TRUE if and only if $n$ is a power of $p$.
#endif
{
	if (n == 0) {
		error("is_prime_power()|n == 0");
		return 0;
		}
	if (n < 0)
		n = -n;
	while (n != 1) {
		if ((n % p) != 0)
			return FALSE;
		n /= p;
		}
	return TRUE;
}

#if TEXDOCU
INT factor_prime_power(INT n, INT *p, INT *e)
#else
Computes $p$ and $e$ with $n=p^e$. 
If $n$ is not a prime power, an error is raised.
#endif
{
	VECTOR_OB vp, ve;

	factor_integer(n, &vp, &ve);
	if (vp.s_li() != 1)
		return error("factor_prime_power() the number is not a prime power");
	*p = vp.s_ii(0);
	*e = ve.s_ii(0);
	return OK;
}

#if TEXDOCU
INT i_power_j(INT i, INT j)
#else
Computes $i^j$ as integer.
There is no checking for overflow.
#endif
{
	INT k, r = 1;

	for (k = 0; k < j; k++)
		r *= i;
	return r;
}

#if TEXDOCU
INT eulerfunc(INT n)
#else
Computes Eulers $\varphi$-function for $n$.
Uses the prime factorization of $n$.
#endif
{
	VECTOR_OB p, e;
	INT i, k, p1, e1;
	
	if (n <= 0)
		return error("eulerfunc(): n <= 0");
	if (n == 1)
		return 1;
	factor_integer(n, &p, &e);
	k = 1;
	for (i = 0; i < p.s_li(); i++) {
		p1 = p.s_ii(i);
		e1 = e.s_ii(i);
		if (e1 > 1)
			k *= i_power_j(p1, e1 - 1);
		k *= (p1 - 1);
		}
	return k;
}

#if TEXDOCU
INT moebius(INT i)
#else
Computes the number-theoretic $\mu$ (= moebius) function of $i$.
#endif
{
	VECTOR_OB vp, ve;
	INT j, a, l;
	
	factor_integer(i, &vp, &ve);
	l = vp.s_li();
	for (j = 0; j < l; j++) {
		a = ve.s_ii(j);
		if (a > 1)
			return 0;
		}
	if (EVEN(l))
		return 1;
	else
		return -1;
}

#if TEXDOCU
INT factor_integer(INT n, VECTOR_OP primes, VECTOR_OP exponents)
#else
Factors the integer $n = \prod_{i=1}^r p_i^{e_i}$. 
The vector primes holds the primes $p_i$,
the vector exponents holds the $e_i$.
#endif
{
	INT d, last_prime = 2;
	INT erg = OK;
	
	if (n == 0)
		return error("factor_integer(): "
		"tried to factor 0");
	if (n == 1) {
		primes->m_il(0);
		exponents->m_il(0);
		return erg;
		}
	d = sp_ge(n, last_prime);
	primes->m_il(1);
	exponents->m_il(1);
	primes->m_ii(0, d);
	exponents->m_ii(0, 1);
	last_prime = d;
	n /= d;
	while (n != 1) {
		d = sp_ge(n, last_prime);
		if (d == last_prime) {
			exponents->s_i(
			primes->s_li() - 1)->inc();
			}
		else {
			primes->inc();
			exponents->inc();
			primes->m_ii(primes->s_li() - 1, d);
			exponents->m_ii(primes->s_li() - 1, 1);
			last_prime = d;
			}
		n /= d;
		}
	return erg;
}

#if TEXDOCU
INT print_factorization(VECTOR_OP primes, VECTOR_OP exponents, BYTE *str)
#else
Prints the factorization into a string.
#endif
{
	INT i, p, e;
	
	for (i = 0; i < primes->s_li(); i++) {
		if (i >= exponents->s_li())
			break;
		p = primes->s_ii(i);
		e = exponents->s_ii(i);
		if (e > 1)
			sprintf(str + strlen(str), "%ld^%ld", p, e);
		else
			sprintf(str + strlen(str), "%ld", p);
		if (i < primes->s_li() - 1)
			strcat(str, " ");
		}
	return OK;
}

#if TEXDOCU
INT sieve_primes(VECTOR_OP Primes, INTEGER_OP nb_primes, INT *first, INT *len)
#else
Computes the prime numbers in the interval $[first, first + len - 1]$ 
using a sive method. 
Primes must hold all smaller primes. 
nb\_primes gives the numerb of used elements in Primes.
#endif
/* Siebe len Zahlen ab first aufsteigend 
 * modulo der in Primes 
 * bereits gespeicherten Primzahlen. 
 * Neu gefundene Primzahlen
 * werden in Primes eingetragen 
 * und nb_primes wird raufgezaehlt. 
 * first wird auf first + len gesetzt */
{
	char *sieve = NIL;
	INTEGER_OB int_ob;
	INT erg = OK, i, j, k, r, prime;

	sieve = (char * ) my_malloc(sizeof(char) * *len, "sieve_primes");
	if (sieve == NIL)
		return error("sieve_primes()|no memory for sieve");
	for (i = 0; i < *len; i++)
		sieve[i] = 0;
	if (*first == 0) {
		Primes->m_il(VECTOR_OVERSIZE);
		nb_primes->m_i(0);
		sieve[0] = 1;
		sieve[1] = 1;
		}
	for (k = 0; k < nb_primes->s_i(); k++) {
		i = Primes->s_ii(k);
		if (*first == 0) /* nie erfuellt ! */
			j = i;
		else {
			r = *first % i;
			if (r)
				j = i - (*first % i);
			else
				j = 0;
			}
		for (; j < *len; j += i)
			sieve[j] = 1;
		}
	for (i = 0; i < *len; i++) {
		if (sieve[i] == 0) {
			prime = i + *first;
			int_ob.m_i(prime);
			Primes->append_element(
			nb_primes, &int_ob);
			for (j = (prime << 1) - *first; 
				j < *len; j += prime)
				sieve[j] = 1;
			}
		}
	*first += *len;
	if (sieve) {
		my_free(sieve);
		sieve = NIL;
		}
	return erg;
}

#if TEXDOCU
INT order_mod_p(INT a, INT p)
#else
Computes the order of $a$ mod $p$, i.~e. the smallest $k$ 
s.~th. $a^k \equiv 1$ mod $p$.
#endif
{
	INT o, b;
	
	if (a < 0)
		return error("order_mod_p() a < 0");
	a %= p;
	if (a == 0)
		return 0;
	if (a == 1)
		return 1;
	o = 1;
	b = a;
	while (b != 1) {
		b *= a;
		b %= p;
		o++;
		}
	return o;
}

#if TEXDOCU
INT primitive_root(INT p)
#else
Computes a primitive element for $\EZ_p$, i.~e. an integer $k$ 
with $2 \le k \le p - 1$ s.~th. the order of $k$ mod $p$ is $p-1$.
#endif
{
	INT i, o;

	if (p < 2)
		return error("primitive_root(): p < 2");
	if (p == 2)
		return 1;
	for (i = 2; i < p; i++) {
		o = order_mod_p(i, p);
		if (o == p - 1) {
			printf("%ld is primitive root mod %ld\n", i, p);
			return i;
			}
		}
	return error("no primitive root found");
#if 0
	VECTOR_OB vp, ve, V;
	INT nb_primes, i, q;
	INT j, jhq, k;
	
	factor_integer(p - 1, &vp, &ve);
	nb_primes = vp.s_li();
	V.m_il_n(p);
	for (i = 0; i < nb_primes; i++) {
		q = vp.s_ii(i);
		for (j = 1; j < q; j++) {
			jhq = 1;
			for (k = 0; k < q; k++) 
				jhq = (jhq * j) % p;
			V.m_ii(jhq, 1);
			}
		}
	for (i = 1; i < p; i++) {
		if (V.s_ii(i) == 0) {
			printf("a primitive root mod %ld is %ld\n", p, i);
			return i;
			}
		}
	printf("p = %ld\n", p); fflush(stdout);
	return error("primitive_root() no primitive root found");
#endif
#if 0
	INT gfp_kind = ik_gfp(p), i, o;
	INTEGER_OB i1;
	
	if (p < 2)
		return error("primitive_root(): p < 2");
	if (p == 2)
		return 1;
	i1.init(gfp_kind);
	for (i = 2; i < p; i++) {
		i1.homo_z(i);
		i1.order(&o);
		if (o == p - 1) {
			printf("%ld is primitive root mod %ld\n", i, p);
			return i;
			}
		}
	return error("no primitive root found");
#endif
	
}

#if TEXDOCU
INT bezout_integer(INT m, INT n, INT *u, INT *v, INT *g)
#else
Computes the gcd of $m$ and $n$ into $g$. 
In addition, the numbers $u$ and $v$ are computed such that 
$g = um + vn$ holds.
#endif
/* Findet den positiven ggT(m, n) mit Vorfaktoren: 
 * g = ggT(m, n) = u * m + v * n */
{
	INT m1, n1, q, r;
	INT u1, u2, u3, v1, v2, v3;
	
	if (u == NIL || v == NIL || g == NIL) {
		Srfs("bezout_integer", "args NIL");
		return(ERROR);
		}
	if (m < 0L) {
		m1 = - m;
		if (bezout_integer(m1, n, u, v, g) != OK) {
			Srff("bezout_integer", "bezout_integer");
			return(ERROR);
			}
		*u = - *u;
		return(OK);
		}
	if (n < 0L) {
		n1 = - n;
		if (bezout_integer(m, n1, u, v, g) != OK) {
			Srff("bezout_integer", "bezout_integer");
			return(ERROR);
			}
		*v = - *v;
		return(OK);
		}
	/* ab jetzt m und n positiv:*/
	if (m < n) {
		return(bezout_integer(n, m, v, u, g));
		}
	/* ab jetzt m, n positiv und n < m */
	if (n == 0L) {
		*u = 1L;
		*v = 0L;
		*g = m;
		return(OK);
		}
	m1 = m; u1 = 1L; v1 = 0L;
	n1 = n; u2 = 0L; v2 = 1L;
	while (TRUE) {
		if (div_rem(m1, n1, &q, &r) != OK) {
			Srff("bezout_integer", "div_rem");
			return(ERROR);
			}
		if (r == 0L) {
			if (n1 < 0L) {
				u2 = - u2;
				v2 = - v2;
				n1 = - n1;
				}
			*u = u2;
			*v = v2;
			*g = n1;
			return(OK);
			}
		u3 = u1 - q * u2;
		v3 = v1 - q * v2;
		m1 = n1; n1 = r;
		u1 = u2; u2 = u3;
		v1 = v2; v2 = v3;
		}
}

#if TEXDOCU
INT asr(INT a, INT m, INT *r)
#else
absolute smallest remainder: Computes $r$ such that 
$a \equiv r$ mod $m$ and $- \frac{m}{2} < r \le \frac{m}{2}$ holds.
#endif
{
	return abs_sm_rem(a, m, r);
}

#if TEXDOCU
INT Asr(INT a, INT m)
#else
returns the absolute smallest remainder of $a$ mod $m$.
#endif
{
	INT r;
	
	abs_sm_rem(a, m, &r);
	return r;
}

#if TEXDOCU
INT abs_sm_rem(INT a, INT m, INT *r)
#else
absolute smallest remainder: Computes $r$ such that 
$a \equiv r$ mod $m$ and $- \frac{m}{2} < r \le \frac{m}{2}$ holds.
#endif
/* Berechnet r als Rest von a mod m: 
 * - m / 2 < r <= m / 2.
 * r kann mit a uebereinstimmen (der Pointer). 
 * Bei m == 0L wird mit Warnung zurueckgekehrt 
 * (Rueckgabe OK) dann: r = a. */
{
	INT q, m0, m1, m_halbe;
	
	if (m == 0) {
		Srfs("abs_sm_rem", "warning: m == 0L");
		*r = a;
		return(OK);
		}
	while (a < 0)
		a += m;
	m0 = m;
	m_halbe = m0 >> 1;
	q = a / m0;
	m1 = a - q * m0;
	if (m1 > m_halbe)
		m1 -= m0;
	if (ODD(m0)) {
		if (m1 < - m_halbe)
			m1 += m0;
		}
	else {
		if (m1 <= - m_halbe)
			m1 += m0;
		}
	*r = m1;
	return(OK);
}

#if TEXDOCU
INT div_rem(INT m, INT n, INT *q, INT *r)
#else
Computes $q$ and $r$ such that 
$m = q n + r$ and $0 \le r < n$ holds.
#endif
/* Division mit Rest von m duch n:
 * m = q * n + r mit 0 <= r < n 
 * Fehler, wenn n == 0 */
{
	INT d;
	
	if (q == NIL || r == NIL) {
		Srfs("div_rem", "args NIL");
		return(ERROR);
		}
	if (n == 0) {
		Srfs("div_rem", "n == 0L");
		return(ERROR);
		}
	d = m / n;
	*q = d;
	*r = m - n * d;
	return(OK);
}

#if TEXDOCU
INT Inverse_mod(INT a, INT p)
#else
Computes $a^{-1}$ mod $p$, i.~e. an integer $b$ such that 
$a b \equiv 1$ mod $p$ holds.
#endif
{
	INT b;

	inverse_mod_integer(a, p, &b);
	return b;
}

#if TEXDOCU
INT inverse_mod_integer(INT a, INT m, INT *b)
#else
Computes $a^{-1}$ mod $p$ into $b$, i.~e. an integer $b$ such that 
$a b \equiv 1$ mod $p$ holds.
#endif
/* Berechnet Inverses von a mod m nach b. 
 * Rueckgabe ERROR, falls a mod m nicht invertierbar. 
 * Wenn m == 0L, so sind nur +-1 invertierbar. */
{
	INT u, v, g;
	BYTE str[256];
	
	if (m == 0L) {
		if (ABS(a) == 1) {
			*b = a;
			}
		else {
			return(ERROR);
			}
		return(OK);
		}
	if (bezout_integer(a, m, &u, &v, &g) != OK) {
		Srff("inverse_mod_integer", "bezout_integer");
		return(ERROR);
		}
	if (g != 1) { /* g != -1, da bezout_integer() normiert */
		sprintf(str, "notinvertible: "
			"a=%ld m=%ld u=%ld v=%ld g=%ld", 
			a, m, u, v, g);
		Srfs("inverse_mod_integer", str);
		return(ERROR);
		}
	*b = u;
	return(OK);
}

#if TEXDOCU
INT ny2(INT p, INT *q, INT *n)
#else
computes ny\_p$(p, 2)$ into $n$. Computes $q := \frac{p}{2^n}$. 
#endif
#if 0
Zaehlt endende Nullen in der 
Binaerdarstellung von $p$ nach $n$;
$q$ enthaelt $p$, nachdem der 2-Anteil abdividiert wurde
(Rechtsshift von $p$ um $n$ Stellen).
$p$ kann auch $< 0$ sein; $q$ ist dann ebenfalls $< 0$. 
q kann mit p uebereinstimmen. 
$p$ muss $\neq 0$ sein.
#endif
{
	INT p1 = p;
	INT n1;
	INT f_negative;
	
	n1 = 0;
	if (p1 == 0) {
		return error("ny2() p == 0");
		}
	if (p1 < 0) {
		p1 = -p1;
		f_negative = TRUE;
		}
	else
		f_negative = FALSE;
	while (TRUE) {
		/* solange p1 kongruent 0 mod 2: */
		if (p1 % 2) /* inkongruent */
			break;
		n1++;
		p1 >>= 1;
		}
	if (f_negative)
		p1 = -p1;
	*q = p1;
	*n = n1;
	return OK;
}

#if TEXDOCU
INT nb_abelian_groups_of_order(INT n)
#else
Computes the number of abelian groups of order $n$.
#endif
{
	VECTOR_OB vp, ve;
	INT k, l, i, j, a;
	VECTOR_OB V;
	
	factor_integer(n, &vp, &ve);
	k = 1;
	l = vp.s_li();
	for (i = 0; i < l; i++) {
		j = ve.s_ii(i);
		vector_of_part(&V, j);
		a = V.s_li();
		k *= a;
		}
	return k;
}

#if TEXDOCU
INT n_Choose_k_first(INT *choice, INT n, INT k)
#else
Computes the first $k$-subset of $\{1,\ldots,n\}$ into choice.
This is always the set $\{0,1,\ldots,k-1\}$.
#endif
{
	INT i;
	
	for (i = 0; i < k; i++) {
		choice[i] = i;
		}
	return TRUE;
}

#if TEXDOCU
INT n_Choose_k_next(INT *choice, INT n, INT k)
#else
Computes the next $k$-subset of $\{1,\ldots,n\}$ into choice.
Returns TRUE if there is another $k$-subset, FALSE otherwise.
#endif
{
	INT i, ii, a;
	
	for (i = 0; i < k; i++) {
		a = choice[k - 1 - i];
		if (a < n - 1 - i) {
			choice[k - 1 - i] = a + 1;
			for (ii = i - 1; ii >= 0; ii--) {
				choice[k - 1 - ii] = choice[k - 1 - ii - 1] + 1;
				}
			return TRUE;
			}
		}
	return FALSE;
}

#if TEXDOCU
INT quadratic_residues(INT p, VECTOR_OP Q, VECTOR_OP N)
#else
Computes the quadratic residues mod $p$ into $Q$, 
the non-quadratic residues into $N$.
#endif
{
	INT ql, nl, i, r;
	
	Q->m_il(0);
	N->m_il(0);
	ql = 0;
	nl = 0;
	for (i = 1; i < p; i++) {
		r = Jacobi(i, p);
		if (r == 1) {
			Q->inc();
			Q->m_ii(ql, i);
			ql++;
			}
		else if (r == -1) {
			N->inc();
			N->m_ii(nl, i);
			nl++;
			}
		else return error("quadratic_residues() r != 1 && r != 0");
		}
	return OK;
}

#if TEXDOCU
INT Jacobi(INT a, INT m)
#else
Computes the Jacobi symbol $\left( \frac{a}{m} \right)$.
#endif
{
	INT a1, m1, ord2, r1;
	INT u, v, g;
	INT f_negative;
	INT t, t1, t2;
	
	a1 = a;
	m1 = m;
	r1 = 1;
	bezout_integer(a1, m1, &u, &v, &g);
	if (ABS(g) != 1) {
		return 0;
		}
	while (TRUE) {
		/* Invariante: 
		 * r1 enthaelt bereits ausgerechnete Faktoren.
		 * ABS(r1) == 1.
		 * Jacobi(a, m) = r1 * Jacobi(a1, m1) und ggT(a1, m1) == 1. */
		if (a1 == 0) {
			return error("Jacobi() a1 == 0");
			}
		a1 = NormRemainder(a1, m1);
		if (a1 < 0)
			f_negative = TRUE;
		else
			f_negative = FALSE;
		ny2(a1, &a1, &ord2);
		/* a1 jetzt immer noch != 0 */
		if (f_negative) {
			t = (m1 - 1) >> 1; /* t := (m1 - 1) / 2 */
			/* Ranmultiplizieren von (-1) hoch t an r1: */
			if (t % 2)
				r1 = -r1; /* Beachte ABS(r1) == 1 */
			/* und a1 wieder positiv machen: */
			a1 = -a1;
			}
		if (ord2 % 2) {
			/* tue nur dann etwas, wenn ord2 ungerade */
			t = (m1 * m1 - 1) >> 3; /* t = (m1 * m1 - 1) / 8 */
			/* Ranmultiplizieren von (-1) hoch t an r1: */
			if (t % 2)
				r1 = -r1; /* Beachte ABS(r1) == 1L */
			}
		if (ABS(a1) <= 1)
			break;
		/* Reziprozitaet: */
		t1 = (m1 - 1) >> 1; /* t1 = (m1 - 1) / 2 */
		t2 = (a1 - 1) >> 1; /* t1 = (a1 - 1) / 2 */
		if ((t1 % 2) && (t2 % 2)) /* t1 und t2 ungerade */
			r1 = -r1; /* Beachte ABS(r1) == 1 */
		t = m1;
		m1 = a1;
		a1 = t;
		}
	if (a1 == 1) {
		return r1;
		}
	if (a1 <= 0) {
		return error("Jacobi() a1 == -1 || a1 == 0");
		}
	return error("Jacobi() wrong termination");
}

#if TEXDOCU
INT NormRemainder(INT a, INT m)
#else
absolute smallest remainder: Computes $r$ such that 
$a \equiv r$ mod $m$ and $- \frac{m}{2} < r \le \frac{m}{2}$ holds.
#endif
{
	INT q, m0, m1, m_halbe;
	
	if (m == 0) {
		return error("NormRemainder() m == 0");
		}
	m0 = m;
	m_halbe = m0 >> 1;
	q = a / m0;
	m1 = a - q * m0;
	if (m1 > m_halbe)
		m1 -= m0;
	if (ODD(m0)) {
		if (m1 < - m_halbe)
			m1 += m0;
		}
	else {
		if (m1 <= - m_halbe)
			m1 += m0;
		}
	return m1;
}

#if TEXDOCU
INT n_choose_k(INT n, INT k)
#else
Computes ${n \choose k}$ as an integer. 
No overflow checking.
#endif
{
	INT a, b;
	
	if (n < 0)
		return error("n_choose_k(): n < 0");
	if (k < 0)
		return error("n_choose_k(): k < 0");
	if (k > n)
		return error("n_choose_k(): k > n");
	if (k == 0)
		return 1;
	if (k == 1)
		return n;
	if (n == 1)
		return 1;
	a = n * n_choose_k(n - 1, k - 1);
	b = a / k;
	if (b * k != a)
		return error("n_choose_k(): b * k != a");
	return b;
}

#if TEXDOCU
INT N_choose_K(SYM_OP n, INT k, SYM_OP res)
#else
Computes ${n \choose k}$ into res as an {\em object}.
This function uses a recursion formula.
#endif
{
	SYM_OB n1, res1, k1, tmp;
	
	if (k < 0)
		return error("N_choose_K(): k < 0");
	k1.m_i_i(k);
	if (k1.gt(n)) {
		res->m_i_i(0);
		return OK;
		}
	if (k == 0) {
		res->m_i_i(1);
		return OK;
		}
	if (k == 1) {
		n->copy(res);
		return OK;
		}
	if (n->einsp()) {
		res->m_i_i(1);
		return OK;
		}
	n->copy(&n1);
	n1.dec();
	N_choose_K(&n1, k - 1, &res1);
	res1.mult(n, &tmp);
	// a = n * n_choose_k(n - 1, k - 1);
	k1.m_i_i(k);
	tmp.ganzdiv(&k1, res);
	// b = a / k;
	return OK;
}

#if TEXDOCU
INT Binomial(INT n, INT k, SYM_OP n_choose_k)
#else
Computes binomial coefficients as {\em objects} 
so that large numbers are no problem. 
This function uses an internal table to remember all 
previously computed values. This may speed up 
those computations where Binomial() is heavily involved !
#endif
{
	static MATRIX_OP T = NIL;
	INT tn, i, j;
	
	if (k < 0) {
		return error("Binomial(): k < 0");
		}
	if (k > n) {
		n_choose_k->m_i_i(0);
		return OK;
		}
	if (n == 0 || k == 0 || k == n) {
		n_choose_k->m_i_i(1);
		return OK;
		}
	if (T == NIL) {
		T = (MATRIX_OP) callocobject("Binomial() table of binomial coefficients");
		T->m_ilih_n(10, 10);
		}
	tn = T->s_hi();
	if (tn < n + 1) {
		MATRIX_OB TT;

#if 0
		printf("reallocating table of binomial coefficients to length %ld\n", n + 1);
		fflush(stdout);
#endif
		TT.m_ilih_n(n + 1, n + 1);
		for (i = 0; i < tn; i++) {
			for (j = 0; j <= i; j++) {
				T->s_ij(i, j)->swap(TT.s_ij(i, j));
				}
			}
		TT.swap(T);
		}
	return Binomial_using_table(n, k, T, n_choose_k);
}

static INT Binomial_using_table(INT n, INT k, MATRIX_OP T, SYM_OP res)
{
	SYM_OB tmp1, tmp2;
	INT tn;
	
	tn = T->s_hi();
	if (tn < n)
		return error("Binomial_using_table: tn < n\n");
	if (k > n) {
		return error("Binomial_using_table: k > n\n");
		}
	if (k > (n >> 1) + 1) {
		Binomial_using_table(n, n - k, T, res);
		T->s_ij(n, n - k)->copy(T->s_ij(n, k));
		return OK;
		}
	if (n == 0) {
		error("Binomial_using_table: n == 0\n");
		return 0;
		}
	if (n < 0) {
		error("Binomial_using_table: n < 0\n");
		return 0;
		}
	if (k < 0) {
		error("Binomial_using_table: k < 0\n");
		return 0;
		}
	if (n == k) {
		T->m_iji(n, k, 1);
		res->m_i_i(1);
		return OK;
		}
	if (k == 0) {
		T->m_iji(n, k, 1);
		res->m_i_i(1);
		return OK;
		}
	if (k == 1) {
		T->m_iji(n, k, n);
		res->m_i_i(n);
		return OK;
		}
	if (T->s_iji(n, k) == 0) {
		Binomial_using_table(n - 1, k - 1, T, &tmp1);
		Binomial_using_table(n - 1, k, T, &tmp2);
		tmp1.add(&tmp2, res);
		res->copy(T->s_ij(n, k));
		return OK;
		}
	T->s_ij(n, k)->copy(res);
	return OK;
}

#if TEXDOCU
INT hamming_bound_q(INT n, INT t, INT q, INTEGER_OP res, INT scale)
#else
Computes the Hamming bound (upper bound) for $GF(q)$-codes:
Given $n, t$ and $q$, we compute the maximal 
information rate $k/n$ which a code can have.
The bound is:
\[
B = \frac{1}{n} \cdot log_q(\frac{q^n}{(\sum_{i=0}^t {n \choose i} )})
\]
The bound is multiplied by scale and then returned in res.
#endif
{
	double d3, d4, d5;
	SYM_OB a, b, c;

	a.m_i_i(i_power_j(q, n));
	binomial_sum(n, t, &b);
	quotient_scaled(&a, &b, &c, scale);
	d3 = (double) c.s_i_i() / (double) scale;
	d4 = log(d3) / log((double)q);
	d4 /= (double) n;
	d5 = d4 * scale;
	res->m_i((INT) d5);
	return OK;

#if 0
	d1 = (double) (1L << n);
	/* printf("(%lf / ", d1); */
	d2 = (double) binomial_sum(n, t);
	/* printf("%lf = ", d2); */
	d3 = d1 / d2;
	/* printf("%lf)", d3); */
	d4 = log(d3) / log(2.);
	return d4 / (double)n;
#endif
}

#if TEXDOCU
INT gilbert_varshamov_bound_q(INT n, INT t, INT q, 
	INTEGER_OP res, INT scale)
#else
Computes the Gilbert-Varshamov lower bound 
for the information-rate $k/n$ of a code with 
given parameters $n, t$ and $q$.
(the bound tells that a code with information rate at least $k/n$ 
must exist). The bound is:
\[
B = \frac{1}{n} \cdot log_q(\frac{q^n}{(\sum_{i=0}^{d-2} {n-1 \choose i} )}) 
\]
where $d = 2t + 1$. ($t$ is the error correction capability).
#endif
{
	double d3, d4, d5;
	SYM_OB a, b, c;
	INT d;

	d = (t << 1) + 1;
	a.m_i_i(i_power_j(q, n));
	binomial_sum(n - 1, d - 2, &b);
	quotient_scaled(&a, &b, &c, scale);
	d3 = (double) c.s_i_i() / (double) scale;
	d4 = log(d3) / log((double)q);
	d4 /= (double) n;
	d5 = d4 * scale;
	res->m_i((INT) d5);
	return OK;
#if 0
	d = (t << 1) + 1;
	d1 = (double) (1L << n);
	/* printf("(%lf / ", d1); */
	d2 = (double) binomial_sum(n - 1, d - 2);
	/* printf("%lf = ", d2); */
	d3 = d1 / d2;
	printf("%lf)", d3);
	d4 = log(d3) / log(2.);
	return d4 / (double) n;
#endif
}

#if TEXDOCU
INT binomial_sum(INT n, INT t, SYM_OP res)
#else
Computes $\sum_{i=0}^t {n \choose i}$.
#endif
{
	INT l = 0;
	INT i;
	SYM_OB tmp1, tmp2, tmp3;
	
	tmp1.m_i_i(0);
	for (i = 0; i <= t; i++) {
		Binomial(n, i, &tmp2);
		tmp1.add(&tmp2, &tmp3);
		tmp3.swap(&tmp1);
		}
	tmp1.swap(res);
	return l;
}

#if TEXDOCU
INT quotient_scaled(SYM_OP a, SYM_OP b, SYM_OP c, INT scale)
#else
#endif
{
	SYM_OB a0, a1;
	
	a0.m_i_i(scale);
	a0.ganzdiv(b, c);
	return OK;
}


#if TEXDOCU
INT inpsl(INT n, INT d)
#else
#endif
{
	if ((n * (d - 1) / d) % 2 == 0)
		return 1;
	else
		return 0;
}

#if TEXDOCU
INT Binomial1(INT n, INT k, SYM_OP res)
#else
#endif
{
	if (n < 0) {
		res->m_i_i(0);
		return OK;
		}
	if (k < 0) {
		res->m_i_i(0);
		return OK;
		}
	return Binomial(n, k, res);
}

#if TEXDOCU
INT nb_orbits_psl_pgl(INT q, INT k, INT f_v, 
	SYM_OP nb_psl_orbits, SYM_OP nb_pgl_orbits)
#else
#endif
{
	SYM_OB sum, sum1, sum2, sum3, sum4;
	SYM_OB rest1, rest2, rest3, rest4;
	SYM_OB psl, psl1, psl2, psl3, psl4;
	SYM_OB pslrest1, pslrest2, pslrest3, pslrest4;
	SYM_OB c1, c2, c3, c4, s1, s2, s3, s4, S, t1, t2, t3, t4, T;
	INTEGER_OB m_ob, q_ob, two_ob;
	SYM_OB a, b, c, aa;
	SYM_OB qq, qqm1, qqmq, qqpq;
	INT r, s, d, m;
	INT i, p, f;
	INT *euler;
	
	factor_prime_power(q, &p, &f);
	if (f_v) {
		printf("q=%ld=%ld^%ld k=%ld\n", q, p, f, k);
		}
	
	euler = (INT *) my_malloc((q + 2) * sizeof(INT), "euler");
	for (i = 1; i < q + 2; i++) 
		euler[i] = eulerfunc(i);
	
	if (k <= 3) {
		nb_psl_orbits->m_i_i(1);
		nb_pgl_orbits->m_i_i(1);
		return OK;
		}
	if (q + 1 < k) {
		nb_psl_orbits->m_i_i(0);
		nb_pgl_orbits->m_i_i(0);
		return OK;
		}

	if (f_v) {
		printf("step1\n"); fflush(stdout);
		}
	Binomial1(q - 2, k - 3, &sum1); // sum1 = choose(q-2,k-3);
	m = k * (k - 1) * (k - 2);
	m_ob.m_i(m);
	q_ob.m_i(q);
	two_ob.m_i(2);
	sum1.quores(&m_ob, &sum, &rest1); // sum = sum1/m; rest1 = sum1%m;
	sum1.mult(&two_ob, &a);
	a.quores(&m_ob, &psl, &aa);  // psl = 2 * sum1/m;
	two_ob.mult(&sum1, &a);
	a.quores(&m_ob, &aa, &pslrest1); //  pslrest1 = 2*sum1%m;
	if (f_v) {
		printf("rest1=");rest1.println();
		printf("pslrest1=");pslrest1.println();
		}
	
	if (k % p == 0)
		r = k / p;
	else 
		r = -1;
	if ((k - 1) % p == 0)
		s = (k - 1) / p;
	else
		s = -1;
	Binomial1(q / p, r, &a);
	Binomial1(q / p, s, &b);
	a.add(&b, &sum2);   // sum2 = ( choose(q/p,r) + choose(q/p,s));
	sum2.quores(&q_ob, &a, &rest2); // rest2 = sum2%q;
	a.add_apply(&sum); // sum = sum + sum2/q;
	two_ob.mult(&sum2, &b);
	b.quores(&q_ob, &c, &aa);
	c.add_apply(&psl); // psl = psl + 2 * sum2/q;
	c.copy(&pslrest2); // pslrest2 = 2* sum2%q;
	if (f_v) {
		printf("rest2=");rest2.println();
		printf("pslrest2=");pslrest2.println();
		printf("sum=");sum.println();
		printf("psl=");psl.println();
		}

	if (f_v) {
		printf("step2\n"); fflush(stdout);
		}
	sum3.m_i_i(0);
	psl3.m_i_i(0);
	for (d = 2; d <= q + 1; d++) {
		if ((q + 1) % d == 0) {
			if (k % d == 0) {
				Binomial1((q + 1) / d, k / d, &a);
				c.m_i_i(euler[d]);
				c.mult(&a, &aa);
				aa.add_apply(&sum3); // sum3 = sum3 + euler[d] *  choose((q+1)/d,k/d);
				if (inpsl(q + 1, d)) {
					aa.add_apply(&psl3);
					}
				}
			} // if
		} // next d
	a.m_i_i(2 * q + 2);
	sum3.quores(&a, &b, &rest3); // rest3 = sum3%(2*q + 2);
	b.add_apply(&sum); // sum = sum + sum3/(2*q + 2);
	a.m_i_i(q + 1);
	psl3.quores(&a, &b, &pslrest3); // pslrest3 = psl3%(q + 1);
	b.add_apply(&psl); // psl = psl + psl3/(q + 1);
	if (f_v) {
		printf("rest3=");rest3.println();
		printf("pslrest3=");pslrest3.println();
		printf("sum=");sum.println();
		printf("psl=");psl.println();
		}
	

	if (f_v) {
		printf("step3\n"); fflush(stdout);
		}
	sum4.m_i_i(0);
	psl4.m_i_i(0);
	for (d = 2; d < q; d++) {
		if ((q - 1) % d == 0) {
			r = (q - 1) / d;
			if (k % d == 0) { // fixpoint which consists only of d-cycles
				Binomial(r, k / d, &a);
				c.m_i_i(euler[d]);
				a.mult(&c, &aa);
				aa.add_apply(&sum4);
				if (inpsl(q - 1, d))
					aa.add_apply(&psl4);
				}
			if ((k - 1) % d == 0) { // fixpoint with one 1-cycle and then only d-cycles
				Binomial(r, (k - 1) / d, &a);
				two_ob.mult_apply(&a); // two choices for the fixpoint 
				c.m_i_i(euler[d]);
				a.mult(&c, &aa);
				aa.add_apply(&sum4);
				if (inpsl(q - 1, d))
					aa.add_apply(&psl4);
				}
			if ((k - 2) % d == 0) { // fixpoint with two 1-cycles and then only d-cycles
				Binomial(r, (k - 2) / d, &a);
				c.m_i_i(euler[d]);
				a.mult(&c, &aa);
				aa.add_apply(&sum4);
				if (inpsl(q - 1, d))
					aa.add_apply(&psl4);
				}
			} // if
		} // next d
	
	if (f_v) {
		printf("step4\n"); fflush(stdout);
		}
	a.m_i_i(2 * q - 2);
	sum4.quores(&a, &b, &rest4); // rest4 = sum4%(2*q - 2);
	b.add_apply(&sum); // sum =  sum + sum4/(2*q - 2);

	a.m_i_i(q - 1);
	psl4.quores(&a, &b, &pslrest4); // pslrest4 = psl4%(q - 1);
	b.add_apply(&psl); // psl =  psl + psl4/(q - 1);
	if (f_v) {
		printf("rest4=");rest4.println();
		printf("pslrest4=");pslrest4.println();
		printf("sum=");sum.println();
		printf("psl=");psl.println();
		}
	
	q_ob.mult(&q_ob, &qq);
	b.m_i_i(-1);
	qq.add(&b, &qqm1);
	b.m_i_i(-q);
	qq.add(&b, &qqmq);
	qq.add(&q_ob, &qqpq);
	b.m_i_i(-q);
	qq.add(&b, &qqmq);
	
	qqm1.mult(&q_ob, &c1);
	two_ob.mult_apply(&c1); // c1 = 2 * q * (q * q - 1)
	
	qqm1.mult(&m_ob, &c2);
	two_ob.mult_apply(&c2); // c2 = 2 * (q * q - 1) * m

	qqmq.mult(&m_ob, &c3); // c3 = (q * q - q) * m

	qqpq.mult(&m_ob, &c4); // c4 = (q * q + q ) * m

	if (f_v) {
		printf("c1=");c1.println();
		printf("c2=");c2.println();
		printf("c3=");c3.println();
		printf("c4=");c4.println();
		}
	
	c1.mult(&rest1, &s1); // s_i = c_i * rest_i
	c2.mult(&rest2, &s2);
	c3.mult(&rest3, &s3);
	c4.mult(&rest4, &s4);
	
	c1.mult(&pslrest1, &t1); // t_i = c_i * pslrest_i
	c2.mult(&pslrest2, &t2);
	c3.mult(&pslrest3, &t3);
	c4.mult(&pslrest4, &t4);
	
	s1.add(&s2, &a);
	a.add(&s3, &b);
	b.add(&s4, &S); // S = sum_{i=1}^4 c_i rest_i
	
	t1.add(&t2, &a);
	a.add(&t3, &b);
	b.add(&t4, &T); // T = sum_{i=1}^4 c_i pslrest_i
	
	S.ganzdiv(&two_ob, &a);
	a.ganzdiv(&q_ob, &b);
	b.ganzdiv(&m_ob, &c);
	c.ganzdiv(&qqm1, &aa);
	aa.add_apply(&sum);
	
	T.ganzdiv(&q_ob, &a);
	a.ganzdiv(&m_ob, &b);
	b.ganzdiv(&qqm1, &aa);
	aa.add_apply(&psl);
	
	sum.copy(nb_pgl_orbits);
	psl.copy(nb_psl_orbits);
	
#if 0
 sum = sum +
  ( 2 * q * (q * q - 1) * rest1 +
    2 * (q * q - 1) * m * rest2 +
    (q * q - q ) * m *    rest3 +
    (q * q + q ) * m *    rest4      )/2/q/m/(q*q - 1);

 psl = psl +
  ( 2 * q * (q * q - 1) * pslrest1 +
    2 * (q * q - 1) * m * pslrest2 +
    (q * q - q ) * m *    pslrest3 +
    (q * q + q ) * m *    pslrest4   )/q/m/(q*q - 1);
    
#endif
	my_free(euler);
	return OK;
}

#if TEXDOCU
INT multinomial(PARTITION_OP p, SYM_OP res, INT f_v)
#else
#endif
{
	INT i, j, l, w;
	SYM_OB a, b, c, d;
	
	l = p->s_li();
	w = p->weight_i();
	d.fakul_int(w);
	// printf("w=%ld d=",w);
	// d.println();
	b.m_i_i(1);
	for (i = 0; i < l; i++) {
		j = p->s_ii(i);
		a.fakul_int(j);
		a.mult(&b, &c);
		c.swap(&b);
		}
	d.ganzdiv_integral(&b, res);
	if (f_v) {
		printf("multinomial ");
		p->print();
		printf(" = ");
		res->println();
		}
	return OK;
}

#if TEXDOCU
INT multinomial_ordered(PARTITION_OP p, SYM_OP res, INT f_v)
#else
#endif
{
	VECTOR_OB m;
	PARTITION_OB q, qq;
	INT i, j, w;
	SYM_OB a, b, c, d;
	
	w = p->weight_i();
	d.fakul_int(w);
	p->t_VECTOR_EXPONENT(&q);
	multinomial(p, &b, f_v);
	m.m_il(0);
	for (i = 1; i <= w; i++) {
		j = q.exp_nb_i_parts(i);
		if (j == 0)
			continue;
		m.inc();
		m.m_ii(m.s_li() - 1, j);
		}
	m.quicksort(m.s_li(), TRUE);
	qq.m_ks(VECTOR, &m);
	if (f_v) {
		printf("exponent multiplicities: ");
		qq.println();
		}
	multinomial(&qq, &a, f_v);
	a.mult(&b, res);
	return OK;
}

#if TEXDOCU
INT stirling_second(INT i, INT j, INT f_ordered, SYM_OP res, INT f_v)
#else
#endif
{
	SYM_OB a, b, c;
	PARTITION_OB part;
	
	if (i == 0) {
		if (j == 0) {
			res->one();
			return OK;
			}
		else {
			res->zero();
			return OK;
			}
		}
	if (j == 0) {
		res->zero();
		return OK;
		}
	if (j > i) {
		res->zero();
		return OK;
		}
	if (f_v) {
		printf("stirling_second");
		if (f_ordered) {
			printf(" ordered ");
			}
		printf("(%ld, %ld)\n", i, j);
		}
	a.m_i_i(0);
	part.first_into_k_parts_VECTOR(i, j);
	do {
		if (f_v) {
			part.println();
			}
		multinomial_ordered(&part, &b, f_v);
		a.add(&b, &c);
		c.swap(&a);
		} while (part.next_into_k_parts_VECTOR(i, j));
	if (!f_ordered) {
		b.fakul_int(j);
		a.ganzdiv_integral(&b, &c);
		c.swap(&a);
		}
	a.swap(res);
	if (f_v) {
		printf("stirling_second");
		if (f_ordered) {
			printf(" ordered ");
			}
		printf("(%ld, %ld) = ", i, j);
		res->println();
		}
	return OK;
}

#if TEXDOCU
INT stirling_first(INT i, INT j, INT f_signless, SYM_OP res, INT f_v)
#else
#endif
{
	INT l, x, ax;
	SYM_OB a, b, c, d;
	PARTITION_OB part;
	INT f_ordered = FALSE;
	
	if (i == 0) {
		if (j == 0) {
			res->one();
			return OK;
			}
		else {
			res->zero();
			return OK;
			}
		}
	if (j == 0) {
		res->zero();
		return OK;
		}
	if (j > i) {
		res->zero();
		return OK;
		}
	if (f_v) {
		printf("stirling_first");
		printf("(%ld, %ld)\n", i, j);
		}
	a.m_i_i(0);
	part.first_into_k_parts_VECTOR(i, j);
	do {
		if (f_v) {
			part.println();
			}
		multinomial_ordered(&part, &b, f_v);
		l = part.s_li();
		for (x = 0; x < l; x++) {
			ax = part.s_ii(x);
			if (ax == 1)
				continue;
			c.fakul_int(ax - 1);
			b.mult(&c, &d);
			d.swap(&b);
			}
		a.add(&b, &c);
		c.swap(&a);
		} while (part.next_into_k_parts_VECTOR(i, j));
	if (!f_ordered) {
		b.fakul_int(j);
		a.ganzdiv_integral(&b, &c);
		c.swap(&a);
		}
	if (!f_signless) {
		if (ODD(i + j))
			a.addinvers_apply();
		}
	a.swap(res);
	if (f_v) {
		printf("stirling_first");
		printf("(%ld, %ld) = ", i, j);
		res->println();
		}
	return OK;
}


