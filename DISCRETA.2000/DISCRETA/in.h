/* in.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef IN_INCLUDED
#define IN_INCLUDED

class integer_ob : public SYM_OB {
public:
	void freeself();
	void m_i(INT i) { 
		freeself(); // new, AB Mai 1997
		ob_kind = INTEGER;
		ob_self.ob_INT = i; };
	void c_i(INT i) { 
		ob_self.ob_INT = i; };
	INT s_i() { return(ob_self.ob_INT); };
	
/* INT freeself(); */
INT sprint(BYTE *s);
INT sprint_latex(BYTE *s);
INT sscan(BYTE *s);
INT invers_apply();
INT addinvers_apply();
INT addinvers(INTEGER_OP b);
INT inc();
INT dec();
INT mult_integer(INTEGER_OP b, INTEGER_OP d);
/* Vorm. mult_integer_integer(a, b, d) */
INT mult(SYM_OP b, SYM_OP d);
/* before: mult_integer(a, b, d) */
INT even();
INT posp();
INT negp();
INT add_integer(INTEGER_OP b, INTEGER_OP c);
/* before: add_integer_integer(OP a, OP b, OP c) */
INT add(SYM_OP b, SYM_OP c);
/* before: add_integer(OP a, OP b, OP c) */
/* das erste object ist vom typ INTEGER, 
 * das ergebnis ist ein leeres object */
INT comp_integer(INTEGER_OP b);
INT compare(SYM_OP b);
/* a ist vom typ INTEGER, b hat unbekannten typ */
INT nullp();
INT einsp();
INT negeinsp();
INT copy(INTEGER_OP b);
INT invers(INTEGER_OP b);
INT add_apply_integer(INTEGER_OP b);
INT add_apply(SYM_OP b);
INT mult_apply_integer(INTEGER_OP b);
INT mult_apply(SYM_OP b);

INT log_10();
/* anzahl stellen */ /* vorm. intlog(OP a) */
INT zero();
INT one();
INT m_one();
INT homo_z(INT z);
INT ganzdiv(INTEGER_OP b, INTEGER_OP c);
INT quores(INTEGER_OP b, INTEGER_OP c, INTEGER_OP d);
INT modulo(INTEGER_OP b, INTEGER_OP c);
/* before: mod_integer(OP a, OP b, OP c) */
};

/* in1.C: */
INT ord_to_order_if_prime(INT ord1, INT *ord);
INT ord_to_order_if_prime_power(
	INT ord1, INT *ord, INT *prime, INT *k);
INT kgv_iipi(INT m, INT n, INT *t);
INT ggt_iipi(INT m, INT n, INT *t);
/* former name: ggt_llpl(), former name: ggt() */
INT smallest_primedivisor(INT n);
INT sp_ge(INT n, INT p_min);
INT ny_p(INT n, INT p);
INT nb_primes(INT n);
/* gives the number of primes 
 * (counting multiplicities) of n. */
INT is_prime(INT n);
INT is_prime_power(INT n, INT p);
INT factor_prime_power(INT n, INT *p, INT *e);
INT i_power_j(INT i, INT j);
INT eulerfunc(INT n);
INT moebius(INT i);
INT factor_integer(INT n, VECTOR_OP primes, VECTOR_OP exponents);
INT print_factorization( VECTOR_OP primes, VECTOR_OP exponents, BYTE *str);
INT sieve_primes(VECTOR_OP Primes, INTEGER_OP nb_primes, INT *first, INT *len);
INT order_mod_p(INT a, INT p);
INT primitive_root(INT p);
INT bezout_integer(INT m, INT n, INT *u, INT *v, INT *g);
/* Findet den positiven ggT(m, n) mit Vorfaktoren: 
 * g = ggT(m, n) = u * m + v * n */
INT asr(INT a, INT m, INT *r);
INT Asr(INT a, INT m);
INT abs_sm_rem(INT a, INT m, INT *r);
/* Berechnet r als Rest von a mod m: 
 * - m / 2 < r <= m / 2.
 * r kann mit a uebereinstimmen (der Pointer). 
 * Bei m == 0L wird mit Warnung zurueckgekehrt 
 * (Rueckgabe OK) dann: r = a. */
INT div_rem(INT m, INT n, INT *q, INT *r);
/* Division mit Rest von m duch n:
 * m = q * n + r mit 0 <= r < n 
 * Fehler, wenn n == 0 */
INT Inverse_mod(INT a, INT p);
INT inverse_mod_integer(INT a, INT m, INT *b);
/* Berechnet Inverses von a mod m nach b. 
 * Rueckgabe ERROR, falls a mod m nicht invertierbar. 
 * Wenn m == 0L, so sind nur +-1 invertierbar. */
INT ny2(INT p, INT *q, INT *n);
/* Zaehlt endende Nullen in der 
 * Binaerdarstellung von p nach n;
 * q enthaelt p, nachdem der 2-Anteil 
 * abdividiert wurde
 * (Rechtsshift von p um n Stellen).
 * p kann auch < 0 sein; q ist dann ebenfalls < 0. 
 * q kann mit p uebereinstimmen. 
 * p muss != 0 sein. */
INT nb_abelian_groups_of_order(INT n);
INT n_Choose_k_first(INT *choice, INT n, INT k);
INT n_Choose_k_next(INT *choice, INT n, INT k);
INT quadratic_residues(INT p, VECTOR_OP Q, VECTOR_OP N);
INT Jacobi(INT a, INT m);
INT NormRemainder(INT a, INT m);
INT n_choose_k(INT n, INT k);
INT N_choose_K(SYM_OP n, INT k, SYM_OP res);
INT Binomial(INT n, INT k, SYM_OP n_choose_k);
INT hamming_bound_q(INT n, INT t, INT q, INTEGER_OP res, INT scale);
INT gilbert_varshamov_bound_q(INT n, INT t, INT q, 
	INTEGER_OP res, INT scale);
INT binomial_sum(INT n, INT t, SYM_OP res);
INT quotient_scaled(SYM_OP a, SYM_OP b, SYM_OP c, INT scale);
INT inpsl(INT n, INT d);
INT Binomial1(INT n, INT k, SYM_OP res);
INT nb_orbits_psl_pgl(INT q, INT k, INT f_v, 
	SYM_OP nb_psl_orbits, SYM_OP nb_pgl_orbits);
INT multinomial(PARTITION_OP p, SYM_OP res, INT f_v);
INT multinomial_ordered(PARTITION_OP p, SYM_OP res, INT f_v);
INT stirling_second(INT i, INT j, INT f_ordered, SYM_OP res, INT f_v);
INT stirling_first(INT i, INT j, INT f_signless, SYM_OP res, INT f_v);

#endif /* IN_INCLUDED */
