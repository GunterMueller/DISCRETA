/* vec.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef VEC_INCLUDED
#define VEC_INCLUDED

#ifndef IN_INCLUDED
#include <DISCRETA/in.h>
#endif

#define VECTOR_SAFE
// range check in s_i() !



class vector_ob : public SYM_OB {
public:
	INT freeself();
	INT freeself_debug(INT print_depth);
	// access vectorstruct:
	// before: S_V_S, C_V_S, S_V_L, C_V_L, S_V_LI
	SYM_OP s_s() { 
		return((SYM_OP)(ob_self.ob_vector->v_self)); };
	void c_s(SYM_OP self) { 
		(ob_self.ob_vector)->v_self = (OP)self; };
	INTEGER_OP s_l() { 
		return((INTEGER_OP)ob_self.ob_vector->v_length); }
	void c_l(INTEGER_OP l) { 
		ob_self.ob_vector->v_length = (OP)l; };
	INT s_li() { return(s_l()->s_i()); };

	// access vector elements:
	// before: S_V_I, S_V_II
	SYM_OP s_i(INT i);
		/* i-tes Element zurueckgeben. */
	INT s_ii(INT i) { return(((INTEGER_OP)s_i(i))->s_i()); };
		/* i-tes Element als Integer zurueckgeben. */
	void c_i(INT i, SYM_OP b) {
		s_i(i)->freeself(); *s_i(i) = *b; };
		/* i-tes Element auf b setzen, 
		 * b wird nicht kopiert, sondern Bestandteil des Vektors. */
	void c_ii(INT i, INT j) { 
		((INTEGER_OP)s_i(i))->c_i(j); };
		/* i-tes Element auf j setzen, ob_kind unveraendert lassen. */
	void m_ii(INT i, INT j) { 
		((INTEGER_OP)s_i(i))->m_i(j); };
		/* i-tes Element auf j setzen, ob_kind auf INTEGER setzen. */
	
	SYM_OB &operator[](int i) { return *s_i(i); };
	
	// vector creation:
	INT m_il(INT l);
	INT m_il_n(INT l);
		/* mit 0 vorbesetzen */
	INT m_l(INTEGER_OP l) { 
		return(m_il(l->s_i())); };
		/* make_length_vector */
	INT m_l_n(INTEGER_OP l) { 
		return(m_il_n(l->s_i())); };
	INT b_ls(INTEGER_OP length, SYM_OP self);
	
	// c_v_i, m_il_v, m_il_nv, m_o_v, b_l_v, b_l_nv, b_o_v, b_ls_v

INT Print();
INT latex(FILE *fp);
INT Latex(FILE *fp);
INT fprint_GAP(FILE *fp);
INT sprint(BYTE *str);
/* appends to str. writes to maximal strlength of 200. */
INT sprint_len(INT len, BYTE *str);
INT realloc_z(INT new_len);
INT realloc(INTEGER_OP new_len);
INT append_i(INT i);
INT append_element(INTEGER_OP len, SYM_OP b);
/* Element b an Vektor der Laenge len anfuegen (len ist Anzahl benutzter 
 * Elemente im Vektor, nicht dessen tatsaechliche Laenge). Es wird also 
 * b zum len-ten Element des Vektors. Gegebenenfalls wird der Vektor 
 * via realloc() verlaengert, und zwar dann gleich um VECTOR_OVERSIZE viele
 * Elemente. Der Zaehler len wird incrementiert. */
INT append_element_itself(INTEGER_OP len, SYM_OP b, INT f_itself);
/* f_itself TRUE: b wird per swap() an seinen Platz im Array bewegt. 
 * f_itself FALSE: via copy(). */
INT insert_at(INTEGER_OP len, INT i, SYM_OP b);
/* Fuegt Element b an i-ter Stelle im Vektor ein. Der Vektor wird 
 * gegebenenfalls verlaengert. Der Zaehler len wird incrementiert. */
/* i <= len notwendig;
 * i == len heisst es wird angefuegt. */
INT insert_at_itself(INTEGER_OP len, INT i, SYM_OP b, INT f_itself);
/* f_itself TRUE: b wird per swap() an seinen Platz im Array bewegt. 
 * f_itself FALSE: via copy(). */
INT delete_ith(INTEGER_OP len, INT i);
/* len wird decrementiert. */
INT v_del_ith2(INT len, INT i);
INT v_insert_at(INT i, SYM_OP b);
INT v_realloc(INT len);
INT v_shorten(INT len);
INT v_minus_v_ip(VECTOR_OP q);
/* in place version (p := p - q, p is this argument). */
INT v_minus_v(VECTOR_OP q, VECTOR_OP p_minus_q)
/* q must be sorted ascendingly (p is this).*/;
INT do_search(
	INT len, INT f_ascending, SYM_OP v, 
	INT *idx, INT *f_found, 
	INT type, void *data);
INT search(INT len, INT f_ascending, SYM_OP v, 
	INT *idx, INT *f_found);
/* Sucht im Array der Laenge len (aufsteigend sortiert falls f_ascending) 
 * nach v. Gibt in idx die Position des ersten Elements
 * zurueck, welches echt groesser v ist (Beachte: idx kann == len sein).
 * f_found, falls dessen Vorgaenger gleich v ist. 
 * Beachte: wenn f_found TRUE ist, so ist des gefundene Element 
 * das (idx-1)-te ! */
INT insert_sorted(INTEGER_OP len, INT f_ascending, SYM_OP v);
INT search_and_insert(INT len, SYM_OP v);
INT search_and_insert_int(INT len, INT i);
INT set_minus(VECTOR_OP B, INTEGER_OP len_A, 
	INTEGER_OP len_B, INT f_ascending);
/* A := A - B (A := this) */
INT set_minus2(SYM_OP B, INTEGER_OP len_A, 
	INT len_B, INT f_ascending);
/* former name: vs_set_minus(VECTOR_OP A, OP B, ... ) */
INT subseteq(VECTOR_OP B, INT lenA, INT lenB, 
	INT f_ascending, INT *f_subseteq);
/* Test, ob A Teilmenge B (A := this). */
INT equal(VECTOR_OP B, INT lenA, INT lenB, 
	INT f_ascending, INT *f_equal);
/* Test, ob A = B. */
/* Annahme: keine wiederholten Eintraege in den Vektoren. */

INT quicksort(INT len, INT f_ascending);
/* in iof.C: */
INT hip();
INT hip1();
INT calc_size_on_file();
INT write_mem(MEM_OP mem, INT debug_depth);
INT read_mem(MEM_OP mem, INT debug_depth);

INT nullp();
INT einsp();
INT vectorp();
INT m_o(SYM_OP ob);
INT b_o(SYM_OP ob);
INT add_apply(VECTOR_OP b);
/* b = b + a */
INT add(VECTOR_OP b, VECTOR_OP c);
INT addinvers(VECTOR_OP erg);
INT addinvers_apply();
INT addtoallelements(SYM_OP zahl, VECTOR_OP ergebnis);
/* before: addtoallvectorelements(OP zahl, OP vector, OP ergebnis) */
INT copy(VECTOR_OP res);
INT compare(VECTOR_OP b);
INT lastof(SYM_OP b);
INT length(INTEGER_OP b);
INT inc();
/* verlaengert den vector um ein leeres object, welches am Ende steht */
/* dabei werden die vector elemente kopiert */
INT dec();
/* kuerzt den vector um 1 */
/* das letzte element wird gestrichen */
INT append(SYM_OP b, VECTOR_OP c);
/* haengt den vector b (bzw. das object b, 
 * welches zum vector konvertiert wird) an den vector a an */
INT append_in_place(SYM_OP b);
INT max_vector(SYM_OP m);
/* kopiert maximales element */
INT sum(INTEGER_OP ergebnis);
/* berechnet die summe der vectorelemente;
 * nur fuer INTEGER vectoren. */
INT mult_scalar(SYM_OP a, VECTOR_OP c);
/* skalarmultiplikation */
/* before: mult_scalar_vector(OP a, OP b, OP c) */
/*         a ist skalar b ist vector c wird vector ACHTUNG ! */
INT mult_matrix(MATRIX_OP b, VECTOR_OP c);
/* before: mult_vector_matrix(OP a, OP b, OP c) */
INT mult_vector(VECTOR_OP b, VECTOR_OP c);
/* componentwise multiplication */
/* before: mult_vector_vector(OP a, OP b, OP c) */
INT scalarproduct(VECTOR_OP b, SYM_OP d);
/* before: scalarproduct_vector(OP a, OP b, OP d) */
INT mult_apply_vector(VECTOR_OP b);
/* componentwise multiplication */
/* before: mult_apply_vector_vector(OP a, OP b) */
INT mult_apply(SYM_OP b);
/* before: mult_apply_vector(OP a, OP b) */
INT search_linear(SYM_OP p);
INT join2(VECTOR_OP V1, VECTOR_OP V2);
INT multiplicities(VECTOR_OP val, VECTOR_OP mult);
INT classify(VECTOR_OP val, VECTOR_OP classes);
INT classify_and_reorder(VECTOR_OP val, VECTOR_OP classes, 
	PERMUTATION_OP p, VECTOR_OP class_lengths);
INT norm(VECTOR_OP mult);
INT sprint_multiplicities(VECTOR_OP mult, BYTE *str);
INT gcd_all_elements(SYM_OP g);
INT gcd_all_elements_and_divide_out(SYM_OP g, VECTOR_OP V);
INT add_apply_elementwise(VECTOR_OP V2, INT sign);
INT divide_out(SYM_OP val);
INT multiply_elementwise(SYM_OP val);
INT apply_perm(PERMUTATION_OP p);
INT apply_perm_to_vector_of_pairs(PERMUTATION_OP p, INT n);





/* in fga/perm.C: */
INT first_lehmercode(INTEGER_OP l);
/* l beleibt erhalten */
/* firstlemercode = 0000...0000 */
INT last_lehmercode(INTEGER_OP l);
/* lastlehmercode = 0123...n-1 */
INT next_lehmercode(VECTOR_OP next);
/* Erzeugt den lexikographisch naechsten Lehmercode nach next. 
 * Rueckgabe LASTLEHMERCODE, falls Ende erreicht 
 * (next ist dann EMPTY). */
INT lehmercode_perm(PERMUTATION_OP b);
/* before: lehmercode_vector(OP vec, OP b) */
/* diese procedure berechnet aus dem lehmercode vec = [v1,....,vn]
die zugehoerige permutation b [e1,...,en] */
INT rz_lehmercode(VECTOR_OP b);
/* bildet die reduzierte zerlegung des lehmercodes
 * bsp this = 321200 dann ist ergebnis 32132354
 * vgl verfahren 1 in diplomarbeit */
INT m_standard_1_n(INT l);
INT m_standard_0_nm1(INT l);
INT m_intvect1(INT a);
INT m_i_times_j(INT i, INT j);
INT power_elementwise(INT k, VECTOR_OP q);
INT power_elementwise_apply(INT k);
INT first_subset(INT n);
INT next_subset(INT n);
INT first_subset_ordered(INT n);
INT next_subset_ordered(INT n);
INT all_k_subsets(INT n, INT k);
INT first_k_subset(INT n, INT k);
INT next_k_subset(INT n, INT k);
INT is_subset(VECTOR_OP set2);
INT is_fix_under(PERMUTATION_OP p);
INT dec_all_entries();
INT inc_all_entries();
INT canonicize_map(LABRA_OP G, LABRA_OP aut, 
	PERMUTATION_OP transporter, INT f_v);
INT Nabla(LABRA_OP G, POLYNOM_OP res, INT n, INT f_v);
INT nabla(POLYNOM_OP res, INT n, INT f_v);
INT compare_images_as_unordered_multisets(VECTOR_OP V2);
};

INT vec_test();
struct vector *callocvectorstruct();
INT VS_search(
	SYM_OP p, INT len, 
	INT f_ascending, SYM_OP v, 
	INT *idx, INT *f_found, 
	INT type, void *data);
/* p is v_self component of a vector */
INT VS_subseteq(SYM_OP A, SYM_OP B, INT lenA, INT lenB, 
	INT f_ascending, INT *f_subseteq);
/* Test, ob A Teilmenge B. */
INT v_println(VECTOR_OP V, INT f_numerated);
INT v_do_println(VECTOR_OP V, 
	INT f_numerated, INT n_on_a_row, 
	INT type, void *data);
INT v_do_fprintln(VECTOR_OP V, FILE *fp, 
	INT f_numerated, INT n_on_a_row, 
	INT type, void *data);

#endif /* VEC_INCLUDED */

