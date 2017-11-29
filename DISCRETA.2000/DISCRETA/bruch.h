/* bruch.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef BRUCH_INCLUDED
#define BRUCH_INCLUDED

class bruch_ob : public SYM_OB {
public:
	INT s_i() { return (ob_self.ob_bruch->b_info); };
	SYM_OP s_o() { return ((SYM_OP) ob_self.ob_bruch->b_oben); };
	INT s_oi() { return (((INTEGER_OP) ob_self.ob_bruch->b_oben)->s_i()); };
	SYM_OP s_u() { return ((SYM_OP) ob_self.ob_bruch->b_unten); };
	INT s_ui() { return (((INTEGER_OP) ob_self.ob_bruch->b_unten)->s_i()); };
	void c_i(INT info) { ob_self.ob_bruch->b_info = info; };
	void c_o(SYM_OP oben) { ob_self.ob_bruch->b_oben = (OP) oben; };
	void c_u(SYM_OP unten) { ob_self.ob_bruch->b_unten = (OP) unten; };

INT m_ioiu(INT oben, INT unten);
INT m_ou(SYM_OP oben, SYM_OP unten);
INT b_ou(SYM_OP oben, SYM_OP unten);
INT m_scalar(SYM_OP a);
INT freeself();
INT copy(BRUCH_OP nach);
INT scalar_it();
INT kuerzen();
INT add(SYM_OP b, SYM_OP c);
INT add_scalar(SYM_OP b, SYM_OP c);
INT add_bruch(BRUCH_OP b, SYM_OP c);
INT addinvers(BRUCH_OP b);
INT addinvers_apply();
INT add_apply(SYM_OP b);
INT add_apply_bruch(BRUCH_OP b);
INT invers(BRUCH_OP b);
INT mult(SYM_OP b, SYM_OP c);
INT mult_integer(SYM_OP b, BRUCH_OP c);
INT mult_bruch(BRUCH_OP b, BRUCH_OP c);
INT mult_apply(SYM_OP b);
INT compare(SYM_OP b);
INT einsp();
INT negeinsp();
INT nullp();
INT sprint(BYTE *str);
INT latex(FILE *fp);
};

INT bruch_anfang();
INT bruch_ende();

#endif /* BRUCH_INCLUDED */

