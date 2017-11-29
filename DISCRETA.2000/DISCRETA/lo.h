/* lo.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef LO_INCLUDED
#define LO_INCLUDED

class longint_ob : public SYM_OB {
public:
INT init();
INT fprint(FILE *f);
INT copy(LONGINT_OP b);
INT freeself();
INT add(SYM_OP b, SYM_OP c);
INT add_integer(INTEGER_OP b, LONGINT_OP c);
INT mult(SYM_OP b, SYM_OP c);
INT mult_longint(LONGINT_OP b, LONGINT_OP c);
INT mult_integer(INTEGER_OP b, LONGINT_OP c);
INT addinvers_apply();
INT addinvers(LONGINT_OP b);
INT add_apply(SYM_OP b);
INT add_apply_integer(INTEGER_OP b);
INT add_apply_longint(LONGINT_OP b);
INT mult_apply(SYM_OP b);
INT mult_apply_longint(LONGINT_OP b);
INT mult_apply_integer(INTEGER_OP b);
#ifdef MATRIXTRUE
INT mult_apply_matrix(MATRIX_OP b);
#endif
INT dec();
INT inc();
INT posp();
INT negp();
INT odd();
INT even();
INT nullp();
INT einsp();
INT m_i(INT i);
INT lo2string(STRING_OP p);
INT string2lo(STRING_OP p);
INT sprint(BYTE *str);
INT sprint_latex(BYTE *str);
INT sscan(BYTE *str);
INT compare(SYM_OP c);
INT mod(SYM_OP e, SYM_OP c);
INT ganzdiv(SYM_OP e, SYM_OP c);
INT quores(SYM_OP e, SYM_OP c, SYM_OP d);
INT invers(LONGINT_OP b);

/* in base/iof.C: */
INT calc_size_on_file();
INT write_mem(MEM_OP mem, INT debug_depth);
INT read_mem(MEM_OP mem, INT debug_depth);
};

INT start_longint(void);
INT t_longint_int(LONGINT_OP a);
INT t_int_longint(INTEGER_OP a, LONGINT_OP c);

#endif /* LO_INCLUDED */

