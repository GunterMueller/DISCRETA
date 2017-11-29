/* in.C: */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>
#ifdef MATRIXTRUE
#include <DISCRETA/ma.h>
#endif
#ifdef LONGINTTRUE
#include <DISCRETA/lo.h>
#endif
#ifdef BRUCHTRUE
#include <DISCRETA/bruch.h>
#endif
#ifdef POLYTRUE
#include <DISCRETA/poly.h> /* for MONOM_OB */
#endif
#ifdef CP_TRUE
#include <DISCRETA/cp.h>
#endif

void INTEGER_OB::freeself()
{
	c_obj_k(EMPTY);
}

INT INTEGER_OB::sprint(BYTE *s)
{
	sprintf(s + strlen(s), "%ld", s_i());
	return(OK);
}

INT INTEGER_OB::sprint_latex(BYTE *s)
{
	sprintf(s + strlen(s), "%ld", s_i());
	return(OK);
}

INT INTEGER_OB::sscan(BYTE *s)
{
	INT eingabe;

	sscanf(s, "%ld", &eingabe);
	m_i(eingabe);
	return OK;
}

INT INTEGER_OB::invers_apply()
{ 
#ifdef BRUCHTRUE
	{
	BRUCH_OB c;
	
	c.m_ioiu(1, s_i());
	c.swap(this);
	return OK;
	}
#endif /* BRUCHTRUE */
	return(0L);
}

INT INTEGER_OB::addinvers_apply()
{ 
	if (s_obj_k() != INTEGER) {
		printf("warning: INT::addinvers_apply() no INT obj");
		return ((SYM_OP)this)->addinvers_apply();
		}
	c_i(- s_i());
	return OK;  
}

INT INTEGER_OB::addinvers(INTEGER_OP b)
{
	if (s_obj_k() != INTEGER) {
		printf("warning: INT::addinvers() no INT obj");
		return ((SYM_OP)this)->addinvers(b);
		}
	b->m_i(- s_i());
	return OK;  
}

INT INTEGER_OB::inc()
{
	INT i;

	if (s_obj_k() != INTEGER) {
		printf("warning: INT::inc() no INT obj");
		return ((SYM_OP)this)->inc();
		}
	i = s_i() + 1;
	c_i(i); 
	return(OK); 
}

INT INTEGER_OB::dec()
{
	INT i;

	if (s_obj_k() != INTEGER) {
		printf("warning: INT::dec() no INT obj");
		return ((SYM_OP)this)->dec();
		}
	i = s_i() - 1;
	c_i(i); 
	return(OK); 
}

INT INTEGER_OB::mult_integer(INTEGER_OP b, INTEGER_OP d)
/* Vorm. mult_integer_integer(a, b, d) */
{
	INT l, i;

	if (s_obj_k() != INTEGER) {
		printf("warning: INT::mult_integer() no INT obj");
		return ((SYM_OP)this)->mult(b, d);
		}
	l = log_10() + b->log_10();
	if (l > 7) {
#ifdef LONGINTTRUE
		LONGINT_OB c, e;

		t_int_longint(b, &e);
		t_int_longint(this, &c);
		c.mult_longint(&e, (LONGINT_OP) d); 
		return(OK);
#else /* LONGINTTRUE */
		printf("INT::mult_integer: %ld * %ld\n", 
			s_i(), b->s_i());
		return error("INT::mult_integer:"
			"no LONGINT");
#endif /* LONGINTTRUE */
		}
	i = s_i() * b->s_i();
	d->m_i(i);
	return(OK);
}

INT INTEGER_OB::mult(SYM_OP b, SYM_OP d)
/* before: mult_integer(a, b, d) */
{
	INT erg = OK;
	
	if (s_obj_k() != INTEGER) {
		printf("warning: INT::mult() no INT obj\n");
		println();
		printobjectkind();
		return ((SYM_OP)this)->mult(b, d);
		}

	switch (b->s_obj_k()) {
#ifdef BRUCHTRUE
		case BRUCH: 
			erg += ((BRUCH_OP)b)->mult(this, d); break;
#endif
		case INTEGER:
			erg += mult_integer((INTEGER_OP)b, (INTEGER_OP) d); break;
#ifdef LONGINTTRUE
		case LONGINT:
			erg += ((LONGINT_OP)b)->mult(this, d); break;
#endif
#ifdef MATRIXTRUE
		case KRANZTYPUS :
		case MATRIX:
			erg += ((MATRIX_OP)b)->mult_scalar(this, (MATRIX_OP)d); break;
#endif
#ifdef MONOMTRUE
		case MONOM:
			erg += ((MONOM_OP)b)->mult_scalar(this, (MONOM_OP) d); break;
#endif
#ifdef POLYTRUE
		/* case POW_SYM:
		case ELM_SYM:
		case HOM_SYM:
		case MONOMIAL:
		case SCHUR:
		case SCHUBERT:
		case GRAL:
		case POLYNOM:
			erg += mult_scalar_polynom(a,b,d);break; !!! */
#endif
#ifdef UNDEF
#ifdef SCHUBERTTRUE
		case SCHUBERT:
			erg += mult_scalar_schubert(a,b,d);break;
#endif
#ifdef SCHURTRUE
		case SCHUR:
			erg += mult_scalar_schur(a,b,d);break;
#endif
#endif
#ifdef CHARTRUE
		case SYMCHAR:
			erg += mult_scalar_symchar(a,b,d);break;
#endif
#ifdef VECTORTRUE
		/* case VECTOR:
			erg += mult_scalar_vector(a,b,d);break; !!! */
#endif
#ifdef PERMTRUE
		case PERMUTATION:
			if (nullp()) 
				((INTEGER_OP)d)->m_i(0); 
			else {
				b->printobjectkind();
				erg += error("INT::mult: wrong second kind");
				}
			break;
#endif
		default:
			{
			b->printobjectkind();
			error("INT::mult:wrong second kind");
			return(ERROR);
			}
		}
	if (erg != OK)
		return error("INT::mult: error in computing");
	return erg;
}

INT INTEGER_OB::even()
{
	if (s_obj_k() != INTEGER) {
		printf("warning: INT::even() no INT obj");
		}
	return (s_i() % 2 == 0); 
}

INT INTEGER_OB::posp()
{
	if (s_obj_k() != INTEGER) {
		printf("warning: INT::posp() no INT obj");
		}
	return (s_i() > 0); 
}

INT INTEGER_OB::negp()
{
	if (s_obj_k() != INTEGER) {
		printf("warning: INT::negp() no INT obj");
		}
	return (s_i() < 0); 
}

INT INTEGER_OB::add_integer(INTEGER_OP b, INTEGER_OP c)
/* before: add_integer_integer(OP a, OP b, OP c) */
{
	INT i;

	if (s_obj_k() != INTEGER) {
		printf("warning: INT::add_integer() no INT obj");
		return ((SYM_OP)this)->add(b, c);
		}
	if ((b->s_i() >1000000L) || 
		(b->s_i() < -1000000L)) {
#ifdef LONGINTTRUE
		LONGINT_OB d;
		
		d.m_i(b->s_i()); 
		add(&d, c); 
		return(OK);
#else
		return error("INT::add_integer:"
		"overflow no LONGINT");
#endif
		};
	i = s_i() + b->s_i();
	c->m_i(i);
	return OK;
}

INT INTEGER_OB::add(SYM_OP b, SYM_OP c)
/* before: add_integer(OP a, OP b, OP c) */
/* das erste object ist vom typ INTEGER, 
 * das ergebnis ist ein leeres object */
{
	INT erg = OK;

	if (s_obj_k() != INTEGER) {
		printf("warning: INT::add() no INT obj");
		return ((SYM_OP)this)->add(b, c);
		}
	switch (b->s_obj_k()) {
#ifdef BRUCHTRUE
		case BRUCH:
			erg += ((BRUCH_OP)b)->add(this, c); break;
#endif
		case INTEGER:
			erg += add_integer((INTEGER_OP)b, (INTEGER_OP)c); break;
#ifdef LONGINTTRUE
		case LONGINT:
			erg += ((LONGINT_OP)b)->add(this, c); break;
#endif
#ifdef POLYTRUE
		/* case POLYNOM:
			erg += add_scalar_polynom(a,b,c); break; !!! */
#endif
		default :
			{
			if (nullp()) 
				return b->copy(c);
			printobjectkind(); 
			b->printobjectkind();
			return error("INT::add: wrong second type");
			};
		}
	if (erg != OK)
		return error("INT::add: (2) error in computation");
	return erg;
}

INT INTEGER_OB::comp_integer(INTEGER_OP b)
{
	INT ai, bi;

	if (s_obj_k() != INTEGER) {
		printf("warning: INT::comp_integer() no INT obj");
		return ((SYM_OP)this)->sym_comp(b);
		}
	ai = s_i();
	if (b->s_obj_k() == BRUCH) {
		return - b->sym_comp(this);
		}
	bi = b->s_i();
	if (ai == bi) return(0);
	if (ai > bi) return(1);
	return(-1);
}

INT INTEGER_OB::compare(SYM_OP b)
/* a ist vom typ INTEGER, b hat unbekannten typ */
{
	if (s_obj_k() != INTEGER) {
		printf("warning: INT::compare() no INT obj");
		return ((SYM_OP)this)->sym_comp(b);
		}
	
	switch (b->s_obj_k()) {
#ifdef BRUCHTRUE
		case BRUCH:
			return (-1) * ((BRUCH_OP)b)->compare(this);
#endif
		case INTEGER:
			return comp_integer((INTEGER_OP)b);
#ifdef LONGINTTRUE
		case LONGINT:
			return (-1) * ((LONGINT_OP)b)->compare(this);
#endif
		default:
			{
			b->printobjectkind();
			error("INT::comp: wrong second type");
			return(ERROR);
			}
		}
} 

INT INTEGER_OB::nullp()
{
	if (s_obj_k() != INTEGER) {
		printf("warning: INT::nullp() no INT obj");
		return ((SYM_OP)this)->nullp();
		}
	return (s_i() == 0);
}

INT INTEGER_OB::einsp()
{
	if (s_obj_k() != INTEGER) {
		printf("warning: INT::einsp() no INT obj");
		return ((SYM_OP)this)->einsp();
		}
	return (s_i() == 1); 
}

INT INTEGER_OB::negeinsp()
{
	if (s_obj_k() != INTEGER) {
		printf("warning: INT::negeinsp() no INT obj");
		return ((SYM_OP)this)->negeinsp();
		}
	return (s_i() == -1); 
}

INT INTEGER_OB::copy(INTEGER_OP b)
{
	if (b == this)
		return error("INT::copy()|b == this");	
	b->freeself();	
	b->c_obj_k(INTEGER);
	b->c_i(s_i());
	return OK;
}
	
INT INTEGER_OB::invers(INTEGER_OP b)
{
	if (s_obj_k() != INTEGER) {
		printf("warning: INT::invers() no INT obj\n");
		return ((SYM_OP)this)->invers(b);
		}

	if (einsp())
		return(copy(b));
	if (negeinsp())
		return(copy(b));
	
#ifdef BRUCHTRUE
	{
	BRUCH_OB c;
	c.m_ioiu(1, s_i());
	c.swap(b);
	return OK;
	}
#else
	error("INT::invers: BRUCH not available");
	printf("%ld\n", s_i());
	return(ERROR);
#endif
}

INT INTEGER_OB::add_apply_integer(INTEGER_OP b)
{ 
	INT i;

	if (s_obj_k() != INTEGER) {
		printf("warning: INT::add_apply_integer() no INT obj");
		return ((SYM_OP)this)->add_apply(b);
		}
	if ( 
		(s_i() >1000000L) ||  (b->s_i() > 1000000L) ||
		(s_i() < -1000000L) ||  (b->s_i() < -1000000L) ) {
#ifdef LONGINTTRUE
		LONGINT_OB c;
		t_int_longint(b, &c);
		return c.add(this, b);
#else /* LONGINTTRUE */
		return error("INT::add_apply_integer:"
			"Overflow no LONGINT");
#endif /* LONGINTTRUE */
		}
	
	i = s_i() + b->s_i();
	b->c_i(i); 
	return OK;
}

INT INTEGER_OB::add_apply(SYM_OP b)
{
	INT erg = OK;
	
	if (s_obj_k() != INTEGER) {
		printf("warning: INT::add_apply() no INT obj");
		return ((SYM_OP)this)->add_apply(b);
		}
	switch (b->s_obj_k()) {
		case INTEGER: 
			erg += add_apply_integer((INTEGER_OP)b); break;
#ifdef LONGINTTRUE 
		case LONGINT: 
			{
			LONGINT_OB c;
			b->swap(&c);
			b->c_obj_k(EMPTY);
			erg += c.add(this, b);
			break;
			}
#endif /* LONGINTTRUE */
		case POLYNOM:
		case SCHUBERT:
		case SCHUR: 
			/* erg += add_apply_scalar_polynom(a, b); !!! */ break;
		default: 
			/*
			b->printobjectkind();
			error("INT::add_apply: wrong second type");
			return(ERROR);
			*/
			{
			SYM_OB c;

			c = *b;
			b->c_obj_k(EMPTY);
			erg += add(&c, b);
			break;
			}
		}
	if (erg != OK)
		error("INT::add_apply: error during computation");
	return erg;
}	

INT INTEGER_OB::mult_apply_integer(INTEGER_OP b)
{ 
	INT l, i;

	if (s_obj_k() != INTEGER) {
		printf("warning: INT::mult_apply_integer() no INT obj");
		return ((SYM_OP)this)->mult_apply(b);
		}
	if ( 
		(s_i() < 10000L) && (s_i() > -10000L) &&
		(b->s_i() < 10000L) && (b->s_i() > -10000L) )  {
		b->m_i(s_i() * b->s_i());
		return( OK ); 
		}
	else
		l = log_10() + b->log_10();
	if ( l > 8L ) {
#ifdef LONGINTTRUE
		LONGINT_OB c;
		t_int_longint(b, &c);
		return c.mult(this, b);
#else /* LONGINTTRUE */
		return 
		error("INT::mult_apply_integer: "
			"LONGINT not available");
#endif /* LONGINTTRUE */
		}
	i = s_i() * b->s_i();
	b->m_i(i);
	return(OK);
}


INT INTEGER_OB::mult_apply(SYM_OP b)
/* b = b * a */
{
	if (s_obj_k() != INTEGER) {
		printf("warning: INT::mult_apply() no INT obj");
		return ((SYM_OP)this)->mult_apply(b);
		}
	switch (b->s_obj_k()) {
		case INTEGER:
			return (mult_apply_integer((INTEGER_OP)b));
#ifdef LONGINTTRUE
		case LONGINT: 
			{
			LONGINT_OB c;
			b->swap(&c);
			return c.mult(this, b);
			}
#endif /* LONGINTTRUE */
		case KRANZTYPUS :
#ifdef MATRIXTRUE
		/* case MATRIX:
			return(mult_apply_integer_matrix(a,b)); !!! */
#endif /* MATRIXTRUE */
#ifdef MONOMTRUE
		/* case MONOM:
			return(mult_apply_scalar_monom(a,b)); !!! */
#endif
		case SCHUR:
		case SCHUBERT:
		case GRAL:
#ifdef POLYNOMTRUE
		case POLYNOM:
			return mult_apply_scalar_polynom(a,b);
#endif
		default: {
			SYM_OP c = (SYM_OP) callocobject("mult_apply_integer");
			INT erg;
			*c = *b; 
			b->c_obj_k(EMPTY);
			erg = mult(c, b);
			if (erg == ERROR) {
				c->printobjectkind();
				error("mult_apply_integer: wrong second type");
				}
			freeall(c); 
			return erg;
			}
		}
}

INT INTEGER_OB::log_10()
/* anzahl stellen */ /* vorm. intlog(OP a) */
{
	INT ai = s_i();
	
	if (ai < 0) ai = -ai;
	if (ai >= 1000000000L) return(10L);
	if (ai >= 100000000L) return(9L);
	if (ai >= 10000000L) return(8L);
	if (ai >= 1000000L) return(7L);
	if (ai >= 100000L) return(6L);
	if (ai >= 10000L) return(5L);
	if (ai >= 1000L) return(4L);
	if (ai >= 100L) return(3L);
	if (ai >= 10L) return(2L);
	return(1);
}

INT INTEGER_OB::zero()
{
	c_i(0);
	return(OK);
}

INT INTEGER_OB::one()
{
	c_i(1);
	return(OK);
}

INT INTEGER_OB::m_one()
{
	c_i(-1);
	return(OK);
}

INT INTEGER_OB::homo_z(INT z)
{
	c_i(z);
	return(OK);
}

INT INTEGER_OB::ganzdiv(INTEGER_OP b, INTEGER_OP c)
{
	INT erg = OK;
	
	switch (b->s_obj_k()) {
	case INTEGER:
		c->m_i(s_i() / b->s_i());
		if ((s_i() < 0) && (b->s_i() > 0)) c->dec();
		if ((s_i() > 0) && (b->s_i() < 0)) c->dec();
		break;
#ifdef LONGINTTRUE
		case LONGINT:
			{
			LONGINT_OB d;
			d.m_i(s_i());
			erg += d.ganzdiv(b, c);
			return erg;
			}
#endif /* LONGINTTRUE */
#ifdef BRUCHTRUE       /* AK 130691 V1.2 */
		case BRUCH: 
#if 0
			if (einsp(S_B_U(b)))
				erg += ganzdiv_integer(
					a,S_B_O(b),c);
			else {
				b->printobjectkind();
				return error("ganzdiv_integer: "
				"wrong bruch as second type");
				}
#endif
				b->printobjectkind();
				return error("INTEGER::ganzdiv: "
				"wrong bruch as second type");
#endif /* BRUCHTRUE */
		default:
			printobjectkind();
			b->printobjectkind();
			return error("INTEGER::ganzdiv: "
			"wrong second type");
		}
	return erg;
}

INT INTEGER_OB::quores(INTEGER_OP b, 
	INTEGER_OP c, INTEGER_OP d)
{
	switch (b->s_obj_k()) {
		case INTEGER: 
			c->m_i(s_i() / b->s_i());
			d->m_i(s_i() % b->s_i());
			return OK;
#ifdef LONGINTTRUE
		case LONGINT:
			{
			LONGINT_OB e;
			e.m_i(s_i()); 
			e.quores(b, c, d);
			return(OK);
			};
#endif /* LONGINTTRUE */
		default:
			{
			b->printobjectkind();
			return error("INTEGER::quores: "
			"wrong second type");
			}
		}
}

INT INTEGER_OB::modulo(INTEGER_OP b, INTEGER_OP c)
/* before: mod_integer(OP a, OP b, OP c) */
{
	switch (b->s_obj_k()) {
		case INTEGER: 
				c->m_i(s_i() % b->s_i());
				return OK;
#ifdef LONGINTTRUE
		case LONGINT: 
			{
			LONGINT_OB d; 
			
			d.m_i(s_i());
			d.mod(b, c);
			return(OK);
			};
#endif /* LONGINTTRUE */
		default:
			{
			b->printobjectkind();
			error("INTEGER::modulo: "
			"wrong second type");
			return(ERROR);
			}
		}
}


