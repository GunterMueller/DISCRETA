/* nu.C: */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef PERMTRUE
#include <DISCRETA/perm.h>
#endif
#ifdef MATRIXTRUE
#include <DISCRETA/ma.h>
#endif
#ifdef CP_TRUE
#include <DISCRETA/cp.h>
#endif
#ifdef LONGINTTRUE
#include <DISCRETA/lo.h>
#endif
#ifdef BRUCHTRUE
#include <DISCRETA/bruch.h>
#endif
#include <DISCRETA/in.h>
#ifdef SOLVABLE_TRUE
#include <DISCRETA/solvable.h>
#endif
#ifdef POLYTRUE
#include <DISCRETA/poly.h>
#endif
#ifdef UNIPOLYTRUE
#include <DISCRETA/unip.h>
#endif

#undef FG_DEBUG

#if TEXDOCU
INT n_choose_k_ob(SYM_OP n, SYM_OP k, SYM_OP a)
#else
computes ${n \choose k}$ with objects. 
k must be a (small) integer ! (integer object). 
#endif
{
	SYM_OB num, denom, tmp;
	INT i, kk;
	
	if (k->s_obj_k() != INTEGER)
		return error("n_choose_k_ob() k must be an integer object");
	kk = k->s_i_i();
	a->m_i_i(1);
	n->copy(&num);
	for (i = 0; i < kk; i++) {
		if (i > 0)
			num.dec();
		// num.m_i_i(n - i);
		
		denom.m_i_i(i + 1);
		a->mult(&num, &tmp);
		tmp.ganzdiv(&denom, a);
		}
	return OK;
}

INT SYM_OB::operator+=(SYM_OB &a)
{
	SYM_OB c;
	
	add(&a, &c);
	swap(&c);
	return OK;
}

INT SYM_OB::operator*=(SYM_OB &a)
{
	SYM_OB c;
	
	mult(&a, &c);
	swap(&c);
	return OK;
}

INT SYM_OB::operator/=(SYM_OB &a)
{
	SYM_OB c;
	
	ganzdiv(&a, &c);
	swap(&c);
	return OK;
}

#if TEXDOCU
INT SYM_OB::m_i_i(INT i)
#endif
{
	INTEGER_OP p = (INTEGER_OP) this;

	freeself();
	p->m_i(i);
	return OK;
}

#if TEXDOCU
INT SYM_OB::s_i_i()
#endif
{
	INTEGER_OP p = (INTEGER_OP) this;
	
	if (s_obj_k() != INTEGER)
		return error("SYM_OB::s_i_i() not an INTEGER object !");
	return p->s_i();
}

#if TEXDOCU
INT SYM_OB::fakul_int(int i)
#endif
{
	INTEGER_OB ii;
	SYM_OB b;
	INT ret;

	ii.m_i(i);
	ret = ii.fakul(&b);
	swap(&b);
	return ret;
}

#if TEXDOCU
INT SYM_OB::order(INT *order)
#endif
{
	SYM_OB q;

	if (s_obj_k() == PERMUTATION) {
#ifdef PERMTRUE
		return(((PERMUTATION_OP)this)->order(order));
#endif
		}
	copy(&q);
	*order = 1;
	while (!q.einsp()) {
		mult(&q, &q);
		(*order)++;
		}
	return(OK);
}

#if TEXDOCU
INT SYM_OB::order_if_prime(INT *ord)
#endif
{
	INT ord1;

	if (s_obj_k() == PERMUTATION) {
#ifdef PERMTRUE
		return(((PERMUTATION_OP)this)->order_if_prime(ord));
#endif
		}
	order(&ord1);
	ord_to_order_if_prime(ord1, ord);
	return(OK);
}

#if TEXDOCU
INT SYM_OB::order_if_prime_power(INT *ord, INT *prime, INT *k)
#endif
{
	INT ord1;
	
	if (s_obj_k() == PERMUTATION) {
#ifdef PERMTRUE
		return(((PERMUTATION_OP)this)->order_if_prime_power(ord, prime, k));
#endif
		}
	order(&ord1);
	ord_to_order_if_prime_power(ord1, ord, prime, k);
	return(OK);
}

#if TEXDOCU
INT SYM_OB::add(SYM_OP b, SYM_OP d)
#else
d := this + b.
#endif
{
	INT erg = OK;
	
	if ((this == d) && (b == d)) {
		SYM_OP c = callocobject("SYM::add()");
		*c = *this; 
		d->c_obj_k(EMPTY);
		erg += c->add(c, d); 
		erg += freeall(c); 
		goto add_ende;
		}
	else if	(this == d) {
		SYM_OP c = callocobject("SYM::add()");
		*c = *this; 
		d->c_obj_k(EMPTY);
		erg += c->add(b, d); 
		erg += freeall(c); 
		goto add_ende;
		}
	else if	(b == d) {
		SYM_OP c = callocobject("SYM::add()");
		*c = *b; 
		d->c_obj_k(EMPTY);
		erg += add(c, d); 
		erg += freeall(c); 
		goto add_ende;
		}
	
	
	else if	(emptyp()) {
		erg += b->copy(d);
		goto add_ende;
		}
	else if	(b->emptyp()) {
		erg +=  copy(d);
		goto add_ende;
		}
	if (! d->emptyp()) 
		if (d->s_obj_k() != INTEGER) 
			erg += d->freeself();

	switch (s_obj_k()) {
	case INTEGER :
		erg += ((INTEGER_OP)this)->add(b, d); break;
#ifdef PARTTRUE
	/* case PARTITION:
		erg += add_partition(a,b,d); break;!!! */
#endif
#ifdef UNIPOLYTRUE
	case UNIPOLY :
		erg += ((UNIPOLY_OP)this)->add((UNIPOLY_OP)b, (UNIPOLY_OP)d); break;
#endif
#ifdef POLYTRUE
	case POLYNOM :
		erg += ((POLYNOM_OP)this)->add(b, (POLYNOM_OP)d); break;
#endif
#ifdef MONOMTRUE
	case MONOM :
		erg += ((MONOM_OP)this)->add(b, (MONOM_OP)d); break;
#endif
#ifdef VECTORTRUE
	case VECTOR :
		erg += ((VECTOR_OP)this)->add((VECTOR_OP)b, (VECTOR_OP)d); break;
#endif
#ifdef LONGINTTRUE
	case LONGINT :
		erg += ((LONGINT_OP)this)->add(b, d); break;
#endif
#ifdef MATRIXTRUE
	case KRANZTYPUS:
	case MATRIX :
		erg += ((MATRIX_OP)this)->add((MATRIX_OP)b, (MATRIX_OP)d); break;
#endif
#ifdef BRUCHTRUE
	case BRUCH :
		erg += ((BRUCH_OP)this)->add(b, d); break;
#endif
	default: 
		{
		if (nullp()) { erg += b->copy(d); break; }
		if (b->nullp()) { erg += copy(d); break; }
			printobjectkind(); 
			b->printobjectkind();
			return error("SYM::add: wrong types");
		}
	};
add_ende:
	if (erg != OK) {
		printobjectkind(); 
		b->printobjectkind();
		return error("SYM::add: error during computation");
		}
	return erg;
}

#if TEXDOCU
INT SYM_OB::add_apply(SYM_OP b)
#else
b := this + b
#endif
{
	INT erg = OK;
	
	if (emptyp()) return OK;
	if (b->emptyp()) {
		erg += copy(b);
		goto add_apply_ende;
		}
	if (nullp()) return OK;
	if (b->nullp()) {
		erg += copy(b);
		goto add_apply_ende;
		}
	switch (s_obj_k()) {
#ifdef BRUCHTRUE
		case BRUCH:
			erg += ((BRUCH_OP)this)->add_apply(b); break;
#endif
		case INTEGER:
			erg += ((INTEGER_OP)this)->add_apply(b); break;
#ifdef LONGINTTRUE
		case LONGINT:
			erg += ((LONGINT_OP)this)->add_apply(b); break;
#endif
#ifdef MATRIXTRUE
		case MATRIX:
			erg += ((MATRIX_OP)this)->add_apply(b); break;
#endif
#ifdef UNIPOLYTRUE
	case UNIPOLY :
			erg += ((UNIPOLY_OP)this)->add_apply((UNIPOLY_OP)b); break;
#endif
#ifdef POLYTRUE
		case POLYNOM:
			erg += ((POLYNOM_OP)this)->add_apply((POLYNOM_OP) b); break;
#endif
		case VECTOR:
			erg += ((VECTOR_OP)this)->add_apply((VECTOR_OP)b); break;
		default: {
			printobjectkind();
			return error("SYM::add_apply: wrong first type");
			}
		}
add_apply_ende:
	if (erg != OK) {
		printobjectkind();
		b->printobjectkind();
		return error("SYM::add_apply: error during computation");
		}
	return erg;
}

#if TEXDOCU
INT SYM_OB::addinvers(SYM_OP res)
#else
res := -this.
#endif
{
	INT erg = OK;

	if (this == res) {
		SYM_OP c = callocobject("SYM::addinvers"); 
		*c = *this; 
		res->c_obj_k(EMPTY);
		erg += c->addinvers(res); 
		erg += freeall(c); 
		if (erg != OK) 
			return error("SYM::addinvers:(1) error in computing");
		return erg;
	};

	if (! res->emptyp()) 
		erg += res->freeself();
	if (emptyp()) return(OK);

	switch (s_obj_k()) {
#ifdef BRUCHTRUE
	case BRUCH :
		erg += ((BRUCH_OP)this)->addinvers((BRUCH_OP) res); break;
#endif
#ifdef INTEGERTRUE
	case INTEGER :
		erg += ((INTEGER_OP)this)->addinvers((INTEGER_OP) res); break;
#endif	
#ifdef LONGINTTRUE
	case LONGINT :
		erg+= ((LONGINT_OP)this)->addinvers((LONGINT_OP) res);break;
#endif
#ifdef MONOMTRUE
	case MONOM :
		erg += ((MONOM_OP)this)->addinvers((MONOM_OP) res); break;
#endif
#ifdef UNIPOLYTRUE
	case UNIPOLY :
		erg += ((UNIPOLY_OP)this)->addinvers((UNIPOLY_OP) res); break;
#endif
#ifdef POLYTRUE
	case POLYNOM :
		erg += ((POLYNOM_OP)this)->addinvers((POLYNOM_OP) res); break;
#endif
#ifdef VECTORTRUE
	case VECTOR :
		erg += ((VECTOR_OP)this)->addinvers((VECTOR_OP)res); break;
#endif
	default: 
		{
			printobjectkind();
			return error("SYM::addinvers: wrong type");
		}
	};
	if (erg != OK) {
		return error("SYM::addinvers: error in computing");
		}
	return erg;
}

#if TEXDOCU
INT SYM_OB::addinvers_apply()
#else
this := -this
#endif
{
	if (emptyp()) 
		return(OK);
	
	switch (s_obj_k()) {
#ifdef BRUCHTRUE
	case BRUCH :
		return(((BRUCH_OP)this)->addinvers_apply());
#endif
#ifdef INTEGERTRUE
	case INTEGER :
		return(((INTEGER_OP)this)->addinvers_apply());
#endif
#ifdef LONGINTTRUE
	case LONGINT :
		return(((LONGINT_OP)this)->addinvers_apply());
#endif
#ifdef MONOMTRUE
	case MONOM :
		return ((MONOM_OP)this)->addinvers_apply();
#endif
#ifdef UNIPOLYTRUE
	case UNIPOLY :
		return ((UNIPOLY_OP)this)->addinvers_apply(); break;
#endif
#ifdef POLYTRUE
	case POLYNOM :
		return ((POLYNOM_OP)this)->addinvers_apply();
#endif
#ifdef VECTORTRUE
	case VECTOR :
		return(((VECTOR_OP)this)->addinvers_apply());
#endif
	default: 
		{
			printobjectkind();
			return error("SYM::addinvers_apply: wrong type");
		}
	};
}

#if TEXDOCU
INT SYM_OB::mult(SYM_OP b, SYM_OP d)
#else
d := this * b.
#endif
{
	SYM_OP c;
	INT erg = OK;
	
	if (b == NULL || d == NULL)
		return error("SYM::mult: NULL object"); 

	if (emptyp())
		return error("SYM::mult: b empty object"); 
	if (b->emptyp()) 
		return error("SYM::mult: b empty object"); 

	/* sonderfaelle bei gleichen variablennamen */
	if	((this == d) && (b == d)) { 
		c = callocobject("SYM::mult");
		*c = *this;
		d->c_obj_k(EMPTY);
		c->mult(c, d);
		freeall(c);
		return(OK); 
		}
	else if	(this == d) { 
		c = callocobject("SYM::mult"); 
		*c = *this; 
		d->c_obj_k(EMPTY);
		c->mult(b, d); 
		freeall(c); 
		return(OK); 
		}
	else if (b == d) { 
		c = callocobject("SYM::mult"); 
		*c = *b ;
		d->c_obj_k(EMPTY);
		mult(c, d); 
		freeall(c); 
		return(OK); 
		}

	/*freigabe des speichers belegt durch d */
	if (! d->emptyp())
		d->freeself();


	switch (s_obj_k()) {
#ifdef BRUCHTRUE
	case BRUCH :
		erg += ((BRUCH_OP)this)->mult(b, d); break;
#endif
#ifdef INTEGERTRUE
	case INTEGER :
		erg += ((INTEGER_OP)this)->mult(b, d); break;
#endif
#ifdef UNIPOLYTRUE
	case UNIPOLY :
		erg += ((UNIPOLY_OP)this)->mult((UNIPOLY_OP) b, (UNIPOLY_OP) d); break;
#endif
#ifdef POLYTRUE
	case POLYNOM :
		erg += ((POLYNOM_OP)this)->mult(b, (POLYNOM_OP) d); break;
#endif
#ifdef MATRIXTRUE
	case MATRIX : 	
		erg += ((MATRIX_OP)this)->mult(b, (MATRIX_OP)d); break;
#endif
#ifdef LONGINTTRUE
	case LONGINT:
		erg += ((LONGINT_OP)this)->mult(b, d); break;
#endif
#ifdef PERMTRUE
	case PERMUTATION: 
		erg += ((PERMUTATION_OP)this)->mult((PERMUTATION_OP)b, 
			(PERMUTATION_OP)d); break;
#endif
#ifdef VECTORTRUE
	case VECTOR :
		switch (b->s_obj_k()) {
#ifdef BRUCHTRUE
		case BRUCH:
#endif
		case LONGINT:
		case INTEGER:
			erg += ((VECTOR_OP)this)->mult_scalar(b, (VECTOR_OP)d); break;
		case VECTOR:
			erg += ((VECTOR_OP)this)->mult_vector(
				(VECTOR_OP)b, (VECTOR_OP)d); break;
#ifdef MATRIXTRUE
		/* case MATRIX:
			erg+=mult_vector_matrix(a,b,d);break; !!! */
#endif
		default: 
				b->printobjectkind();
				error("SYM::mult_vector:wrong second type");
				return ERROR;
		};
		break;
#endif
		
	default: 
		{
			printobjectkind(); 
			b->printobjectkind();
			return error("SYM::mult:wrong types");
		}
	}
	if (erg != OK)
		return error("SYM::mult: error in computation");
	return erg;
}

#if TEXDOCU
INT SYM_OB::mult_apply(SYM_OP b)
#else
b := this * b.
#endif
{
	INT erg = OK;
	
	if (this == b) {
		SYM_OP c;
		c = callocobject("SYM::mult_apply");
		copy(c); 
		c->mult_apply(b); 
		freeall(c); 
		return(OK);
		}
	if (b->emptyp()) 
		return OK;
	if (emptyp()) 
		return(b->freeself());

	switch (s_obj_k()) {
#ifdef BRUCHTRUE
		case BRUCH:
			erg += ((BRUCH_OP)this)->mult_apply(b); break;
#endif
		case INTEGER:
			erg += ((INTEGER_OP)this)->mult_apply(b); break;
#ifdef LONGINTTRUE
		case LONGINT:
			erg += ((LONGINT_OP)this)->mult_apply(b); break;
#endif
#ifdef MATRIXTRUE
		case MATRIX:
			erg += ((MATRIX_OP)this)->mult_apply(b); break;
#endif
#ifdef UNIPOLYTRUE
		case UNIPOLY :
			erg += ((UNIPOLY_OP)this)->mult_apply((UNIPOLY_OP) b); break;
#endif
#ifdef POLYTRUE
		case POLYNOM:
			erg += ((POLYNOM_OP)this)->mult_apply(b); break;
#endif
		case VECTOR:
			erg += ((VECTOR_OP)this)->mult_apply(b); break;
		default: {
			printobjectkind();
			error("SYM::mult_apply: wrong first type");
			return ERROR;
			}
		}
	if (erg != OK)
		error("SYM::mult_apply: "
			"error during computation");
	return erg;
}

#if TEXDOCU
INT SYM_OB::invers(SYM_OP b)
#else
b := this$^{-1}$.
#endif
{
	SYM_OP c;
	INT erg = OK;
	
	/* sonderfaelle bei gleichen variablennamen */
	if (this == b) { 
		c = callocobject("SYM::invers");
		*c = *this; 
		b->c_obj_k(EMPTY); 
		c->invers(b); 
		freeall(c); 
		return(OK); 
		}

	if (! b->emptyp()) 
		erg += b->freeself();
	switch (s_obj_k()) {
#ifdef BRUCHTRUE
	case BRUCH :
		erg += ((BRUCH_OP)this)->invers((BRUCH_OP) b); break;
#endif
#ifdef INTEGERTRUE
	case INTEGER :
		erg += ((INTEGER_OP)this)->invers((INTEGER_OP)b); break;
#endif
	case LONGINT :
		erg += ((LONGINT_OP)this)->invers((LONGINT_OP)b); break;
#ifdef MATRIXTRUE
	case MATRIX : 
			erg += ((MATRIX_OP)this)->invers((MATRIX_OP)b); break;
#endif
#ifdef PERMTRUE
	case PERMUTATION : 
			erg += ((PERMUTATION_OP)this)->invers((PERMUTATION_OP)b); break;
#endif
	default:
			printobjectkind();
			return error("SYM::invers: wrong type");
	};
	if (erg != OK) {
		printobjectkind();
		error("SYM::invers: error during computation");
		}
	return erg;
}

#if TEXDOCU
INT SYM_OB::invers_apply()
#endif
{
	INT erg = OK;
	
	if (emptyp()) 
		return(OK);
	switch (s_obj_k()) {
#ifdef INTEGERTRUE
	case INTEGER :
		erg += ((INTEGER_OP)this)->invers_apply(); break;
#endif /* INTEGERTRUE */
	default: {
		SYM_OP c = callocobject("SYM::invers_apply");
		erg += copy(c);
		erg += c->invers(this);
		erg += freeall(c);
		}
	}
	if (erg != OK) {
		error("SYM::invers_apply: error during computation");
		}
	return erg;
}

#if TEXDOCU
INT SYM_OB::inc()
#endif
{
	switch (s_obj_k()) {
#ifdef INTEGERTRUE
	case INTEGER :
		return(((INTEGER_OP)this)->inc());
#endif
#ifdef LONGINTTRUE
	case LONGINT :
		return(((LONGINT_OP)this)->inc());
#endif
#ifdef MATRIXTRUE
	/* case MATRIX :
		return inc_matrix(a); !!! */
#endif
#ifdef PARTTRUE
	/* case PARTITION :
		return(INC_PARTITION(a));!!! */
#endif
#ifdef PERMTRUE
	/* case PERMUTATION :
		return(inc_permutation(a)); !!! */
#endif
#ifdef VECTORTRUE
	case VECTOR :
		return(((VECTOR_OP)this)->inc());
#endif
	default:
		{
			printobjectkind();
			return error("SYM::inc: wrong type");
		}
	};
}

#if TEXDOCU
INT SYM_OB::dec()
#endif
{
	switch (s_obj_k()) {
#ifdef INTEGERTRUE
	case INTEGER :
		return(((INTEGER_OP)this)->dec());
#endif
#ifdef LONGINTTRUE
	case LONGINT :
		return(((LONGINT_OP)this)->dec());
#endif
#ifdef PARTTRUE
	/* case PARTITION :
		return(dec_partition(a));!!! */
#endif
#ifdef PERMTRUE
	/* case PERMUTATION :
		return(dec_permutation(a)); !!! */
#endif
#ifdef VECTORTRUE
	case VECTOR :
		return(((VECTOR_OP)this)->dec());
#endif
	default: { 
		printobjectkind();
		return error("SYM::dec: wrong type"); }
	};
}

#if TEXDOCU
INT SYM_OB::zero()
#endif
{
	switch (s_obj_k()) {
	case INTEGER:
		return ((INTEGER_OP)this)->zero();
#ifdef LONGINTTRUE
	case LONGINT:
		return ((LONGINT_OP)this)->m_i(0);
#endif
#ifdef MATRIXTRUE
	case MATRIX:
		return ((MATRIX_OP)this)->zero();
#endif
#ifdef UNIPOLYTRUE
	case UNIPOLY:
		return ((UNIPOLY_OP)this)->zero();
#endif
	default:
		printobjectkind();
		return error("SYM::zero()|not yet implemented");
	}
	return(OK);
}

#if TEXDOCU
INT SYM_OB::one()
#endif
{
	switch (s_obj_k()) {
	case INTEGER:
		return ((INTEGER_OP)this)->one();
#ifdef LONGINTTRUE
	case LONGINT:
		return ((LONGINT_OP)this)->m_i(1);
#endif
#ifdef MATRIXTRUE
	case MATRIX:
		return ((MATRIX_OP)this)->one();
#endif
#ifdef UNIPOLYTRUE
	case UNIPOLY:
		return ((UNIPOLY_OP)this)->one();
#endif
#ifdef PERMTRUE
	case PERMUTATION:
		return ((PERMUTATION_OP)this)->one();
#endif
	default:
		printobjectkind();
		return error("SYM::one()|not yet implemented");
	}
	return(OK);
}

#if TEXDOCU
INT SYM_OB::m_one()
#endif
{
	switch (s_obj_k()) {
	case INTEGER:
		return ((INTEGER_OP)this)->m_one();
#ifdef LONGINTTRUE
	case LONGINT:
		return ((LONGINT_OP)this)->m_i(-1);
#endif
#ifdef MATRIXTRUE
	case MATRIX:
		return ((MATRIX_OP)this)->m_one();
#endif
#ifdef UNIPOLYTRUE
	case UNIPOLY:
		return ((UNIPOLY_OP)this)->m_one();
#endif
	default:
		printobjectkind();
		return error("SYM::m_one()|not yet implemented");
	}
	return(OK);
}

#if TEXDOCU
INT SYM_OB::homo_z(INT z)
#endif
{
	switch (s_obj_k()) {
	case INTEGER:
		return ((INTEGER_OP)this)->homo_z(z);
#ifdef LONGINTTRUE
	case LONGINT:
		return ((LONGINT_OP)this)->m_i(z);
#endif
#ifdef MATRIXTRUE
	case MATRIX:
		return ((MATRIX_OP)this)->homo_z(z);
#endif
	default:
		printobjectkind();
		return error("SYM::homo_z()|not yet implemented");
	}
	return(OK);
}

#if TEXDOCU
INT SYM_OB::sym_div(SYM_OP b, SYM_OP d)
#else
$d := this \cdot b^{-1}$
#endif
{
	/* AK 031286 als invers*mult */
	INT erg = OK;
	
	SYM_OP c = callocobject("SYM::sym_div");
	erg += b->invers(c); 
	erg += mult(c, d); 
	erg += freeall(c); 
	if (erg != OK) {
		error("sym_div: error during computation");
		}
	return erg;
}

#if TEXDOCU
INT SYM_OB::quores(SYM_OP b, SYM_OP c, SYM_OP d)
#else
#endif
/* c = ganzdiv(a,b)  d = mod(a,b) */
{
	SYM_OP e; 
	INT erg;
	
	if (c == d)
		return error("quores: c == d");
	if (this == c) {
		e = callocobject("SYM::quores");
		*e = *this;
		c_obj_k(EMPTY);
		erg = e->quores(b, c, d); 
		freeall(e);
		return(erg);
		}
	if (this == d) {
		e = callocobject("SYM::quores");
		*e = *this;
		d->c_obj_k(EMPTY);
		erg = e->quores(b,c,d); 
		freeall(e);
		return(erg);
		}
	if (b == c) {
		e =callocobject("SYM::quores");
		*e = *b;
		c->c_obj_k(EMPTY);
		erg = quores(e, c, d); 
		freeall(e);
		return(erg);
		}
	if (b == d) {
		e = callocobject("SYM::quores");
		*e = *b;
		d->c_obj_k(EMPTY);
		erg = quores(e, c, d); 
		freeall(e);
		return(erg);
		}
	if (! d->emptyp())
		d->freeself();
	if (! c->emptyp())
		c->freeself();
	if (emptyp() || b->emptyp())
		return OK;
	if (b->nullp())
		return error("quores: null division");
	if (b->einsp()) {
		/* printf("quores: b is one, just copying\n"); */
		copy(c);
		d->init(s_obj_k());
		d->zero();
		}

	switch (s_obj_k()) {
#ifdef INTEGERTRUE
	case INTEGER:
		return ((INTEGER_OP)this)->quores((INTEGER_OP) b, 
			(INTEGER_OP) c, (INTEGER_OP) d);
#endif
#ifdef LONGINTTRUE
	case LONGINT:
		return(((LONGINT_OP)this)->quores(b, c, d));
#endif
	default:
		{
			printobjectkind();
			return error("quores:wrong first type");
		}
	}
}

#if TEXDOCU
INT SYM_OB::ganzdiv_integral(SYM_OP b, SYM_OP d)
#else
raises an error if the division is not integral.
#endif
{
	SYM_OB r;
	
	quores(b, d, &r);
	if (!r.nullp()) {
		printf("ganzdiv_integral: remainder != 0\n");
		println();
		b->println();
		d->println();
		r.println();
		return ERROR;
		}
	return OK;
}

#if TEXDOCU
INT SYM_OB::ganzdiv(SYM_OP b, SYM_OP d)
#else
$d := this / b$.
#endif
{
	SYM_OP c;
	INT erg = OK;
	
	/* sonderfaelle bei gleichen variablennamen */
	if (this == d) {
		c = callocobject("SYM::ganzdiv");
		*c = *this;
		c_obj_k(EMPTY); 
		erg += c->ganzdiv(b, d); 
		erg += freeall(c); 
		return(erg); 
		}
	if (b == d) {
		c = callocobject("SYM::ganzdiv");
		*c = *b;
		d->c_obj_k(EMPTY); 
		erg += ganzdiv(c,d); 
		erg += freeall(c); 
		return(erg); 
		}

	if (! d->emptyp()) 
		d->freeself();

	/*falls beides leere objecte => d auch leer  */
	if (emptyp() || b->emptyp()) 
		return OK;

	if (b->nullp())	
		return error("ganzdiv: null division");
	if (b->einsp()) 
		return copy(d);



	switch (s_obj_k()) {
#ifdef INTEGERTRUE
	case INTEGER : 
		erg += ((INTEGER_OP)this)->ganzdiv(
			(INTEGER_OP) b, 
			(INTEGER_OP) d); break;
#endif	
#ifdef LONGINTTRUE
	case LONGINT : erg += ((LONGINT_OP)this)->ganzdiv(b, d); break;
#endif
#ifdef BRUCHTRUE
	case BRUCH:
		((BRUCH_OP)this)->scalar_it();
		if (s_obj_k() != BRUCH)
			ganzdiv(b, d);
		break;
#endif
	default: {
			printobjectkind(); 
			return error("SYM::ganzdiv: wrong first type");
			}
		}
	if (erg != OK) {
		error("ganzdiv: error during computation");
		}

	return erg;
}

#if TEXDOCU
INT SYM_OB::fakul(SYM_OP d)
#endif
{
	INT i, j;
	INTEGER_OP a;
	SYM_OB e, f;

	if (s_obj_k() != INTEGER) {
		error("fakul: no INTEGER");
		return ERROR;
		}
	a = (INTEGER_OP) this;
	if (a->s_i() < 0) {
		error("fakul: negativ INTEGER");
		return ERROR;
		}

	if (this == d) {
		SYM_OB c;
		 
		c = *this; 
		c_obj_k(EMPTY);
		c.fakul(d); 
		return(OK);
		}

	if (!d->emptyp()) 
		d->freeself();
	j = a->s_i();
	i = 2;
	((INTEGER_OP)d)->m_i(1);
	for (i = 2; i <= j; i++) {
		((INTEGER_OP)&e)->m_i(i);
		d->mult(&e, &f);
		f.swap(d);
		}
	return OK;
}

#if TEXDOCU
INT SYM_OB::sub(SYM_OP b, SYM_OP c)
#else
$c := this - b$.
#endif
{
	INT erg = OK;
	SYM_OP d;
	INTEGER_OP ai, bi, ci;
	
	if (s_obj_k() == INTEGER &&
		b->s_obj_k() == INTEGER) {
		ai = (INTEGER_OP) this;
		bi = (INTEGER_OP) b;
		ci = (INTEGER_OP) c;
		if ((ai->s_i() < 10000000L) && 
			(ai->s_i() > -10000000L) &&
			(bi->s_i() < 10000000L) && 
			(bi->s_i() > -10000000L)) {
			ci->m_i(ai->s_i() - bi->s_i());
			return OK;
			}
		}
	d = callocobject("SYM::sub");
	erg += b->addinvers(d); 
	erg += add(d, c); 
	erg += freeall(d);
	if (erg != OK) {
		printobjectkind();
		b->printobjectkind();
		error("sub: errors in computing");
		}
	return erg;
}

#undef DEBUG_BEZOUT

#if TEXDOCU
INT bezout(SYM_OP m, SYM_OP n, SYM_OP u, SYM_OP v, SYM_OP g)
#else
$g = \gcd(m,n) = u \cdot m + v \cdot n$.
#endif
{
	SYM_OB m1, n1, q, r;
	SYM_OB u1, u2, u3, v1, v2, v3, tmp1, tmp2;
	INT erg = OK;
	
	u->freeself();
	v->freeself();
	g->freeself();
#ifdef DEBUG_BEZOUT
	printf("m->s_obj_k() = %ld\n", m->s_obj_k());
#endif
	u->init(m->s_obj_k());
	v->init(m->s_obj_k());
	g->init(m->s_obj_k());
	u1.init(m->s_obj_k());
	u2.init(m->s_obj_k());
	u3.init(m->s_obj_k());
	v1.init(m->s_obj_k());
	v2.init(m->s_obj_k());
	v3.init(m->s_obj_k());
	tmp1.init(m->s_obj_k());
	tmp2.init(m->s_obj_k());
	if (m->nullp()) {
		erg += n->copy(g);
		goto l_norm;
		}
	if (n->nullp()) {
		erg += m->copy(g);
		goto l_norm;
		}
	/* from now on m != 0, n != 0 */
	
#ifdef DEBUG_BEZOUT
	printf("bezout: m = ");
	m->println();
	printf("n = ");
	n->println();
#endif
	erg += m->quores(n, &q, &r);
#ifdef DEBUG_BEZOUT
	printf("q = ");
	q.println();
	printf("r = ");
	r.println();
#endif
	if (q.nullp()) {
		erg += bezout(n, m, v, u, g);
		goto l_norm;
		}
	/* from now on: deg(m) >= deg(n)  */
	if (n->nullp()) {
		if (m->nullp()) {
			erg += g->zero();
			erg += u->zero();
			erg += v->zero();
			goto l_exit;
			}
		erg += m->copy(g);
		erg += u->one();
		erg += v->zero();
		goto l_exit;
		}
	/* m1 = m; u1 = 1; v1 = 0;
	 * n1 = n; u2 = 0; v2 = 1; */
	erg += m->copy(&m1);
	erg += u1.one();
	erg += v1.zero();
	erg += n->copy(&n1);
	erg += u2.zero();
	erg += v2.one();
	while (TRUE) {
#ifdef DEBUG_BEZOUT
		printf("bezout: m1 = ");
		m1.println();
		printf("n1 = ");
		n1.println();
#endif
		erg += m1.quores(&n1, &q, &r);
#ifdef DEBUG_BEZOUT
		printf("q = ");
		q.println();
		printf("r = ");
		r.println();
#endif
		if (r.nullp()) {
			/* u = u2;
			 * v = v2;
			 * g = n1; */
			erg += u2.copy(u);
			erg += v2.copy(v);
			erg += n1.copy(g);
			break;
			}
		/* u3 = u1 - q * u2;
		 * v3 = v1 - q * v2; */
		erg += q.mult(&u2, &tmp1);
		erg += q.mult(&v2, &tmp2);
		erg += tmp1.addinvers_apply();
		erg += tmp2.addinvers_apply();
		erg += u1.add(&tmp1, &u3);
		erg += v1.add(&tmp2, &v3);
		/* m1 = n1; n1 = r;
		 * u1 = u2; u2 = u3;
		 * v1 = v2; v2 = v3; */
		erg += n1.copy(&m1);
		erg += r.copy(&n1);
		erg += u2.copy(&u1);
		erg += u3.copy(&u2);
		erg += v2.copy(&v1);
		erg += v3.copy(&v2);
		}
l_norm:
	if (!g->nullp()) {
		tmp1.init(m->s_obj_k());
		erg += tmp1.one();
#ifdef DEBUG_BEZOUT
		printf("g != 0, normiere:");
		printf("tmp1 = 1 = ");
		tmp1.println();
		printf("g = ");
		g->println();
#endif
		erg += tmp1.quores(g, &q, &r);
#ifdef DEBUG_BEZOUT
		printf("q = ");
		q.println();
		printf("r = ");
		r.println();
#endif
		if (r.nullp()) {
			/* g is a unit and q it's inverse */
			erg += g->mult(&q, &tmp1);
				/* make g to one */
			tmp1.copy(g);
			erg += u->mult(&q, &tmp1);
			tmp1.copy(u);
			erg += v->mult(&q, &tmp1);
			tmp1.copy(v);
			}
		}
l_exit:
	return erg;
}

#if TEXDOCU
INT inverse_mod(SYM_OP a, SYM_OP b, SYM_OP m, INT *f_notinvertible)
#endif
{
	INT erg = OK;
	SYM_OB u, v, g;
	
	erg += bezout(a, m, &u, &v, &g);
	if (!g.einsp())
		*f_notinvertible = TRUE;
	else {
		*f_notinvertible = FALSE;
		/* 1 = u * a + v * m */
		erg += u.copy(b);
		}
	return erg;
}

#undef DEBUG_PPP

#if TEXDOCU
INT prime_power_parts(SYM_OP g, VECTOR_OP gpp, 
	INT type, void *data)
#endif
/* es wird an den Vektor gpp angehaengt. */
{
	INT order, i, j, k, k1, len, ii;
	SYM_OB g1;
	VECTOR_OB op, oe;
	
	do_order(g, &order, type, data);
#ifdef DEBUG_PPP
	printf("prime_power_parts(): ");
	do_print(g, type, data);
	printf(" order = %ld\n", order);
#endif
	factor_integer(order, &op, &oe);
	len = op.s_li();
	if (len == 1) {
		gpp->inc();
		g->copy(gpp->s_i(gpp->s_li() - 1));
		return OK;
		}
	for (i = 0; i < len; i++) {
		k = 1;
		for (j = 0; j < len; j++) {
			if (j == i)
				continue;
			k1 = i_power_j(op.s_ii(j), oe.s_ii(j));
			k *= k1;
			}
		do_power(g, k, &g1, type, data);
#ifdef DEBUG_PPP
		printf("prime_power_parts(): k = %ld g1 = ", k);
		do_print(&g1, type, data);
		printf("\n");
#endif
		for (ii = 0; ii < gpp->s_li(); ii++) {
			if (do_comp(gpp->s_i(ii), &g1, type, data) == 0)
				goto weiter;
			}
		gpp->inc();
		g1.swap(gpp->s_i(gpp->s_li() - 1));
weiter:
			;
		}
	return OK;
}

#if TEXDOCU
INT vec_to_vec_pp(VECTOR_OP V, VECTOR_OP Vpp, 
	INT type, void *data)
#endif
{
	INT i, len, erg = OK;
	
	len = V->s_li();
	Vpp->m_il(0);
	for (i = 0; i < len; i++) {
		erg += prime_power_parts(
			V->s_i(i), Vpp, type, data);
		}
	return erg;
}

/* type 0: ordinary SYMMETRICA mult()
 * (1: reverse mult (PERMUTATIONs) - no longer needed)
 * 2: ENUM mult (needs #define ENUM_TRUE) 
 * 3: iperm (data is CHUNK_MEMH pointer)
 * 4: FG_OB::gt_*
 */

#if TEXDOCU
INT do_mult(SYM_OP a, SYM_OP b, SYM_OP c, INT type, void *data)
#endif
{
	INT erg = OK;
	
	if (type == DO_TYPE_SYM || 
		type == DO_TYPE_PERM) {
		erg += a->mult(b, c);
		}

#ifdef SOLVABLE_TRUE
	else if (type == DO_TYPE_FG) {
		FG_OP G;
		
#ifdef FG_DEBUG
		printf("FG: do_mult()"); fflush(stdout);
#endif
		G = (FG_OP) data;
		G->mult_op((INTEGER_OP) a, (INTEGER_OP) b, (INTEGER_OP) c);
		}
#endif
	else {
		Srfs("do_mult", "nyi");
		printf("type = %ld\n", type);
		return ERROR;
		}
	return erg;
}

#if TEXDOCU
INT do_invers(SYM_OP a, SYM_OP b, INT type, void *data)
#endif
{
	INT erg = OK;
	
	if (type == DO_TYPE_SYM || 
		type == DO_TYPE_PERM) {
		erg += a->invers(b);
		}
#ifdef SOLVABLE_TRUE
	else if (type == DO_TYPE_FG) {
		FG_OP G;
		
#ifdef FG_DEBUG
		printf("FG: do_invers()"); fflush(stdout);
#endif
		G = (FG_OP) data;
		G->inv_op((INTEGER_OP) a, (INTEGER_OP) b);
		}
#endif
	else {
		Srfs("do_invers", "nyi");
		printf("type = %ld\n", type);
		return ERROR;
		}
	return erg;
}

#if TEXDOCU
INT do_einsp(SYM_OP a, INT type, void *data)
#endif
{
	INT res;
	
	if (type == DO_TYPE_SYM || 
		type == DO_TYPE_PERM) {
		res = a->einsp();
		}
#ifdef SOLVABLE_TRUE
	else if (type == DO_TYPE_FG) {
		FG_OP G;
		
#ifdef FG_DEBUG
		printf("FG: do_einsp()"); fflush(stdout);
#endif
		G = (FG_OP) data;
		return G->onep_op((INTEGER_OP) a);
		}
#endif
	else {
		Srfs("do_einsp", "nyi");
		printf("type = %ld\n", type);
		return ERROR;
		}
	return res;
}

#if TEXDOCU
INT do_one(SYM_OP a, INT type, void *data)
#endif
{
	INT erg = OK;
	
	if (type == DO_TYPE_SYM || 
		type == DO_TYPE_PERM) {
		erg += a->one();
		return erg;
		}
#ifdef SOLVABLE_TRUE
	else if (type == DO_TYPE_FG) {
		FG_OP G;
		
#ifdef FG_DEBUG
		printf("FG: do_one()");
#endif
		G = (FG_OP) data;
		G->one_op((INTEGER_OP) a);
#ifdef FG_DEBUG
		printf("a = %ld\n", ((INTEGER_OP) a)->s_i()); fflush(stdout);
#endif
		}
#endif
	else {
		Srfs("do_one", "nyi");
		printf("type = %ld\n", type);
		return ERROR;
		}
	return erg;
}

#if TEXDOCU
INT do_order(SYM_OP a, INT *order, INT type, void *data)
#endif
{
	INT erg = OK;
	
	if (type == DO_TYPE_SYM || 
		type == DO_TYPE_PERM) {
		erg += a->order(order);
		return erg;
		}
#ifdef SOLVABLE_TRUE
	else if (type == DO_TYPE_FG) {
		FG_OP G;
		
#ifdef FG_DEBUG
		printf("FG: do_order()"); fflush(stdout);
#endif
		G = (FG_OP) data;
		G->gt_order(((INTEGER_OP) a)->s_i(), order);
		}
#endif
	else {
		Srfs("do_order", "nyi");
		printf("type = %ld\n", type);
		return ERROR;
		}
	return erg;
}

#if TEXDOCU
INT do_order_if_prime(SYM_OP a, INT *ord, INT type, void *data)
#endif
{
	INT erg = OK;
	
	if (type == DO_TYPE_SYM || 
		type == DO_TYPE_PERM) {
		erg += a->order_if_prime(ord);
		return erg;
		}
#ifdef SOLVABLE_TRUE
	else if (type == DO_TYPE_FG) {
		FG_OP G;
		
#ifdef FG_DEBUG
		printf("FG: do_order_if_prime()"); fflush(stdout);
#endif
		G = (FG_OP) data;
		G->gt_order_if_prime(((INTEGER_OP) a)->s_i(), ord);
		}
#endif
	else {
		Srfs("do_order_if_prime", "nyi");
		printf("type = %ld\n", type);
		return ERROR;
		}
	return erg;
}

#if TEXDOCU
INT do_order_if_prime_power(SYM_OP a, 
	INT *ord, INT *prime, INT *k, 
	INT type, void *data)
#endif
{
	INT erg = OK;
	
	if (type == DO_TYPE_SYM || 
		type == DO_TYPE_PERM) {
		erg += a->order_if_prime_power(
			ord, prime, k);
		return erg;
		}
#ifdef SOLVABLE_TRUE
	else if (type == DO_TYPE_FG) {
		FG_OP G;
		
#ifdef FG_DEBUG
		printf("FG: do_order_if_prime_power()"); fflush(stdout);
#endif
		G = (FG_OP) data;
		G->gt_order_if_prime_power(((INTEGER_OP) a)->s_i(), ord, prime, k);
		}
#endif
	else {
		Srfs("do_order_if_prime_power", "nyi");
		printf("type = %ld\n", type);
		return ERROR;
		}
	return erg;
}

#if TEXDOCU
INT do_power(SYM_OP a, INT i, SYM_OP b, INT type, void *data)
#endif
{
	INT j;
	SYM_OB c, d;

	a->copy(&c);
	if (i == 0) {
		do_one(&c, type, data);
		}
	else {
		j = 1;
		while (j < i) {
			do_mult(a, &c, &d, type, data);
			d.swap(&c);
			j++;
			}
		}
	c.swap(b);
	return OK;
}

#if TEXDOCU
INT do_print(SYM_OP a, INT type, void *data)
#endif
{
	INT erg = OK;
	
	if (type == DO_TYPE_SYM || 
		type == DO_TYPE_PERM || 
		type == DO_TYPE_FG) {
		erg += a->print();
		return erg;
		}
	Srfs("do_print", "nyi");
	printf("type = %ld\n", type);
	return ERROR;
}

#if TEXDOCU
INT do_sprint(SYM_OP a, BYTE *str, INT type, void *data)
#endif
{
	INT erg = OK;
	
	if (type == DO_TYPE_SYM || 
		type == DO_TYPE_PERM || 
		type == DO_TYPE_FG) {
		erg += a->sprint(str);
		return erg;
		}
	Srfs("do_sprint", "nyi");
	printf("type = %ld\n", type);
	return ERROR;
}

#if TEXDOCU
INT do_latex(SYM_OP a, FILE *fp, INT type, void *data)
#endif
{
	INT erg = OK;
	
	if (type == DO_TYPE_SYM || 
		type == DO_TYPE_PERM || 
		type == DO_TYPE_FG) {
		erg += a->latex(fp);
		return erg;
		}
	Srfs("do_latex", "nyi");
	printf("type = %ld\n", type);
	return ERROR;
}


#if TEXDOCU
INT do_image_of(SYM_OP a, INT i, INT type, void *data)
#endif
{
	INT res = -1;
	
	if (type == DO_TYPE_SYM || 
		type == DO_TYPE_PERM) {
#ifdef PERMTRUE
		PERMUTATION_OP p;
		
		if (a->s_obj_k() != PERMUTATION) {
			printf("do_image_of(): "
				"not a permutation\n");
			return -1;
			}
		p = (PERMUTATION_OP) a;
		res = p->s_ii(i) - 1;
#endif
		return res;
		}
	Srfs("do_image_of", "nyi");
	printf("type = %ld\n", type);
	return -1;
}

#if TEXDOCU
INT do_comp(SYM_OP a, SYM_OP b, INT type, void *data)
#endif
{
	INT res;
	
	if (type == DO_TYPE_SYM || 
		type == DO_TYPE_PERM || 
		type == DO_TYPE_FG) {
		res = a->sym_comp(b);
		return res;
		}
	Srfs("do_comp", "nyi");
	printf("type = %ld\n", type);
	return -1;
}

#if TEXDOCU
INT do_copy(SYM_OP a, SYM_OP b, INT type, void *data)
#endif
{
	INT erg = OK;
	
	if (type == DO_TYPE_SYM || 
		type == DO_TYPE_PERM || 
		type == DO_TYPE_FG) {
		erg += a->copy(b);
		return erg;
		}
	Srfs("do_copy", "nyi");
	printf("type = %ld\n", type);
	return ERROR;
}


