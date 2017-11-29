/* poly.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef POLYTRUE

#include <DISCRETA/poly.h>
#include <DISCRETA/bruch.h>
#include <DISCRETA/part.h>
#include <DISCRETA/lb.h>
#include <DISCRETA/fga.h> // for reduce_generators_labra

#if TEXDOCU
INT POLYNOM_OB::Print()
#else
#endif
{
	POLYNOM_OP q;
	
	q = this;
	while (q != NIL) {
		if (q->s_obj_k() == EMPTY) {
			printf("#\n");
			}
		else {
			q->s_mo()->println();
			}
		q = q->s_n();
		}
	return OK;
}

#if TEXDOCU
INT POLYNOM_OB::b_mn(MONOM_OP m, POLYNOM_OP next)
#else
build\_monom\_next\_polynom: put the monomial and the next pointer 
together in a polynomial.
#endif
{
	INT erg;
	
	erg = ((LIST_OP)this)->b_sn(m, (LIST_OP) next);
	c_obj_k(POLYNOM);
	return erg;
}

#if TEXDOCU
INT POLYNOM_OB::b_skn(VECTOR_OP self, SYM_OP koeff, POLYNOM_OP next)
#else
build\_self\_coefficient\_next\_polynom: 
make a new monomial containing koeff and self and 
place this monomial together with the next pointer into a poynomial.
#endif
{
	INT erg = OK;
	
	erg += b_mn((MONOM_OP) callocobject("POLY::b_skn"), next);
	erg += s_mo()->b_sk(self, koeff);
	return erg;
}

#if TEXDOCU
INT POLYNOM_OB::b_skin(VECTOR_OP self, INT koeff, POLYNOM_OP next)
#else
like b\_skn but here the coeffitient is an integer.
#endif
{
	INT erg = OK;
	
	erg += b_mn((MONOM_OP) callocobject("POLY::b_skin"), next);
	erg += s_mo()->b_ski(self, koeff);
	return erg;
}

#if TEXDOCU
INT POLYNOM_OB::m_skn(VECTOR_OP self, SYM_OP koeff, POLYNOM_OP next)
#else
make\_self\_coefficient\_next\_polynomial: make a new polynomial 
using a copy of self, koeff. 
#endif
{
	INT erg = OK;
	
	erg += b_mn((MONOM_OP) callocobject("POLY::m_skn"), next);
	erg += s_mo()->m_sk(self, koeff);
	return erg;
}

#if TEXDOCU
INT POLYNOM_OB::m_skin(VECTOR_OP self, INT koeff, POLYNOM_OP next)
#else
like m\_skin but here koeff is an integer.
#endif
{
	INT erg = OK;
	
	erg += b_mn((MONOM_OP) callocobject("POLY::m_skin"), next);
	erg += s_mo()->m_ski(self, koeff);
	return erg;
}

#if TEXDOCU
INT POLYNOM_OB::m_iindex(INT i)
#else
generates the polynomial $1 \cdot x^i$.
#endif
{
	return m_iindex_iexponent(i, 1);
}

#if TEXDOCU
INT POLYNOM_OB::m_iindex_iexponent(INT i, INT j)
#else
generates the polynomial $1 \cdot x^j$.
#endif
{
	if (i < 0)
		return error("POLY::m_iindex_iexponent()| i < 0");
	freeself();
	b_skn((VECTOR_OP) callocobject("POLY::m_iindex_iexponent"), callocobject("POLY::m_iindex_iexponent"), NIL);
	s_s()->m_il_n(i + 1);
	((INTEGER_OP) s_k())->m_i(1);
	s_s()->m_ii(i, j);
	return OK;
}

#if TEXDOCU
INT POLYNOM_OB::nullp()
#else
Returns TRUE if the polynomial is the zero polynomial.
#endif
{
	POLYNOM_OP z = this;
	
	if (((LIST_OP)z)->emptyp()) 
		return TRUE;
	while (z != NULL) {
		if (! z->s_k()->nullp())
			return FALSE;
		z = z->s_n();
		}
	return TRUE;
}


#if TEXDOCU
INT POLYNOM_OB::einsp()
#else
Returns TRUE if the polynomial is the one polynomial.
#endif
{
	POLYNOM_OP z = this;
	VECTOR_OP s;
	INT i;
	
	if (((LIST_OP)z)->emptyp())
		return FALSE;
	if (! z->s_k()->einsp())
		return FALSE;
	/* now check whether self == 000000 */
	s = s_s();
	for (i = 0; i < s->s_li(); i++)
		if (s->s_ii(i) != 0)
			return FALSE;
	z = z->s_n();
	if (z != NULL)
		return FALSE;
	return TRUE;
}


#if TEXDOCU
INT POLYNOM_OB::m_scalar(SYM_OP a)
#else
make the polynomial $a \cdot x^0$.
#endif
{
	VECTOR_OP s;

	b_skn((VECTOR_OP) callocobject("POLY::m_scalar"), callocobject("POLY::m_scalar"), NULL);
	a->copy(s_k());
	s = s_s();
	s->m_il_n(1);
	return OK;
}

#if TEXDOCU
INT POLYNOM_OB::mult_scalar(SYM_OP a, POLYNOM_OP erg)
#else
multiply all coefficients by a.
#endif
/* before: mult_scalar_polynom(OP a, OP poly, OP erg) */
{
#if FALSE
	if ((nullp(a)) || (nullp()))
		return m_i_i(0L, erg);
#endif
	/* return trans2formlist(a,poly,erg,mult); */
	return ((LIST_OP)this)->transform_by(a, (LIST_OP)erg, mult_abc);
}

#undef DEBUG_MULT_POLYNOM

#if TEXDOCU
INT POLYNOM_OB::mult_polynom(POLYNOM_OP zwei, POLYNOM_OP c)
#else
c := this * zwei. zwei must be a polynomial.
#endif
/* mult_polynom_polynom(OP eins, OP zwei, OP c) */
{
	POLYNOM_OP z, ez, zz;
	
#ifdef DEBUG_MULT_POLYNOM
	printf("in mult_polynom():\n");
	println();
	zwei->println();
	fflush(stdout);
#endif
	/* c freigeben ! */
	zz = zwei;
	while (zz != NULL) {
		z = (POLYNOM_OP) callocobject("POLY::mult_polynom");
		copy(z);

		ez = z;
		while (ez != NULL) {
#ifdef DEBUG_MULT_POLYNOM
			printf("ez->s = ");
			ez->s_s()->println();
			printf("zz->s = ");
			zz->s_s()->println();
			fflush(stdout);
#endif
			ez->s_s()->add_apply(zz->s_s());
#ifdef DEBUG_MULT_POLYNOM
			printf("ez->s = ");
			ez->s_s()->println();
			fflush(stdout);
#endif
			zz->s_k()->mult_apply(ez->s_k());
			ez = ez->s_n();
			}
#ifdef DEBUG_MULT_POLYNOM
			printf("vor insert_list() z = ");
			z->println();
			printf("c = ");
			c->println();
			fflush(stdout);
#endif
		((LIST_OP)c)->insert_list((LIST_OP)z, add_koeff, comp_monomvector_monomvector);
#ifdef DEBUG_MULT_POLYNOM
		printf("nach insert_list() c = ");
		c->println();
		fflush(stdout);
#endif
		zz = zz->s_n();
		}
#ifdef DEBUG_MULT_POLYNOM
	printf("ende mult_polynom():\n");
	c->println();
	fflush(stdout);
#endif
	return(OK);
}

#if TEXDOCU
INT POLYNOM_OB::mult(SYM_OP b, POLYNOM_OP d)
#else
d := this * b. b can be of various types.
#endif
{
	INT erg = OK;
	
#if FALSE
	if (nullp(a))
		return m_i_i(0L,d);
#endif

	switch (b->s_obj_k()) {
	case INTEGER:
	case LONGINT:
#ifdef BRUCHTRUE
	case BRUCH: 
#endif /* BRUCHTRUE */
		erg += mult_scalar(b, (POLYNOM_OP) d);
		break;
	case POLYNOM: erg += mult_polynom((POLYNOM_OP)b, (POLYNOM_OP)d);
		break;
#ifdef SCHUBERTTRUE
	case SCHUBERT: erg+=mult_schubert_polynom(b,a,d);
		break;
#endif /* SCHUBERTTRUE */
#ifdef SCHURTRUE
	case SCHUR: erg += mult_scalar_schur(a,b,d); break;
#endif /* SCHURTRUE */
#ifdef MATRIXTRUE
	/* case MATRIX: erg += mult_scalar_matrix(a,b,d);!!! */
		break;
#endif /* MATRIXTRUE */
#ifdef MONOMTRUE
	/* case MONOM: erg += mult_scalar(a, d); break;!!! */
#endif /* MONOMTRUE */
	default: 
		b->printobjectkind();
		erg += error("POLYNOM::mult: wrong second type"); 
	}
	return erg;
}

#if TEXDOCU
INT POLYNOM_OB::mult_apply_scalar(SYM_OP a)
#endif
/* before: mult_apply_scalar_polynom(OP a, OP b) */
{
	POLYNOM_OP z = this;
	
#if FALSE
	if (nullp(a))
		{
		return m_i_i(0L,b);
		}
#endif
	while (z != NULL) {
		s_k()->mult_apply(a);
		z = z->s_n();
		}
	return(OK);
}

#if TEXDOCU
INT POLYNOM_OB::mult_apply_polynom(POLYNOM_OP b)
#else
b := this * b. b must be a polynomial. 
#endif
/* before: mult_apply_polynom_polynom(OP a, OP b) */
{
	POLYNOM_OP c = (POLYNOM_OP) callocobject("POLY::mult_apply_polynom");
	INT erg = OK;
	
	*c = *b;
	b->c_obj_k(EMPTY);
	erg += mult_polynom(c, b); 
	erg += freeall(c);
	if (erg != OK) {
		error("POLY::mult_apply_polynom: error during computation");
		}
	return erg;
}

#if TEXDOCU
INT POLYNOM_OB::mult_apply(SYM_OP b)
#else
b := this * b. b can be of various types.
#endif
/* before: mult_apply_polynom(OP a, OP b) */
{
	INT erg = OK;
	
#if FALSE
	if (nullp(b))
		return m_i_i(0L,b);
#endif
	switch (b->s_obj_k()) {
		case BRUCH:
		case LONGINT:
		case INTEGER:  /* a ist scalar, b ist poly */
			/* erg+= mult_apply_polynom_scalar(a,b); break; */
			erg += mult_apply_scalar(b); break;
		case POLYNOM: 
			erg+= mult_apply_polynom((POLYNOM_OP) b); break;
		default:
			{
			SYM_OP c = callocobject("POLY::mult_apply");
			*c = *b;
			b->c_obj_k(EMPTY);
			erg = mult(c, (POLYNOM_OP) b);
			if (erg == ERROR) {
				c->printobjectkind();
				error("POLYNOM::mult_apply:wrong second type");
				}
			freeall(c); return erg;
			}
		}
	if (erg != OK) {
		error("POLYNOM::mult_apply: error during computation");
		}
	return erg;
}


#if TEXDOCU
INT POLYNOM_OB::addinvers(POLYNOM_OP b)
#else
b := -this.
#endif
/* before: addinvers_polynom(OP a, OP b) */
{ 
	return(((LIST_OP)this)->transform((LIST_OP)b, addinvers_ab)); 
}

#if TEXDOCU
INT POLYNOM_OB::addinvers_apply()
#else
this := - this.
#endif
/* before: addinvers_apply_polynom(OP a) */
{ 
	return(((LIST_OP)this)->transform_apply(addinvers_apply_a)); 
}

#if TEXDOCU
INT POLYNOM_OB::add_polynom(POLYNOM_OP b, POLYNOM_OP c)
#endif
/* before: add_polynom_polynom(OP a, OP b, OP c) */
{
	INT erg = OK;
	POLYNOM_OP d = (POLYNOM_OP) callocobject("POLY::add_polynom");
	
	if (! c->emptyp()) 
		erg += c->freeself();
	erg += copy(d); 
	erg += b->copy(c);
	

	((LIST_OP)c)->insert((LIST_OP)d, add_koeff, comp_monomvector_monomvector);

#if FALSE
	if (c->nullp())
		if (c->s_mo() != NULL) {
				b_sn_l(NULL,NULL,c);
				C_O_K(c,S_O_K(a));
				}
#endif
	return erg;
}
 
#if TEXDOCU
INT POLYNOM_OB::add_scalar(SYM_OP a, POLYNOM_OP c)
#endif
/* before: add_scalar_polynom(OP a, OP b, OP c) */
/* a scalar, b polynom, c result */
{
	INT erg = OK;
	POLYNOM_OP d;

	if (! c->emptyp()) 
		erg += c->freeself();
	
	if (a->nullp())
		return copy(c);

	d = (POLYNOM_OP) callocobject("POLY::add_scalar");
	erg += d->m_scalar(a);
	erg += add_polynom(d, c);
	erg += freeall(d);
	return erg;
}

#if TEXDOCU
INT POLYNOM_OB::add(SYM_OP b, POLYNOM_OP c)
#endif
/* before: add_polynom(OP a, OP b, OP c) */
{
	INT erg = OK;
	
	if (nullp())
		return b->copy(c);
	if (b->nullp())
		return copy(c);

	switch (b->s_obj_k()) {
	case INTEGER:
	case LONGINT:
	case BRUCH: erg += add_scalar(b, c); break;
	case GRAL:
	case POLYNOM: erg += add_polynom((POLYNOM_OP) b, c); break;
	default: {
		b->printobjectkind();
		return error("POLYNOM::add: wrong second type");
		}
	}
	if (erg != OK) {
		error("POLYNOM::add: error during computation");
		}	
	return erg;
}

#if TEXDOCU
INT POLYNOM_OB::add_apply_polynom(POLYNOM_OP b)
#endif
/* before: add_apply_polynom_polynom(OP a, OP b) */
{
	INT erg = OK;
	POLYNOM_OP c;
	
	c = (POLYNOM_OP) callocobject("POLY::add_apply_polynom");
	erg += copy(c);
	((LIST_OP)b)->insert_list((LIST_OP)c, add_koeff, comp_monomvector_monomvector);
	return erg;
}

#if TEXDOCU
INT POLYNOM_OB::add_apply_scalar(SYM_OP a)
#endif
/* before: add_apply_polynom_scalar(OP a, OP b) 
 * a polynom, b scalar. berechnet b := a + b. (?) !!! */
{
	POLYNOM_OP c = (POLYNOM_OP) callocobject("POLY::add_apply_scalar");
	INT erg = OK;
	
	erg += copy(c);
	erg += c->add(a, this);
	erg += freeall(c);
	return erg;
}

#if TEXDOCU
INT POLYNOM_OB::add_apply(SYM_OP b)
#endif
{
	INT erg = OK;
	
#if FALSE
	if ((b->emptyp)) || b->nullp())
		return(copy_polynom(a,b));
#endif
	switch (b->s_obj_k()) {
		case BRUCH:
		case LONGINT:
		case INTEGER: erg += add_apply_scalar(b); break;
		case POLYNOM:  erg += add_apply_polynom((POLYNOM_OP) b); break;
#ifdef SCHUBERTTRUE
		case SCHUBERT:  erg += add_apply_polynom_schubert(a,b);break;
#endif /* SCHUBERTTRUE */
		default:
			b->printobjectkind();
			return error("POLYNOM::add_apply: wrong second type");
		}
	if (erg != OK)
		error("POLYNOM::add_apply: error during computation");
	return erg;
}

#if TEXDOCU
#else
\subsection{Cycle-Indices}
The cycle index of a finite permutation group $G$ acting on a set $X$ of $n$ elements is 
defined as 
\[
C(G, \, X) \; := \; \frac{1}{|G|} \sum_{g \in G\atop c(g) = (a_1,\ldots, a_n)} 
\prod_{i=1}^n Z_i^{a_i},
\]
where $c(g) \, = \, (a_1, a_2, \ldots, a_n)$ 
is the {\em cycle type} of the permutation $g$. 
As the cycle type is constant on conjugacy classes, 
one gets the following refinement of this formula. 
Here $\cC = \{g_1,\ldots,g_c\}$ is a transversal of the classes of $G$.
\[
C(G, \, X) \; := \; \frac{1}{|G|} \sum_{i=1\atop c(g_i) = (a_1,\ldots, a_n)}^c 
|C^G(g_i)| \prod_{i=1}^n Z_i^{a_i},
\]
#endif

#if TEXDOCU
INT POLYNOM_OB::cycle_ind_Cn(INT n)
#else
Computes the cycle index of $C_n$, the cyclic group of order $n$.


Let now $G = C_n = \la g \ra$ with $g^n = 1$.
Then $G$ is abelian and therefore its classes 
consist of one element each. The cycle structure of 
any $g^i$, $1 \le i \le n$ consists of $d$ cycles 
of length $\frac{n}{d}$ with $d =  \gcd (n, i)$. 
The element order of $g^i$ is thus $\frac{n}{d}$.
The number of integers $i \le n$ which have $\gcd (n,i) = d$ 
is $\varphi(d)$ where $\varphi$ is Eulers function. 
So
\begin{align*}
C(C_n, \, X)  
=& \; \sum_{d \mid n} \varphi(d) Z_{d}^{n/d}
\end{align*}
#endif
{
	if (n < 1)
		return error("POLY::cycle_ind_Cn() n < 1");
	freeself();
	if (n == 1)
		return m_iindex(0);
	init(POLYNOM);
#ifdef BRUCHTRUE
	{
		POLYNOM_OB p;
		BRUCH_OB b;
		VECTOR_OB V;
		INT d;

		for (d = 1; d <= n; d++) {
			if (n % d)
				continue;
			/* now d divides n */
			b.m_ioiu(eulerfunc(d), n);
			p.m_skn(&V, &b, NIL);
			p.s_s()->m_il_n(n);
			p.s_s()->m_ii(d - 1, n / d);
			p.add_apply(this);
			}
		
	}
#else
	return error("cycle_ind_Cn(): BRUCH not available");
#endif
	return OK;
}

#if TEXDOCU
INT POLYNOM_OB::cycle_ind_Dn(INT n)
#else
In order to compute the cycle index of $D_n$, the dihedral group 
on $n$ elements (of order $2n$) we note that $C_n$ is a subgroup 
of index two. $\frac{1}{2}C(C_n, X)$ therefore counts 
the element cycles of the elements lying in $C_n$. 
For the elements of $D_n \backslash C_n$ we have to divide 
two cases:
\begin{enumerate}
\item $2 \mid n$: \\
In this case, $n / 4$ elements consist of $n / 2$ cycles of length 2
(forming one conjugacy class), the other $n/4$ elements 
have 2 fixpoints and $(n-2)/2$ 2--cycles (again forming one 
conjugacy class). This gives
\begin{align*}
C(D_n, \, X)  
=& \; \frac{1}{2} C(C_n, \, X) + 
\frac{1}{4} Z_2^{n/2} + \frac{1}{4} Z_1^2 Z_2^{(n-2)/2}.
\end{align*}
\item $2 \nmid  n$:\\
In this case, $n/2$ elements of $D_n \backslash C_n$ 
have 1 fixpoint and $(n-1)/2$ cycles of length two 
(forming one conjugacy class). 
Thus
\begin{align*}
C(D_n, \, X) 
=& \; 
\frac{1}{2} C(C_n, \, X) + 
\frac{1}{2} Z_1 Z_2^{(n-1)/2}
\end{align*}
for $2 \nmid n$.
\end{enumerate}
#endif
{
	if (n < 1)
		return error("POLY::cycle_ind_Dn() n < 1");
	freeself();
	if (n == 1)
		return m_iindex(0);
	init(POLYNOM);
#ifdef BRUCHTRUE
	{
		BRUCH_OB halb;
		POLYNOM_OB p, q;
		BRUCH_OB b;
		VECTOR_OB V;
		
		/* Cn is of index 2 in Dn, so 
		 * multiply by 1/2: */
		cycle_ind_Cn(n);
		halb.m_ioiu(1, 2);
		((SYM_OP) &halb)->mult(this, &q);

		if (EVEN(n)) {
			b.m_ioiu(1, 4);
			p.m_skn(&V, &b, NIL);
			p.s_s()->m_il_n(n);
			p.s_s()->m_ii(1, (n >> 1));
				/* n/2 2-cycles in this element */
			p.add_apply(&q);
			
			b.m_ioiu(1, 4);
			p.m_skn(&V, &b, NIL);
			p.s_s()->m_il_n(n);
			p.s_s()->m_ii(0, 2);
				/* 2 fix points and */
			p.s_s()->m_ii(1, ((n - 2) >> 1));
				/* (n-2)/2 2-cycles in this element */
			p.add_apply(&q);
			}
		else {
			b.m_ioiu(1, 2);
			p.m_skn(&V, &b, NIL);
			p.s_s()->m_il_n(n);
			p.s_s()->m_ii(0, 1);
				/* 1 fix point and */
			p.s_s()->m_ii(1, ((n - 1) >> 1));
				/* (n-1)/2 2-cycles in this element */
			p.add_apply(&q);
			}
		swap(&q);
	}
#else
	return error("cycle_ind_Dn(): BRUCH not available");
#endif
	return OK;
}

#if TEXDOCU
INT POLYNOM_OB::cycle_ind_Sn(INT n)
#else
The conjugacy classes of $S_n$ are labelled by partitions. 
For each partition $\alpha \vdash n$, let $\pi_\alpha$ be an element 
of cycle type $\alpha$. 
Thus
\begin{align*}
C(S_n, \, X) 
=& \;
\frac{1}{n!} 
\sum_{\alpha = (a_1,\ldots, a_n) \vdash n} 
|C^{S_n}(\pi_\alpha)| \prod_{i=1}^n Z_i^{a_i} \\
=& \;
\sum_{\alpha = (a_1,\ldots, a_n) \vdash n} 
\frac{1}{|C_{S_n}(\pi_\alpha)|}
\prod_{i=1}^n Z_i^{a_i} 
\end{align*}
#endif
{
	if (n < 1)
		return error("POLY::cycle_ind_Sn() n < 1");
	freeself();
	if (n == 1)
		return m_iindex(0);
	init(POLYNOM);
	{
		PARTITION_OB type_v, type_e;
		POLYNOM_OB p;
		SYM_OB z, z1;
		INT i, j, k;
		
		type_v.first(n);
		do {
			type_v.t_VECTOR_EXPONENT(&type_e);
			p.b_skn((VECTOR_OP) callocobject("POLY::cycle_ind_Sn"), 
				callocobject("POLY::cycle_ind_Sn"), NIL);
			type_e.s_s()->copy(p.s_s());

			/* the coefficient:
			 * class_len / group_order = 
			 * 1 / centralizer order. */
			
			/* compute centralizer order into z: */
			z.m_i_i(1);
			for (i = 1; i <= n; i++) {
				j = p.s_s_iii(i - 1);
				/* i = cycle length, 
				 * j = multiplicity. */
				z1.fakul_int(j);
				z1.mult_apply(&z);
				z1.m_i_i(i);
				for (k = 0; k < j; k++)
					z1.mult_apply(&z);
				}
			z.invers(p.s_k());
			// p.println();	
			p.add_apply(this);
			} while (type_v.next_VECTOR(&type_v));
	}
	return OK;
}

#if TEXDOCU
INT POLYNOM_OB::cycle_ind_An(INT n)
#else
The conjugacy classes of $A_n$ are either conjugacy classes 
of even elements of $S_n$ or pairs of classes (of equal size)
which together form one conjugacy class of even elements of $S_n$.
In is well known that a conjugacy class of even elements of $S_n$ 
splits into two classes of equal size in $A_n$ if and only if the 
cycle type of its elements consists only of odd cycles 
and all cycles haved different length. 
For example, $S_5$ has one conjugacy class of 5-cycles $(1,2,3,4,5)$ (with 24 elements) 
but $A_5$ has 2 classes of 5-cycles, each class containing 12 elements.
Note, however, 
that in order to compute the cycle type one only studies cycle types 
and not conjugacy classes, therefore the splitting of classes does not matter.

In order to obtain all elements of $A_n$ we loop over all 
partitions $\alpha \vdash n$ with $\pi_\alpha$ even. 
Note that $\pi_\alpha$ is even if and only if $\sum_{j=2 \atop 2 \mid j}^n a_j$ 
is an even integer.
So we get
\begin{align*}
C(A_n, \, X) 
=& \; 
\frac{2}{n!} 
\sum_{\alpha = (a_1,\ldots, a_n) \vdash n \atop 2 \mid (\sum_{j=2 \atop 2 \mid j}^n a_j)} 
|C^{S_n}(\pi_\alpha)| \prod_{i=1}^n Z_i^{a_i} \\
=& \;
\sum_{\alpha = (a_1,\ldots, a_n) \vdash n \atop 2 \mid (\sum_{j=2 \atop 2 \mid j}^n a_j)} 
\frac{2}{|C_{S_n}(\pi_\alpha)|} 
\prod_{i=1}^n Z_i^{a_i}.
\end{align*}
#endif
{
	if (n < 1)
		return error("POLY::cycle_ind_An() n < 1");
	freeself();
	if (n == 1)
		return m_iindex(0);
	init(POLYNOM);
	{
		PARTITION_OB type_v, type_e;
		POLYNOM_OB p;
		INTEGER_OB z, z1;
		INT i, j, k;
		
		type_v.first(n);
		do {
			type_v.t_VECTOR_EXPONENT(&type_e);
			j = 0;
			for (i = 1; i < n; i += 2)
				j += type_e.s_ii(i);
			if (ODD(j))
				continue;
			p.m_skn((VECTOR_OP) callocobject("POLY::cycle_ind_An"), 
				callocobject("POLY::cycle_ind_An"), NIL);
			type_e.s_s()->copy(p.s_s());

			/* the coefficient:
			 * we have class_len in Sn == 
			 * class_len in An !
			 * class_len / (group_order Sn / 2) = 
			 * 2 / centralizer order. */
			
			/* compute centralizer order into z: */
			z.m_i(1);
			for (i = 1; i <= n; i++) {
				j = p.s_s_iii(i - 1);
				/* i = cycle length, 
				 * j = multiplicity. */
				p.s_s_i(i - 1)->fakul(&z1);
				z1.mult_apply(&z);
				z1.m_i(i);
				for (k = 0; k < j; k++)
					z1.mult_apply(&z);
				}
			z1.m_i(2);
			((BRUCH_OP) p.s_k())->m_ou(&z1, &z);
			((BRUCH_OP) p.s_k())->kuerzen();
			p.add_apply(this);
			} while (type_v.next_VECTOR(&type_v));
	}
	return OK;
}

#if TEXDOCU
INT POLYNOM_OB::cycle_index_direct_product(POLYNOM_OP p, POLYNOM_OP q)
#else
Given cycle indices $p$ and $q$ for group actions $_GX$ and $_HY$, 
this function computes the cycle index for the action of $G \times H$ 
on $X \times Y$.
The result is in the this object. 
#endif
{
	POLYNOM_OP pp, qq;
	POLYNOM_OB rr;
	VECTOR_OP vp, vq;
	VECTOR_OB v;
	INT lp, lq, l, i, j, ai, aj, c, d, g, k;
	SYM_OB a;
	
	init(POLYNOM);
	pp = p;
	while (pp != NIL){
		qq = q;
		while (qq != NIL){
			qq->s_k()->mult(pp->s_k(), &a);
			vp = pp->s_s();
			vq = qq->s_s();
			lp = vp->s_li();
			lq = vq->s_li();
			l = lp * lq;
			v.m_il_n(l);
			for (i=0; i<lp ; i++){
				ai = vp->s_ii(i);
				if (ai == 0) continue;
				for (j = 0; j < lq; j++){
					aj = vq->s_ii(j);
					if (aj == 0) continue;
					ggt_iipi(i+1, j+1, &g);
					k = (i+1) * (j+1) / g;
					c = ai * aj * g;
					d = v.s_ii(k-1);
					d += c;
					v.m_ii(k-1, d);
				}
			}
			rr.m_skn(&v, &a, NIL);
			rr.add_apply(this);
			qq = qq->s_n();

		}
		
		pp = pp->s_n();
	}
	return OK;
}

#if TEXDOCU
INT POLYNOM_OB::cycle_index_direct_sum(POLYNOM_OP p, POLYNOM_OP q)
#else
Given cycle indices $p$ and $q$ for group actions $_GX$ and $_HY$, 
this function computes the cycle index for the action of $G \times H$ 
on $X \cup Y$.
The result is in the this object. 
#endif
{
	POLYNOM_OP pp, qq;
	POLYNOM_OB rr;
	VECTOR_OP vp, vq;
	VECTOR_OB v;
	INT lp, lq, l, i, ai, c;
	SYM_OB a;
	
	init(POLYNOM);
	pp = p;
	while (pp != NIL){
		qq = q;
		while (qq != NIL){
			qq->s_k()->mult(pp->s_k(), &a);
			vp = pp->s_s();
			vq = qq->s_s();
			lp = vp->s_li();
			lq = vq->s_li();
			l = MAXIMUM(lp, lq);
			v.m_il_n(l);
			for (i = 0; i < l; i++){
				if (i < lp){
					ai = vp->s_ii(i);
					c = v.s_ii(i);
					c += ai;
					v.m_ii(i, c);
				}
				if (i < lq){
					ai = vq->s_ii(i);
					c = v.s_ii(i);
					c += ai;
					v.m_ii(i, c);
				}
			}
			rr.m_skn(&v, &a, NIL);
			rr.add_apply(this);
			qq = qq->s_n();

		}
		pp = pp->s_n();
	}
	return OK;
}

#if TEXDOCU
INT POLYNOM_OB::cycle_index_Sk_Snmk(INT n, INT k)
#else
Computes the cycle index of the group $S_k \times S_{n-k}$. 
This function uses cycle\_ind\_Sn and cycle\_index\_direct\_product.
#endif
{
	POLYNOM_OB p, q;
	
	p.cycle_ind_Sn(k);
	q.cycle_ind_Sn(n-k);
	cycle_index_direct_product(&p, &q);
	return OK;
}

#if TEXDOCU
INT number_of_double_cosets_A_Sn_B(INT n, POLYNOM_OP ca, POLYNOM_OP cb, SYM_OP res)
#else
This function computes the number of double cosets $A \backslash S_n / B$ 
with $A, B \le S_n$. 
The cycle indices of $A$ and $B$ are ca and cb. The result is put into res.
The following formula is used, it is an application of the Lemma of Cauchy-Frobenius:
\begin{align*}
|A \backslash S_n / B| 
=& \; 
\frac{n!}{|A| \cdot |B|} 
\sum_{\alpha \vdash n} 
\frac{|C^{S_n}(\pi_\alpha) \cap A| \cdot 
|C^{S_n}(\pi_\alpha) \cap B|}{|C^{S_n}(\pi_\alpha)|} \\
=& \; 
\frac{1}{|A| \cdot |B|} 
\sum_{\alpha \vdash n} 
|C^{S_n}(\pi_\alpha) \cap A| \cdot 
|C^{S_n}(\pi_\alpha) \cap B| \cdot |C_{S_n}(\pi_\alpha)| \\
\end{align*}
#endif
{
	PARTITION_OB part, part2;
	POLYNOM_OP p, q;
	SYM_OB a, b, c, d;
	
	res->m_i_i(0);
	part2.first(n);
	do 
	{
		part2.swap(&part);
		part.t_VECTOR_EXPONENT(&part2);
		p = ca;
		while (p != NIL) {
			if (part2.s_s()->sym_comp(p->s_s()) == 0) {
				break;
			}
			p = p->s_n();
		}
		q = cb;
		while (q != NIL) {
			if (part2.s_s()->sym_comp(q->s_s()) == 0) {
				break;
			}
			q = q->s_n();
		}
		if (p && q) {
			p->s_k()->mult(q->s_k(), &a);
			part2.centralizer_order_Sn(&b);
			a.mult(&b, &c);
			res->add(&c, &d);
			d.swap(res);
			}		
	}
	while (part.next_VECTOR(&part2));
	return OK;
}		

#if TEXDOCU
INT POLYNOM_OB::cycle_index_arbitrary(VECTOR_OP gen)
#else
Computes the cycle index of the group generated by the elements in gen.
Uses labelled branchings.
#endif
{
	LABRA_OB L;
	SYM_OB go;
	
	reduce_generators_labra(gen, &go, FALSE /* f_verbose */, &L);
	cycle_index_labra(&L);
	return OK;
}


#if TEXDOCU
INT POLYNOM_OB::cycle_index_labra(LABRA_OP L)
#else
Computes the cycle index of the group L.
#endif
{
	VECTOR_OB path;
	PERMUTATION_OB p;
	VECTOR_OB type;
	POLYNOM_OB rr;
	SYM_OB a;
	
	L->group_order(&a);
	a.invers_apply();
	init(POLYNOM);
	L->elements_first(&path, &p);
	do {
		p.cycle_type(&type);
		
		// p.print(); printf("type: "); type.println();
		rr.m_skn(&type, &a, NIL);
		rr.add_apply(this);
		
	} while (L->elements_next(&path, &p));
	return OK;
}



// #include "mw_poly.C"
/* mw_poly.C */
/* Martin Wiesend Januar 1996 */

#include <DISCRETA/ma.h>

static void mw_mult_apply_polynom(POLYNOM_OP p, POLYNOM_OP q);
static INT mw_mult_polynom(POLYNOM_OP p, POLYNOM_OP q, POLYNOM_OP r);
static INT PolyaSubst_in_Monom(SYM_OP newKoeff, VECTOR_OP Monom, POLYNOM_OP r);
static INT cycle_type_in_Sm_times_Sn(VECTOR_OP v1, VECTOR_OP v2, VECTOR_OP v3);
static void Kopf(FILE *fp);
static void Fuss(FILE *fp);
static void mw_Latex_Anzahlmatrix(MATRIX_OP M, INT m, INT n);

static void mw_mult_apply_polynom(POLYNOM_OP p, POLYNOM_OP q)
/* p = p * q */
{
   POLYNOM_OB r;
   r.init(POLYNOM);

   mw_mult_polynom(p, q, &r);
   p->freeself();
   p->init(POLYNOM);
   r.add_apply(p);
}
static INT mw_mult_polynom(POLYNOM_OP p, POLYNOM_OP q, POLYNOM_OP r)
{
   r->init(POLYNOM);

   POLYNOM_OB h;
   VECTOR_OB  self;
   INT        koeff = 0, len_p, len_q, len;
   POLYNOM_OP zp, zq;

   len_p = p->s_s()->s_li();
   len_q = q->s_s()->s_li();
   len = (len_p > len_q) ? len_p : len_q; // p habe 3 Var, q habe 5 ==> r hat 5
   h.m_skin(&self, koeff, NULL);
   h.s_s()->m_il_n(len);

   zp = p;
   while(zp != NULL) {
      zq = q;
      while(zq != NULL) {
         /*
         koeff = zp->s_k_ii() * zq->s_k_ii();
         h.c_k_i(koeff);
         for(i=0; i<len_p; i++)
            h.s_s()->c_ii(i,zp->s_s_iii(i));
         for(i=len_p; i<len; i++)
            h.s_s()->c_ii(i,0);
         for(i=0; i<len_q; i++) {
            pot = h.s_s_iii(i) + zq->s_s_iii(i);
            h.s_s()->c_ii(i,pot);
         } // for i
         */
         mult_abc(zp->s_k(), zq->s_k(), h.s_k());
         zp->s_s()->add(zq->s_s(), h.s_s());
         h.add_apply(r);
         zq = zq->s_n();
      } // while zq
      zp = zp->s_n();
   } // while zp

   return OK;
}

INT POLYNOM_OB::sum_of_coefficients(SYM_OP sum)
{
   POLYNOM_OP zeiger = this;

   sum->init(INTEGER);

   while (zeiger != NULL) { 
      zeiger->s_k()->add_apply(sum); 
      zeiger = zeiger->s_n();
   } 

   return OK;
}

static INT PolyaSubst_in_Monom(SYM_OP newKoeff, VECTOR_OP Monom, POLYNOM_OP r)
{
   if(newKoeff->s_obj_k() != BRUCH) {
      return error("*** warning 1 in POLYNOM_OB::PolyaSubst");
   }
   INT           i, j;
   INT           len = Monom->s_li();
   VECTOR_OB     self;
   INT           koeff = 1;
   POLYNOM_OP  h, p;  // p = 1 + y^i
   VECTOR_OB     newSelf;
   POLYNOM_OP  q;
   /* POLYNOM_OP    r; */

   r->init(POLYNOM);

   //*************** q := newKoeff ***************
   q = (POLYNOM_OP) callocobject("PolyaSubst_in_Monom");
   q->init(POLYNOM);
   q->m_skn(&newSelf, newKoeff, NULL);
   q->s_s()->m_il_n(1); 
   q->s_s()->c_ii(0,0); 

   for (i = 0; i < len; i++) {
         for (j = 0; j < Monom->s_ii(i); j++) {
            //************* p := 1 + y^i *************
            p = (POLYNOM_OP) callocobject("PolyaSubst_in_Monom");
            p->init(POLYNOM);
            h = (POLYNOM_OP) callocobject("PolyaSubst_in_Monom");
            h->init(POLYNOM);
            h->m_skin(&self, koeff, NULL);
            h->s_s()->m_il_n(1);   // einzige Variable y
            h->s_s()->c_ii(0,0);
            h->add_apply(p);
            h->s_s()->c_ii(0, i + 1); // wg. Indizierung i+1
            h->add_apply(p);
            freeall(h);

						mw_mult_apply_polynom(q, p);
            /*   *q *= *p;   */
            freeall(p);
         } // for j
   } // for i

   q->copy(r);
   freeall(q);
   return OK;
}

INT POLYNOM_OB::PolyaSubst(POLYNOM_OP erzFkt)
{
   erzFkt->init(POLYNOM);

   POLYNOM_OB  q;
   POLYNOM_OP  zeiger = this;  // zum Durchlaufen des Zykelindex

   while (zeiger != NULL) {   // das naechste Monom im Zykelindex
      //*************** q = 1/12 * (1 + y^6)^3 * (1 + y^17)
      PolyaSubst_in_Monom(zeiger->s_k(), zeiger->s_s(), &q);
      q.add_apply(erzFkt);           // erzFkt += q
      q.freeself();

      zeiger = zeiger->s_n();
   } // while

   return OK;
}

static INT cycle_type_in_Sm_times_Sn(VECTOR_OP v1, VECTOR_OP v2, VECTOR_OP v3)
{
   INT i, j, a_i, b_j, I, J;
   INT l1 = v1->s_li();
   INT l2 = v2->s_li();
   INT ggT, kgV, ZwSumme;

   for (i = 0; i < l1; i++) {
      a_i = v1->s_ii(i);
      if (a_i > 0) {
         I = i + 1; 
         for (j = 0; j < l2; j++) {
            b_j =  v2->s_ii(j);
            if (b_j > 0) {
               J = j + 1; 
               ggt_iipi(I, J, &ggT);
               kgV = (I * J / ggT) - 1;
               ZwSumme = v3->s_ii(kgV) + ggT * a_i * b_j;
               v3->c_ii(kgV, ZwSumme); 
            } // if b_j
         } // for j
      } // if a_i
   } // for i

   return OK;
}

#include <DISCRETA/lo.h>

INT POLYNOM_OB::cycle_ind_Sm_times_Sn(INT m, INT n)
{
   init(POLYNOM);
#ifdef BRUCHTRUE
   {
      POLYNOM_OB cycSm, cycSn;
      cycSm.cycle_ind_Sn(m);
      cycSn.cycle_ind_Sn(n);

      BRUCH_OB newKoeff;
      LONGINT_OB Eins, Nenner, N1, N2;
      Eins.m_i(1);
      VECTOR_OB newSelf;
      POLYNOM_OB p;
      POLYNOM_OP zeiger1, zeiger2;  // zum Durchlaufen der Zykelindices

      //++++++++++++++ durchlaufe alle Paare von ZykelTYPEN
      // POLYNOM ==> LISTTRUE ist definiert, s. also LIST_OP->sprint
      zeiger1 = &cycSm;
      while (zeiger1 != NULL) {
      zeiger2 = &cycSn;
      while (zeiger2 != NULL) {

      //++++++++++++++

         //*********** Koeffizient (ein Bruch) des Monoms p
         // Nenner := Nenner(Koeff(i.Monom(cycSm)))
         //           * Nenner(Koeff(i.Monom(cycSm)))
         if (zeiger1->s_k()) {
            if(m == 1)
               N1.m_i(1);  // der Koeff von zeiger1 ist kein Bruch!
            else
               N1.m_i(((BRUCH_OP)zeiger1->s_k())->s_ui()); 
         }
         else 
            return error("*** warning 1 in POLYNOM_OB::cycle_ind_Sm_times_Sn");
         if (zeiger2->s_k()) {
            if(n == 1)
               N2.m_i(1);
            else
               N2.m_i(((BRUCH_OP)zeiger2->s_k())->s_ui()); 
         }
         else
            return error("*** warning 2 in POLYNOM_OB::cycle_ind_Sm_times_Sn");
         // der Koeffizient des durch cycle_ind_Sn kreierten Polynoms 
         // hat ob_kind 4, ist also ein BRUCH, s. stypes.h
         Nenner.m_i(1);
         N1.mult_apply_longint(&Nenner);
         N2.mult_apply_longint(&Nenner);
         newKoeff.m_ou(&Eins, &Nenner);

         //*********** der Vektor der Potenzen des Monoms p
         if (!zeiger1->s_s())
            return error("*** warning 3 in POLYNOM_OB::cycle_ind_Sm_times_Sn");
         if (!zeiger2->s_s())
            return error("*** warning 4 in POLYNOM_OB::cycle_ind_Sm_times_Sn");
         p.m_skn(&newSelf, &newKoeff, NIL); // Koeff = der Bruch newKoeff
         p.s_s()->m_il_n(m*n);              // Vektor der Laenge m*n, init. mit 0en
         cycle_type_in_Sm_times_Sn(zeiger1->s_s(), zeiger2->s_s(), p.s_s());
                                             // den Vektor fuellen

         //*********** das Monom p zum Polynom this dazuaddieren 
         p.add_apply(this); 
      //++++++++++++++
      zeiger2 = zeiger2->s_n();
      }
      zeiger1 = zeiger1->s_n();
      }
      //++++++++++++++ Ende while zeiger...
   }
#else
   return error("cycle_ind_Sm_times_Sn(): BRUCH not available");
#endif

   return OK;  
}

//*************************************************************************
// Sei c == 'x'.
// indiziere == TRUE  ==> Variablen heissen x_1, x_2, x_3, ...
//           == FALSE ==>                   x, y, z; 
// if Koeff-objectkind = BRUCH_OB oder INT oder LONGINT oder ...
// Das ausgegebene Polynom formatieren (Zeilenumbrueche, linker Rand, ...)
// muss der Benutzer selbst nach seinen individuellen Vorstellungen!
//
// Achtung: es wird nicht kontrolliert, ob negative Koeffizienten dabei sind.
//*************************************************************************

INT POLYNOM_OB::LaTeX(FILE *outLatex, char c, INT indiziere)
{
   if((c < 'A') || (c > 'Z') && (c < 'a') || (c > 'z')) {
      return error("*** warning 1 in POLYNOM_OB::LaTeX");
   }
 
   POLYNOM_OP zeiger = this;  // zum Durchlaufen des Zykelindex
   SYM_OP Koeff, Zaehler, Nenner;
   VECTOR_OP Monom;
   INT len, pot;
   char i;
   char str[2], str2[64], str3[64];

   fprintf(outLatex, "   ");  // Zeilen umbrechen mit \\&&
   while (zeiger != NULL) {
      Koeff = zeiger->s_k();
      switch(Koeff->s_obj_k()) {
         case BRUCH:
            Zaehler = ((BRUCH_OP) Koeff)->s_o();
            Nenner  = ((BRUCH_OP) Koeff)->s_u();
            str2[0] = '\0';
            Zaehler->sprint(str2);
            str3[0] = '\0';
            Nenner->sprint(str3);
            fprintf(outLatex, "\\frac{%s}{%s}", str2, str3);
           break;
         case INTEGER:
            fprintf(outLatex, "%ld", ((INTEGER_OP)Koeff)->s_i());
            break;
         case LONGINT:
            str2[0] = '\0';
            Koeff->sprint(str2);
            fprintf(outLatex, "%s", str2);
            break;
         default:
            return error("*** warning 3 in POLYNOM_OB::LaTeX");
      };  // switch Koeff 

      Monom = zeiger->s_s();
      len = Monom->s_li();
      if (indiziere == TRUE) {
         str[0] = c; str[1] = '\0';
         for(i = 0; i < len; i++) {
            pot = Monom->s_ii(i);
            if (pot > 0) {
            	fprintf(outLatex, " %s_{%ld}", str, (INT) (i + 1)); 
               if(pot >= 2) 
               	fprintf(outLatex, "^{%ld}", (INT) pot);
            }
         }
      } 
      else // indiziere == FALSE
         if((c + len > 'z') || (c + len > 'Z') && (c + len < 'a')) {
            return error("*** warning 2 in POLYNOM_OB::LaTeX");
         }
         else 
            for (i = 0; i < len; i++) {
               pot = Monom->s_ii(i);
               if (pot > 0) {
                  str[0] = c + i;
                  str[1] = '\0';
            			fprintf(outLatex, " %s", str); 
                 if(pot >= 2)
            				fprintf(outLatex, "^{%ld}", (INT) pot); 
               }
            }
      zeiger = zeiger->s_n();
      if (zeiger != NULL)
         fprintf(outLatex, "\n + "); 
      else
         fprintf(outLatex, ".\n"); 
   } // while

  return OK;
}

//*** eine m x n - Matrix mit LONGINT-Werten wird im Latex-Format in
//*** das File mwAnz.tex geschrieben

static void mw_Latex_Anzahlmatrix(MATRIX_OP M, INT m, INT n)
{
  BYTE str[64];
  INT i, j;
  FILE *fp;
  
  fp = fopen("mwAnz.tex", "w");
  Kopf(fp);
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabulator}{c||");
  for(i=0; i < n; i++)
     fprintf(fp, "r");
  fprintf(fp, "}\n");
  fprintf(fp, "$|G|\\backslash |M|$");
  for(i=0; i < n; i++)
     fprintf(fp, "&%ld", i+1);
  fprintf(fp, "\\\n");
  fprintf(fp, "\\hline\\hline\n");

  for(i=0; i < m; i++) {
     fprintf(fp, "%ld & ", i + 1);
     for(j=0; j < n; j++) {
        str[0] = 0;
        M->s_ij(i, j)->sprint(str);
        fprintf(fp, "%s", str);
        if(j < n-1)
        	fprintf(fp, " & ");
     }
     if(i < m-1)
     	fprintf(fp, "\\\\");
    	fprintf(fp, "\n");
  } 

  fprintf(fp, "\\end{tabulator}\n");
  fprintf(fp, "\\end{center}\n");
  Fuss(fp);
  fclose(fp);
}

static void Kopf(FILE *fp)
{
  fprintf(fp, "\\documentstyle[a4,german]{report}\n");
  fprintf(fp, "\\begin{document}\n\n");
}

static void Fuss(FILE *fp)
{
  fprintf(fp, "\\end{document}\n\n");
}

#include <DISCRETA/divs.h> /* for STRING_OB */

void mw_Latex_ZykInd_und_erzFkt(POLYNOM_OP ZykInd, POLYNOM_OP erzFkt, INT m, INT n)
{
  BYTE *cstr = new BYTE[3]; 
  STRING_OB str1, str2;
  str1.init("mw_");
  cstr[0] = 0;
  sprintf(cstr, "%ld_%ld", m, n);
  str2.init(cstr);
  str1.append(&str2); 
  str2.init(".tex");
  str1.append(&str2);
	FILE *fp;
	
  fp = fopen(str1.s_str(), "w");

  Kopf(fp);
  fprintf(fp, "\\begin{eqnarray*}\n");
  fprintf(fp, "\\lefteqn{ C(S_{\\underline{%ld}}\\times S_{\\underline{%ld}},", m, n);
  fprintf(fp, "\\underline{%ld}\\times\\underline{%ld}) = }\\\\", m, n);
  ZykInd->LaTeX(fp, 'z', TRUE);
  fprintf(fp, "\\end{eqnarray*}\n\n");

  fprintf(fp, "\\begin{eqnarray*}\n");
  fprintf(fp, "\\lefteqn{ C(S_{\\underline{%ld}}\\times S_{\\underline{%ld}},", m, n);
  fprintf(fp, "\\underline{%ld}\\times\\underline{%ld}\\mid 1+y) = }\\\\", m, n);
  erzFkt->LaTeX(fp, 'y', FALSE);
  fprintf(fp, "\\end{eqnarray*}\n\n");
  Fuss(fp);

  fclose(fp);
}

void mw_Zykelindizes_und_Anzahlen(INT m, INT n)
{
	MATRIX_OB M;
   INT Max, i, j; 
   if (n > m) { Max = n; n = m; m = Max; } // ==>   m >= n

   POLYNOM_OP ZykInd, erzFkt;
   M.m_ilih(n, m);

   for(j = 1; j <= n; j++)
   for(i = j; i <= m; i++)
   {
      printf("i = %ld  j = %ld\n", i, j);
      fflush(stdout);
      erzFkt = (POLYNOM_OP) callocobject("mw_Zykelindizes_und_Anzahlen");
      erzFkt->init(POLYNOM);
      ZykInd = (POLYNOM_OP) callocobject("mw_Zykelindizes_und_Anzahlen");
      ZykInd->init(POLYNOM);
      ZykInd->cycle_ind_Sm_times_Sn(i,j);
      
      ZykInd->PolyaSubst(erzFkt);
      
      mw_Latex_ZykInd_und_erzFkt(ZykInd, erzFkt, i, j);
      erzFkt->sum_of_coefficients(M.s_ij(i-1, j-1));
      if (i > j)
      	M.s_ij(i-1, j-1)->copy(M.s_ij(j-1,i-1));
      freeall(ZykInd);
      freeall(erzFkt);
   }
   mw_Latex_Anzahlmatrix(&M, m, n);
}


// end of mw_poly.C 


#if TEXDOCU
INT POLYNOM_OB::nabla(POLYNOM_OP res, LABRA_OP G, INT n, INT f_v)
#else
#endif
{
	POLYNOM_OP p;
	POLYNOM_OB q, q1, q2, q3;
	VECTOR_OB V;
	INT i = 0;

	p = this;
	while (p != NIL) {
		q.freeself();
		if (f_v) {
			printf("p=");
			p->s_mo()->println();
			printf("\n");
			}
		p->s_s()->copy(&V);
		V.println();
		V.Nabla(G, &q, n, f_v);
		if (f_v) {
			printf("Nabla=");
			q.Print();
			printf("\n");
			}
		q1.freeself();
		q.mult(p->s_k(), &q1);
		if (i == 0)
			q1.swap(&q2);
		else {
			q3.freeself();
			q1.add(&q2, &q3);
			q3.swap(&q2);
			}
		p = p->s_n();
		i++;
		}
	q2.swap(res);
	return OK;
}



#endif /* POLYTRUE */

