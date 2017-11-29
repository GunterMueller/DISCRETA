/* bruch.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef BRUCHTRUE

#include <DISCRETA/bruch.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/poly.h>
#include <DISCRETA/vec.h>

static struct bruch * callocbruch();

INT bruch_anfang()
{
	return OK;
}

INT bruch_ende()
{
	INT erg = OK;
	return erg;
}

static struct bruch * callocbruch()
{
	struct bruch * ergebnis;
	
	ergebnis = (struct bruch *) my_malloc(sizeof(struct bruch), "callocbruch");
	if (ergebnis == NULL) {
		error("callocbruch(): no memory");
		return NIL;
		}
	return ergebnis;
}

INT BRUCH_OB::m_ioiu(INT oben, INT unten)
/* m_ioiu_b(oben,unten,ergebnis) */
/* ein bruch mit einem integer eintrag im zaehler und einem
 * integer eintrag im nenner z.b. oben = 3 unten = 5 --> 3/5 */
{
	INT erg = OK;
	
	erg += b_ou(callocobject("BRUCH"), callocobject("BRUCH"));
	((INTEGER_OP) s_o())->m_i(oben);
	((INTEGER_OP) s_u())->m_i(unten);
	return erg;
}

INT BRUCH_OB::m_ou(SYM_OP oben, SYM_OP unten)
/* m_ou_b(oben, unten, ergebnis) */
{
	INT erg = OK;
	
	erg += b_ou(callocobject("BRUCH"), callocobject("BRUCH"));
	erg += oben->copy(s_o());
	erg += unten->copy(s_u());
	return erg;
}

INT BRUCH_OB::b_ou(SYM_OP oben, SYM_OP unten)
/* b_ou_b(oben,unten,ergebnis) */
{
	INT erg = OK;
	OBJECTSELF d;

	d.ob_bruch = callocbruch();
	b_ks(BRUCH, d);
	c_o(oben);
	c_u(unten);
	c_i(NGEKUERZT);
	return erg;
}

INT BRUCH_OB::m_scalar(SYM_OP a)
/* m_scalar_bruch(a, b) */
/* macht aus scalar bruch */
/* die integerzahl 5 wird z.B. 5/1 */
{
	SYM_OB a1;
	
	if (a->s_obj_k() == BRUCH)
		return error("BRUCH_OB::m_scalar(): a is bruch");
	b_ou(callocobject("BRUCH"), callocobject("BRUCH"));
	a->copy(s_o());
	a->copy(s_u());
	s_u()->one();
	return OK;
}

INT BRUCH_OB::freeself()
{
	OBJECTSELF d;
	
	d = s_obj_s();
	freeall(s_o());
	freeall(s_u());
	my_free(d.ob_bruch);
	c_obj_k(EMPTY);
	return OK;
}

INT BRUCH_OB::copy(BRUCH_OP nach)
{
	INT erg = OK;
	
	if (this == nach) 
		return error("BRUCH::copy() this == nach");
	erg += nach->m_ou(s_o(), s_u());
	return erg;
}

INT BRUCH_OB::scalar_it()
{
	SYM_OB a;
	
	if (s_o()->nullp() || s_u()->einsp()) {
		s_o()->copy(&a);
		swap(&a);
		return OK;
		}
	return OK;
}

INT BRUCH_OB::kuerzen()
{
	
#if 0
	if (kuerzen_yn == 1L) 
		return OK; /* d.h. nicht kuerzen */
#endif

	if (s_obj_k() != BRUCH)
		return(OK);

	if (s_i() == GEKUERZT)
		return OK;

	if (s_o()->nullp()) {
		SYM_OB a;
		
		s_o()->copy(&a);
		swap(&a);
		return OK;
		}
	if (s_u()->nullp()) {
		return error("BRUCH:: s_u() is 0");
		}

	if ((s_o()->s_obj_k() == INTEGER) &&
	    (s_u()->s_obj_k() == INTEGER) ) {
		INT oi = s_oi();
		INT ui = s_ui();
		INT g;
		
		if (ui == 0)
			return error("BRUCH:: s_u() is 0 (integer)");
		ggt_iipi(oi, ui, &g);
		oi /= g;
		ui /= g;
		if (ui < 0) {
			ui *= -1;
			oi *= -1;
			}
		((INTEGER_OP)s_o())->m_i(oi);
		((INTEGER_OP)s_u())->m_i(ui);
		c_i(GEKUERZT);
		scalar_it();
		return(OK);
	}
	{
	SYM_OB u, v, g, tmp1, q, r;
	
#if 0
	printf("BRUCH_OB:: kuerzen() calling bezout()\n");
	fflush(stdout);
#endif
	bezout(s_o(), s_u(), &u, &v, &g);
#if 0
	printf("g, o, u=\n");
	g.println();
	fflush(stdout);
#endif
	
	/* kuerzen: */
	s_o()->quores(&g, &u, &v);
	u.swap(s_o());
	s_u()->quores(&g, &u, &v);
	u.swap(s_u());
#if 0
	s_o()->println();
	s_u()->println();
	fflush(stdout);
#endif
	
	/* Nenner ist Einheit ? */
	tmp1.init(s_u()->s_obj_k());
	tmp1.one();
#if 0
	printf("BRUCH_OB:: kuerzen() calling quores()\n");
	fflush(stdout);
#endif
	tmp1.quores(s_u(), &q, &r);
	if (r.nullp()) { /* s_u() is a unit and q it's inverse */
		s_u()->mult(&q, &u);
		u.swap(s_u());
		s_o()->mult(&q, &v);
		v.swap(s_o());
#if 0
		printf("o, u=\n");
		s_o()->println();
		s_u()->println();
		fflush(stdout);
#endif
		}
	c_i(GEKUERZT);
	scalar_it();
	return(OK);
	}
#if 0
	SYM_OB a;
	OP ggterg, moderg;
	INT erg=OK;
	INT ggtierg,vorzeichen=1L;
	moderg = callocobject("BRUCH");
	if (negp(S_B_U(bruch))) {
		vorzeichen *= -1L;
		erg += addinvers_apply(S_B_U(bruch));
	}
	if (negp(S_B_O(bruch))) {
		vorzeichen *= -1L;
		erg += addinvers_apply(S_B_O(bruch));
	}


	ggterg = callocobject("BRUCH");

	erg += ggt(S_B_O(bruch),S_B_U(bruch),ggterg);

	erg += ganzdiv(S_B_O(bruch),ggterg,S_B_O(bruch));
	erg += ganzdiv(S_B_U(bruch),ggterg,S_B_U(bruch));
	erg += freeall(ggterg);

	if (einsp(S_B_U(bruch)) )
	{
		erg += copy(S_B_O(bruch),moderg);
		erg += freeself(bruch);
		*bruch = *moderg;
		if (vorzeichen == -1L)
			erg += addinvers_apply(bruch);
		C_O_K(moderg,EMPTY);
		freeall(moderg);
		if (S_O_K(bruch) == BRUCH)
			C_B_I(bruch,GEKUERZT);
		return OK;
	}

	erg += freeall(moderg);
	if (vorzeichen == -1L)
		erg += addinvers_apply(bruch);
	if (erg != OK)
		return EDC("kuerzen");
	if (S_O_K(bruch) == BRUCH)
		C_B_I(bruch,GEKUERZT);

	return erg;
#endif
}

INT BRUCH_OB::add(SYM_OP b, SYM_OP c)
/* add_bruch(a,b,c) */
{
	INT erg = OK;
	
	if (s_obj_k() != BRUCH) {
		printf("BRUCH::add: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->add(b, c);
		}
	switch (b->s_obj_k()) {
	case INTEGER:
	case LONGINT: 
		erg += add_scalar(b, c);
		break;
	case BRUCH: 
		erg += add_bruch((BRUCH_OP) b, c); 
		break;
#ifdef POLYTRUE
	case POLYNOM:
		erg += ((POLYNOM_OP) b)->add_scalar(this, (POLYNOM_OP) c);
		break;
#endif /* POLYTRUE */
	default :
		return error("BRUCH::add wrong type");
	}
	return erg;
}

INT BRUCH_OB::add_scalar(SYM_OP b, SYM_OP c)
{
	INT erg = OK;
	BRUCH_OB d;
	
	if (s_obj_k() != BRUCH) {
		printf("BRUCH::add_scalar: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->add(b, c);
		}
	d.m_scalar(b);
	erg += add_bruch(&d, c);
	return erg;
}

INT BRUCH_OB::add_bruch(BRUCH_OP b, SYM_OP c)
/* add_bruch_bruch(a,b,c) */
{
	SYM_OB a1, b1, d, e;
	BRUCH_OB f;
	INT erg = OK;

	if (s_obj_k() != BRUCH || b->s_obj_k() != BRUCH) {
		printf("BRUCH::add_bruch: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->add(b, c);
		}
	erg += s_u()->mult(b->s_u(), &d);
	erg += s_o()->mult(b->s_u(), &a1); 
	erg += s_u()->mult(b->s_o(), &b1); 
	erg += a1.add(&b1, &e);
	f.m_ou(&e, &d);
	f.kuerzen();
	f.swap(c);
	return erg;
}

INT BRUCH_OB::addinvers(BRUCH_OP b)
{
	INT erg = OK;
	
	if (s_obj_k() != BRUCH) {
		printf("BRUCH::addinvers: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->addinvers(b);
		}
	erg += b->b_ou(callocobject("BRUCH"), callocobject("BRUCH"));
	erg += s_o()->addinvers(b->s_o());
	erg += s_u()->copy(b->s_u());
	return erg;
}

INT BRUCH_OB::addinvers_apply()
{ 
	if (s_obj_k() != BRUCH) {
		printf("BRUCH::addinvers_apply: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->addinvers_apply();
		}
	return(s_o()->addinvers_apply()); 
}

INT BRUCH_OB::add_apply(SYM_OP b)
/* add_apply_bruch(a,b) */
{
	INT erg = OK;

	if (s_obj_k() != BRUCH) {
		printf("BRUCH::add_apply: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->add_apply(b);
		}
	switch (b->s_obj_k()) {
	case BRUCH:
		erg += add_apply_bruch((BRUCH_OP) b); 
		break;
	default:
		{
			SYM_OB c;
			c = *b;
			b->c_obj_k(EMPTY);
			erg += add(&c, b);
			break;
		}
	}
	return erg;
}

INT BRUCH_OB::add_apply_bruch(BRUCH_OP b)
/* add_apply_bruch_bruch(a,b) */
/* b = b + a */
{
	INT erg = OK;
	BRUCH_OB c;

	if (s_obj_k() != BRUCH || b->s_obj_k() != BRUCH) {
		printf("BRUCH::add_apply_bruch: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->add_apply(b);
		}
	c = *b;
	b->c_obj_k(EMPTY);
	erg += add_bruch(&c, b); /* kuerzt */
	return erg;
}

INT BRUCH_OB::invers(BRUCH_OP b)
/* invers_bruch(a,b) */
{
	INT erg = OK;
	
	if (s_obj_k() != BRUCH) {
		printf("BRUCH::invers: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->invers(b);
		}
	erg += b->b_ou(callocobject("BRUCH"), callocobject("BRUCH"));
	if (s_o()->nullp())
		return error("BRUCH::invers: s_o() is 0");
	s_u()->copy(b->s_o());
	s_o()->copy(b->s_u());
	return erg;
}

INT BRUCH_OB::mult(SYM_OP b, SYM_OP c)
/* mult_bruch(a,b,c) */
{
	INT erg = OK;
	
	if (s_obj_k() != BRUCH) {
		printf("BRUCH::mult: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->mult(b, c);
		}
	switch (b->s_obj_k()) {
	case BRUCH:  
		erg += mult_bruch((BRUCH_OP) b, (BRUCH_OP) c);
		break;
#ifdef INTEGERTRUE
	case LONGINT:
	case INTEGER: 
		erg += mult_integer(b, (BRUCH_OP) c);
		break;
#endif /* INTEGERTRUE */
#ifdef MATRIXTRUE
	case MATRIX: 
		erg += ((MATRIX_OP)b)->mult_scalar(this, (MATRIX_OP) c);
		break;
#endif /* MATRIXTRUE */
#ifdef MONOMTRUE
	case MONOM: 
		erg += ((MONOM_OP)b)->mult_scalar(this, (MONOM_OP) c);
		break;
#endif /* MONOMTRUE */
#ifdef POLYTRUE
	case GRAL:
	case POLYNOM: 
		erg += ((POLYNOM_OP)b)->mult_scalar(this, (POLYNOM_OP) c);
		break;
#endif /* POLYTRUE */
#ifdef SCHURTRUE
	case SCHUR: 
		erg += mult_scalar_schur(a,b,c);
		break;
#endif /* SCHURTRUE */
#ifdef CHARTRUE
	case SYMCHAR: 
		erg += mult_scalar_symchar(a,b,c);
		break;
#endif /* CHARTRUE */
#ifdef VECTORTRUE
	case VECTOR:  
		erg += ((VECTOR_OP)b)->mult_scalar(this, (VECTOR_OP) c);
		break;
#endif /* VECTORTRUE */
	default:
		b->printobjectkind();
		error("BRUCH::mult: wrong second type");
		return(ERROR);
	}
	if (erg != OK)
		return error("BRUCH::mult: error in computing");
	return erg;
}

INT BRUCH_OB::mult_integer(SYM_OP b, BRUCH_OP c)
/* mult_bruch_integer(a,b,c) */
{
	INT erg = OK;
	/* AK fuer integer und longint */

	if (s_obj_k() != BRUCH) {
		printf("BRUCH::mult_integer: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->mult(b, c);
		}
	erg += copy(c);
	erg += b->mult(s_o(), c->s_o());
	erg += c->kuerzen();
	return erg;
}

INT BRUCH_OB::mult_bruch(BRUCH_OP b, BRUCH_OP c)
{
	INT erg = OK;

	if (s_obj_k() != BRUCH || b->s_obj_k() != BRUCH) {
		printf("BRUCH::mult_bruch_bruch: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->mult(b, c);
		}
	erg += c->b_ou(callocobject("BRUCH"), callocobject("BRUCH"));
	erg += s_o()->mult(b->s_o(), c->s_o());
	erg += s_u()->mult(b->s_u(), c->s_u());
	erg += c->kuerzen();
	return erg;
}

INT BRUCH_OB::mult_apply(SYM_OP b)
/* mult_apply_bruch(a,b) */
/* a is BRUCHobject */
{
	INT erg = OK;
	
	if (s_obj_k() != BRUCH) {
		printf("BRUCH::mult_apply: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->mult_apply(b);
		}
	switch (b->s_obj_k()) {
	case BRUCH: 
		erg += s_o()->mult_apply(((BRUCH_OP)b)->s_o());
		erg += s_u()->mult_apply(((BRUCH_OP)b)->s_u());
		((BRUCH_OP)b)->c_i(NGEKUERZT);
		erg += ((BRUCH_OP)b)->kuerzen();
		break;
	default:
		{
		SYM_OB c;
		erg += mult(b, &c);
		c.swap(b);
		}
	}
	return erg;
}

INT BRUCH_OB::compare(SYM_OP b)
/* comp_bruch(a,b) */
{
	BRUCH_OB b1;
	
	if (b->s_obj_k() == BRUCH) {
		/* a/b < c/d   <==>  ad < cb */
		INT erg;
		SYM_OB c, d;

		s_o()->mult(((BRUCH_OP)b)->s_u(), &c);
		((BRUCH_OP)b)->s_o()->mult(s_u(), &d);
		erg = c.sym_comp(&d);
#if 0
		if 	(
		    (negp(S_B_U(a)) && negp(S_B_U(b)))
		    || 
		    (posp(S_B_U(a)) && posp(S_B_U(b)))
		    )
			erg = comp(c,d);
		else 
			erg = comp(d,c);
#endif
		return(erg);
	}
	else if (b->s_obj_k() == INTEGER) {
		b1.m_scalar(b);
		return compare(&b1);
		}
	else	{
		b->printobjectkind();
		return error("comp_bruch: wrong second type");
	}
}

INT BRUCH_OB::einsp()
{
	if (s_obj_k() != BRUCH) {
		printf("BRUCH::einsp: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->einsp();
		}
	kuerzen();
	return s_o()->einsp();
}

INT BRUCH_OB::negeinsp()
{
	if (s_obj_k() != BRUCH) {
		printf("BRUCH::negeinsp: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->negeinsp();
		}
	kuerzen();
	return s_o()->negeinsp();
}

INT BRUCH_OB::nullp()
{
	if (s_obj_k() != BRUCH) {
		printf("BRUCH::nullp: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->nullp();
		}
#if 0
	printf("BRUCH_OB::nullp() kuerzen\n");
	fflush(stdout);
#endif
	kuerzen();
	// now the object may have changed to integer !
	if (s_obj_k() != BRUCH) {
		printf("BRUCH::nullp: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->nullp();
		}
#if 0
	printf("BRUCH_OB::nullp() s_o=\n");
	fflush(stdout);
	s_o()->println();
#endif
	return s_o()->nullp();
}

INT BRUCH_OB::sprint(BYTE *str)
{
	BYTE str1[256], str2[256];
	
	if (s_obj_k() != BRUCH) {
		printf("BRUCH::sprint: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->sprint(str);
		}
	if (strlen(str) > 200)
		return OK;
	str1[0] = 0;
	str2[0] = 0;
	s_o()->sprint(str1);
	s_u()->sprint(str2);
	sprintf(str + strlen(str), "%s/%s", str1, str2);
	return OK;
}

INT BRUCH_OB::latex(FILE *fp)
{
	if (s_obj_k() != BRUCH) {
		printf("BRUCH::latex: warning: not of type BRUCH\n");
		return ((SYM_OP)this)->latex(fp);
		}
	s_o()->latex(fp);
	fprintf(fp, "/");
	s_u()->latex(fp);
	return OK;
}

#endif /* BRUCHTRUE */

