/* mon.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef MONOMTRUE

#include <DISCRETA/poly.h>
#include <DISCRETA/ma.h>

static struct monom * callocmonom(void);

static struct monom * callocmonom()
{
	struct monom *c = (struct monom *) my_malloc(sizeof(struct monom), "callocmonom()");
	if (c == NULL) error("callocmonom:no mem");
	return(c);
}

#if TEXDOCU
INT MONOM_OB::m_sk(VECTOR_OP self, SYM_OP koeff)
#endif
{
	b_sk((VECTOR_OP) callocobject("MONOM::m_sk()"), callocobject("MONOM::m_sk"));
	((SYM_OP)self)->copy(s_s());  
	koeff->copy(s_k());  
	return OK;
}

#if TEXDOCU
INT MONOM_OB::m_ski(VECTOR_OP self, INT koeff)
#endif
{
	b_sk((VECTOR_OP) callocobject("MONOM::m_ski"), callocobject("MONOM::m_ski"));
	((SYM_OP)self)->copy(s_s());  
	m_k_i(koeff);
	return OK;
}

#if TEXDOCU
INT MONOM_OB::m_siki(INT self, INT koeff)
#endif
/* z.B. fuer MONOPOLYs verwendet. */
{
	b_sk((VECTOR_OP) callocobject("MONOM::m_siki"), callocobject("MONOM::m_siki"));
	((INTEGER_OP)s_s())->m_i(self);
	m_k_i(koeff);
	return OK;
}

#if TEXDOCU
INT MONOM_OB::b_sk(VECTOR_OP self, SYM_OP koeff)
#endif
/* build_self koeff_monom */
{
	OBJECTSELF d;

	d.ob_monom = callocmonom();
	
	b_ks(MONOM, d);
	c_s(self);
	c_k(koeff);
	return OK;
}

#if TEXDOCU
INT MONOM_OB::b_ski(VECTOR_OP self, INT koeff)
#endif
{
	OBJECTSELF d;
	INTEGER_OP k;

	d.ob_monom = callocmonom();
	k = (INTEGER_OP) callocobject("MONOM::b_ski");
	
	b_ks(MONOM, d);
	c_s(self);
	c_k(k);
	return OK;
}

#if TEXDOCU
INT MONOM_OB::fprint(FILE *f)
#endif
{
	s_k()->fprint(f); 
	fprintf(f, " "); 
	((SYM_OP)s_s())->fprint(f);
	return OK;
}

#if TEXDOCU
INT MONOM_OB::sprint(BYTE *str)
#endif
/* appends to str. writes to maximal strlength of 200. */
{
	BYTE str1[256];
	
	if (s_k()) {
		str1[0] = 0;
		((SYM_OP)s_k())->sprint(str1);
		}
	else {
		strcpy(str1, "*");
		}
	if (strlen(str) + strlen(str1) < 200) {
		strcat(str, str1);
		strcat(str, " ");
		}
	else
		return OK;
	if (s_s()) {
		str1[0] = 0;
		((SYM_OP)s_s())->sprint(str1);
		}
	else {
		strcpy(str1, "*");
		}
	if (strlen(str) + strlen(str1) < 200)
		strcat(str, str1);
	return OK;
}

#if TEXDOCU
INT MONOM_OB::freeself()
#endif
{
	OBJECTSELF d; 
	INT erg = OK;
	
	d = s_obj_s(); 
	erg += freeall(s_s()); 

	erg += freeall(s_k()); 

	my_free(d.ob_monom);
	c_obj_k(EMPTY);
	return(OK);
}

#if TEXDOCU
INT MONOM_OB::compare(MONOM_OP b)
#endif
/* Bei gleichheit von self wird koeff verglichen */ 
{
	INT erg;
	
	erg = ((SYM_OP)s_s())->sym_comp(b->s_s());
	if (erg != 0)
		return(erg);
	return s_k()->sym_comp(b->s_k());
}

#if TEXDOCU
INT MONOM_OB::copy(MONOM_OP b)
#endif
{
	return b->m_sk(s_s(), s_k());
}

#if TEXDOCU
INT MONOM_OB::mult_scalar(SYM_OP a, MONOM_OP c)
#endif
/* a ist skalar */
/* before: mult_scalar_monom(OP a, OP b, OP c) */
{
	if (! c->emptyp()) 
		c->freeself();
	c->b_sk((VECTOR_OP) callocobject("MONOM::mult_scalar"), callocobject("MONOM::mult_scalar"));
	((SYM_OP)s_s())->copy(c->s_s()); 
	return s_k()->mult(a, c->s_k());
}

#if TEXDOCU
INT MONOM_OB::mult_apply_scalar(SYM_OP a)
#endif
/* before: mult_apply_scalar_monom(OP a, OP b) */
{
	return a->mult_apply(s_k());
}

#if TEXDOCU
INT MONOM_OB::add(SYM_OP b, MONOM_OP c)
#endif
{
	INT erg = OK;
	
	switch (b->s_obj_k()) {
#ifdef HOMSYMTRUE
		case HOM_SYM: erg += add_monom_homsym(a,b,c);break;
#endif /* HOMSYMTRUE */
#ifdef SCHURTRUE
		case SCHUR: erg += add_monom_schur(a,b,c);break;
#endif /* SCHURTRUE */
		default: return(error("MONOM::add:"));
		}
	if (erg != OK) 
		error("MONOM::add: error during computation");
	return erg;
}

#if TEXDOCU
INT MONOM_OB::addinvers(MONOM_OP b)
#endif
{
	b->b_sk((VECTOR_OP)callocobject("MONOM::addinvers"), callocobject("MONOM::addinvers"));
	((SYM_OP)s_s())->copy(b->s_s());
	s_k()->addinvers(b->s_k());
	return(OK);
}

#if TEXDOCU
INT MONOM_OB::addinvers_apply()
#endif
{ 
	return(s_k()->addinvers_apply()); 
}

#if TEXDOCU
INT add_koeff(SYM_OP a, SYM_OP b)
#endif
/* eqhandle bei insert bei polynomen, monopoly.
 * Addiert die beiden koeffizienten der monome a und b
 * nach dem monom b */
{
	MONOM_OP a1, b1;
	
	a1 = (MONOM_OP) a;
	b1 = (MONOM_OP) b;
	if (a->s_obj_k() != MONOM) 
			return error("add_koeff: first != MONOM");
	if (b->s_obj_k() != MONOM) 
			return error("add_koeff: second != MONOM");
	a1->s_k()->add_apply(b1->s_k());
	a1->freeself();
	if (b1->s_k()->nullp()) 
		b1->freeself();
	return(OK);
}

#if TEXDOCU
INT comp_monomvector_monomvector(SYM_OP a, SYM_OP b)
#endif
/* vergleicht zwei monome von vectorform */
/* wenn keine vectorform dann der normale vergleich */
{
	INT i, j, erg;
	VECTOR_OP as, bs;
	MONOM_OP a1, b1;
	
	a1 = (MONOM_OP) a;
	b1 = (MONOM_OP) b;
	if (a->s_obj_k() != MONOM) {
		a->printobjectkind();
		return error("comp_monomvector_monomvector: a != MONOM");
		}
	if (b->s_obj_k() != MONOM) {
		b->printobjectkind();
		return error("comp_monomvector_monomvector: b != MONOM");
		}
	as = a1->s_s(); 
	bs = b1->s_s();

#ifdef SCHUBERTTRUE
	if ((as->s_obj_k() == PERMUTATION) &&
		(bs->s_obj_k() == PERMUTATION)) { /* if schubert polynom */
		erg = 1L;
		if (S_P_LI(bs) > S_P_LI(as)) {
			a=bs;
			bs=as;
			as=a;
			erg= -1L;
		}
		/* as ist laenger als bs */
		for (i=0L; i<S_P_LI(as); i++)
			{
			if (i < S_P_LI(bs))
				{
				if (S_P_II(as,i) > S_P_II(bs,i)) return erg*1L;
				if (S_P_II(as,i) < S_P_II(bs,i)) return erg*-1L;
				}
			else {
				if (S_P_II(as,i) < i+1L) return erg*-1L;
				if (S_P_II(as,i) > i+1L) return erg*1L;
				}
			}
		return 0L;
		}
#endif
#ifdef MATRIXTRUE
	if ((as->s_obj_k() == MATRIX) &&
		(bs->s_obj_k() == MATRIX)) {
		
		MATRIX_OP A, B;
		INT h, l;
		
		A = (MATRIX_OP)as;
		B = (MATRIX_OP)bs;
		h = (A->s_hi() > B->s_hi() ? A->s_hi() : B->s_hi());
		l = (A->s_li() > B->s_li() ? A->s_li() : B->s_li());
		for (i = 0; i < h; i++)
		for (j = 0; j < l; j++) {
			if ((i < A->s_hi()) &&
				(i < B->s_hi()) &&
				(j < A->s_li()) &&
				(j < B->s_li())) 
				if ((erg = A->s_ij(i, j)->sym_comp(B->s_ij(i, j))) != 0) return erg;
			if ((i >= A->s_hi()) &&
				(j < B->s_li()))
				if (! B->s_ij(i, j)->nullp())
					return -1;
			if ((i >= B->s_hi()) &&
				(j < A->s_li()))
				if (! A->s_ij(i, j)->nullp())
					return 1;
			if ((i < A->s_hi()) &&
				(j >= B->s_li()))
				if (! A->s_ij(i, j)->nullp())
					return 1;
			if ((i < B->s_hi()) &&
				(j >= A->s_li()))
				if (! B->s_ij(i, j)->nullp())
					return -1;
			}
		return 0;
		}
#endif
	if ((as->s_obj_k() != VECTOR) || 
		(bs->s_obj_k() != VECTOR))
		return ((SYM_OP)as)->sym_comp((SYM_OP)bs);
	
	for (i = 0; i < as->s_li(); i++) {
		if (as->s_i(i)->s_obj_k() != INTEGER)
			error("comp_monomvector_monomvector: as no INTEGERVECTOR");
		else if (i >= bs->s_li()) {
			if (as->s_ii(i) != 0)
				return 1; 
			}
		else if (bs->s_i(i)->s_obj_k() != INTEGER)
			error("comp_monomvector_monomvector: bs no INTEGERVECTOR");
		else if (as->s_ii(i) > bs->s_ii(i))
			return 1;
		else if (as->s_ii(i) < bs->s_ii(i))
			return -1;
		}

	for (j = i; j < bs->s_li(); j++) {
		if (bs->s_i(j)->s_obj_k() != INTEGER) 
			error("comp_monomvector_monomvector:bs no INTEGERVECTOR");
		else if (bs->s_ii(j) != 0) 
			return -1; 
		}
	return 0;
}

#endif /* MONOMTRUE */


