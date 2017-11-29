/* discreta.C: */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef MATRIXTRUE
#include <DISCRETA/ma.h>
#endif
#ifdef PERMTRUE
#include <DISCRETA/perm.h>
#endif
#ifdef PARTTRUE
#include <DISCRETA/part.h>
#endif
#ifdef LONGINTTRUE
#include <DISCRETA/lo.h>
#endif
#ifdef LISTTRUE
#include <DISCRETA/list.h>
#endif
#ifdef BRUCHTRUE
#include <DISCRETA/bruch.h>
#endif
#ifdef DIVS_TRUE
#include <DISCRETA/divs.h>
#endif
#ifdef POLYTRUE
#include <DISCRETA/poly.h>
#endif
#ifdef UNIPOLYTRUE
#include <DISCRETA/unip.h>
#endif
#ifdef DB_TRUE
#include <DISCRETA/db.h>
#endif
#include <stdlib.h> // for getenv

SYM_OP *speicher;       /* global variable for callocobject/freeall */
INT speicherposition;

INT base_kind[SYM_MAX_KIND];

static INT init_base_kind();
BYTE *discreta_home;
BYTE *discreta_arch;

#if TEXDOCU
INT discreta_init()
#else
initialization of DISCRETA
must be called at the beginning of every 
DISCRETA program.
#endif
{
	BYTE str[1024];

	if (sizeof(INT4) != 4) {
		printf("warning: sizeof(INT4) != 4");
		}
	/* time((time_t *)&l);
	srand((unsigned long)l*l); */
	fflush(stdout); fflush(stderr);
	init_base_kind();
	speicher = (SYM_OP *) my_malloc(SPEICHERSIZE * sizeof(SYM_OP), "discreta_init()");
	if (speicher == NULL)
		error("discreta_init:no mem");
	speicherposition = -1;
	test_swap();
#ifdef LONGINTTRUE
	start_longint();
#endif
#ifdef DB_TRUE
	database_init();
#endif
	discreta_home = getenv("DISCRETA_HOME");
	if (discreta_home == NULL) {
		printf("discreta_init(): WARNING: $DISCRETA_HOME not set !\n");
		fflush(stdout);
		}
	discreta_arch = getenv("DISCRETA_ARCH");
	if (discreta_arch == NULL) {
		printf("discreta_init(): WARNING: $DISCRETA_ARCH not set !\n");
		fflush(stdout);
		}
	if (discreta_home) {
		sprintf(str, "%s/lib/this_is", discreta_home);
		if (file_size(str) <= 0) {
			printf("discreta_init(): WARNING: "
				"can't find my library (DISCRETA_HOME/lib) !\n");
			fflush(stdout);
			}
		}
	return OK;
}

#if TEXDOCU
INT discreta_exit()
#else
freeing local data of DISCRETA
should be called before exiting a DISCRETA program
#endif
{
	INT i;
	INT erg = OK;
	
#ifdef PARTTRUE
	part_ende();
#endif
#ifdef DB_TRUE
	database_exit();
#endif
	for (i = speicherposition; i >= 0; i--)
		my_free(speicher[i]);
	my_free(speicher);
	fflush(stdout);
	fflush(stderr);
	return erg;
}

#if TEXDOCU
INT SYM_OB::swap(SYM_OP b)
#else
swaps the objects this and b
#endif
/* a becomes b and b becomes a */
{
	struct object c;
	if (this == b)
		return error("swap:identical");
	c = *(OP)this; 
	*(OP)this = *(OP)b; 
	*(OP)b = c; 
	return(OK);
}

#if TEXDOCU
INT freeall(SYM_OP a)
#else
frees the self part of object a
and afterwards deallocates the menory 
for the object itself. 
This is the correct way to destroy objects 
which were previously allocated via callocobject()
#endif
{ 
	INT erg = OK;
	
	if (! a->emptyp()) 
		if (a->s_obj_k() != INTEGER)
			erg += a->freeself(); 

	if (speicherposition + 1 < SPEICHERSIZE)
		speicher[++speicherposition] = a;
	else
		my_free(a); 
	if (erg != OK)
		return error("freeall: error in computing");
	return erg;
}


#if TEXDOCU
INT SYM_OB::freeself()
#else
frees the self part of the object.
The object gets the objectkind EMPTY
#endif
{
	INT erg = OK;
	
	if (emptyp()) 
		return OK;
	
	switch (s_obj_bk()) {
		case INTEGER :
			((INTEGER_OP)this)->freeself();
			break;
		case MEM:
			erg += ((MEM_OP)this)->freeself();
			break;
		case STRING:
			erg += ((STRING_OP)this)->freeself();
			break;
		case VECTOR: 
		case GROUP_SELECTION_KIND :
			erg += ((VECTOR_OP)this)->freeself();
			break;
		case MATRIX :
			erg += ((MATRIX_OP)this)->freeself();
			break;
#ifdef UNIPOLYTRUE
		case UNIPOLY :
			erg += ((UNIPOLY_OP)this)->freeself();
			break;
#endif
#ifdef BRUCHTRUE
		case BRUCH :
			erg += ((BRUCH_OP)this)->freeself();
			break;
#endif
#ifdef LISTTRUE
		case POLYNOM:
		case LIST:
			erg += ((LIST_OP)this)->freeself();
			break; 
#endif
#ifdef LONGINTTRUE
		case LONGINT :
			erg += ((LONGINT_OP)this)->freeself();
			break; 
#endif
#ifdef MONOMTRUE
		case MONOM :
			erg += ((MONOM_OP)this)->freeself();
			break;
#endif
#ifdef PARTTRUE
		case AUG_PART : 
		case PARTITION :
			erg += ((PARTITION_OP)this)->freeself();
			break;
#endif
#ifdef PERMTRUE
		case PERMUTATION :
			erg += ((PERMUTATION_OP)this)->freeself();
			break;
#endif
		default:
			{
			printobjectkind();
			{
			INT i = s_obj_k();
			printf("obj_k = %ld\n", i);
			return error("freeself:wrong type");
			}
			}
		};
	c_obj_k(EMPTY);
	if (erg != OK) 
		return error("freeself: error in computing");
	return erg;
}

#if TEXDOCU
INT SYM_OB::freeself_debug(INT print_depth)
#else
frees the self part of the object.
The object gets the objectkind EMPTY
#endif
{
	INT erg = OK;
	
	if (emptyp()) 
		return OK;
	if (print_depth > 0) {
		BYTE str[1000];
		
		str[0] = 0;
		sprintobjectkind(str);
		printf("SYM_OB::freeself_debug for: %s object", str);
		printf("\n");
		fflush(stdout);
		}
	switch (s_obj_bk()) {
		case INTEGER :
			((INTEGER_OP)this)->freeself();
			break;
		case MEM:
			erg += ((MEM_OP)this)->freeself();
			break;
		case STRING:
			erg += ((STRING_OP)this)->freeself();
			break;
		case VECTOR: 
		case GROUP_SELECTION_KIND :
			erg += ((VECTOR_OP)this)->freeself_debug(print_depth);
			break;
		case MATRIX :
			erg += ((MATRIX_OP)this)->freeself();
			break;
#ifdef UNIPOLYTRUE
		case UNIPOLY :
			erg += ((UNIPOLY_OP)this)->freeself();
			break;
#endif
#ifdef BRUCHTRUE
		case BRUCH :
			erg += ((BRUCH_OP)this)->freeself();
			break;
#endif
#ifdef LISTTRUE
		case POLYNOM:
		case LIST:
			erg += ((LIST_OP)this)->freeself();
			break; 
#endif
#ifdef LONGINTTRUE
		case LONGINT :
			erg += ((LONGINT_OP)this)->freeself();
			break; 
#endif
#ifdef MONOMTRUE
		case MONOM :
			erg += ((MONOM_OP)this)->freeself();
			break;
#endif
#ifdef PARTTRUE
		case AUG_PART : 
		case PARTITION :
			erg += ((PARTITION_OP)this)->freeself();
			break;
#endif
#ifdef PERMTRUE
		case PERMUTATION :
			erg += ((PERMUTATION_OP)this)->freeself();
			break;
#endif
		default:
			{
			printobjectkind();
			{
			INT i = s_obj_k();
			printf("obj_k = %ld\n", i);
			return error("freeself:wrong type");
			}
			}
		};
	c_obj_k(EMPTY);
	if (erg != OK) 
		return error("freeself: error in computing");
	return erg;
}

#if TEXDOCU
INT SYM_OB::copy(SYM_OP b)
#else
copies this into b.
b is overwritten 
(it is emptied with a freeself before the copy)
#endif
{
	INT erg = OK;
	
	if (this == b) 
		return(OK);
	if (b == NULL) 
		return error("copy:b == NULL");

	if (!b->emptyp()) 
		erg += b->freeself();
	if (emptyp())	
		return(OK);
	switch (s_obj_bk()) {
		case MEM:
			erg += ((MEM_OP)this)->copy((MEM_OP)b);
			break;
		case STRING:
			erg += ((STRING_OP)this)->copy((STRING_OP)b);
			break;
#ifdef BRUCHTRUE
		case BRUCH :
			erg += ((BRUCH_OP)this)->copy((BRUCH_OP)b);
			break;
#endif
#ifdef INTEGERTRUE
		case INTEGER :
			erg += ((INTEGER_OP)this)->copy((INTEGER_OP)b);
			break;	
#endif
#ifdef LISTTRUE
		case POLYNOM:
		case LIST : 
			erg += ((LIST_OP)this)->copy((LIST_OP)b);
			break;
#endif
#ifdef LONGINTTRUE
		case LONGINT :
			erg += ((LONGINT_OP)this)->copy((LONGINT_OP)b);
			break;
#endif
#ifdef MATRIXTRUE
		case MATRIX :
			erg += ((MATRIX_OP)this)->copy((MATRIX_OP)b);
			break;
#endif
#ifdef UNIPOLYTRUE
		case UNIPOLY :
			erg += ((UNIPOLY_OP)this)->copy((UNIPOLY_OP) b);
			break;
#endif
#ifdef MONOMTRUE
		case MONOM :
			erg += ((MONOM_OP)this)->copy((MONOM_OP)b);
			break;
#endif
#ifdef PARTTRUE
		case AUG_PART :
			erg += ((PARTITION_OP)this)->copy((PARTITION_OP) b);
			b->c_obj_k(AUG_PART);
			break;
		case PARTITION :
			erg += ((PARTITION_OP)this)->copy((PARTITION_OP) b);
			break;
#endif
#ifdef PERMTRUE
		case PERMUTATION : 
			erg += ((PERMUTATION_OP)this)->copy((PERMUTATION_OP)b);
			break;
#endif
#ifdef VECTORTRUE
		case VECTOR :   
			erg += ((VECTOR_OP)this)->copy((VECTOR_OP)b);
			break;
#endif
		default:	
			printobjectkind();
			return error("SYM::copy: wrong type");
		};

	if (erg != OK) {
		error("SYM::copy: error during computation");
		}
	return erg;
}

#if TEXDOCU
INT SYM_OB::nullp()
#else
TRUE if 0, FALSE otherwise
#endif
{
	switch (s_obj_k()) {
		case INTEGER:
			return (((INTEGER_OP)this)->nullp()); 
#ifdef LONGINTTRUE
		case LONGINT:
			return ((LONGINT_OP)this)->nullp();
#endif
#ifdef BRUCHTRUE
		case BRUCH:
			return ((BRUCH_OP)this)->nullp();
#endif
#ifdef MATRIXTRUE
		case MATRIX:
			return ((MATRIX_OP)this)->nullp();
#endif
#ifdef UNIPOLYTRUE
		case UNIPOLY :
			return ((UNIPOLY_OP)this)->nullp();
#endif
#ifdef POLYTRUE
		case POLYNOM:
			return ((POLYNOM_OP)this)->nullp();
#endif
#ifdef VECTORTRUE
		case VECTOR:
			return (((VECTOR_OP)this)->nullp()); 
#endif
		default:
			return(FALSE);
		}
}

#if TEXDOCU
INT SYM_OB::einsp()
#else
TRUE if 1, FALSE otherwise
#endif
{
	switch (s_obj_k()) {
#ifdef BRUCHTRUE
		case BRUCH: 
			return (((BRUCH_OP)this)->einsp()); 
#endif
		case INTEGER:
			return (((INTEGER_OP)this)->einsp()); 
#ifdef LONGINTTRUE
		case LONGINT:
			return ((LONGINT_OP)this)->einsp();
#endif
#ifdef MATRIXTRUE
		case MATRIX:
			return ((MATRIX_OP)this)->einsp();
#endif
#ifdef UNIPOLYTRUE
		case UNIPOLY :
			return ((UNIPOLY_OP)this)->einsp();
#endif
#ifdef PERMTRUE
		case PERMUTATION:
			return (((PERMUTATION_OP)this)->einsp()); 
#endif
#ifdef POLYTRUE
		case POLYNOM:
			return ((POLYNOM_OP)this)->einsp();
#endif
		case VECTOR:
			return (((VECTOR_OP)this)->einsp()); 
		default:
			return(FALSE);
		}
}

#if TEXDOCU
INT SYM_OB::negeinsp()
#else
TRUE if -1, FALSE otherwise
#endif
{
	switch (s_obj_k()) {
#ifdef INTEGERTRUE
		case INTEGER:
			return (((INTEGER_OP)this)->negeinsp()); 
#endif
#ifdef BRUCHTRUE
		case BRUCH:
			return (((BRUCH_OP)this)->negeinsp()); 
#endif
		default:
			return(FALSE);
		};
	}

#if TEXDOCU
INT SYM_OB::sym_comp(SYM_OP b)
#else
compares this and b.
returns -1 / 0 / 1
#endif
{
	if (emptyp() && b->emptyp())
		return(0);
	else if (emptyp()) return(-1);
	else if (b->emptyp()) return(1);
	else switch(s_obj_bk()) {
		case STRING:
			return(strcmp(((STRING_OP)this)->s_str(), 
				((STRING_OP)b)->s_str()));
#ifdef BRUCHTRUE
		case BRUCH : 
			return ((BRUCH_OP)this)->compare(b);
#endif
#ifdef INTEGERTRUE
		case INTEGER : 
			if (b->s_obj_k() == INTEGER)
				return ((INTEGER_OP)this)->comp_integer((INTEGER_OP)b);
			else
				return ((INTEGER_OP)this)->compare(b);
#endif
#ifdef LONGINTTRUE
		case LONGINT : 
			return ((LONGINT_OP)this)->compare(b);
#endif
#ifdef MATRIXTRUE
		case MATRIX : 
			return ((MATRIX_OP)this)->compare((MATRIX_OP)b);
#endif
#ifdef UNIPOLYTRUE
		case UNIPOLY :
			return error("compare not yet implemented for UNIPOLYs !");
#endif
#ifdef MONOMTRUE
		case MONOM :
			return ((MONOM_OP)this)->compare((MONOM_OP)b);
#endif
#ifdef LISTTRUE
		case POLYNOM:
		case LIST :
			return ((LIST_OP)this)->compare((LIST_OP)b);
#endif
#ifdef PARTTRUE
		case PARTITION:
			return ((PARTITION_OP)this)->compare((PARTITION_OP)b);
#endif
#ifdef PERMTRUE
		case PERMUTATION: 
			return ((PERMUTATION_OP)this)->compare((PERMUTATION_OP)b);
#endif
#ifdef VECTORTRUE
		case VECTOR : 
			return ((VECTOR_OP)this)->compare((VECTOR_OP)b);
#endif
		default:{
			printobjectkind();
			error("SYM::sym_comp:not possible");
			return(ERROR);
			}
		}
}

#if TEXDOCU
INT SYM_OB::eq(SYM_OP b)
#endif
{
	if (sym_comp(b) == 0) return(TRUE);
	return(FALSE);
}

#if TEXDOCU
INT SYM_OB::neq(SYM_OP b)
#endif
{
	if (sym_comp(b) != 0) return(TRUE);
	return(FALSE);
}

#if TEXDOCU
INT SYM_OB::ge(SYM_OP b)
#endif
{
	if (sym_comp(b) >= 0) return(TRUE);
	return(FALSE);
}

#if TEXDOCU
INT SYM_OB::gt(SYM_OP b)
#endif
{
	if (sym_comp(b) > 0) return(TRUE);
	return(FALSE);
}

#if TEXDOCU
INT SYM_OB::le(SYM_OP b)
#endif
{
	if (sym_comp(b) <= 0) return(FALSE);
	return(TRUE);
}

#if TEXDOCU
INT SYM_OB::lt(SYM_OP b)
#endif
{
	if (sym_comp(b) < 0) return(TRUE);
	return(FALSE);
}

#if TEXDOCU
INT SYM_OB::init(OBJECTKIND kind)
#endif
{
	INT erg = OK;
	
	if (kind == EMPTY) 
		return error("SYM::init:EMPTY");
	if (!emptyp()) 
		erg += freeself();
	switch (kind) {
		case MEM:
			erg += ((MEM_OP)this)->init(0, NIL); break;
		case STRING:
			erg += ((STRING_OP)this)->init(""); break;
#ifdef BRUCHTRUE
		case BRUCH: 
			erg += ((BRUCH_OP)this)->b_ou(callocobject("SYM::init(BRUCH)"), callocobject("SYM::init(BRUCH)"));
				break;
#endif
		case INTEGER: 
			c_obj_k(kind);
			erg += ((INTEGER_OP)this)->zero();
			break;
#ifdef LONGINTTRUE
		case LONGINT:
			erg += ((LONGINT_OP)this)->init();
			break;
#endif
#ifdef MATRIXTRUE
		case MATRIX:
			return error("can't initialize MATRIX");
			break;
#endif
#ifdef MONOMTRUE
		case MONOM:
			erg += ((MONOM_OP)this)->b_sk((VECTOR_OP) callocobject("SYM::init(MONOM)"), 
				callocobject("SYM::init(MONOM)"));
				break;
#endif
#ifdef PARTTRUE
		case PARTITION: 
			erg += ((PARTITION_OP)this)->b_ks(
				VECTOR, (VECTOR_OP) callocobject("SYM::init(PARTITION)"));
				break;
#endif
#ifdef PERMTRUE
		case PERMUTATION: 
			return error("can't initialize PERMUTATION");
#endif
#ifdef LISTTRUE
		case POLYNOM:
		case LIST: 
			erg += ((LIST_OP)this)->b_sn(NULL, NULL);
			c_obj_k(kind);
			break;
#endif
#ifdef VECTORTRUE  
		case VECTOR:
		case COMP: 
			erg += ((VECTOR_OP)this)->b_ls(
				(INTEGER_OP)callocobject("SYM::init(VECTOR"), NULL); 
			((VECTOR_OP)this)->s_l()->m_i(0);
			c_obj_k(kind);
			break;
#endif
		default: 
			fprintf(stderr,"kind = %ld\n",(INT) kind);
			return error("SYM::init: wrong kind");
		}
	
	if (erg != OK) {
		error("SYM::init: error during computation");
		}
	return erg;
}

#if TEXDOCU
INT SYM_OB::listp()
#endif
{
	OBJECTKIND kind = s_obj_k();
	
	if (kind == LIST || 
		kind == POLYNOM || 
		kind == MONOMIAL)
		return(TRUE);
	else
		return(FALSE);
}
	    
#if TEXDOCU
INT copy_ab(SYM_OP a, SYM_OP b)
#endif
{
	return a->copy(b);
}

#if TEXDOCU
INT compare_ab(SYM_OP a, SYM_OP b)
#endif
{
	return a->sym_comp(b);
}

#if TEXDOCU
INT mult_abc(SYM_OP a, SYM_OP b, SYM_OP c)
#endif
{
	return a->mult(b, c);
}

#if TEXDOCU
INT addinvers_ab(SYM_OP a, SYM_OP b)
#endif
{
	return a->addinvers(b);
}

#if TEXDOCU
INT addinvers_apply_a(SYM_OP a)
#endif
{
	return a->addinvers_apply();
}

#if TEXDOCU
SYM_OP callocobject(BYTE *memory_label)
#else
allocation routine for free-memory objects. 
for example, a very short 
DISCRETA program could look like this:
//CODE
main()
{
discreta_init();
{
SYM_OP a = callocobject("object a");
// do something with a
freeall(a);
}
discreta_exit();
}
///CODE

but !
normally, in C++ objects are automatically allocated 
and destroyed. So, the previous example could be 
rewritten as
//CODE
main()
{
discreta_init();
{
SYM_OB a;
// do something with a
}
discreta_exit();
}
///CODE
#endif
{
	SYM_OP c;

	if (speicherposition >= 0)
		c = speicher[speicherposition--];
	else
		c = (SYM_OP) my_malloc((unsigned long) sizeof(struct object), memory_label);
	if (c == NULL) error("callocobject: no memory");
	c->c_obj_k(EMPTY);
	return(c);
}

#if TEXDOCU
INT SYM_OB::power_int(SYM_OP res, INT l)
#else
Lueneburg~\cite{Lueneburg87a}:
$RussPower(a, b, l) := a * b^l;$. 
\begin{enumerate}
\item
$RussPower(a, b, - l) = RussPower(a, 1 / b, l);$
\item
$RussPower(a, b, l + 1) = a * b^{l + 1} = a * b * b^l = RussPower(a * b, b, l);$
\item
$RussPower(a, b, 2 * l) = a * b^{2 * l} = a * {(b^2)}^l = RussPower(a, b^2, l);$
\item
$RussPower(a, b, 0) = a * b^0 = a;$
\end{enumerate}
this algorithm is sometimes also called
\lq repeated squaring and multiplying\rq
#endif
{
	SYM_OB a, b, c, d;

	if (l == 0) {
		copy(&b);
		b.one();
		b.swap(res);
		return OK;
		}
	copy(&a);
	a.one();
	copy(&b);
	if (l <= 0) {
		return error("SYM_OB::power_int(): l <= 0");
		}
	while (l != 0) {
		if (ODD(l)) {
			/* Rule ii): 
			 * a := a * b; l--; */
			a.mult(&b, &c);
			c.swap(&a);
			l--;
			continue; /* l kann 0 geworden sein. */
			}
		/* now: EVEN(l) and l != 0 */
		/* Rule iii): 
		 * b := b * b; l := l / 2; */
		b.mult(&b, &c);
		c.swap(&b);
		l >>= 1;
		}
	a.swap(res);
	return OK;
}

#if TEXDOCU
INT SYM_OB::power_int_apply(INT l)
#endif
{
	SYM_OB res;
	
	power_int(&res, l);
	swap(&res);
	return OK;
}




#if TEXDOCU
static INT init_base_kind()
#endif
{
	INT l;
	
	for (l = 0; l < SYM_MAX_KIND; l++) {
		base_kind[l] = EMPTY;
		}
	base_kind[EMPTY]                   = EMPTY;         /* 0 */
	base_kind[INTEGER]                 = INTEGER;       /* 1 */
	base_kind[VECTOR]                  = VECTOR;        /* 2 */
	base_kind[PARTITION]               = PARTITION;     /* 3 */
	base_kind[BRUCH]                   = BRUCH;         /* 4 */
	;
	base_kind[PERMUTATION]             = PERMUTATION;   /* 6 */
	;
	base_kind[TABLEAUX]                = TABLEAUX;      /* 8 */
	base_kind[POLYNOM]                 = POLYNOM;       /* 9 */
	base_kind[SCHUR]                   = SCHUR;         /* 10 */
	base_kind[MATRIX]                  = MATRIX;        /* 11 */
	base_kind[AUG_PART]                = AUG_PART;      /* 12 */
	base_kind[HOM_SYM]                 = HOM_SYM;       /* 13 */
	base_kind[SCHUBERT]                = SCHUBERT;      /* 14 */
	base_kind[INTEGERVECTOR]           = INTEGERVECTOR; /* 15 */
	base_kind[KOSTKA]                  = KOSTKA;        /* 16 */
	;
	base_kind[SYMCHAR]                 = SYMCHAR;       /* 18 */
	base_kind[WORD]                    = WORD;          /* 19 */
	base_kind[LIST]                    = LIST;          /* 20 */
	base_kind[MONOM]                   = MONOM;         /* 21 */
	base_kind[LONGINT]                 = LONGINT;       /* 22 */
	;
	base_kind[BINTREE]                 = BINTREE;       /* 24 */
	base_kind[GRAPH]                   = GRAPH;         /* 25 */
	base_kind[COMP]                    = COMP;          /* 26 */
	base_kind[KRANZTYPUS]              = KRANZTYPUS;    /* 27 */
	base_kind[POW_SYM]                 = POW_SYM;       /* 28 */
	base_kind[MONOMIAL]                = MONOMIAL;      /* 29 */
	base_kind[BTREE]                   = BTREE;         /* 30 */
	base_kind[KRANZ]                   = KRANZ;         /* 31 */
	base_kind[GRAL]                    = GRAL;          /* 32 */
	base_kind[ELM_SYM]                 = ELM_SYM;       /* 33 */
	;
	base_kind[FF]                      = FF;            /* 35 */
	base_kind[SGL]                     = VECTOR;        /* 36 */
	base_kind[SGO]                     = VECTOR;        /* 37 */
	;
	base_kind[MEM]                     = MEM;           /* 39 */
	base_kind[DCY]                     = VECTOR;        /* 40 */
	base_kind[41]                      = VECTOR;        /* 41 */ // ZECH
	;
	;
	base_kind[STRING]                  = STRING;        /* 44 */
	;
	;
	;
	;
	;
	base_kind[DATABASE]                = VECTOR;        /* 50 */
	base_kind[BAYERTREE]               = VECTOR;        /* 51 */
	base_kind[ZE]                      = VECTOR;        /* 52 */
	base_kind[CONTI]                   = VECTOR;        /* 53 */
	base_kind[FG]                      = VECTOR;        /* 54 */
	base_kind[55]                      = VECTOR;        /* 55 */
	base_kind[LABRA_KIND]              = VECTOR;        /* 56 */
	;
	base_kind[CODE]                    = VECTOR;        /* 58 */
	;
	;
	base_kind[BT_KEY_KIND]             = VECTOR;        /* 61 */
	base_kind[CODE_ESSENTIALS_KIND]         = VECTOR;   /* 62 */
	base_kind[GED_KIND]                = VECTOR;        /* 63 */
	;
	;
	;
	;
	;
	;
	base_kind[DESIGN_PARAMETER_KIND]          = VECTOR;    /* 70 */
	base_kind[KONTEXT_KIND]                   = VECTOR;    /* 71 */
	base_kind[CLASS_REP_KIND]                 = VECTOR;    /* 72 */
	base_kind[GROUP_CANONIC_FORM_KIND]                 = VECTOR;    /* 73 */
	;
	;
	;
	base_kind[BITVEC_KIND]                     = VECTOR;    /* 77 */
	base_kind[GROUP_SELECTION_KIND]            = VECTOR;    /* 78 */
	base_kind[UNIPOLY]                         = UNIPOLY;   /* 79 */
	base_kind[TREE_KIND]                       = VECTOR;    /* 80 */
	base_kind[GEO_BY_BASE_BLOCKS_KIND]         = VECTOR;    /* 81 */
	base_kind[GENERATORS_KIND]                 = VECTOR;    /* 82 */
	base_kind[DESIGN_PARAMETER_SOURCE_KIND]    = VECTOR;    /* 83 */
	base_kind[SOLID_KIND]                      = VECTOR;    /* 84 */
	base_kind[M_V_KONTEXT_KIND]                = VECTOR;    /* 85 */
	base_kind[L_F_KONTEXT_KIND]                = VECTOR;    /* 86 */
	;
	;
	;
	;
	
	/* < 128; see iof.c */
	
	return OK;
}

