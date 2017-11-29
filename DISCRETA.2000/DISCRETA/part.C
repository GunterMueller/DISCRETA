/* part.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef PARTTRUE

#include <DISCRETA/part.h>
#include <DISCRETA/ma.h>

static struct partition * callocpartition(void);
static INT rh_kostka_kostka(INT *um, INT *ziel, INT *inh, 
	INT k, INT i, INT zahl, INT st, INT len, INT n, INT *zaehler);

static struct partition *callocpartition()
{
	struct partition *erg;

	erg = (struct partition *) my_malloc(sizeof(struct partition), "callocpartition");

	if (erg == NULL) 
		error("callocpartition: no memory");

	return(erg);
}

INT part_ende()
{
	return(OK);
}

INT PARTITION_OB::freeself()
{
	OBJECTSELF d;

	if (s_obj_k() != PARTITION) {
		((SYM_OP) this)->freeself();
		return OK;
		}
	freeall(s_s());
	d = s_obj_s();
	my_free(d.ob_partition);
	c_obj_k(EMPTY);
	return OK;
}

#if TEXDOCU
INT PARTITION_OB::m_ks(OBJECTKIND kind, VECTOR_OP self)
#else
make\_kind.self\_partition: allocates a partition of the given kind 
(which can be VECTOR or EXPONENT). 
The vector self is copied into the new partition object.
#endif
{
	VECTOR_OP s = NULL;
	
	if (self != NULL) { 
		s = (VECTOR_OP) callocobject("PART::m_sk");
		self->copy(s);
		}
	return b_ks(kind, s);
}

#if TEXDOCU
INT PARTITION_OB::b_ks(OBJECTKIND kind, VECTOR_OP self)
#else
build\_kind\_self\_partition: allocates a partition of the given kind 
(which can be VECTOR or EXPONENT). 
The vector self becomes part of the new partition object.
#endif
{
	OBJECTSELF d;
	
	if (! emptyp())
		freeself();
	d.ob_partition = callocpartition();
	c_obj_k(PARTITION);
	c_obj_s(d);
	c_k(kind);
	c_s(self);
	return OK;
}

#if TEXDOCU
INT PARTITION_OB::m_kli(OBJECTKIND kind, INT l)
#else
make kind length integer: allocates a new partition object of the 
given kind. A new vector is allocated for the partition components.
#endif
{
	INT erg = OK;
	
	erg += b_ks(kind, (VECTOR_OP) callocobject("PART::m_kli"));
	erg += s_s()->m_il(l);
	return erg;
}

#if TEXDOCU
INT PARTITION_OB::b_i(INTEGER_OP l)
#else
Allocates a new partition of kind VECTOR of length 1. 
The integer value of $l$ becomes the part of the partition.
Example: $5 \mapsto [5]$.
#endif
{
	INT erg = OK;
	
#ifdef UNDEF
	/* ist auch erlaubt z.B. Hall Littlewood */
	if (negp(integer)) 
		return error("b_i_pa: negativ entry");
	if (nullp(integer)) 
		return error("b_i_pa: null entry");
#endif /* UNDEF */
	erg += m_kli(VECTOR, 1);
	m_ii(0, l->s_i());
	return erg;
}

#if TEXDOCU
INT PARTITION_OB::partitionp()
#else
TRUE if the object is of type PARTITION and all entries are $> 0$.
If the partition is of type VECTOR, all entries must be strictly increasing.
#endif
{
	INT i;
	
	if (s_obj_k() != PARTITION)
		return FALSE;
	if (s_k() == VECTOR) {
		INT m = 1;
		for (i = 0; i < s_li(); i++) {
			if (s_i(i)->s_obj_k() != INTEGER)
				return FALSE;
			if (s_ii(i) < m)
				return FALSE;
			m = s_ii(i);
			}
		return TRUE;
		}
	if (s_k() == EXPONENT)
		for (i = 0; i < s_li(); i++)
			if (s_i(i)->s_obj_k() != INTEGER)
				return FALSE;
	return FALSE;
}

#if TEXDOCU
INT PARTITION_OB::weight_i()
#else
Returns the weight of the partition, i.~e. 
$n$ if $p\vdash n$ is a partition of $n$.
#endif
{
	INT i, w = 0;
	
	if (s_obj_k() != PARTITION)
		return ERROR;
	if (s_k() == VECTOR) {
		for (i = s_li() - 1; i >= 0; i--)
			w += s_ii(i);
		return(w);
		}
	if (s_k() == EXPONENT) {
		for (i = s_li() - 1; i >= 0; i--)
			w += (i + 1) * s_ii(i);
		return(w); 
		}
	error("PART::weight_i: wrong kind of part");
	return -1;
}

#if TEXDOCU
INT PARTITION_OB::weight_augmented_i()
#else
#endif
/* before: weight_augpart() */
{
	INT i, k = 0;
	
	if (s_obj_k() != AUG_PART)
		return ERROR;
	for (i = s_li() - 1; i >= 0; i--)
		k = k + s_ii(i) - i;

	return k;
}

#if TEXDOCU
INT PARTITION_OB::length_i()
#else
Gives the number of parts of the partition.
#endif
{
	INT i, l = 0;
	
	if (s_obj_k() != PARTITION)
		return ERROR;
	if (s_k() == VECTOR)
		return s_li();
	if (s_k() == EXPONENT) {
		for (i = 0; i < s_li(); i++)
			l += s_ii(i);
		return l;
		}
	error("PART::length_i: wrong kind of part");
	return -1;
}

#if TEXDOCU
INT PARTITION_OB::strictp()
#else
true if no equal parts
#endif
{
	INT i;

	if (s_obj_k() != PARTITION)
		return error("PART::strictp: wrong type of object");
	if (s_k() == VECTOR) {
		for (i = 1; i < s_li(); i++)
			if (s_ii(i) == s_ii(i - 1))
				return FALSE;
		return TRUE;
		}
	else if (s_k() == EXPONENT) {
		for (i = 0; i < s_li(); i++)
			if (s_ii(i) > 1)
				return FALSE;
		return TRUE;
		}
	else
		return error("PART::strictp: wrong type of part");
}

INT PARTITION_OB::copy(PARTITION_OP b)
{
	INT erg = OK;
	
	erg += b->b_ks(s_k(), (VECTOR_OP) callocobject("PART::copy"));
	erg += s_s()->copy(b->s_s());
	return erg;
}

INT PARTITION_OB::fprint(FILE *f)
{
	INT i;
	
	if (s_k() == FROBENIUS)
		return s_s()->fprint(f);
	for (i = 0; i < s_li(); i++)
		if (s_ii(i) < 10) {
			/* Partitionsteile kleiner 10 
			 * werden als Zahlen geschrieben */
			fprintf(f, "%ld", s_ii(i));
			}
		else if (s_ii(i) < 16) {
			/* partitionsteile von 10 bis 15 werden als 
			 * A,B,C,D,E,F geschrieben */
			fprintf(f, "%c", (BYTE)(s_ii(i) + 55));
			}
		else {
			/* sonst werden die Teile als Zahl mit 
			 * abschliessenden senkrechten 
			 * Strich geschrieben */
			fprintf(f, "%c%ld", '|', s_ii(i));
			};
	return(OK);
}

INT PARTITION_OB::sprint(BYTE *str)
/* appends to str. 
 * writes to maximal strlength of 200. */
{
	INT i;
	BYTE str1[256];
	
	strcat(str, "[");
	for (i = 0; i < s_li(); i++) {
		str1[0] = 0;
		sprintf(str1, "%ld", s_ii(i));
		if (strlen(str) + strlen(str1) < 200)
			strcat(str, str1);
		else
			return OK;
		if (i < s_li() - 1)
			strcat(str, ", ");
		}
	strcat(str, "]");
	return(OK);
}

INT PARTITION_OB::sprint_latex(BYTE *str)
{
	INT i, ai, f_need_space;
	BYTE str1[256];
	
	strcat(str, "[");
	if (s_k() == VECTOR) {
		for (i = 0; i < s_li(); i++) {
			str1[0] = 0;
			sprintf(str1, "%ld", s_ii(i));
			if (strlen(str) + strlen(str1) < 200)
				strcat(str, str1);
			else
				return OK;
			if (i < s_li() - 1)
				strcat(str, ", \\,");
			}
		}
	else {
		f_need_space = FALSE;
		for (i = 1; i <= s_li(); i++) {
			ai = s_ii(i - 1);
			if (ai == 0)
				continue;
			str1[0] = 0;
			if (f_need_space)
				strcpy(str1, ", \\,");
			sprintf(str1 + strlen(str1), "%ld^{%ld}", i, ai);
			f_need_space = TRUE;
			if (strlen(str) + strlen(str1) < 200)
				strcat(str, str1);
			else
				return OK;
			}
		}
	strcat(str, "]");
	return(OK);
}

INT PARTITION_OB::latex(FILE *fp)
{
	BYTE str[1024];

	str[0] = 0;
	sprint_latex(str);
	fprintf(fp, "%s", str);
	return OK;
}

INT PARTITION_OB::compare(PARTITION_OP b)
{
	INT i;
	INT erg;
	/* char *ac, *bc; */
	
	if (s_obj_k() != PARTITION)
		return error("PART::compare: a not PARTITION");
	if (b->s_obj_k() != PARTITION)
		return error("PART::compare: b not PARTITION");
		
	if (s_k() != b->s_k()) 
		return error("PART::compare: different kind of partitions");

	if (s_k() == VECTOR ) {
		for (i = 0; i < s_li(); i++) {
			if (i >= b->s_li())
				return 1;
			if (s_ii(i) > b->s_ii(i))
				return 1;
			if (s_ii(i) < b->s_ii(i))
				return -1;
			}
		if (i < b->s_li())
			return -1;

#if FALSE
		ac = (char *) S_V_S(S_PA_S(a));
		bc = (char *) S_V_S(S_PA_S(b));
		if (S_PA_LI(a) == S_PA_LI(b))
			{
/*
			printf("fall1: ");
			for (j=0;
				j<S_PA_LI(a)*(sizeof(struct object));j++) 
				printf("%d%d",ac[j],bc[j]);
*/
			erg =  (INT)memcmp(ac,bc,
				( sizeof(struct object) * S_PA_LI(a) ));
			goto cpende;
			}
		if (S_PA_LI(a) < S_PA_LI(b))
			{
/*
			printf("fall2: ");
			for (j=0;
				j<S_PA_LI(a)*(sizeof(struct object));j++) 
				printf("%d%d",ac[j],bc[j]);
*/
			erg = (INT) memcmp(ac,bc,
				(sizeof(struct object) * S_PA_LI(a) ));
			if (erg == 0L)  erg = -1L;
			goto cpende;
			}
		if (S_PA_LI(a) > S_PA_LI(b))
			{
/*
			printf("fall3: ");
			for (j=0;
				j<S_PA_LI(b)*(sizeof(struct object));j++) 
				printf("%d%d",ac[j],bc[j]);
*/
			erg = (INT)memcmp(ac,bc,
				(sizeof(struct object) * S_PA_LI(b) ));
			if (erg == 0L)  erg = 1L;
			goto cpende;
			}
#endif

		}
	else if (s_k() == EXPONENT) {
		for (i = 0; i < s_li(); i++) {
			if (i >= b->s_li()) {
				if (s_ii(i) != 0) {
					erg = 1;
					goto cpende;
					}
				}
			else if (s_ii(i) > b->s_ii(i)) {
				erg = 1;
				goto cpende;
				}
			else if (s_ii(i) < b->s_ii(i)) {
				erg = -1;
				goto cpende;
				}
			}
		
		for (; i < b->s_li(); i++) {
			if (b->s_ii(i) != 0) {
				erg = -1; 
				goto cpende;
				}
			}
		}
	erg = 0;
	goto cpende;
cpende:
/*
	print(a); println(b);
	printf("cpende: %ld\n",erg);
*/
	return erg;
}

#if TEXDOCU
INT PARTITION_OB::conjugate(PARTITION_OP b)
#else
Computes the conjugate partition into $b$. 
Algorithm according to MacDonald. 
Not much faster than conjugatepartition.
#endif
{
	INT i, j, k = 0, m;
	/* k ist die Adresse an der in b geschrieben wird. */

	if (s_obj_k() != PARTITION)
		return ERROR;
#if FALSE /* !!! */
	if (s_k() == EXPONENT) {
		PARTITION_OP c = (PARTITION_OP) callocobject("PART::conjugate");
		erg += t_EXPONENT_VECTOR(c);
		erg += c->conjugate(b);
		erg += freeall(c);
		erg += b->t_VECTOR_EXPONENT(b);
		return erg;
		}
#endif
	if (s_k() != VECTOR)
		return ERROR;

	b->b_ks(VECTOR, (VECTOR_OP) callocobject("PART::conjugate"));
	j = s_li() - 1;
	b->s_s()->m_il(s_ii(j));

	/* dies sind die Adressen 
	 * in den beiden Partitionen */
	m = b->s_li() + s_li() + 1;
	/* dies ist die Laenge der Permutation + 1 */
	for (i = m - 1; i > 0; i--) {
		if (j >= 0) {
			if (i == s_ii(j) + j + 1)
				j--;
			else {
				b->m_ii(k, m - i - k - 1);
				k++;
				}
			}
		else {
			b->m_ii(k, m - i - k - 1);
			k++;
			}
		}
	return OK;
}

#if TEXDOCU
INT PARTITION_OB::equal_parts(INTEGER_OP b)
#else
return TRUE if the partition $a$ has $\ge b$ equal parts
#endif
{
	INT i, j = 0, k = 0;
	
	for (i = 0; i < s_li(); i++) {
		if (s_ii(i) == k)
			j++;
		else {
			k = s_ii(i);
			j = 1;
			}
		if (j == b->s_i())
			return TRUE;
		}
	return FALSE;
}

#if TEXDOCU
INT PARTITION_OB::dimension(INTEGER_OP b)
#else
Dimension of the representation of the $S_n$ labelled by the 
given partition. Hook-fomula.
#endif
{
	INTEGER_OB zaehler, nenner, zw;
	INT i, j, k;
	INT erg = OK;

	if (s_obj_k() != PARTITION)
		return ERROR;
	if (s_k() == EXPONENT) {
		return error("PARTITION::dimension(): s_k() == EXPONENT");
#if FALSE /* !!! */
		PARTITION_OP c;
		
		c = (PARTITION_OP) callocobject("PART::dimension");
		erg += t_EXPONENT_VECTOR(c);
		erg += c->dimension(b);
		erg += freeall(c);
		return erg;
#endif
		}
	if (s_k() != VECTOR)
		return ERROR;
	zw.m_i(weight_i());

	erg += zw.fakul(&zaehler);
	nenner.m_i(1);
	for (i = 0; i < s_li(); i++) {
		for (j = 0; j < s_ii(s_li() - 1 - i); j++) {
			k = hook_length_i(i, j);
			if (k != 1) {
				zw.m_i(k);
				erg += zw.mult_apply(&nenner);
				}
			}
		}
	erg += zaehler.ganzdiv(&nenner, b);
		/* statt div AK 170889 */
	if (erg != OK)
		error("PART::dimension: "
		"error during computation");
	return erg;
}

#if TEXDOCU
INT PARTITION_OB::dimension_augmented(INTEGER_OP b)
#else
#endif
/* before: dimension_augpart() */
/* b becomes the Dimension of the 
 * corresponding irred. representation */
{
	INTEGER_OB zaehler, nenner, zw;
	
	INT i, j, k, erg = OK;


	if (! b->emptyp())
		b->freeself();
	if (s_li() == 1) {
		b->m_i(1);
		return OK;
		}
	if (s_ii(s_li() - 1) == s_li()) { /* 1^n */
		b->m_i(1);
		return OK;
		}
	if (s_ii(s_li() - 2) == s_li() - 2) { /* n */
		b->m_i(1);
		return OK;
		}
	if (s_li() == 2) {
		if (s_ii(0) == 1) {
			b->m_i(s_ii(1) - 1);
			return OK;
			}
		}
	zw.m_i(weight_augmented_i());
	erg += zw.fakul(&zaehler);
	nenner.m_i(1);
	for (i = 0; i < s_li(); i++) {
		for (j = 0; j < s_ii(i) - i; j++) {
			k = hook_length_augmented_i(i, j);
			if (k != 1) {
				zw.m_i(k);
				erg += zw.mult_apply(&nenner);
				}
			}
		}

	erg += zaehler.ganzdiv(&nenner, b);
		/* statt div AK 170889 */
	return erg;
}

#if TEXDOCU
INT PARTITION_OB::hook_length_i(INT i, INT j)
#else
Hook length as integer.
#endif
{
	INT e, k;
	
	if (s_obj_k() != PARTITION) {
		return error("PART::hook_length: no part");
		}
	if (s_k() == EXPONENT) {
		PARTITION_OB c;
		
		t_EXPONENT_VECTOR(&c);
		e = c.hook_length_i(i, j);
		return e;
		}
	if (s_k() != VECTOR)
		return ERROR;

	if (i >= s_li())
		return 0;
	e = s_ii(s_li() - 1 - i);
	if (e <= j)
		return 0;
	e -= j;
	/* nun noch die zeilen dazu */
	for (k = i + 1; k < s_li(); k++) 
		if (s_ii(s_li() - 1 - k) -1 >= j)
			e++;
		else
			break;
	return e;
}

#if TEXDOCU
INT PARTITION_OB::hook_length_augmented_i(INT i, INT j)
#else
#endif
/* before: hook_length_augpart */
{
	INT e, k;
	
	if (i >= s_li())
		return 0;
	if (j >= s_ii(i) - i)
		return 0;
	e = s_ii(i) - j - i;
	/* nun noch die zeilen dazu */
	for (k = i - 1; k >= 0; k--)
		if (s_ii(k) - 1 - k >= j) 
			e++;
		else
			break;
	return e;
}

#if TEXDOCU
INT PARTITION_OB::t_VECTOR_EXPONENT(PARTITION_OP nach)
#else
Converts a partition from the VECTOR to the EXPONENT format.
Example: $[1,2,2,3] \mapsto [1,2,1,0,0,0,0,0]\simeq [1^1, 2^2, 3^1, 4^0,\ldots, 8^0]$
#endif
{
	INT i, w;
	INT erg = OK;
	
	if (this == nach) {
		PARTITION_OB l;
		l = *this;
		c_obj_k(EMPTY);
		erg += l.t_VECTOR_EXPONENT(nach);
		return erg;
		}
	if (! nach->emptyp()) 
		erg += nach->freeself();
	if (emptyp()) {
		erg += nach->b_ks(EXPONENT, (VECTOR_OP) callocobject("PART::t_VECTOR_EXPONENT"));
		nach->s_s()->m_il_n(1);
		return erg;
		}
	w = weight_i();
	nach->b_ks(EXPONENT, (VECTOR_OP) callocobject("PART::t_VECTOR_EXPONENT"));
	nach->s_s()->m_il_n(w);

	for (i = 0; i < s_li(); i++)
		nach->s_i(s_ii(i) - 1)->inc();

	return erg;
}

#if TEXDOCU
INT PARTITION_OB::t_EXPONENT_VECTOR(PARTITION_OP b)
#else
Converts from EXPONENT type to VECTOR type.
#endif
{
	INT i, j, z;

	if (this == b) {
		PARTITION_OB l;
		
		copy(&l);
		l.t_EXPONENT_VECTOR(b);
		return OK;
		}
	/* sum_vector(S_PA_S(a),l); */
	j = 0;
	for (i = 0; i < s_li(); i++)
		j += s_ii(i);
	if (! b->emptyp())
		b->freeself();
	b->b_ks(VECTOR, (VECTOR_OP) callocobject("PART::t_EXPONENT_VECTOR"));
	b->s_s()->m_il_n(j);
	z = 0;
	for (i = 0; i < s_li(); i++) {
		for (j = 0; j < s_ii(i); j++) {
			b->m_ii(z, i + 1);
			z++;
			}
		}
	return OK;
}

#undef DEBUG_INDUCE2

#if TEXDOCU
INT PARTITION_OB::induce2(PARTITION_OP nach)
#else
Assume $\pi \in S_n$ with cycle\_type($\pi) = p$, the given partition. 
The routine computes the cycle type of $\pi$ acting on the two-subsets 
of $n$ into nach.
#endif
{
	INT n, i, j, l, n2, ai, aj, k, g, m, m1, i2, ai2;

	if (s_k() != EXPONENT) {
		return error("PART::induce2 only for type EXPONENT");
		}
	n = weight_i();
	l = s_li();
	if (l != n)
		return error("PART::induce2 l != n");
	n2 = (n * (n - 1)) >> 1;
	nach->b_ks(EXPONENT, (VECTOR_OP) callocobject("PART::induce2"));
	nach->s_s()->m_il_n(n2);
	for (i = 1; i <= l - 1; i++) {
		ai = s_ii(i - 1);
		if (ai == 0)
			continue;
		for (j = i + 1; j <= l; j++) {
			aj = s_ii(j - 1);
			if (aj == 0)
				continue;
			ggt_iipi(i, j, &g);
			m = ai * aj * g;
			k = i * j / g;
#ifdef DEBUG_INDUCE2
			printf("i = %ld j = %ld ai = %ld aj = %ld g = %ld k = %ld\n", 
				i, j, ai, aj, g, k);
#endif
			if (k > n2)
				return error("PART::induce2(): k > n2");
			m1 = nach->s_ii(k - 1);
			nach->m_ii(k - 1, m1 + m);
			}
		}
	for (i = 1; i <= n2; i++) {
		if (i <= l) {
			ai = s_ii(i - 1);
			if (ai != 0) {
				if (ODD(i)) {
					m = (ai * (i * ai - 1) ) >> 1;
					}
				else {
					m = (ai * (i * ai - 2) ) >> 1;
					}
				m1 = nach->s_ii(i - 1);
				nach->m_ii(i - 1, m1 + m);
				}
			}
		i2 = i << 1;
		if (i2 <= l) {
			ai2 = s_ii(i2 - 1);
#ifdef DEBUG_INDUCE2
			printf("i = %ld i2 = %ld ai2 = %ld\n", i, i2, ai2);
#endif
			if (ai2 != 0) {
				m1 = nach->s_ii(i - 1);
				nach->m_ii(i - 1, m1 + ai2);
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT PARTITION_OB::augment()
#else
adds the staircase $0,1,2,3,\ldots$. Example: $[1,1,1,3] \mapsto [1,2,3,6]$
#endif
/* before: augpart */
{
	INT i;
	
	if (s_obj_k() != PARTITION)
		return error("augment(): s_obj_k() != PARTITION");
	c_obj_k(AUG_PART);
	for (i = 0; i < s_li(); i++)
		c_ii(i, s_ii(i) + i);
	return OK;
}

#if TEXDOCU
INT PARTITION_OB::de_augment()
#else
Subtracts the staircase
#endif
{
	INT i;
	
	if (s_obj_k() != AUG_PART)
		return error("de_augment(): s_obj_k() != PARTITION");
	c_obj_k(PARTITION);
	for (i = 0; i < s_li(); i++)
		c_ii(i, s_ii(i) - i);
	return OK;
}

static INT stripexistp(PARTITION_OP part, INT length, INT i);
static INT addstrip(PARTITION_OP part, INT k, INT i, INT hi);
static INT removestrip(PARTITION_OP part, INT k, INT i);
static INT calculate(INT sign, PARTITION_OP rep, 
	PARTITION_OP part, INTEGER_OP res);

static INT stripexistp(PARTITION_OP part, INT length, INT i)
{
	INTEGER_OP z = part->s_i(i);
	INT h2;

	h2 = z->s_i();

	for (; i >= 0; i--, z--)
		if ( (z->s_i() + length) == h2) 
			return FALSE;
	return TRUE;
}

static INT addstrip(PARTITION_OP part, INT k, INT i, INT hi)
{
	i = i - hi;
	
	/* in l wird angesetzt */
	while ((k--) > 0) {
		if (i == part->s_li() - 1) {
			part->m_ii(i, part->s_ii(i) + k + 1);
			return OK; 
			}
		else if (part->s_ii(i) < part->s_ii(i + 1))
			part->s_i(i)->inc();
		else if (part->s_ii(i) == part->s_ii(i + 1))
			part->s_i(++i)->inc();
		else
			error("addstrip:");
		}
	return OK;
}

static INT removestrip(PARTITION_OP part, INT k, INT i)
/* erzeugt neue Partition part in der ab der Zeile i ein
 * Streifen der Laenge length entfernt wurde .
 * Ergebnis ist die Hakenlaenge */
{
	INT l;
	
	l = i;
	while ((k--) > 0) {
		if (i == 0) 
			part->s_i(0)->dec();
		else if (part->s_ii(i) > part->s_ii(i - 1))
			part->s_i(i)->dec();
		else 	
			part->s_i(--i)->dec();
		};
	return (l - i);
}

static INT calculate(INT sign, PARTITION_OP rep, 
	PARTITION_OP part, INTEGER_OP res)
{
	INT i, hooklength, l;
	INT erg = OK; 

	if (part->s_li() == 0) { 
		if (sign == 1) 
			erg += res->inc(); 
		else if (sign == -1) 
			erg += res->dec();
		else 
			erg += ERROR;
		return erg;
		}
	if (part->s_li() == 1) {
		/* Robinson Lemma 4.11 */
		if (rep->s_li() == 1) {
				res->m_i(1);
				return OK;
				}
		if (rep->s_ii(rep->s_li() - 2) > 
				rep->s_li() - 1 )
				return OK;
		/* rep is haken */
		for (i = 0; i < rep->s_li(); i++)
			if (rep->s_ii(i) > i)
				break;
		i = rep->s_li() - i;
		/* i is laenge der part */
		if (sign == 1)
			if (i % 2 == 0)
				return res->dec();
			else 
				return res->inc();
		else
			if (i % 2 == 0)
				return res->inc();
			else 
				return res->dec();
		}
	if (part->s_ii(part->s_li() - 1) == 1) {
		/* dimension */
		INTEGER_OB newrep;

		erg += rep->dimension_augmented(&newrep);
		if (sign == -1) 
			erg += newrep.addinvers_apply();
		erg += newrep.add_apply(res); 
		return erg;
		}
	l = part->s_li() - 1;
	for (i = rep->s_li() - 1; i >= 0; i--) {
		if (part->s_ii(l) <= rep->s_ii(i)) {
			if 	(stripexistp(rep, part->s_ii(l), i)) { 
				hooklength = removestrip(
					rep, part->s_ii(l), i); 
				erg += part->s_l()->dec();
				erg += calculate( ((hooklength % 2 == 0) ? 
					sign : - sign),
					rep, part, res); 
				erg += part->s_l()->inc();
				erg += addstrip(rep, 
				part->s_ii(l), i, hooklength);
				}
			}
		}
	return erg;
}

#if TEXDOCU
INT charvalue(PARTITION_OP rep, PARTITION_OP part, INTEGER_OP res)
#else
Computes the value of the irreducible representation of $S_n$ labelled by $part$ 
for the conjugacy class of elements labelled by $rep$ into $res$. 
The result is an integer.
#endif
{
	if (rep == part) {
		PARTITION_OB part2;
		
		part->copy(&part2);
		return charvalue(rep, &part2, res);
		}
	rep->augment();
	res->m_i(0);
	calculate(1, rep, part, res);
	rep->de_augment();
	return OK;
}

#if TEXDOCU
INT charvalue_ij(VECTOR_OP V, INT rep_i, INT part_i, INTEGER_OP res)
#else
Computes charvalue(V[i], V[j]).
#endif
{
	PARTITION_OP rep = (PARTITION_OP) V->s_i(rep_i);
	PARTITION_OP part = (PARTITION_OP) V->s_i(part_i);
	return charvalue(rep, part, res);
}

#if TEXDOCU
INT chartafel(INT n, MATRIX_OP M)
#else
Computes the character table of $S_n$.
#endif
/* before: chartafel(OP a, OP res) */
{
	VECTOR_OB V; /* the partitions */
	
	vector_of_part(&V, n);
	chartafel_vop(n, M, &V);
	return OK;
}

INT chartafel_vop(INT n, MATRIX_OP M, VECTOR_OP V)
{
	INT dim, i, j, index;
	PARTITION_OB conj_part;
	PARTITION_OP pi, pj;

	dim = V->s_li();
	M->m_ilih(dim, dim);
	/* alternating character: */
	i = dim - 1;
	for (j = 0; j < dim; j++) {
		charvalue_ij(V, i, j, 
			(INTEGER_OP) M->s_ij(i, j));
		}
	/* one character: */
	for (j = 0; j < dim; j++) {
		M->m_iji(0, j, 1);
		}
	for (i = 0; i < dim; i++) {
		if (!M->s_ij(i, 0)->emptyp())
			/* already computed */
			continue;
		pi = (PARTITION_OP) V->s_i(i);
		for (j = 0; j < dim; j++) {
			pj = (PARTITION_OP) V->s_i(j);
			/* JK Cor 2.4.9: */
			if (pi->s_li() - 1 + pi->s_ii(pi->s_li() - 1) >= 
				pj->s_ii(pj->s_li() - 1)) {
				charvalue_ij(V, i, j, 
					(INTEGER_OP) M->s_ij(i, j));
				}
			else
				M->m_iji(i, j, 0);
			}
		/* associated character: */
		pi->conjugate(&conj_part);
		for (index = i + 1; index < dim; index++)
			if (((SYM_OP) &conj_part)->sym_comp(
				V->s_i(index)) == 0)
				break;
		if (index < dim) {
			/* char[index] := 
			 *      char[i] * alternating character: */
			for (j = 0; j < dim; j++) {
				M->m_iji(index, j, 
					M->s_iji(i, j) * M->s_iji(dim - 1, j));
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT PARTITION_OB::first(INT n)
#else
Computes the first partition of $n$.
#endif
/* before: first_partition(OP n, OP part) */
{
	INTEGER_OB int_ob;
	
	if (n <= 0)
		return error("PARTITION::first: n <= 0");
	int_ob.m_i(n);
	b_i(&int_ob);
	return OK;
}

#if TEXDOCU
INT PARTITION_OB::next_VECTOR_apply()
#else
Computes the next partition in VECTOR type.
Returns TRUE if a new partition is found.
#endif
{
	PARTITION_OB n;
	INT ret;
	
	ret = next_VECTOR(&n);
	swap(&n);
	return ret;
}

#if TEXDOCU
INT PARTITION_OB::next_VECTOR(PARTITION_OP next)
#else
Computes the next partition in VECTOR type into next.
Returns TRUE if a new partition is found.
#endif
/* before: next_part_VECTOR(OP n, OP part) */
/* Nijenhuis ch. 9 */
/* returns TRUE if a new permutation was found. */
{
	INT i;
	
	if (next == this) {
		PARTITION_OB next1;
		INT r;
		
		r = next_VECTOR(&next1);
		swap(&next1);
		return r;
		}
	if (s_ii(0) > 1) {
		next->m_kli(VECTOR, s_li() + 1);
		next->s_i(0)->m_i(1);
		next->s_i(1)->m_i(s_ii(0) - 1);
		for (i = 1; i < s_li(); i++)
			next->s_i(i + 1)->m_i(s_ii(i));
		return TRUE;
		}
	for (i = 1; i < s_li(); i++)
		if (s_ii(i) > 1)
			break;
	if (i == s_li())
		return FALSE;
	{
		INT k, m, n, j, o;
		
		k = s_li() - i;
		m = s_ii(i);
		n = m - 1;
		j = (i + m) / n;
		o = (i + m) % n;
		if (o == 0)
			j--;
		next->m_kli(VECTOR, j + k);
		if (o != 0) {
			next->s_i(0)->m_i(o);
			o = 1;
			}
		for (m = o; m <= j; m++)
			next->s_i(m)->m_i(n);
		for ( ; m < j + k; m++, i++)
			next->s_i(m)->m_i(s_ii(i + 1));
	}
	return TRUE;
}

#if TEXDOCU
INT PARTITION_OB::first_into_k_parts_VECTOR(INT n, INT k)
#else
Computes the first partition of $n$ into exactly $k$ parts.
#endif
{
	if (n <= 0)
		return error("PARTITION::first_into_k_parts_VECTOR: n <= 0");
	if (k <= 0)
		return error("PARTITION::first_into_k_parts_VECTOR: k <= 0");
	if (k > n)
		return error("PARTITION::first_into_k_parts_VECTOR: k > n");
	first(n);
	while (s_li() != k) {
		if (!next_VECTOR_apply())
			return FALSE;
		}
	return TRUE;
}

#if TEXDOCU
INT PARTITION_OB::next_into_k_parts_VECTOR(INT n, INT k)
#else
Computes the next partition of $n$ into exactly $k$ parts.
#endif
{
	PARTITION_OB p;
	
	do {
		if (!next_VECTOR_apply())
			return FALSE;
		} while (s_li() != k);
	return TRUE;
}

#if TEXDOCU
INT vector_of_part(VECTOR_OP V, INT n)
#else
Computes a vector holding all partitions of $n$. 
The partitions are generated via first and next\_VECTOR.
#endif
/* before: makevectorofpart(OP n, OP vec) */
{
	PARTITION_OB next;
	INT i;
	
	if (n < 0)
		return error("vector_of_part(): n < 0");
	if (n == 0)
		return V->m_il(0);
	V->m_il(1);
	next.first(n);
	((SYM_OP) &next)->swap(V->s_i(0));
	i = 0;
	while (((PARTITION_OP) V->s_i(i))->next_VECTOR(&next)) {
		i++;
		V->inc();
		((SYM_OP) &next)->swap(V->s_i(i));
		}
	return OK;
}

#if TEXDOCU
INT PARTITION_OB::class_length_Sn(SYM_OP res)
#else
Computes the length of the conjugacy class in $S_n$ 
of elements labelled by the cycle type which is the given partition.
#endif
{
	SYM_OB a, b;
	
	a.fakul_int(weight_i());
	centralizer_order_Sn(&b);
	a.ganzdiv_integral(&b, res);
	return OK;	
}

#if TEXDOCU
INT PARTITION_OB::centralizer_order_Sn(SYM_OP res)
#else
Computes the centralizer order of the element $\pi \in S_n$ 
which is labelled by the cycle type which is the given partition.
#endif
{
	PARTITION_OB q;
	SYM_OB a, b, c, d, e;
	INT i, j, n;
	
	n = weight_i();
	if (s_k() == VECTOR){
		t_VECTOR_EXPONENT(&q);
		}
	else
		copy(&q);
	a.m_i_i(1);
	for (i = 0; i < n; i++){
		j = q.s_ii(i);
		b.m_i_i(i + 1);
		b.power_int_apply(j);
		c.fakul_int(j);
		b.mult(&c, &d);
		a.mult(&d, &e);
		e.swap(&a);
		}
	a.swap(res);
	return OK;
}

#if 0
#if TEXDOCU
INT PARTITION_OB::class_len(SYM_OP l)
#else
Computes the length of the conjugacy class in $S_n$ 
of elements labelled by the cycle type which is the given partition.
#endif
{
	INT i, g = 1, l = 1, n, sp = 1;
	PARTITION_OB p;

	if (s_k() == VECTOR)
		copy(&p);
	else {
		t_VECTOR_EXPONENT(&q);
		}
	if (s_k() != VECTOR)
		return error("PART::class_len(): wrong kind");
	n = weight_i();
	/* centralizer order: */
	for (i = 0; i < s_li(); i++) {
		if (i > 0) {
			if (s_ii(i) == s_ii(i - 1)) {
				sp++;
				l *= sp;
				}
			else
				sp = 1;
			}
		l *= s_ii(i);
		}
	/* class_len = group_order / centralizer order: */
	for (i = n; i >= 2; i--) {
		if ((l % i) == 0)
			l /= i;
		else
			g *= i;
		}
	g /= l;
	return g;
}
#endif

INT PARTITION_OB::exp_nb_i_parts(INT i)
{
	INT l, a;
	
	l = s_li();
	if (i > l)
		return 0;
	a = s_ii(i - 1);
	return a;
}

#if TEXDOCU
INT kostka_number(PARTITION_OP content, PARTITION_OP umriss, INTEGER_OP res)
#else
#endif
/* Ralf Hager */
{
	INT *um,*hilf,*ziel,*inh;
	INT	i;
	INT	zaehler = (INT)0;
	INT	len,n, wt, wt2;
	PARTITION_OB d;

	res->freeself();
	wt = umriss->weight_i();
	wt2 = content->weight_i();
	if (wt != wt2) 
		return error("kostka_number weight differs");
	um = (INT *) my_malloc((wt + 1) * sizeof(INT), "kostka_number()");
	hilf = (INT *) my_malloc((wt + 1) * sizeof(INT), "kostka_number()");
	ziel = (INT *) my_malloc((wt + 1) * sizeof(INT), "kostka_number()");
	inh = (INT *) my_malloc((wt + 1) * sizeof(INT), "kostka_number()");

	umriss->conjugate(&d);
	for (i = 0; i <= wt; i++) { 
		hilf[i]=0; 
		um[i]=0; 
		ziel[i]=0; 
		inh[i]=0; 
		}
	n= d.s_li();
	for (i = 1; i <= n; ++i) { 
		um[i] = d.s_ii(n - i);
		hilf[i] = um[i];
		}
	um[0] = -1;
	for (i = 0; i < n; ++i) {
		um[i+1] = um[1] + i;
		ziel[i+1] = um[i+1] - hilf[i+1];
		}
	len = content->s_li();
	for (i = 1; i <= len; ++i)
		inh[i] = content->s_ii(len - i);
	my_free(hilf);
		
	rh_kostka_kostka(um,ziel,inh,0,0,inh[1], 1, len, n, &zaehler);
	res->m_i(zaehler);
	my_free(um);
	my_free(inh);
	my_free(ziel);
	return OK;
}

static INT rh_kostka_kostka(INT *um, INT *ziel, INT *inh, 
	INT k, INT i, INT zahl, INT st, INT len, INT n, INT *zaehler)
/* Ralf Hager */
{
	INT	l;

	if (i == zahl)	
		{
	     if (st == len) ++(*zaehler);
	     else rh_kostka_kostka(um, ziel, inh, 0, 0, inh[st+1], st+1, len, n, zaehler);
		}
	else
		{
		for (l = k+1; l <= n; ++l)
			if (((um[l]-1) > um[l-1]) && (um[l] > ziel[l]))
			 {
			 um[l]--;
			 rh_kostka_kostka(um,ziel,inh,l,i+1L,zahl,st,len,n,zaehler);
			 um[l]++;
			 }
		}
	return OK;
}

#if TEXDOCU
INT PARTITION_OB::equal_parts(VECTOR_OP v)
#endif
{
	INT i, j, l, a, b, c;
	VECTOR_OB w;
	
	l = v->s_li();
	if (l == 0) {
		m_kli(VECTOR, 0);
		return OK;
		}
	w.m_il(l + 1);
	j = 0;
	w.m_ii(j++, 0);
	for (i = 1; i < l; i++) {
		if (v->s_i(i)->sym_comp(v->s_i(i - 1)) != 0) {
			w.m_ii(j++, i);
			}
		}
	w.m_ii(j, i);
	m_kli(VECTOR, j);
	for (i = 0; i < j; i++) {
		a = w.s_ii(i);
		b = w.s_ii(i + 1);
		c = b - a;
		m_ii(i, c);
		}
	return OK;
}


#endif /* PARTTRUE */
