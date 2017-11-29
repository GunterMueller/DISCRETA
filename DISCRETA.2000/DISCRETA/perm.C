/* perm.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef PERMTRUE

#include <DISCRETA/perm.h>


#define PERM_DEBUG

/*
 * PERM
 */

static struct permutation * callocpermutation();

INT test_perm()
{
	PERMUTATION_OB p, q, r;
	
	p.m_il(6); p.one();
	q.m_il(6); q.one();
	r.m_il(6); r.one();
	p.Add2Cycle(1, 2);
	q.Add2Cycle(2, 3);
	p.println();
	q.println();
	p.mult(&q, &r);
	r.println();
	return OK;
}

static struct permutation * callocpermutation()
{
	struct permutation *
	b = (struct permutation *) my_malloc(sizeof(struct permutation), "callocpermutation");
	if (b == NULL)
		error("callocpermutation:no mem");
	return b;
}

INT PERMUTATION_OB::freeself()
{
	OBJECTSELF d;
	
	if (!s_s()->emptyp()) {
		s_s()->freeself();
		}
	freeall(s_s()); 
	d = s_obj_s(); 
	my_free(d.ob_permutation); 
	c_obj_k(EMPTY);
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::m_l(INTEGER_OP len)
#else
make length: allocates a permutation of the given length.
length is specified as an integer object.
#endif
{
	return m_il(len->s_i());
}

#if TEXDOCU
INT PERMUTATION_OB::m_il(INT l)
#else
make length: allocates a permutation of the given length.
length is specified as an integer.
#endif
{
	INT erg = OK;
	
	erg += b_ks((OBJECTKIND) VECTOR, 
		(VECTOR_OP)callocobject("PERM::m_il s"));
	if (erg != OK)
		return error("PERMUTATION::m_il(): error in b_ks()");
	erg += s_s()->m_il(l);
	if (erg != OK)
		return error("PERMUTATION::m_il(): error in VECTOR::m_il()");
	return erg;
}

INT PERMUTATION_OB::m_ks(OBJECTKIND kind, VECTOR_OP self)
{
	INT erg = OK;
	
	erg += b_ks(kind, (VECTOR_OP)callocobject("PERM::m_ks s"));
	erg += self->copy(s_s());
	return erg;
}

INT PERMUTATION_OB::b_ks(OBJECTKIND kind, VECTOR_OP self)
{
	OBJECTSELF b;

	((SYM_OP) this)->freeself();
	b.ob_permutation = callocpermutation();
	if (b.ob_permutation == NIL)
		return ERROR;
	((SYM_OP)this)->b_ks(PERMUTATION, b);
	c_s(self); 
	c_k(kind); 
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::one()
#else
Computes the identity permutation.
The permutation in this must already be allocated. 
#endif
{
	INT N, l;
	
	N = s_li();
	for (l = 1; l <= N; l++) {
		m_ii(l - 1, l);
		}
	return(OK);
}

#if TEXDOCU
INT PERMUTATION_OB::print_list()
#else
Prints a permutation in list notation.
#endif
{
	INT i, j, l;
	
	l = s_li();
	for (i = 0; i < l; i++) {
		j = s_ii(i);
		printf("%ld ", j);
		}
	printf("\n");
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::sprint_list(BYTE *str)
#else
Prints a permutation in list notation into a string.
#endif
{
	INT i, j, l;
	
	l = s_li();
	sprintf(str + strlen(str), "[");
	for (i = 0; i < l; i++) {
		j = s_ii(i);
		sprintf(str + strlen(str), "%ld", j);
		if (i < l - 1)
			sprintf(str + strlen(str), ",");
		}
	sprintf(str + strlen(str), "]");
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::fprint_list(FILE *fp)
#else
Prints a permutation in list notation to FILE.
#endif
{
	INT i, j, l;
	
	l = s_li();
	for (i = 0; i < l; i++) {
		j = s_ii(i);
		fprintf(fp, "%ld ", j);
		}
	fprintf(fp, "\n");
	return OK;
}

INT f_perm_print_start_with_zero = FALSE;

INT PERMUTATION_OB::sprint(BYTE *str)
{
	BYTE have_seen[MAX_PERM_N + 1];
	INT l, l1, first, next, len, N, ii;
	INT f_nothing_printed_at_all = TRUE;
	
	if (str == NIL) {
		Srfs("PERM::sprint", "args NIL");
		return(ERROR);
		}
	N = s_li();
	if (N > MAX_PERM_N) {
		Srfs("PERM::sprint", "N > MAX_PERM_N");
		return(ERROR);
		}
	for (l = 1; l <= N; l++) {
		have_seen[l] = FALSE;
		}
	l = 1L;
	while (TRUE) {
		if (l > N) {
			if (f_nothing_printed_at_all) {
				sprintf(Eostr(str), "id");
				}
			return(OK);
			}
		if (have_seen[l]) {
			l++;
			continue;
			}
		/* Bearbeite Zyklus, beginnend mit l: */
		first = l;
		l1 = l;
		len = 1;
		while (TRUE) {
			have_seen[l1] = TRUE;
			next = s_ii(l1 - 1);
			if (next > MAX_PERM_N) {
				Srfs("PERM::sprint", 
				"next > MAX_PERM_N");
				printf("l1 = %ld next = %ld\n", 
				l1, next);
				return(ERROR);
				}
			if (next == first) {
				break;
				}
			if (have_seen[next]) {
				printf("%s\n", str);
				printf("(");
				for (ii = 0; ii < N; ii++) {
					printf("%2ld ", ii);
					}
				printf(")\n");
				printf("(");
				for (ii = 0; ii < N; ii++) {
					printf("%2ld ", s_ii(ii) - 1);
					}
				printf(")\n");
				fflush(stdout);
				Srfs("PERM::sprint", 
				"have_seen[next]");
				return(ERROR);
				}
			l1 = next;
			len++;
			}
		if (len == 1)
			continue;
		f_nothing_printed_at_all = FALSE;
		/* Drucke Zyklus, beginnend mit first: */
		l1 = first;
		sprintf(Eostr(str), "(");
		while (TRUE) {
			if (f_perm_print_start_with_zero)
				sprintf(Eostr(str), "%ld", l1 - 1);
			else
				sprintf(Eostr(str), "%ld", l1);
			next = s_ii(l1 - 1);
			if (next == first) {
				break;
				}
			sprintf(Eostr(str), " ");
			l1 = next;
			}
		sprintf(Eostr(str), ")");
		}
	return(OK);
}

INT PERMUTATION_OB::fprint_GAP(FILE *fp)
{
	BYTE str[10000];

	str[0] = 0;
	sprint_GAP(str);
	fprintf(fp, "%s", str);
	return OK;
}

INT PERMUTATION_OB::sprint_GAP(BYTE *str)
{
	BYTE have_seen[MAX_PERM_N + 1];
	INT l, l1, first, next, len, N, ii;
	INT f_nothing_printed_at_all = TRUE;
	
	if (str == NIL) {
		Srfs("PERM::sprint", "args NIL");
		return(ERROR);
		}
	N = s_li();
	if (N > MAX_PERM_N) {
		Srfs("PERM::sprint", "N > MAX_PERM_N");
		return(ERROR);
		}
	for (l = 1; l <= N; l++) {
		have_seen[l] = FALSE;
		}
	l = 1L;
	while (TRUE) {
		if (l > N) {
			if (f_nothing_printed_at_all) {
				sprintf(Eostr(str), "()");
				}
			return(OK);
			}
		if (have_seen[l]) {
			l++;
			continue;
			}
		/* Bearbeite Zyklus, beginnend mit l: */
		first = l;
		l1 = l;
		len = 1;
		while (TRUE) {
			have_seen[l1] = TRUE;
			next = s_ii(l1 - 1);
			if (next > MAX_PERM_N) {
				Srfs("PERM::sprint", 
				"next > MAX_PERM_N");
				printf("l1 = %ld next = %ld\n", 
				l1, next);
				return(ERROR);
				}
			if (next == first) {
				break;
				}
			if (have_seen[next]) {
				printf("%s\n", str);
				printf("(");
				for (ii = 0; ii < N; ii++) {
					printf("%2ld ", ii);
					}
				printf(")\n");
				printf("(");
				for (ii = 0; ii < N; ii++) {
					printf("%2ld", s_ii(ii) - 1);
					if (ii < N - 1)
						printf(", ");
					}
				printf(")\n");
				fflush(stdout);
				Srfs("PERM::sprint", 
				"have_seen[next]");
				return(ERROR);
				}
			l1 = next;
			len++;
			}
		if (len == 1)
			continue;
		f_nothing_printed_at_all = FALSE;
		/* Drucke Zyklus, beginnend mit first: */
		l1 = first;
		sprintf(Eostr(str), "(");
		while (TRUE) {
			if (f_perm_print_start_with_zero)
				sprintf(Eostr(str), "%ld", l1 - 1);
			else
				sprintf(Eostr(str), "%ld", l1);
			next = s_ii(l1 - 1);
			if (next == first) {
				break;
				}
			sprintf(Eostr(str), ", ");
			l1 = next;
			}
		sprintf(Eostr(str), ")");
		}
	return(OK);
}

INT PERMUTATION_OB::latex(FILE *fp)
{
	BYTE str[1024];

	str[0] = 0;
	sprint_latex(str);
	fprintf(fp, "%s", str);
	return OK;
}

INT PERMUTATION_OB::sprint_latex(BYTE *s)
{
	BYTE have_seen[MAX_PERM_N + 1];
	INT l, l1, first, next, len, N;
	INT f_nothing_printed_at_all = TRUE;
	
	N = s_li();
	if (N > MAX_PERM_N) {
		Srfs("PERM::sprint_latex", "N > MAX_PERM_N");
		return(ERROR);
		}
	for (l = 1; l <= N; l++) {
		have_seen[l] = FALSE;
		}
	l = 1L;
	while (TRUE) {
		if (l > N) {
			if (f_nothing_printed_at_all) {
				strcat(s, "id");
				}
			return(OK);
			}
		if (have_seen[l]) {
			l++;
			continue;
			}
		/* Bearbeite Zyklus, beginnend mit l: */
		first = l;
		l1 = l;
		len = 1;
		while (TRUE) {
			have_seen[l1] = TRUE;
			next = s_ii(l1 - 1);
			if (next > MAX_PERM_N) {
				Srfs("PERM::sprint_latex", 
				"next > MAX_PERM_N");
				return(ERROR);
				}
			if (next == first) {
				break;
				}
			if (have_seen[next]) {
				Srfs("PERM::sprint_latex", 
				"have_seen[next]");
				return(ERROR);
				}
			l1 = next;
			len++;
			}
		if (len == 1)
			continue;
		f_nothing_printed_at_all = FALSE;
		/* Drucke Zyklus, beginnend mit first: */
		l1 = first;
		if (strlen(s) > 200)
			return OK;
		sprintf(s + strlen(s), "(");
		while (TRUE) {
			sprintf(s + strlen(s), "%ld", l1);
			next = s_ii(l1 - 1);
			if (next == first) {
				break;
				}
			sprintf(s + strlen(s), " \\,");
			l1 = next;
			}
		sprintf(s + strlen(s), ")");
		}
	return(OK);
}

INT PERMUTATION_OB::Add2Cycle(INT i0, INT i1)
{
	INT N;
	
	N = s_li();
	if (i0 < 1L || i0 > N || i1 < 1L || i1 > N) {
		Srfs("PERM::Add2Cycle", "i? < 1 || i? > N");
		return(ERROR);
		}
	m_ii(i0 - 1, i1);
	m_ii(i1 - 1, i0);
	return(OK);
}

INT PERMUTATION_OB::Add3Cycle(INT i0, INT i1, INT i2)
{
	INT N;
	
	N = s_li();
	if (i0 < 1L || i0 > N || i1 < 1L || 
		i1 > N || i2 < 1L || i2 > N) {
		Srfs("PERM::Add3Cycle", "i? < 1 || i? > N");
		return(ERROR);
		}
	m_ii(i0 - 1, i1);
	m_ii(i1 - 1, i2);
	m_ii(i2 - 1, i0);
	return(OK);
}

INT PERMUTATION_OB::Add4Cycle(INT i0, INT i1, INT i2, 
	INT i3)
{
	INT N;
	
	N = s_li();
	if (i0 < 1L || i0 > N || i1 < 1L || 
		i1 > N || i2 < 1L || i2 > N || i3 < 1L || i3 > N) {
		Srfs("PERM::Add4Cycle", "i? < 1 || i? > N");
		return(ERROR);
		}
	m_ii(i0 - 1, i1);
	m_ii(i1 - 1, i2);
	m_ii(i2 - 1, i3);
	m_ii(i3 - 1, i0);
	return(OK);
}

INT PERMUTATION_OB::Add5Cycle(INT i0, INT i1, INT i2, 
	INT i3, INT i4)
{
	INT N;
	
	N = s_li();
	if (i0 < 1L || i0 > N || i1 < 1L || 
		i1 > N || i2 < 1L || i2 > N || 
		i3 < 1L || i3 > N || i4 < 1L || i4 > N) {
		Srfs("PERM::Add5Cycle", "i? < 1 || i? > N");
		return(ERROR);
		}
	m_ii(i0 - 1, i1);
	m_ii(i1 - 1, i2);
	m_ii(i2 - 1, i3);
	m_ii(i3 - 1, i4);
	m_ii(i4 - 1, i0);
	return(OK);
}

INT PERMUTATION_OB::Add6Cycle(INT i0, INT i1, INT i2, 
	INT i3, INT i4, INT i5)
{
	INT N;
	
	N = s_li();
	if (i0 < 1L || i0 > N || i1 < 1L || 
		i1 > N || i2 < 1L || i2 > N || 
		i3 < 1L || i3 > N || i4 < 1L || 
		i4 > N || i5 < 1L || i5 > N) {
		Srfs("PERM::Add6Cycle", "i? < 1 || i? > N");
		return(ERROR);
		}
	m_ii(i0 - 1, i1);
	m_ii(i1 - 1, i2);
	m_ii(i2 - 1, i3);
	m_ii(i3 - 1, i4);
	m_ii(i4 - 1, i5);
	m_ii(i5 - 1, i0);
	return(OK);
}

INT PERMUTATION_OB::Add7Cycle(INT i0, INT i1, INT i2, 
	INT i3, INT i4, INT i5, INT i6)
{
	INT N;
	
	N = s_li();
	if (i0 < 1L || i0 > N || i1 < 1L || 
		i1 > N || i2 < 1L || i2 > N || 
		i3 < 1L || i3 > N || i4 < 1L || 
		i4 > N || i5 < 1L || i5 > N || 
		i6 < 1L || i6 > N) {
		Srfs("PERM::Add7Cycle", "i? < 1 || i? > N");
		return(ERROR);
		}
	m_ii(i0 - 1, i1);
	m_ii(i1 - 1, i2);
	m_ii(i2 - 1, i3);
	m_ii(i3 - 1, i4);
	m_ii(i4 - 1, i5);
	m_ii(i5 - 1, i6);
	m_ii(i6 - 1, i0);
	return(OK);
}

INT PERMUTATION_OB::Add8Cycle(INT i0, INT i1, INT i2, 
	INT i3, INT i4, INT i5, 
	INT i6, INT i7)
{
	INT N;
	
	N = s_li();
	if (i0 < 1L || i0 > N || i1 < 1L || 
		i1 > N || i2 < 1L || i2 > N || 
		i3 < 1L || i3 > N || i4 < 1L || 
		i4 > N || i5 < 1L || i5 > N || 
		i6 < 1L || i6 > N || i7 < 1L || i7 > N) {
		Srfs("PERM::Add8Cycle", "i? < 1 || i? > N");
		return(ERROR);
		}
	m_ii(i0 - 1, i1);
	m_ii(i1 - 1, i2);
	m_ii(i2 - 1, i3);
	m_ii(i3 - 1, i4);
	m_ii(i4 - 1, i5);
	m_ii(i5 - 1, i6);
	m_ii(i6 - 1, i7);
	m_ii(i7 - 1, i0);
	return(OK);
}

INT PERMUTATION_OB::AddNCycle(INT first, INT len)
{
	INT N, i;
	
	N = s_li();
	if (first < 1 || first + len - 1 > N) {
		Srfs("PERM::AddNCycle", 
		"first < 1 || first + len - 1 > N");
		return(ERROR);
		}
	for (i = 0; i < len; i++) {
		if (i == len - 1)
			m_ii(first + i - 1, first);
		else
			m_ii(first + i - 1, first + i + 1);
		}
	return(OK);
}

INT PERMUTATION_OB::AddNCycleOffset(INT first, INT len, INT offset)
/* len kann auch 1 sein ! */
{
	INT j;
	
	for (j = 0; j < len; j++) {
		if (j < len - 1)
			m_ii(first + j * offset - 1, 
				first + (j + 1) * offset);
		else
			m_ii(first + j * offset - 1, first);
		}
	return(OK);
}

INT PERMUTATION_OB::is_full_cycle()
{
	INT i0, i1, i, l;

	l = s_li();	
	i0 = 0;
	for (i = 0; i < l; i++) {
		i1 = s_ii(i0) - 1;
		if (i1 == 0)
			break;
		i0 = i1;
		}
	if (i == l - 1)
		return TRUE;
	else
		return FALSE;
}

INT PERMUTATION_OB::order(INT *order)
/* Rueckgabe: order == 1 falls p == id. 
 * Sonst order = Elementordnung(p). */
{
	INT N, order1, g;
	BYTE have_seen[MAX_PERM_N + 1L];
	INT l, l1, first, next, len;
	
	if (order == NIL) {
		Srfs("PERM::order", "args NIL");
		return(ERROR);
		}
	N = s_li();
	if (N > MAX_PERM_N) {
		Srfs("PERM::order", "N > MAX_PERM_N");
		return(ERROR);
		}
	for (l = 1; l <= N; l++) {
		have_seen[l] = FALSE;
		}
	order1 = 1;
	l = 1L;
	while (TRUE) {
		if (l > N) {
			*order = order1;
			return(OK);
			}
		if (have_seen[l]) {
			l++;
			continue;
			}
		/* Bearbeite Zyklus, beginnend mit l: */
		first = l;
		l1 = l;
		len = 1;
		while (TRUE) {
			have_seen[l1] = TRUE;
			next = s_ii(l1 - 1);
			if (next >= MAX_PERM_N) {
				Srfs("PERM::order", 
				"next >= MAX_PERM_N");
				return(ERROR);
				}
			if (next == first) {
				break;
				}
			if (have_seen[next]) {
				Srfs("PERM::order", 
				"have_seen[next]");
				return(ERROR);
				}
			l1 = next;
			len++;
			}
		/* found a cycle; length in len */
		if (len == 1)
			continue;
		/* order1 = kgv(order1, len) */
		if (ggt_iipi(order1, len, &g) != OK) {
			Srff("PERM::order", "ggt_iipi");
			return(ERROR);
			}
		len /= g;
		order1 *= len;
		}
	return(OK);
}

INT PERMUTATION_OB::order_if_prime(INT *order)
/* Rueckgabe: order == 0, 
 * wenn Elementordnung (p) nicht prim und nicht eins;
 * sonst in order die entsprechende Primzahl 
 * bzw eins (wenn p == id). */
{
	INT N, order1, prime;
	BYTE have_seen[MAX_PERM_N + 1L];
	INT l, l1, first, next, len;
	
	if (order == NIL) {
		Srfs("PERM::order_if_prime", "args NIL");
		return(ERROR);
		}
	N = s_li();
	if (N > MAX_PERM_N) {
		Srfs("p_OrderIfPrime", "N > MAX_PERM_N");
		return(ERROR);
		}
	for (l = 1; l <= N; l++) {
		have_seen[l] = FALSE;
		}
	order1 = 1;
	l = 1L;
	while (TRUE) {
		if (l > N) {
			*order = order1;
			return(OK);
			}
		if (have_seen[l]) {
			l++;
			continue;
			}
		/* Bearbeite Zyklus, beginnend mit l: */
		first = l;
		l1 = l;
		len = 1;
		while (TRUE) {
			have_seen[l1] = TRUE;
			next = s_ii(l1 - 1);
			if (next >= MAX_PERM_N) {
				Srfs("PERM::order_if_prime", 
				"next >= MAX_PERM_N");
				return(ERROR);
				}
			if (next == first) {
				break;
				}
			if (have_seen[next]) {
				Srfs("PERM::order_if_prime", 
				"have_seen[next]");
				return(ERROR);
				}
			l1 = next;
			len++;
			}
		/* found a cycle; length in len */
		if (len == 1)
			continue;
		prime = smallest_primedivisor(len);
		if (len != prime) {
			/* Zykellaenge nicht prim -> Abbruch. */
			*order = 0L;
			return(OK);
			}
		if (order1 != 1 && prime != order1) {
			/* Verschiedene Primzahlen 
			 * als Zykellaengen: 
			 * zusammengesetzte Elementordnung 
			 * -> Abbruch */
			*order = 0L;
			return(OK);
			}
		else {
			/* entweder: order1 != 1L && 
			 *    prime_of_order1 == order1
			 * oder:     order1 == 1
			 * in beiden Faellen: */
			order1 = prime;
			}
		}
}

INT PERMUTATION_OB::order_if_prime_power(
	INT *order, INT *prime, INT *k)
/* Rueckgabe: order == 0, 
 * wenn Elementordnung keine Primpotenz;
 * sonst in order die entsprechende 
 * Primpotenz (order == prime ^ k). 
 * order == 1 falls p == id 
 * (und prime == 1 und k == 0). */
{
	INT N, order1, prime1, k1, prime2, prime3, k2;
	/* order1 = prime1 ^ k1 */
	BYTE have_seen[MAX_PERM_N + 1L];
	INT l, l1, first, next, len, len1;
	
	if (order == NIL || prime == NIL || k == NIL) {
		Srfs("PERM::order_if_prime_power", "args NIL");
		return(ERROR);
		}
	N = s_li();
	if (N > MAX_PERM_N) {
		Srfs("PERM::order_if_prime_power", 
		"N > MAX_PERM_N");
		return(ERROR);
		}
	for (l = 1; l <= N; l++) {
		have_seen[l] = FALSE;
		}
	order1 = 1;
	prime1 = 1;
	k1 = 0;
	l = 1;
	while (TRUE) {
		if (l > N) {
			*order = order1;
			*prime = prime1;
			*k = k1;
			return(OK);
			}
		if (have_seen[l]) {
			l++;
			continue;
			}
		/* Bearbeite Zyklus, beginnend mit l: */
		first = l;
		l1 = l;
		len = 1;
		while (TRUE) {
			have_seen[l1] = TRUE;
			next = s_ii(l1 - 1);
			if (next >= MAX_PERM_N) {
				Srfs("PERM::order_if_prime_power", 
				"next >= MAX_PERM_N");
				return(ERROR);
				}
			if (next == first) {
				break;
				}
			if (have_seen[next]) {
				Srfs("PERM::order_if_prime_power", 
				"have_seen[next]");
				return(ERROR);
				}
			l1 = next;
			len++;
			}
		/* found a cycle; length in len */
		if (len == 1)
			continue;
		len1 = len;
		prime2 = smallest_primedivisor(len);
		if (prime1 != 1 && 
			prime1 != prime2) {
			/* Es taucht eine andere Primzahl 
			 * als Zykellaenge auf -> Abbruch */
			*order = 0;
			return(OK);
			}
		k2 = 1;
		len /= prime2;
		while (len != 1) {
			prime3 = smallest_primedivisor(len);
			if (prime3 != prime2) {
				/* Zykellaenge nicht von 
				 * Primpotenz -> Abbruch */
				*order = 0;
				return(OK);
				}
			len /= prime2;
			k2++;
			}
		/* Es wurde len zerlegt in: len1 == prime2 ^ k2 */
		if (prime1 == 1) { /* erster nichttrivialer Zyklus: */
			order1 = len1;
			prime1 = prime2;
			k1 = k2;
			}
		else { /* order1 = kgv(order1, len1): */
			k1 = MAX(k1, k2);
			order1 = MAX(order1, len1);
			}
		}
	return(OK);
}

INT PERMUTATION_OB::compare(PERMUTATION_OP b)
{
	return s_s()->sym_comp(b->s_s());
}


INT PERMUTATION_OB::invers(PERMUTATION_OP b)
/* before: invers_permutation(OP perm, OP b) */
{
	INT i, erg = OK;
	VECTOR_OP self;

	if (s_k() != VECTOR)
		return error("PERM::invers: wrong perm type");
	if (!b->emptyp()) 
		erg += ((SYM_OP) b)->freeself();

	self = (VECTOR_OP) callocobject("PERM::invers: s");
	erg += self->m_il(s_li());
	for (i = 0; i < s_li(); i++)
		self->m_ii(s_ii(i) - 1, i + 1);
	erg += b->b_ks(VECTOR, self);
	b->c_obj_k(s_obj_k());
	return erg;
}

INT PERMUTATION_OB::mult(PERMUTATION_OP b, PERMUTATION_OP c)
/* before: mult_permutation(OP a, OP b, OP c) */
/* before: c := a(b(i)) - (erst b, dann a). */
/* now: first this, then b into c */
{
	INT i;
	INT erg = OK;
	/* if (b->nullp())
		return m_i_i(0L,c); */
	if (b->s_obj_k() != PERMUTATION)
		return error("PERM::mult: wrong second type");
	if ((s_k() != VECTOR) || (b->s_k() != VECTOR))
		return error("PERM::mult: only for VECTOR type");
	erg += copy(c);
	for (i = 0; i < c->s_li(); i++)
		c->m_ii(i, b->s_ii(s_ii(i) - 1));
	return erg;
}

INT PERMUTATION_OB::copy(PERMUTATION_OP b)
{
	INT erg;
	
	((SYM_OP) b)->freeself();
	erg = b->b_ks(s_k(), (VECTOR_OP)callocobject("PERM::copy s"));
	erg += s_s()->copy(b->s_s());
	b->c_obj_k(s_obj_k());
	if (erg != OK)
		return erg;
	return OK;
}

INT PERMUTATION_OB::einsp()
{
	INT i;
	
	for (i = s_li() - 1; i >= 0; i--) {
		if (s_ii(i) != (i + 1))
			return(FALSE);
		}
	return(TRUE);
}

#if TEXDOCU
INT VECTOR_OB::first_lehmercode(INTEGER_OP l)
#else
Computes the first lehmercode for a permutation of given degree $l$.
This is the vector $0,\ldots,0$.
#endif
{
	INT i;
	
	m_il(l->s_i());
	for (i = 0; i < s_li(); i++)
		m_ii(i, 0);
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::last_lehmercode(INTEGER_OP l)
#else
Computes the last lehmercode for a permutation of given degree $l$.
This is the vector $0,1,\ldots,n-1$.
#endif
{
	INT i, j = l->s_i() - 1;

	m_il(l->s_i());
	for (i = 0; i < l->s_i(); i++, j--)
		m_ii(i, j);
	return(OK);
}

#if TEXDOCU
INT VECTOR_OB::next_lehmercode(VECTOR_OP next)
#else
Computes the lexicographically next lehmercode.
Returns LASTLEHMERCODE if there is no next lehmercode.
#endif
/* Erzeugt den lexikographisch naechsten Lehmercode nach next. 
 * Rueckgabe LASTLEHMERCODE, falls Ende erreicht 
 * (next ist dann EMPTY). */
{
	INT i, j;
	
	copy(next);
	for (i = next->s_li() - 1, j = 0; i >= 0; i--, j++) {
		if (next->s_ii(i) < j)
			return(next->s_i(i)->inc());
		else
			next->c_ii(i, 0);
		}
	next->freeself();
	return(LASTLEHMERCODE);
}

#if TEXDOCU
INT VECTOR_OB::lehmercode_perm(PERMUTATION_OP b)
#else
Computes the permutation $b$ defined by its lehmercode (this).
#endif
/* before: lehmercode_vector(OP vec, OP b) */
/* Berechnet aus dem Lehmercode vec = [v1,....,vn]
 * die zugehoerige Permutation b [e1,...,en] */
{
	INT i, j, k;
	VECTOR_OP self, liste;
	
	self = (VECTOR_OP) callocobject("PERM::lehmercode_perm self");
	liste = (VECTOR_OP) callocobject("PERM::lehmercode_perm liste");

	self->m_il(s_li());
	liste->m_il(s_li());
	/* initialisierung zweier vektoren fuer
	 * eine Liste und fuer die zu berechnende Permutation */
	for (i = 0; i < liste->s_li(); i++)
		liste->m_ii(i, i + 1);
	/* liste ist jetzt ein vector [1,2,3,....,n] */
	
	for (i = 0; i < s_li(); i++) {
		k = s_ii(i);
		/* k ist ist das i-te Element aus vec, also vi */
		self->m_ii(i, liste->s_ii(k));
		/* daher ist ei = k-te Element aus der aktuellen Liste*/
		
		for (j = k; j < (s_li() - 1) - i; j++)
			liste->c_ii(j, liste->s_ii(j + 1));
			/* in der liste wird das k-te Element gestrichen.
			und von rechts aufgefuellt */
		}
	freeall(liste);

	b->b_ks(VECTOR, self);
	/* bildung einer Permutation aus dem vector */
	return(OK);
}

#if TEXDOCU
INT PERMUTATION_OB::first_permutation(INTEGER_OP l)
#else
Computes the first permutation of degree $l$.
#endif
/* l bleibt erhalten */
{
	INT erg = OK;
	
	erg += m_il(l->s_i());
	one();
#if FALSE
	for (i = 0; i < s_li(); i++)
		m_ii(i, i + 1);
#endif
	return erg;
}

#if TEXDOCU
INT PERMUTATION_OB::next_permutation_lex(PERMUTATION_OP next_perm)
#else
Computes the next permutation using an elementary algorithm.
#endif
/* Fischer Krause */
{
	INT r, s, i, j, erg;
	PERMUTATION_OP c;
	
	if (this == next_perm) {
		c = (PERMUTATION_OP) callocobject("PERM::next_permutation_lex c");
		*c = *this;
		next_perm->c_obj_k(EMPTY);
		erg = c->next_permutation_lex(next_perm);
		freeall(c);
		return erg;
		}
	copy(next_perm);
	for (r = next_perm->s_li() - 2; r >= 0; r--)
		if (next_perm->s_ii(r) < next_perm->s_ii(r + 1))
			break;
	if (r == -1)
		return LASTPERMUTATION;
	for (s = 0; s < next_perm->s_li() - r - 1; s++)
		if (next_perm->s_ii(r) > next_perm->s_ii(r + s + 1))
			break;
	next_perm->s_i(r)->swap(next_perm->s_i(r + s));
	for (i = r + 1, j = next_perm->s_li() - 1; i < j; i++, j--)
		next_perm->s_i(i)->swap(next_perm->s_i(j));
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::perm_lehmercode(VECTOR_OP vec)
#else
Computes the lehmercode of a permutation.
#endif
/* before: lehmercode_permutation(
 *   OP perm, OP vec) */
/* Berechnet zur Permutation perm = [p1,....,pn]
 *  den zugehoerigen Lehmercode vec [v1,...,vn] */
{
	INT i,j,k;
	INT erg = OK;
	
	if (s_k() == ZYKEL) {
		 /* erg += t_ZYKEL_VECTOR(perm); !!! */
		erg += perm_lehmercode(vec);
		return erg;
		}

	erg += vec->m_il(s_li());
	/* erzeugt ein Vectorobject */
	for (i = 0; i < s_li(); i++) {
		k = 0;
		for (j = i + 1; j < s_li(); j++)
			if (s_ii(j) < s_ii(i))
				k++;
		/* k ist die Anzahl der Permutationselemente
		 * rechts von pi, die kleiner sind */
		vec->m_ii(i, k);
		/* k wird an der richtigen Stelle im
				Vector notiert */
		}
	return erg;
}

#if TEXDOCU
INT PERMUTATION_OB::perm_lehmercode2(VECTOR_OP vec)
#else
#endif
/* before: lehmercode2_permutation(
 *   OP perm, OP vec) */
{
	INT i, j, k;
	
	s_s()->copy(vec);
	for (i = 0; i < vec->s_li(); ) {
		k = vec->s_ii(i) - 1;
		vec->m_ii(i, k);
		i++;
		for (j = i; j < vec->s_li(); j++)
			if (vec->s_ii(j) > k)
				vec->s_i(j)->dec();
		}
	return(OK);
}

#if TEXDOCU
INT PERMUTATION_OB::next_permutation(PERMUTATION_OP next)
#else
Computes the next permutation using lehmercodes.
#endif
/* before: next_permutation(OP start, OP n) */
/* Erzeuge naechste Permutation via Lehmercode nach n. 
 * Rueckgabe LASTPERMUTATION, wenn Ende. n ist dann 
 * unveraendert. */
{
	VECTOR_OB zwa, zwb;
	
	perm_lehmercode(&zwa);
	if (zwa.next_lehmercode(&zwb) == LASTLEHMERCODE) {
		return(LASTPERMUTATION);
		}
	zwb.lehmercode_perm(next);
	return(OK);
}

#if TEXDOCU
INT PERMUTATION_OB::signum()
#else
+1 if the permutation is even, -1 if it is odd.
#endif
/* before: signum_permutation(OP perm, OP b) */
{
	VECTOR_OB zwischen;
	INT i, summe = 0;
	
	perm_lehmercode(&zwischen);
	for (i = 0; i < zwischen.s_li(); i++)
		summe += zwischen.s_ii(i);
	if (EVEN(summe))
		return (1);
	else
		return (-1);
}

#if TEXDOCU
INT PERMUTATION_OB::numberof_inversions()
#else
Computes the number of inversions of the permutation. 
The number is computed as the sum of the entries in the lehmercode.
#endif
/* before: numberof_inversionen(OP a, OP b) */
/* Rueckgabe ist die Anzahl der 
 * Inversionen in der permutation a */
{
	INT summe = 0, i;
	VECTOR_OB vec;

	if (s_obj_k() != PERMUTATION)
		return error("PERM::numberof_inversions: wrong type");
	if ((s_k() != VECTOR) && (s_k() != ZYKEL))
		return error("PERM::numberof_inversions: wrong perm type");

	perm_lehmercode(&vec);
	for (i = 0; i < vec.s_li(); i++)
		summe += vec.s_ii(i);
	return summe;
}

#if TEXDOCU
INT PERMUTATION_OB::rz(VECTOR_OP c)
#else
reduced decomposition.
#endif
{
	INT erg = OK;
	VECTOR_OP lc = (VECTOR_OP) callocobject("PERM::rz lc");

	erg += perm_lehmercode(lc);
	erg += lc->rz_lehmercode(c);
	erg += freeall(lc);
	return erg;
}

#if TEXDOCU
INT VECTOR_OB::rz_lehmercode(VECTOR_OP b)
#else
Reduzierte Zerlegung des Lehmercodes. 
bsp 321200 liefert 32132354
#endif
{
	INT i = s_li();
	INT k, j, erg = OK;
	
	INTEGER_OP zw = (INTEGER_OP) callocobject("PERM::rz_lehmercode zw");

	if (b == NULL)
		return error("rz_lehmercode: b == NULL");

	erg += sum(zw); 
	if (zw->nullp()) {
		erg += b->m_il(0);
		erg += freeall(zw);
		goto ende;
		}
	k = zw->s_i();
	erg += b->m_l(zw);
	/* die Laenge der reduzierten Zerlegung 
	 * ist die summe des lehmercodes */
	while (i-- > 0) {
		if (s_ii(i) > 0) {
			for (j = 0; j < s_ii(i); j++) {
				k--;
				if (k < 0)
					return error(
					"rz_lehmercode: k < 0");
				b->m_ii(k, i + 1 + j);
				}
			}
		}
	freeall(zw);
ende:
	return erg;
}

#define MAX_PERM_LEN 10000

#if TEXDOCU
INT PERMUTATION_OB::sscan(BYTE **str, INT f_v)
#else
Scans a permutation from a string.
#endif
{
	BYTE *p = *str;
	INT perm[MAX_PERM_LEN];
	INT cycle[MAX_PERM_LEN];
	INT i, a_last, a, ma, dig, ci;
	BYTE s[256], c;
	INT si;

	for (i = 0; i < MAX_PERM_LEN; i++)
		perm[i] = i;
	ma = 1;
	m_il(1);
	one();
	while (TRUE) {
		ci = 0;
		if (*p != '(') {
			*str = p;
			break;
			}
		if (f_v) {
			printf("opening parenthesis\n"); fflush(stdout);
			}
		p++;
		while (TRUE) {
			while (*p == ' ' || *p == '\t')
				p++;
			si = 0;
			while (*p == ' ' || *p == '\t')
				p++;

			// read first digit:
			c = *p;
			if (c >= '0' && c <= '9') {
				if (f_v) {
					printf("character %c\n", c); fflush(stdout);
					}
				s[si++] = c;
				p++;
				}
			else {
				*p = 0;
				printf("sscan(): character unexpected in %s", *str);
				fflush(stdout);
				if (c != 0)
					printf("%s", p + 1);
				printf("\n");
				fflush(stdout);
				*p = c;
				return ERROR;
				}
			
			// read following digits:
			while (TRUE) {
				c = *p;
				if (c >= '0' && c <= '9') {
					if (f_v) {
						printf("character %c\n", c); fflush(stdout);
						}
					s[si++] = c;
					p++;
					}
				else {
					break;
					}
				} // reading digits
			
			while (*p == ' ' || *p == '\t')
				p++;
			if (*p == ',')
				p++;
			while (*p == ' ' || *p == '\t')
				p++;
#if 0
			if (*p != ')' && *p != ',') {
				c = *p;
				*p = 0;
				printf("sscan(): character unexpected in %s", *str);
				if (c != 0)
					printf("%s\n", p + 1);
				*p = c;
				return ERROR;
				}
#endif
			s[si] = 0;
			sscanf(s, "%ld", &dig);
			if (f_v) {
				printf("digit %ld\n", dig);
				fflush(stdout);
				}
			si = 0;
			if (dig == 0) {
				printf("permutation permutes the elements 1,2,... (0 not allowed)\n");
				return ERROR;
				}
			cycle[ci++] = dig;
			if (*p == ')') {
				if (f_v) {
					printf("closing parenthesis, cycle: ");
					for (i = 0; i < ci; i++)
						printf("%ld ", cycle[i]);
					printf("\n");
					fflush(stdout);
					}
				for (i = 0; i < ci; i++)
					ma = MAXIMUM(ma, cycle[i]);
				for (i = 1; i < ci; i++) {
					a_last = cycle[i - 1];
					a = cycle[i];
					perm[a_last - 1] = a - 1;
					}
				if (ci > 1) {
					a_last = cycle[ci - 1];
					a = cycle[0];
					perm[a_last - 1] = a - 1;
					}
				if (f_v) {
					printf("perm (max.deg = %ld): ", ma);
					for (i = 0; i < ma; i++)
						printf("%ld ", perm[i] + 1);
					printf("\n");
					fflush(stdout);
					}
				ci = 0;
				p++;
				break;
				}
#if 0
			if (f_v) {
				printf("reading ','\n");
				}
			p++;
#endif
			} // loop for one cycle
		while (*p == ' ' || *p == '\t')
			p++;
		ci = 0;
		} // end of loop over all cycles
	if (f_v) {
		printf("get perm\n");
		fflush(stdout);
		}
	m_il(ma);
	for (i = 1; i <= ma; i++) {
		a = perm[i - 1] + 1;
		m_ii(i - 1, a);
		}
	if (f_v) {
		println();
		}
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::cycle_type(VECTOR_OP v)
#else
Computes the cycle type of the permutation into the vector $v$.
$v[i]$ contains the number of $i+1$ cycles of the permutation.
#endif
{
	INT i, j, n, next, l;
	VECTOR_OB v1;
	
	n = s_li();
	v->m_il_n(n);
	v1.m_il_n(n);
	for (i = 0; i < n; i++) {
		if (v1.s_ii(i) == 1)
			continue;
		v1.m_ii(i, 1);
		j = i;
		l = 1;
		while (TRUE) {
			next = s_ii(j) - 1;
			if (v1.s_ii(next) == 1)
				break;
			v1.m_ii(next, 1);
			j = next;
			l++;
			}
		v->s_i(l - 1)->inc();
		}
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::embed(INT n)
#else
Embeds the permutation into $S_n$.
#endif
{
	INT i, j, l;
	PERMUTATION_OB p1;
	
	l = s_li();
	if (n < l) {
		printf("PERM::embed: %ld = n < l = %ld\n", n, l);
		return error("PERM::embed");
		}
	if (n == l)
		return OK;
	p1.m_il(n);
	for (i = 0; i < l; i++) {
		j = s_ii(i);
		p1.m_ii(i, j);
		}
	for (i = l; i < n; i++) {
		p1.m_ii(i, i + 1);
		}
	swap(&p1);
	return(OK);
}

#if TEXDOCU
INT PERMUTATION_OB::conjugate(PERMUTATION_OP kappa, PERMUTATION_OP p_kappa)
#else
#endif
{
	PERMUTATION_OB kappa_inv, a, b;
	
	kappa->invers(&kappa_inv);
	kappa_inv.mult(this, &a);
	a.mult(kappa, &b);
	b.swap(p_kappa);
	return OK;
}




#endif /* PERMTRUE */


