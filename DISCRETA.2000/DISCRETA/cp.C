/* cp.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef CP_TRUE

#include <DISCRETA/cp.h>

void sp_nil(SPERM *p)
{
	p->l = 0;
	p->a = NIL;
}

INT sp_int(SPERM *p, INT l, BYTE *where)
{
	sp_free(p);
	p->a = (UINT2 *) my_malloc(l * sizeof(UINT2), where);
	if (p->a == NIL) {
		Srfs("sp_int", "no memory");
		return ERROR;
		}
	p->l = l;
	return OK;
}

INT sp_free(SPERM *p)
{
	if (p->a) {
		my_free(p->a);
		p->a = NIL;
		}
	p->l = 0;
	return OK;
}

void sp_free_it(SPERM *p)
{
	sp_free(p);
	my_free(p);
}

INT sp_mv(SPERM *p, SPERM *q)
{
	INT i;
	
	if (p->l != q->l) {
		Srfs("sp_mv", "p->l != q->l");
		return ERROR;
		}
	for (i = 0; i < p->l; i++)
		q->a[i] = p->a[i];
	return OK;
}

INT sp_id(SPERM *p)
{
	INT i;
	
	for (i = 0; i < p->l; i++)
		p->a[i] = i;
	return OK;
}

INT sp_mult(SPERM *a, SPERM *b, SPERM *c)
/* erst a, dann b; Ergebnis nach c */
{
	INT i, j;
	
	if (a->l != b->l) {
		Srfs("sp_mult", "a->l != b->l");
		return ERROR;
		}
	if (c->l != a->l) {
		Srfs("sp_mult", "c->l != a->l");
		return ERROR;
		}
	/* sp_int(c, a->l, "sp_int sp_mult"); */
	for (i = 0; i < a->l; i++) {
		j = a->a[i];
		c->a[i] = b->a[j];
		}
	return OK;
}

INT sp_inv(SPERM *a, SPERM *b)
/* b:= a^-1 */
{
	INT i, j;
	
	if (b->l != a->l) {
		Srfs("sp_inv", "b->l != a->l");
		return ERROR;
		}
	/* sp_int(b, a->l, "sp_int sp_inv"); */
	for (i = 0; i < a->l; i++) {
		j = a->a[i];
		b->a[j] = i;
		}
	return OK;
}

INT sp_power(SPERM *a, SPERM *res, INT exp)
{
	SPERM *b = NIL;
	SPERM *c = NIL;
	SPERM *d = NIL;
	INT erg = ERROR;
	INT len;
	
	len = a->l;
	sp_int(res, len, "sp_int sp_power");
	b = (SPERM *) my_malloc(sizeof(SPERM), "sp_power");
	c = (SPERM *) my_malloc(sizeof(SPERM), "sp_power");
	d = (SPERM *) my_malloc(sizeof(SPERM), "sp_power");
	sp_nil(b);
	sp_nil(c);
	sp_nil(d);
	sp_int(b, len, "sp_int sp_power");
	sp_int(c, len, "sp_int sp_power");
	sp_int(d, len, "sp_int sp_power");
	
	sp_id(c);
	sp_mv(a, b);
	while (exp > 0) {
		/* res = b^exp * c */
		if (ODD(exp)) {
			erg += sp_mult(b, c, d);
			erg += sp_mv(d, c);
			exp--;
			continue; /* exp == 0 possible */
			}
		if (EVEN(exp)) {
			erg += sp_mult(b, b, d);
			erg += sp_mv(d, b);
			exp >>= 1;
			}
		}
	erg += sp_mv(c, res);

	erg = OK;
	if (b) {
		sp_free_it(b);
		b = NIL;
		}
	if (c) {
		sp_free_it(c);
		c = NIL;
		}
	if (d) {
		sp_free_it(d);
		d = NIL;
		}
	return erg;
}

INT sp_mult_apply_tau_r(SPERM *a, INT i, INT j)
/* a := a (i j). */
{
	INT k, i1, j1;
	
	for (k = 0; k < a->l; k++) {
		if (a->a[k] == i)
			i1 = k;
		if (a->a[k] == j)
			j1 = k;
		}
	/* now: i1 a == i, j1 a == j */
	a->a[i1] = j;
	a->a[j1] = i;
	return OK;
}

INT sp_mult_apply_tau_l(SPERM *a, INT i, INT j)
/* a := (i j) a. */
{
	INT i1, j1;
	
	i1 = a->a[i];
	j1 = a->a[j];
	/* now: i -> j -> j1; j -> i -> i1 */
	a->a[i] = j1;
	a->a[j] = i1;
	return OK;
}

INT sp_mult_apply_forwc_r(SPERM *a, INT i, INT l)
/* a := a (i i+1 ... i+l-1). */
{
	INT t[SPERM_MAX_N], m, j, k;
	
	if (l > SPERM_MAX_N) {
		Srfs("sp_mult_apply_forwc_r", "l > SPERM_MAX_N");
		return ERROR;
		}
	for (m = 0; m < a->l; m++) {
		j = a->a[m];
		if (j >= i && (k = j - i) < l)
			t[k] = m;
		}
	/* now: t[k] -> i+k -> i+k+1 for 0 < k < l-1
	 *      t[l-1] -> i+l-1 -> i */
	for (k = 0; k < l - 1; k++)
		a->a[t[k]] = (i + k + 1);
	a->a[t[l - 1]] = i;
	return OK;
}

INT sp_mult_apply_backwc_r(SPERM *a, INT i, INT l)
/* a := a (i+l-1 i+l-2 ... i+1 i). */
{
	INT t[SPERM_MAX_N], m, j, k;
	
	if (l > SPERM_MAX_N) {
		Srfs("sp_mult_apply_backwc_r", "l > SPERM_MAX_N");
		return ERROR;
		}
	for (m = 0; m < a->l; m++) {
		j = a->a[m];
		if (j >= i && (k = j - i) < l)
			t[k] = m;
		}
	/* now: t[k] -> i+k -> i+k-1 for 1 < k < l
	 *      t[0] -> i -> i+l-1 */
	for (k = 1; k < l; k++)
		a->a[t[k]] = (i + k - 1);
	a->a[t[0]] = (i + l - 1);
	return OK;
}

INT sp_mult_apply_forwc_l(SPERM *a, INT i, INT l)
/* a := (i i+1 ... i+l-1) a. */
{
	INT t, m;
	
	/* i+m -> i+m+1 -> a[i+m+1]  for m < l - 1 
	 * i+l-1 -> i -> a[i] */
	t = a->a[i];
	for (m = 0; m < l - 1; m++)
		a->a[i + m] = a->a[i + m + 1];
	a->a[i+l-1] = t;
	return OK;
}

INT sp_mult_apply_backwc_l(SPERM *a, INT i, INT l)
/* a := (i+l-1 i+l-2 ... i+1 i) a. */
{
	INT t, m;
	
	/* i+m -> i+m-1 -> a[i+m-1]  for 1 < m < l 
	 * i -> i+l-1 -> a[i+l-1] */
	t = a->a[i + l - 1];
	for (m = l - 1; m > 0; m--)
		a->a[i + m] = a->a[i + m - 1];
	a->a[i] = t;
	return OK;
}

INT sp_onep(SPERM *p)
{
	INT i;
	
	for (i = 0; i < p->l; i++)
		if (p->a[i] != i)
			return FALSE;
	return TRUE;
}

INT sp_cmp(SPERM *a, SPERM *b)
{
	INT i;
	
	if (a->l != b->l) {
		Srfs("sp_cmp", "a->l != b->l");
		return 0;
		}
	for (i = 0; i < a->l; i++) {
		if (a->a[i] < b->a[i])
			return -1;
		if (a->a[i] > b->a[i])
			return 1;
		}
	return 0;
}

void sp_print(SPERM *p)
{
	BYTE s[10000];
	
	s[0] = 0;
	sp_sprint(p, s);
	printf("%s\n", s);
}

INT sp_sprint(SPERM *p, BYTE *s)
/* haengt an s an. */
{
	INT *have_seen = NIL;
	INT l, l1, first, next, len;
	INT f_nothing_printed_at_all = TRUE;
	BYTE str1[256];
	BYTE str2[256];
	
	if (p == NIL || s == NIL) {
		Srfs("sp_sprint", "args NIL");
		return ERROR;
		}
	str1[0] = 0;
	have_seen = (INT *) my_malloc(p->l * sizeof(INT), "sp_sprint");
	if (have_seen == NIL) {
		Srfs("sp_sprint", "no memory");
		return ERROR;
		}
	for (l = 0; l < p->l; l++) {
		have_seen[l] = FALSE;
		}
	l = 0;
	while (TRUE) {
		if (l >= p->l) {
			if (f_nothing_printed_at_all) {
				strcat(s, "id");
				}
			else {
				strcat(s, str1);
				}
			if (have_seen)
				my_free(have_seen);
			return OK;
			}
		if (have_seen[l]) {
			l++;
			continue;
			}
		/* cycle starting at l: */
		first = l;
		l1 = l;
		len = 1;
		while (TRUE) {
			have_seen[l1] = TRUE;
			next = p->a[l1];
			if (next == first) {
				break;
				}
			l1 = next;
			len ++;
			}
		if (len == 1) {
			l++;
			continue;
			}
		f_nothing_printed_at_all = FALSE;
		/* print cycle, starting at first: */
		l1 = first;
		strcat(str1, "(");
		while (TRUE) {
			sprintf(str2, "%ld", l1);
			strcat(str1, str2);
			next = p->a[l1];
			if (next == first) {
				break;
				}
			strcat(str1, " ");
			l1 = next;
			}
		strcat(str1, ")");
		}
}

INT sp_latex(SPERM *p, FILE *fp)
/* haengt an s an. */
{
	INT *have_seen = NIL;
	INT l, l1, first, next, len;
	INT f_nothing_printed_at_all = TRUE;
	BYTE str1[256];
	BYTE str2[256];
	
	str1[0] = 0;
	have_seen = (INT *) my_malloc(p->l * sizeof(INT), "sp_latex");
	if (have_seen == NIL) {
		Srfs("sp_latex", "no memory");
		return ERROR;
		}
	for (l = 0; l < p->l; l++) {
		have_seen[l] = FALSE;
		}
	l = 0;
	while (TRUE) {
		if (l >= p->l) {
			if (f_nothing_printed_at_all) {
				fprintf(fp, "id");
				}
			else {
				fprintf(fp, "%s", str1);
				}
			if (have_seen)
				my_free(have_seen);
			return OK;
			}
		if (have_seen[l]) {
			l++;
			continue;
			}
		/* cycle starting at l: */
		first = l;
		l1 = l;
		len = 1;
		while (TRUE) {
			have_seen[l1] = TRUE;
			next = p->a[l1];
			if (next == first) {
				break;
				}
			l1 = next;
			len ++;
			}
		if (len == 1) {
			l++;
			continue;
			}
		f_nothing_printed_at_all = FALSE;
		/* print cycle, starting at first: */
		l1 = first;
		strcat(str1, "(");
		while (TRUE) {
			sprintf(str2, "%ld", l1);
			strcat(str1, str2);
			next = p->a[l1];
			if (next == first) {
				break;
				}
			strcat(str1, " \\, ");
			l1 = next;
			}
		strcat(str1, ")");
		}
}

INT sp_test(void)
{
	SPERM p, q, r;
	BYTE s[256];
	
	s[0] = 0;
	sp_nil(&p);
	sp_nil(&q);
	sp_nil(&r);
	sp_int(&p, 5, "sp_test");
	sp_id(&p);
	sp_mult_apply_backwc_r(&p, 1, 3);
	sp_mult_apply_backwc_l(&p, 2, 3);
	/* sp_mult_apply_tau_r(&p, 0, 1);
	sp_mult_apply_tau_r(&p, 1, 2); */
	sp_inv(&p, &q);
	sp_mult(&p, &q, &r);
	strcat(s, "p = ");
	sp_sprint(&p, s);
	strcat(s, "; q = ");
	sp_sprint(&q, s);
	strcat(s, "; r = ");
	sp_sprint(&r, s);
	printf("\n%s\n", s);
	sp_free(&p);
	sp_free(&q);
	sp_free(&r);
	return OK;
}


void cp_nil(CPERM *p)
{
	p->l = 0;
	p->a = NIL;
}

INT cp_int(CPERM *p, INT l)
{
	cp_free(p);
	p->a = (UBYTE *) my_malloc(l, "cp_int");
	if (p->a == NIL) {
		Srfs("cp_int", "no memory");
		return ERROR;
		}
	p->l = l;
	return OK;
}

INT cp_free(CPERM *p)
{
	if (p->a) {
		my_free(p->a);
		p->a = NIL;
		}
	p->l = 0;
	return OK;
}

void cp_free_it(CPERM *p)
{
	cp_free(p);
	my_free(p);
}

INT cp_mv(CPERM *p, CPERM *q)
{
	INT i;
	
	if (p->l != q->l) {
		Srfs("cp_mv", "p->l != q->l");
		return ERROR;
		}
	for (i = 0; i < p->l; i++)
		q->a[i] = p->a[i];
	return OK;
}

INT cp_id(CPERM *p)
{
	INT i;
	
	for (i = 0; i < p->l; i++)
		p->a[i] = (char)i;
	return OK;
}

INT cp_mult(CPERM *a, CPERM *b, CPERM *c)
/* erst a, dann b; Ergebnis nach c */
{
	INT i, j;
	
	if (a->l != b->l) {
		Srfs("cp_mult", "a->l != b->l");
		return ERROR;
		}
	cp_int(c, a->l);
	for (i = 0; i < a->l; i++) {
		j = a->a[i];
		c->a[i] = b->a[j];
		}
	return OK;
}

INT cp_inv(CPERM *a, CPERM *b)
/* b:= a^-1 */
{
	INT i, j;
	
	cp_int(b, a->l);
	for (i = 0; i < a->l; i++) {
		j = a->a[i];
		b->a[j] = (UBYTE) i;
		}
	return OK;
}

INT cp_power(CPERM *a, CPERM *res, INT exp)
{
	CPERM *b = NIL;
	CPERM *c = NIL;
	CPERM *d = NIL;
	INT erg = ERROR;
	INT len;
	
	len = a->l;
	cp_int(res, len);
	b = (CPERM *) my_malloc(sizeof(CPERM), "cp_power");
	c = (CPERM *) my_malloc(sizeof(CPERM), "cp_power");
	d = (CPERM *) my_malloc(sizeof(CPERM), "cp_power");
	cp_nil(b);
	cp_nil(c);
	cp_nil(d);
	cp_int(b, len);
	cp_int(c, len);
	cp_int(d, len);
	
	cp_id(c);
	cp_mv(a, b);
	while (exp > 0) {
		/* res = b^exp * c */
		if (ODD(exp)) {
			erg += cp_mult(b, c, d);
			erg += cp_mv(d, c);
			exp--;
			continue; /* exp == 0 possible */
			}
		if (EVEN(exp)) {
			erg += cp_mult(b, b, d);
			erg += cp_mv(d, b);
			exp >>= 1;
			}
		}
	erg += cp_mv(c, res);

	erg = OK;
	if (b) {
		cp_free_it(b);
		b = NIL;
		}
	if (c) {
		cp_free_it(c);
		c = NIL;
		}
	if (d) {
		cp_free_it(d);
		d = NIL;
		}
	return erg;
}

INT cp_mult_apply_tau_r(CPERM *a, INT i, INT j)
/* a := a (i j). */
{
	INT k, i1, j1;
	
	for (k = 0; k < a->l; k++) {
		if (a->a[k] == i)
			i1 = k;
		if (a->a[k] == j)
			j1 = k;
		}
	/* now: i1 a == i, j1 a == j */
	a->a[i1] = (UBYTE) j;
	a->a[j1] = (UBYTE) i;
	return OK;
}

INT cp_mult_apply_tau_l(CPERM *a, INT i, INT j)
/* a := (i j) a. */
{
	INT i1, j1;
	
	i1 = a->a[i];
	j1 = a->a[j];
	/* now: i -> j -> j1; j -> i -> i1 */
	a->a[i] = (UBYTE) j1;
	a->a[j] = (UBYTE) i1;
	return OK;
}

INT cp_mult_apply_forwc_r(CPERM *a, INT i, INT l)
/* a := a (i i+1 ... i+l-1). */
{
	INT t[256], m, j, k;
	
	if (l > 256) {
		Srfs("cp_mult_apply_forwc_r", "l > 256");
		return ERROR;
		}
	for (m = 0; m < a->l; m++) {
		j = a->a[m];
		if (j >= i && (k = j - i) < l)
			t[k] = m;
		}
	/* now: t[k] -> i+k -> i+k+1 for 0 < k < l-1
	 *      t[l-1] -> i+l-1 -> i */
	for (k = 0; k < l - 1; k++)
		a->a[t[k]] = (UBYTE) (i + k + 1);
	a->a[t[l - 1]] = (UBYTE) i;
	return OK;
}

INT cp_mult_apply_backwc_r(CPERM *a, INT i, INT l)
/* a := a (i+l-1 i+l-2 ... i+1 i). */
{
	INT t[256], m, j, k;
	
	if (l > 256) {
		Srfs("cp_mult_apply_backwc_r", "l > 256");
		return ERROR;
		}
	for (m = 0; m < a->l; m++) {
		j = a->a[m];
		if (j >= i && (k = j - i) < l)
			t[k] = m;
		}
	/* now: t[k] -> i+k -> i+k-1 for 1 < k < l
	 *      t[0] -> i -> i+l-1 */
	for (k = 1; k < l; k++)
		a->a[t[k]] = (UBYTE) (i + k - 1);
	a->a[t[0]] = (UBYTE) (i + l - 1);
	return OK;
}

INT cp_mult_apply_forwc_l(CPERM *a, INT i, INT l)
/* a := (i i+1 ... i+l-1) a. */
{
	INT t, m;
	
	/* i+m -> i+m+1 -> a[i+m+1]  for m < l - 1 
	 * i+l-1 -> i -> a[i] */
	t = a->a[i];
	for (m = 0; m < l - 1; m++)
		a->a[i + m] = a->a[i + m + 1];
	a->a[i+l-1] = (UBYTE) t;
	return OK;
}

INT cp_mult_apply_backwc_l(CPERM *a, INT i, INT l)
/* a := (i+l-1 i+l-2 ... i+1 i) a. */
{
	INT t, m;
	
	/* i+m -> i+m-1 -> a[i+m-1]  for 1 < m < l 
	 * i -> i+l-1 -> a[i+l-1] */
	t = a->a[i + l - 1];
	for (m = l - 1; m > 0; m--)
		a->a[i + m] = a->a[i + m - 1];
	a->a[i] = (UBYTE) t;
	return OK;
}

INT cp_onep(CPERM *p)
{
	INT i;
	
	for (i = 0; i < p->l; i++)
		if (p->a[i] != i)
			return FALSE;
	return TRUE;
}

INT cp_cmp(CPERM *a, CPERM *b)
{
	INT i;
	
	if (a->l != b->l) {
		Srfs("cp_cmp", "a->l != b->l");
		return 0;
		}
	for (i = 0; i < a->l; i++) {
		if (a->a[i] < b->a[i])
			return -1;
		if (a->a[i] > b->a[i])
			return 1;
		}
	return 0;
}

void cp_print(CPERM *p)
{
	BYTE s[256];
	
	s[0] = 0;
	cp_sprint(p, s);
	printf("%s\n", s);
}

INT cp_sprint(CPERM *p, BYTE *s)
/* haengt an s an. */
{
	INT *have_seen = NIL;
	INT l, l1, first, next, len;
	INT f_nothing_printed_at_all = TRUE;
	BYTE str1[256];
	BYTE str2[256];
	
	if (p == NIL || s == NIL) {
		Srfs("cp_sprint", "args NIL");
		return ERROR;
		}
	str1[0] = 0;
	have_seen = (INT *) my_malloc(p->l * sizeof(INT), "cp_sprint");
	if (have_seen == NIL) {
		Srfs("cp_sprint", "no memory");
		return ERROR;
		}
	for (l = 0; l < p->l; l++) {
		have_seen[l] = FALSE;
		}
	l = 0;
	while (TRUE) {
		if (l >= p->l) {
			if (f_nothing_printed_at_all) {
				strcat(s, "id");
				}
			else {
				strcat(s, str1);
				}
			if (have_seen)
				my_free(have_seen);
			return OK;
			}
		if (have_seen[l]) {
			l++;
			continue;
			}
		/* cycle starting at l: */
		first = l;
		l1 = l;
		len = 1;
		while (TRUE) {
			have_seen[l1] = TRUE;
			next = p->a[l1];
			if (next == first) {
				break;
				}
			l1 = next;
			len ++;
			}
		if (len == 1) {
			l++;
			continue;
			}
		f_nothing_printed_at_all = FALSE;
		/* print cycle, starting at first: */
		l1 = first;
		strcat(str1, "(");
		while (TRUE) {
			sprintf(str2, "%ld", l1);
			strcat(str1, str2);
			next = p->a[l1];
			if (next == first) {
				break;
				}
			strcat(str1, " ");
			l1 = next;
			}
		strcat(str1, ")");
		}
}

INT cp_latex(CPERM *p, FILE *fp)
/* haengt an s an. */
{
	INT *have_seen = NIL;
	INT l, l1, first, next, len;
	INT f_nothing_printed_at_all = TRUE;
	BYTE str1[256];
	BYTE str2[256];
	
	str1[0] = 0;
	have_seen = (INT *) my_malloc(p->l * sizeof(INT), "cp_latex");
	if (have_seen == NIL) {
		Srfs("cp_latex", "no memory");
		return ERROR;
		}
	for (l = 0; l < p->l; l++) {
		have_seen[l] = FALSE;
		}
	l = 0;
	while (TRUE) {
		if (l >= p->l) {
			if (f_nothing_printed_at_all) {
				fprintf(fp, "id");
				}
			else {
				fprintf(fp, "%s", str1);
				}
			if (have_seen)
				my_free(have_seen);
			return OK;
			}
		if (have_seen[l]) {
			l++;
			continue;
			}
		/* cycle starting at l: */
		first = l;
		l1 = l;
		len = 1;
		while (TRUE) {
			have_seen[l1] = TRUE;
			next = p->a[l1];
			if (next == first) {
				break;
				}
			l1 = next;
			len ++;
			}
		if (len == 1) {
			l++;
			continue;
			}
		f_nothing_printed_at_all = FALSE;
		/* print cycle, starting at first: */
		l1 = first;
		strcat(str1, "(");
		while (TRUE) {
			sprintf(str2, "%ld", l1);
			strcat(str1, str2);
			next = p->a[l1];
			if (next == first) {
				break;
				}
			strcat(str1, " \\, ");
			l1 = next;
			}
		strcat(str1, ")");
		}
}

INT cp_test(void)
{
	CPERM p, q, r;
	BYTE s[256];
	
	s[0] = 0;
	cp_nil(&p);
	cp_nil(&q);
	cp_nil(&r);
	cp_int(&p, 5);
	cp_id(&p);
	cp_mult_apply_backwc_r(&p, 1, 3);
	cp_mult_apply_backwc_l(&p, 2, 3);
	/* cp_mult_apply_tau_r(&p, 0, 1);
	cp_mult_apply_tau_r(&p, 1, 2); */
	cp_inv(&p, &q);
	cp_mult(&p, &q, &r);
	strcat(s, "p = ");
	cp_sprint(&p, s);
	strcat(s, "; q = ");
	cp_sprint(&q, s);
	strcat(s, "; r = ");
	cp_sprint(&r, s);
	printf("\n%s\n", s);
	cp_free(&p);
	cp_free(&q);
	cp_free(&r);
	return OK;
}

/* 
 * IPERMs: permutations using chunk memory
 * (and INT handles therefore)
 */

INT iperm_test(void)
{
	INT deg = 5;
	CPERM p, q, r;
	INT ip, iq, ir;
	CHUNK_MEMH cm(
		deg /* entry_size */, 
		64 /* chunk_size */, 
		TRUE /* f_verbose */);
	
	cp_nil(&p);
	cp_nil(&q);
	cp_nil(&r);
	cp_int(&p, deg);
	cp_id(&p);
	cp_mult_apply_backwc_r(&p, 1, 3);
	cp_mult_apply_backwc_l(&p, 2, 3);
	ip = iperm_alloc(&cm);
	iq = iperm_alloc(&cm);
	ir = iperm_alloc(&cm);
	cm.print_info();

	cperm2iperm(&cm, &p, ip);
	printf("p = ");
	iperm_print(&cm, ip);
	
	iperm_invers(&cm, ip, iq);
	iperm_mult(&cm, ip, iq, ir);
	printf(" q = ");
	iperm_print(&cm, iq);
	printf(" r = ");
	iperm_print(&cm, ir);
	printf("\n");
	
	iperm_free(&cm, ip);
	iperm_free(&cm, iq);
	iperm_free(&cm, ir);
	cm.print_info();
	
	cp_free(&p);
	cp_free(&q);
	cp_free(&r);
	return OK;
}

INT iperm_sprint(CHUNK_MEMH *cm, INT a, BYTE *s)
{
	CPERM cp;
	
	cp_nil(&cp);
	cp.l = cm->entry_size;
	cp.a = (UBYTE *) cm->hdl2ptr(a);
	cp_sprint(&cp, s);
	cp_nil(&cp);
	return OK;
}

INT iperm_latex(CHUNK_MEMH *cm, INT a, FILE *fp)
{
	CPERM cp;
	
	cp_nil(&cp);
	cp.l = cm->entry_size;
	cp.a = (UBYTE *) cm->hdl2ptr(a);
	cp_latex(&cp, fp);
	cp_nil(&cp);
	return OK;
}

INT iperm_print(CHUNK_MEMH *cm, INT a)
{
	BYTE s[256];
	
	s[0] = 0;
	iperm_sprint(cm, a, s);
	printf("%s", s);
	return OK;
}

INT iperm_println(CHUNK_MEMH *cm, INT a)
{
	iperm_print(cm, a);
	printf("\n");
	return OK;
}

INT iperm_print_vec(CHUNK_MEMH *cm, VECTOR_OP V)
{
	INT i, l, a;
	
	l = V->s_li();
	for (i = 0; i < l; i++) {
		a = V->s_ii(i);
		iperm_println(cm, a);
		}
	return OK;
}

INT iperm_alloc(CHUNK_MEMH *cm)
{
	INT a;
	
	a = cm->chunk_alloc_hdl();
	if (a == -1) {
		printf("iperm_alloc(): no memory\n");
		return -1;
		}
	return a;
}

INT iperm_free(CHUNK_MEMH *cm, INT a)
{
	cm->chunk_free_hdl(a);
	return OK;
}

INT iperm_free_vec(CHUNK_MEMH *cm, VECTOR_OP V)
{
	INT i, l, a;
	
	l = V->s_li();
	for (i = 0; i < i; i++) {
		a = V->s_ii(i);
		iperm_free(cm, a);
		}
	return OK;
}

INT iperm_id(CHUNK_MEMH *cm, INT a)
{
	UBYTE *pa;
	INT i, l;

	pa = (UBYTE *) cm->hdl2ptr(a);
	l = cm->entry_size;
	for (i = 0; i < l; i++) {
		pa[i] = (UBYTE) i;
		}
	return OK;
}

INT iperm_mv(CHUNK_MEMH *cm, INT a, INT b)
/* b := a */
{
	UBYTE *pa, *pb;
	INT i, l;

	pa = (UBYTE *) cm->hdl2ptr(a);
	pb = (UBYTE *) cm->hdl2ptr(b);
	l = cm->entry_size;
	for (i = 0; i < l; i++) {
		pb[i] = pa[i];
		}
	return OK;
}

INT iperm_mult(CHUNK_MEMH *cm, INT a, INT b, INT c)
/* erst a, dann b; Ergebnis nach c */
{
	UBYTE *pa, *pb, *pc;
	INT i, j, l;

	pa = (UBYTE *) cm->hdl2ptr(a);
	pb = (UBYTE *) cm->hdl2ptr(b);
	pc = (UBYTE *) cm->hdl2ptr(c);
	l = cm->entry_size;
	for (i = 0; i < l; i++) {
		j = pa[i];
		pc[i] = pb[j];
		}
	return OK;
}

INT iperm_invers(CHUNK_MEMH *cm, INT a, INT b)
/* b := a^{-1} */
{
	UBYTE *pa, *pb;
	INT i, j, l;

	pa = (UBYTE *) cm->hdl2ptr(a);
	pb = (UBYTE *) cm->hdl2ptr(b);
	l = cm->entry_size;
	for (i = 0; i < l; i++) {
		j = pa[i];
		pb[j] = (UBYTE) i;
		}
	return OK;
}

INT iperm_invers_apply(CHUNK_MEMH *cm, INT *p_a)
/* a := a^{-1} 
 * (the handle will change). */
{
	UBYTE *pa, *pb;
	INT a, b, i, j, l;
	
	a = *p_a;
	b = iperm_alloc(cm);
	pa = (UBYTE *) cm->hdl2ptr(a);
	pb = (UBYTE *) cm->hdl2ptr(b);
	l = cm->entry_size;
	for (i = 0; i < l; i++) {
		j = pa[i];
		pb[j] = (UBYTE) i;
		}
	iperm_free(cm, a);
	*p_a = b;
	return OK;
}

INT iperm_onep(CHUNK_MEMH *cm, INT a)
{
	UBYTE *pa;
	INT i, l;
	
	pa = (UBYTE *) cm->hdl2ptr(a);
	l = cm->entry_size;
	for (i = 0; i < l; i++)
		if (pa[i] != i)
			return FALSE;
	return TRUE;
}

INT iperm_cmp(CHUNK_MEMH *cm, INT a, INT b)
{
	UBYTE *pa, *pb;
	INT i, l;
	
	pa = (UBYTE *) cm->hdl2ptr(a);
	pb = (UBYTE *) cm->hdl2ptr(b);
	l = cm->entry_size;
	for (i = 0; i < l; i++) {
		if (pa[i] < pb[i])
			return -1;
		if (pa[i] > pb[i])
			return 1;
		}
	return 0;
}

INT iperm_image_of(CHUNK_MEMH *cm, INT a, INT i)
{
	UBYTE *pa;
	
	pa = (UBYTE *) cm->hdl2ptr(a);
	return (INT) pa[i];
}

INT cperm2iperm(CHUNK_MEMH *cm, CPERM *cp, INT a)
{
	INT i, l;
	UBYTE *pa;
	
	l = cp->l;
	if (l != cm->entry_size) {
		printf("cperm2iperm(): wrong degree\n");
		return ERROR;
		}
	pa = (UBYTE *) cm->hdl2ptr(a);
	for (i = 0; i < l; i++)
		pa[i] = cp->a[i];
	return OK;
}

INT perm2iperm(CHUNK_MEMH *cm, PERMUTATION_OP p, INT a)
{
	INT i, l;
	UBYTE *pa;
	
	l = p->s_li();
	if (l != cm->entry_size) {
		printf("perm2iperm(): wrong degree\n");
		return ERROR;
		}
	pa = (UBYTE *) cm->hdl2ptr(a);
	for (i = 0; i < l; i++)
		pa[i] = (UBYTE) p->s_ii(i) - 1;
	return OK;
}

INT iperm2perm(CHUNK_MEMH *cm, INT a, PERMUTATION_OP p)
{
	INT i, l;
	UBYTE *pa;
	
	l = cm->entry_size;
	p->m_il(l);
	pa = (UBYTE *) cm->hdl2ptr(a);
	for (i = 0; i < l; i++)
		p->m_ii(i, (INT) pa[i] + 1);
	return OK;
}

INT perm2iperm_vec(CHUNK_MEMH *cm, VECTOR_OP V, VECTOR_OP W)
{
	PERMUTATION_OP p;
	INT i, l, a;
	
	l = V->s_li();
	W->m_il(l);
	for (i = 0; i < l; i++) {
		p = (PERMUTATION_OP) V->s_i(i);
		a = iperm_alloc(cm);
		perm2iperm(cm, p, a);
		W->m_ii(i, a);
		}
	return OK;
}

INT iperm2perm_vec(CHUNK_MEMH *cm, VECTOR_OP V, VECTOR_OP W)
{
	PERMUTATION_OB p;
	INT i, l, a;
	
	l = V->s_li();
	W->m_il(l);
	for (i = 0; i < l; i++) {
		a = V->s_ii(i);
		iperm2perm(cm, a, &p);
		p.swap(W->s_i(i));
		}
	return OK;
}

/* 
 * CHUNK MEMORY
 * handle 0 not used.
 */

chunk_memh::chunk_memh()
{
	entries_used = -1;
	entry_size = 0;
	entry_size_INT = 0;
	free_hdl = -1;
	mem = NIL;
	times = 0;
	chunk_size = 0;
	f_verbose = TRUE;
	printf("chunk_memh(): "
	"entry_size = %ld entry_size_INT = %ld\n"
	"chunk_size = %ld f_verbose = %ld\n", 
	this->entry_size, this->entry_size_INT, 
	this->chunk_size, this->f_verbose);
}

chunk_memh::chunk_memh(INT entry_size, INT chunk_size, INT f_verbose)
{
	INT es;
	
	if (entry_size <= 0) {
		error("chunk_memh(): entry_size <= 0");
		return;
		}
	entries_used = -1;
	es = entry_size;
	while ((es % sizeof(INT)) != 0)
		es++;
	this->entry_size = entry_size;
	entry_size_INT = es / sizeof(INT);
	free_hdl = -1;
	mem = NIL;
	times = 0;
	this->chunk_size = chunk_size;
	this->f_verbose = f_verbose;
	printf("chunk_memh(1, 2, 3): "
	"entry_size = %ld entry_size_INT = %ld\n"
	"chunk_size = %ld f_verbose = %ld\n", 
	this->entry_size, this->entry_size_INT, 
	this->chunk_size, this->f_verbose);
}

chunk_memh::~chunk_memh()
{
	INT i;

	for (i = 0; i < times; i++)
		my_free(mem_tbl[i]);
	times = 0;
	free_hdl = -1;
	entries_used = -1;
	entry_size = 0;
	entry_size_INT = 0;
	chunk_size = 0;
	mem = NIL;
}

void CHUNK_MEMH::print_info(void)
{
	INT l, hdl;
	INT *p;
	
	printf("CHUNK_MEMH::info: "
		"entry_size = %ld entry_size_INT = %ld\n"
		"chunk_size = %ld ", 
		entry_size, entry_size_INT, 
		chunk_size);
	fflush(stdout);
	l = 0;
	hdl = free_hdl;
	while (hdl != -1) {
		l++;
		p = (INT *) hdl2ptr(hdl);
		hdl = *p;
		}
	printf("free handles: %ld used chunks: %ld\n", l, times);
}

void *CHUNK_MEMH::hdl2ptr(INT hdl)
{
	INT i, j, *p;
	
	if (!hdl_check(hdl)) {
		printf("CHUNK_MEMH::hdl2ptr() handle illegal.\n");
		return NIL;
		}
	i = hdl / chunk_size; /* chunk number */
	j = hdl - i * chunk_size; /* entry number */
	p = &mem_tbl[i][j * entry_size_INT];
	return (void *) p;
}

INT CHUNK_MEMH::hdl_check(INT hdl)
{
	if (!is_legal_hdl(hdl)) {
		printf("illegal handle: %ld times = %ld chunk_size = %ld\n", 
			hdl, times, chunk_size);
		return FALSE;
		}
	return TRUE;
}

INT CHUNK_MEMH::is_legal_hdl(INT hdl)
{
	if (hdl <= 0 || hdl >= (times - 1) * chunk_size + entries_used)
		return FALSE;
	return TRUE;
}

INT CHUNK_MEMH::chunk_alloc_hdl(void)
{
	INT *p, hdl;
	
	if (free_hdl != -1) {
		hdl = free_hdl;
		p = (INT *) hdl2ptr(hdl);
		free_hdl = *p;
			/* maybe -1 if last 
			 * free entry is used again. */
		return hdl;
		}
	if (entries_used >= chunk_size || entries_used == -1) {
		if (!alloc_new_chunk()) {
			printf("CHUNK_MEMH::chunk_alloc(): no memory for new chunk");
			return -1;
			}
		if (times == 1)
			entries_used++;
				/* do not use hdl 0 */
		}
	hdl = (times - 1) * 
		chunk_size + entries_used;
	entries_used++;
	return hdl;
}

void CHUNK_MEMH::chunk_free_hdl(INT hdl)
{
	INT *p;
	
	if (!hdl_check(hdl)) {
		printf("CHUNK_MEMH::chunk_free_hdl() handle illegal.\n");
		return;
		}
	p = (INT *) hdl2ptr(hdl);
	*p = free_hdl;
	free_hdl = hdl;
}

INT CHUNK_MEMH::alloc_new_chunk(void)
{
	INT size;
	
	if (entry_size == 0) {
		printf("CHUNK_MEMH::alloc_new_chunk(): entry_size not initialized");
		return FALSE;
		}
	if (entry_size_INT == 0) {
		printf("CHUNK_MEMH::alloc_new_chunk(): entry_size_INT not initialized");
		return FALSE;
		}
	if (chunk_size == 0) {
		printf("CHUNK_MEMH::alloc_new_chunk(): chunk_size not initialized");
		return FALSE;
		}
	size = chunk_size * 
		entry_size_INT * sizeof(INT);
	mem = (INT *) my_malloc(size, "alloc_new_chunk");
	if (mem == NIL) {
		printf("CHUNK_MEMH::"
		"alloc_new_chunk(): "
		"no memory");
		return FALSE;
		}
	mem_tbl[times++] = mem;
	if (f_verbose) {
		printf("\n*** alloc_new_chunk(): %ld structs allocated (%ld times); %ld bytes\n", 
			chunk_size, times, size);
		fflush(stdout);
		}
	entries_used = 0;
	return TRUE;
}

#endif /* CP_TRUE */


