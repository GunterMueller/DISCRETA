/* iso_sub.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef GEO_TRUE

#include <DISCRETA/geo.h>

#define DEBUG_ISO

#define MAX_ISO 250000

#define MAX_TYPE 200
	/* at least 2 * MAX_VB + 1 */

#define MAX_SET_SIZE MAX_VB

typedef struct ordered_set ORDERED_SET;
typedef struct sub_iso SUB_ISO;

struct ordered_set {
	INT a[MAX_SET_SIZE];
	INT size;
};

struct sub_iso {
	ISO_INFO *info;
	ORDERED_SET *E; /* MAX_VB */
	ISO_GEO_DATA A, B;
};

static void subiso_print_beginning(SUB_ISO *iso, INT k);
static INT subiso_do_level(SUB_ISO *iso, INT i);
static void subiso_nil(SUB_ISO *iso);
static void subiso_exit(SUB_ISO *iso);
static INT subiso_print_Ei(SUB_ISO *iso, INT i);
static INT subiso_int(
	SUB_ISO *iso, ISO_INFO *info);


void igd_nil(ISO_GEO_DATA *geo)
{
	geo->f_pc = FALSE;
	geo->theX = NIL;
	geo->pctheX = NIL;
	geo->dim_n = 0;
	geo->R = NIL;
	geo->f_R_allocated = FALSE;
	geo->f_use_ddp = FALSE;
	geo->f_use_ddb = FALSE;
	geo->ddp = NIL;
	geo->ddb = NIL;
	geo->nb_bcol = 0;
	geo->nb_pcol = 0;
	geo->bcol = NIL;
	geo->pcol = NIL;
	cp_nil(&geo->p);
	cp_nil(&geo->pv);
	cp_nil(&geo->q);
	cp_nil(&geo->qv);
	geo->f_tdo = FALSE;
	geo->f_colors = FALSE;
	geo->f_transpose_it = FALSE;
}

void igd_exit(ISO_GEO_DATA *geo)
{
	if (geo->f_R_allocated)
		my_free(geo->R);
	geo->R = NIL;
	if (geo->theX)
		my_free(geo->theX);
	geo->theX = NIL;
	cp_free(&geo->p);
	cp_free(&geo->pv);
	cp_free(&geo->q);
	cp_free(&geo->qv);
}

void igd_init_theX(
	ISO_GEO_DATA *geo, INT *theX, INT *R, INT dim_n, 
	INT v, INT b)
{
	geo->f_pc = FALSE;
	geo->theX = theX;
	geo->R = R;
	geo->dim_n = dim_n;
	geo->v = v;
	geo->b = b;
}

void igd_init_pctheX(
	ISO_GEO_DATA *geo, BYTE *theX, INT *R, INT dim_n, 
	INT v, INT b)
{
	geo->f_pc = TRUE;
	geo->pctheX = theX;
	geo->R = R;
	geo->dim_n = dim_n;
	geo->v = v;
	geo->b = b;
}

INT igd_init_dd(
	ISO_GEO_DATA *geo, 
	INT f_ddp, SHORT *ddp, 
	INT f_ddb, SHORT *ddb)
{
	geo->f_use_ddp = f_ddp;
	if (f_ddp) {
		if (ddp == NIL)
			return error("igd_init_dd()"
			"ddp == NIL");
		geo->ddp = ddp;
		}
	geo->f_use_ddb = f_ddb;
	if (f_ddb) {
		if (ddb == NIL)
			return error("igd_init_dd()" 
			"ddb == NIL");
		geo->ddb = ddb;
		}
	return OK;
}

void igd_init_color(ISO_GEO_DATA *p, 
	INT nb_pcol, INT nb_bcol, INT *pcol, INT *bcol)
{
	p->nb_pcol = nb_pcol;
	p->nb_bcol = nb_bcol;
	p->pcol = pcol;
	p->bcol = bcol;
	p->f_colors = TRUE;
}

void igd_init_tdo(
	ISO_GEO_DATA *p, TDOSS *tdoss)
{
	INT i, j, len;

	len = tdoss_m(tdoss);
	p->tdo_m = len;
	for (i = 0; i < len; i++)
		p->tdo_V[i] = tdoss_Vi(tdoss, i);
	len = tdoss_n(tdoss);
	p->tdo_n = len;
	for (j = 0; j < len; j++)
		p->tdo_B[j] = tdoss_Bj(tdoss, j);
	p->f_tdo = TRUE;
}

void igd_init_tdo_V_B(ISO_GEO_DATA *p, 
	INT V, INT B, INT *Vi, INT *Bj)
{
	INT i, j;

	p->tdo_m = V;
	for (i = 0; i < V; i++)
		p->tdo_V[i] = Vi[i];
	p->tdo_n = B;
	for (j = 0; j < B; j++)
		p->tdo_B[j] = Bj[j];
	p->f_tdo = TRUE;
}

INT igd_init_igd(
	ISO_GEO_DATA *geo, ISO_GEO_DATA *geo_info)
{
	INT i, j, k;
	INT size;
	INT *theX1, *R1, dim_n1;

	size = geo_info->v * geo_info->dim_n * sizeof(INT);
	geo->theX = (INT *) my_malloc(size, "iso_sub");
	if (geo->theX == NIL)
		return error("igd_init_igd(): no memory");
	geo->dim_n = geo_info->dim_n;
	for (i = 0; i < geo_info->v; i++)
		for (j = 0; j < geo_info->dim_n; j++) {
			if (j < geo_info->R[i]) {
				if (geo_info->f_pc) {
					k = (INT) geo_info->pctheX
						[i * geo_info->dim_n + j];
					}
				else {
					k = geo_info->theX
						[i * geo_info->dim_n + j];
					}
				geo->theX[i * geo->dim_n + j] = k;
				}
			else {
				geo->theX[i * geo->dim_n + j] = 0;
				}
			}
	if (geo->f_transpose_it) {
#ifdef DEBUG_ISO
		printf("igd_init_igd(): vor transpose:\n");
		fflush(stdout);
#endif
#ifdef DM_TRUE
		if (!inc_transpose(geo_info->R, geo->theX, 
			FALSE /* f_full */, 
			geo->dim_n, geo_info->v, geo_info->b, 
			&theX1, &dim_n1, &R1))
			return error("idg_init_idg() error in "
				"inc_transpose()");
#endif
#ifdef DEBUG_ISO
		printf("igd_init_igd(): nach transpose\n");
		fflush(stdout);
#endif
		my_free(geo->theX);
		geo->theX = theX1;
		geo->dim_n = dim_n1;
		geo->R = R1;
		geo->v = geo_info->b;
		geo->b = geo_info->v;
		geo->f_R_allocated = TRUE;
		geo->tdo_m = geo_info->tdo_n;
		geo->tdo_n = geo_info->tdo_m;
		for (i = 0; i < geo_info->tdo_n; i++)
			geo->tdo_V[i] = geo_info->tdo_B[i];
		for (i = 0; i < geo_info->tdo_m; i++)
			geo->tdo_B[i] = geo_info->tdo_V[i];
		geo->f_use_ddp = geo_info->f_use_ddb;
		geo->ddp = geo_info->ddb;
		geo->f_use_ddb = geo_info->f_use_ddp;
		geo->ddb = geo_info->ddp;
		} /* if (f_transpose_it) */
	else {
		geo->R = geo_info->R;
		geo->f_R_allocated = FALSE;
		geo->v = geo_info->v;
		geo->b = geo_info->b;
		geo->tdo_m = geo_info->tdo_m;
		geo->tdo_n = geo_info->tdo_n;
		for (i = 0; i < geo_info->tdo_m; i++)
			geo->tdo_V[i] = geo_info->tdo_V[i];
		for (i = 0; i < geo_info->tdo_n; i++)
			geo->tdo_B[i] = geo_info->tdo_B[i];
		geo->f_use_ddp = geo_info->f_use_ddp;
		geo->ddp = geo_info->ddp;
		geo->f_use_ddb = geo_info->f_use_ddb;
		geo->ddb = geo_info->ddb;
		}
	return OK;
}

INT igd_init_hbar(ISO_GEO_DATA *p)
{
	INT i;

	if (p->f_colors)
		return error("igd_init_hbar() "
			"colors not supported yet");
	if (p->f_tdo)
		return error("igd_init_hbar() "
			"tdo in igd not supported yet");
	/* ein einziger hbar Bereich: */
	for (i = 0; i < p->V; i++) {
		p->hbar[i] = p->B;
		p->hlen[i] = 0;
		p->grid_entry[i] = 0;
		}
	p->hbar[0] = -1;
	p->hlen[0] = p->V;
	p->hlen1[0] = 0;
	p->hlen01[0] = p->v;
	p->G_max = 1;
	return OK;
}

void igd_init_cperms(ISO_GEO_DATA *p)
{
	cp_int(&p->p, p->V);
	cp_int(&p->pv, p->V);
	cp_int(&p->q, p->B);
	cp_int(&p->qv, p->B);
	cp_id(&p->p);
	cp_id(&p->pv);
	cp_id(&p->q);
	cp_id(&p->qv);
}

INT igd_find_incidence(ISO_GEO_DATA *p, 
	INT i, INT j)
{
	INT r, R, a;

	R = p->R[i];
	for (r = 0; r < R; r++) {
		a = p->theX[i * p->dim_n + r];
		if (a == j)
			return TRUE;
		if (a > j)
			return FALSE;
		}
	return FALSE;
}

INT igd_print_theX(ISO_GEO_DATA *p)
{
	INT len, len01;
	INT ge, i, ii, i0, j, j0, f_inc;
	INT v = p->V;
	INT b = p->B;

	for (j = 0; j < b; j++) {
		j0 = p->qv.a[j];
		printf("%ld ", j0);
		}
	i = 0;
	for (ge = 0; ge < p->G_max; ge++) {
#ifdef DEBUG_ISO
		printf("igd_print_theX()\n");
		printf("ge = %ld i = %ld hbar[i] = %ld\n", 
			ge, i, p->hbar[i]);
		printf("hlen[i] = %ld hlen01[i] = %ld hlen1[i] = %ld\n", 
			p->hlen[i], p->hlen01[i], p->hlen1[i]);
		fflush(stdout);
#endif
		len = p->hlen[i];
		len01 = p->hlen01[i];
		for (ii = 0; ii < len01; ii++) {
			i0 = p->pv.a[i + ii];
			printf("%2ld ", i0);
			for (j = 0; j < b; j++) {
				j0 = p->qv.a[j];
				f_inc = igd_find_incidence(p, i0, j0);
				if (f_inc)
					printf("X");
				else
					printf(".");
				}
			printf("\n");
			}
		for (ii = len01; ii < len; ii++) {
			i0 = p->pv.a[i + ii];
			printf("%2ld ", i0);
			printf("\n");
			}
		printf("-\n");
		i += len;
		}
	if (i != p->V)
		return error("igd_print_theX(): i != p->V");
	printf("\n");
	return OK;
}

INT igd_add_point(ISO_GEO_DATA *p, 
	INT k, INT j0, INT f_verbose)
{
	INT j, first, i, ii, ge, len, len01, len1;
	INT i0, f_inc;

	if (f_verbose) {
		printf("igd_add_point(%ld, %ld)\n", k, j0);
		fflush(stdout);
		}

	j = p->q.a[j0];
	if (j != k) {
		if (j < k) {
			return error("igd_add_point(): j < k");
			}
		/* Spaltenpermutation notwendig: */
		cp_mult_apply_tau_r(&p->q, k, j);
		cp_mult_apply_tau_l(&p->qv, k, j);
		}
	first = 0;
	for (i = 0; i < p->G_max; i++) {
		len = p->hlen[first];
		len01 = p->hlen01[first];
		len1 = 0;
		/* untersuche den Bereich der Zeilen, 
		 * welche echt vorhanden sind, 
		 * auf Einsen (keine Einschubzeilen): */
		for (ii = 0; ii < len01; ii++) {
			i0 = p->pv.a[first + ii];
			f_inc = igd_find_incidence(p, i0, j0);
			if (!f_inc)
				continue;
			if (len1 != ii) {
				/* Zeilenpermutation notwendig: */
				cp_mult_apply_tau_r(&p->p, 
					first + len1, first + ii);
				cp_mult_apply_tau_l(&p->pv, 
					first + len1, first + ii);
				}
			len1++;
			}
		p->hlen1[first] = len1;
			/* = Anzahl der Einsen im i-ten Bereich */

		first += len;
		}
	if (first != p->V)
		return error("igd_add_point(): first != p->V");
	return OK;
}

INT igd_join(ISO_GEO_DATA *A, ISO_GEO_DATA *B, 
	INT k, INT f_verbose)
{
	INT first, i, ii, ge;
	INT Alen, Alen01, Alen1;
	INT Blen, Blen01, Blen1;
	INT i0, f_inc, l;
	INT low_new_len, high_new_len;
	INT low_new_len01, high_new_len01;
	INT first1;

	if (f_verbose) {
		printf("igd_join()\n");
		fflush(stdout);
		}

	/* teste zuerst, ob join moeglich ist: */
	first = 0;
	for (i = 0; i < A->G_max; i++) {
		if (A->hbar[first] >= k) {
			igd_print_theX(A);
			igd_print_theX(B);
			error("igd_join(): A->hbar[first] >= k");
			return FALSE;
			}
		if (B->hbar[first] >= k) {
			error("igd_join(): B->hbar[first] >= k");
			return FALSE;
			}
		Alen = A->hlen[first];
		Blen = B->hlen[first];
		if (Alen != Blen) {
			error("igd_join(): Alen != Blen");
			return FALSE;
			}
		Alen01 = A->hlen01[first];
		Blen01 = B->hlen01[first];
		Alen1 = A->hlen1[first];
		Blen1 = B->hlen1[first];

		/* Sonderfaelle, wo nicht gesplittet 
		 * werden muss: */
		if (Blen1 == 0) {
			if (Alen1 > 0)
				return FALSE;
			goto loop1;
			}
		if (Blen1 == Blen) {
			if (Alen1 < Alen01)
				return FALSE;
			goto loop1;
			}

		if (Alen1 > Blen1)
			return FALSE;
				/* zu viele Kreuze in A */
		if (Alen - Alen01 < Blen1 - Alen1)
			return FALSE;
				/* zu wenig Fuellzeilen */
loop1:
		first += Alen;
		}

	first = 0;
	for (i = 0; i < A->G_max; i++) {
		Alen = A->hlen[first];
		Blen = B->hlen[first];
		Alen01 = A->hlen01[first];
		Blen01 = B->hlen01[first];
		Alen1 = A->hlen1[first];
		Blen1 = B->hlen1[first];
		/* Sonderfaelle, wo nicht gesplittet 
		 * werden muss: */
		if (Blen1 == 0)
			goto loop2;
		if (Blen1 == Blen)
			goto loop2;

		if ((Alen1 != Alen01) && (Alen1 < Blen1)) {
			/* es muss vorgeschoben werden: */
			l = Blen1 - Alen1;
			for (ii = 0; ii < l; ii++) {
				cp_mult_apply_tau_r(&A->p, 
					first + Alen1 + ii, 
					first + Alen01 + ii);
				cp_mult_apply_tau_l(&A->pv, 
					first + Alen1 + ii, 
					first + Alen01 + ii);
				}
			}
		low_new_len = Blen1;
		high_new_len = Alen - low_new_len;
		low_new_len01 = Alen1;
		high_new_len01 = Alen01 - Alen1;
		first1 = first + low_new_len;
		A->hlen[first] = low_new_len;
		B->hlen[first] = low_new_len;
		A->hlen01[first] = low_new_len01;
		B->hlen01[first] = low_new_len;
		A->hbar[first1] = k;
		B->hbar[first1] = k;
		A->hlen[first1] = high_new_len;
		B->hlen[first1] = high_new_len;
		A->hlen01[first1] = high_new_len01;
		B->hlen01[first1] = high_new_len;

loop2:
		first += Alen;
		}
	/* grid_entry setzen: */
	igd_set_ge(A);
	igd_set_ge(B);
	return TRUE;
}

INT igd_del_point(ISO_GEO_DATA *p, 
	INT k, INT f_A, INT f_verbose)
/* hiernach ist len1 nicht mehr gueltig */
/* Loescht alle hbars > k */
{
	INT i, first, first1, k1, ii, l;
	INT low_len, low_len01;
	INT high_len, high_len01;

	if (f_verbose) {
		printf("igd_del_point(): ");
		if (f_A)
			printf("A ");
		else
			printf("B ");
		printf("%ld\n", k);
		fflush(stdout);
		}

	first = 0;
	while (first < p->V) {
		low_len = p->hlen[first];
		low_len01 = p->hlen01[first];
		first1 = first + low_len;
		if (first1 >= p->V)
			break;
		k1 = p->hbar[first1];
		if (k1 >= k) {
			/* hbar bei first1 loeschen: */
			high_len = p->hlen[first1];
			high_len01 = p->hlen01[first1];
			if (f_A && (low_len01 < low_len) && 
				(high_len01 > 0)) {
				/* es muss geschoben werden: */
				l = low_len - low_len01;
				for (ii = 0; ii < l; ii++) {
					cp_mult_apply_tau_r(&p->p, 
						first + low_len - ii - 1, 
						first1 + high_len01 - ii - 1);
					cp_mult_apply_tau_l(&p->pv, 
						first + low_len - ii - 1, 
						first1 + high_len01 - ii - 1);
					}
				}
			p->hbar[first1] = p->B;
			p->hlen[first1] = 0;
			p->hlen01[first1] = 0;
			p->hlen[first] = low_len + high_len;
			p->hlen01[first] = low_len01 + high_len01;
			/* Eintrag wurde geloescht;
			 * folgender hbar kann nicht ebenfalls 
			 * in derselben Stufe entstanden sein: */
			first = first1 + high_len;

			}
		else
			first += low_len;
		}
	igd_set_ge(p);
	return OK;
}

void igd_set_ge(ISO_GEO_DATA *p)
{
	INT i, first, len, ii;

	first = 0;
	i = 0;
	while (first < p->V) {
		len = p->hlen[first];
		for (ii = 0; ii < len; ii++)
			p->grid_entry[first + ii] = i;

		first += len;
		i++;
		}
	p->G_max = i;
}

static void subiso_print_beginning(SUB_ISO *iso, INT k)
{
	INT j, j0;

	for (j = 0; j < k; j++) {
		j0 = iso->B.qv.a[j];
		printf("%ld ", j0);
		}
	printf("|");
	for (j = k; j < iso->B.B; j++) {
		j0 = iso->B.qv.a[j];
		printf("%ld ", j0);
		}
	printf("\n");
}

static INT subiso_do_level(SUB_ISO *iso, INT i)
/* Rueckgabe FALSE, 
 * falls wegen f_break_after_first
 * abgebrochen wurde. */
{
	INT j, j0, l;
	BYTE s[256];
	
	if (iso->info->f_very_verbose) {
		printf("*** subiso: level %ld ***\n", i);
		}

	if (i == iso->A.b) {
		iso->info->nb_isomorphisms++;
		if (iso->info->nb_isomorphisms > MAX_ISO)
			return FALSE;
		if (iso->info->f_verbose) {
			printf("*** found isomorphism nr.%ld !\n", 
				iso->info->nb_isomorphisms);
			subiso_print_beginning(iso, i);
			/* s[0] = 0;
			cp_sprint(&iso->B.qv, s);
			printf("%s\n", s); */
			for (j = 0; j < iso->A.b; j++) {
				j0 = iso->B.qv.a[j];
				printf("%ld ", j0);
				}
			printf("\n");
			if (iso->info->f_very_verbose) {
				igd_print_theX(&iso->A);
				igd_print_theX(&iso->B);
				}
			}
		
		/* if (iso->info->f_verbose) {
			igd_print_theX(&iso->A);
			igd_print_theX(&iso->B);
			} */

		if (iso->info->f_break_after_fst)
			return(FALSE);
				/* terminate program */
		else
			return(TRUE);
				/* proceed further */
		}

	l = 0;
	for (j = i; j < iso->B.B; j++) {
		j0 = iso->B.qv.a[j];
		iso->E[i].a[l] = j0;
		l++;
		}
	iso->E[i].size = l;

	/* if (iso->info->f_verbose) {
		subiso_print_Ei(iso, i);
		} */
	
	if (iso->E[i].size == 0) {
		return(TRUE);
		}
	
	while (iso->E[i].size > 0) {
	
		if (iso->info->f_very_verbose) {
			printf("*** subiso: level %ld *** \n", i);
			subiso_print_beginning(iso, i);
			subiso_print_Ei(iso, i);
			}

		igd_add_point(&iso->A, i, i, 
			iso->info->f_very_verbose);
		igd_add_point(&iso->B, i, iso->E[i].a[0], 
			iso->info->f_very_verbose);

		if (!igd_join(&iso->A, &iso->B, i, 
			iso->info->f_very_verbose))
			goto loop;
		
		if (!subiso_do_level(iso, i + 1)) {
			return FALSE;
			}
		
loop:
		/* delete first of E[i]: */
		for (j = 1; j < iso->E[i].size; j++) {
			iso->E[i].a[j - 1] = 
			iso->E[i].a[j];
			}
		iso->E[i].size--;

		igd_del_point(&iso->A, i, TRUE /* f_A */, 
			iso->info->f_very_verbose);
		igd_del_point(&iso->B, i, FALSE /* f_A */, 
			iso->info->f_very_verbose);

		
		}

	return TRUE;
}

static void subiso_nil(SUB_ISO *iso)
{
	iso->info = NIL;
	iso->E = NIL;
	igd_nil(&iso->A);
	igd_nil(&iso->B);
}

static void subiso_exit(SUB_ISO *iso)
{
	my_free(iso->E);
	iso->E = NIL;
	igd_exit(&iso->A);
	igd_exit(&iso->B);
}

static INT subiso_print_Ei(SUB_ISO *iso, INT i)
{
	INT j;
	ORDERED_SET *p;
	
	p = &iso->E[i];
	printf("E[%ld] = { ", i);
	for (j = 0; j < p->size; j++) {
		printf("%2ld ", p->a[j]);
		}
	printf("}\n");
	return TRUE;
}

static INT subiso_int(
	SUB_ISO *iso, ISO_INFO *info)
{
	INT i, j, k, I, first;
	INT last_p1, len, size, dim_n;
	INT *theX1, *R1, dim_n1;
	
	iso->info = info;

	igd_init_igd(&iso->A, &iso->info->A);
	igd_init_igd(&iso->B, &iso->info->B);
	iso->A.V = iso->B.v;
	iso->A.B = iso->B.b;
	iso->B.V = iso->B.v;
	iso->B.B = iso->B.b;
	igd_init_cperms(&iso->A);
	igd_init_cperms(&iso->B);
	
	iso->E = (ORDERED_SET *) 
		my_malloc(sizeof(ORDERED_SET) * 
			iso->A.B, "iso_sub");
	if (iso->E == NIL) {
		Srfs("iso2_int", 
			"no memory for iso->E");
		return FALSE;
		}

#if 0
	/* ein einziger hbar Bereich: */
	for (i = 0; i < info->v; i++) {
		iso->hbar[i] = info->v;
		iso->hlen[i] = 0;
		iso->grid_entry[i] = 0;
		}
	iso->hbar[0] = -1;
	iso->hlen[0] = info->v;
	iso->G_max = 1;
#endif
	/* init grid: */

	igd_init_hbar(&iso->A);
	igd_init_hbar(&iso->B);

#if 0
	igd_add_point(&iso->A, 0, 0);
	igd_add_point(&iso->B, 0, 0);
	i = igd_join(&iso->A, &iso->B, 0);
	printf("join = %ld\n", i);
	igd_print_theX(&iso->A);
	igd_print_theX(&iso->B);

	igd_add_point(&iso->A, 1, 1);
	igd_add_point(&iso->B, 1, 1);
	i = igd_join(&iso->A, &iso->B, 1);
	printf("join = %ld\n", i);
	igd_print_theX(&iso->A);
	igd_print_theX(&iso->B);

	igd_del_point(&iso->A, 1, TRUE /* f_A */);
	igd_del_point(&iso->B, 1, FALSE /* f_A */);
	igd_print_theX(&iso->A);
	igd_print_theX(&iso->B);
	fflush(stdout);
#endif

#if 0
	/* use the tdo info in iso->A (not iso->B): */
	iso->G_max = 0;
	for (i = 0; i < iso->A.V; i++) {
		iso->hbar[i] = iso->A.B;
		iso->hlen[i] = 0;
		}
	iso->G_max = iso->A.tdo_m;
	first = 0;
	for (I = 0; I < iso->G_max; I++) {
		len = iso->A.tdo_V[I];
		iso->hbar[first] = -1;
		iso->hlen[first] = len;
		last_p1 = first + len;
		for (i = first; i < last_p1; i++) {
			iso->grid_entry[i] = I;
			}
		first = last_p1;
		}
#endif
	return TRUE;
}

INT subiso_test(ISO_INFO *info)
{
	SUB_ISO *iso = NIL;
	
	iso = (SUB_ISO *) my_malloc(sizeof(SUB_ISO), "iso_sub");
	if (iso == NIL)
		return error("subiso_test() no memory for iso");
	subiso_nil(iso);
	if (!subiso_int(iso, info))
		return error("subiso_test() error in subiso_int");
	
	subiso_do_level(iso, 0);
	
	subiso_exit(iso);
	my_free(iso);
	return OK;
}

#endif /* GEO_TRUE */


