/* iso.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef GEO_TRUE

#include <DISCRETA/geo.h>

#undef DEBUG_HBAR

#undef DEBUG_ISO

#define MAX_ISO 250000

#define MAX_TYPE 200
	/* at least 2 * MAX_VB + 1 */

#define MAX_SET_SIZE MAX_VB

typedef struct ordered_set ORDERED_SET;
typedef struct iso2 ISO2;
typedef struct iso_grid ISO_GRID;

struct ordered_set {
	INT a[MAX_SET_SIZE];
	INT size;
};

struct iso2 {
	ISO_INFO *info;
	INT Adim_n;
	INT Bdim_n;
	INT *AtheX;
	INT *BtheX;
	INT *A_R;
	INT *B_R;
	INT f_A_R_allocated;
	INT f_B_R_allocated;
	INT v_t;
	INT b_t;
	INT tdo_m_t;
	INT tdo_n_t;
	INT *tdo_V_t;
	INT *tdo_B_t;
	ORDERED_SET *E; /* MAX_VB */
	
	INT f_use_ddp;
	INT f_use_ddb;

	SHORT *Addp; /* (v \atop 2) entries */
	SHORT *Bddp;
	SHORT *Addb; /* (b \atop 2) entries */
	SHORT *Bddb;

	/* for the partitioning of the rows: */
	INT hbar[MAX_VB];
	/* column that is responsible for 
	 * a hbar before this row. 
	 * if equal to b, 
	 * then no hbar is present. 
	 * initially hbar[i] = -1 
	 * for all i in i_hbar[]. */
	INT hlen[MAX_VB];
	INT grid_entry[MAX_VB];
	/* the number of hbars 
	 * that lie over this row;
	 * or: number into the 
	 * hbar fields for this row. 
	 * this will be used as 
	 * an index into type[][]. */
	INT G_max;
	/* number of grid_entries. */
	CPERM Ap;
		/* row permutation for A; 
		 * degree iso->v_t */
	CPERM Apv; /* Ap^-1 */
	CPERM Bp;
		/* row permutation for B; 
		 * degree iso->v_t */
	CPERM Bpv; /* Bp^-1 */
		/* given AtheX and BtheX, 
		 * applying Ap /Bp to the rows, 
		 * q to the columns, You get the 
		 * matrix C of the isomorphism. 
		 * q is inside ISO_GRID and depends on 
		 * the actual level k. */
	
	/* for the partitioning of the columns: */
	/* type_len = G_max
	 *   + 1 if derivative (Ad, Bd) is present
	 *   + k if (Add, Bdd) is present:
	 *         	(k is the actual level)
	 *         	for all the columns as yet in C. */
	ISO_GRID *G1; /* for A */
	ISO_GRID *G2; /* for B */

	/* ISO_GEO_DATA A, B;
	INT hlen1[MAX_VB]; */
};

struct iso_grid {
	INT m; /* = iso->b_t */
	INT n; /* = type_len */
	CPERM q;
		/* column permutation; 
		 * degree m */
	CPERM qv; /* q^-1 */
	INT type[MAX_VB][MAX_TYPE];
	INT G_max;
	INT first[MAX_GRID];
	INT len[MAX_GRID];
	INT type_idx[MAX_GRID];
	INT grid_entry[MAX_GRID];
};

static void iso2_nil(ISO2 *iso);
static INT iso2_int(ISO2 *iso, ISO_INFO *info);
static INT iso2_ext(ISO2 *iso);
static INT iso2_do_level(ISO2 *iso, INT i);
/* Rueckgabe FALSE, 
 * falls wegen f_break_after_first
 * abgebrochen wurde. */
static INT iso2_B_add_point(
	ISO2 *iso, INT k, INT b0);
static INT iso2_A_add_point(ISO2 *iso, INT k);
static INT iso2_A_del_point(ISO2 *iso, INT k);
static INT check_hbar(ISO2 *iso, INT k);
static INT iso2_check_iso(
	ISO2 *iso, INT k, ORDERED_SET *E);
static INT iso2_print_Ei(ISO2 *iso, INT i);
static INT iso2_calc_grid(
	ISO2 *iso, ISO_GRID *G, INT k, 
	INT *theX, INT dim_n, INT *d, 
	SHORT *dd, CPERM *p);
/* bislang: max_r konstant, statt r Array
 * type_len nur ohne 
 * zweite Ableitungen (d und dd). 
 * iso->grid_entry[] 
 * muss bereits berechnet sein. */
static INT iso2_radix_sort(
	ISO2 *iso, ISO_GRID *G, INT radix, 
	INT first, INT last);
static INT iso2_insert_idx(
	ISO_GRID *G, INT first, 
	INT len, INT radix, 
	INT search_this, INT *idx);
static INT print_hbar(ISO2 *iso);
static INT iso2_print_grid(ISO_GRID *G);


/*
 *
 */

INT iso_info_init_A_INT(
	ISO_INFO *info, INT *theA, INT f_full)
{
	info->pcAtheX = NIL;
	info->AtheX = theA;
	info->Af_full = f_full;
	info->f_pcA = FALSE;
	return TRUE;
}

INT iso_info_init_A_BYTE(
	ISO_INFO *info, BYTE *theA, INT f_full)
{
	info->pcAtheX = theA;
	info->AtheX = NIL;
	info->Af_full = f_full;
	info->f_pcA = TRUE;
	return TRUE;
}

INT iso_info_init_B_INT(
	ISO_INFO *info, INT *theB, INT f_full)
{
	info->pcBtheX = NIL;
	info->BtheX = theB;
	info->Bf_full = f_full;
	info->f_pcB = FALSE;
	return TRUE;
}

INT iso_info_init_B_BYTE(
	ISO_INFO *info, BYTE *theB, INT f_full)
{
	info->pcBtheX = theB;
	info->BtheX = NIL;
	info->Bf_full = f_full;
	info->f_pcB = TRUE;
	return TRUE;
}

INT iso_info_init_ddp(
	ISO_INFO *info, INT f_ddp, 
	SHORT *Addp, SHORT *Bddp)
{
	info->f_use_ddp = FALSE;
	if (f_ddp) {
		info->f_use_ddp = TRUE;
		if (Addp == NIL) {
			Srfs("iso_info_init_ddp", 
			"Addp == NIL");
			return FALSE;
			}
		info->Addp = Addp;
		if (Bddp == NIL) {
			Srfs("iso_info_init_ddp", 
			"Bddp == NIL");
			return FALSE;
			}
		info->Bddp = Bddp;
		}
	return TRUE;
}

INT iso_info_init_ddb(
	ISO_INFO *info, INT f_ddb, 
	SHORT *Addb, SHORT *Bddb)
{
	info->f_use_ddb = FALSE;
	if (f_ddb) {
		info->f_use_ddb = TRUE;
		if (Addb == NIL) {
			Srfs("iso_info_init_ddb", 
			"Addb == NIL");
			return FALSE;
			}
		info->Addb = Addb;
		if (Bddb == NIL) {
			Srfs("iso_info_init_ddb", 
			"Bddb == NIL");
			return FALSE;
			}
		info->Bddb = Bddb;
		}
	return TRUE;
}

INT iso_info_init_tdo(
	ISO_INFO *info, TDOSS *tdoss)
{
	INT i, j, len;

	len = tdoss_m(tdoss);
	info->tdo_m = len;
	for (i = 0; i < len; i++)
		info->tdo_V[i] = tdoss_Vi(tdoss, i);
	len = tdoss_n(tdoss);
	info->tdo_n = len;
	for (j = 0; j < len; j++)
		info->tdo_B[j] = tdoss_Bj(tdoss, j);
	return TRUE;
}

INT iso_info_init_tdo_V_B(ISO_INFO *info, 
	INT V, INT B, INT *Vi, INT *Bj)
{
	INT i, j;

	info->tdo_m = V;
	for (i = 0; i < V; i++)
		info->tdo_V[i] = Vi[i];
	info->tdo_n = B;
	for (j = 0; j < B; j++)
		info->tdo_B[j] = Bj[j];
	return TRUE;
}

INT init_ISO2()
{
	printf("size ISO_INFO = %ld\n", 
		(INT)sizeof(ISO_INFO));
	printf("size ISO2 = %ld\n", 
		(INT)sizeof(ISO2));
	printf("size ISO_GRID = %ld\n", 
		(INT)sizeof(ISO_GRID));
	return TRUE;
}

static void iso2_nil(ISO2 *iso)
{
	iso->info = NIL;
	iso->AtheX = NIL;
	iso->BtheX = NIL;
	iso->A_R = NIL;
	iso->B_R = NIL;
	iso->f_A_R_allocated = FALSE;
	iso->f_B_R_allocated = FALSE;
	iso->E = NIL;
	cp_nil(&iso->Ap);
	cp_nil(&iso->Apv);
	cp_nil(&iso->Bp);
	cp_nil(&iso->Bpv);
	iso->G1 = NIL;
	iso->G2 = NIL;

	/* igd_nil(&iso->A);
	igd_nil(&iso->B); */
}

static INT iso2_int(
	ISO2 *iso, ISO_INFO *info)
{
	INT i, j, k, I, first;
	INT last_p1, len, size, dim_n;
	/* TDOSS *tdoss; */
	INT *theX1, *R1, dim_n1;
	
	iso->info = info;

	size = info->v * 
		info->max_r * sizeof(INT);
	iso->AtheX = (INT *) my_malloc(size, "iso.C");
	iso->BtheX = (INT *) my_malloc(size, "iso.C");
	if (iso->AtheX == NIL || 
		iso->BtheX == NIL) {
		Srfs("iso2_int", 
		"no memory for iso->AtheX, BtheX");
		return FALSE;
		}
	if (info->Af_full)
		dim_n = MAX_R;
	else
		dim_n = info->max_r;
	for (i = 0; i < info->v; i++) 
		for (j = 0; j < info->max_r; j++) {
			if (info->f_pcA) {
				k = (INT) 
				info->pcAtheX[i * dim_n + j];
				}
			else {
				k = info->AtheX[i * dim_n + j];
				}
			iso->AtheX[i * info->max_r + j] = k;
			}
	iso->Adim_n = info->max_r;
	if (info->f_transpose_it) {
#ifdef DEBUG_ISO
		printf("iso - vor transpose\n");
		fflush(stdout);
#endif
#ifdef DM_TRUE
		if (!inc_transpose(
			info->R, iso->AtheX, 
			FALSE /* f_full */ , 
			info->max_r, info->v, info->b, 
			&theX1, &dim_n1, &R1)) {
			printf("iso2_int() "
			"error in inc_transpose()\n");
			if (info->f_pcA) {
				inc_theX_print_INT_BYTE(
					info->R, 
					info->pcAtheX, 
					FALSE /* f_INT */, 
					FALSE /* f_full */, 
					info->max_r, info->v);
				}
			else {
				inc_theX_print_INT_BYTE(
					info->R, info->AtheX, 
					TRUE /* f_INT */, 
					FALSE /* f_full */, 
					info->max_r, info->v);
				}
			error("iso2_int() stop");
			return FALSE;
			}
#endif
#ifdef DEBUG_ISO
		printf("iso - nach transpose\n");
		fflush(stdout);
#endif
		my_free(iso->AtheX);
		iso->AtheX = theX1;
		iso->Adim_n = dim_n1;
		iso->A_R = R1;
		iso->v_t = info->b;
		iso->b_t = info->v;
#ifdef DEBUG_ISO
		printf("info->b = %ld "
			"info->v = %ld\n", 
			info->b, info->v);
			fflush(stdout);
#ifdef DM_TRUE
		inc_theX_print_INT_BYTE(iso->A_R, 
			iso->AtheX, 
			TRUE /* f_INT */, 
			FALSE /* f_full */, 
			iso->Adim_n, iso->v_t);
		fflush(stdout);
#endif
#endif
		iso->f_A_R_allocated = TRUE;
		iso->tdo_m_t = info->tdo_n;
		iso->tdo_n_t = info->tdo_m;
		iso->tdo_V_t = info->tdo_B;
		iso->tdo_B_t = info->tdo_V;
		iso->f_use_ddp = info->f_use_ddb;
		iso->Addp = info->Addb;
		iso->f_use_ddb = info->f_use_ddp;
		iso->Addb = info->Addp;
		}
	else {
		iso->A_R = info->R;
		iso->f_A_R_allocated = FALSE;
		iso->v_t = info->v;
		iso->b_t = info->b;
		iso->tdo_m_t = info->tdo_m;
		iso->tdo_n_t = info->tdo_n;
		iso->tdo_V_t = info->tdo_V;
		iso->tdo_B_t = info->tdo_B;
		iso->f_use_ddp = info->f_use_ddp;
		iso->Addp = info->Addp;
		iso->f_use_ddb = info->f_use_ddb;
		iso->Addb = info->Addb;
		}

	if (info->Bf_full)
		dim_n = MAX_R;
	else
		dim_n = info->max_r;
	for (i = 0; i < info->v; i++) 
		for (j = 0; j < info->max_r; j++) {
			if (info->f_pcB) {
				k = (INT)
				info->pcBtheX[i * dim_n + j];
				}
			else {
				k = info->BtheX[i * dim_n + j];
				}
			iso->BtheX[i * info->max_r + j] = k;
			}
	iso->Bdim_n = info->max_r;
	if (info->f_transpose_it) {
#ifdef DM_TRUE
		inc_transpose(
			info->R, iso->BtheX, 
			FALSE /* f_full */ , 
			info->max_r, info->v, info->b, 
			&theX1, &dim_n1, &R1);
		my_free(iso->BtheX);
		iso->BtheX = theX1;
		iso->Bdim_n = dim_n1;
		iso->B_R = R1;
		iso->f_B_R_allocated = TRUE;
		iso->Bddp = info->Bddb;
		iso->Bddb = info->Bddp;
#endif
		}
	else {
		iso->B_R = info->R;
		iso->f_B_R_allocated = FALSE;
		iso->Bddp = info->Bddp;
		iso->Bddb = info->Bddb;
		}
	if (info->f_transpose_it) {
		for (i = 0; i < iso->v_t; i++) {
			if (iso->A_R[i] != iso->B_R[i]) {
				Srfs("iso2_int", 
				"iso->A_R[i] != iso->B_R[i]");
				return FALSE;
				}
			}
		}
	cp_int(&iso->Ap, iso->v_t);
	cp_int(&iso->Apv, iso->v_t);
	cp_int(&iso->Bp, iso->v_t);
	cp_int(&iso->Bpv, iso->v_t);
	cp_id(&iso->Ap);
	cp_id(&iso->Apv);
	cp_id(&iso->Bp);
	cp_id(&iso->Bpv);
	
	iso->E = (ORDERED_SET *) 
		my_malloc(sizeof(ORDERED_SET) * 
			iso->b_t, "iso.C");
	if (iso->E == NIL) {
		Srfs("iso2_int", 
			"no memory for iso->E");
		return FALSE;
		}
	iso->G1 = (ISO_GRID *) 
		my_malloc(sizeof(ISO_GRID), "iso.C");
	iso->G2 = (ISO_GRID *) 
		my_malloc(sizeof(ISO_GRID), "iso.C");
	if (iso->G1 == NIL || iso->G2 == NIL) {
		Srfs("iso2_int", "no memory for iso->G");
		return FALSE;
		}
	cp_nil(&iso->G1->q);
	cp_nil(&iso->G1->qv);
	cp_nil(&iso->G2->q);
	cp_nil(&iso->G2->qv);
	cp_int(&iso->G1->q, iso->b_t);
	cp_int(&iso->G1->qv, iso->b_t);
	cp_int(&iso->G2->q, iso->b_t);
	cp_int(&iso->G2->qv, iso->b_t);
	cp_id(&iso->G1->q);
	cp_id(&iso->G1->qv);
	cp_id(&iso->G2->q);
	cp_id(&iso->G2->qv);
	
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

	iso->G_max = 0;
	for (i = 0; i < iso->v_t; i++) {
		iso->hbar[i] = iso->b_t;
		iso->hlen[i] = 0;
		}
	iso->G_max = iso->tdo_m_t;
	first = 0;
	for (I = 0; I < iso->G_max; I++) {
		len = iso->tdo_V_t[I];
		iso->hbar[first] = -1;
		iso->hlen[first] = len;
		last_p1 = first + len;
		for (i = first; i < last_p1; i++) {
			iso->grid_entry[i] = I;
			}
		first = last_p1;
		}
	return TRUE;
}

static INT iso2_ext(ISO2 *iso)
{
	if (iso->f_A_R_allocated) {
		my_free(iso->A_R);
		iso->A_R = NIL;
		iso->f_A_R_allocated = FALSE;
		}
	if (iso->f_B_R_allocated) {
		my_free(iso->B_R);
		iso->B_R = NIL;
		iso->f_B_R_allocated = FALSE;
		}
	if (iso->AtheX) {
		my_free(iso->AtheX);
		iso->AtheX = NIL;
		}
	if (iso->BtheX) {
		my_free(iso->BtheX);
		iso->BtheX = NIL;
		}
	if (iso->E) {
		my_free(iso->E);
		iso->E = NIL;
		}
	cp_free(&iso->Ap);
	cp_free(&iso->Apv);
	cp_free(&iso->Bp);
	cp_free(&iso->Bpv);
	if (iso->G1) {
		cp_free(&iso->G1->q);
		cp_free(&iso->G1->qv);
		my_free(iso->G1);
		iso->G1 = NIL;
		}
	if (iso->G2) {
		cp_free(&iso->G2->q);
		cp_free(&iso->G2->qv);
		my_free(iso->G2);
		iso->G2 = NIL;
		}
	return TRUE;
}

INT iso_test(ISO_INFO *info)
{
	ISO2 *iso = NIL;
	
	iso = (ISO2 *) my_malloc(sizeof(ISO2), "iso.C");
	if (iso == NIL) {
		Srfs("iso_test", "no memory for iso");
		return FALSE;
		}
	iso2_nil(iso);
	if (!iso2_int(iso, info)) {
		Srff("iso_test", "iso2_int");
		return FALSE;
		}
	
	iso2_do_level(iso, 0);
	
	iso2_ext(iso);
	my_free(iso);
	return TRUE;
}


static INT iso2_do_level(ISO2 *iso, INT i)
/* Rueckgabe FALSE, 
 * falls wegen f_break_after_first
 * abgebrochen wurde. */
{
	INT j;
	BYTE s[256];
	
	if (iso->info->f_verbose) {
		printf("*** ISO2: level %ld ***\n", i);
		}

	if (i == iso->b_t - 1) {
		iso->info->nb_isomorphisms++;
		if (iso->info->nb_isomorphisms > MAX_ISO)
			return FALSE;
		if (iso->info->f_verbose) {
			printf("*** found isomorphism nr.%ld !\n", 
				iso->info->nb_isomorphisms);
			s[0] = 0;
			cp_sprint(&iso->Bp, s);
			printf("%s\n", s);
			}
		
		/* if (iso->info->f_verbose) {
			print_iso(A, isoA);
			print_iso(B, isoB);
			} */

		if (iso->info->f_break_after_fst)
			return(FALSE);
				/* terminate program */
		else
			return(TRUE);
				/* proceed further */
		}
		
	iso2_calc_grid(iso, iso->G1, i, iso->AtheX, 
		iso->Adim_n, iso->info->Ad, 
		iso->Addb, &iso->Ap);
	
	if (iso->info->f_verbose) {
		printf("A:\n");
		print_theX_pq(iso->AtheX, iso->Adim_n, 
			iso->v_t, iso->b_t, iso->A_R, 
			&iso->Apv, &iso->G1->qv);
		iso2_print_grid(iso->G1);
		}

	iso2_calc_grid(iso, iso->G2, i, 
		iso->BtheX, 
		iso->Bdim_n, iso->info->Bd, 
		iso->Bddb, &iso->Bp);
	
	if (iso->info->f_verbose) {
		printf("B:\n");
		print_theX_pq(iso->BtheX, iso->Bdim_n, 
			iso->v_t, iso->b_t, iso->B_R, 
			&iso->Bpv, &iso->G2->qv);
		iso2_print_grid(iso->G2);
		}
	
	iso2_check_iso(iso, i, &iso->E[i]);
	
	if (iso->info->f_verbose) {
		iso2_print_Ei(iso, i);
		}
	
	if (iso->E[i].size == 0L) {
		return(TRUE);
		}
	
	if (!iso2_A_add_point(iso, i)) {
		Srff("iso2_do_level", "iso2_A_add_point");
		printf("A: (i = %ld)\n", i);
		print_theX_pq(iso->AtheX, iso->Adim_n, 
			iso->v_t, iso->b_t, iso->A_R, 
			&iso->Apv, &iso->G1->qv);
		return FALSE;
		}
	
	if (iso->info->f_verbose) {
		print_hbar(iso);
		}

	while (iso->E[i].size > 0L) {
	
		if (!iso2_B_add_point(
			iso, i, iso->E[i].a[0])) {
			Srff("iso2_do_level", 
			"iso2_B_add_point");
			return FALSE;
			}
		
		if (!iso2_do_level(iso, i + 1L)) {
			return(FALSE);
			}
		
		/* delete first of E[i]: */
		for (j = 1; j < iso->E[i].size; j++) {
			iso->E[i].a[j - 1] = 
			iso->E[i].a[j];
			}
		iso->E[i].size--;
		if (iso->info->f_verbose) {
			iso2_print_Ei(iso, i);
			}
		
		}

	if (!iso2_A_del_point(iso, i)) {
		Srff("iso2_do_level", 
		"iso2_A_del_point");
		return FALSE;
		}
	
	if (iso->info->f_verbose) {
		printf("nach A_del_points():\n");
		print_hbar(iso);
		}

	return(TRUE);
}

static INT iso2_B_add_point(
	ISO2 *iso, INT k, INT b0)
{
	INT i, first, len, k0;
	INT b, last_len, j, l, j1, o;
	
	b = iso->G2->q.a[b0];
	if (b != k) {
		cp_mult_apply_tau_r(
			&iso->G2->q, k, b);
		cp_mult_apply_tau_l(
			&iso->G2->qv, k, b);
		}
	k0 = iso->G2->qv.a[k];
	if (k0 != b0) {
		Srfs("iso2_B_add_point", "k0 != b0");
		return FALSE;
		}
	first = 0; /* Zeile */
	i = 0; /* horizontal grid_entry */
	last_len = -1;
	while (TRUE) {
		/* we must have a hbar at first: */
		if (iso->hbar[first] > k) {
			Srfs("iso2_B_add_point", 
			"iso->hbar[first] > k");
			return FALSE;
			}
		len = iso->hlen[first];
		if (iso->hbar[first] == k) {
			if (last_len == -1) {
				Srfs("iso2_B_add_point", 
				"last_len == -1");
				return FALSE;
				}
			l = 0;
			for (j = -last_len; j < len; j++) {
				j1 = iso->Bpv.a[first + j];
				for (o = 0; 
					o < iso->B_R[j1]; o++) {
					/* before: max_r */
					if (iso->BtheX[(j1 * iso->Bdim_n + o)] 
						== k0)
						break;
					}
				if (o < iso->B_R[j1]) {
					/* this means: (j1, k0) is in BtheX. 
					 * swap row first - last_len + l 
					 * with first + j. */
					if (- last_len + l != j) {
						cp_mult_apply_tau_r(&iso->Bp, 
						first - last_len + l, first + j);
						cp_mult_apply_tau_l(&iso->Bpv, 
						first - last_len + l, first + j);
						}
					l++;
					if (l > last_len) {
						Srfs("iso2_B_add_point", 
						"l > last_len");
						return FALSE;
						}
					}
				else {
					/* (j1, k0) not in BtheX. */
					}
				} /* next j */
			} /* if (iso->hbar[first] == k) */
		i++;
		first += len;
		last_len = len;
		if (i >= iso->G_max)
			break;
		}
	if (first != iso->v_t) {
		Srfs("iso2_B_add_point", 
		"first != iso->v_t");
		return FALSE;
		}
	return TRUE;
}

static INT iso2_A_add_point(
	ISO2 *iso, INT k)
{
	INT i, first, len, k1;
	INT k_type_idx, new_len;
	INT ge, j, l, j1, o;
	
	if (!check_hbar(iso, k)) {
		Srff("iso2_A_add_point", "check_hbar");
		return FALSE;
		}
	k1 = iso->G1->q.a[k];
	k_type_idx = iso->G1->type_idx[k1];
	if (k != k1) {
		cp_mult_apply_tau_r(
			&iso->G1->q, k, k1);
		cp_mult_apply_tau_l(
			&iso->G1->qv, k, k1);
		}
	first = 0; /* Zeile */
	i = 0; /* h grid_entry */
	ge = 0; /* new h grid_entry */
	while (TRUE) {
		/* we must have a hbar at first: */
		if (iso->hbar[first] >= k) {
			Srfs("iso2_A_add_point", 
			"iso->hbar[first] >= k");
			return FALSE;
			}
		len = iso->hlen[first];
		new_len = iso->G1->type[k_type_idx][i];
		/* printf("i = %ld ge = %ld "
			"len = %ld new_len = %ld\n", 
			i, ge, len, new_len); */
		if (new_len && new_len < len) {
			/* add a new hbar: */
			iso->hbar[first + new_len] = k;
			iso->hlen[first] = new_len;
			iso->hlen[first + new_len] = 
				len - new_len;
			for (j = 0; j < new_len; j++)
				iso->grid_entry[first + j] = ge;
			ge++;
			for (j = new_len; j < len; j++)
				iso->grid_entry[first + j] = ge;
			
			/* resort the rows: */
			l = 0;
				/* the number of rows 
				 * with the k-th block;
				 * 0 <= l <= new_len. */
			for (j = 0; j < len; j++) {
				j1 = iso->Apv.a[first + j];
				for (o = 0; 
					o < iso->A_R[j1]; o++) {
					/* before: max_r */
					if (iso->AtheX
						[(j1 * iso->Adim_n + o)] 
						== k)
						break;
					}
				if (o < iso->A_R[j1]) {
					/* this means: (j1, k) is in AtheX. 
					 * swap row first + l 
					 * with first + j. */
					if (j > l) {
						cp_mult_apply_tau_r(
						&iso->Ap, first + l, 
						first + j);
						cp_mult_apply_tau_l(
						&iso->Apv, first + l, 
						first + j);
						}
					l++;
					if (l > new_len) {
						Srfs("iso2_A_add_point", 
						"l > new_len");
						return FALSE;
						}
					}
				else {
					/* (j1, k) not in AtheX. */
					}
				} /* next j */
			} /* if (new_len && new_len < len) */
		else {
			for (j = 0; j < len; j++)
				iso->grid_entry[first + j] = ge;
			}
		ge++;
		i++;
		first += len;
		if (i >= iso->G_max)
			break;
		}
	iso->G_max = ge;
#ifdef DEBUG_HBAR
	if (!check_hbar(iso, k + 1)) {
		Srff("iso2_A_add_point", 
		"check_hbar(k + 1)");
		return FALSE;
		}
#endif
	return TRUE;
}

static INT iso2_A_del_point(
	ISO2 *iso, INT k)
{
	INT i, first, new_len, last_len, ge, j;
	
	first = 0; /* Zeile */
	i = 0; /* new horizontal grid_entry */
	ge = 0; /* old horizontal grid_entry */
	last_len = -1;
	while (TRUE) {
		/* we must have a hbar at first: */
		if (iso->hbar[first] > k) {
			Srfs("iso2_A_del_point", 
			"iso->hbar[first] > k");
			return FALSE;
			}
		new_len = iso->hlen[first];
		if (iso->hbar[first] == k) {
			if (last_len == -1) {
				Srfs("iso2_A_del_point", 
				"last_len == -1");
				return FALSE;
				}
			iso->hbar[first] = iso->b_t;
				/* no hbar */
			iso->hlen[first - last_len] += new_len;
			ge--;
			for (j = 0; j < new_len; j++)
				iso->grid_entry[first + j] = ge;
			ge++;
			last_len = last_len + new_len;
			}
		else {
			for (j = 0; j < new_len; j++)
				iso->grid_entry[first + j] = ge;
			ge++;
			last_len = new_len;
			}
		first += new_len;
		i++;
		if (i >= iso->G_max)
			break;
		}
	if (first != iso->v_t) {
		Srfs("iso2_A_del_point", 
		"first != iso->v_t");
		return FALSE;
		}
	iso->G_max = ge;
#ifdef DEBUG_HBAR
	if (!check_hbar(iso, k)) {
		Srff("iso2_A_del_point", 
		"check_hbar");
		return FALSE;
		}
#endif
	return TRUE;
}

static INT check_hbar(ISO2 *iso, INT k)
{
	INT i, j, len, first;
	
	first = 0;
	i = 0;
	while (TRUE) {
		if (iso->hbar[first] >= k) {
			Srfs("check_hbar", 
			"iso->hbar[first] >= k");
			print_hbar(iso);
			return FALSE;
			}
		len = iso->hlen[first];
		for (j = 0; j < len; j++) {
			if (iso->grid_entry
				[first + j] != i) {
				Srfs("check_hbar", 
				"iso->grid_entry"
				"[first + j] != i");
				print_hbar(iso);
				return FALSE;
				}
			}
		i++;
		first += len;
		if (i >= iso->G_max)
			break;
		}
	if (first != iso->v_t) {
		Srfs("check_hbar", 
			"first != iso->v_t");
		print_hbar(iso);
		return FALSE;
		}
	return TRUE;
}

static INT iso2_check_iso(
	ISO2 *iso, INT k, ORDERED_SET *E)
{
	INT i, j, first, i1, i2, k1, ge;
	
	E->size = 0L;
	if (iso->G1->G_max != iso->G2->G_max)
		return TRUE;
	for (i = 0; i < iso->G1->G_max; i++) {
		if (iso->G1->first[i] != 
			iso->G2->first[i])
			return TRUE;
		if (iso->G1->len[i] != 
			iso->G2->len[i])
			return TRUE;
		}
	for (i = 0; i < iso->G1->G_max; i++) {
		first = iso->G1->first[i];
		i1 = iso->G1->type_idx[first];
		i2 = iso->G2->type_idx[first];
		for (j = 0; j < iso->G1->n; j++) {
			if (iso->G1->type[i1][j] != 
				iso->G2->type[i2][j])
				return TRUE;
			}
		}
	/* now: a bijection exists. */
	
	if (k != iso->G2->first[0]) {
		Srfs("iso2_check_iso", 
		"k != iso->G2->first[0]");
		return FALSE;
		}
	k1 = iso->G1->q.a[k];
		/* der k-te Block liegt momentan 
		 * in der k1-ten Spalte. */
	ge = iso->G1->grid_entry[k1];
	/* links (A): nimm den Block, 
	 *    in dem k liegt. 
	 * rechts (B): alle Bloecke 
	 *   des Grideintrags ge
	 * nach E[k]. */
	first = iso->G2->first[ge];
	for (i = 0; i < iso->G2->len[ge]; i++) {
		i1 = first + i;
		i2 = iso->G2->qv.a[i1];
		/* die urspruengliche Blocknummer. */
		E->a[E->size++] = i2;
		}
	return TRUE;
}

static INT iso2_print_Ei(ISO2 *iso, INT i)
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

static INT iso2_calc_grid(
	ISO2 *iso, ISO_GRID *G, INT k, 
	INT *theX, INT dim_n, INT *d, 
	SHORT *dd, CPERM *p)
/* (bislang: max_r konstant, 
 * statt R Array) jetzt: R[]
 * type_len nur ohne 
 *   zweite Ableitungen (d und dd). 
 * iso->grid_entry[] muss 
 *   bereits berechnet sein. */
{
	INT len, i, j, l, x, i1, j1, x1, ge;
	INT begin_d, begin_dd, d1;
	/* TDOSS *tdoss; */
	INT first, tdo_n, Bj;
	
	G->m = iso->b_t;
	G->n = iso->G_max;
	begin_d = G->n;
	if (iso->info->f_use_d) {
		G->n++;
		begin_dd = begin_d + 1;
		}
	else
		begin_dd = begin_d + 1;
	if (iso->f_use_ddb) {
		G->n += k;
		}
	len = G->m - k;
	for (i = 0; i < len; i++)
		for (j = 0; j < G->n; j++)
			G->type[i][j] = 0;
	for (i = 0; i < G->m; i++) {
		if (i >= k)
			G->type_idx[i] = i - k;
		else
			G->type_idx[i] = - 1;
		}
	for (i = 0; i < iso->v_t; i++) {
		for (j = 0; j < iso->A_R[i]; j++) {
			/* before: max_r */
			x = theX[i * dim_n + j];
			i1 = p->a[i];
			x1 = G->q.a[x];
			if (x1 >= k) {
				ge = iso->grid_entry[i1];
				G->type[x1 - k][ge]++;
				}
			}
		}
	G->G_max = 0;
	G->first[0] = k;
	
	if (iso->info->f_use_d) {
		for (i = 0; i < G->m; i++) {
			d1 = d[i];
			i1 = G->q.a[i];
			if (i1 >= k) {
				G->type[i1 - k][begin_d] = d1;
				}
			}
		}
	if (iso->f_use_ddb) {
		for (i = 0; i < G->m; i++) {
			for (j = 0; j < G->m; j++) {
				if (j == i)
					continue;
				l = ij2k(i, j, G->m);
				d1 = dd[l];
				i1 = G->q.a[i];
				j1 = G->q.a[j];
				if (i1 >= k && j1 < k) {
					G->type[i1 - k]
					[begin_dd + j1] = d1;
					}
				}
			}
		}

	if (iso->info->f_very_verbose) {
		printf("in iso2_calc_grid() (k = %ld):\n", k);
		for (i = 0; i < len; i++) {
			for (j = 0; j < G->n; j++) {
				printf("%ld ", G->type[i][j]);
				}
			printf("\n");
			}
		}
	
	/* tdoss = iso->info->tdo; */
	/* tdo_n = tdoss_n(tdoss); */
	tdo_n = iso->tdo_n_t;
	first = 0;
	for (j = 0; j < tdo_n; j++) {
		/* Bj = tdoss_Bj(tdoss, j); */
		Bj = iso->tdo_B_t[j];
		if (k < first + Bj) {
			iso2_radix_sort(
			iso, G, 0, k, first + Bj - 1);
			first += Bj;
			break;
			}
		first += Bj;
		}
	for (j++; j < tdo_n; j++) {
		/* Bj = tdoss_Bj(tdoss, j); */
		Bj = iso->tdo_B_t[j];
		iso2_radix_sort(
		iso, G, 0, first, first + Bj - 1);
		first += Bj;
		}
	return TRUE;
}

static INT iso2_radix_sort(
	ISO2 *iso, ISO_GRID *G, INT radix, 
	INT first, INT last)
{
	INT f_found, idx, k, l, t, i1, j;
	INT first0, first1, res, k1;
	CPERM *perm, *perm_inv;

	if (first == last || radix == G->n) {
		/* Berech first .. last 
		 * als neuen grid entry nach G eintragen. 
		 * grid_entry[first..last] setzen. */
		k = G->G_max;
		if (G->first[k] != first) {
			Srfs("iso2_radix_sort", 
			"G->first[k] != first");
			return TRUE;
			}
		for (l = first; l <= last; l++)
			G->grid_entry[l] = k;
		G->len[k] = last - first + 1;
		G->first[k + 1] = last + 1;
		G->G_max++;
		if (iso->info->f_very_verbose) {
			printf("iso2_radix_sort()|"
			"new entry: first = %ld len = %ld : ", 
				(INT)first, (INT)G->len[k]);
			i1 = G->type_idx[first];
			for (j = 0; j < G->n; j++) {
				printf("%ld ", (INT)G->type[i1][j]);
				}
			printf("\n");
			}
		return TRUE;
		}
	for (k = first; k <= last; k++) {
		f_found = iso2_insert_idx(
			G, first, k - first, 
			radix, k, &idx);
		if (idx != k) {
			/* s = (idx idx+1 ... k) 
			 * auf den aktuellen Stand 
			 * der Matrix anwenden. 
			 *   q := q * s, qv := s^-1 * qv */
			perm = &G->q;
			perm_inv = &G->qv;
			cp_mult_apply_forwc_r(
			perm, idx, k - idx + 1);
			cp_mult_apply_backwc_l(
			perm_inv, idx, k - idx + 1);
			t = G->type_idx[k];
			for (l = k; l > idx; l--)
				G->type_idx[l] = G->type_idx[l - 1];
			G->type_idx[idx] = t;
			/* grid_entry ist noch 
			 * nicht gesetzt, 
			 * braucht nicht mitge-
			 * schoben zu werden. */
			}
		}
	first0 = first;
	first1 = G->type_idx[first0];
	for (k = first + 1; k <= last; k++) {
		k1 = G->type_idx[k];
		res = G->type[k1][radix] - 
			G->type[first1][radix];
		if (res > 0) {
			Srfs("iso2_radix_sort", 
			"not descending");
			return FALSE;
			}
		if (res < 0) {
			iso2_radix_sort(
			iso, G, radix + 1, first0, k - 1);
			first0 = k;
			first1 = G->type_idx[first0];
			}
		if (k == last) {
			iso2_radix_sort(
			iso, G, radix + 1, first0, k);
			}
		}
	return (TRUE);
}

static INT iso2_insert_idx(
	ISO_GRID *G, INT first, INT len, 
	INT radix, INT search_this, 
	INT *idx)
{
	INT i, st1, cur, cur1, res;
	INT f_found;
	
	st1 = G->type_idx[search_this];
	f_found = FALSE;
	for (i = 0; i < len; i++) {
		cur = first + i;
		cur1 = G->type_idx[cur];
		res = G->type[cur1][radix] - 
		G->type[st1][radix];
		if (res == 0)
			f_found = TRUE;
		if (res < 0) {
			*idx = cur;
			return f_found;
			}
		}
	*idx = first + len;
	return f_found;
}

INT print_theX_pq(
	INT *theX, INT dim_n, 
	INT v, INT b, INT *R, 
	CPERM *pv, CPERM *qv)
/* before: INT r */
{
	INT i, j, i1, j1, o;
	
	printf("   ");
	for (j = 0; j < b; j++) {
		printf("%2ld ", (INT)qv->a[j]);
		}
	printf("\n");
	for (i = 0; i < v; i++) {
		i1 = pv->a[i];
		printf("%2ld ", i1);
		for (j = 0; j < b; j++) {
			j1 = qv->a[j];
			for (o = 0; o < R[i1]; o++) {
				/* before: r */
				if (theX[i1 * dim_n + o] == j1)
					break;
				}
			if (o < R[i1])
				printf(" X ");
			else
				printf(" . ");
			}
		printf("\n");
		}
	printf("\n");
	return TRUE;
}

INT print_theX(INT *theX, 
	INT dim_n, INT v, INT b, INT *R)
{
	INT i, j, i1, j1, o;
	
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			for (o = 0; o < R[i]; o++) {
				/* before: r */
				if (theX[i * dim_n + o] == j)
					break;
				}
			if (o < R[i])
				printf("X");
			else
				printf(".");
			}
		printf("\n");
		}
	printf("\n");
	return TRUE;
}

static INT print_hbar(ISO2 *iso)
{
	INT i, j, len, first;
	
	for (i = 0; i < iso->info->v; i++)
		printf("%ld %ld\n", 
		iso->hbar[i], iso->hlen[i]);
	first = 0;
	i = 0;
	while (TRUE) {
		len = iso->hlen[first];
		printf("i = %ld: at %ld len %ld\n", 
		i, first, len);
		i++;
		first += len;
		if (i >= iso->G_max)
			break;
		}
	return TRUE;
}

static INT iso2_print_grid(ISO_GRID *G)
{
	INT i, j, first, i1;
	
	for (i = 0; i < G->G_max; i++) {
		printf("at %ld %ld x (", 
			G->first[i], G->len[i]);
		first = G->first[i];
		i1 = G->type_idx[first];
		for (j = 0; j < G->n; j++) {
			printf("%ld", G->type[i1][j]);
			if (j < G->n)
				printf(" ");
			}
		printf(")\n");
		}
	printf("total: %ld\n", 
		G->first[G->G_max]);
	return TRUE;
}

#endif /* GEO_TRUE */

