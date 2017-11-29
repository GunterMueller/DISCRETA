/* dc_draw.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef LADDER_TRUE

#include <DISCRETA/perm.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/lb.h>
#include <DISCRETA/cp.h>
#include <DISCRETA/divs.h>
#include <DISCRETA/ladder.h>
#include <DISCRETA/fga.h>

#include <DISCRETA/gfq.h> /* for PSU_3(q^2) and Sz(q) */

#define MAX_STEP 64

#if TEXDOCU
DC_LATTICE *dc_graphs(INT deg, BYTE *path)
#else
/* path without trailing '/' 
 * produces path/g-%ld.ps and the copy 
 * path/g.ps
 */
#endif
{
	VECTOR_OB G, A2;
	SYM_OB id;
	DCY_OP dc; /* [MAX_STEP] */
	INT erg = OK, i;
	INT deg2, up_to_k, up_to_step;
	INT type;
	void *data;
	INT t0, t1, user_time;
	BYTE str[256];
	DC_LATTICE *dcl;
	
	t0 = os_ticks();
	
	dc = (DCY_OP) my_malloc(MAX_STEP * sizeof(DCY_OB), "dc_graphs");
	if (dc == NIL) {
		printf("dc_graphs(): no memory for dc\n");
		return NIL;
		}
		
	deg2 = (deg * (deg - 1)) >> 1;
	up_to_k = deg2 - 1;
	up_to_step = 2 * up_to_k + 1;
	if (up_to_step >= MAX_STEP) {
		printf("dc_graphs(): "
		"up_to_step too large\n");
		return NIL;
		}
	for (i = 0; i <= up_to_step; i++) {
		dc[i].c_obj_k(EMPTY);
		}
	
	/* 
	 * the group G: 
	 */
	printf("initializing G:\n"); fflush(stdout);
	{
	VECTOR_OB G1;
	
	erg += symmetric_generators(&G1, deg2);
	type = DO_TYPE_SYM;
	data = NIL;
	G1.swap(&G);
	v_do_println(&G, 
		FALSE /* f_numerated */, 
		1 /* n_on_a_row */, 
		type, data);
	}

	/* 
	 * the group A: 
	 */
	printf("initializing A:\n"); fflush(stdout);
	{
	VECTOR_OB A, A1;
	
	erg += symmetric_generators(&A, deg);
	erg += A1.m_il(A.s_li());
	for (i = 0; i < A.s_li(); i++)
		((PERMUTATION_OP) A.s_i(i))->induce2(
			(PERMUTATION_OP) A1.s_i(i), deg);
	A1.swap(&A2);
	}
	A2.println();
		
	/*
	 * subgroup ladder: 
	 */
	printf("initializing subgroup ladder:\n");
	fflush(stdout);
	{
	VECTOR_OB T;
	
	for (i = 0; i <= up_to_step; i++) {
		/* printf("step %ld:\n", i);
		fflush(stdout); */
		dc[i].initialize_Young(&T, 
			deg2, i /* step */, 
			0 /* type */, NIL /* data */);
		/* printf("T = ");
		T.println(); */
		T.swap(dc[i].s_T());
		}
	}
		
	printf("initializing id, D, Ad:\n");
	fflush(stdout);
	do_copy(A2.s_i(0), &id, type, data);
	do_one(&id, type, data);
	dc[0].s_D()->m_il(1);
	do_copy(&id, dc[0].s_D_i(0), type, data);
	
	dc[0].s_Ad()->m_il(1);
	((SYM_OP) &A2)->swap(
		dc[0].s_Ad()->s_i(0));
	
	printf("D = ");
	dc[0].s_D()->println();
	printf("Ad = ");
	dc[0].s_Ad()->println();

	printf("computing subgroup ladder:\n");
	fflush(stdout);
	for (i = 1; i <= up_to_step; i++)
		erg += dc_do_step(&dc[0], i, &id, 
		deg2 /* deg */, 
		TRUE /* f_verbose */, 
		type, data);
	
	{
	BYTE s[256], path_name[256], file_name[256];
	
	sprintf(file_name, "g_%ld", deg);

#if 0
#ifdef WWW_REQUEST
	strcpy(path_name, "/source/www/anton/DC/");

	sprintf(s, "rm %s%s.ps", 
		path_name, file_name);
	call_system(s);

	sprintf(s, "rm %sg.ps", 
		path_name);
	call_system(s);

	sprintf(s, "rm %sg.xbm", 
		path_name);
	call_system(s);
#else
	strcpy(path_name, "");
#endif
#endif

	strcpy(path_name, path);
	strcat(path_name, "/");
	sprintf(s, "%s%s.ps", path_name, file_name);
	
	dcl = open_dcl(dc, deg /* n */, deg2 /* n2 */, 
		up_to_step, 
		FALSE /* f_with_perm */, 
		TRUE /* f_verbose */, 
		1 /* width_in_pages */, 
		s, type, data);

#if 0
#ifdef WWW_REQUEST
	sprintf(s, "chmod ugo+w %s%s.ps", 
		path_name, file_name);
	call_system(s);
	
	sprintf(s, "cp %s%s.ps %sg.ps", 
		path_name, file_name, 
		path_name);
	call_system(s);

	sprintf(s, "chmod ugo+w %sg.ps", 
		path_name);
	call_system(s);
	
	fflush(stdout);
	sprintf(s, "cd %s; /usr/local/bin/pstoxbm"
		" g.ps g.xbm", 
		path_name);
	call_system(s);
	
#endif
#endif


	sprintf(s, "cp %s%s.ps %sg.ps", 
		path_name, file_name, 
		path_name);
	call_system(s);

	}
	
	t1 = os_ticks();
	user_time = t1 - t0;
	strcpy(str, "Running time for dc_graphs(): ");
	str[0] = 0;
	print_delta_time(user_time, str);
	printf("%s\n", str);
	
	return dcl;
}

#if TEXDOCU
DC_LATTICE *open_dcl(DCY_OP dc0, INT n, INT n2, 
	INT up_to_step, INT f_with_perm, 
	INT f_verbose, INT width_in_pages, 
	BYTE *file_name, INT type, void *data)
#endif
{
	DC_LATTICE *p;
	MATRIX_OP Ainf_t = (MATRIX_OP) callocobject("open_dcl");
	VECTOR_OP nl = (VECTOR_OP) callocobject("open_dcl");
	VECTOR_OP orbit_size = (VECTOR_OP)callocobject("open_dcl");
	MEM_OP plaz = (MEM_OP) callocobject("open_dcl");
	MEM_OP o_dx = (MEM_OP) callocobject("open_dcl");
	MATRIX_OB Acover;
	SYM_OB go;
	INT i, j, k, nb_d, step, dc_idx, ol;
	DCY_OP L;
	// VECTOR_OB layer;
	
	p = (DC_LATTICE *) my_malloc(sizeof(DC_LATTICE), "open_dcl");
	if (p == NIL)
		return NIL;
	dc_calc_go(dc0, &go);
	p->type = type;
	p->data = data;
	p->f_with_perm = f_with_perm;
	p->dc = dc0;
	p->n = n;
	p->n2 = n2;
	p->up_to_step = up_to_step;
	p->nb_d = dc_calc_nb_d(dc0, up_to_step);
	k = (up_to_step + 1) / 2;
	
	dc_calc_Ainf_t(dc0, up_to_step, Ainf_t, /* &layer, */ f_verbose);
	// printf("layer vector: ");
	// layer.println();
	Ainf_t->Asup2Acover(&Acover);
	Acover.Acover2nl(nl);
	p->Ainf_t = Ainf_t;
	p->nl = nl;
	p->orbit_size = orbit_size;

	orbit_size->m_il(p->nb_d);
	for (i = 0; i < p->nb_d; i++) {
		dc_dc_no_to_dc_idx(dc0, i, &step, &dc_idx);
		ol = dc_orbit_size(dc0, &go, step, dc_idx, FALSE /* f_v */);
		printf("dc no %ld in step %ld dc_idx = %ld: orbit size = %ld\n", 
			i, step, dc_idx, ol);
		orbit_size->m_ii(i, ol);
		}
	{
	INT x_max = DC_X_PIX * width_in_pages;
	INT y_max = DC_Y_PIX 
		- (INT)((k + 1) * DC_GRAPH_RAD * 2.6) // space for graphs 
		- (INT)(DC_GRAPH_RAD * 3.); // space for legende 
	
	printf("calling vbp(): k = %ld x_max = %ld y_max = %ld\n", k, x_max, y_max);

	vbp(nl, orbit_size, x_max, y_max, plaz, o_dx, TRUE);
	}
	p->plaz = plaz;
	p->o_dx = o_dx;
	p->extrema[0] = 0.;
	p->extrema[1] = (double) DC_X_PIX * width_in_pages;
	p->extrema[2] = 0.;
	p->extrema[3] = (double) DC_Y_PIX * width_in_pages;
	p->extrema[4] = -1.;
	p->extrema[5] = 1.;
	
	/* enlarge the CO-system to shorten the graphics output ! */
	
	{
	INT *plazierung;
	INT x, y, la;
	
	plazierung = (INT *) p->plaz->ob_self.ob_charpointer;
	for (i = 0; i < p->nb_d; i++) {
		x = plazierung[2 * i];
		y = plazierung[2 * i + 1];
		dc_dc_no_to_dc_idx(dc0, i, &step, &dc_idx);
		la = step;
		// la = layer.s_ii(i);
		// y += DC_GRAPH_RAD * 3; // space for legende
		x += DC_GRAPH_RAD * 4 * width_in_pages; // space for legende
		y += (INT)(DC_GRAPH_RAD * 1.3);
		y += la * (INT)(DC_GRAPH_RAD * 2.6);
		plazierung[2 * i] = x;
		plazierung[2 * i + 1] = y;
		}
	}

	nb_d = 1;
	p->d_layer[0] = 0;
	p->d_layer_idx[0] = 0;
	for (i = 1; i <= up_to_step; i++) {
		if (ODD(i)) {
			L = dc0 + i;
			k = L->s_D()->s_li();
			for (j = 0; j < k; j++) {
				p->d_layer[nb_d + j] = i;
				p->d_layer_idx[nb_d + j] = j;
				}
			nb_d += k;
			}
		}
	if (nb_d != p->nb_d) {
		printf("open_dcl(): nb_d != p->nb_d");
		return NIL;
		}
	
	draw_PS_file(p, &p->vdev, DC_X_PIX, DC_Y_PIX, 
		file_name, dcl_draw_func);

	return p;
}

#if TEXDOCU
INT free_dcl(DC_LATTICE *p)
#endif
{
	INT i, l;
	
	if (p == NIL)
		return OK;
	if (p->dc) {
		l = p->up_to_step;
		for (i = 0; i <= l; i++) {
			p->dc[i].freeself();
			}
		// my_free(p->dc);
		p->dc = NIL;
		}
	if (p->plaz) {
		freeall(p->plaz);
		p->plaz = NIL;
		}
	if (p->o_dx) {
		freeall(p->o_dx);
		p->o_dx = NIL;
		}
	if (p->Ainf_t) {
		freeall(p->Ainf_t);
		p->Ainf_t = NIL;
		}
	if (p->nl) {
		freeall(p->nl);
		p->nl = NIL;
		}
	if (p->orbit_size) {
		freeall(p->orbit_size);
		p->orbit_size = NIL;
		}
	my_free(p);
	return OK;
}

#define DRAW_LEGENDE

#if TEXDOCU
INT dcl_draw_func(void *dc_lattice, INT x, INT y, INT w, INT h, VDEVICE *vdev)
#endif
{
	DC_LATTICE *p = (DC_LATTICE *) dc_lattice;
	INT i, j, nachbar, a_ij;
	INT n, graph_n, graph_n2;
	double x1, y1, z1;
	double x2, y2, z2;
	INT pix_x0, pix_y0;
	INT pix_x1, pix_y1;
	INT pix_x2, pix_y2;
	INT dx, dy;
	BYTE str[256];
	INT gra_x[DC_MAX_GRAPH_N];
	INT gra_y[DC_MAX_GRAPH_N];
	INT *plazierung, *o_dx;
	INT width_in_pages;
	INT at_x, at_y, odx;
	
	Vswr_mode(vdev, 2);
	plazierung = (INT *) p->plaz->ob_self.ob_charpointer;
	o_dx = (INT *) p->o_dx->ob_self.ob_charpointer;
	n = p->nl->s_ii(0) - 2;
	graph_n = p->n;
	graph_n2 = p->n2;
	if (!dcl_calc_xy(graph_n, gra_x, gra_y)) {
		Srff("dcl_draw_func", "dcl_calc_xy");
		return(FALSE);
		}
	width_in_pages = (INT)(p->extrema[1] / DC_X_PIX);
	printf("width_in_pages = %ld\n", width_in_pages);
	fflush(stdout);

#ifdef DRAW_LEGENDE
	// draw legende:
	dcl_draw_legende(p, width_in_pages, gra_x, gra_y, vdev);
#endif

	for (i = 0; i < p->nl->s_ii(0) - 1; i++) {
		at_x = plazierung[2 * i];
		at_y = plazierung[2 * i + 1];
		odx = o_dx[i];
		dcl_draw_orbit(p, i, at_x, at_y, odx, 
			width_in_pages, gra_x, gra_y, vdev);

		x1 = (double) plazierung[2 * i];
		y1 = (double) plazierung[2 * i + 1] + DC_GRAPH_RAD * 2.5;
		z1 = 0.;
		if (!user2vdev(vdev, x1, y1, z1, 
			&pix_x0, &pix_y0, p->extrema)) {
			Srff("dcl_draw_func", "user2vdev");
			return(FALSE);
			}
		for (j = p->nl->s_ii(i); 
			j < p->nl->s_ii(i + 1); j++) {
			nachbar = p->nl->s_ii(j);
			x2 = (double) plazierung[2 * nachbar];
			y2 = (double) plazierung[2 * nachbar + 1] - DC_GRAPH_RAD * 2.5;
			z2 = 0.;
			if (!user2vdev(vdev, x2, y2, z2, 
				&pix_x1, &pix_y1, p->extrema)) {
				Srff("dcl_draw_func", "user2vdev");
				return(FALSE);
				}
			dcl_line(vdev, pix_x0, pix_y0, 
				pix_x1, pix_y1);
			a_ij = p->Ainf_t->s_iji(i, nachbar);
			sprintf(str, "%ld", a_ij);
			dx = pix_x1 - pix_x0;
			dy = pix_y1 - pix_y0;
#if FALSE
			norm_vector = sqrt(dx * dx + dy * dy);
			pix_x3 = pix_x2 + (LONG)((double)dx * 
				DC_GRAPH_OFFSET_AIJ / norm_vector);
			pix_y3 = pix_y2 + (LONG)((double)dy * 
				DC_GRAPH_OFFSET_AIJ / norm_vector);
#endif
			pix_x2 = pix_x0 + (INT) ((double)dx * 0.3);
			pix_y2 = pix_y0 + (INT) ((double)dy * 0.3);
			V_gtext(vdev, pix_x2, pix_y2, str);
			}
		}
	return(TRUE);
}

#if TEXDOCU
INT dcl_draw_orbit(DC_LATTICE *p, INT orbit_idx, INT at_x, INT at_y, INT dx, 
	INT width_in_pages, INT *gra_x, INT *gra_y, VDEVICE *vdev)
#endif
{
	double ddx;
	double x1, y1, z1;
	INT pix_x0, pix_y0;
	BYTE str[256];
	INT i, l, dc_k;
	INT layer, layer_idx;
	DCY_OP Li;
	PERMUTATION_OB d;
	VECTOR_OB R, O;
	VECTOR_OP A, pR;
	
	A = (VECTOR_OP) p->dc->s_Ad_i(0);
	printf("orbit_idx = %ld\n", orbit_idx); fflush(stdout);
	x1 = (double) at_x;
	y1 = (double) at_y;
	z1 = 0.;
	if (!user2vdev(vdev, x1, y1, z1, 
		&pix_x0, &pix_y0, p->extrema)) {
		Srff("dcl_draw_orbit", "user2vdev");
		return(FALSE);
		}
	layer = p->d_layer[orbit_idx];
	layer_idx = p->d_layer_idx[orbit_idx];
	// printf("layer = %ld\n", layer); fflush(stdout);
	// printf("layer_idx = %ld\n", layer_idx); fflush(stdout);
	Li = p->dc + layer;
	Li->s_D_i(layer_idx)->copy(&d);
	dc_k = Li->s_k_i();
	// printf("dc_k = %ld\n", dc_k); fflush(stdout);
	if (layer == 1)
		dc_k++;
	sprintf(str, "%ld/%ld ", layer, layer_idx);
	d.sprint(str + strlen(str));
	if (p->f_with_perm)
		;
	else
		str[0] = 0;
	// dcl_3text(vdev, pix_x0, pix_y0, str, "", "");
	dc_get_k_set(&d, &R, dc_k, FALSE /* f_v */);
	R.println(); fflush(stdout);
	dc_calc_set_orbit(A, &R, &O, TRUE /* f_v */, FALSE /* f_vv */);
	l = O.s_li();
	ddx = (double)dx / l;
	ddx *= .7;
	
	for (i = 0; i < l; i++) {
		x1 = (double) at_x - (l - 1) * ddx * .5 + i * ddx;
		y1 = (double) at_y;
		z1 = 0.;
		if (!user2vdev(vdev, x1, y1, z1, 
			&pix_x0, &pix_y0, p->extrema)) {
			Srff("dcl_draw_orbit", "user2vdev");
			return(FALSE);
			}
	
		pR = (VECTOR_OP) O.s_i(i);
		dcl_print_graph(p, pR, pix_x0, pix_y0, 
			width_in_pages, gra_x, gra_y, vdev);
		}

	
	return OK;
}

#if TEXDOCU
INT dcl_print_graph(DC_LATTICE *p, VECTOR_OP R, INT pix_x, INT pix_y, 
	INT width_in_pages, INT *gra_x, INT *gra_y, VDEVICE *vdev)
#endif
{
	INT k, i, l, i1, j1;
	INT pix_x2, pix_y2;
	INT pix_x3, pix_y3;
	
	// draw the points of the graph:
	// pix_x1 = pix_x0 + DC_GRAPH_OFFSET_X;
	// pix_y1 = pix_y0;
	for (k = 0; k < p->n; k++) {
		pix_x2 = pix_x + gra_x[k];
		pix_y2 = pix_y + gra_y[k];
		dcl_halb(vdev, pix_x2, pix_y2, "");
		}
	
	// draw the edges of the graph:
	l = R->s_li();
	for (i = 0; i < l; i++) {
		k = R->s_ii(i);
		k2ij(k, &i1, &j1, p->n);
		pix_x2 = pix_x + gra_x[i1];
		pix_y2 = pix_y + gra_y[i1];
		pix_x3 = pix_x + gra_x[j1];
		pix_y3 = pix_y + gra_y[j1];
		dcl_line(vdev, 
			pix_x2, pix_y2, 
			pix_x3, pix_y3);
		}
	return OK;
}

#if TEXDOCU
INT dcl_draw_legende(DC_LATTICE *p, 
	INT width_in_pages, INT *gra_x, INT *gra_y, VDEVICE *vdev)
#endif
{
	double x1, y1, z1;
	INT pix_x0, pix_y0;
	INT pix_x1, pix_y1;
	INT pix_x2, pix_y2;
	INT pix_x3, pix_y3;
	BYTE str[256];
	INT kk, i, k, j;
	
	x1 = DC_GRAPH_RAD * 3 * width_in_pages; // x_max
	y1 = DC_GRAPH_RAD * 3; // y_max
	z1 = 0.;
	kk = 0;
	for (i = 0; i < p->n - 1; i++) {
		if (!user2vdev(vdev, x1, y1, z1, 
			&pix_x0, &pix_y0, p->extrema)) {
			Srff("dcl_draw_legende", "user2vdev");
			return(FALSE);
			}
		pix_x1 = pix_x0 + DC_GRAPH_OFFSET_X;
		pix_y1 = pix_y0;
		for (k = 0; k < p->n; k++) {
			pix_x2 = pix_x1 + gra_x[k];
			pix_y2 = pix_y1 + gra_y[k];
			dcl_halb(vdev, pix_x2, pix_y2, "");
			}
		str[0] = 0;
		for (j = i + 1; j < p->n; j++) {
			pix_x2 = pix_x1 + gra_x[i];
			pix_y2 = pix_y1 + gra_y[i];
			pix_x3 = pix_x1 + gra_x[j];
			pix_y3 = pix_y1 + gra_y[j];
			dcl_line(vdev, 
				pix_x2, pix_y2, 
				pix_x3, pix_y3);	
			sprintf(str + strlen(str), "%ld", kk);
			if (j < p->n - 1)
				strcat(str, ",");
#if 0
			dx = pix_x3 - pix_x2;
			dy = pix_y3 - pix_y2;
			pix_x3 = pix_x2 + (INT) ((double)dx * 0.3);
			pix_y3 = pix_y2 + (INT) ((double)dy * 0.3);
			V_gtext(vdev, pix_x3, pix_y3, str);
#endif
			kk++;
			} // next j
		pix_x2 = pix_x1 + gra_x[0];
		pix_y2 = pix_y1 + gra_y[0];
		V_gtext(vdev, pix_x2, pix_y2, str);
		
		y1 += DC_GRAPH_RAD * 4;
		}
	return OK;
}

#if TEXDOCU
INT dcl_calc_xy(INT n, INT *gra_x, INT *gra_y)
#endif
{
	INT l;
	double phi;
	
	phi = M_2PI / (double) n;
	for (l = 0; l < n; l++) {
		gra_x[l] = (INT)(cos(phi * (double)l) * 
			DC_GRAPH_RAD);
		gra_y[l] = (INT)(sin(phi * (double)l) * 
			DC_GRAPH_RAD);
		}
	return TRUE;
}

#if TEXDOCU
void dcl_3text(VDEVICE *vdev, INT x, INT y, 
	BYTE *text1, BYTE *text2, BYTE *text3)
#endif
{
	INT rad, off_x, off_y1, off_y2, off_y3;
	
#ifdef SYSTEMUNIX
	rad = 5;
	off_x = 15;
	off_y1 = - 25;
	off_y2 = - 5;
	off_y3 = 15;
#else
	rad = 4;
	off_x = 5;
	off_y1 = -12;
	off_y2 = -2;
	off_y2 = 8;
#endif
	V_pie(vdev, x, y, rad, 0, 3600);
	V_gtext(vdev, x + off_x, y + off_y1, text1);
	V_gtext(vdev, x + off_x, y + off_y2, text2);
	V_gtext(vdev, x + off_x, y + off_y3, text3);
}

#if TEXDOCU
void dcl_halb(VDEVICE *vdev, INT x, INT y, BYTE *text)
#endif
{
	INT rad, off_x, off_y;
	
#ifdef SYSTEMUNIX
	rad = 2;
	off_x = 15;
	off_y = 20;
#else
	rad = 2;
	off_x = 5;
	off_y = 8;
#endif
	V_pie(vdev, x, y, rad, 0, 3600);
	V_gtext(vdev, x + off_x, y + off_y, text);
}

#if TEXDOCU
INT dcl_line(VDEVICE *vdev, INT x0, INT y0, INT x1, INT y1)
#endif
{
	INT pxy[4];
	
	pxy[0] = x0;
	pxy[1] = y0;
	pxy[2] = x1;
	pxy[3] = y1;
	V_pline(vdev, 2, pxy);
	return TRUE;
}

#endif /* LADDER_TRUE */

/* end of dc_draw.C */
