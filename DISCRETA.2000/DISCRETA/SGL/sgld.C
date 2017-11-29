/* sgld.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SGL_TRUE

#include <DISCRETA/ma.h>
#include <DISCRETA/perm.h>
#include <DISCRETA/dimino.h>
#include <DISCRETA/divs.h>
#include <DISCRETA/sgl.h>
#include <DISCRETA/lb.h>

INT sgll_draw_sims = TRUE;

INT sgll_draw_labra(void *sgll, 
	INT x, INT y, INT w, INT h, VDEVICE *vdev)
{
	SG_LATTICE_LOCAL *p = (SG_LATTICE_LOCAL *) sgll;
	VECTOR_OP Gen, gen;
	LABRA_OB LA;
	INT i, l, deg;
	INT f_v = FALSE;
	INT f_vv = FALSE;
	PERMUTATION_OP perm;
	double x1, y1, z1;
	double x2, y2, z2;
	double x3, y3, z3;
	INT pix_x1, pix_y1;
	INT pix_x2, pix_y2;
	INT pix_x3, pix_y3;
	BYTE str[1024];
	INT j, k, prev;
	INT unit_x = 40000;
	INT unit_y = 30000;
	INT x0 = 0, y0 = vdev->co[3];
	
	printf("vdev->co: %ld %ld %ld %ld\n", vdev->co[0], vdev->co[1], vdev->co[2], vdev->co[3]);
	printf("y0 = %ld\n", y0);
	Gen = p->generators;
	l = Gen->s_li();
	for (i = 0; i < l; i++) {
		printf("Orbit %ld: generators:\n", i);
		fflush(stdout);
		gen = (VECTOR_OP) Gen->s_i(i);
		deg = 0;
		if (gen->s_li()) {
			perm = (PERMUTATION_OP) gen->s_i(0);
			deg = perm->s_li();
			LA.Init(deg, gen, gen->s_li(), f_v, f_vv);
			LA.jerrum(FALSE);
			// LA.Print();
			// LA.Print_T(f_v);
			// LA.print_group_order();
			for (j = 0; j < deg; j++) {
				if (sgll_draw_sims) {
					for (k = 0; k < deg; k++) {
						if (LA.path_exists(j, k)) {
							x1 = (double) (x0 + j * unit_x);
							y1 = (double) (y0 - j * unit_y);
							z1 = 0.;
							if (!user2vdev(vdev, x1, y1, z1, 
								&pix_x1, &pix_y1, p->extrema)) {
								Srff("sgll_draw_func", "user2vdev");
								return(FALSE);
								}
							x1 = (double) (x0 + k * unit_x);
							y1 = (double) y0 - (double)((j + 1) * unit_y) * 0.8;
							z1 = 0.;
							if (!user2vdev(vdev, x1, y1, z1, 
								&pix_x2, &pix_y2, p->extrema)) {
								Srff("sgll_draw_func", "user2vdev");
								return(FALSE);
								}
							bind_mol(vdev, pix_x1, pix_y1, pix_x2, pix_y2, 1);
							}
						}
					}
				else { // if (!sgll_draw_sims) 
					prev = LA.s_V_ii(j);
					perm = LA.s_KM_i(j);
					str[0] = 0;
					perm->sprint(str);
					if (prev != j) {
						x1 = (double) (x0 + (prev + 1) * unit_x);
						y1 = (double) (y0 - (prev + 1) * unit_y);
						z1 = 0.;
						// printf("y0 = %ld prev = %ld unit_y = %ld y1 = %lf\n", y0, prev, unit_y, y1);
						if (!user2vdev(vdev, x1, y1, z1, 
							&pix_x1, &pix_y1, p->extrema)) {
							Srff("sgll_draw_func", "user2vdev");
							return(FALSE);
							}
						x1 = (double) (x0 + (j + 1) * unit_x);
						y1 = (double) (y0 - (prev + 1 + 1) * unit_y);
						z1 = 0.;
						if (!user2vdev(vdev, x1, y1, z1, 
							&pix_x2, &pix_y2, p->extrema)) {
							Srff("sgll_draw_func", "user2vdev");
							return(FALSE);
							}
						bind_mol(vdev, pix_x1, pix_y1, pix_x2, pix_y2, 1);
						}
					else {
						x1 = (double) (x0 + (j + 1) * unit_x);
						y1 = (double) (y0 - (j + 1) * unit_y);
						z1 = 0.;
						if (!user2vdev(vdev, x1, y1, z1, 
							&pix_x1, &pix_y1, p->extrema)) {
							Srff("sgll_draw_func", "user2vdev");
							return(FALSE);
							}
						x2 = (double) (x0 + (j + 1) * unit_x);
						y2 = (double) (y0 - (j + 1 + 1) * unit_y);
						z2 = 0.;
						if (!user2vdev(vdev, x2, y2, z2, 
							&pix_x2, &pix_y2, p->extrema)) {
							Srff("sgll_draw_func", "user2vdev");
							return(FALSE);
							}
						bind_mol(vdev, pix_x1, pix_y1, pix_x2, pix_y2, 1);
						}
					} // else sgll_draw_sims
				if (!sgll_draw_sims) {
					x3 = (x1 + x2) / 2;
					y3 = (y1 + y2) / 2;
					y3 -= (double) j * (unit_y / 3);
					z3 = 0.;
					if (!user2vdev(vdev, x3, y3, z3, 
						&pix_x3, &pix_y3, p->extrema)) {
						Srff("sgll_draw_func", "user2vdev");
						return(FALSE);
						}
					draw_mol(vdev, pix_x3, pix_y3, str);
					}
				}
			}
		printf("\n");
		x0 += unit_x * deg;
		x0 += unit_x * 1;
		if (i % 4 == 0) {
			x0 = 0;
			y0 -= unit_y * deg;
			y0 -= unit_y * 1;
			}
		} /* next i */
	return OK;
}

static INT sgl_Line(VDEVICE *vdev, INT *x, INT *y, INT a, INT b) 
{
	Bind_mol(vdev, 1, x[a], y[a], x[b], y[b], 1);
	return OK;
}


INT sgll_draw_boxed_func(void *sgll, 
	INT x, INT y, INT w, INT h, VDEVICE *vdev)
{
	SG_LATTICE_LOCAL *p = (SG_LATTICE_LOCAL *) sgll;
	INT Px[20], Py[20];
	INT px[20], py[20];
	INT i;
	
	sgll_draw_func(sgll, x, y, w, h, vdev);
	Px[0] = 0;
	Py[0] = 0;
	Px[1] = SGL_VBP_X_PIX;
	Py[1] = 0;
	Px[2] = SGL_VBP_X_PIX;
	Py[2] = SGL_VBP_Y_PIX;
	Px[3] = 0;
	Py[3] = SGL_VBP_Y_PIX;
	for (i = 0; i <= 3; i++) {
		user2vdev(vdev, Px[i], Py[i], 0, &px[i], &py[i], p->extrema);
		}
	sgl_Line(vdev, px, py, 0, 1);
	sgl_Line(vdev, px, py, 1, 2);
	sgl_Line(vdev, px, py, 2, 3);
	sgl_Line(vdev, px, py, 3, 0);

	return TRUE;
}

INT sgll_draw_func(void *sgll, 
	INT x, INT y, INT w, INT h, VDEVICE *vdev)
{
	SG_LATTICE_LOCAL *p = (SG_LATTICE_LOCAL *) sgll;
	INT i, j, nachbar;
	INT n;
	INT *plazierung;
	double x1, y1, z1;
	INT ix0, iy0;
	INT ix1, iy1;
	INT pix_x1, pix_y1;
	INT rad;
	
	Vswr_mode(vdev, 2);
	n = p->nl->s_ii(0) - 2;
	plazierung = (INT *) p->plaz->ob_self.ob_charpointer;


	user2vdev(vdev, 0, 0, 0, &ix0, &iy0, p->extrema);
	user2vdev(vdev, SGL_RADIUS, 0, 0, &ix1, &iy1, p->extrema);
	rad = ABS(ix1 - ix0);
	printf("rad = %ld\n", rad);

	for (i = 0; i < p->nl->s_ii(0) - 1; i++) {
		for (j = p->nl->s_ii(i); 
			j < p->nl->s_ii(i + 1); j++) {
			nachbar = p->nl->s_ii(j);
			sgll_d_neighbour(p, i, nachbar, vdev);
			}
		}
	
	for (i = 0; i < p->nl->s_ii(0) - 1; i++) {
		x1 = (double) plazierung[2 * i];
		y1 = (double) plazierung[2 * i + 1];
		z1 = 0.;
		if (!user2vdev(vdev, x1, y1, z1, 
			&pix_x1, &pix_y1, p->extrema)) {
			Srff("sgll_draw_func", "user2vdev");
			return(FALSE);
			}
		sgll_d_orbit(p, i, rad, vdev);
		}
	return(TRUE);
}

INT sgll_d_orbit(SG_LATTICE_LOCAL *sgll, INT orbit, INT rad, VDEVICE *vdev)
{
	double orbit_dx, dx1, dx1_halbe;
	double x1, y1, z1, x1o, x1oo;
	INT pix_x0, pix_y0, pix_x1, pix_y1;
	INT j, l;
	INT *plazierung;
	INT *o_dx;
	INT ix1, iy1;
	
	plazierung = (INT *) sgll->plaz->ob_self.ob_charpointer;
	o_dx = (INT *) sgll->o_dx->ob_self.ob_charpointer;
	x1 = (double) plazierung[2 * orbit];
	y1 = (double) plazierung[2 * orbit + 1];
	z1 = 0.;
	if (sgll->draw_lines_type == 3) {
		sgll_get_ko_center(sgll, orbit, &ix1, &iy1, vdev);
		sgll_draw_dot(vdev, &ix1, &iy1, 0, rad);
		}
	else {
		l = sgll->orbit_size->s_ii(orbit);
		if (l > 1) {
			sgll_get_ko(sgll, orbit, 0, &pix_x0, &pix_y0, vdev);
			sgll_get_ko(sgll, orbit, l - 1, &pix_x1, &pix_y1, vdev);
			bind_mol(vdev, pix_x0, pix_y0, pix_x1, pix_y1, 1);
			}
		orbit_dx = o_dx[orbit] * 0.7;
		dx1 = orbit_dx / (double) sgll->orbit_size->s_ii(orbit);
		dx1_halbe = dx1 * .5;
		x1o = x1 - orbit_dx * .5;
		for (j = 0; j < l; j++) {
			x1oo = x1o + dx1_halbe;
			if (!user2vdev(vdev, x1oo, y1, z1, 
				&pix_x1, &pix_y1, sgll->extrema)) {
				Srff("sgll_d_orbit", "user2vdev");
				return FALSE;
				}
			if (j == 0) {
				pix_x0 = pix_x1;
				pix_y0 = pix_y1;
				}
			sgll_draw_dot(vdev, &pix_x1, &pix_y1, 0, rad);
			// draw_mol(vdev, pix_x1, pix_y1, "");
			x1o += dx1;
			}
		}
	return TRUE;
}

INT sgll_d_neighbour(SG_LATTICE_LOCAL *sgll, 
	INT orbit, INT neighbour, VDEVICE *vdev)
{
	INT i, j, l1, o1, l2, o2;
	INT x1, x2, y1, y2;
	
	sgll->L->orbit2lo(orbit, &l1, &o1);
	sgll->L->orbit2lo(neighbour, &l2, &o2);
	if (sgll->draw_lines_type == 0) { /* Asup */
		j = 0;
		for (i = 0; 
			i < sgll->orbit_size->s_ii(orbit); 
			i++) {
			sgll_d_neighbour2(sgll, 
				orbit, neighbour, 
				l1, o1, i, l2, o2, j, vdev);
			}
		}
	else if (sgll->draw_lines_type == 1) { /* Ainf */
		i = 0;
		for (j = 0; 
			j < sgll->orbit_size->s_ii(neighbour); 
			j++) {
			sgll_d_neighbour2(sgll, 
				orbit, neighbour, 
				l1, o1, i, l2, o2, j, vdev);
			}
		}
	else if (sgll->draw_lines_type == 2) { /* all */
		for (i = 0; 
			i < sgll->orbit_size->s_ii(orbit); 
			i++) {
			for (j = 0; 
				j < sgll->orbit_size->s_ii(neighbour); 
				j++) {
				sgll_d_neighbour2(sgll, 
					orbit, neighbour, 
					l1, o1, i, l2, o2, j, vdev);
				}
			}
		}
	else if (sgll->draw_lines_type == 3) { // only poset of orbits
		sgll_get_ko_center(sgll, orbit, &x1, &y1, vdev);
		sgll_get_ko_center(sgll, neighbour, &x2, &y2, vdev);
		Bind_mol(vdev, 1, x1, y1, x2, y2, 1);
		}
	return TRUE;
}

INT sgll_d_neighbour2(SG_LATTICE_LOCAL *sgll, 
	INT orbit, INT neighbour, 
	INT l1, INT o1, INT rep1, 
	INT l2, INT o2, INT rep2, 
	VDEVICE *vdev)
{
	INT f_is_subgroup;
	INT pix_x0, pix_y0, pix_x1, pix_y1;
	
	if (sgll->L->IsSubgroup(l1, o1, rep1, 
		l2, o2, rep2, &f_is_subgroup) != OK) {
		Srff("sgll_d_neighbour2", 
		"SG_LATTICE::IsSubgroup");
		return FALSE;
		}
	if (f_is_subgroup) {
		sgll_get_ko(sgll, orbit, rep1, 
			&pix_x0, &pix_y0, vdev);
		sgll_get_ko(sgll, neighbour, rep2, 
			&pix_x1, &pix_y1, vdev);
		bind_mol(vdev, pix_x0, 
			pix_y0, pix_x1, pix_y1, 1);
		}
	return OK;
}

INT sgll_get_ko_center(SG_LATTICE_LOCAL *sgll, 
	INT orbit, INT *pix_x1, INT *pix_y1, VDEVICE *vdev)
{
	double x1, y1, z1;
	INT *plazierung;
	INT *o_dx;
	
	plazierung = (INT *) sgll->plaz->ob_self.ob_charpointer;
	o_dx = (INT *) sgll->o_dx->ob_self.ob_charpointer;
	x1 = (double) plazierung[2 * orbit];
	y1 = (double) plazierung[2 * orbit + 1];
	z1 = 0.;
	if (!user2vdev(vdev, x1, y1, z1, 
		pix_x1, pix_y1, sgll->extrema)) {
		Srff("sgll_get_ko_center", "user2vdev");
		return FALSE;
		}
	return TRUE;
}

INT sgll_get_ko(SG_LATTICE_LOCAL *sgll, 
	INT orbit, INT rep, 
	INT *pix_x1, INT *pix_y1, VDEVICE *vdev)
{
	double orbit_dx, dx1, dx1_halbe;
	double x1, y1, z1, x1o, x1oo;
	INT *plazierung;
	INT *o_dx;
	
	plazierung = (INT *) sgll->plaz->ob_self.ob_charpointer;
	o_dx = (INT *) sgll->o_dx->ob_self.ob_charpointer;
	x1 = (double) plazierung[2 * orbit];
	y1 = (double) plazierung[2 * orbit + 1];
	z1 = 0.;
	orbit_dx = o_dx[orbit] * 0.7;
	dx1 = orbit_dx / (double) sgll->orbit_size->s_ii(orbit);
	dx1_halbe = dx1 * .5;
	x1o = x1 - orbit_dx * .5;
	x1o += rep * dx1;
	x1oo = x1o + dx1_halbe;
	if (!user2vdev(vdev, x1oo, y1, z1, 
		pix_x1, pix_y1, sgll->extrema)) {
		Srff("sgll_get_ko", "user2vdev");
		return FALSE;
		}
	return TRUE;
}

#define SGLL_GREY_INTERIOR 70
#define SGLL_GREY_COLOR 0

INT sgll_draw_dot(VDEVICE *vdev, INT *x, INT *y, INT i, INT rad)
{
	Vsf_interior(vdev, SGLL_GREY_INTERIOR);
	Vsf_color(vdev, SGLL_GREY_COLOR);
	V_pie(vdev, x[i], y[i], rad, 0, 0);
	Vsf_interior(vdev, 0);
	V_arc(vdev, x[i], y[i], rad, 0 /* begang */, 3600 /* endang */);
	return OK;
}

#if 0
INT gpl_d_generators(GP_LATTICE *gpl, 
	INT orbit, VDEVICE *vdev)
{
	INT l1, o1, index_in_normalizer, k, s_int;
	FGO_OP fgo;
	FGG_OP H;
	SYM_OP s1; /* in G1 */
	SYM_OP s; /* in G */
	BYTE str1[256];
	BYTE str_tmp[256];
	INT pix_x0, pix_y0;
	
	gpl->L->orbit2lo(orbit, &l1, &o1);
	gpl->L->IndexInNormalizer(l1, o1, 
		&index_in_normalizer);
	fgo = gpl->L->s_theOrbits_ij(l1, o1);
	H = fgo->s_H();
	strcpy(str1, "|<");
	for (k = 0; k < H->s_Ssize_i(); k++) {
		s1 = H->s_S_i(k);
		if (gpl->f_enumerated_group) {
			s_int = ((INTEGER_OP)s1)->s_i();
			s = gpl->G->s_Hsorted_i(s_int);
			}
		else {
			s = s1;
			}
		str_tmp[0] = 0;
		s->sprint(str_tmp);
		strcat(str1, str_tmp);
		if (k < H->s_Ssize_i() - 1)
			strcat(str1, ", ");
		else {
			sprintf(str_tmp, ">| = %ld; %ld", 
				H->s_Hsize_i(), index_in_normalizer);
			strcat(str1, str_tmp);
			}
		gpl_get_ko(gpl, orbit, 0, 
			&pix_x0, &pix_y0, vdev);
		V_gtext(vdev, 
			pix_x0, pix_y0 + 20 + k * 10, str1);
		str1[0] = 0;
		}
	return TRUE;
}
#endif

#endif /* SGL_TRUE */

