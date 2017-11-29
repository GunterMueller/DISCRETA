// t140e.C
// (modified t140.C)
//
// 26.10.1999
//
// draws a picture of the subgroup lattice
// for the isomorphism classification according to the 
// 8-design paper
//
// Evi Haberberger
//

#include <DISCRETA/discreta.h>
#include <DISCRETA/graphics.h>

#include <stdio.h>		// printf
#include <string.h>		// memset
#include <time.h>		// clock

#include <stdlib.h>

typedef struct draw_local DRAW_LOCAL;

struct draw_local 
{
	double extrema[6];
	BYTE *finput;
	BYTE *fname;
	INT pic;
	INT draw_lines_type;	// 0 = Asup
				// 1 = Ainf
				// 2 = all
	MATRIX_OP info, incma, Ainf, Acover;
	VECTOR_OP nl, orbit_size, basics;
	VECTOR_OP Px, Py, Os, x, y;
	INT nb_con, nb_gr;
};

#define BUFSIZE 50000

#define RADIUS (6)
#define RADIUS1 (180)

#define MAX_POINTS 500

#define DASH 30
#define EPSILON 0.0001

#define FNAME_MASK "%s.%s"
#define FACTOR 0.9


INT draw_pic(void *draw_data, INT x, INT y, INT w, INT h, VDEVICE *vdev);
INT draw_orbit(VDEVICE *vdev, DRAW_LOCAL *p, INT *Px, INT *Py, INT *x, INT *y, INT orbit, INT rad);
INT draw_neighbour(VDEVICE *vdev, DRAW_LOCAL *p, INT *Px, INT *Py, INT *x, INT *y, 
	INT orbit, INT neighbour);
INT draw_neighbour2(VDEVICE *vdev, DRAW_LOCAL *p, INT *Px, INT *Py, 
	INT o1, INT rep1, INT o2, INT rep2);
INT get_ko(VDEVICE *vdev, DRAW_LOCAL *p, 
	INT orbit, INT rep, 
	INT *Px, INT *Py, 
	INT *pix_x1, INT *pix_y1);
INT lincomb3(INT *Px, INT *Py, INT a1, INT i1, INT a2, INT i2, INT a3, INT i3, INT idx);
INT find_midpoint(INT *Px, INT *Py, INT i1, INT i2, INT idx);
INT insert_lattice(DRAW_LOCAL *latt);

INT do_it(BYTE *finput_base, BYTE *fname_base, INT pic)
{
	BYTE fname[1024];
	BYTE finput[1024];
	VDEVICE vdev;
	DRAW_LOCAL *data = NIL;

	data = (DRAW_LOCAL *) my_malloc(sizeof(DRAW_LOCAL), "do_it");
	data->extrema[0] = -500;
	data->extrema[1] = 500;
	data->extrema[2] = 0;
	data->extrema[3] = 1000;
	data->extrema[4] = 0;
	data->extrema[5] = 0;
	data->finput = finput;
	data->fname = fname;
	data->pic = pic;
	data->draw_lines_type = 2;
	
		
	sprintf(fname, FNAME_MASK, fname_base, "mp");
	sprintf(finput, FNAME_MASK, finput_base, "txt");

	VECTOR_OB orbit_size, nl, Os, Px, Py, basics;
	MATRIX_OB info, incma, Ainf, Acover;

	data->info = &info;
	data->incma = &incma;
	data->Ainf = &Ainf;
	data->Acover = &Acover;
	data->nl = &nl;
	data->orbit_size = &orbit_size;
	data->basics = &basics;
		
	data->Px = &Px;
	data->Py = &Py;
	data->Os = &Os;

	insert_lattice(data);
		
	printf("draw .mp:\n"); fflush(stdout);
	draw_mp_file(data, &vdev, fname, 1000 /* factor 1000 */, draw_pic);
	printf("finished.\n"); fflush(stdout);

	return OK;
}

int main(int argc, char **argv)
{
	INT t0, t1, user_time;

	if ( argc < 3 ) 
	{
		printf("usage: t140.out [options] finput fname pic \n");
		printf("draws a picture of the group lattice \n");
		printf("from finput.txt to fname.mp\n");
		printf("options:\n");
		return 1;
		}
	discreta_init();
	{
		// INT t0, t1, user_time;
		INT i;
		BYTE s[256];
		BYTE *fname;
		BYTE *finput;
		INT f_v = FALSE;
		INT pic;
		
		t0 = os_ticks();
		{
			for (i = 1; i < argc - 1; i++) 
			{
				if (strcmp(argv[i], "-v") == 0) 
				{
					f_v = TRUE;
				}
				else
					break;
			}
			finput = argv[i++];
			fname = argv[i++];
			sscanf(argv[i++], "%ld", &pic);
			do_it(finput, fname, pic);
		}
	
		t1 = os_ticks();
		user_time = t1 - t0;
		s[0] = 0;
		print_delta_time_100(user_time, s);
		printf("total computing time: %s\n", s);
		fflush(stdout);
	}
	discreta_exit();
	return 0;
}

INT draw_pic(void *draw_data, INT x_, INT y_, INT w_, INT h_, VDEVICE *vdev)
{
	DRAW_LOCAL *p = (DRAW_LOCAL *) draw_data;
	INT *Px, *Py, *Pz;
	INT x0, y0, x1, y1;
	INT *x, *y;
	INT *pts;
	INT rad, rad1;
	INT i, j;
	INT nb_nodes = 300;
	
	// printf("vdev->co: %ld %ld %ld %ld\n", vdev->co[0], vdev->co[1], vdev->co[2], vdev->co[3]);

	user2vdev(vdev, 0, 0, 0, &x0, &y0, p->extrema);
	user2vdev(vdev, RADIUS, 0, 0, &x1, &y1, p->extrema);
	rad = ABS(x1 - x0);
	printf("rad = %ld\n", rad);
	
	user2vdev(vdev, 0, 0, 0, &x0, &y0, p->extrema);
	user2vdev(vdev, RADIUS, 0, 0, &x1, &y1, p->extrema);
	rad1 = ABS(x1 - x0);
	printf("rad1 = %ld\n", rad1);
	
	Px = (INT *) my_malloc(nb_nodes * sizeof(INT), "draw_pic Px");
	Py = (INT *) my_malloc(nb_nodes * sizeof(INT), "draw_pic Py");
	Pz = (INT *) my_malloc(nb_nodes * sizeof(INT), "draw_pic Pz");
	x = (INT *) my_malloc(nb_nodes * sizeof(INT), "draw_pic x");
	y = (INT *) my_malloc(nb_nodes * sizeof(INT), "draw_pic y");
	pts = (INT *) my_malloc(nb_nodes * sizeof(INT), "draw_pic pts");
	
	Px[0] = 0;
	Py[0] = 0;
	for (i = 0; i < nb_nodes; i++) 
	{
		Pz[i] = 0;
	}
	
	on_circle_int(Px, Py, 1, 45, 200);
	on_circle_int(Px, Py, 2, 90, 200);
	on_circle_int(Px, Py, 3, 135, 200);
	
	lincomb3(Px, Py, 0, 1, 1, 2, 0, 3, 4);		/* P */
	lincomb3(Px, Py, 0, 1, 2, 2, 0, 3, 5);		/* N_A(P) */
	lincomb3(Px, Py, 1, 1, 2, 2, 0, 3, 6);		/* A */
	lincomb3(Px, Py, 0, 1, 2, 2, 1, 3, 7);		/* N_{N_G(A)}(P) */
	lincomb3(Px, Py, 1, 1, 2, 2, 1, 3, 8);		/* N_G(A) */
	lincomb3(Px, Py, 0, 1, 2, 2, 2, 3, 9);		/* N_G(P) */
	lincomb3(Px, Py, 1, 1, 2, 2, 2, 3, 10);		 
	lincomb3(Px, Py, 3, 1, 2, 2, 2, 3, 11);		/* G */
	
	find_midpoint(Px, Py, 0, 4, 12);
	find_midpoint(Px, Py, 4, 5, 13);
	find_midpoint(Px, Py, 5, 6, 14);
	find_midpoint(Px, Py, 5, 7, 15);
	find_midpoint(Px, Py, 7, 9, 16);
	
	printf("Py[6] = %ld\n", Py[6]);
	printf("Py[11] = %ld\n", Py[11]);
	
	for (i = 0; i < p->nb_con; i++) {
		Px[100 + i] = p->Px->s_ii(i) + Px[6];
		Py[100 + i] = p->Py->s_ii(i) + Py[6];
		}
	for (i = 0; i < p->nb_con; i++) {
		Px[200 + i] = p->Os->s_ii(i);
		Py[200 + i] = 0;
		}

	for (i = 0; i <= 16; i++) 
	{
		user2vdev(vdev, Px[i], Py[i], 0, &x[i], &y[i], p->extrema);
	}
	for (i = 0; i < p->nb_con; i++) 
	{
		user2vdev(vdev, Px[100 + i], Py[100 + i], 0, &x[100 + i], &y[100 + i], p->extrema);
	}
	for (i = 0; i < p->nb_con; i++) 
	{
		user2vdev(vdev, Px[200 + i], Py[200 + i], 0, &x[200 + i], &y[200 + i], p->extrema);
	}
	
	draw_line(vdev, x, y, 0, 5);
	draw_line(vdev, x, y, 5, 6);
	draw_line(vdev, x, y, 7, 8);
	draw_line(vdev, x, y, 9, 11);
	draw_line(vdev, x, y, 5, 9);
	draw_line(vdev, x, y, 6, 10);

	draw_dot(vdev, x, y, 0, rad);
	draw_dot(vdev, x, y, 4, rad);
	draw_dot(vdev, x, y, 5, rad);
	draw_dot(vdev, x, y, 6, rad);
	draw_dot(vdev, x, y, 7, rad);
	draw_dot(vdev, x, y, 8, rad);
	draw_dot(vdev, x, y, 9, rad);
	draw_dot(vdev, x, y, 11, rad);

	draw_label(vdev, "1", "tl", x, y, 0);
	draw_label(vdev, "$P$", "tl", x, y, 4);
	draw_label(vdev, "$N_A(P)$", "tl", x, y, 5);
	draw_label(vdev, "$N_{N_G(A)}(P)$", "r", x, y, 7);
	draw_label(vdev, "$N_G(A)$", "r", x, y, 8);
	draw_label(vdev, "$N_G(P)$", "r", x, y, 9);
	


	BYTE str[1000];
	
	for (i = 0; i < 5; i++)
	{
		sprintf(str, "%ld", p->basics->s_ii(i));
		draw_label(vdev, str, "l", x, y, 12+i);
	}

	if (p->Ainf->s_li()>2)
	{
		printf("draw_neighbours:\n"); fflush(stdout);
		for (i = 0; i < p->nb_con; i++) 
		{
			for (j = p->nl->s_ii(i); 
				j < p->nl->s_ii(i + 1); j++) 
			{
				INT nachbar = p->nl->s_ii(j);
				printf("draw_neighbour: %ld -> %ld\n", i, nachbar); 
				fflush(stdout);
				draw_neighbour(vdev, p, Px, Py, x, y, i, nachbar);
			}
		}
		for (i = 0; i < p->nb_con; i++) 
		{
			printf("draw_orbit %ld:\n", i); fflush(stdout);
			draw_orbit(vdev, p, Px, Py, x, y, i, rad);
		}
	
		Vst_boxed(vdev, FALSE);
		Vst_overwrite(vdev, FALSE);
#if 0
		for (i = (p->nb_con)-2; i > 0; i--)
		{
#endif			
		for (i = 1; i < (p->nb_con)-1; i++) 
		{
			sprintf(str, "$%ld \\times B_{%ld}$", 
				p->info->s_iji(i, 2), 
				i);
			get_ko(vdev, p, i, 0, 
				Px, Py, 
				&x0, &y0);
			// draw_label(vdev, str, "tr", &x0, &y0, 0);
			x0 = x[100 + i];
			y0 = y[100 + i] - 7 * rad;
			draw_label(vdev, str, "t", &x0, &y0, 0);
		}
	}	
	Vst_boxed(vdev, FALSE);
	Vst_overwrite(vdev, TRUE);

	sprintf(str, "$A$(%ld): %ld", 
			p->info->s_iji(0, 1),
			p->info->s_iji(0, 3));
	draw_label(vdev, str, "l", x, y, 6);

	sprintf(str, "$G$");
	draw_label(vdev, str, "l", x, y, 11);

#if 0
	for (i = 1; i < (p->nb_con)-1; i++) {
		draw_label(vdev, ((MATRIX_OP)p->info)->s_iji(i, 3), "br", x, y, 200+i);
		}
#endif

#if 0
	for (i = 0; i < p->nb_con; i++) {
		draw_dot(vdev, x, y, 100 + i, rad);
		}
#endif

	my_free(Px);
	my_free(Py);
	my_free(Pz);
	my_free(x);
	my_free(y);
	my_free(pts);

	return OK;
}

INT draw_rep(VDEVICE *vdev, DRAW_LOCAL *p, INT orbit, INT rep, INT x, INT y, INT rad)
{
	draw_dot(vdev, &x, &y, 0, rad);
	
	return OK;
}

INT draw_orbit(VDEVICE *vdev, DRAW_LOCAL *p, INT *Px, INT *Py, INT *x, INT *y, INT orbit, INT rad)
{
	double x1, y1, z1;
	INT pix_x0, pix_y0, pix_x1, pix_y1;
	INT j, l;
	INT ix1, iy1;
	
	x1 = (double) x[100 + orbit];
	y1 = (double) y[100 + orbit];
	z1 = 0.;
	if (p->draw_lines_type == 3) {
		ix1 = (INT) x1;
		iy1 = (INT) y1;
		draw_rep(vdev, p, orbit, 0, ix1, iy1, rad);
		}
	else {
		l = p->orbit_size->s_ii(orbit);
		if (l > 1) {
			get_ko(vdev, p, orbit, 0, Px, Py, &pix_x0, &pix_y0);
			get_ko(vdev, p, orbit, l - 1, Px, Py, &pix_x1, &pix_y1);
			bind_mol(vdev, pix_x0, pix_y0, pix_x1, pix_y1, 1);
			}
		for (j = 0; j < l; j++) {
			get_ko(vdev, p, orbit, j, Px, Py, &pix_x0, &pix_y0);
			// draw_mol(vdev, pix_x1, pix_y1, "");
			draw_rep(vdev, p, orbit, j, pix_x0, pix_y0, rad);
			}
		}
	return TRUE;
}

INT draw_neighbour(VDEVICE *vdev, DRAW_LOCAL *p, INT *Px, INT *Py, INT *x, INT *y, 
	INT orbit, INT neighbour)
{
	INT i, j, o1, o2;
	INT x1, x2, y1, y2;
	
	o1 = orbit;
	o2 = neighbour;
	if (p->draw_lines_type == 0) { /* Asup */
		j = 0;
		for (i = 0; 
			i < p->orbit_size->s_ii(orbit); 
			i++) {
			draw_neighbour2(vdev, p, Px, Py, o1, i, o2, j);
			}
		}
	else if (p->draw_lines_type == 1) { /* Ainf */
		i = 0;
		for (j = 0; 
			j < p->orbit_size->s_ii(neighbour); 
			j++) {
			draw_neighbour2(vdev, p, Px, Py, o2, j, o1, i);
			}
		}
	else if (p->draw_lines_type == 2) { /* all */
		for (i = 0; 
			i < p->orbit_size->s_ii(orbit); 
			i++) {
			for (j = 0; 
				j < p->orbit_size->s_ii(neighbour); 
				j++) {
				draw_neighbour2(vdev, p, Px, Py, o1, i, o2, j);
				}
			}
		}
	else if (p->draw_lines_type == 3) { // only poset of orbits
		x1 = x[100 + orbit];
		y1 = y[100 + orbit];
		x2 = x[100 + neighbour];
		y2 = y[100 + neighbour];
		Bind_mol(vdev, 1, x1, y1, x2, y2, 1);
		}
	return TRUE;
}

INT draw_neighbour2(VDEVICE *vdev, DRAW_LOCAL *p, INT *Px, INT *Py, 
	INT o1, INT rep1, INT o2, INT rep2)
{
	MATRIX_OP I;
	INT pix_x0, pix_y0, pix_x1, pix_y1;
	
	printf("draw_neighbour2 o1 = %ld rep1 = %ld o2 = %ld rep2 = %ld\n", 
		o1, rep1, o2, rep2);
	if (o1 < o2) { 
		I = (MATRIX_OP) p->incma->s_ij(o1, o2);
		I->Print();
		if (I->s_iji(rep1, rep2)) {
			get_ko(vdev, p, o1, rep1, Px, Py, &pix_x0, &pix_y0);
			get_ko(vdev, p, o2, rep2, Px, Py, &pix_x1, &pix_y1);
			bind_mol(vdev, pix_x0, pix_y0, pix_x1, pix_y1, 1);
			}
		}	
	else {
		I = (MATRIX_OP) p->incma->s_ij(o2, o1);
		I->Print();
		if (I->s_iji(rep2, rep1)) {
			get_ko(vdev, p, o1, rep1, Px, Py, &pix_x0, &pix_y0);
			get_ko(vdev, p, o2, rep2, Px, Py, &pix_x1, &pix_y1);
			bind_mol(vdev, pix_x0, pix_y0, pix_x1, pix_y1, 1);
			}
		}	
	return OK;
}

#define ORBIT_SCALE 0.87

INT get_ko(VDEVICE *vdev, DRAW_LOCAL *p, 
	INT orbit, INT rep, 
	INT *Px, INT *Py, 
	INT *pix_x1, INT *pix_y1)
{
	double orbit_dx, dx1, dx1_halbe;
	double x1, y1, z1, x1o, x1oo;

	x1 = (double) Px[100 + orbit];
	y1 = (double) Py[100 + orbit];
	z1 = 0.;
	orbit_dx = (double) p->Os->s_ii(orbit) * ORBIT_SCALE;
	dx1 = orbit_dx / (double) p->orbit_size->s_ii(orbit);
	dx1_halbe = dx1 * .5;
	x1o = x1 - orbit_dx * .5;
	x1o += rep * dx1;
	x1oo = x1o + dx1_halbe;
	if (!user2vdev(vdev, x1oo, y1, z1, 
		pix_x1, pix_y1, p->extrema)) {
		Srff("get_ko", "user2vdev");
		return FALSE;
		}
	return TRUE;
}

INT lincomb3(INT *Px, INT *Py, INT a1, INT i1, INT a2, INT i2, INT a3, INT i3, INT idx)
{
	INT x, y;
	
	x = a1 * Px[i1] + a2 * Px[i2] + a3 * Px[i3];
	y = a1 * Py[i1] + a2 * Py[i2] + a3 * Py[i3];
	Px[idx] = x;
	Py[idx] = y;
	return OK;
}

INT find_midpoint(INT *Px, INT *Py, INT i1, INT i2, INT idx)
{
	INT x, y;
	
	x = (Px[i1] + Px[i2])/2;
	y = (Py[i1] + Py[i2])/2;
	Px[idx] = x;
	Py[idx] = y;
	
	return OK;
}

INT insert_lattice(DRAW_LOCAL *latt)
{
	INT i, j;
	INT nb_con, nb_gr;
	INT a, b, l1, l2, ii, jj;
	FILE *fp_in;
	INT size_x = 400;
	INT size_y = 564;
	MATRIX_OB Asup, I;
	
	printf("read input file:\n"); fflush(stdout);
	fp_in = fopen(latt->finput, "r");
	printf("is open...\n"); fflush(stdout);
	fscanf(fp_in, "%ld %ld\n", &nb_con, &nb_gr);
	printf("nb of conjugates: %ld, groups total: %ld\n", nb_con, nb_gr);
	fflush(stdout);
	
	latt->nb_con = nb_con;
	latt->nb_gr = nb_gr;
	
	latt->info->m_ilih_n(5, nb_con);
	latt->Ainf->m_ilih_n(nb_con, nb_con);
	latt->incma->m_ilih(nb_con, nb_con);
	latt->basics->m_il(5);
	
	printf("save info as matrix\n"); fflush(stdout);
	for (i = 0; i < nb_con; i++)
	{
		for (j = 0; j < 5; j++)
		{
			fscanf(fp_in, "%ld", &a);
			latt->info->m_iji(i, j, a);
		}
		
	} 
	printf("get Ainf\n"); fflush(stdout);
	for (i = 0; i < nb_con; i++)
	{
		for (j = 0; j < nb_con; j++)
		{
			fscanf(fp_in, "%ld", &a);
			latt->Ainf->m_iji(i, j, a);
		}
		
	}
	printf("get incma\n"); fflush(stdout);
	for (i = 0; i < nb_con; i++)
	{
		l1 = latt->info->s_iji(i, 2);
		for (j = i + 1; j < nb_con; j++)
		{
			l2 = latt->info->s_iji(j, 2);
			fscanf(fp_in, "%ld %ld", &a, &b);
			if (a != i + 1)
				return error("file format error");
			if (b != j + 1)
				return error("file format error");
			I.m_ilih(l2, l1);
			for (ii = 0; ii < l1; ii++) {
				for (jj = 0; jj < l2; jj++) {
					fscanf(fp_in, "%ld", &a);
					I.m_iji(ii, jj, a);
				}
			}
			printf("i=%ld j=%ld\n", i, j);
			I.Print();
			I.swap((MATRIX_OP) latt->incma->s_ij(i, j));	
		}	
	}
	printf("get basics\n"); fflush(stdout);
	for (i = 0; i < 5; i++)
	{
		fscanf(fp_in, "%ld", &a);
		latt->basics->m_ii(i, a);
	}
	fclose(fp_in);
	
	printf("group info:\n");
	latt->info->fprint_raw(stdout);
	printf("Ainf=\n");
	latt->Ainf->fprint_raw(stdout);
	
	latt->Ainf->Ainf2Asup(&Asup);
	printf("Asup=\n");
	Asup.fprint_raw(stdout);
	
	Asup.Asup2Acover(latt->Acover);
	printf("Acover=\n");
	latt->Acover->fprint_raw(stdout);
	
	latt->Acover->Acover2nl(latt->nl);
	printf("neighbour list:\n"); 
	latt->nl->println();
	
	latt->orbit_size->m_il(nb_con);
	for (i = 0; i < nb_con; i++)
	{
		latt->orbit_size->m_ii(i, latt->info->s_iji(i, 2));
	}
	printf("orbit sizes:\n");
	latt->orbit_size->println();
	printf("\n"); fflush(stdout);
	
	place_lattice(latt->nl, latt->orbit_size, size_x, size_y,
		latt->Px, latt->Py, latt->Os, TRUE /* f_upsidedown */, TRUE /* f_v */);
	
	INT x0 = latt->Px->s_ii(0);
	INT y0 = latt->Py->s_ii(0);
	for (i = 0; i < nb_con; i++) {
		latt->Px->m_ii(i, latt->Px->s_ii(i) - x0);
		latt->Py->m_ii(i, latt->Py->s_ii(i) - y0);
	}
	latt->Px->m_ii(nb_con - 1, 0);
	latt->Py->m_ii(nb_con - 1, size_y);
	 
	printf("\nPx=\n");
	latt->Px->println();
	printf("\nPy=\n");
	latt->Py->println();

	return OK;
}
