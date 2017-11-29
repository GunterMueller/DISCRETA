/* ged.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef GRAPHICS_TRUE

#include <DISCRETA/graphics.h>
#include <DISCRETA/divs.h>

/* 
 * GED
 */

#define GED_X_PIX 1000
#define GED_Y_PIX 1000


#define BUF_SIZE 16000

INT GED_OB::read_graph_from_file(INT f_verbose, BYTE *fname)
{
	BYTE buf[BUF_SIZE];
	BYTE str1[256];
	BYTE *p;
	FILE *fp;
	INT line, f_break, f_graph = FALSE, rem_level = 0;
	INT max_x = -1, max_y = -1, x, y, i, l, a, b, f_not_added;
	INT nb_points = 0;
	
	fp = fopen(fname, "r");
	
	
	Init();
	line = 0;
	while (TRUE) {
		buf[0] = 0;
		line++;
		if (fgets(buf, BUF_SIZE, fp) == NULL) {
			if (rem_level) {
				printf("warning: end of file "
				"occured inside REM\n");
				}
			break;
			}
		buf[strlen(buf) - 1] = 0; /* clear new-line */
		if (f_verbose) {
			printf("%ld: %s\n", line, buf);
			fflush(stdout);
			}
		p = buf;
		f_break = FALSE;
		
		if (!s_scan_token(&p, str1))
			continue;

		if (rem_level) {
			if (strcmp(str1, "MER") == 0) {
				rem_level--;
				}
			else if (strcmp(str1, "REM") == 0) {
				rem_level++;
				}
			continue;
			}
		
		if (strcmp(str1, "END") == 0) {
			if (rem_level) {
				printf("warning: END "
				"occured inside REM\n");
				}
			break;
			}
		else if (strcmp(str1, "REM") == 0) {
			rem_level++;
			continue;
			}
		else if (strcmp(str1, "GRAPH") == 0) {
			if (f_graph) {
				printf("error: GRAPH seen twice\n");
				f_break = TRUE;
				}
			f_graph = TRUE;
			if (!s_scan_int(&p, &nb_points)) {
				printf("GRAPH: number of points not specified\n");
				f_break = TRUE;
				}
			if (nb_points <= 0) {
				printf("GRAPH: number of points <= 0\n");
				f_break = TRUE;
				}
			s_KO()->m_ilih(KO_LENGTH, nb_points);
			s_cur_id()->m_i(nb_points);
			for (i = 0; i < nb_points; i++) {
				init_KO(i, i);
				}
			}
		else if (strcmp(str1, "COORDS_RANGE") == 0) {
			if (!f_graph) {
				printf("error: GRAPH keyword not yet seen\n");
				f_break = TRUE;
				}
			if (!s_scan_int(&p, &max_x)) {
				printf("GRAPH: error reading max_x\n");
				f_break = TRUE;
				}
			if (!s_scan_int(&p, &max_y)) {
				printf("GRAPH: error reading max_y\n");
				f_break = TRUE;
				}
			s_max_x()->m_i(max_x);
			s_max_y()->m_i(max_y);
			}
		else if (strcmp(str1, "COORDS") == 0) {
			if (max_x < 0) {
				printf("error: COORDS_RANGE keyword not yet seen\n");
				f_break = TRUE;
				}
			for (i = 0; i < nb_points; i++) {
				fscanf(fp, "%ld %ld", &x, &y);
			
#if 0
				if (!s_scan_int(&p, &x)) {
					printf("GRAPH: error reading x\n");
					f_break = TRUE;
					}
				if (!s_scan_int(&p, &y)) {
					printf("GRAPH: error reading y\n");
					f_break = TRUE;
					}
#endif

				init_KO_xy(i, x, y);
				}
			}
		else if (strcmp(str1, "TEXT") == 0) {
			if (max_x < 0) {
				printf("error: COORDS_RANGE keyword not yet seen\n");
				f_break = TRUE;
				}
			
			if (!s_scan_int(&p, &x)) {
				printf("GRAPH: error reading x\n");
				f_break = TRUE;
				}
			if (!s_scan_int(&p, &y)) {
				printf("GRAPH: error reading y\n");
				f_break = TRUE;
				}
			i = s_KO()->s_hi();
			add_KO(x, y, FALSE /* f_visible */);
			init_KO_string(i, p);
			}
		else if (strcmp(str1, "ADJ") == 0) {
			if (!f_graph) {
				printf("error: GRAPH keyword not yet seen\n");
				f_break = TRUE;
				}
			if (!s_scan_int(&p, &a)) {
				printf("GRAPH: error reading a\n");
				f_break = TRUE;
				}
			if (!s_scan_int(&p, &l)) {
				printf("GRAPH: error reading l\n");
				f_break = TRUE;
				}
			for (i = 0; i < l; i++) {
				if (!s_scan_int(&p, &b)) {
					printf("GRAPH: error reading b\n");
					f_break = TRUE;
					}
				add_line(a, b, &f_not_added, FALSE);
				}
			}
		else {
			printf("unknown command: %s\n", str1);
			}
			
		if (f_break)
			break;
		}
			
	
	fclose(fp);
	return OK;
}

INT GED_OB::field_name(INT i, INT j, BYTE *str)
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "max_x"; break;
	case 1: s = "max_y"; break;
	case 2: s = "cur_id"; break;
	case 3: s = "sel_idx"; break;
	case 4: s = "old_sel_idx"; break;
	case 5: s = "KO"; break;
	case 6: s = "lines"; break;
	default:
		return error("GED::field_name()|"
		"i out of range");
	}
	strcpy(str, s);
	return OK;
}

INT GED_OB::Init()
{
	INT erg = OK;
	
	erg += m_il(7);
	c_obj_k(GED_KIND);
	s_max_x()->m_i(GED_X_PIX);
	s_max_y()->m_i(GED_Y_PIX);
	s_cur_id()->m_i(0);
	s_sel_idx()->m_i(-1);
	s_old_sel_idx()->m_i(-1);
	s_KO()->m_ilih_n(KO_LENGTH, 0);
	s_lines()->m_ilih_n(LINES_LENGTH, 0);
	
	return(OK);
}

INT GED_OB::sprint(BYTE *s)
{
	BYTE str1[256];
	
	sprintf(str1, "GED");
	if (strlen(s) + strlen(str1) < 200)
		strcat(s, str1);
	return OK;
}

INT GED_OB::click(GED_LOCAL *gl, 
	VDEVICE *vdev, INT x, INT y, 
	INT click_count, 
	INT f_shift, INT f_control, 
	INT button, INT f_points, INT f_lines)
{
	double dx, dy;
	INT ix, iy;
	INT new_KO_idx, new_line_idx;
	INT idx_nearest, f_found;
	INT f_not_added;
	
	vdev2user(vdev, x, y, 
		&dx, &dy, gl->extrema);
	ix = (INT) dx;
	iy = (INT) dy;
	search_nearest(gl, ix, iy, 
		&idx_nearest, &f_found);
	if (f_found)
		printf("idx_nearest = %ld "
			"id = %ld\n", idx_nearest, 
			s_KO_iji(idx_nearest, 0));

	if (click_count == 1) {
		if (f_found && f_control) {
			/* 1 x Klick auf Knoten mit control:
			 * Loeschen des Knotens: */

			draw_KO_full(gl, vdev, idx_nearest);
				/* einmal Zeichnen zum Loeschen */
			delete_KO_idx_full(idx_nearest);
			return OK;
			}
		if (f_found && idx_nearest != s_sel_idx_i()) {
			/* neue Selektion: */

			new_selection(gl, vdev, idx_nearest);
			/* uebertrage sel_idx nach old_sel_idx 
			 * um evtl. alte Selektion z.B. beim 
			 * generieren von Linien 
			 * zur Verfuegung zu haben. */
			}
		else if (f_found) {
			/* Deselektion: */

			new_selection(gl, vdev, -1);
			}
		/* wenn nicht gefunden, so ignoriere;
		 * es kann evtl. ein Doppelklick werden. 
		 * Eine eventuelle Selection bleibt in sel_idx 
		 * (wird nicht nach old_sel_idx uebertragen !) */
		}
	else if (click_count == 2) {
		if (!f_found) {
			/* neuen Knoten generieren: */
			add_KO(ix, iy, f_points /* f_visible */);
			new_KO_idx = s_KO()->s_hi() - 1;
			draw_ko(gl, vdev, new_KO_idx);

			if (f_lines && s_sel_idx_i() != -1) {
				/* Erinnerung: falls nicht auf 
				 * einen Knoten geklickt wurde, 
				 * so bleibt sel_idx erhalten. 
				 * Hier wird jetzt eine Linie 
				 * zum zuletzt selektierten Knoten 
				 * generiert: */
				add_line_idx(s_sel_idx_i(), 
					new_KO_idx, &f_not_added, FALSE);
				if (f_not_added) {
					/* Eine solche Linie 
					 * gibt es bereits */
					return OK;
					}
				new_line_idx = s_lines()->s_hi() - 1;
				draw_line(gl, vdev, new_line_idx);
				new_selection(gl, vdev, new_KO_idx);
				}
			}
		else {
			if (f_shift) {
				/* Doppelklick mit shift:
				 * Editieren des Knotentextes: */
				draw_ko(gl, vdev, idx_nearest);
				KO_edit_text(idx_nearest);
				draw_ko(gl, vdev, idx_nearest);
				return OK;
				}
			if (f_lines && s_old_sel_idx_i() != -1) {
				/* Linie generieren, 
				 * Endpunkt existiert bereits. */
				printf("s_old_sel_idx_i() = %ld "
					"idx_nearest = %ld\n", 
					s_old_sel_idx_i(), idx_nearest);
				add_line_idx(s_old_sel_idx_i(), 
					idx_nearest, &f_not_added, FALSE);
				if (f_not_added)
					return OK;
				new_line_idx = s_lines()->s_hi() - 1;
				draw_line(gl, vdev, new_line_idx);
				new_selection(gl, vdev, idx_nearest);
				}
			}
		}
	return OK;
}

INT GED_OB::track(GED_LOCAL *gl, VDEVICE *vdev, 
	INT x0, INT y0, INT x1, INT y1, 
	INT f_shift, INT f_control, 
	INT button, INT f_points, INT f_lines, 
	INT *f_have_tracked)
{
	double dx, dy;
	INT ix0, iy0, ix1, iy1;
	INT idx_nearest0, idx_nearest1, f_found;
	
	*f_have_tracked = FALSE;
	vdev2user(vdev, x0, y0, 
		&dx, &dy, gl->extrema);
	ix0 = (INT) dx;
	iy0 = (INT) dy;
	search_nearest(gl, ix0, iy0, 
		&idx_nearest0, &f_found);
	if (!f_found)
		return OK;
	vdev2user(vdev, x1, y1, 
		&dx, &dy, gl->extrema);
	ix1 = (INT) dx;
	iy1 = (INT) dy;
	search_nearest(gl, ix1, iy1, 
		&idx_nearest1, &f_found);
	if (f_found)
		return OK;
	if (draw_KO_full(gl, vdev, idx_nearest0) != OK)
		return error("track: error in draw_KO_full");
	s_KO_ij(idx_nearest0, 1)->m_i(ix1);
	s_KO_ij(idx_nearest0, 2)->m_i(iy1);
	draw_KO_full(gl, vdev, idx_nearest0);
	*f_have_tracked = TRUE;
	return OK;
}

INT GED_OB::track_smooth(
	GED_LOCAL *gl, VDEVICE *vdev, 
	INT x0, INT y0, INT x1, INT y1, 
	INT f_shift, INT f_control, 
	INT button, INT f_points, INT f_lines, 
	INT *f_have_tracked)
{
	double dx, dy;
	INT ix0, iy0, ix1, iy1;
	INT idx_nearest0, idx_nearest1, f_found;
	INT f_not_added, new_line_idx;
	
	*f_have_tracked = FALSE;
	vdev2user(vdev, x0, y0, &dx, &dy, gl->extrema);
	ix0 = (INT) dx;
	iy0 = (INT) dy;
	search_nearest(gl, ix0, iy0, 
		&idx_nearest0, &f_found);
	if (!f_found) {
		add_KO(ix0, iy0, f_points /* f_visible */);
		idx_nearest0 = s_KO()->s_hi() - 1;
		draw_ko(gl, vdev, idx_nearest0);
		}
	vdev2user(vdev, x1, y1, 
		&dx, &dy, gl->extrema);
	ix1 = (INT) dx;
	iy1 = (INT) dy;
	search_nearest(gl, ix1, iy1, 
		&idx_nearest1, &f_found);
	if (!f_found) {
		add_KO(ix1, iy1, f_points /* f_visible */);
		idx_nearest1 = s_KO()->s_hi() - 1;
		draw_ko(gl, vdev, idx_nearest1);
		}

	if (f_lines) {
		add_line_idx(idx_nearest0, 
			idx_nearest1, &f_not_added, FALSE);
		if (f_not_added) {
			/* *f_have_tracked = TRUE; */
			return OK;
			}
		new_line_idx = s_lines()->s_hi() - 1;
		draw_line(gl, vdev, new_line_idx);
		}
	*f_have_tracked = TRUE;
	return OK;
}

INT GED_OB::delete_KO_idx_full(INT i)
{
	INT id, j, ID0, ID1;

	id = s_KO_iji(i, 0);
	if (s_old_sel_idx_i() == i)
		s_old_sel_idx()->m_i(-1);
	if (s_sel_idx_i() == i)
		s_sel_idx()->m_i(-1);
	for (j = 0; j < s_lines()->s_hi(); j++) {
		ID0 = s_lines_iji(j, 1);
		ID1 = s_lines_iji(j, 2);
		if (ID0 != id && ID1 != id)
			continue;
		del_line(j);
		j--;
		}
	del_KO(i);
	return OK;
}

static INT st_InputINT(BYTE *text, 
	INT old_value, INT *new_value);
static INT st_InputSTRING(BYTE *text, 
	BYTE old_value[256], BYTE new_value[256]);

INT GED_OB::KO_edit_text(INT i)
{
	STRING_OP s;
	BYTE str1[256], str2[256];
	INT font_nr, font_size;

	s = (STRING_OP) s_KO()->s_ij(i, 4);
	strcpy(str1, s->s_str());
	st_InputSTRING("text", str1, str2);
	s->init(str2);
	st_InputINT("which font", 
		s_KO_iji(i, 5), &font_nr);
	st_InputINT("size in points", 
		s_KO_iji(i, 6), &font_size);
	s_KO_ij(i, 5)->m_i(font_nr);
	s_KO_ij(i, 6)->m_i(font_size);
	return OK;
}

static INT st_InputINT(BYTE *text, 
	INT old_value, INT *new_value)
{
	BYTE str[256], c;
	INT i;
	
	printf("%s (%ld): ", text, old_value);
	fflush(stdout);
	i = 0;
	str[0] = 0;
	while (TRUE) {
		c = (BYTE) getchar();
		if (c == '\n')
			break;
		str[i] = c;
		str[i + 1] = 0;
		i++;
		}
	printf("%s\n", str);
	/* scanf("%s\n", str); */
	if (strlen(str) == 0)
		*new_value = old_value;
	else
		sscanf(str, "%ld", new_value);
	return TRUE;
}

static INT st_InputSTRING(BYTE *text, 
	BYTE old_value[256], BYTE new_value[256])
{
	BYTE str[256], c;
	INT i;
	
	printf("%s (%s): ", text, old_value);
	fflush(stdout);
	i = 0;
	str[0] = 0;
	while (TRUE) {
		c = (BYTE) getchar();
		if (c == '\n')
			break;
		str[i] = c;
		str[i + 1] = 0;
		i++;
		}
	printf("%s\n", str);
	if (strlen(str) == 0)
		strcpy(new_value, old_value);
	else
		strcpy(new_value, str);
	return TRUE;
}

INT GED_OB::new_selection(GED_LOCAL *gl, 
	VDEVICE *vdev, INT new_sel)
{
	if (s_sel_idx_i() != -1)
		draw_ko_selection(gl, vdev, s_sel_idx_i());
	s_old_sel_idx()->m_i(s_sel_idx_i());
	s_sel_idx()->m_i(new_sel);
	if (s_sel_idx_i() != -1)
		draw_ko_selection(gl, vdev, s_sel_idx_i());
	return OK;
}

INT GED_OB::add_line_idx(INT idx1, INT idx2, INT *f_not_added, INT f_v)
{
	INT id1, id2;
	
	if (idx1 >= s_KO()->s_hi() || idx1 < 0) {
		printf("add_line_idx()|idx1 >= s_hi()\n");
		return OK;
		}
	if (idx2 >= s_KO()->s_hi() || idx2 < 0) {
		printf("add_line_idx()|idx2 >= s_hi()\n");
		return OK;
		}
	id1 = s_KO_iji(idx1, 0);
	id2 = s_KO_iji(idx2, 0);
	add_line(id1, id2, f_not_added, f_v);
	return OK;
}

INT GED_OB::add_line(INT id0, INT id1, INT *f_not_added, INT f_v)
{
	INT i, id, idx, f_found;
	
	*f_not_added = FALSE;
	if (id0 == id1) {
		*f_not_added = TRUE;
		return OK;
		}
	search_line_id0_id1(id0, id1, &idx, &f_found);
	if (f_found) {
		*f_not_added = TRUE;
		return OK;
		}
	s_lines()->inc_row();
	i = s_lines()->s_hi() - 1;
	id = s_cur_id_i();
	s_cur_id()->inc();
	s_lines_ij(i, 0)->m_i(id);
	s_lines_ij(i, 1)->m_i(id0);
	s_lines_ij(i, 2)->m_i(id1);
	if (f_v) {
		printf("new line %ld: %ld %ld\n", id, id0, id1);
		}
	return OK;
}

INT GED_OB::del_line(INT i)
{
	INT j;

	for (j = i + 1; j < s_lines()->s_hi(); j++) {
		s_lines_ij(j - 1, 0)->m_i(s_lines_iji(j, 0));
		s_lines_ij(j - 1, 1)->m_i(s_lines_iji(j, 1));
		s_lines_ij(j - 1, 2)->m_i(s_lines_iji(j, 2));
		}
	s_lines()->s_h()->dec();
	return OK;
}

INT GED_OB::search_line_id0_id1(INT id0, INT id1, 
	INT *idx, INT *f_found)
{
	INT i, l, ID0, ID1;
	
	l = s_lines()->s_hi();
	for (i = 0; i < l; i++) {
		ID0 = s_lines_iji(i, 1);
		ID1 = s_lines_iji(i, 2);
		if ((ID0 == id0 && ID1 == id1) || 
			(ID0 == id1 && ID1 == id0)) {
			*idx = i;
			*f_found = TRUE;
			return OK;
			}
		}
	*f_found = FALSE;
	return OK;
}

INT GED_OB::add_KO(INT x, INT y, INT f_visible)
{
	INT idx, id;
	STRING_OB s;
	
	s_KO()->inc_row();
	idx = s_KO()->s_hi() - 1;
	id = s_cur_id_i();
	s_cur_id()->inc();
	s_KO_ij(idx, 0)->m_i(id);
	s_KO_ij(idx, 1)->m_i(x);
	s_KO_ij(idx, 2)->m_i(y);
	s_KO_ij(idx, 3)->m_i(f_visible);
	s.init("");
	s.swap(s_KO()->s_ij(idx, 4));
	s_KO_ij(idx, 5)->m_i(8);
	/* 0: Times-Roman 
	 * 4: Helvetica
	 * 8: Courier
	 * 12: Symbol */
	s_KO_ij(idx, 6)->m_i(12); /* point */
	return OK;
}

INT GED_OB::init_KO(INT i, INT id)
{
	STRING_OB s;
	
	s_KO_ij(i, 0)->m_i(id);
	s_KO_ij(i, 1)->m_i(0 /* x */ );
	s_KO_ij(i, 2)->m_i(0 /* y */ );
	s_KO_ij(i, 3)->m_i(TRUE /* f_visible */ );
	s.init("");
	s.swap(s_KO()->s_ij(i, 4));
	s_KO_ij(i, 5)->m_i(8);
	/* 0: Times-Roman 
	 * 4: Helvetica
	 * 8: Courier
	 * 12: Symbol */
	s_KO_ij(i, 6)->m_i(12); /* point */
	return OK;
}

INT GED_OB::init_KO_xy(INT i, INT x, INT y)
{
	s_KO_ij(i, 1)->m_i(x);
	s_KO_ij(i, 2)->m_i(y);
	return OK;
}

INT GED_OB::init_KO_string(INT i, BYTE *str)
{
	STRING_OB s;
	
	s.init(str);
	s.swap(s_KO()->s_ij(i, 4));
	return OK;
}

INT GED_OB::del_KO(INT i)
{
	INT j;
	
	s_KO()->s_ij(i, 4)->freeself();
	for (j = i + 1; j < s_KO()->s_hi(); j++) {
		s_KO_ij(j - 1, 0)->m_i(s_KO_iji(j, 0));
		s_KO_ij(j - 1, 1)->m_i(s_KO_iji(j, 1));
		s_KO_ij(j - 1, 2)->m_i(s_KO_iji(j, 2));
		s_KO_ij(j - 1, 3)->m_i(s_KO_iji(j, 3));
		s_KO_ij(j - 1, 4)->swap(s_KO_ij(j, 4));
		s_KO_ij(j - 1, 5)->m_i(s_KO_iji(j, 5));
		s_KO_ij(j - 1, 6)->m_i(s_KO_iji(j, 6));
		}
	s_KO()->s_h()->dec();
	return OK;
}

INT GED_OB::search_KO_id(
	INT id, INT *idx, INT *f_found)
{
	INT i, l, id1;
	
	l = s_KO()->s_hi();
	for (i = 0; i < l; i++) {
		id1 = s_KO_iji(i, 0);
		if (id1 == id) {
			*idx = i;
			*f_found = TRUE;
			return OK;
			}
		}
	*f_found = FALSE;
	return OK;
}

INT GED_OB::search_nearest(GED_LOCAL *gl, 
	INT x, INT y, INT *idx, INT *f_found)
{
	INT i, l, x1, y1, d, d_min, i_min;
	
	d_min = s_max_x_i() / 50;
	i_min = -1;
	l = s_KO()->s_hi();
	for (i = 0; i < l; i++) {
		x1 = s_KO_iji(i, 1);
		y1 = s_KO_iji(i, 2);
		d = (INT) sqrt((double) 
			((x1 - x) * (x1 - x) + (y1 - y) * (y1 - y) ));
		if (d < d_min)
			i_min = i;
		}
	if (i_min == -1) {
		*f_found = FALSE;
		return OK;
		}
	*f_found = TRUE;
	*idx = i_min;
	return OK;
}

INT GED_OB::draw(GED_LOCAL *gl, VDEVICE *vdev)
{
	INT i;
	
	if (s_KO()->s_li() == 5) {
		s_KO()->realloc_column(KO_LENGTH);
		for (i = 0; i < s_KO()->s_hi(); i++) {
			s_KO_ij(i, 5)->m_i(8);
			s_KO_ij(i, 6)->m_i(12);
			}
		}
	Vswr_mode(vdev, 1); /* 1 = Copy 2 = XOR */
	for (i = 0; i < s_KO()->s_hi(); i++) {
		draw_ko(gl, vdev, i);
		if (i == s_sel_idx_i())
			draw_ko_selection(gl, vdev, i);
		}
	for (i = 0; i < s_lines()->s_hi(); i++) {
		draw_line(gl, vdev, i);
		}
	return OK;
}

static void ged_draw_ko(VDEVICE *vdev, 
	INT x, INT y, INT f_visible, 
	BYTE *text, INT font_nr, INT font_size);
static void ged_draw_selection(
	VDEVICE *vdev, INT x, INT y);
static void ged_line(VDEVICE *vdev, 
	INT x0, INT y0, INT x1, INT y1);

INT GED_OB::draw_KO_full(
	GED_LOCAL *gl, VDEVICE *vdev, INT i)
{
	INT id, id1, id2, j;

	if (i >= s_KO()->s_hi() || i < 0) {
		printf("draw_KO_full()|i >= s_hi()\n");
		return OK;
		}
	if (i == s_sel_idx_i())
		draw_ko_selection(gl, vdev, i);
	draw_ko(gl, vdev, i);
	id = s_KO_iji(i, 0);
	for (j = 0; j < s_lines()->s_hi(); j++) {
		id1 = s_lines_iji(j, 1);
		id2 = s_lines_iji(j, 2);
		if (id1 == id || id2 == id)
			draw_line(gl, vdev, j);
		}
	return OK;
}

INT GED_OB::draw_ko(
	GED_LOCAL *gl, VDEVICE *vdev, INT i)
{
	double x, y, z;
	INT id, pix_x1, pix_y1;
	INT f_visible, font_nr, font_size;
	STRING_OP s;
	BYTE str[256];
	
	if (i >= s_KO()->s_hi() || i < 0) {
		printf("draw_ko()|i >= s_hi()\n");
		return OK;
		}
	Vswr_mode(vdev, 1); /* 1 = Copy 2 = XOR */
	id = s_KO_iji(i, 0);
	x = (double) s_KO_iji(i, 1);
	y = (double) s_KO_iji(i, 2);
	z = 0.;
	if (!user2vdev(vdev, x, y, z, 
		&pix_x1, &pix_y1, gl->extrema)) {
		Srff("GED::draw_ko", "user2vdev");
		return(FALSE);
		}
	// printf("GED::draw_ko x=%ld y=%ld -> %ld,%ld\n", 
	//	s_KO_iji(i, 1), s_KO_iji(i, 2), pix_x1, pix_y1);
	sprintf(str, "%ld", id);
	f_visible = s_KO_iji(i, 3);
	s = (STRING_OP) s_KO()->s_ij(i, 4);
	font_nr = s_KO_iji(i, 5);
	font_size = s_KO_iji(i, 6);
	ged_draw_ko(vdev, pix_x1, pix_y1, f_visible, 
		s->s_str(), font_nr, font_size);
	return OK;
}

INT GED_OB::draw_ko_selection(
	GED_LOCAL *gl, VDEVICE *vdev, INT i)
{
	double x, y, z;
	INT id, pix_x1, pix_y1;
	
	if (i >= s_KO()->s_hi() || i < 0) {
		printf("draw_ko_selection()|i >= s_hi()\n");
		return OK;
		}
	Vswr_mode(vdev, 2); /* 1 = Copy 2 = XOR */
	id = s_KO_iji(i, 0);
	x = (double) s_KO_iji(i, 1);
	y = (double) s_KO_iji(i, 2);
	z = 0.;
	if (!user2vdev(vdev, x, y, z, 
		&pix_x1, &pix_y1, gl->extrema)) {
		Srff("GED::draw_ko", "user2vdev");
		return(FALSE);
		}
	ged_draw_selection(vdev, pix_x1, pix_y1);
	return OK;
}

INT GED_OB::draw_line(
	GED_LOCAL *gl, VDEVICE *vdev, INT i)
{
	double x, y, z;
	INT id, id1, id2, idx1, idx2;
	INT f_found1, f_found2;
	INT pix_x1, pix_y1, pix_x2, pix_y2;
	
	if (i >= s_lines()->s_hi() || i < 0) {
		printf("draw_line()|i >= s_hi()\n");
		return OK;
		}
	Vswr_mode(vdev, 1); /* 1 = Copy 2 = XOR */
	id = s_lines_iji(i, 0);
	id1 = s_lines_iji(i, 1);
	id2 = s_lines_iji(i, 2);
	search_KO_id(id1, &idx1, &f_found1);
	search_KO_id(id2, &idx2, &f_found2);
	if (!f_found1 || !f_found2) {
		Srfs("GED::draw_line", 
			"!f_found1 || !f_found2");
		printf("id1 = %ld idx1 = %ld "
			"f_found1 = %ld\n", id1, idx1, f_found1);
		printf("id2 = %ld idx2 = %ld "
			"f_found2 = %ld\n", id1, idx1, f_found1);
		return(FALSE);
		}
	x = (double) s_KO_iji(idx1, 1);
	y = (double) s_KO_iji(idx1, 2);
	z = 0.;
	if (!user2vdev(vdev, x, y, z, 
		&pix_x1, &pix_y1, gl->extrema)) {
		Srff("GED::draw_line", "user2vdev");
		return(FALSE);
		}
	x = (double) s_KO_iji(idx2, 1);
	y = (double) s_KO_iji(idx2, 2);
	z = 0.;
	if (!user2vdev(vdev, x, y, z, 
		&pix_x2, &pix_y2, gl->extrema)) {
		Srff("GED::draw_line", "user2vdev");
		return(FALSE);
		}
	ged_line(vdev, pix_x1, pix_y1, pix_x2, pix_y2);
#if 0
	sprintf(str, "%ld %ld", id1, id2);
	pix_x3 = (pix_x1 + pix_x2) >> 1;
	pix_y3 = (pix_y1 + pix_y2) >> 1;
	V_gtext(vdev, pix_x3, pix_y3, str);
#endif
	return OK;
}

static void ged_draw_ko(VDEVICE *vdev, 
	INT x, INT y, INT f_visible, 
	BYTE *text, INT font_nr, INT font_size)
{
	INT rad, off_x, off_y;
	INT char_width, char_height;
	INT cell_width, cell_height;
	INT hor_in, hor_out;
	INT vert_in, vert_out;
	
#ifdef SYSTEMUNIX
	rad = 5;
	off_x = 15;
	off_y = 20;
#else
	rad = 4;
	off_x = 5;
	off_y = 8;
#endif
	if (f_visible) {
		if (strlen(text)) {
			Vst_font(vdev, font_nr);
			Vst_point(vdev, font_size, 
				&char_width, &char_height, 
				&cell_width, &cell_height); 
			V_pie_text(vdev, x, y, rad, 0, 3600, text);
			// printf("GED:: V_pie_tex() x = %ld y = %ld\n", x, y);
			return;
			// printf("ged.C visible KO, V_pie_text(%s)\n", text);
			}
		else {
			V_pie(vdev, x, y, rad, 0, 3600);
			}
		}
	else {
		if (vdev->type != 0 && vdev->type != 3 && vdev->type != 4)
			/* 0 = PS, 3 = epic file output 4 = MetaPost */
			V_arc(vdev, x, y, rad, 0, 3600);
		}
	if (strlen(text)) {
		Vst_font(vdev, font_nr);
		Vst_point(vdev, font_size, 
			&char_width, &char_height, 
			&cell_width, &cell_height); 
		hor_in = 1; // centered
		vert_in = 2; // up (hanging text)
		Vst_alignment(vdev, hor_in, vert_in, &hor_out, &vert_out);
		V_gtext(vdev, x + off_x, y + off_y, text);
		}
}

static void ged_draw_selection(
	VDEVICE *vdev, INT x, INT y)
{
	INT rad, off_x, off_y, sox, soy;
	/* INT pxy[32]; */
	
#ifdef SYSTEMUNIX
	rad = 5;
	off_x = 15;
	off_y = 20;
#else
	rad = 4;
	off_x = 5;
	off_y = 8;
#endif
	sox = rad << 1;
	soy = rad << 1;
	if (vdev->type != 0 && vdev->type != 3)
		/* not to .ps, (epic) .tex */
		draw_box(vdev, x - sox, y - soy, 
			x + sox, y + soy);
}

static void ged_line(VDEVICE *vdev, 
	INT x0, INT y0, INT x1, INT y1)
{
	INT pxy[4];
	
	pxy[0] = x0;
	pxy[1] = y0;
	pxy[2] = x1;
	pxy[3] = y1;
	V_pline(vdev, 2, pxy);
}

#endif /* GRAPHICS_TRUE */

