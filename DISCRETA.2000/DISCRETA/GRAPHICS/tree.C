/* tree.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1997
 */


#include <DISCRETA/discreta.h>

#ifdef GRAPHICS_TRUE

#include <DISCRETA/divs.h>
#include <DISCRETA/graphics.h>

#define BUF_SIZE 16000

static INT add_to_tree(TREE_OP p, VECTOR_OP V, INT *len, 
	INT first_V, INT len_V, INT *idx, INT depth, INT f_v);
static void my_place(TREE_OP p, INT x_unit, INT y_unit, INT idx, INT depth);
static void write_out_ged(TREE_OP p, GED_OP G, INT idx, INT depth, INT f_v);
static void my_ged_label(TREE_OP p, GED_OP G, INT idx, INT depth);
static void my_ged_user_label(TREE_OP p, GED_OP G, INT idx, INT depth, INT f_v);
static void print_node(MATRIX_OP p, INT i);
static INT get_value(VECTOR_OP V, INT i, INT depth);
static INT get_de(VECTOR_OP V, INT i);

INT TREE_OB::field_name(INT i, INT j, BYTE *str)
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "depth"; break;
	case 1: s = "nb_leaves"; break;
	case 2: s = "nb_nodes"; break;
	case 3: s = "V"; break;
	case 4: s = "Tree"; break;
	case 5: s = "X"; break;
	case 6: s = "Y"; break;
	case 7: s = "max_x"; break;
	case 8: s = "max_y"; break;
	default:
		return error("TREE::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

INT TREE_OB::Init()
{
	INT erg = OK;
	
	erg += m_il(9);
	c_obj_k(TREE_KIND);
	
	return(OK);
}

INT TREE_OB::sprint(BYTE *s)
{
	BYTE str1[256];
	
	sprintf(str1, "TREE");
	if (strlen(s) + strlen(str1) < 200)
		strcat(s, str1);
	return OK;
}

INT TREE_OB::read_tree_from_file(INT f_verbose, BYTE *fname)
{
	BYTE buf[BUF_SIZE];
	BYTE str1[256];
	BYTE str2[256];
	BYTE str3[256];
	BYTE *p, *pp;
	FILE *fp;
	INT line, f_break, f_tree = FALSE, rem_level = 0;
	INT max_x = -1, max_y = -1;
	INT x, i, l, a;
	INT max_depth, de;
	VECTOR_OB V, v;
	STRING_OB label;
	INTEGER_OB V_len;
	
	fp = fopen(fname, "r");
	
	
	Init();
	line = 0;
	while (TRUE) {
		buf[0] = 0;
		line++;
		if (fgets(buf, BUF_SIZE, fp) == NULL) {
			if (rem_level) {
				printf("warning: end of file occured inside REM\n");
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

		if (strncmp(str1, "TREE", 4) == 0) {
			if (f_tree) {
				printf("error: TREE seen twice\n");
				f_break = TRUE;
				break;
				}
			if (max_x == -1) {
				printf("error: COORDS_RANGE not seen !\n");
				f_break = TRUE;
				break;
				}
			f_tree = TRUE;
#if 0
			if (!s_scan_int(&p, &depth)) {
				printf("TREE: depth not specified\n");
				f_break = TRUE;
				}
			if (depth <= 0) {
				printf("TREE: depth <= 0\n");
				f_break = TRUE;
				}
			s_depth()->m_i(depth);
#endif
			max_depth = 0;
			V.m_il(VECTOR_OVERSIZE);
			V_len.m_i(0);
			while (TRUE) {
				buf[0] = 0;
				line++;
				if (fgets(buf, BUF_SIZE, fp) == NULL) {
					printf("warning: end of file occured inside TREE\n");
					}
				buf[strlen(buf) - 1] = 0; /* clear new-line */
				if (f_verbose) {
					printf("%ld: %s\n", line, buf);
					fflush(stdout);
					}
				p = buf;
				if (strncmp(buf, "TREE_END", 8) == 0)
					break;

				if (rem_level) {
					if (strncmp(buf, "MER", 3) == 0) {
						rem_level--;
						}
					else if (strncmp(buf, "REM", 3) == 0) {
						rem_level++;
						}
					continue;
					}
				
				if (strncmp(buf, "END", 3) == 0) {
					if (rem_level) {
						printf("warning: END occured inside REM\n");
						}
					break;
					}
				else if (strncmp(buf, "REM", 3) == 0) {
					rem_level++;
					continue;
					}
				
				if (!s_scan_int(&p, &de)) {
					printf("TREE: error reading depth\n");
					f_break = TRUE;
					}
				max_depth = MAXIMUM(de, max_depth);
				v.m_il(de + 1);
				for (i = 0; i < de; i++) {
					if (!s_scan_int(&p, &x)) {
						printf("TREE: error reading x in line %ld\n", line);
						f_break = TRUE;
						}
					v.m_ii(i, x);
					}
				pp = p;
				while (*pp != 0) {
					if (pp[0] == '%' && pp[1] == '%' && pp[2] == '%')
						*pp = 0;
					else
						pp++;
					}
				label.init(p);
				label.swap((STRING_OP) v.s_i(de));
				V.append_element(&V_len, &v);
				} // while
			s_depth()->m_i(max_depth);
			V.realloc_z(V_len.s_i());
			if (f_verbose) {
				l = V.s_li();
				printf("read TREE of length %ld max_depth = %ld\n", l, max_depth);
				V.Print();
				fflush(stdout);
				}
			build_tree(&V, f_verbose);
			place(f_verbose);
			V.swap(s_V());
			}
		else if (strncmp(str1, "COORDS_RANGE", 12) == 0) {
			if (!s_scan_int(&p, &max_x)) {
				printf("TREE: error reading max_x\n");
				f_break = TRUE;
				}
			if (!s_scan_int(&p, &max_y)) {
				printf("TREE: error reading max_y\n");
				f_break = TRUE;
				}
			s_max_x()->m_i(max_x);
			s_max_y()->m_i(max_y);
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

INT TREE_OB::build_tree(VECTOR_OP V, INT f_v)
{
	INT h, l = V->s_li();
	INT first, len, first_V, len_V, i, idx;

	h = 1 + (s_depth_i() + 1) * l;
	s_Tree()->m_ilih(5, h);
	if (f_v) {
		printf("allocated Tree of length %ld\n", h);
		}

	// add the root:
	first = 1;
	len = -1;
	first_V = 0;
	len_V = l;
	s_Tree()->m_iji(0, 0, first);
	s_Tree()->m_iji(0, 1, len); // unknown
	s_Tree()->m_iji(0, 2, first_V);
	s_Tree()->m_iji(0, 3, len_V);
	s_Tree()->m_iji(0, 4, -1); // root has no label
	if (f_v) {
		print_node(s_Tree(), 0);
		}
	len = 0;
	idx = first;
	add_to_tree(this, V, &len, first_V, len_V, &idx, 0 /* depth */, f_v);
	s_Tree()->m_iji(0, 1, len);
	
	if (f_v) {
		printf("the tree:\n");
		printf("idx:first_son,num_suns,first_V,len_V,label\n");
		for (i = 0; i < idx; i++) {
			print_node(s_Tree(), i);
			}
		}
	if (idx > h) {
		printf("warning: idx = %ld > h = %ld\n", idx, h);
		}
	s_nb_nodes()->m_i(idx);
	s_nb_leaves()->m_i(l);

	return OK;
}

static INT add_to_tree(TREE_OP p, VECTOR_OP V, INT *len, 
	INT first_V, INT len_V, INT *idx, INT depth, INT f_v)
{
	INT i = 0, i0 = 0, val, val_old;
	INT f, l, fV, lV, nb_different_values;
	INT old_idx = *idx, local_idx, k;

	if (len_V <= 0)
		return error("add_to_tree() len_V <= 0");
	i = 0;
	i0 = 0;
	val_old = get_value(V, first_V + i, depth);
	nb_different_values = 1;
	while (i < len_V) {
		val = get_value(V, first_V + i, depth);
		if (val != val_old) {
			if (f_v && 0) {
				printf("from %ld len %ld in V value %ld at depth %ld\n", 
					first_V + i0, i - i0, val_old, depth);
				}
			nb_different_values++;
			i0 = i;
			}
		val_old = val;
		i++;
		}
	if (f_v && 0) {
		printf("from %ld len %ld in V value %ld at depth %ld\n", 
			first_V + i0, i - i0, val_old, depth);
		}
	if (f_v) {
		printf("add_to_tree() depth %ld, nb_different_values=%ld\n", 
			depth, nb_different_values);
		fflush(stdout);
		}
	(*idx) += nb_different_values;
	
	i = 0;
	i0 = 0;
	k = 0;
	val_old = get_value(V, first_V + i, depth);
	while (i < len_V) {
		val = get_value(V, first_V + i, depth);
		if (val != val_old) {
			local_idx = old_idx + k;
			f = *idx;
			l = 0;
			fV = first_V + i0;
			lV = i - i0;
			p->s_Tree()->m_iji(local_idx, 0, *idx);
			p->s_Tree()->m_iji(local_idx, 1, l); // unknown
			p->s_Tree()->m_iji(local_idx, 2, fV);
			p->s_Tree()->m_iji(local_idx, 3, lV);
			p->s_Tree()->m_iji(local_idx, 4, val_old); // root has no label
			if (val_old != -1) {
				if (f_v) {
					print_node(p->s_Tree(), local_idx);
					}
				add_to_tree(p, V, &l, fV, lV, idx, depth + 1, f_v);
				p->s_Tree()->m_iji(local_idx, 1, l);
				}
			else {
				p->s_Tree()->m_iji(local_idx, 0, -1);
				}
			k++;
			i0 = i;
			}
		val_old = val;
		i++;
		}
	
	local_idx = old_idx + k;
	f = *idx;
	l = 0;
	fV = first_V + i0;
	lV = i - i0;
	p->s_Tree()->m_iji(local_idx, 0, *idx);
	p->s_Tree()->m_iji(local_idx, 1, l); // unknown
	p->s_Tree()->m_iji(local_idx, 2, fV);
	p->s_Tree()->m_iji(local_idx, 3, lV);
	p->s_Tree()->m_iji(local_idx, 4, val_old); // root has no label
	if (val_old != -1) {
		if (f_v) {
			print_node(p->s_Tree(), local_idx);
			}
		add_to_tree(p, V, &l, fV, lV, idx, depth + 1, f_v);
		p->s_Tree()->m_iji(local_idx, 1, l);
		}
	else {
		p->s_Tree()->m_iji(local_idx, 0, -1);
		}

	*len = nb_different_values;
	return OK;
}

INT TREE_OB::place(INT f_v)
{
	INT nb_nodes, nb_leaves;
	INT depth, max_x, max_y;
	INT x, y, y_unit, x_unit;
	INT i;
	double d;
	
	nb_nodes = s_nb_nodes_i();
	nb_leaves = s_nb_leaves_i();
	depth = s_depth_i() + 2; // one for the root, one for the lowest leaves !
	max_x = s_max_x_i();
	max_y = s_max_y_i();
	s_X()->m_il(nb_nodes + 1);
	s_Y()->m_il(nb_nodes + 1);
	y_unit = (INT)((double)max_y / ((double)depth + 2.5));
	x_unit = (INT)((double)max_x / ((double)nb_leaves + 1.5));
	x = s_max_x_i() / 2;
	// d = (double) s_depth_i() + .3; // position of leaves
	d = (double) s_depth_i() + .6;
	y = (INT)(((double)d + .5) * (double)y_unit);
	s_X_i(nb_nodes)->m_i(x);
	s_Y_i(nb_nodes)->m_i(y);
	if (f_v) {
		printf("nb_nodes = %ld\n", s_nb_nodes_i());
		printf("nb_leaves = %ld\n", s_nb_leaves_i());
		printf("depth = %ld\n", depth);
		printf("max_x = %ld\n", max_x);
		printf("max_y = %ld\n", max_y);
		printf("y_unit = %ld\n", y_unit);
		printf("x_unit = %ld\n", x_unit);
		fflush(stdout);
		}

	my_place(this, x_unit, y_unit, 0, 0);
	for (i = 0; i < nb_nodes; i++) {
		// x = s_X_ii(i);
		y = s_Y_ii(i);
		y = max_y - y;
		s_Y_i(i)->m_i(y);
		}
	if (f_v) {
		printf("placement:\n");
		for (i = 0; i < nb_nodes; i++) {
			x = s_X_ii(i);
			y = s_Y_ii(i);
			printf("%ld: %ld,%ld\n", i, x, y);
			}
		printf("the invisible node:\n");
		i = nb_nodes;
		x = s_X_ii(i);
		y = s_Y_ii(i);
		printf("%ld: %ld,%ld\n", i, x, y);
		}

	return OK;
}

static void my_place(TREE_OP p, INT x_unit, INT y_unit, INT idx, INT depth)
{
	INT f, l, fV, lV;
	INT x, y, i;
	double d;

	f = p->s_Tree()->s_iji(idx, 0);
	l = p->s_Tree()->s_iji(idx, 1);
	fV = p->s_Tree()->s_iji(idx, 2);
	lV = p->s_Tree()->s_iji(idx, 3);
	d = (double) depth;
	if (f == -1)
		d = (double) p->s_depth_i() + .3;
	x = (INT)(((double)fV + ((double) lV * .5)  + 0.5) * (double) x_unit);
	y = (INT)(((double)d + .5) * (double)y_unit);
	p->s_X_i(idx)->m_i(x);
	p->s_Y_i(idx)->m_i(y);
	for (i = 0; i < l; i++) {
		if (l) {
			my_place(p, x_unit, y_unit, f + i, depth + 1);
			}
		}
}

INT TREE_OB::init_ged(GED_OP G, INT f_v)
{   
	INT nb_leaves, nb_points, depth;
	INT i, j, k, ll, F, L, f, l, f0, a, a0, next, f_not_added;
	INT x, y, first_i, first_ip1, n, n0;
	VECTOR_OP If, Il, In;
	
	depth = s_depth_i();
	l = s_X()->s_li();
	if (f_v) {
		printf("TREE::init_ged() found %ld nodes\n", l);
		fflush(stdout);
		}
	nb_points = l;
	G->Init();
	G->s_max_x()->m_i(s_max_x_i());
	G->s_max_y()->m_i(s_max_y_i());
	G->s_KO()->m_ilih(KO_LENGTH, nb_points);
	G->s_cur_id()->m_i(nb_points);
	for (i = 0; i < nb_points; i++) {
		G->init_KO(i, i);
		}
	
	// make last node invisible:
	i = nb_points - 1;
	G->s_KO_ij(i, 3)->m_i(FALSE); // invisible
	G->init_KO_string(i, " ");
	
	for (i = 0; i < l; i++) {
		x = s_X_ii(i);
		y = s_Y_ii(i);
		G->init_KO_xy(i, x, y);
		if (f_v) {
			printf("KO %ld = (%ld,%ld)\n", i, x, y);
			fflush(stdout);
			}
		}
	write_out_ged(this, G, 0, 0, f_v);
	return OK;
}   
    
static void write_out_ged(TREE_OP p, GED_OP G, INT idx, INT depth, INT f_v)
{
	INT f, l, fV, lV;
	INT x0, y0, x, y, i;
	INT f_not_added;

	f = p->s_Tree()->s_iji(idx, 0);
	l = p->s_Tree()->s_iji(idx, 1);
	fV = p->s_Tree()->s_iji(idx, 2);
	lV = p->s_Tree()->s_iji(idx, 3);
	for (i = 0; i < l; i++) {
		if (p->s_Tree()->s_iji(f + i, 0) != -1)
			G->add_line(idx, f + i, &f_not_added, f_v);
		if (l) {
			write_out_ged(p, G, f + i, depth + 1, f_v);
			}
		}
}

INT TREE_OB::init_ged_labels(GED_OP G, INT f_v)
{
	INT depth, l;
	
	depth = s_depth_i();
	l = s_X()->s_li();
	if (f_v) {
		printf("TREE::init_ged_labels() found %ld nodes\n", l);
		fflush(stdout);
		}

	my_ged_label(this, G, 0, 0);
	return OK;
}

static void my_ged_label(TREE_OP p, GED_OP G, INT idx, INT depth)
{
	INT f, l, fV, lV;
	INT x0, y0, x, y, i;
	BYTE str[1024];
	INT a;

	f = p->s_Tree()->s_iji(idx, 0);
	l = p->s_Tree()->s_iji(idx, 1);
	fV = p->s_Tree()->s_iji(idx, 2);
	lV = p->s_Tree()->s_iji(idx, 3);
	for (i = 0; i < l; i++) {
		if (p->s_Tree()->s_iji(f + i, 0) != -1) {
			a = p->s_Tree()->s_iji(f + i, 4);
			sprintf(str, "%ld", a);
			G->init_KO_string(f + i, str);
			}
		if (l) {
			my_ged_label(p, G, f + i, depth + 1);
			}
		}
}

INT TREE_OB::init_ged_user_labels(GED_OP G, INT f_v)
{
	INT depth, l;
	
	depth = s_depth_i();
	l = s_X()->s_li();
	if (f_v) {
		printf("TREE::init_ged_user_labels() found %ld nodes\n", l);
		fflush(stdout);
		}

	my_ged_user_label(this, G, 0, 0, f_v);
	return OK;
}

static void my_ged_user_label(TREE_OP p, GED_OP G, INT idx, INT depth, INT f_v)
{
	VECTOR_OP v;
	INT f, l, fV, lV;
	INT x0, y0, x, y, i;
	BYTE str[1024], *pp;
	INT a, de;
	STRING_OP s;

	f = p->s_Tree()->s_iji(idx, 0);
	l = p->s_Tree()->s_iji(idx, 1);
	fV = p->s_Tree()->s_iji(idx, 2);
	lV = p->s_Tree()->s_iji(idx, 3);
	for (i = 0; i < l; i++) {
		if (p->s_Tree()->s_iji(f + i, 0) != -1) {
			}
		else {
			v = p->s_V_i(fV);
			de = get_de(p->s_V(), fV);
			s = (STRING_OP) v->s_i(de);
			if (s->s_obj_k() == STRING)
				pp = s->s_str();
			else
				pp = "";
				// leading $ means no translation in epic mode !
				// disabled now in epic.C !

			sprintf(str, "%s", pp);
			
			if (f_v) {
				printf("leave %ld gets label: %s\n", fV, str);
				fflush(stdout);
				}
			G->init_KO_string(f + i, str);
			G->s_KO_ij(f + i, 3)->m_i(FALSE); // invisible
			}
		if (l) {
			my_ged_user_label(p, G, f + i, depth + 1, f_v);
			}
		}
}

INT TREE_OB::leaves_in_array()
// embeds each leave-text into 
// $\makebox(0,0)[t]{$\begin{array}[t]{l}..\end{array}$}
// where .. stands for the original leave label.
{
	STRING_OP s;
	VECTOR_OP v;
	INT i, l, de;
	BYTE str[10000];

	l = s_V()->s_li();
	for (i = 0; i < l; i++) {
		v = s_V_i(i);
		de = v->s_li() - 1;
		s = (STRING_OP) v->s_i(de);
		sprintf(str, "\\makebox(0,0)[t]{"
			"$\\begin{array}[t]{l}"	
				"%s\\\\"
				"\\end{array}$}", 
			s->s_str());
		s->init(str);
		}
	return OK;
}

static void print_node(MATRIX_OP p, INT i)
{
	printf("%ld: first=%ld len=%ld first_V=%ld len_V=%ld label=%ld\n", 
		i, 
		p->s_iji(i, 0), 
		p->s_iji(i, 1), 
		p->s_iji(i, 2), 
		p->s_iji(i, 3), 
		p->s_iji(i, 4));
	fflush(stdout);

}

static INT get_value(VECTOR_OP V, INT i, INT depth)
{
	VECTOR_OP v;

	v = (VECTOR_OP) V->s_i(i);
	if (v->s_li() - 1 > depth)
		return v->s_ii(depth);
	return -1;
}

static INT get_de(VECTOR_OP V, INT i)
{
	VECTOR_OP v;

	v = (VECTOR_OP) V->s_i(i);
	return v->s_li() - 1;
}


#endif /* GRAPHICS_TRUE */
