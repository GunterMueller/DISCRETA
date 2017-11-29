/* gt_col_util.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SOLVABLE_TRUE

#include <DISCRETA/solvable.h>

#if TEXDOCU
INT color_join(VECTOR_OP col1, VECTOR_OP col2, 
	VECTOR_OP col3, MATRIX_OP Col, INT f_v, INT f_vv)
#endif
{
	VECTOR_OB V;
	INTEGER_OB l1_ob;
	INT i, ii, j, l, l1 = 0, a, b, idx, f_found;
	
	if (f_vv) {
		printf("color_join: old coloring (length %ld %ld):\n", 
			col1->s_li(), col2->s_li());
		col1->println();
		col2->println();
		fflush(stdout);
		}
	l = col1->s_li();
	if (col2->s_li() != l)
		return error("color_join() different length !");
	Col->m_ilih(2, l1);
	V.m_il(2);
	
	/* at first, we determine all colors i.e.
	 * combinations (a,b), a \in col1, b \in col2
	 */
	for (i = 0; i < l; i++) {
		a = col1->s_ii(i);
		b = col2->s_ii(i);
		V.m_ii(0, a);
		V.m_ii(1, b);
		Col->search_row(l1, 2 /* nb_keys */, 
			TRUE /* f_ascending */, &V, &idx, &f_found);
		if (f_found) {
			idx--;
			}
		else {
			if (f_vv) {
				printf("found new color %ld/%ld\n", a, b);
				fflush(stdout);
				}
			l1_ob.m_i(l1);
			Col->insert_row_at(&l1_ob, idx, &V, 2 /* V_len */);
			l1++;
			Col->realloc_row(l1);
			if (f_vv) {
				Col->Print();
				fflush(stdout);
				}
			}
		}
	/* now we create the new vector of colors (col3) */
	col3->m_il(l);
	for (i = 0; i < l; i++) {
		a = col1->s_ii(i);
		b = col2->s_ii(i);
		V.m_ii(0, a);
		V.m_ii(1, b);
		Col->search_row(l1, 2 /* nb_keys */, 
			TRUE /* f_ascending */, &V, &idx, &f_found);
		if (!f_found)
			return error("color_join(): not found");
		idx--;
		col3->m_ii(i, idx);
		}
	if (f_vv) {
		printf("new vector of colors: (length %ld)\n", col3->s_li());
		col3->println();
		fflush(stdout);
		}
	if (f_v) {
		printf("color_join(): found %ld colors\n", Col->s_hi());
		fflush(stdout);
		}
	return OK;
}

#if TEXDOCU
INT color_v_join(VECTOR_OP col_0, VECTOR_OP col_v, 
	VECTOR_OP new_col, MATRIX_OP new_Col, 
	INT f_v, INT f_vv)
#endif
{
	VECTOR_OP col_i;
	VECTOR_OB V;
	INTEGER_OB l1_ob;
	INT i, ii, j, l, l1 = 0, a, b, idx, f_found;
	INT nt, Nt;
	
	if (f_vv) {
		printf("color_v_join()\n");
		printf("old coloring col_0 (length %ld):\n", col_0->s_li());
		col_0->println();
		}
	l = col_0->s_li();
	nt = col_v->s_li();
	if (f_vv) {
		printf("new colorings:\n");
		}
	for (i = 0; i < nt; i++) {
		col_i = (VECTOR_OP) col_v->s_i(i);
		if (col_i->s_li() != l) {
			printf("i = %ld col_i->s_li() = %ld l = %ld\n", i, col_i->s_li(), l);
			fflush(stdout);
			return error("color_v_join() different length !");
			}
		if (f_vv) {
			printf("col_%ld:\n", i);
			col_i->println();
			fflush(stdout);
			}
		}
	Nt = 1 + nt;
	new_Col->m_ilih(Nt, l1);
	V.m_il(Nt);
	
	/* at first, we determine all colors i.e.
	 * combinations (a,b_1,b_2,\ldots,b_nt), a \in col_0, b_i \in col_i
	 */
	for (i = 0; i < l; i++) {
		a = col_0->s_ii(i);
		V.m_ii(0, a);
		for (j = 0; j < nt; j++) {
			col_i = (VECTOR_OP) col_v->s_i(j);
			b = col_i->s_ii(i);
			V.m_ii(1 + j, b);
			}
		new_Col->search_row(l1, Nt /* nb_keys */, 
			TRUE /* f_ascending */, &V, &idx, &f_found);
		if (f_found) {
			idx--;
			}
		else {
			if (f_vv) {
				printf("found new color:");
				V.println();
				fflush(stdout);
				}
			l1_ob.m_i(l1);
			new_Col->insert_row_at(&l1_ob, idx, &V, Nt /* V_len */);
			l1++;
			new_Col->realloc_row(l1);
			if (f_vv) {
				new_Col->Print();
				fflush(stdout);
				}
			}
		}
	/* now we create the new vector of colors (new_col) */
	new_col->m_il(l);
	for (i = 0; i < l; i++) {
		a = col_0->s_ii(i);
		V.m_ii(0, a);
		for (j = 0; j < nt; j++) {
			col_i = (VECTOR_OP) col_v->s_i(j);
			b = col_i->s_ii(i);
			V.m_ii(1 + j, b);
			}
		new_Col->search_row(l1, Nt /* nb_keys */, 
			TRUE /* f_ascending */, &V, &idx, &f_found);
		if (!f_found)
			return error("color_v_join(): not found");
		idx--;
		new_col->m_ii(i, idx);
		}
	if (f_vv) {
		printf("new vector of colors (length %ld):\n", new_col->s_li());
		new_col->println();
		fflush(stdout);
		}
	if (f_v) {
		printf("color_v_join(): found %ld colors\n", new_Col->s_hi());
		fflush(stdout);
		}
	return OK;
}

#if TEXDOCU
INT number_of_colors(MATRIX_OP Col, INT f_v)
#endif
{
	INT nc;
	
	nc = Col->s_hi();
	if (f_v)
		printf("number of colors: %ld\n", nc);
	return nc;
}

#if TEXDOCU
INT prefer_seldom_colors(VECTOR_OP col, VECTOR_OP new_col, 
	PERMUTATION_OP col_perm_new_old, PERMUTATION_OP col_perm_old_new)
#endif
{
	VECTOR_OB val, mult;
	INT i, j, l, j0, a, a0, b, n;
	
	n = col->s_li();
	col->multiplicities(&val, &mult);
	l = mult.s_li();
	col_perm_new_old->m_il(l);
	
	/* lets determine the i-th important color: */
	for (i = 0; i < l; i++) {
		j0 = -1; /* no minimum yet */
		for (j = l - 1; j >= 0; j--) {
			a = mult.s_ii(j);
			if (a == 0)
				continue; /* this color is already chosen */
			if (j0 == -1) { /* lets make the first choice */
				j0 = j;
				a0 = a;
				continue;
				}
			if (a < a0) { /* we found a better color (better = more seldom) */
				j0 = j;
				a0 = a;
				}
			}
		if (j0 < 0)
			return error("prefer_seldom_colors() no color found");
		col_perm_new_old->m_ii(i, j0 + 1);
		mult.m_ii(j0, 0); /* do not choose it again ! */
		} /* next i */
	col_perm_new_old->invers(col_perm_old_new);
	new_col->m_il(n);
	for (i = 0; i < n; i++) {
		a = col->s_ii(i);
		b = col_perm_old_new->s_ii(a) - 1;
		new_col->m_ii(i, b);
		}
	return OK;
}

#if TEXDOCU
INT get_ci_ordered(VECTOR_OP cl, VECTOR_OP eo, 
	VECTOR_OP val, VECTOR_OP mult, INT f_v, INT f_vv)
#endif
{
	VECTOR_OB type;

	type.join2(cl, eo);
	if (f_v) {
		printf("type=");
		type.println();
		fflush(stdout);
		}
	type.multiplicities(val, mult);
	if (f_v) {
		printf("val=");
		val->println();
		fflush(stdout);
		printf("mult=");
		mult->println();
		fflush(stdout);
		}
	return OK;
}

#if TEXDOCU
INT sprint_class_structure(VECTOR_OP class_len, 
	VECTOR_OP class_el_ord, BYTE *str, INT f_v)
#endif
{
	/* VECTOR_OB type; */
	VECTOR_OB val, mult;
	VECTOR_OB col3;
	MATRIX_OB Col;
	INT n, i, l, a, b, c, m;

	color_join(class_len, class_el_ord, &col3, &Col, FALSE /* f_v */, FALSE /* f_vv */);

	col3.multiplicities(&val, &mult);
	
	/* val.sprint_multiplicities(&mult, str); */
	
	l = val.s_li();
	for (i = 0; i < l; i++) {
		c = val.s_ii(i);
		/* b = c % n;
		c -= b;
		a = c / n; */
		a = Col.s_iji(c, 0);
		b = Col.s_iji(c, 1);
		m = mult.s_ii(i);
		sprintf(str + strlen(str), "%ld x %ld/%ld", m, a, b);
		if (i < l - 1)
			strcat(str, ", ");
		}
	if (f_v) {
		printf("class structure: mult x cass-length "
			"/ element-order: \n%s\n", str);
		}
	return OK;
}

#endif /* SOLVABLE_TRUE */

