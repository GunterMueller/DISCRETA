/* dimino.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef DIMINO_TRUE

#include <DISCRETA/dimino.h>

INT gruppen_elemente(VECTOR_OP G_gen, VECTOR_OP G, INT type, void *data)
/* mit Reduktion des Erzeugendensystems */
{
	VECTOR_OB Ggen;
	INT i;
	
	Ggen.m_il(0);
	G->m_il(0);
	for (i = 0; i < G_gen->s_li(); i++) {
		dimino_extend(G, &Ggen, G_gen->s_i(i), NIL, 0, NIL, type, data);
		}
	Ggen.swap(G_gen);
	return OK;
}

INT gruppen_elemente1(VECTOR_OP G, INT type, void *data)
{
	VECTOR_OB G1, G2;
	INT erg = OK, i;
	
	erg += G->copy(&G1);
	G2.m_il(0);
	G->m_il(0);
	for (i = 0; i < G1.s_li(); i++)
		erg += dimino_extend(G, &G2, G1.s_i(i), NIL, 0, NIL, type, data);
	return erg;
}

INT dimino_extend(VECTOR_OP G, VECTOR_OP G_gen, 
	SYM_OP g, VECTOR_OP bad_list, INT bad_list_len, 
	INT *f_bad_element_found, 
	INT type, void *data)
{
	VECTOR_OB G_tmp, reps;
	SYM_OP g1;
	SYM_OB id, h, h1;
	INT erg = OK, nb_G;
	INT j, k, l, i1, l1, idx, f_found;
	
	if (bad_list)
		*f_bad_element_found = FALSE;
	g->copy(&id);
	do_one(&id, type, data);
	if (G->emptyp() || G->s_li() == 0) {
		G->m_il(1);
		id.copy(G->s_i(0));
		}
	erg += G->search(G->s_li(), TRUE /* f_ascending */, g, &idx, &f_found);
	if (f_found)
		return OK;
	if (G_gen->emptyp())
		G_gen->m_il(1);
	else
		G_gen->inc();
	g->copy(G_gen->s_i(G_gen->s_li() - 1));
	erg += G->copy(&G_tmp);
	reps.m_il(1);
	erg += id.copy(reps.s_i(0));
	/* the dimino induction step: 
	 * apply all generators 
	 * to all coset reps. */
	for (j = 0; j < reps.s_li(); j++) {
		for (i1 = 0; 
			i1 < G_gen->s_li(); i1++) {
			if (j == 0 && 
				i1 < G_gen->s_li() - 1)
				continue;
			g1 = G_gen->s_i(i1);
			do_mult(reps.s_i(j), g1, &h, type, data);
			l1 = G->s_li();
			erg += G->search(l1, TRUE /* f_ascending */, &h, &idx, &f_found);
			if (f_found)
				continue;
			reps.inc();
			erg += h.copy(reps.s_i(reps.s_li() - 1));
			erg += G->v_realloc(l1 + G_tmp.s_li());
			for (k = 0; k < G_tmp.s_li(); k++) {
				do_mult(G_tmp.s_i(k), &h, &h1, type, data);
				if (bad_list) {
					bad_list->search(bad_list_len, TRUE, 
						&h, &idx, &f_found);
					if (f_found) {
						*f_bad_element_found = TRUE;
						return OK;
						}
					}
				/* G hat jetzt l1 + k 
				 * benuetzte Elemente */
				erg += G->search(l1 + k, TRUE /* f_ascending */, 
					&h1, &idx, &f_found);
				if (f_found) {
					printf("dimino_extend(): h1 found\n");
					printf("k = %ld\n", k);
					printf("G_tmp[k] = ");
					G_tmp.s_i(k)->println();
					G_tmp.println();
					printf("h = ");
					h.println();
					printf("h1 = ");
					h1.println();
					G->println();
					return error("..");
					}
				for (l = l1 + k; l > idx; l--) {
					G->s_i(l - 1)->swap(G->s_i(l));
					}
				erg += h1.copy(G->s_i(idx));
				} /* next k */
			} /* next i1 */
		} /* next j */
	return OK;
}

INT dimino_extend_normal(VECTOR_OP G, VECTOR_OP G_gen, SYM_OP g, INT type, void *data)
/* before a10_dimino_extend_normal() */
{
	VECTOR_OB G_tmp;
	SYM_OB g1, g2, id, h;
	INT erg = OK, nb_G, k;
	INT l, l1, idx, f_found;
	
	g->copy(&id);
	do_one(&id, type, data);
	if (G->emptyp() || G->s_li() == 0) {
		G->m_il(1);
		id.copy(G->s_i(0));
		}
	erg += G->search(G->s_li(), TRUE /* f_ascending */, g, &idx, &f_found);
	if (f_found)
		return OK;
	if (G_gen->emptyp() || G_gen->s_li() == 0)
		G_gen->m_il(1);
	else
		G_gen->inc();
	g->copy(G_gen->s_i(G_gen->s_li() - 1));
	erg += G->copy(&G_tmp);
	g->copy(&g1);
	while (TRUE) {
		if (do_einsp(&g1, type, data))
			break;
		if (G_tmp.search(G_tmp.s_li(), TRUE, 
			&g1, &idx, &f_found) != OK)
			return error("dimino_extend_normal(): error in VEC::search()");
		if (f_found)
			break;
		l1 = G->s_li();
		G->v_realloc(l1 + G_tmp.s_li());
		for (k = 0; k < G_tmp.s_li(); k++) {
			do_mult(G_tmp.s_i(k), &g1, &h, type, data);
			/* G hat jetzt l1 + k 
			 * benuetzte Elemente */
			erg += G->search(l1 + k, TRUE /* f_ascending */, &h, &idx, &f_found);
			if (f_found)
				return error("dimino_extend_normal(): h1 found\n");
			for (l = l1 + k; l > idx; l--) {
				G->s_i(l - 1)->swap(G->s_i(l));
				}
			erg += h.copy(G->s_i(idx));
			} /* next k */
		if (do_mult(&g1, g, &g2, type, data) != OK)
			return error("dimino_extend_normal(): error in do_mult()");
		g2.swap(&g1);
		}
	return OK;
}

INT RelativeOrder(VECTOR_OP H, SYM_OP g, INT *relative_order, INT type, void *data)
{
	SYM_OB g1, g2;
	INT l, idx, f_found;
	
	g->copy(&g1);
	l = 1;
	while (TRUE) {
		if (do_einsp(&g1, type, data))
			break;
		if (H->search(H->s_li(), TRUE, &g1, &idx, &f_found) != OK)
			return error("RelativeOrder(): error in VEC::search()");
		if (f_found)
			break;
		if (do_mult(g, &g1, &g2, type, data) != OK)
			return error("RelativeOrder(): error in do_mult()");
		g2.swap(&g1);
		l++;
		}
	*relative_order = l;
	return OK;
}

#endif /* DIMINO_TRUE */

