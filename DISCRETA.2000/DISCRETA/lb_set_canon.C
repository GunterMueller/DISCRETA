// lb_set_canon.C


typedef struct set_canon SET_CANON;

struct set_canon {
	LABRA_OP G;
	INT nb_x;
	INT deg;
	
	INT *X;
	INT *chi;
	
	INT base_len;
	INT *base_idx;
	INT **base_orbit;
	INT *base_orbit_len;
	INT max_orbit_length;
	
	INT *path;
	INT *path0;
	INT *X0;
	INT *chi0;
	
	INT *E;
	INT *e_len;
	INT *e_idx;
	
	INT nb_aut;
	INT f_v;
	INT f_vv;
};

INT calc_canonical_form_of_set(LABRA_OP G, INT nb_x, INT *X, INT *Y, INT f_v, INT f_vv);
void free_set_canon(SET_CANON *sc);
INT init_set_canon(SET_CANON **p_sc, LABRA_OP G, INT nb_x, INT *X, INT f_v, INT f_vv);
INT set_canon_apply_sigma_ij(SET_CANON *sc, INT i, INT j, INT f_v);
INT set_canon_apply_sigma_ij_inv(SET_CANON *sc, INT i, INT j, INT f_v);
INT set_canon_apply(SET_CANON *sc, INT i, INT j, PERMUTATION_OP p, INT f_v);
INT set_canon_classify(SET_CANON *sc, INT i, INT f_v, INT f_vv);
INT set_canon_backtrack(SET_CANON *sc, INT i);


INT calc_canonical_form_of_set(LABRA_OP G, INT nb_x, INT *X, INT *Y, INT f_v, INT f_vv)
{
	SET_CANON *sc;
	INT s = f_perm_print_start_with_zero;
	INT i, ago;
	
	f_perm_print_start_with_zero = TRUE;
	init_set_canon(&sc, G, nb_x, X, f_v, f_vv);
	set_canon_backtrack(sc, 0);
	
	ago = sc->nb_aut;
	for (i = 0; i < sc->nb_x; i++) {
		Y[i] = sc->X0[i];
		}
	f_perm_print_start_with_zero = s;
	if (f_v) {
		printf("the canonical set is: {");
		for (i = 0; i < sc->nb_x; i++) {
			printf("%ld ", sc->X0[i]);
			}
		printf("}\n");
		printf("nb_aut=%ld\n", sc->nb_aut);
		}

	free_set_canon(sc);
	return ago;
}

void free_set_canon(SET_CANON *sc)
{
	INT i;
	
	my_free(sc->X);
	my_free(sc->chi);
	
	for (i = 0; i < sc->base_len; i++)
		my_free(sc->base_orbit[i]);
	my_free(sc->base_idx);
	my_free(sc->base_orbit);
	my_free(sc->base_orbit_len);
	
	my_free(sc->path);
	my_free(sc->path0);
	my_free(sc->X0);
	my_free(sc->chi0);
	
	my_free(sc->E);
	my_free(sc->e_len);
	my_free(sc->e_idx);

	my_free(sc);
}

INT init_set_canon(SET_CANON **p_sc, LABRA_OP G, INT nb_x, INT *X, INT f_v, INT f_vv)
{
	SET_CANON *sc = NULL;
	INT i, j, t, a;
	
	sc = (SET_CANON *) my_malloc(sizeof(SET_CANON), "SET_CANON");
	sc->G = G;
	sc->deg = G->s_degree_i();
	sc->nb_x = nb_x;
	// sc->X = X;
	sc->base_len = G->s_tidx()->s_li();

	sc->X = (INT *) my_malloc((sc->base_len + 1) * sc->nb_x * sizeof(INT), "X");
	sc->chi = (INT *) my_malloc((sc->base_len + 1) * sc->deg * sizeof(INT), "X");
	for (i = 1; i < sc->nb_x; i++) {
		if (X[i - 1] >= X[i])
			return error("init_set_canon(): X[i - 1] >= X[i]");
		}
	for (i = 0; i < sc->nb_x; i++) {
		sc->X[0 * sc->nb_x + i] = X[i];
		}
	for (i = 0; i < sc->deg; i++) {
		sc->chi[i] = 0;
		}
	for (i = 0; i < sc->nb_x; i++) {
		sc->chi[X[i]] = 1;
		}


	sc->base_idx = (INT *) my_malloc((sc->base_len + 1) * sizeof(INT), "base_idx");
	sc->base_orbit = (INT **) my_malloc(sc->base_len * sizeof(INT *), "base_orbit");
	sc->base_orbit_len = (INT *) my_malloc(sc->base_len * sizeof(INT), "base_orbit_len");
	for (i = 0; i < sc->base_len; i++) {
		sc->base_idx[i] = t = G->s_tidx_ii(i);
		sc->base_orbit_len[i] = G->s_T_i(t)->s_li();
		sc->base_orbit[i] = (INT *) my_malloc(sc->base_orbit_len[i] * sizeof(INT *), "base_orbit[i]");
		for (j = 0; j < sc->base_orbit_len[i]; j++) {
			sc->base_orbit[i][j] = G->s_T_iji(t, j);
			}
		if (sc->base_orbit[i][0] != i)
			return error("init_set_canon() sc->base_orbit[i][0] != i");
		}
	if (sc->base_idx[0] != 0)
		return error("init_set_canon() sc->base_idx[0] != 0");
	sc->base_idx[sc->base_len] = sc->deg;

	sc->max_orbit_length = 0;
	for (i = 0; i < sc->base_len; i++) {
		j = sc->base_orbit_len[i];
		sc->max_orbit_length = MAX(sc->max_orbit_length, j);
		}

	sc->path = (INT *) my_malloc(sc->base_len * sizeof(INT), "path");
	sc->path0 = (INT *) my_malloc(sc->base_len * sizeof(INT), "path0");
	sc->X0 = (INT *) my_malloc(sc->nb_x * sizeof(INT), "X0");
	sc->chi0 = (INT *) my_malloc(sc->deg * sizeof(INT), "chi0");
	for (i = 0; i < sc->base_len; i++) {
		sc->path0[i] = 0;
		}
	for (i = 0; i < sc->deg; i++) {
		sc->chi0[i] = 0;
		}
	for (i = 0; i < sc->nb_x; i++) {
		a = X[i];
		sc->X0[i] = a;
		sc->chi0[a] = 1;
		}
	
	sc->E = (INT *) my_malloc(sc->max_orbit_length * sc->base_len * sizeof(INT), "E");
	sc->e_len = (INT *) my_malloc(sc->base_len * sizeof(INT), "e_len");
	sc->e_idx = (INT *) my_malloc(sc->base_len * sizeof(INT), "e_idx");
	
	sc->nb_aut = 0;
	sc->f_v = f_v;
	sc->f_vv = f_vv;
	
	*p_sc = sc;
	return OK;
}

INT set_canon_apply_sigma_ij(SET_CANON *sc, INT i, INT j, INT f_v)
{
	PERMUTATION_OB p;
	
	sc->G->rep_ij(i, j, &p);
	set_canon_apply(sc, i, j, &p, f_v);
	return OK;
}

INT set_canon_apply_sigma_ij_inv(SET_CANON *sc, INT i, INT j, INT f_v)
{
	PERMUTATION_OB p, q;
	
	sc->G->rep_ij(i, j, &p);
	p.invers(&q);
	set_canon_apply(sc, i, j, &q, f_v);
	return OK;
}

INT set_canon_apply(SET_CANON *sc, INT i, INT j, PERMUTATION_OP p, INT f_v)
{
	INT jj, a, b, i1, i2, i3, k;
	
	i1 = i * sc->nb_x;
	i2 = i1 + sc->nb_x;
	i3 = (i + 1) * sc->deg;
	for (jj = 0; jj < sc->deg; jj++) {
		sc->chi[i3 + jj] = 0;
		}
	if (f_v) {
		printf("set_canon_apply() i=%ld: j=%ld p=", i, j);
		p->println();
		printf("X={");
		for (jj = 0; jj < sc->nb_x; jj++) {
			a = sc->X[i1 + jj];
			printf("%ld ", a);
			}
		printf("} -> {");
		}
	for (jj = 0; jj < sc->nb_x; jj++) {
		a = sc->X[i1 + jj];
		b = p->s_ii(a) - 1;
		sc->chi[i3 + b] = 1;
		}
	k = 0;
	for (jj = 0; jj < sc->deg; jj++) {
		if (sc->chi[i3 + jj]) {
			sc->X[i2 + k++] = jj;
			}
		}
	if (f_v) {
		for (jj = 0; jj < sc->nb_x; jj++) {
			a = sc->X[i2 + jj];
			printf("%ld ", a);
			}
		printf("}\n");
		fflush(stdout);
		}
	if (k != sc->nb_x) {
		printf("k=%ld nb_x=%ld\n", k, sc->nb_x);
		return error("set_canon_apply() k != sc->nb_x");
		}
	return OK;
}

INT set_canon_classify(SET_CANON *sc, INT i, INT f_v, INT f_vv)
{
	INT len, j, k, k1, idx1, idx2, ii, a, b, t_len, kk;
	VECTOR_OB S, T, Tidx, Type, Idx, Type0;
	VECTOR_OP p_idx;
	INT idxx, f_found;
	INT f_equal_beginning, cmp;
	
	if (f_v) {
		printf("set_canon_classify() i=%ld\n", i);
		}
	idx1 = sc->base_idx[i];
	idx2 = sc->base_idx[i + 1];
	S.m_il(0);
	len = sc->base_orbit_len[i];
	for (k = 0; k < len; k++) {
		j = sc->base_orbit[i][k];
		if (sc->chi[i * sc->deg + j])
			S.append_i(k);
		}
	if (S.s_li() == 0) {
		S.m_il(len);
		for (k = 0; k < len; k++)
			S.m_ii(k, k);
		}
	len = S.s_li();
	if (f_v) {
		printf("S=");
		S.println();
		}

	
	f_equal_beginning = TRUE;
	if (f_v) {
		printf("beginning:\n");
		printf("chi0: ");
		for (ii = 0; ii < idx1; ii++) {
			a = sc->chi0[ii];
			printf("%ld", a);
			}
		printf("\n");
		printf("chi:  ");
		for (ii = 0; ii < idx1; ii++) {
			b = sc->chi[i * sc->deg + ii];
			printf("%ld", b);
			}
		printf("\n");
		}
	for (ii = 0; ii < idx1; ii++) {
		a = sc->chi0[ii];
		b = sc->chi[i * sc->deg + ii];
		if (b < a) {
			return error("set_canon_classify() b < a");
			}
		if (b > a) {
			f_equal_beginning = FALSE;
			break;
			}
		}
	if (f_vv) {
		printf("f_equal_beginning=%ld\n", f_equal_beginning);
		}
	if (f_equal_beginning) {
		Type0.m_il(idx2 - idx1);
		for (ii = idx1; ii < idx2; ii++) {
			a = sc->chi0[ii];
			Type0.m_ii(ii - idx1, a);
			}
		if (f_vv) {
			printf("type0=");
			Type0.println();
			}
		}



	T.m_il(0);
	Tidx.m_il(0);
	t_len = 0;
	for (k = 0; k < len; k++) {
		k1 = S.s_ii(k);
		j = sc->base_orbit[i][k1];
		set_canon_apply_sigma_ij_inv(sc, i, j, f_vv);
		if (f_equal_beginning) {
			cmp = 0;
			for (ii = idx1; ii < idx2; ii++) {
				a = sc->chi[(i + 1) * sc->deg + ii];
				b = sc->chi0[ii];
				if (a < b) {
					cmp = -1;
					break;
					}
				if (a > b) {
					cmp = 1;
					break;
					}
				}
			if (cmp < 0)
				continue;
			}
		Type.m_il(idx2 - idx1);
		for (ii = idx1; ii < idx2; ii++) {
			a = sc->chi[(i + 1) * sc->deg + ii];
			Type.m_ii(ii - idx1, a);
			}
		if (f_vv) {
			printf("%ld: type=", k);
			Type.println();
			}
		if (f_equal_beginning) {
			a = Type.sym_comp(&Type0);
			if (f_vv) {
				printf("compare with Type0 gives %ld\n", a);
				}
			if (a < 0) {
				continue;
				}
			}
					
		T.search(t_len, TRUE /* f_ascending */, &Type, 
			&idxx, &f_found);
		if (f_vv) {
			printf("f_found=%ld idxx=%ld\n", f_found, idxx);
			}
		if (f_found) {
			idxx--;
			((VECTOR_OP) Tidx.s_i(idxx))->append_i(k1);
			}
		else {
			T.inc();
			Tidx.inc();
			for (kk = t_len; kk > idxx; kk--) {
				T.s_i(kk)->swap(T.s_i(kk - 1));
				Tidx.s_i(kk)->swap(Tidx.s_i(kk - 1));
				}
			Type.swap((VECTOR_OP) T.s_i(idxx));
			((VECTOR_OP) Tidx.s_i(idxx))->m_il(1);
			((VECTOR_OP) Tidx.s_i(idxx))->m_ii(0, k1);
			t_len++;
			}
		}
	if (f_vv) {
		printf("sorted type vectors:\n");
		for (ii = 0; ii < t_len; ii++) {
			T.s_i(ii)->print();
			printf(" for: ");
			Tidx.s_i(ii)->println();
			}
		}
	if (t_len == 0)
		return FALSE;
	p_idx = (VECTOR_OP) Tidx.s_i(t_len - 1);
	sc->e_len[i] = len = p_idx->s_li();
	for (k = 0; k < len; k++) {
		sc->E[i * sc->max_orbit_length + k] = p_idx->s_ii(k);
		}
	sc->e_idx[i] = 0;
	if (f_v) {
		printf("E_%ld = {", i);
		for (ii = 0; ii < len; ii++) {
			k = sc->E[i * sc->max_orbit_length + ii];
			j = sc->base_orbit[i][k];
			printf("%ld ", j);
			}
		printf("}\n");
		}
	return TRUE;
}

INT set_canon_backtrack(SET_CANON *sc, INT i)
{
	INT ii, j, k, a, b, f_equal;
	
	if (i == sc->base_len) {
		f_equal = TRUE;
		for (ii = 0; ii < sc->deg; ii++) {
			a = sc->chi0[ii];
			b = sc->chi[sc->base_len * sc->deg + ii];
			if (b < a)
				return error("set_canon_backtrack() b < a");
			if (a < b) {
				f_equal = FALSE;
				break;
				}
			}
		if (!f_equal) {
			for (ii = 0; ii < i; ii++) {
				a = sc->e_idx[ii];
				sc->path0[ii] = a;
				}
			for (ii = 0; ii < sc->deg; ii++) {
				sc->chi0[ii] = sc->chi[sc->base_len * sc->deg + ii];
				}
			for (ii = 0; ii < sc->nb_x; ii++) {
				sc->X0[ii] = sc->X[sc->base_len * sc->nb_x + ii];
				}
			sc->nb_aut = 1;
			}
		else {
			sc->nb_aut++;
			if (sc->f_v) {
				printf("found automorphism no %ld\n", sc->nb_aut);
				}
			}
		return OK;
		}
	if (!set_canon_classify(sc, i, sc->f_v, sc->f_vv))
		return OK;
	for (sc->e_idx[i] = 0; sc->e_idx[i] < sc->e_len[i]; sc->e_idx[i]++) {
		k = sc->E[i * sc->max_orbit_length + sc->e_idx[i]];
		j = sc->base_orbit[i][k];
		if (sc->f_v) {
			printf("backtrack i=%ld j=%ld\n", i, j);fflush(stdout);
			}
		set_canon_apply_sigma_ij_inv(sc, i, j, sc->f_v);
		set_canon_backtrack(sc, i + 1);
		}
	return OK;
}


