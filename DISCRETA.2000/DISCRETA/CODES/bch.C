/* bch.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1997
 */


#include <DISCRETA/discreta.h>

#ifdef CODES_TRUE

#include <DISCRETA/unip.h>
#include <DISCRETA/gfq.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/codes.h>
#include <DISCRETA/lb.h>
#ifndef DB_INCLUDED
#include <DISCRETA/db.h>
#endif
#include <DISCRETA/fga.h>

#if TEXDOCU
INT unip_mult_gfq(ZECH_DATA *Z, UNIPOLY_OP a, UNIPOLY_OP b, UNIPOLY_OP c)
#else
Multiplies two $GF(q)$ polynomials: c := a * b.
#endif
{
	INT d1, d2, d3, i, j, k, tmp;

	d1 = a->degree();
	d2 = b->degree();
	d3 = d1 + d2;
	c->m_il_n(d3 + 1);
	for (i = 0; i <= d1; i++) {
		for (j = 0; j <= d2; j++) {
			k = i + j;
			tmp = z_mult_num(Z, a->s_ii(i), b->s_ii(j));
			tmp = z_add_num(Z, tmp, c->s_ii(k));
			c->m_ii(k, tmp);
			}
		}
	return OK;
}

#if TEXDOCU
INT unip_init_Xma(ZECH_DATA *Z, INT a, UNIPOLY_OP p)
#else
Initializes the polynomial $p := (X - a)$ where $a$ is an $GF(q)$ element.
Polinomials of this type are used by the minpoly routine.
#endif
{
	INT tmp;
	
	p->m_il(2);
	p->m_ii(1, 1);
	tmp = z_negate_num(Z, a);
	p->m_ii(0, tmp);
	return OK;
}

#define DEBUG_MINPOLY

#if TEXDOCU
INT unip_minpoly(ZECH_DATA *Z, INT a, UNIPOLY_OP m, INT f_v)
#else
Computes the minimum polynomial of the $GF(q)$ element $a$ 
into $m$. The minimumpolynomial is 
$\prod (X - a^\gamma)$ where $\gamma$ ranges over all 
$Gal(GF(q) / GF(p))$ elements. 
This is in fact the set of conjugates, i.e. the orbit under the 
``frobenius'' map $x \mapsto x^p$. This map is stored in 
the ZECH\_DATA structure.
#endif
{
	UNIPOLY_OB u, v;
	INT a_save, a1, i, d, x;

	if (f_v) {
		printf("unip_minpoly() q=%ld q=%ld\n", Z->q, a);
		fflush(stdout);
		}
	if (a == 0) {
		return error("unip_minpoly() a == 0");
		}
	a_save = a;
	unip_init_Xma(Z, a, &u);
	if (f_v) {
		printf("X-%ld = ", a);
		u.println();
		fflush(stdout);
		}
	u.copy(m);
	a1 = a;
	while (TRUE) {
		a1 = z_apply_frob_num(Z, a1);
		if (a1 == a)
			break;
		unip_init_Xma(Z, a1, &u);
		if (f_v) {
			printf("X-%ld = ", a1);
			u.println();
			fflush(stdout);
			}
		unip_mult_gfq(Z, m, &u, &v);
		v.swap(m);
		if (f_v) {
			printf("m=");
			m->println();
			fflush(stdout);
			}
		}

	// test if all coefficients lie in the ground field:
	d = m->degree();
	for (i = 0; i <= d; i++) {
		x = m->s_ii(i);
		if (x >= Z->chi) {
			printf("unip_minpoly() i = %ld x = %ld chi=%ld\n", i, x, Z->chi);
			fflush(stdout);
			return error("stop");
			}
		}
#ifdef DEBUG_MINPOLY
	if (f_v) {
		printf("minpoly of %ld: ", a_save);
		m->println();
		fflush(stdout);
		}
#endif
	return OK;
}

#define MAX_DEG 20

#if TEXDOCU
INT unip_QR_code_generator_polynomial(UNIPOLY_OP g, INT n, INT p, INT f_v, INT f_vv)
#else
#endif
{
	UNIPOLY_OB u, v;
	ZECH_DATA *Z;
	INT m, q, i0, i;
	VECTOR_OB Q, N;
	INT f_first, f_zech_v, l, j, jj, jjj;
	
	if (f_v) {
		printf("unip_QR_code_generator_polynomial() n = %ld p = %ld\n", 
			n, p);
		}
	m = order_mod_p(p, n);
	q = i_power_j(p, m);
	if (m >= MAX_DEG) {
		return error("unip_QR_code_generator_polynomial() m == MAX_DEG");
		}
	i0 = (q - 1) / n;
	if (f_v) {
		printf("unip_QR_code_generator_polynomial() choosing "
			"q = p^m = %ld = %ld^%ld i0 = %ld\n", q, p, m, i0);
		fflush(stdout);
		}
	f_zech_v = f_v;
	if (q > 50)
		f_zech_v = FALSE;
	Z = zech_open(m, p, f_zech_v);
	quadratic_residues(n, &Q, &N);
	if (f_v) {
		printf("quadratic residues mod %ld:\n", n);
		Q.println();
		printf("quadratic non-residues mod %ld:\n", n);
		N.println();
		fflush(stdout);
		}
	f_first = TRUE;
	l = Q.s_li();
	for (i = 0; i < l; i++) {
		j = Q.s_ii(i);
		jj = (j * i0) % (q - 1);
		jjj = Z->Num_inv[jj];
		if (f_v) {
			printf("j = %ld (j * i0) %% (q - 1) = %ld jjj = %ld\n", j, jj, jjj);
			}
		unip_init_Xma(Z, jjj, &u);
		if (f_v) {
			printf("(X-\\alpha^%ld) = (X-(\\alpha^%ld)^%ld) = ", jj, i0, j);
			fflush(stdout);
			u.println();
			fflush(stdout);
			}
		if (f_first) {
			u.swap(g);
			f_first = FALSE;
			}
		else {
			unip_mult_gfq(Z, g, &u, &v);
			v.swap(g);
			}
		if (f_v) {
			printf("g = ");
			g->println();
			fflush(stdout);
			}
		}
	
	zech_free(Z);
	return OK;
}

#if TEXDOCU
INT unip_BCH_generator_polynomial(UNIPOLY_OP g, INT n, INT designed_distance, 
	INT *bose_distance, INT p, INT f_v, INT f_vv)
#else
Constructs a generator polynomial for a (narrow sense) BCH code.
The polynomial is (with $\delta:=$ designed\_distance)
\[
lcm(m_b, m_{b+1}, \ldots, m_{b+\delta-2}) =: g
\]
(designed\_distance -1 terms). Here, $m_i$ denotes the 
minimum polynomial of the element $\alpha^i \in GF(q) = GF(p^m)$.
$m$ is so chosen such that $n$ divides $q^m - 1$ i.e. 
$GF(p^m)$ contains $n$-th roots of unity.
$b=1$ leads to narrow sense BCH codes.
We know the following:
\begin{enumerate}
\item
$d \ge \delta$
\item
$k \ge n - m(\delta - 1)$
\end{enumerate}
We remind the reader that the generator polynomial is in fact a 
polynomial over $GF(p)$ (and not $GF(q)$). 
The extension field is used only for the construction process. 
Afterwards, all elements lie in $GF(p)$. 
Because we use numeric numbering of $GF(q)$ elements, 
the numbers $0,1,\ldots,p-1$ stand for the $GF(p)$ elements. 
This way, we get a polynomial with coefficients 
only in $0,1,\ldots,p-1$. 
The bose distance is returned. This is the largest 
$d \ge \delta$ such that BCH(n,d,p) gives the same generator 
polynomial (and the same code). So, the true minimum distance is 
$\ge d$.
#endif
{
	UNIPOLY_OB u, v;
	ZECH_DATA *Z;
	INT m, q, i, ai, bi;
	VECTOR_OB sel, sel1;
	INT f_first;
	
	if (f_v) {
		printf("unip_BCH_generator_polynomial() n = %ld delta = %ld p = %ld\n", 
			n, designed_distance, p);
		}
#if 0
	for (m = 1; m < MAX_DEG; m++) {
		q = i_power_j(p, m);
		if ((q - 1) % n == 0)
			break;
		}
#endif
	m = order_mod_p(p, n);
	q = i_power_j(p, m);
	if (m >= MAX_DEG) {
		return error("unip_BCH_generator_polynomial() m == MAX_DEG");
		}
	if (f_v) {
		printf("unip_BCH_generator_polynomial() choosing "
			"q = p^m = %ld = %ld^%ld\n", q, p, m);
		fflush(stdout);
		}
	if (1 + designed_distance - 2 >= q - 1)
		return error("unip_BCH_generator_polynomial() "
			"1 + designed_distance - 2 >= q - 1");
	Z = zech_open(m, p, f_vv);
	sel.m_il_n(q);
	sel1.m_il_n(q);
	for (i = 1; i <= 1 + designed_distance - 2; i++) {
		ai = Z->Num_inv[i];
		if (f_v) {
			printf("i = %ld, ai = %ld\n", i, ai);
			fflush(stdout);
			}
		if (sel.s_ii(ai))
			continue;
		if (f_v) {
			printf("no yet used\n");
			}
		sel.m_ii(ai, 1);
		sel1.m_ii(ai, 1);
		bi = z_apply_frob_num(Z, ai);
		if (f_v) {
			printf("orbit of conjugate elements:\n");
			printf("{ %ld ", ai);
			}
		while (bi != ai) {
			if (f_v) {
				printf("%ld ", bi);
				}
			sel.m_ii(bi, 1);
			bi = z_apply_frob_num(Z, bi);
			}
		if (f_v) {
			printf("}\n");
			}
		}
	for (; i <= q - 2; i++) {
		ai = Z->Num_inv[i];
		if (sel.s_ii(ai)) 
			continue;
		break;
		}
	i--;
	// now: i is the last consecutive element contained 
	// in the roots. 
	// So, we have $a^{1+0},a^{1+1},\ldots,a^{1+(i+1)-2}=a^{i}$
	// in the roots of the code. 
	// this gives (i+1)-2 = bose_distance - 2
	*bose_distance = i + 1;
	if (f_v) {
		printf("bose_distance = %ld\n", *bose_distance);
		fflush(stdout);
		}
	f_first = TRUE;
	for (i = 1; i <= q - 1; i++) {
		if (sel1.s_ii(i) == 0)
			continue;
		unip_minpoly(Z, i, &u, f_vv);
		if (f_v) {
			printf("minpoly %ld = ", i);
			fflush(stdout);
			u.println();
			fflush(stdout);
			}
		if (f_first) {
			u.swap(g);
			f_first = FALSE;
			}
		else {
			unip_mult_gfq(Z, g, &u, &v);
			v.swap(g);
			}
		if (f_v) {
			printf("g = ");
			g->println();
			fflush(stdout);
			}
		}
	zech_free(Z);
	return OK;
}

#if TEXDOCU
INT generator_mat_cyclic_code(UNIPOLY_OP g, INT n, MATRIX_OP G)
#else
Computes the generator matrix for the cyclic code 
of length n generated by the polynomial $g$. This is 
a $k \times n$ matrix where $k:=n-deg(g)$.
It contains the coefficients of $g$ in its rows (shifted cyclically).
#endif
{
	INT d, k, i, j, a;
	
	d = g->degree();
	k = n - d;
	G->m_ilih_n(n, k);
	for (i = 0; i < k; i++) {
		for (j = 0; j <= d; j++) {
			a = g->s_ii(j);
			G->m_iji(i, i + j, a);
			}
		}
	return OK;
}

/* 
 * CODE
 */

INT CODE_OB::field_name(INT i, INT j, BYTE *str)
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "id"; break;
	case 1: s = "n"; break;
	case 2: s = "k"; break;
	case 3: s = "q"; break;
	case 4: s = "d"; break;
	case 5: s = "ago"; break;
	case 6: s = "ago_text"; break;
	case 7: s = "M"; break;
	case 8: s = "W"; break;
	case 9: s = "T"; break;
	default:
		return error("CODE::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

INT CODE_OB::Init(INT n, INT k, INT q)
{
	INT erg = OK;
	
	erg += m_il(10);
	c_obj_k(CODE);
	
	s_id()->m_i(-1);
	s_n()->m_i(n);
	s_k()->m_i(k);
	s_q()->m_i(q);
	s_d()->m_i(0);
	s_ago()->m_i_i(0);
	s_ago_text()->init("                    ");
	s_W()->m_il(n + 1);
	s_T()->m_il(0);
	return erg;
}

INT CODE_OB::Init_M(MATRIX_OP M)
{
	INT erg = OK;
	
	erg += M->copy(s_M());
	return erg;
}

INT CODE_OB::sprint(BYTE *s)
{
	BYTE str1[256];
	BYTE str2[256];
	
	sprintf(str1, "CODE: n = %ld "
		"k = %ld q = %ld d = %ld id=%ld", 
		s_n_i(), s_k_i(), s_q_i(), s_d_i(), s_id_i());
	if (!s_ago()->emptyp()) {
		str2[0] = 0;
		s_ago()->sprint(str2);
		sprintf(str1 + strlen(str1), " ago=");
		strcat(str1, str2);
		}
	if (strlen(s) + strlen(str1) < 200)
		strcat(s, str1);
	return OK;
}
 
#if TEXDOCU
INT CODE_OB::Print()
#endif
{
	println();
	s_M()->fprint_raw(stdout);
	return OK;
}

#if TEXDOCU
INT CODE_OB::fprint(FILE *fp)
#endif
{
	fprintf(fp, "%ld %ld %ld %ld\n", s_n_i(), s_k_i(), s_q_i(), s_d_i());
	s_M()->fprint_raw(fp);
	return OK;
}

#if TEXDOCU
INT CODE_OB::set_ago(VECTOR_OP aut_gen)
#endif
{
	SYM_OB go;
	LABRA_OB L;
	BYTE str[1024], str1[1024];
	INT i, l;
	
	reduce_generators_labra(aut_gen, &go, FALSE /* f_verbose */, &L);
	str[0] = 0;
	str1[0] = 0;
	go.sprint(str);
	l = strlen(str);
	for (i = l; i < 20; i++) {
		strcat(str1, " ");
		}
	strcat(str1, str);
	printf("CODE_OB::set_ago: |Aut|=%s\n", str1);
	printf("$|Aut(C)| = %s$\\\\\n", str1);
	
	s_ago_text()->init(str1);
	go.copy(s_ago());
	return OK;
}

#if TEXDOCU
INT CODE_OB::init_map(INT n, INT k, INT q, INT *map, PERM_REP *P)
#endif
{
   	INT *vec1;
	INT i, j;
	MATRIX_OB M;
	
	if (n <= 0)
		return error("CODE_OB::init_map() n <= 0");
	if (k <= 0)
		return error("CODE_OB::init_map() k <= 0");
	if (n != s_n_i())
		return error("CODE_OB::init_map() n is wrong");
	if (k != s_k_i())
		return error("CODE_OB::init_map() k is wrong");
	if (q != s_q_i())
		return error("CODE_OB::init_map() q is wrong");
	vec1 = (INT *) my_malloc(sizeof(INT) * k, "CODE_OB::init_map");
	M.m_ilih(n, k);
	for (j = 0; j < n; j++) {
		perm_rep_i2vec(P, map[j], vec1);
		for (i = 0; i < k; i++) {
			M.m_iji(i, j, vec1[i]);
			}
		}
	M.swap(s_M());
	my_free(vec1);
	return OK;
}


#if TEXDOCU
INT CODE_OB::mindist()
#endif
{
	INT d;
	
	d = code_Mindist_of_matrix(s_M(), s_q_i());
	s_d()->m_i(d);
	return OK;
}

#if TEXDOCU
INT i_db_codes(DATABASE_OP *p_db, BYTE *path)
#else
Creates a database object for the database of codes.
The btree access paths are: \\
btree 0: id, n, k, d, ago \\
btree 1: n, k, d, ago, id \\
btree 2: n, d, k, ago, id \\
btree 3: d, k, n, ago, id \\
btree 4: k, d, n, ago, id \\
btree 5: ago, n, k, d, id \\
path with trailing stash
#endif
{
	DATABASE_OP db;
	BAYERTREE_OP bt;
	BYTE file_db[1024];
	BYTE file_idx0[1024];
	BYTE file_idx1[1024];
	BYTE file_idx2[1024];
	BYTE file_idx3[1024];
	BYTE file_idx4[1024];
	BYTE file_idx5[1024];

	strcpy(file_db, path);
	strcat(file_db, "codes.db");
	strcpy(file_idx0, path);
	strcat(file_idx0, "codes0.idx");
	strcpy(file_idx1, path);
	strcat(file_idx1, "codes1.idx");
	strcpy(file_idx2, path);
	strcat(file_idx2, "codes2.idx");
	strcpy(file_idx3, path);
	strcat(file_idx3, "codes3.idx");
	strcpy(file_idx4, path);
	strcat(file_idx4, "codes4.idx");
	strcpy(file_idx5, path);
	strcat(file_idx5, "codes5.idx");
	db = (DATABASE_OP) callocobject("i_db_codes()");
	db->init(file_db, CODE);

	db->add_btree(file_idx0, 
		TRUE /* f_duplicatekeys */, 
		0 /* btree_idx */ );
	bt = db->s_btree_i(0);
	bt->add_key_INT4(CODE_INDEX_ID, 0);
	bt->add_key_INT4(CODE_INDEX_N, 0);
	bt->add_key_INT4(CODE_INDEX_K, 0);
	bt->add_key_INT4(CODE_INDEX_D, 0);
	bt->add_key_STRING(20, CODE_INDEX_AGO, 0);

	db->add_btree(file_idx1, 
		TRUE /* f_duplicatekeys */, 
		1 /* btree_idx */ );
	bt = db->s_btree_i(1);
	bt->add_key_INT4(CODE_INDEX_N, 0);
	bt->add_key_INT4(CODE_INDEX_K, 0);
	bt->add_key_INT4(CODE_INDEX_D, 0);
	bt->add_key_STRING(20, CODE_INDEX_AGO, 0);
	bt->add_key_INT4(CODE_INDEX_ID, 0);

	db->add_btree(file_idx2, 
		TRUE /* f_duplicatekeys */, 
		2 /* btree_idx */ );
	bt = db->s_btree_i(2);
	bt->add_key_INT4(CODE_INDEX_N, 0);
	bt->add_key_INT4(CODE_INDEX_D, 0);
	bt->add_key_INT4(CODE_INDEX_K, 0);
	bt->add_key_STRING(20, CODE_INDEX_AGO, 0);
	bt->add_key_INT4(CODE_INDEX_ID, 0);

	db->add_btree(file_idx3, 
		TRUE /* f_duplicatekeys */, 
		3 /* btree_idx */ );
	bt = db->s_btree_i(3);
	bt->add_key_INT4(CODE_INDEX_D, 0);
	bt->add_key_INT4(CODE_INDEX_K, 0);
	bt->add_key_INT4(CODE_INDEX_N, 0);
	bt->add_key_STRING(20, CODE_INDEX_AGO, 0);
	bt->add_key_INT4(CODE_INDEX_ID, 0);

	db->add_btree(file_idx4, 
		TRUE /* f_duplicatekeys */, 
		4 /* btree_idx */ );
	bt = db->s_btree_i(4);
	bt->add_key_INT4(CODE_INDEX_K, 0);
	bt->add_key_INT4(CODE_INDEX_D, 0);
	bt->add_key_INT4(CODE_INDEX_N, 0);
	bt->add_key_STRING(20, CODE_INDEX_AGO, 0);
	bt->add_key_INT4(CODE_INDEX_ID, 0);

	db->add_btree(file_idx5, 
		TRUE /* f_duplicatekeys */, 
		5 /* btree_idx */ );
	bt = db->s_btree_i(5);
	bt->add_key_STRING(20, CODE_INDEX_AGO, 0);
	bt->add_key_INT4(CODE_INDEX_N, 0);
	bt->add_key_INT4(CODE_INDEX_K, 0);
	bt->add_key_INT4(CODE_INDEX_D, 0);
	bt->add_key_INT4(CODE_INDEX_ID, 0);

	*p_db = db;
	return OK;
}

#if TEXDOCU
void e_db_codes(DATABASE_OP db)
#else
Frees the object.
#endif
{
	freeall(db);
}

#if TEXDOCU
INT db_codes_create_and_close(BYTE *path)
#else
Creates an empty database and closes it.
#endif
{
	DATABASE_OP db;
	
	if (i_db_codes(&db, path) != OK) {
		return error("db_codes_create_and_close() i_db_codes()");
		}
	if (db->create() != OK) {
		return error("db_codes_create_and_close() create()");
		}
	if (db->close() != OK) {
		return error("db_codes_create_and_close() close()");
		}
	e_db_codes(db);
	return OK;
}

#if TEXDOCU
INT db_codes_info(BYTE *path, INT btree_idx)
#else
Calls the btree function \lq print\_pages\rq for the btree index \lq btree\_idx\rq.
#endif
{
	BAYERTREE_OP btree;
	DATABASE_OP db;
	
	if (i_db_codes(&db, path) != OK) {
		return error("db_codes_info() i_db_codes()");
		}
	/* da_print(db, 0); */
	btree = db->s_btree_i(btree_idx);
	btree->print_pages();
	e_db_codes(db);
	return OK;
}

#if TEXDOCU
INT db_codes_dump(BYTE *path, INT btree_idx, BYTE *fname)
#endif
{
	BAYERTREE_OP btree;
	DATABASE_OP db;
	INT i, l;
	KEYTYPE key;
	DATATYPE data;
	CODE_OB C;
	FILE *fp = NIL;
	
	if (fname) {
		fp = fopen(fname, "w");
		}
	if (i_db_codes(&db, path) != OK) {
		return error("db_codes_dump() i_db_codes()");
		}
	btree = db->s_btree_i(btree_idx);
	//btree->print_pages();
	
	if (db->open() != OK) {
		return error("db_codes_dump() error in DB::open()");
		}
	if (btree->len(&l) != OK) {
		return error("db_codes_dump() error in BT::len()");
		}
	for (i = 0; i < l; i++) {
		if (btree->ith(i, &key, &data) != OK) {
			return error("db_codes_dump() error in BT::ith()");
			}
		if (db->get_op(&data, &C) != OK) {
			return error("db_codes_dump() error in DB::get_op()");
			}
		printf("%ld: ", i);
		C.println();
		fflush(stdout);
		if (fp)
			C.fprint(fp);
		} // next i
	
	if (db->close() != OK) {
		return error("db_codes_dump() error in DB::close()");
		}
	e_db_codes(db);
	if (fp) {
		fclose(fp);
		}
	return OK;
}

#if TEXDOCU
INT db_codes_highest_id(DATABASE_OP db, INT *highest_id)
#else
returns the highest id which is used in the database.
#endif
{
	BAYERTREE_OP btree;
	INT btree_idx;
	KEYTYPE key;
	INT4 *pi;
	INT id, len;
	DATATYPE data;

	btree_idx = 0;
	btree = db->s_btree_i(btree_idx);
	if (btree->open() != OK) {
		return error("db_codes_highest_id() error in BT::open()");
		}
	if (btree->len(&len) != OK) {
		return error("db_codes_highest_id() error in BT::len()");
		}
	if (len == 0) {
		btree->close();
		*highest_id = -1;
		return OK;
		}
	if (btree->ith(len - 1, &key, &data) != OK) {
		return error("db_codes_highest_id() error in BT::ith()");
		}
	pi = (INT4 *) &key;
	id = pi[0];
	*highest_id = id;
	if (btree->close() != OK) {
		return error("db_codes_highest_id() error in BT::close()");
		}
	return OK;
}

#if TEXDOCU
INT db_codes_add(DATABASE_OP db, CODE_OP C)
#endif
{
	INT highest_id;
	SYM_OB go;
	INT f_v = FALSE;
	INT f_vv = FALSE;
	
	db_codes_highest_id(db, &highest_id);
	C->s_id()->m_i(highest_id + 1);

	printf("adding code to database:\n");
	C->println();
	fflush(stdout);

	if (db->open_DB() != OK) {
		return error("db_codes_add() error in DB::open_DB()");
		}
	printf("database opened\n");
	fflush(stdout);
	if (db->add_op1(C, f_v, f_vv) != OK) {
		return error("db_codes_add() error in DB::add_op()");
		}
	printf("code added\n");
	fflush(stdout);
	if (db->close_DB() != OK) {
		return error("db_codes_add() error in DB::close()");
		}
	printf("database closed\n");
	fflush(stdout);
	return OK;
}

#endif /* CODES_TRUE */

