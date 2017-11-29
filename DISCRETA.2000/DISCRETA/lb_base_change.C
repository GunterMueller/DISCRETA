/* lb_base_change.C: */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef LABRA_TRUE

#include <DISCRETA/lb.h>
#include <DISCRETA/ma.h> // for calc_transversal_matrix()

typedef struct cycle_local {
	INT r, s, n;
	LABRA_OP A, B; // A: the input labra, B the output labra
	VECTOR_OP root;
	VECTOR_OP bottom_merke;

	VECTOR_OP cv; // childrens vector (vector of vectors)
	VECTOR_OP orbit; // orbit of \pi_s under G^(j), r < j <= s
	VECTOR_OP delta; // not allocated, only a pointer into cv !
	VECTOR_OP cosetrep;
} CYCLE_LOCAL;

static INT cycle_top(CYCLE_LOCAL *cl);
static INT cycle_bottom(CYCLE_LOCAL *cl);
static INT cycle_middle(CYCLE_LOCAL *cl, INT f_v);
static INT cycle_middle_init_data(CYCLE_LOCAL *cl);
INT cycle_middle_update_orbit(CYCLE_LOCAL *cl, INT j);

#undef DEBUG_CYCLE

#if TEXDOCU
INT LABRA_OB::cycle(INT r, INT s, INT f_v, INT f_vv)
#else
//PRE
// base change of labelled branchings:
//
// input:
// a labelled branching for the base 
// $\pi_1$, $\pi_2$, ..., $\pi_{r-1}$,  $|$  $\pi_r$, ... , $\pi_s$,          $|$ $\pi_{s+1}$, ..., $\pi_n$
//
// output:
// a labelled branching for the base 
// $\pi_1$, $\pi_2$, ..., $\pi_r-1$,  $|$  $\pi_s$, $\pi_r$, ... , $\pi_{s-1}$, $|$ $\pi_{s+1}$, ..., $\pi_n$
//
// (so we pre-multiply the cycle (r r+1 r+2 .. s) to the current base).
//
// the three parts [1..r-1], [r..s]  and [s+1..n] (separated by '$|$' in the above lists)
// are called the top, middle and bottom parts respectively.
// they are treated in cycle-bottom, cycle-middle and cycle-top in this order.
//
///PRE
#endif
{
	INT i;
	PERMUTATION_OP perm;
	CYCLE_LOCAL *cl;

	cl = (CYCLE_LOCAL *) my_malloc(sizeof(CYCLE_LOCAL), "LABRA::cycle() CYCLE_LOCAL");
	cl->n = s_degree_i();
	cl->r = r;
	cl->s = s;
	cl->A = this;
	cl->B = (LABRA_OP) callocobject("LABRA::cycle() labra B");
	cl->B->Init_no_generators(cl->n, s_f_verbose_i(), s_f_very_verbose_i());
	for (i = 0; i < cl->n; i++) {
		perm = (PERMUTATION_OP) cl->B->s_KM_i(i);
		perm->m_il(cl->n);
		perm->one();
		perm = (PERMUTATION_OP) cl->B->s_EM_i(i);
		perm->m_il(cl->n);
		perm->one();
		}
	cl->root = (VECTOR_OP) callocobject("LABRA::cycle() root");
	cl->bottom_merke = (VECTOR_OP) callocobject("LABRA::cycle() bottom_merke");
	cl->cv = (VECTOR_OP) callocobject("LABRA::cycle() cv");
	cl->orbit = (VECTOR_OP) callocobject("LABRA::cycle() orbit");
	cl->cosetrep = (VECTOR_OP) callocobject("LABRA::cycle() cosetrep");

	if (f_v) {
		printf("LABRA::cycle() r = %ld s = %ld\n", r, s);
		fflush(stdout);
		}
	if (f_vv) {
		printf("before cycle_bottom()\n");
		cl->A->Print();
		printf("the old base:\n");
		cl->A->s_base()->print_list();
		printf("the old inverse base:\n");
		cl->A->s_inv_base()->print_list();
		printf("the old V vector:\n");
		cl->A->s_V()->println();
		fflush(stdout);
		}
	cycle_bottom(cl);
	if (f_vv) {
		printf("after cycle_bottom()\n");
		printf("the new base:\n");
		cl->B->s_base()->print_list();
		printf("the new inverse base:\n");
		cl->B->s_inv_base()->print_list();
		printf("cl->root:\n");
		cl->root->println();
		printf("cl->bottom_merke:\n");
		cl->bottom_merke->println();
		// Print();
		fflush(stdout);
		}
	
	cycle_middle_init_data(cl);
	// orbit of \pi_s under G^(s)
	if (f_vv) {
		printf("after cycle_middle_init_data()\n");
		printf("child vectors:\n");
		for (i = 0; i < cl->n; i++) {
			printf("%ld: ", i);
			cl->cv->s_i(i)->println();
			}
		printf("orbit=child_vector[s]:\n");
		cl->orbit->println();
		printf("cosetrep:\n");
		for (i = 0; i < cl->n; i++) {
			if (cl->orbit->s_ii(i)) {
				printf("%ld: ", i);
				cl->cosetrep->s_i(i)->println();
				}
			}
		printf("calling cycle_middle:\n");
		}
	
	cycle_middle(cl, f_vv);
	if (f_vv) {
		printf("vor cycle_top()\n");
		// Print();
		fflush(stdout);
		}
	cycle_top(cl);
	cl->B->update_EM();
	cl->B->check_consistency();
#ifdef ANGST
	printf("consistency check OK\n");
#endif
	if (f_vv) {
		printf("nach cycle_top()\n");
		cl->B->Print();
		fflush(stdout);
		}
	
	// s_base()->copy(cl->B->s_base());   
	// s_inv_base()->copy(cl->B->s_inv_base());   

	cl->B->swap(cl->A);
	
	freeall(cl->B);
	freeall(cl->root);
	freeall(cl->bottom_merke);
	freeall(cl->cv);
	freeall(cl->orbit);
	freeall(cl->cosetrep);
	cl->B = NIL;
	my_free(cl);
	return OK;
}

static INT cycle_top(CYCLE_LOCAL *cl)
{
	LABRA_OP A = cl->A;
	LABRA_OP B = cl->B;
	INT i, i_el, i_el_old, i_prev, i_prev_el, j, j1, j_el, m;
	PERMUTATION_OB per;
	
#ifdef DEBUG_CYCLE
	printf("cycle_top()\n");
	fflush(stdout);
#endif
	for (i = cl->r; i < cl->n; i++) {
		i_el = B->s_base_ii(i) - 1;
		if (cl->root->s_ii(i)) {
			m = A->s_inv_base_ii(i_el) - 1;
#ifdef DEBUG_CYCLE
			printf("root[i=%ld] = 1, i_el = %ld, m = %ld\n", i, i_el, m);
#endif
			// now: \pi_m = \pi'_i

			// search the predecessor \pi_j of \pi_m with j < r 
			// if there is such an element:	
			j = m;
			while ((j1 = A->s_V_ii(j)) != j && j >= cl->r)
				j = j1;
#ifdef DEBUG_CYCLE
			printf("j = %ld\n", j);
#endif
			if (j < cl->r) {
				// if such an element \pi_j exists:
				B->s_V_i(i)->m_i(j);
#ifdef ANGST
				if (A->s_V_ii(m) == m) {
					if (!A->s_EM_i(i_el)->einsp())
						return error("cycle_top(): !A->s_EM_i(i_el)->einsp()");
					}
#endif
				A->s_EM_i(i_el)->swap(B->s_EM_i(i_el));
				if (A->s_V_ii(m) == j) {
					// m is an immediate descendant of j
					A->s_KM_i(i_el)->swap(B->s_KM_i(i_el));
					}
				else {
					j_el = A->s_base_ii(j) - 1;
					A->s_EM_i(j_el)->invers(&per);
					per.mult(B->s_EM_i(i_el), B->s_KM_i(i_el));
					}
				} // if i < cl->r
			} // if root[i]
		else {
			// root[i] == 0
			i_prev = B->s_V_ii(i);
			i_prev_el = B->s_base_ii(i_prev) - 1;
#ifdef DEBUG_CYCLE
			printf("root[i=%ld] = 0, i_el = %ld\n", i, i_el);
			printf("i_prev = %ld, i_prev_el = %ld\n", i_prev, i_prev_el);
#endif
			B->s_EM_i(i_prev_el)->copy(&per);
			per.mult(B->s_KM_i(i_el), B->s_EM_i(i_el));
			}
		} // next i

	for (i = 0; i < cl->r; i++) {
		i_el = B->s_base_ii(i) - 1;
		i_el_old = A->s_base_ii(i) - 1;
		B->s_V_i(i)->m_i(A->s_V_ii(i));
		A->s_EM_i(i_el_old)->swap(B->s_EM_i(i_el));
		A->s_KM_i(i_el_old)->swap(B->s_KM_i(i_el));
		} // next i
	
#ifdef DEBUG_CYCLE
	printf("end of cycle_top()\n");
#endif
	return OK;
}

static INT cycle_bottom(CYCLE_LOCAL *cl)
{
	LABRA_OP A = cl->A;
	LABRA_OP B = cl->B;
	INT i, j, nb_bm;
	
	cl->root->m_il_n(cl->n);
	for (i = 0; i < cl->n; i++)
		cl->root->m_ii(i, 1);

	for (i = 0; i < cl->n; i++)
		B->s_V_i(i)->m_i(i);
	
	cl->bottom_merke->m_il(cl->n);
	nb_bm = 0;
	
	// the new base:
	// pre-multiply the cycle (r r+1 r+2 ... s)
	A->s_base()->copy(B->s_base());
	B->s_base_i(cl->r)->m_i(A->s_base_ii(cl->s));
	for (i = cl->r + 1; i <= cl->s; i++) {
		B->s_base_i(i)->m_i(A->s_base_ii(i - 1));
		}
	B->s_base()->invers(B->s_inv_base());
	
	for (j = cl->s + 1; j < cl->n; j++) {
		i = A->s_V_ii(j);
		if (i > cl->s && i != j) {
			B->s_V_i(j)->m_i(i);
			cl->root->m_ii(j, 0);
			cl->bottom_merke->m_ii(nb_bm++, j);
			}
		} // next j
	
	cl->bottom_merke->realloc_z(nb_bm);
	return OK;
}

static INT cycle_middle(CYCLE_LOCAL *cl, INT f_v)
{		
	PERMUTATION_OB alpha, alpha_inv, per, per1;
	LABRA_OP A = cl->A;
	LABRA_OP B = cl->B;
	INT i, j, jm1_el, m, m_el, n, p, q, p_el, q_el, nb_bm, j_el, j_el_new, s_el;
	
	n = cl->n;
	
	// loop-invariant:
	// orbit and cosetrep describe the orbit of \pi_s under G^(j)
	// for j = s down to r + 1
	
	for (j = cl->s; j > cl->r; j--) {
		if (f_v) {
			printf("\ncycle_middle() j = %ld\n\n", j);
			fflush(stdout);
			}
		
		jm1_el = A->s_base_ii(j - 1) - 1;
		cl->delta = (VECTOR_OP) cl->cv->s_i(j - 1);
		// delta is the orbit of \pi_j-1.
		// delta[m] = 1 iff \pi_m is a child of \pi_j-1 (in the branching A)
		// this orbit is shrunken to the new orbit of \pi_j in B
		// this is due to the fact that the group in B has to fix \pi_s

		
		// loop over the orbit:
		for (m = j - 1; m < n; m++) {
			m_el = A->s_base_ii(m) - 1;
			// we test if \pi_m lies in the old orbit of \pi_j-1:
			
			if (cl->delta->s_ii(m)) { // if \pi_m is in delta
				if (f_v) {
					printf("m = %ld m_el = %ld in \\Delta=\\pi_{j-1}^{G^{(j-1)}}\n", m, m_el);
					}
				
				A->s_EM_i(jm1_el)->invers(&per);
				per.mult(A->s_EM_i(m_el), &alpha);
				alpha.invers(&alpha_inv);
				// \alpha := EM(\pi_j-1)^{-1} * EM(\pi_m)  
				// (= rep_ij(j-1,m) )
				// maps \pi_j-1 upon \pi_m
#ifdef ANGST
				if (alpha.s_ii(jm1_el) - 1 != m_el)
					return error("cycle_middle() alpha.s_ii(jm1_el) - 1 != m_el");
#endif

				
				// later on we will use the equality
				// \pi'_j^\alpha = \pi_j-1^\alpha = \pi_m = \pi'_q
				// \alpha can be pre-multiplied by an element in G^(j) 
				// without changing equality.
				// this will be done by cosetrep[] later on.

				// now we have to check if \pi_m is a child of \pi_j {\em in B}.
				
				s_el = A->s_base_ii(cl->s) - 1;
				p_el = alpha_inv.s_ii(s_el) - 1;
				// we are going to compute 
				// p_el = s_el^{\alpha^{-1}} 
				//      = s_el^{EM(\pi_m)^{-1} * EM(\pi_j-1)}
				
				p = A->s_inv_base_ii(p_el) - 1;
				// \pi_s^{\alpha^{-1}} = \pi_p
				if (f_v) {
					printf("\\pi_p = \\pi_%ld = p_el = %ld = \\pi_s^\\alpha^-1 ", p, p_el);
					}
				
				if (cl->orbit->s_ii(p)) {
					if (f_v) {
						printf("lies in \\pi_s orbit\n");
						}
					// if \pi_p is in the orbit of \pi_s
					
					// now we know:
					// \pi_s^{\alpha^{-1}} = \pi_p = \pi_s^cosetrep[p]
					// =>
					// \pi_s^{cosetrep[p]\alpha} = \pi_s
					// => 
					// cosetrep[p]\alpha fixes \pi_s 
					// and 
					// cosetrep[p] lies in G^(j).
					//
					// so, cosetrep[p]\alpha maps \pi_j-1 upon \pi_m.
#ifdef ANGST
					if (alpha_inv.s_ii(s_el) - 1 != p_el)
						return error("cycle_middle() alpha_inv.s_ii(s_el) - 1 != p_el");
#endif
				
					q = B->s_inv_base_ii(m_el) - 1;
					// we are going to prepare an edge 
					// \pi'_j -> \pi'_q in B
					if (q != j && cl->root->s_ii(q) == 1) {
						cl->root->m_ii(q, 0);
						if (f_v) {
							printf("new edge \\pi'j -> \\pi'_q, j = %ld, q = %ld\n", j, q);
							printf("root[%ld] = 0\n", q);
							}
						B->s_V_i(q)->m_i(j);
						cl->cosetrep->s_i(p)->mult(&alpha, B->s_KM_i(m_el));
						} // if
					
					} // if \pi_p is in orbit of \pi_s
				else {
					if (f_v) {
						printf("is not in \\pi_s orbit \n");
						}
					}

				} // if \pi_m is in delta
			} // next m
	
		if (f_v) {
			printf("root vector (in B): \n");
			cl->root->println();
			printf("calling update_orbit(%ld)\n", j);
			}
	
		cycle_middle_update_orbit(cl, j);
			// compute orbit / cosetreps for \pi_s under G^(j-1)
		
		if (f_v) {
			printf("\\pi_s orbit (in A): \n");
			cl->orbit->println();
			}
	
		} // next j


	
	cl->delta = cl->orbit;
	// cl->delta = cl->orbitlist;
	if (f_v) {
		printf("we found the following \\pi_s orbit:\n");
		printf("with \\pi_m = m_el = \\pi'_q \n");
		for (m = 0; m < n; m++) {
			if (cl->delta->s_ii(m) == 0)
				continue;
			m_el = A->s_base_ii(m) - 1;
			q = B->s_inv_base_ii(m_el) - 1;
			printf("m = %ld m_el = %ld \\mapsto q = %ld root[q] = %ld\n", 
				m, m_el, q, cl->root->s_ii(q));
			}
		for (q = 0; q < n; q++) {
			if (cl->root->s_ii(q) == 0)
				continue;
			q_el = B->s_base_ii(q) - 1;
			m = A->s_inv_base_ii(q_el) - 1;
			printf("q = %ld q_el = %ld \\mapsto m = %ld orbit[m] = %ld\n", 
				q, q_el, m, cl->orbit->s_ii(m));
			}
		}
	
	
	for (m = 0; m < n; m++) {
		if (cl->delta->s_ii(m) == 0)
			continue;
		m_el = A->s_base_ii(m) - 1;
		q = B->s_inv_base_ii(m_el) - 1;
		if (f_v) {
			printf("\\pi_s orbit element m = %ld m_el = %ld q = %ld root[q] = %ld\n", 
				m, m_el, q, cl->root->s_ii(q));
			}
		if (q != cl->r && cl->root->s_ii(q) == 1) {
			cl->root->m_ii(q, 0);
			B->s_V_i(q)->m_i(cl->r);
			q_el = B->s_base_ii(q) - 1;
			cl->cosetrep->s_i(m)->swap(B->s_KM_i(q_el));
			}
		} // next m

	nb_bm = cl->bottom_merke->s_li();
	for (i = 0; i < nb_bm; i++) {
		j = cl->bottom_merke->s_ii(i);
		j_el = A->s_base_ii(j) - 1;
		j_el_new = B->s_base_ii(j) - 1;
		A->s_KM_i(j_el)->swap(B->s_KM_i(j_el_new));
		} // next i

	if (f_v) {
		printf("end of cycle_middle()\n");
		printf("root vector (in B): \n");
		cl->root->println();
		}
	return OK;
}

static INT cycle_middle_init_data(CYCLE_LOCAL *cl)
{
	PERMUTATION_OP perm, perm1;
	LABRA_OP A = cl->A;
	// LABRA_OP B = cl->B;
	INT i, i_el, j, n;
	
	n = cl->n;
	A->child_vectors(cl->cv);
	// cv[i][j] = 1 iff \pi_j is a child of \pi_i

	// determine the orbit of \pi_s under G^(s):
	cl->cv->s_i(cl->s)->copy(cl->orbit);
	// orbit[i] = 1 iff \pi_i is a child of \pi_s
	
	cl->cosetrep->m_il(n);
	// if orbit[i] = 1 cosetrep[i] should be an element in G^(s) 
	// which maps \pi_s upon \pi_i
	
	perm = (PERMUTATION_OP) cl->cosetrep->s_i(cl->s);
	perm->m_il(n);
	perm->one();
	// cosetrep[s] := id
	
	for (i = cl->s + 1; i < cl->n; i++) {
		i_el = A->s_base_ii(i) - 1;
		if (cl->orbit->s_ii(i)) {
			j = A->s_V_ii(i);
				// cosetrep[j] already computed !
			if (j != i) {
				perm = (PERMUTATION_OP) cl->cosetrep->s_i(j);
				perm1 = (PERMUTATION_OP) cl->cosetrep->s_i(i);
				perm->mult(A->s_KM_i(i_el), perm1);
					// post-multiply KM(i_el) !
				}
			}
		} // next i
	// cl->cv->s_i(cl->s)->copy(cl->orbitlist);
	return OK;
}

INT cycle_middle_update_orbit(CYCLE_LOCAL *cl, INT j)
// extends orbit and cosetrep
// according to the (j-1)-th stabilizer
// afterwards: orbit and cosetrep describe 
// the orbit of \pi_s under G^(j-1)
{
	VECTOR_OB gens, Q;
	PERMUTATION_OB per, per1;
	LABRA_OP A = cl->A;
	// LABRA_OP B = cl->B;
	INT n, p, q, p_el, q_el, nb_Q, i;
	INT idx, idx_el, jj;

	n = cl->n;
	Q.m_il(n);
	nb_Q = 0;
	for (i = 0; i < n; i++) {
		if (cl->orbit->s_ii(i)) {
			Q.m_ii(nb_Q, i);
			nb_Q++;
			}
		}
#ifdef DEBUG_CYCLE
	printf("calling all_stab_generators(%ld)\n", j - 2);
#endif
	A->all_stab_generators(j - 2, &gens);
	while (nb_Q) {
		p = Q.s_ii(0);
		for (i = 1; i < nb_Q; i++) {
			Q.m_ii(i - 1, Q.s_ii(i));
			}
		nb_Q--;
		p_el = A->s_base_ii(p) - 1;
		for (jj = 0; jj < gens.s_li(); jj++) {
			idx = gens.s_ii(jj);	
			idx_el = A->s_base_ii(idx) - 1;
			A->the_son(idx, &per);
#if 0
#ifdef ANGST
			if (per.s_ii(j - 1) - 1 != j - 1) {
				printf("j = %ld idx = %ld base[idx] = %ld\n", j, idx, idx_el);
				printf("per = ");
				per.println();
				return error("cycle_middle_update_orbit() "
					"per.s_ii(j - 1) - 1 != j - 1");
				}
#endif
#endif
			q_el = per.s_ii(p_el) - 1;
			q = A->s_inv_base_ii(q_el) - 1;
			if (cl->orbit->s_ii(q) == 0) {
				
#ifdef DEBUG_CYCLE
				printf("\\pi_s orbit extended by q = %ld\n", q);
				fflush(stdout);
#endif
				cl->orbit->m_ii(q, 1);
				Q.m_ii(nb_Q, q);
				nb_Q++;
				cl->cosetrep->s_i(p)->mult(&per, cl->cosetrep->s_i(q));
#ifdef ANGST
				{
				INT s_el = A->s_base_ii(cl->s) - 1;
				PERMUTATION_OP per1;

				per1 = (PERMUTATION_OP) cl->cosetrep->s_i(q);
				if (per1->s_ii(s_el) - 1 != q_el)
					return error("cycle_middle_update_orbit() per1->s_ii(s_el) - 1) != q_el");
				}
#endif
				}
			} // next jj
		} // while
	return OK;
}

#if TEXDOCU
INT LABRA_OB::cycle_id()
#endif
{
	INT i, i_el, n;

	n = s_degree_i();
	for (i = 0; i < n; i++) {
		i_el = s_base_ii(i) - 1;
		if (i_el != i) {
			cycle(i, s_inv_base_ii(i) - 1, FALSE, FALSE);
			}
		}
	return OK;
}

#endif /* LABRA_TRUE */

