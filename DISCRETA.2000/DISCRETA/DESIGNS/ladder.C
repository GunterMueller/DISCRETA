/* ladder.C */

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

#undef FREE_STAB
#define MAX_STEP 64

#if TEXDOCU
INT DCY_OB::field_name(INT i, INT j, BYTE *str)
#endif
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "step"; break;
	case 1: s = "type"; break;
	case 2: s = "k"; break;
	case 3: s = "fDown"; break;
	case 4: s = "reduce_generators_mode"; break;
	case 5: s = "T"; break;
	case 6: s = "TDidx"; break;
	case 7: s = "TDfusel"; break;
	case 8: s = "D"; break;
	case 9: s = "Ad"; break;
	case 10: s = "omega"; break;
	case 11: s = "oiti"; break;
	default:
		return error("DCY::field_name()|"
		"i out of range");
	}
	strcpy(str, s);
	return OK;
}

#if TEXDOCU
INT DCY_OB::initialize_Young(VECTOR_OP T, INT n, INT i, INT type, void *data)
#else
This function initializes the DCY object for a transversal between Young subgroups 
of $S_n$. 
type is set to 0 to indicate a transversal between Young subgroups.
$i$ becomes the step number of the ladder. 
$k$ is computed and indicated the moved point by the transversals. 
$k$ is 0 for $i=0$ and $i=1$, 
$i/2$ for even $i$ and $(i+1)/2$ otherwise.
fDown indicates a downstep. it is true for $i=0$ and even $i$, 
false for $i=1$ and odd $i$.
The vectors D and Ad are allocated with length 0. 
In case of a downstep, we use the transversal $(1, 2, \ldots, n-k)^j$ 
for $j=0,\ldots,n-k-1$. 
For an upstep, we use $(n-k+1,\ldots,n)^j$ for $j=0,\ldots,k-1$.
The transversals are initialized only in the case $i > 0$ ($i=o$ means that there is 
nothing to compute).
#endif
{
	PERMUTATION_OB p;
	INT erg = OK, k, fDown;

	erg += m_il(12);
	c_obj_k(DCY);
	s_step()->m_i(i);
	if (i == 0 || i == 1)
		k = 0;
	else if (EVEN(i))
		k = i >> 1;
	else
		k = (i + 1) >> 1;
	s_k()->m_i(k);
	s_type()->m_i(0); /* type = 0: young */
	if (i == 0)
		fDown = FALSE;
	else if (i == 1)
		fDown = TRUE;
	else if (EVEN(i))
		fDown = TRUE;
	else
		fDown = FALSE;
	s_fDown()->m_i(fDown);
	s_reduce_generators_mode()->m_i(1);
	s_D()->m_il(0);
	s_Ad()->m_il(0);
	p.m_il(n);
	p.one();
	if (fDown) {
		/* (1 2 .. n - k)^j for j = 0 to n - k - 1 */
		erg += p.AddNCycle(1 /* first */, n - k /* len */);
		/* printf("p = ");
		p.println(); fflush(stdout); */
		erg += vec_multiples(T, &p);
		}
	else if (i > 0) {
		/* (n - k + 1 ... n)^j for j = 0 to k - 1 */
		erg += p.AddNCycle(n - k + 1 /* first */, k /* len */);
		/* printf("p = ");
		p.println(); fflush(stdout); */
		erg += vec_multiples(T, &p);
		}
	else
		T->m_il(0);
	/* if (i > 0)
		T->println(); */
	return erg;
}

#if TEXDOCU
INT DCY_OB::initialize_Wreath(VECTOR_OP T, INT n, INT i, INT type, void *data)
#else
This initializes transversals in wreath product $S_2 \wr S_n$ according to Weinrich.
type is set to 1 to indicate a wreath product type transversal.
#endif
{
	PERMUTATION_OB p;
	INT erg = OK, k, fDown;

	erg += m_il(12);
	c_obj_k(DCY);
	s_step()->m_i(i);
	k = i >> 2;
	s_k()->m_i(k);
	s_type()->m_i(1); /* type = 1: wreath */
	if ((i % 4) == 0)
		fDown = FALSE;
	else if ((i % 4) == 1)
		fDown = TRUE;
	else if ((i % 4) == 2)
		fDown = TRUE;
	else if ((i % 4) == 3)
		fDown = FALSE;
	s_fDown()->m_i(fDown);
	s_reduce_generators_mode()->m_i(1);
	
	s_D()->m_il(0);
	s_Ad()->m_il(0);
	p.m_il(n);
	p.one();
	if (i > 0 && (i % 4) == 1) {
		/* (1 2 ... (n - 2 * k))^j 
		 * for j = 0 to (n - 2 * k) */
		erg += p.AddNCycle(1 /* first */, n - 2 * k /* len */);
		erg += vec_multiples(T, &p);
		}
	else if (i > 0 && (i % 4) == 2) {
		/* (1 2 ... (n - 2 * k - 1))^j 
		 * for j = 0 to (n - 2 * k - 1) */
		erg += p.AddNCycle(1 /* first */, n - 2 * k - 1 /* len */);
		erg += vec_multiples(T, &p);
		}
	else if (i > 0 && (i % 4) == 3) {
		/* ((n - 2 * k - 1) (n - 2 * k))^j 
		 * for j = 0 to 1 */
		erg += p.AddNCycle(n - 2 * k - 1 /* first */, 2 /* len */);
		erg += vec_multiples(T, &p);
		}
	else if (i > 0 && (i % 4) == 0) {
		/* (((n - 2 * k + 1) (n - 2 * k + 3) (n - 1)) 
			((n - 2 * k + 2) (n - 2 * k + 4) (n)))^j 
		 * for j = 0 to k */
		erg += p.AddNCycleOffset(n - 2 * k + 1, k /* len */, 2 /* offset */);
		erg += p.AddNCycleOffset(n - 2 * k + 2, k /* len */, 2 /* offset */);
		erg += vec_multiples(T, &p);
		}
	else
		T->m_il(0);
	/* if (i > 0) {
		printf("%ld:\n", i);
		T->println();
		} */
	return erg;
}

#if TEXDOCU
INT DCY_OB::initialize_arbitrary(VECTOR_OP T, INT n, INT i, 
	INT omega, VECTOR_OP oiti, INT type, void *data)
#else
type is set to 2 to indicate an arbitrary transversal type (derived from 
a Young transversal).
$k$ and fDown are the same as in the Young transversal.
omega and oiti are copied into the DCY object.
#endif
{
	INT erg = OK, k, fDown;

	erg += m_il(12);
	c_obj_k(DCY);
	s_step()->m_i(i);
	if (i == 0 || i == 1)
		k = 0;
	else if (EVEN(i))
		k = i >> 1;
	else
		k = (i + 1) >> 1;
	s_k()->m_i(k);
	s_type()->m_i(2); /* type = 2: arbitrary (using omega and oiti) */
	if (i == 0)
		fDown = FALSE;
	else if (i == 1)
		fDown = TRUE;
	else if (EVEN(i))
		fDown = TRUE;
	else
		fDown = FALSE;
	s_fDown()->m_i(fDown);
	s_reduce_generators_mode()->m_i(1);
	if (i == 0) {
		s_omega()->m_i(-1);
		}
	else {
		s_omega()->m_i(omega);
		oiti->copy(s_oiti());
		T->copy(s_T());
		}
	s_D()->m_il(0);
	s_Ad()->m_il(0);
	return erg;
}

#if TEXDOCU
INT dc_young_downstep(SYM_OP m, SYM_OP dv, INT n, INT k, INT *idx, 
	INT type, void *data)
#else
This function is used during the downstep of the ladder algorithm.
The image $a= dv(m(n-k))$ is evaluated ($n-k-1$ if counting elements from 0)
and the index of the corresponding 
transversal element in $T$ is returned.
We assume that we are working with the transversal 
$(1,2,\ldots,n-k)^j$ for $j=0,1,\ldots,n-k-1$ as from initialize\_Young.
#endif
{
	INT i0, i1;
	
	i0 = do_image_of(m, n - k - 1, type, data);
	i1 = do_image_of(dv, i0, type, data);
	if (i1 >= n - k) {
		printf("dc_young_downstep()|i1 >= n - k\n");
		return ERROR;
		}
	if (i1 == n - k - 1)
		*idx = 0; /* id f"ur j = 0 */
	else
		*idx = i1 + 1; /* j = i1 + 1 */
	return OK;
}

#if TEXDOCU
INT dc_in_test_young(SYM_OP pi, INT n, INT k, INT fDown, 
	INT type, void *data)
#endif
/* Rueckgabe TRUE/FALSE */
{
	INT j, o;
	
	if (fDown) {
		o = do_image_of(pi, n - k - 1, type, data);
		if (o != n - k - 1)
			return FALSE;
		}
	for (j = 0; j <  k; j++) {
		o = do_image_of(pi, n - j - 1, type, data);
		if (o < n - k)
			return FALSE;
		}
	return TRUE;
}

#if TEXDOCU
INT DCY_OB::do_dc_young_downstep(SYM_OP m, SYM_OP dv, INT n, INT *idx, 
	INT type, void *data)
#endif
{
	INT k1, i, ti;

	if (s_type_i() == 0)
		dc_young_downstep(m, dv, n, s_k_i(), idx, type, data);			
	else if (s_type_i() == 1) {
		if (s_step_i() % 4 == 1)
			k1 = 2 * s_k_i();
		else if (s_step_i() % 4 == 2)
			k1 = 2 * s_k_i() + 1;
		else
			return error("do_dc_young_downstep(): step not cong 1 or 2 mod 4");			
		dc_young_downstep(m, dv, n, k1, idx, type, data);
		if (*idx >= s_T()->s_li())
			return error("do_dc_downstep(): *idx >= s_T()->s_li()");	
		}			
	else if (s_type_i() == 2) {
		i = s_omega_i();
		i = do_image_of(m, i, type, data);
		i = do_image_of(dv, i, type, data);
		if (i < 0 || i >= s_oiti()->s_li())
			return error("do_dc_young_downstep() i out of range");
		ti = s_oiti_ii(i);
		if (ti < 0 || ti >= s_T()->s_li())
			return error("do_dc_young_downstep() ti out of range");
		*idx = ti;
		}
	else
		return error("do_dc_young_downstep(): unknown type");
	return OK;			
}

#if TEXDOCU
INT dc_do_step(DCY_OP dc0, 
	INT step, SYM_OP id, INT n, INT f_verbose, 
	INT type, void *data)
#endif
{
	INT fDown, erg = OK;
	
	if (step <= 0)
		return error("dc_do_step()|step <= 0");
	fDown = dc0[step].s_fDown_i();
	if (fDown)
		erg += dc0[step].do_downstep(dc0 + (step - 1), id, n, step, f_verbose, type, data);
	else
		erg += dc0[step].do_upstep(dc0, id, n, step, f_verbose, type, data);
	if (erg)
		return error("error in dc_do_step()");
	return erg;
}

#define STEP_VERBOSE 4

#if TEXDOCU
INT DCY_OB::do_downstep(DCY_OP dc_last, SYM_OP id, 
	INT deg, INT step, INT f_verbose, INT type, void *data)
#endif
{
	VECTOR_OB TD, Ad;
	SYM_OB m1, n, u, dv, m;
	INT erg = OK, i, j, k, l, ii;
	INT ol, cur, idx, f_found, index, nb_D, Glen;
	INT *queue = NIL;
	INT nb_queue = 0;
	SYM_OB group_order;
	VECTOR_OP G;
	
	if (f_verbose) {
		printf("step %ld (down):", step);
		fflush(stdout);
		}
	index = s_T()->s_li();
	if (f_verbose) {
		printf(" index = %ld", index);
		fflush(stdout);
		}
	nb_D = dc_last->s_D()->s_li();
	if (f_verbose) {
		printf(" nb_D = %ld\n", nb_D);
		fflush(stdout);
		}
	queue = (INT *) my_malloc(index * sizeof(INT), "DCY::do_downstep");
	if (queue == NIL)
		return error("DCY::do_downstep: no memory for queue");
	erg += s_TDidx()->m_ilih(nb_D, index);
	erg += s_TDfusel()->m_ilih(nb_D, index);
	if (erg != OK)
		return error("DC::do_downstep: no memory for TDidx / TDfusel");
	for (i = 0; i < index; i++) {
		for (j = 0; j < nb_D; j++) {
			s_TDidx()->m_iji(i, j, -1);
			}
		}
	erg += TD.m_il(index);
	if (erg != OK)
		return error("DC::do_downstep: no memory for TD");

	s_D()->m_il(0);
	s_Ad()->m_il(0);
	for (j = 0; j < nb_D; j++) {
		if (f_verbose) {
			printf("%ld ", j);
			if ((j + 1) % 10 == 0)
				printf("\n");
			fflush(stdout);
			}
		G = dc_last->s_Ad_i(j);
		Glen = G->s_li();
		do_invers(dc_last->s_D_i(j), &dv, type, data);
		/* dc_last->s_D_i(j)->invers(&dv); */
		vec_translate(s_T(), dc_last->s_D_i(j), &TD, type, data);
		for (i = 0; i < index; i++) {
			if (s_TDidx_iji(i, j) != -1)
				continue;
			k = s_D()->s_li();
			s_TDidx()->m_iji(i, j, k);
			do_copy(id, s_TDfusel_ij(i, j), type, data);
			
			nb_queue = 0;
			queue[nb_queue++] = i;
			if (f_verbose && step <= STEP_VERBOSE) {
				printf("%ld/%ld: ", i, j);
				/* printf("\n");
				TD->s_i(i)->println(); */
				fflush(stdout);
				}
			ol = 0;
			while (nb_queue) {
				cur = queue[0];
				for (l = 1; l < nb_queue; l++)
					queue[l - 1] = queue[l];
				nb_queue--;
				for (l = 0; l < Glen; l++) {
					do_mult(TD.s_i(cur), G->s_i(l), &m, type, data);
					/* G->s_i(l)->mult(TD.s_i(cur), &m); */
					do_dc_young_downstep(&m, &dv, deg, &idx, type, data);			
					if (s_TDidx_iji(idx, j) == -1) {
						if (f_verbose && step <= STEP_VERBOSE) {
							printf("%ld ", idx);
							fflush(stdout);
							ol++;
							if (ol % 20 == 0) {
								printf("\n");
								fflush(stdout);
								}
							}
						do_mult(s_TDfusel_ij(cur, j), G->s_i(l), s_TDfusel_ij(idx, j), type, data);
						s_TDidx()->m_iji(idx, j, k);
						queue[nb_queue++] = idx;
						}
					} /* for l */
				} /* while */
			/* a new double coset is found:
			 * add it to D: */
			s_D()->inc();
			do_copy(TD.s_i(i), s_D_i(s_D()->s_li() - 1), type, data);
			
			/* the stabilizer: */
			Ad.m_il(1);
			do_copy(id, Ad.s_i(0), type, data);
			
			for (cur = 0; cur < index; cur++) {
				if (s_TDidx_iji(cur, j) != k)
					continue;
				for (l = 0; l < Glen; l++) {
					do_mult(TD.s_i(cur), G->s_i(l), &m, type, data);
					do_dc_young_downstep(&m, &dv, deg, &idx, type, data);
					
					/* form the schreier generator: */
					do_mult(s_TDfusel()->s_ij(cur, j), G->s_i(l), &m1, type, data);
					do_invers(s_TDfusel()->s_ij(idx, j), &n, type, data);
					do_mult(&m1, &n, &u, type, data);
					
					if (!do_einsp(&u, type, data)) {
						/* empty stabilizer generators 
						 * causes problems in VS_search */
						Ad.do_search(Ad.s_li(), TRUE, &u, &idx, &f_found, type, data);
						if (!f_found) {
							Ad.inc();
							for (ii = Ad.s_li() - 1; ii > idx; ii--)
								Ad.s_i(ii)->swap(Ad.s_i(ii - 1));
							do_copy(&u, Ad.s_i(idx), type, data);
							}
						}
					} /* for l */
				} /* for cur */
			if (f_verbose && step <= STEP_VERBOSE) {
				printf(" - # %ld", Ad.s_li());
				fflush(stdout);
				}
			if (s_reduce_generators_mode_i() == 2) {
				LABRA_OB lab;
				
				reduce_generators_labra(&Ad, &group_order, 
					FALSE /* f_verbose */, &lab);
				}
			else if (s_reduce_generators_mode_i() == 1) {
				reduce_generators(&Ad, &group_order, 
					FALSE /* f_verbose */);
				}
			if (s_reduce_generators_mode_i() != 0) {
				if (f_verbose && step <= STEP_VERBOSE) {
					printf(" reduced to %ld go ", 
					Ad.s_li());
					group_order.println();
					}
				}
			if (f_verbose && step <= STEP_VERBOSE) {
				printf("\n");
				fflush(stdout);
				}
			s_Ad()->inc();
			Ad.swap(s_Ad_i(s_Ad()->s_li() - 1));
			} /* for i */
#ifdef FREE_STAB
		G->freeself();
#endif
		} /* for j */
#ifdef FREE_STAB
	dc_last->s_Ad()->freeself();
#endif
	if (f_verbose) {
		printf("\n");
		fflush(stdout);
		}
	printf("step %ld (down); found %ld double cosets\n", step, s_D()->s_li());
	fflush(stdout);
	if (queue) {
		my_free(queue);
		queue = NIL;
		}
	return erg;
}

#if TEXDOCU
INT DCY_OB::do_upstep(DCY_OP dc0, SYM_OP id, 
	INT deg, INT step, INT f_verbose, INT type, void *data)
#endif
{
	INT erg = OK, i, j, k, ol, kk;
	INT dc_idx, idx, f_found, index, nb_D, ii;
	VECTOR_OB TD, Ad;
	SYM_OB d, fusel;

	index = s_T()->s_li();
	nb_D = dc0[step - 1].s_D()->s_li();
	if (f_verbose) {
		printf("step %ld (up): index = %ld nb_D = %ld\n", step, index, nb_D);
		fflush(stdout);
		}
	erg += s_TDidx()->m_ilih(nb_D, 1);
	erg += s_TDfusel()->m_ilih(nb_D, 1);
	for (j = 0; j < nb_D; j++)
		s_TDidx()->m_iji(0, j, -1);
	if (erg != OK)
		return error("DCY::do_upstep: no memory for TDidx / TDfusel");
	erg += TD.m_il(index);
	if (erg != OK)
		return error("DCY::do_upstep: no memory for TD");
	s_D()->m_il(0);
	s_Ad()->m_il(0);
	kk = 0;
	for (i = 0; i < nb_D; i++) {
		if (s_TDidx_iji(0, i) != -1)
			continue;
		vec_translate(s_T(), dc0[step - 1].s_D_i(i), &TD, type, data);
		k = s_D()->s_li();
		if (f_verbose) {
			printf("%ld ", i);
			if ((kk + 1) % 10 == 0)
				printf("\n");
			kk++;
			fflush(stdout);
			}
		if (f_verbose && step <= STEP_VERBOSE) {
			printf("%ld: ", i);
			fflush(stdout);
			}
		ol = 0;
		do_copy(id, 
		s_TDfusel_ij(0, i), type, data);
		s_TDidx()->m_iji(0, i, k);
		
		/* add the double coset: */
		erg += s_D()->inc();
		do_copy(dc0[step - 1].s_D_i(i), &d, type, data);
		do_copy(&d, s_D_i(s_D()->s_li() - 1), type, data);
		
		/* copy the old stabilizer: */
#if 0
		dc0[step -1].s_Ad_i(i)->copy(&Ad); /* !!! */
#endif
		{
		VECTOR_OP Ad_old;
		INT ii, l;
		
		Ad_old = dc0[step - 1].s_Ad_i(i);
		l = Ad_old->s_li();
		Ad.m_il(l);
		for (ii = 0; ii < l; ii++) {
			do_copy(Ad_old->s_i(ii), 
			Ad.s_i(ii), type, data);
			}
		}
		
		for (j = 0; j < index; j++) {
			erg += dc_trace_coset(dc0, TD.s_i(j), step - 1, 
				&dc_idx, &fusel, type, data);
			if (erg)
				return error("DCY::do_upstep: error in DCY::trace_coset");
			if (dc_idx == i) {
				if (!do_einsp(&fusel, type, data)) {
					Ad.do_search(Ad.s_li(), TRUE, &fusel, &idx, &f_found, type, data);
					if (!f_found) {
						Ad.inc();
						for (ii = Ad.s_li() - 1; ii > idx; ii--)
							Ad.s_i(ii)->swap(Ad.s_i(ii - 1));
						do_copy(&fusel, Ad.s_i(idx), type, data);
						}
					}
				}
			else { /* merge two double cosets */
				if (f_verbose && step <= STEP_VERBOSE) {
					printf("%ld ", dc_idx);
					fflush(stdout);
					ol++;
					if (ol % 20 == 0) {
						printf("\n");
						fflush(stdout);
						}
					}
				s_TDidx()->m_iji(0, dc_idx, k);
				do_invers(&fusel, s_TDfusel_ij(0, dc_idx), type, data);
				}
			} /* for j */
		erg += s_Ad()->inc();
		Ad.swap(s_Ad_i(s_Ad()->s_li() - 1));
		
		if (f_verbose && step <= STEP_VERBOSE) {
			printf("\n");
			fflush(stdout);
			}
		} /* for i */
	if (f_verbose) {
		printf("\n");
		fflush(stdout);
		}
#ifdef FREE_STAB
	dc0[step -1].s_Ad()->freeself();
#endif
	printf("step %ld (up); found %ld double cosets\n", step, s_D()->s_li());
	fflush(stdout);
	return erg;
}

#undef TRACE_COSET_DEBUG

#if TEXDOCU
INT dc_trace_coset(DCY_OP dc0, SYM_OP g, INT i, INT *dc_idx, SYM_OP fusel, 
	INT type, void *data)
#endif
{
	SYM_OB d, dv, av;
	SYM_OB tmp, tmp1, fusel1;
	DCY_OP L;
	INT j, d_idx_jm1, d_idx;
	INT idx, deg;
	INT erg = OK;
	
	L = dc0;
	d_idx = 0; /* L[0].D[0] = id */
	do_copy(L->s_D_i(d_idx), &d, type, data);
	do_copy(&d, fusel, type, data);
	deg = ((PERMUTATION_OP) &d)->s_li();
	do_one(fusel, type, data);
	/* fusel->one(); */
	for (j = 1; j <= i; j++) {
		L = dc0 + j;
		d_idx_jm1 = d_idx;
		/* now Ljm1 d fusel = Ljm1 g */
		if (L->s_fDown_i()) {
			do_invers(fusel, &av, type, data);
			do_invers(&d, &dv, type, data);
			do_mult(g, &av, &tmp, type, data);
			
			/* now tmp1 = g ajm1_inv djm1_inv, 
			 * and is element of Ljm1 */
			/* find it's coset of Lj in Ljm1: */
			L->do_dc_young_downstep(&tmp, &dv, deg, &idx, type, data);			
			d_idx = L->s_TDidx_iji(idx, d_idx_jm1);
			if (d_idx < 0)
				return error("trace_coset (down): d_idx < 0");
			do_mult(L->s_TDfusel_ij(idx, d_idx_jm1), 
			fusel, &fusel1, type, data);
			do_copy(&fusel1, fusel, type, data);
			do_copy(L->s_D_i(d_idx), &d, type, data);
			/* now Lj d fusel = Lj g */
			
#ifdef TRACE_COSET_DEBUG
			/* we check it: */
			erg += fusel->invers(a_inv);
			erg += d->invers(d_inv);
			erg += a_inv->mult(g, tmp);
			erg += d_inv->mult(tmp, tmp1);
			if (!in_test_young(tmp1, L->s_n_i(), L->s_k_i(), TRUE))
				return error("DCC::trace_coset()|!in_test_young fDown TRUE");
#endif
			}
		else {
			/* g fusel^-1 is element of Ljm1 d 
			 * this coset falls into 
			 * the coset TDidx[d_idx_jm1]
			 * a new fusing element a 
			 * with Lj d a contains Ljm1 d 
			 * is TDfusel[d_idx_jm1] */
			d_idx = L->s_TDidx_iji(0, d_idx_jm1);
			if (d_idx < 0)
				return error("trace_coset (up): d_idx < 0");
			do_mult(L->s_TDfusel_ij(0, d_idx_jm1), fusel, &fusel1, type, data);
			do_copy(&fusel1, fusel, type, data);
			do_copy(L->s_D_i(d_idx), &d, type, data);
			/* now Lj d fusel = Lj g */

#ifdef TRACE_COSET_DEBUG
			/* we check it: */
			erg += fusel->invers(a_inv);
			erg += d->invers(d_inv);
			erg += a_inv->mult(g, tmp);
			erg += d_inv->mult(tmp, tmp1);
			if (!in_test_young(tmp1, L->s_n_i(), L->s_k_i(), FALSE))
				return error("DCC::trace_coset()|!in_test_young fDown FALSE");
#endif
			}
		}
	*dc_idx = d_idx;
	return erg;
}

#if TEXDOCU
INT vec_translate(VECTOR_OP V, SYM_OP x, VECTOR_OP Vx, 
	INT type, void *data)
#endif
{
	INT erg = OK, i;
	
	for (i = 0; i < V->s_li(); i++)
		do_mult(V->s_i(i), x, Vx->s_i(i), type, data);
	return erg;
}

#if TEXDOCU
INT vec_multiples(VECTOR_OP V, SYM_OP p)
#endif
/* V = { p^0, p^1, ... p^k-1 }, wenn p^k = id */
{
	SYM_OB pi, pip1;

	p->copy(&pi);
	pi.one();
	V->m_il(1);
	pi.swap(V->s_i(0));
	p->copy(&pi);
	while (!pi.einsp()) {
		V->inc();
		pi.copy(V->s_i(V->s_li() - 1));
		p->mult(&pi, &pip1);
		pip1.swap(&pi);
		}
	return OK;
}

#endif /* LADDER_TRUE */
