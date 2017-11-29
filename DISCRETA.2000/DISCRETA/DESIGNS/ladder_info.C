/* ladder_info.C */

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
#include <DISCRETA/fga.h> /* for group types */


#if TEXDOCU
#else
The ladder game (snakes and ladders) has been invented by 
Schmalz~\cite{Schmalz90},~\cite{Schmalz92b},~\cite{Schmalz93} 
in order to compute double cosets. 
#endif


#if TEXDOCU
LADDER_INFO *init_ladder_info(
	VECTOR_OP generators, SYM_OP go, INT deg, 
	BYTE *g_label, BYTE *g_label_tex, 
	INT t, INT k, INT lambda, 
	INT f_verbose, 
	BYTE *generators_fname)
#else
This function opens a LADDER\_INFO structure (allocates the memory) 
and initializes the fields: 
generators, deg, g\_label, g\_label\_tex, t, k, lambda, f\_verbose. 
The vector dc is allocated; it will be used later on to hold the ladder 
for the algorithm Leiterspiel. 
An object holding the identity group element id is allocated and initialized.
If generators\_fname is given, it is copied into the structure.
up\_to\_step is computed as $2\cdot k - 1$.
The dc vector is allocated to hold up\_to\_step $+ 1$ elements.
type is set to DO\_TYPE\_SYM and data is set to NIL.
#endif
{
	LADDER_INFO *li;

	li = (LADDER_INFO *) my_malloc(sizeof(LADDER_INFO), "init_ladder_info");
	if (li == NIL) {
		error("init_ladder_info: no memory");
		return NIL;
		}
	li->generators = generators;
	li->go = go;
	li->deg = deg;
	li->g_label = g_label;
	li->g_label_tex = g_label_tex;
	
	li->t = t;
	li->k = k;
	li->lambda = lambda;
	li->f_verbose = f_verbose;

	li->dc = (VECTOR_OP) callocobject("li->dc");
	li->id = callocobject("li->dc");
	generators->s_i(0)->copy(li->id);
	li->id->one();

	li->generators_fname[0] = 0;
	if (generators_fname)
		strcpy(li->generators_fname, generators_fname);
	
	li->up_to_step = 2 * li->k - 1;
	printf("up_to_step = %ld\n", li->up_to_step);
	fflush(stdout);
	li->dc->m_il(li->up_to_step + 1);
	li->type = DO_TYPE_SYM;
	li->data = NIL;
	return li;
}

#if TEXDOCU
void free_ladder_info(LADDER_INFO *li)
#else
This function deallocates the structure LADDER\_INFO and all its 
allocated objects (dc and id).
#endif
{
	if (li->dc)
		freeall(li->dc);
	if (li->id)
		freeall(li->id);
	my_free(li);
	printf("LADDER_INFO closed\n");
}

#if TEXDOCU
INT li_init_file_names(LADDER_INFO *li)
#else
Initializes txt\_out to hold the file name 
$KM\_\langle group \rangle\_t\langle t \rangle\_k \langle k \rangle.txt$. 
bin\_out holds the same name but with .bin in the end instead.
A file \lq km\_fname\rq is written which holds the file name in txt\_out 
of the KM-file to be created later on.
#endif
{
	sprintf(li->txt_out, "KM_");
	strcat(li->txt_out, li->g_label);
	sprintf(li->txt_out + strlen(li->txt_out), "_t%ld_k%ld", li->t, li->k);
	strcpy(li->str, li->txt_out);
	strcat(li->txt_out, ".txt");

	strcpy(li->bin_out, li->str);
	strcat(li->bin_out, ".bin");
	
	li->fp = fopen("km_fname", "w");
	fputs(li->txt_out, li->fp);
	fclose(li->fp);

	printf("%s\n", li->txt_out); fflush(stdout);
	return OK;
}

#if TEXDOCU
void li_message(LADDER_INFO *li)
#else
This function displays the DISCRETA banner and prints the values of 
g\_label, t, k, lambda, f\_verbose.
#endif
{
	printf(discreta_copyright_text);

	printf("LADDER_INFO opened:\n");
	printf("group: %s\n", li->g_label);
	printf("t = %ld k = %ld\n", li->t, li->k);
	printf("lambda = %ld "
		"f_verbose = %ld\n", 
		li->lambda, li->f_verbose);
}

#if TEXDOCU
INT li_init_transversals(LADDER_INFO *li)
#else
Prepares the transversals for the Leiterspiel algorithm. 
Calls the function DCY::initialize\_Young() for all steps of the ladder.
It also initializes the first (0-th) ladder step: 
There is only one double coset for $G \backslash G / A$, namely 
the identity (for instance). The stabilizer Stab$_A(G id)=A$ is 
initialized in the 0-th entry of the Ad vector.
#endif
{
	INT i;
	VECTOR_OP A = li->generators;
	DCY_OP dcy = (DCY_OP) li->dc->s_s();
	VECTOR_OB T;
	
	printf("li_init_transversals() up_to_step = %ld\n", li->up_to_step);
	fflush(stdout);
	do_copy(A->s_i(0), li->id, li->type, li->data);
	do_one(li->id, li->type, li->data);
	for (i = 0; i <= li->up_to_step; i++) {
		printf("init transversal i = %ld\n", i);
		fflush(stdout);
		dcy[i].initialize_Young(&T, 
			li->deg, i /* step */, 
			0 /* type */, NIL /* data */);
		T.swap(dcy[i].s_T());
		if (i <= 6 || TRUE)
			dcy[i].
			s_reduce_generators_mode()->m_i(2);
				/* 0: no reduction 
				 * 1: dimino 
				 * 2: labra */
		}
		
	dcy[0].s_D()->m_il(1);
	do_copy(li->id, dcy[0].s_D_i(0), li->type, li->data);
	
	dcy[0].s_Ad()->m_il(1);
	((SYM_OP) A)->swap(dcy[0].s_Ad()->s_i(0));
	/* dcy[0].s_D()->println();
	dcy[0].s_Ad()->println(); */
	printf("finished with li_init_transversals()\n");
	fflush(stdout);
	return OK;
}

#if TEXDOCU
INT li_leiterspiel(LADDER_INFO *li)
#else
This function calls dc\_do\_step() for all ladder steps 
$1,\ldots,$up\_to\_step.
Afterwards, the number of double cosets on each step is printed.
#endif
{
	INT i, nb_D;
	DCY_OP dcy = (DCY_OP) li->dc->s_s();

	for (i = 1; i <= li->up_to_step; i++) {
		printf("li_leiterspiel(): before step %ld\n", i);
		fflush(stdout);
		dc_do_step(&dcy[0], i, li->id, 
		li->deg /* deg */, 
		li->f_verbose, 
		li->type, li->data);
		}
	printf("nach li_leiterspiel()\n");
	for (i = 1; i <= li->up_to_step; i++) {
		nb_D = dcy[i].s_D()->s_li();
		printf("number of double cosets on step %ld: %ld\n", i, nb_D);
		}
	fflush(stdout);
	return OK;
}

#if TEXDOCU
INT li_Mtk(LADDER_INFO *li, MATRIX_OP M, INT t, INT k)
#else
This function calls dc\_Mtk() with the appropriate DCY pointer.
#endif
{
	DCY_OP dcy = (DCY_OP) li->dc->s_s();

	dc_Mtk(dcy, t, k, M);
	return OK;
}

#if TEXDOCU
INT li_stab_go(LADDER_INFO *li, VECTOR_OP V, INT t)
#else
This function computes the stabilizer orders i.e. the orders 
of all groups in the Ad-vector for $t$-orbit representatives.
The functions reduce\_generators\_labra() and/or reduce\_generators() 
are called and the result is the vector V.
#endif
{
	DCY_OP dcy = (DCY_OP) li->dc->s_s();
	DCY_OP dc2;
	INT l, i, t1;

	if (t == 0)
		t1 = 0;
	else if (t == 1)
		t1 = 1;
	else
		t1 = 2 * t - 1;
	/* stabilizer order t-sets: */
	printf("computing stabilizer order %ld sets\n", t);
	fflush(stdout);
	dc2 = dcy + t1;
	l = dc2->s_D()->s_li();
	V->m_il_n(l);
	for (i = 0; i < l; i++) {
		if (i <= 2 || TRUE) {
			LABRA_OB L;
			
			reduce_generators_labra(
				dc2->s_Ad_i(i), 
				V->s_i(i), 
				FALSE /* f_verbose */, &L);
			}
		else {
			reduce_generators(
				dc2->s_Ad_i(i), 
				V->s_i(i), 
				FALSE /* f_verbose */);
			}
		/* V->s_i(i)->println(); */
		}
	/* printf("\n"); */
	return OK;
}

#define BUFSIZE_ 1024

#if TEXDOCU
INT li_perm_rep(LADDER_INFO *li, INT t)
#else
This function is too specialized at the moment. 
It reads generators of a group from the file \lq normperm\rq 
which must normalize the group A of the ladder. 
In this case, normperm induces a permutation group on the 
double cosets in each step and here we are interested 
in its action on the $t$-orbit representatives 
(or better the $t$-orbits).
Here, normperm consists of exactly one element, call it $n$.
The permutation representation $r$ of $n$ on $t$-orbits is calculated.
Then a file \lq 14\rq is opened which holds solution vectors 
of designs with blocks of size $t$ (so $t$ should better be called $k$). 
These solution vectors are then mapped under $r$, 
i.e. $r$ is applied to each entry 1 in these vectors 
and the resulting vectors are printed on the screen.
This function can be used to determine isomorphic solutions of designs.
#endif
{
	VECTOR_OB P;
	PERMUTATION_OP p;
	PERMUTATION_OB q, fusel, r;
	DCY_OP dcy = (DCY_OP) li->dc->s_s();
	DCY_OP L;
	PERMUTATION_OP d;
	PERMUTATION_OB dv;
	VECTOR_OB R;
	INT l, i, j, t1, kk, ii;
	FILE *fp;
	char str[BUFSIZE_];
	char str1[BUFSIZE_];
	

	if (read_file_of_generators(&P, "normperm") == ERROR)
		return OK;
	/* perm_vec_gen_ik(&P, p_deg, p_per_kind); */
	p = (PERMUTATION_OP) P.s_i(0);
	
	if (t == 0)
		t1 = 0;
	else if (t == 1)
		t1 = 1;
	else
		t1 = 2 * t - 1;
	printf("computing representation on %ld - orbits\n", t);
	fflush(stdout);
	L = dcy + t1;
	l = L->s_D()->s_li();
	r.m_il(l);
	for (i = 0; i < l; i++) {
		d = (PERMUTATION_OP) L->s_D_i(i);
		d->mult(p, &q); 
		dc_trace_coset(
			dcy, &q, 
			t1, &ii /* dc_idx */, &fusel, 
			0 /* type */, NIL /* data */);
		r.m_ii(i, ii + 1);
		printf("%ld ", ii);
		fflush(stdout);
		}
	printf("\n");
	
	fp = fopen("14", "r");
	if (fp == NULL)
		return OK;
	printf("applying permutation to the 14 solutions...\n");
	for (i = 0; i < 14; i++) {
		fgets(str, BUFSIZE_, fp);
		for (j = 0; j < l; j++)
			str1[j] = '0';
		str1[l] = 0;
		for (j = 0; j < l; j++) {
			if (str[j] == '0')
				continue;
			kk = r.s_ii(j) - 1;
			str1[kk] = '1';
			}
		printf("%s\n", str1);
		}
	fclose(fp);
	return OK;
}

#if TEXDOCU
INT li_set_reps(LADDER_INFO *li, VECTOR_OP V, INT t)
#else
This function computes orbit-representatives of all $t$-orbits. 
A vector $V$ is created. The orbit representatives hold elements 
$1,\ldots,n$ when $n$ is the degree of the group (so we start 
with 1 here as we do it with the images of permutations).
#endif
{
	DCY_OP dcy = (DCY_OP) li->dc->s_s();
	DCY_OP L;
	PERMUTATION_OP d;
	PERMUTATION_OB dv;
	VECTOR_OB R;
	INT l, i, j, j0, j1, t1, kk;

	if (t == 0)
		t1 = 0;
	else if (t == 1)
		t1 = 1;
	else
		t1 = 2 * t - 1;
	printf("computing representatives %ld sets\n", t);
	fflush(stdout);
	L = dcy + t1;
	l = L->s_D()->s_li();
	V->m_il(l);
	for (i = 0; i < l; i++) {
		d = (PERMUTATION_OP) L->s_D_i(i);
		d->invers(&dv);
		R.m_il(t);
		j0 = li->deg - (t + 1) + 1;
		for (j = 0; j < t; j++) {
			j1 = j0 + j;
			kk = d->s_ii(j1) - 1;
			R.m_ii(j, kk + 1);
			/* + 1 to have the same numbering 
			 * of elements as in the permutations ! */
			}
		R.swap(V->s_i(i));

		/* V->s_i(i)->println(); */
		}
	/* printf("\n"); */
	return OK;
}

#if TEXDOCU
INT li_print_dc_info(LADDER_INFO *li, 
	VECTOR_OP stab_go, VECTOR_OP reps, 
	INT t, FILE *fp, INT f_verbose)
#else
#endif
{
	DCY_OP dcy = (DCY_OP) li->dc->s_s();
	DCY_OP L;
	SYM_OB length;
	VECTOR_OP R, V;
	INT i, j, t1, l, len;
	BYTE str[256];
	BYTE str1[256];

	if (t == 0)
		t1 = 0;
	else if (t == 1)
		t1 = 1;
	else
		t1 = 2 * t - 1;
	L = dcy + t1;
	l = L->s_D()->s_li();
	fprintf(fp, "dc info %ld-sets:\n", t);
	for (i = 0; i < l; i++) {
		if (f_verbose) {
			fprintf(fp, "%ld: stab_go = ", i);
			str[0] = 0;
			stab_go->s_i(i)->sprint(str);
			fprintf(fp, "%s ", str);
			}
		li->go->ganzdiv(stab_go->s_i(i), &length);
		R = (VECTOR_OP) reps->s_i(i);
		fprintf(fp, "{");
		len = R->s_li();
		for (j = 0; j < len; j++) {
			fprintf(fp, " %ld", R->s_ii(j));
			if (j < len - 1)
				fprintf(fp, ",");
			}
		fprintf(fp, " }_");
		str[0] = 0;
		str1[0] = 0;
		length.sprint(str);
		stab_go->s_i(i)->sprint(str1);
		fprintf(fp, "%s,%s\n", str, str1);

		if (f_verbose) {
			V = L->s_Ad_i(i);
			if (V->s_li()) {
				fprintf(fp, "stab = ");
				v_do_fprintln(V, fp, 
				FALSE /* f_numerated */, 
				1 /* n_on_a_row */, 
				0 /* type */, NIL /* data */ );
				}
			}
		fflush(fp);
		}
	return OK;
}


#if TEXDOCU
INT li_print_M_asc(LADDER_INFO *li, 
	MATRIX_OP M, INT t, INT k, INT f_k2, INT k2)
#else
Calls km\_print\_M\_asc().
#endif
{
	DCY_OP dcy = (DCY_OP) li->dc->s_s();
	VECTOR_OP G_gen;

	G_gen = dcy[0].s_Ad_i(0);
	km_print_M_asc(li->g_label, li->g_label_tex, li->txt_out, 
		G_gen, li->go, li->deg, 
		M, t, k, f_k2, k2);
	return OK;
}


#endif /* LADDER_TRUE */

