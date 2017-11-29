/* list.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef LISTTRUE

#include <DISCRETA/list.h>

INT test_list()
{
	INT i, j;
	LIST_OP list = (LIST_OP) callocobject("test_list");
	INTEGER_OP z;
	
	list->init(LIST);
	j = 7; /* hoffentlich primitiv mod 61 (?! nein !) */
	for (i = 0; i < 20; i++) {
		j = (j * j) % 61;
		z = (INTEGER_OP) callocobject("test_list");
		z->m_i(j);
		list->insert(z, NULL, list_cmp_func);
		z = NIL;
		}
	list->fprint(stdout);
	list->freeself();
	printf("nach freeself\n");
	freeall(list);
	return OK;
}

INT list_cmp_func(SYM_OP a, SYM_OP b)
{
	INT i;
	
	i = a->sym_comp(b);
	if (i != 0)
		return i;
	return 1;
}

static struct list * calloclist(void);

static struct list * calloclist()
{
	struct list *a = 
		(struct list *) my_malloc(sizeof(struct list), "calloclist()");
	if (a == NULL)
		error("kein Speicher in calloclist");

	return(a);
}

INT LIST_OB::m_sn(SYM_OP self, LIST_OP next)
/* before: m_sn_l(OP self, OP nx, OP a)*/
{
	SYM_OP s = NULL;
	LIST_OP n = NULL;
	
	if (self != NULL) { 
		s = callocobject("LIST");
		self->copy(s);
		}
	if (next != NULL) {
		n = (LIST_OP) callocobject("LIST");
		next->copy(n);
		}
	return b_sn(s, n);
}

INT LIST_OB::b_sn(SYM_OP self, LIST_OP next)
/* build_self next_list */
/* before: b_sn_l(OP self, OP nx, OP a) */
{
	OBJECTSELF d;

	d.ob_list = calloclist();
	b_ks(LIST, d); 
	c_s(self); 
	c_n(next);
	return OK;
}

INT LIST_OB::emptyp()
/* true falls es sich um eine leere liste handelt, d.h. self == NULL */
{ 
	if (!listp()) 
		return FALSE;
	if (s_s() != NULL) 
		return FALSE;
	return TRUE;
}

INT LIST_OB::lastp()
{
	return (s_n() == NULL);
}

INT LIST_OB::length()
/* before: length_list(OP list, OP erg) */
{
	INT l = 0;
	LIST_OP zeiger = this;
	
	if (emptyp()) 
		return 0;
	while (zeiger != NULL) {
		l++;
		zeiger = zeiger->s_n();
		}
	return l;
}

INT LIST_OB::fprint(FILE *f)
/* ausgabe eines list-objects
 * ausgabe bis einschliesslich next == NULL */
{
	LIST_OP zeiger = this;
	OBJECTSELF d;
	
	d = s_obj_s();
	if (d.ob_list == NULL)
		return error("LIST::fprint: s_obj_s == NULL");

	if ((s_s() == NULL) && (s_s() == NULL))
		/* ein list object wird initialisiert mit b_sn_l(NULL,NULL,obj) */
		fprintf(f, "empty list");
	else
		while (zeiger != NULL) {
			zeiger->s_s()->fprint(f);
			fprintf(f, "  ");
			zeiger = zeiger->s_n();
		}
	return(OK);
}

INT LIST_OB::sprint(BYTE *str)
/* appends to str. writes to maximal strlength of 200. */
{
	BYTE str1[256];
	LIST_OP zeiger = this, zeiger1;
	OBJECTSELF d;
	
	d = s_obj_s();
	if (d.ob_list == NULL)
		return error("LIST::sprint: s_obj_s == NULL");
	if ((s_s() == NULL) && (s_n() == NULL)) {
		strcat(str, "empty list");
		return OK;
		}
	strcat(str, "(");
	while (zeiger != NULL) {
		str1[0] = 0;
		if (zeiger->s_s())
			zeiger->s_s()->sprint(str1);
		else
			strcpy(str1, "***");
		if (strlen(str) + strlen(str1) < 200)
			strcat(str, str1);
		else
			return OK;
		zeiger1 = zeiger->s_n();
		if (zeiger1)
			strcat(str, ", ");
		zeiger = zeiger1;
		}
	strcat(str, ")");
	return OK;
}

INT LIST_OB::freeself()
{
	OBJECTSELF d; 
	INT erg = OK;
	LIST_OP z = this, za = NULL;
	
	if (s_obj_k() == EMPTY)
		return OK;
	while (z != NULL) {
		if (z->s_s() != NULL) {
			erg += freeall(z->s_s());
			z->c_s(NIL);
			}
		za = z;
		z = z->s_n();
		za->c_n(NIL);
		if (za != this) 
			freeall(za);
		else {
			d = za->s_obj_s();
			my_free(d.ob_list);
			}
		}
	c_obj_k(EMPTY);
	return erg;
}

INT LIST_OB::compare(LIST_OP b)
/* before: comp_list(OP a,OP b); */
/* vergleich zweier listen, z.b. 1,1,3  < 1,2,2 z.b. 2,2,3  > 2/3 */
{
	INT erg = s_s()->sym_comp(b->s_s());
	
	if (erg == 0) {
		if ((s_n() == NULL) && (b->s_n() == NULL))
			return(0); /* a == b */
		else if (s_n() == NULL)
			return(-1); /* a < b */
		else if (b->s_n() == NULL)
			return(1); /* a > b */
		else
			return s_n()->compare(b->s_n()); /* rest ist wieder liste */
		}
	else return erg;
}

INT LIST_OB::transform_apply(INT (*tf)(SYM_OP self))
/* before: transform_apply_list(OP von, INT (*tf)(OP zeiger)) */
{
	LIST_OP zeiger = this;
	
	while (zeiger != NULL) {
		(*tf)(zeiger->s_s());
		zeiger = zeiger->s_n();
		}
	return(OK);
}

INT LIST_OB::transform(LIST_OP nach, INT (*tf)(SYM_OP self_von, SYM_OP self_nach))
/* before: transformlist(OP von, OP nach, INT (*tf)(OP zeiger, OP nachzeiger)) */
/* nach wird vor dem Ueberschreiben freigegeben. */
{
	LIST_OP zeiger = this;
	LIST_OP nachzeiger = nach;
	OBJECTSELF d;
	
	while (zeiger != NULL) {
		d = zeiger->s_obj_s();
		if (d.ob_list == NULL)
			return error("LIST::transform: sos = NULL");
		if (zeiger->s_s() != NULL) {
			nachzeiger->b_sn(callocobject("LIST"), NULL);
			nachzeiger->c_obj_k(zeiger->s_obj_k());
			(*tf)(zeiger->s_s(), nachzeiger->s_s());
			}
		else {
			nachzeiger->b_sn(NULL, NULL);
			nachzeiger->c_obj_k(zeiger->s_obj_k());
			}
		if (! zeiger->lastp())
			nachzeiger->c_n((LIST_OP) callocobject("LIST"));
	
		zeiger = zeiger->s_n();
		nachzeiger = nachzeiger->s_n();
		}
	return(OK);
}


INT LIST_OB::copy(LIST_OP nach)
{
	OBJECTSELF d;
	INT erg = OK;
	
	d = s_obj_s();
	if (d.ob_list == NULL)
		return error("LIST::copy:sos = NULL");
	erg += transform(nach, copy_ab);
	nach->c_obj_k(s_obj_k());
	return erg;
}

INT LIST_OB::transform_by(SYM_OP ve, LIST_OP to,
	INT (*tf)(SYM_OP ve, SYM_OP self_from, SYM_OP self_to))
/* before: trans2formlist(OP ve, OP vz, OP nach,INT (*tf)(OP ve, OP zeiger, OP nachzeiger)) */
/* ve ist konstante , vz ist liste */
{
	LIST_OP zeiger = this;
	LIST_OP nachzeiger = to;
	INT erg;
	
	if (! to->emptyp())
		to->freeself();
	while (zeiger != NULL) {
		nachzeiger->b_sn(callocobject("LIST"), NULL);
		nachzeiger->c_obj_k(zeiger->s_obj_k());
		erg = (*tf)(ve, zeiger->s_s(), nachzeiger->s_s());
		if (erg == ERROR)
			return error("transform_by: function returns error");
		if (! zeiger->lastp()) {
			nachzeiger->c_n((LIST_OP) callocobject("LIST"));
			nachzeiger = nachzeiger->s_n(); 
			}
		zeiger = zeiger->s_n();
		}
	return(OK);
}

INT LIST_OB::insert(SYM_OP insert_this, 
	INT (*equality_handler)(SYM_OP self_from, SYM_OP self_to), 
	INT (*compare_func)(SYM_OP from, SYM_OP to))
/* before: insert_list(OP von, OP nach, INT (*eh)(OP von, OP nn), INT (*cf)(OP von, OP nn)) */
/* fuegt das object von in die liste nach ein */
/* von ist keine liste */
{
	LIST_OP c;
	
	if (insert_this->listp()) 
		return(insert_list((LIST_OP)insert_this, equality_handler, compare_func));
	c = (LIST_OP) callocobject("LIST");
#ifdef POLYNOM_TRUE
	if (s_ind_k() == POLYNOM) {
		if (scalarp(insert_this)) {
			b_skn_po(callocobject("LIST"),insert_this,NULL,c);
			m_il_v(1L,S_PO_S(c));
			m_i_i(0L,S_PO_SI(c,0L));
			}
		else { 
			b_sn_l(insert_this,NULL,c); 
			C_O_K(c,S_O_K(nach)); 
			}
		goto insert_c;
		}
#endif /* POLYNOMTRUE */
	c->b_sn(insert_this, NULL);
	c->c_obj_k(s_obj_k());
	return insert_list(c, equality_handler, compare_func);
}

#undef DEBUG_INSERT_LIST

INT LIST_OB::insert_list(LIST_OP insert_this_list, 
	INT (*equality_handler)(SYM_OP self_from, SYM_OP self_to), 
	INT (*compare_func)(SYM_OP from, SYM_OP to))
/* before: insert_list_list_2(OP von, OP nach, 
	INT (*eh)(OP von, OP nn), INT (*cf)(OP von, OP nn)) */
/* Das object insert_this_list wird verbraucht */
/* compare_func NULL means compare_ab() */
/* ersatz fuer insert_list_list programmiert nach
christopher J. van Wyk : Data  structures and c programs */
{
	LIST_OB dummy; /* an empty list object, the beginning of the merged list. */
	struct list dummy_list;
	LIST_OP p; /* the least element of the (partially) merged list */
	INT erg = OK;
	OBJECTSELF d;
	OBJECTKIND kind = insert_this_list->s_obj_k();
	LIST_OP nn, altnext;
	
#ifdef DEBUG_INSERT_LIST
	printf("in insert_list() (ANFANG) this =\n");
	println();
	printobjectkind();
	print_kind(s_obj_k());
	print_kind(s_ind_k());
	printf("emptyp = %ld\n", (INT) emptyp());
	fflush(stdout);
#endif

	if (emptyp()) {
#ifdef DEBUG_INSERT_LIST
		printf("in insert_list() empty object, calling init()\n");
		fflush(stdout);
#endif
		init(kind);
		}
	
	if (insert_this_list == NIL)
		return error("LIST::insert_list(): insert_this_list == NIL");
	
	if (s_s() == NULL) {
#ifdef DEBUG_INSERT_LIST
		printf("in insert_list() if (s_s() == NULL)\n");
		printf("insert_this_list=\n");
		insert_this_list->println();
		fflush(stdout);
#endif
		c_s(insert_this_list->s_s());
		c_n(insert_this_list->s_n());
		insert_this_list->c_s(NIL);
		insert_this_list->c_n(NIL);
		freeall(insert_this_list);
		/* d = insert_this_list->s_obj_s(); 
		my_free(d.ob_list); 
		my_free(insert_this_list); */
#ifdef DEBUG_INSERT_LIST
		printf("this = \n");
		println();
		fflush(stdout);
#endif
		if (s_s()->emptyp()) 
			freeself();
			/* first list element empty means whole list empty */
#ifdef DEBUG_INSERT_LIST
		printf("this = \n");
		println();
		fflush(stdout);
#endif
		return(OK);
		}
	if (insert_this_list->s_s() == NULL) {
		freeall(insert_this_list);
		return(OK);
		}


	if (s_s()->emptyp())
		return error("LIST::insert_list: result is a LIST with empty self");

#ifdef DEBUG_INSERT_LIST
	printf("in insert_list() this =\n");
	println();
	fflush(stdout);
#endif

	nn = (LIST_OP) callocobject("LIST");
	*nn = *this;
	this->c_obj_k(EMPTY); 
	
#ifdef DEBUG_INSERT_LIST
	printf("in insert_list() nn =\n");
	nn->println();
	fflush(stdout);
#endif

	p = &dummy;
	d.ob_list = &dummy_list;
	p->c_obj_s(d);
	p->c_obj_k(LIST);

	if (compare_func == NULL)
		compare_func = compare_ab;
	
	while ((insert_this_list != NULL) && (nn != NULL)) {
#if 0
#ifdef DEBUG_INSERT_LIST
		printf("in insert_list() vor compare:\n");
		printf("insert_this_list->s_s()=\n");
		insert_this_list->s_s()->println();
		printf("nn->s_s()=\n");
		nn->s_s()->println();
		fflush(stdout);
#endif
#endif
		erg = (*compare_func)(insert_this_list->s_s(), nn->s_s());
#if 0
#ifdef DEBUG_INSERT_LIST
		printf("erg = %ld\n", erg);
		fflush(stdout);
#endif
#endif
		if (erg < 0) {
			p->c_n(insert_this_list);
			insert_this_list = insert_this_list->s_n();
			p = p->s_n();
			}
		else if (erg > 0) {
			p->c_n(nn);
			nn = nn->s_n();
			p = p->s_n();
			}
		else {
			if (equality_handler == NULL) ;
			else {
#if 0
#ifdef DEBUG_INSERT_LIST
				printf("in insert_list() vor equality_handler:\n");
				printf("insert_this_list->s_s()=\n");
				insert_this_list->s_s()->println();
				printf("nn->s_s()=\n");
				nn->s_s()->println();
				fflush(stdout);
#endif
#endif
				(*equality_handler)(insert_this_list->s_s(), nn->s_s());
#if 0
#ifdef DEBUG_INSERT_LIST
				printf("in insert_list() nach equality_handler:\n");
				printf("insert_this_list->s_s()=\n");
				insert_this_list->s_s()->println();
				printf("nn->s_s()=\n");
				nn->s_s()->println();
				fflush(stdout);
#endif
#endif
				}
			if (! nn->s_s()->emptyp()) {
				/* equality_handler hat nicht geloescht, 
				 * nn wird uebernommen: */
				p->c_n(nn);
				p = p->s_n();
				nn = nn->s_n();
				}
			else {
				/* Element nn loeschen, nn auf Nachfolger setzen: */
				freeall(nn->s_s());
				altnext = nn->s_n();
				nn->c_s(NIL);
				nn->c_n(NIL);
				freeall(nn);
				/* d = nn->s_obj_s(); 
				free(d.ob_list);
				free(nn); */
				nn = altnext;
				}
			
			/* Element insert_this_list loeschen, 
			 * insert_this_list auf Nachfolger setzen: */
			freeall(insert_this_list->s_s());
			altnext = insert_this_list->s_n();
			insert_this_list->c_s(NIL);
			insert_this_list->c_n(NIL);
			freeall(insert_this_list);
			/* d = insert_this_list->s_obj_s(); 
			free(d.ob_list);
			free(insert_this_list); */
			insert_this_list = altnext;
			}
		} /* while */
	
	p->c_n(NULL); /* terminate the merged list. */
	
	/* now at least one of insert_this_list and nn is NULL: */
	if (insert_this_list == NULL) 
		insert_this_list = nn;
	if (insert_this_list != NULL)
		p->c_n(insert_this_list); /* append a remainder list */
	
	/* remember this is EMPTY at this point */
	if (dummy.s_n() == NULL) {
		init(kind); /* an empty list */
		}
	else { 
		*this = *(dummy.s_n());
		dummy.s_n()->c_obj_k(EMPTY);
		freeall(dummy.s_n());
		dummy.c_n(NIL);
		/* free(dummy.s_n()); */
		}
	dummy.c_obj_k(EMPTY); /* do not free anything in the destructor */
	return(OK);
}

#endif /* LISTTRUE */


