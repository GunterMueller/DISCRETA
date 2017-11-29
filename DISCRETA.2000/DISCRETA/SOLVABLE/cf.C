/* cf.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SOLVABLE_TRUE

#include <DISCRETA/solvable.h>
#include <DISCRETA/lb.h>

#undef DEBUG_CALC_CF

/* 
 * GROUP_CANONIC_FORM_OB
 */

#if TEXDOCU
INT GROUP_CANONIC_FORM_OB::field_name(INT i, INT j, BYTE *str)
#endif
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "nb_gen"; break;
	case 1: s = "p0"; break;
	case 2: s = "p0v"; break;
	case 3: s = "g"; break;
	case 4: s = "go"; break;
	case 5: s = "ro"; break;
	case 6: s = "re"; break;
	case 7: s = "GG"; break;
	case 8: s = "ago_trunc"; break;
	case 9: s = "f_has_auts"; break;
	case 10: s = "f_auts_on_canonic_form"; break;
	case 11: s = "auts"; break;
	case 12: s = "first_moved"; break;
	default:
		return error("GROUP_CANONIC_FORM::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

#if TEXDOCU
INT GROUP_CANONIC_FORM_OB::init()
#endif
{
	
	m_il(13);
	c_obj_k(GROUP_CANONIC_FORM_KIND);
	s_nb_gen()->m_i(0);
	s_f_has_auts()->m_i(FALSE);
	s_ago_trunc()->m_i(0);
	return OK;
}

#if TEXDOCU
INT GROUP_CANONIC_FORM_OB::sprint(BYTE *s)
#endif
{
	BYTE str[512];
	
	sprintf(str, "GROUP_CANONIC_FORM: ");
	if (strlen(s) + strlen(str) < 200)
		strcat(s, str);
	return OK;
}

static INT f_calc_hash_first_time = TRUE;

#if TEXDOCU
INT cf_read_and_print_hash(VECTOR_OP hash)
#endif
{
	INT size, nb_ints;
	INT nb_gen, i, j, k, go, ago;
	INT ago_trunc;
	UBYTE *p_bytes = NIL, *p;
	INT4 *pi4, w4;
	UINT4 *ints = NIL;
	UINT w;

	nb_ints = hash->s_li();
	p_bytes = (UBYTE *) my_malloc(nb_ints * sizeof(INT4) + 
		16 /* for security */ , "cf.C: cf_read_and_print_hash");

	pi4 = (INT4 *) p_bytes;
	for (i = 0; i < nb_ints; i++) {
		w4 = hash->s_ii(i);
		pi4[i] = w4;
		}
	block_swap_bytes((SCHAR *) p_bytes + 8, sizeof(INT4), nb_ints - 2);
		/* we swap only the byte data !!! */
	go = pi4[0];
	printf("go = %ld\n", go);
	ago = pi4[1];
	printf("ago = %ld\n", ago);
	p = p_bytes + 8;
	nb_gen = *p++;
	printf("nb_gen = %ld\n", nb_gen);
	
	printf("go: ");
	for (i = 0; i < nb_gen; i++) {
		k = *p++;
		printf("%ld ", k);
		}
	printf("\n");
	
	printf("ro: ");
	for (i = 0; i < nb_gen; i++) {
		k = *p++;
		printf("%ld ", k);
		}
	printf("\n");
	
	printf("re: ");
	for (i = 0; i < nb_gen; i++) {
		k = *p++;
		printf("%ld ", k);
		}
	printf("\n");
	
	printf("GG: \n");
	for (i = 0; i < nb_gen; i++) {
		for (j = 0; j < nb_gen; j++) {
			k = *p++;
			printf("%ld ", k);
			}
		printf("\n");
		}
	printf("\n");
	fflush(stdout);
	my_free(p_bytes);
	return OK;
}

#if TEXDOCU
INT GROUP_CANONIC_FORM_OB::calc_hash(FILE *fp_txt, 
	VECTOR_OP hash, INT n, INT f_v, INT f_vv)
#endif
{
	INT size, nb_ints;
	INT nb_gen, i, j, k;
	INT ago_trunc;
	UBYTE *p_bytes = NIL, *p;
	UINT4 *pi4, w4;
	UINT4 *ints = NIL;
	UINT w;
	
#ifdef DEBUG_CALC_CF
	f_v = TRUE;
#endif

	ago_trunc = s_ago_trunc_i();
	nb_gen = s_nb_gen_i();
	size = 4 + 4 + 1 + nb_gen * 3 + nb_gen * nb_gen;
	nb_ints = (size >> 2) + 1;
	if (f_vv) {
		fprintf(fp_txt, "GROUP_CANONIC_FORM::calc hash(): n = %ld "
			"hash of size %ld bytes, %ld ints\n", n, size, nb_ints);
		fflush(fp_txt);
		}
	p_bytes = (UBYTE *) my_malloc(nb_ints * sizeof(INT4) + 
		16 /* for security */ , "cf.C: calc_hash");
	pi4 = (UINT4 *) p_bytes;
	*pi4++ = (UINT4) n;
	if (f_v)
		printf("go = %ld\n", n);
	*pi4++ = (UINT4) ago_trunc;
	if (f_v)
		printf("ago_trunc = %ld\n", ago_trunc);
	p = p_bytes + 8;
	*p++ = (BYTE) nb_gen;
	if (f_v)
		printf("nb_gen = %ld\n", nb_gen);
	for (i = 0; i < nb_gen; i++) {
		j = s_go_ii(i + 1);
		*p++ = (UBYTE)j;
		if (f_v)
			printf("%ld ", j);
		}
	if (f_v)
		printf("\n");
	for (i = 0; i < nb_gen; i++) {
		j = s_ro_ii(i);
		*p++ = (UBYTE)j;
		if (f_v)
			printf("%ld ", j);
		}
	if (f_v)
		printf("\n");
	for (i = 0; i < nb_gen; i++) {
		j = s_re_ii(i);
		*p++ = (UBYTE)j;
		if (f_v)
			printf("%ld ", j);
		}
	if (f_v)
		printf("\n");
	for (i = 0; i < nb_gen; i++) {
		for (j = 0; j < nb_gen; j++) {
			k = s_GG_iji(i, j);
			*p++ = (UBYTE)k;
			if (f_v)
				printf("%ld ", k);
			}
		if (f_v)
			printf("\n");
		}
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	if (f_v && f_calc_hash_first_time) {
		fprintf(fp_txt, "GROUP_CANONIC_FORM::calc_hash() \n"
		"first time, BEFORE swap: hash = \n");
		for (i = 0; i < size; i++) {
			j = p_bytes[i];
			fprintf(fp_txt, "%ld ", j);
			}
		fprintf(fp_txt, "\n");
		fflush(fp_txt);
		f_calc_hash_first_time = FALSE;
		}
	block_swap_bytes((SCHAR *) p_bytes + 8, sizeof(INT4), nb_ints - 2);
		/* we swap only the byte data !!! */
	if (f_v) {
		fprintf(fp_txt, "GROUP_CANONIC_FORM::calc_hash() hash = \n");
		for (i = 0; i < size; i++) {
			j = p_bytes[i];
			fprintf(fp_txt, "%ld ", j);
			}
		fprintf(fp_txt, "\n");
		fflush(fp_txt);
		}
	pi4 = (UINT4 *) p_bytes;
	hash->m_il(nb_ints);
	for (i = 0; i < nb_ints; i++) {
		w4 = pi4[i];
		hash->m_ii(i, w4);
		}
	if (f_v) {
		cf_read_and_print_hash(hash);
		}
	if (f_vv) {
		fprintf(fp_txt, "GROUP_CANONIC_FORM::calc hash(): hash =\n");
		hash->fprintln(fp_txt);
		fflush(stdout);
		}
	if (p_bytes) {
		my_free(p_bytes);
		p_bytes = NIL;
		}
	return OK;
}

#endif /* SOLVABLE_TRUE */

