/* cp.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef CP_INCLUDED
#define CP_INCLUDED

#ifndef PERM_INCLUDED
#include <DISCRETA/perm.h>
#endif

typedef struct sperm SPERM;

struct sperm {
	INT l;
	UINT2 *a;
		/* a permutation of 
		 * { 0, 1 ... l - 1 } */
};

#define SPERM_MAX_N 4096

void sp_nil(SPERM *p);
INT sp_int(SPERM *p, INT l, BYTE *where);
INT sp_free(SPERM *p);
void sp_free_it(SPERM *p);
INT sp_mv(SPERM *p, SPERM *q);
INT sp_id(SPERM *p);
INT sp_mult(SPERM *a, SPERM *b, SPERM *c);
INT sp_inv(SPERM *a, SPERM *b);
INT sp_power(SPERM *a, SPERM *res, INT exp);
INT sp_mult_apply_tau_r(SPERM *a, INT i, INT j);
INT sp_mult_apply_tau_l(SPERM *a, INT i, INT j);
INT sp_mult_apply_forwc_r(SPERM *a, INT i, INT l);
INT sp_mult_apply_backwc_r(SPERM *a, INT i, INT l);
INT sp_mult_apply_forwc_l(SPERM *a, INT i, INT l);
INT sp_mult_apply_backwc_l(SPERM *a, INT i, INT l);
INT sp_onep(SPERM *p);
INT sp_cmp(SPERM *a, SPERM *b);
void sp_print(SPERM *p);
INT sp_sprint(SPERM *p, BYTE *s);
INT sp_latex(SPERM *p, FILE *fp);
INT sp_test(void);

typedef struct cperm CPERM;

struct cperm {
	INT l;
	UBYTE *a;
		/* a permutation of 
		 * { 0, 1 ... l - 1 } */
};

void cp_nil(CPERM *p);
INT cp_int(CPERM *p, INT l);
INT cp_free(CPERM *p);
void cp_free_it(CPERM *p);
INT cp_mv(CPERM *p, CPERM *q);
INT cp_id(CPERM *p);
INT cp_mult(CPERM *a, CPERM *b, CPERM *c);
/* erst a, dann b; Ergebnis nach c */
INT cp_inv(CPERM *a, CPERM *b);
/* b:= a^-1 */
INT cp_power(CPERM *a, CPERM *res, INT exp);
INT cp_mult_apply_tau_r(CPERM *a, INT i, INT j);
/* a := a (i j). */
INT cp_mult_apply_tau_l(CPERM *a, INT i, INT j);
/* a := (i j) a. */
INT cp_mult_apply_forwc_r(CPERM *a, INT i, INT l);
/* a := a (i i+1 ... i+l-1). */
INT cp_mult_apply_backwc_r(CPERM *a, INT i, INT l);
/* a := a (i+l-1 i+l-2 ... i+1 i). */
INT cp_mult_apply_forwc_l(CPERM *a, INT i, INT l);
/* a := (i i+1 ... i+l-1) a. */
INT cp_mult_apply_backwc_l(CPERM *a, INT i, INT l);
/* a := (i+l-1 i+l-2 ... i+1 i) a. */
INT cp_onep(CPERM *p);
INT cp_cmp(CPERM *a, CPERM *b);
void cp_print(CPERM *p);
INT cp_sprint(CPERM *p, BYTE *s);
/* haengt an s an. */
INT cp_latex(CPERM *p, FILE *fp);
INT cp_test(void);

class chunk_memh {
public:
	INT entries_used;
		/* -1 for: mem not initialized */
		/* next free entry in current chunk */
	INT entry_size; /* in BYTEs */
	INT entry_size_INT;
	INT chunk_size; /* # of entries per chunk */
	INT free_hdl; /* -1 for: no free handles */
	INT *mem; /* current chunk, = mem_tbl[times - 1]. */
	INT times; /* # of currently allocated chunks. */
	INT *mem_tbl[100000];
	INT f_verbose;

	chunk_memh();
	chunk_memh(INT entry_size, 
		INT chunk_size, INT f_verbose);
	~chunk_memh();
	void print_info(void);
	void *hdl2ptr(INT hdl);
	INT hdl_check(INT hdl);
	INT is_legal_hdl(INT hdl);
	INT chunk_alloc_hdl(void);
		/* returns -1 if no memory left. */
	void chunk_free_hdl(INT hdl);
	INT alloc_new_chunk(void);
};

/* 
 * IPERMs: permutations using chunk memory
 * (and INT handles therefore)
 */

INT iperm_test(void);
INT iperm_sprint(CHUNK_MEMH *cm, INT a, BYTE *s);
INT iperm_latex(CHUNK_MEMH *cm, INT a, FILE *fp);
INT iperm_print(CHUNK_MEMH *cm, INT a);
INT iperm_println(CHUNK_MEMH *cm, INT a);
INT iperm_print_vec(CHUNK_MEMH *cm, VECTOR_OP V);
INT iperm_alloc(CHUNK_MEMH *cm);
INT iperm_free(CHUNK_MEMH *cm, INT a);
INT iperm_free_vec(CHUNK_MEMH *cm, VECTOR_OP V);
INT iperm_id(CHUNK_MEMH *cm, INT a);
INT iperm_mv(CHUNK_MEMH *cm, INT a, INT b);
/* b := a */
INT iperm_mult(CHUNK_MEMH *cm, INT a, INT b, INT c);
/* erst a, dann b; Ergebnis nach c */
INT iperm_invers(CHUNK_MEMH *cm, INT a, INT b);
/* b := a^{-1} */
INT iperm_invers_apply(CHUNK_MEMH *cm, INT *p_a);
/* a := a^{-1} 
 * (the handle will change). */
INT iperm_onep(CHUNK_MEMH *cm, INT a);
INT iperm_cmp(CHUNK_MEMH *cm, INT a, INT b);
INT iperm_image_of(CHUNK_MEMH *cm, INT a, INT i);
INT cperm2iperm(CHUNK_MEMH *cm, CPERM *cp, INT a);
INT perm2iperm(CHUNK_MEMH *cm, PERMUTATION_OP p, INT a);
INT iperm2perm(CHUNK_MEMH *cm, INT a, PERMUTATION_OP p);
INT perm2iperm_vec(CHUNK_MEMH *cm, VECTOR_OP V, VECTOR_OP W);
INT iperm2perm_vec(CHUNK_MEMH *cm, VECTOR_OP V, VECTOR_OP W);

#endif /* CP_INCLUDED */
