/* generators.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef GENERATORS_INCLUDED
#define GENERATORS_INCLUDED

#ifndef VDI_INCLUDED
#include <DISCRETA/graphics.h>
#endif
#ifndef MA_INCLUDED
#include <DISCRETA/ma.h>
#endif
#ifndef FGA_INCLUDED
#include <DISCRETA/fga.h>
#endif
#ifndef LABRA_INCLUDED
#include <DISCRETA/lb.h>
#endif
#ifndef DIVS_INCLUDED
#include <DISCRETA/divs.h> /* for STRING ->s_str() */
#endif
#ifndef DB_INCLUDED
#include <DISCRETA/db.h>
#endif

/* generators.C: */

class generators_ob : public VECTOR_OB {
public:
	INTEGER_OP s_id() { 
		return((INTEGER_OP) s_i(0)); };
	INT s_id_i() { 
		return(s_id()->s_i()); };
	INTEGER_OP s_deg() { 
		return((INTEGER_OP) s_i(1)); };
	INT s_deg_i() { 
		return(s_deg()->s_i()); };
	VECTOR_OP s_gen() { 
		return((VECTOR_OP) s_i(2)); };
	PERMUTATION_OP s_gen_i(INT i) { 
		return((PERMUTATION_OP) s_gen()->s_i(i)); };
	STRING_OP s_label() { 
		return((STRING_OP) s_i(3)); };
	BYTE *s_label_s() { 
		return(s_label()->s_str()); };
	STRING_OP s_label_tex() { 
		return((STRING_OP) s_i(4)); };
	BYTE *s_label_tex_s() { 
		return(s_label_tex()->s_str()); };
	STRING_OP s_order() { 
		return((STRING_OP) s_i(5)); };
	BYTE *s_order_s() { 
		return(s_order()->s_str()); };
	SYM_OP s_go() { 
		return((INTEGER_OP) s_i(6)); };

	INT field_name(INT i, INT j, BYTE *str);
	INT init();
	INT sprint(BYTE *s);
};

INT i_db_generators(DATABASE_OP *p_db, BYTE *path);
void e_db_generators(DATABASE_OP db);
INT do_db_generators_create_and_close(BYTE *path);
INT do_db_generators_info(BYTE *path, INT btree_idx);
INT do_db_generators_dump(BYTE *path, INT btree_idx);
INT db_generators_highest_id(DATABASE_OP db, INT *highest_id);
INT db_generators_load_id(DATABASE_OP db, INT id, GENERATORS_OP p);
INT db_generators_add(DATABASE_OP db, VECTOR_OP gen, 
	BYTE *label, BYTE *label_tex);


#endif /* GENERATORS_INCLUDEDE */

