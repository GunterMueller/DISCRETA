/* dimino.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef DIMINO_INCLUDED
#define DIMINO_INCLUDED

INT gruppen_elemente(VECTOR_OP G_gen, VECTOR_OP G, INT type, void *data);
INT gruppen_elemente1(VECTOR_OP G, INT type, void *data);
INT dimino_extend(VECTOR_OP G, VECTOR_OP G_gen, 
	SYM_OP g, VECTOR_OP bad_list, INT bad_list_len, 
	INT *f_bad_element_found, 
	INT type, void *data);
INT dimino_extend_normal(VECTOR_OP G, VECTOR_OP G_gen, SYM_OP g, INT type, void *data);
INT RelativeOrder(VECTOR_OP H, SYM_OP g, INT *relative_order, INT type, void *data);

#endif /* DIMINO_INCLUDED */
