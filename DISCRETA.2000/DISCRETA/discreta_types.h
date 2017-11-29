/* discreta_types.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef DISCRETA_TYPES_INCLUDED
#define DISCRETA_TYPES_INCLUDED

#define TEXDOCU 1

#define INT_HAS_4_BYTES
// define this if sizeof(int) is 4
// undef this if sizeof(long) is 4

#define USE_OLD_MALLOC
/* used in iof.C: my_malloc() / my_free();
 * if define, olf malloc() and free() are used for memory allocation;
 * new and delete  otherwise */

/* the default base-configuration: */
#define SYM_BASE
#define SYM_DATA
#define SYM_SCALAR
#define SYM_FGA
#define SYM_DCC
#define SYM_DB
#define SYM_GEO
#define GEO_TRUE
#define GRAPHICS_TRUE
#define SOLVABLE_TRUE
#define CODES_TRUE
#define DIMINO_TRUE			/* dimino.C */
#define SGL_TRUE			/* SGL */


#ifdef CONF_99

#undef RUCKDESCHEL_TRUE
#else

#undef RUCKDESCHEL_TRUE
#endif


typedef long INT;
typedef unsigned long UINT;
typedef long LONG;
typedef unsigned long ULONG;
typedef short SHORT;
typedef char BYTE;
typedef unsigned char UBYTE;
typedef char SCHAR;
typedef float FLOAT;
typedef BYTE TSTRING;

#ifdef INT_HAS_4_BYTES
typedef int INT4;
typedef unsigned int UINT4;
#else
typedef long INT4;
typedef unsigned long UINT4;
#endif

typedef short INT2;
typedef unsigned short UINT2;

#ifdef SYSTEMUNIX
#endif

#ifdef SYSTEMMSDOS
#define MSDOS
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>
#ifdef SYSTEMMAC
#include <stdlib.h>
#include <console.h>
#endif
#ifdef SYSTEMMSDOS
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <io.h>
#endif

/* to avoid malloc/free: */
/* #define malloc(a)	("a*a";++--)
#define free(a)     ("a*a";++--) */

#define MINIMUM(x, y)   ( ((x) < (y)) ?  (x) : (y) )
#define MAXIMUM(x, y)   ( ((x) > (y)) ?  (x) : (y) )
#define MIN(x, y)   ( ((x) < (y)) ?  (x) : (y) )
#define MAX(x, y)   ( ((x) > (y)) ?  (x) : (y) )
#define ABS(x)      ( ((x) <  0 ) ? (-(x)) : (x) )
#define EVEN(x)     ( ((x) % 2) == 0 )
#define ODD(x)      ( ((x) % 2) == 1 )

#define MAX_PERM_N 1024
#define VECTOR_OVERSIZE 16L


/* ANSI C: */
/* #define _(a) a */

#define NIL 0L

#ifndef TRUE
#define TRUE 1L
#endif
#ifndef FALSE
#define FALSE 0L
#endif

#ifndef OK
#define OK 0L
#endif
#ifndef ERROR
#define ERROR -1L
#endif
#define IMPROPER 1001L /* MD */

#define EQUAL 300792L 
	/* return value of check_equal_ */


/*
 * 
 */


#ifdef SYM_BASE
#define INTEGERTRUE 1			/* in.C, in1.C */
#endif

#ifdef SYM_DATA
#define DIVS_TRUE 1			/* divs.C: STRING_OB, MEM_OB, CONTI_OB */
#define LISTTRUE 1			/* list.C */
#define MATRIXTRUE 1			/* ma.C */
#define VECTORTRUE 1			/* vec.C */
#define UNIPOLYTRUE 1			/* unip.C */
#endif

#ifdef SYM_SCALAR
#define BRUCHTRUE 1			/* bruch.C */
#define GFQ_TRUE 1			/* gfq_psu.C, gfq_sz.C, gfq_zech.C, perm_rep.C */
#define LONGINTTRUE 1			/* lo.C */
#define SINGER_TRUE			/* singer.C */
#endif

#ifdef SYM_FGA
#define CP_TRUE				/* cp.C */
#define FGA_TRUE			/* fga.C, fga_gen.C */
#undef KONTEXT_TRUE			/* kontext.C */
#define LABRA_TRUE			/* lb.C */
#define MONOMTRUE 1			/* mon.C */
#define PARTTRUE 1			/* part.C */
#define PERMTRUE 1			/* perm.C */
#define POLYTRUE 1			/* poly.C */
#endif

#ifdef SYM_DCC
#define LADDER_TRUE			/* LADDER */
#define SOLID_TRUE			/* solid.C */
#endif


#ifdef SYM_DB
#define DB_TRUE			/* BT, BT_KEY, DB */
#define GENERATORS_TRUE			/* generators.C */
#endif







/* definitionen fuer object.c */

typedef struct object * OP;
typedef long OBJECTKIND;  

struct loc { 
	INT w2, w1, w0;
	struct loc *nloc;
};

struct longint {
	struct loc *floc;
	signed char signum; /* -1,0,+1 */
	INT laenge;
};


struct ganzdaten {
	INT basis, basislaenge, 
		vorrat, speicher, 
		auspos, auslaenge, auszz;
	unsigned char folgezeichen;
	struct loc *floc;
};

struct zahldaten {
	/* signed */ char ziffer[13];
		/* AB 941024 */
	INT mehr;
	INT ziffernzahl;
	struct loc  *fdez;
};

struct memory { 
	OP mem_length; 
	BYTE *mem_self; 
};

struct vector { 
	OP v_length; 
	OP v_self; 
};

struct list { 
	OP l_self; 
	OP l_next; 
};

struct partition { 
	OBJECTKIND pa_kind; 
	OP pa_self; 
};
#define LASTPARTITION 1234
#define EXPONENT (OBJECTKIND)88
#define FROBENIUS (OBJECTKIND)92
#define LASTCOMP 1234L

struct permutation { 
	OBJECTKIND p_kind; 
	OP p_self; 
};
#define LASTLEHMERCODE 12L
#define LASTPERMUTATION 13L
#define LAST_PERMUTATION 13L
#define LASTSHUFFLE 12048802L
#define ZYKEL (OBJECTKIND)40888	
	/* fuer zykelschreibweise */
#define BAR (OBJECTKIND)25		
	/* AK 260292 fuer barred perm */
#define BARCYCLE (OBJECTKIND)26		
	/* AK 260292 fuer barred perm */

struct monom { 
	OP mo_self;
	OP mo_koeff;  
};

struct bruch { 
	OP b_oben; 
	OP b_unten; 
	INT b_info; 
};
#define GEKUERZT 40892L
#define NGEKUERZT 408921L

struct sym_matrix { 
	OP m_length; 
	OP m_height; 
	OP m_self; 
};
#define SINGULAER 2903884


/* return value fuer insert */
#define INSERTEQ 301288  /* falls eq */
#define INSERTOK 3012881  /* falls insert */

typedef union {
	INT ob_INT;
	char *ob_INTpointer;
	char *ob_charpointer;
	struct bruch *ob_bruch;
	struct list *ob_list;
	struct longint *ob_longint;
	struct sym_matrix *ob_matrix;
	struct monom *ob_monom;
	struct partition *ob_partition;
	struct permutation *ob_permutation;
	struct vector *ob_vector;
	struct memory *ob_memory;
} OBJECTSELF;

struct object { 
	OBJECTKIND ob_kind; 
	OBJECTSELF ob_self;
};

#define SYM_MAX_KIND 512

/* the SYMMETRICA types 
 * smaller than 128 !!!!
 * see iof.C
 */

#define EMPTY 0         /* 290590 */
#define INTEGER 1
#define VECTOR 2
#define PARTITION 3
#define BRUCH 4
#define PERMUTATION 6
#define SKEWPARTITION 7  /* 020488 */
#define TABLEAUX 8       /* 020488 */
#define POLYNOM 9

#define SCHUR 10
#define MATRIX 11
#define AUG_PART 12
#define HOM_SYM 13
#define HOMSYM 13
#define SCHUBERT 14
#define INTEGERVECTOR 15
#define KOSTKA 16

#define SYMCHAR 18
#define WORD 19

#define LIST 20       /* 210688 */
#define MONOM 21      /*230688*/
#define LONGINT 22    /* 170888 */

#define BINTREE 24    /* 291288 */
#define GRAPH 25      /* 210889 */
#define COMP 26       /* 300889 */
#define KRANZTYPUS 27 /* 280390 */
#define POW_SYM 28
#define MONOMIAL 29   /* 090992 */

#define BTREE 30
#define KRANZ 31
#define GRAL 32         /* 200691 */
#define ELM_SYM 33      /* 090992 */
#define FINITEFIELD  35 /* 250193 */
#define FF  35 		/* 250193 */
#define SGL 36			/* AB 311294 */
#define SGO 37			/* AB 311294 */
#define MEM 39          		/* AB 110893, 171293 */

#define DCY 40 	        /* AB 311294 */

#define STRING 44       /* AB 081093 */





#define DATABASE 50                    /* AB 161293 */
#define BAYERTREE 51                   /* AB 161293 */
#define ZE 52                          /* AB 251293 */
#define CONTI 53                       /* AB 251293 */
#define FG 54                          /* AB 261293 */

#define LABRA_KIND 56                  /* AB 120494 */

#define CODE 58                        /* AB 160494 */



#define BT_KEY_KIND 61                 /* AB 290494 */
#define CODE_ESSENTIALS_KIND 62        /* AB 020594 */
#define GED_KIND 63                    /* AB 120694 */



#define DESIGN_PARAMETER_KIND 70            /* AB 050196 */
#define KONTEXT_KIND 71                     /* AB 140296 */
#define CLASS_REP_KIND 72                   /* AB 260296 */
#define GROUP_CANONIC_FORM_KIND 73          /* AB 090396 */



#define BITVEC_KIND 77                      /* AB 310596 */
#define GROUP_SELECTION_KIND 78             /* AB 161196 */
#define UNIPOLY 79                          /* AB 100197 */

#define TREE_KIND 80                        /* AB 130297 */
#define GEO_BY_BASE_BLOCKS_KIND 81          /* AB 220697 */
#define GENERATORS_KIND 82                  /* AB 280797 */
#define DESIGN_PARAMETER_SOURCE_KIND 83     /* AB 290998 */
#define SOLID_KIND 84				/* EH 310399 */
#define M_V_KONTEXT_KIND 85                 /* BW 101098 */
#define L_F_KONTEXT_KIND 86                 /* BW 281298*/


typedef class sym_ob SYM_OB;
typedef SYM_OB *SYM_OP;
typedef class integer_ob INTEGER_OB;
typedef INTEGER_OB *INTEGER_OP;
typedef class vector_ob VECTOR_OB;
typedef VECTOR_OB *VECTOR_OP;
typedef class matrix_ob MATRIX_OB;
typedef MATRIX_OB *MATRIX_OP;
typedef class permutation_ob PERMUTATION_OB;
typedef PERMUTATION_OB *PERMUTATION_OP;


typedef class string_ob STRING_OB;
typedef STRING_OB *STRING_OP;
typedef class mem_ob MEM_OB;
typedef MEM_OB *MEM_OP;
typedef class conti_ob CONTI_OB;
typedef CONTI_OB *CONTI_OP;
typedef class ze_ob ZE_OB;
typedef ZE_OB *ZE_OP;
typedef class fg_ob FG_OB;
typedef FG_OB *FG_OP;
typedef class bt_ob BAYERTREE_OB;
typedef BAYERTREE_OB *BAYERTREE_OP;
typedef class db_ob DATABASE_OB;
typedef DATABASE_OB *DATABASE_OP;
typedef class list_ob LIST_OB;
typedef LIST_OB *LIST_OP;
typedef class monom_ob MONOM_OB;
typedef MONOM_OB *MONOM_OP;
typedef class polynom_ob POLYNOM_OB;
typedef POLYNOM_OB *POLYNOM_OP;
typedef class partition_ob PARTITION_OB;
typedef PARTITION_OB *PARTITION_OP;
typedef class labra_ob LABRA_OB;
typedef LABRA_OB *LABRA_OP;
typedef class code_ob CODE_OB;
typedef CODE_OB *CODE_OP;
typedef class bt_key_ob BT_KEY_OB;
typedef BT_KEY_OB *BT_KEY_OP;
typedef class code_essentials_ob CODE_ESSENTIALS_OB;
typedef CODE_ESSENTIALS_OB *CODE_ESSENTIALS_OP;
typedef class ged_ob GED_OB;
typedef GED_OB *GED_OP;
typedef class longint_ob LONGINT_OB;
typedef LONGINT_OB *LONGINT_OP;
typedef class bruch_ob BRUCH_OB;
typedef BRUCH_OB *BRUCH_OP;
typedef class sgo_ob SGO_OB;
typedef SGO_OB *SGO_OP;
typedef class sgl_ob SGL_OB;
typedef SGL_OB *SGL_OP;
typedef class dcy_ob DCY_OB;
typedef DCY_OB *DCY_OP;
typedef class km_data_ob KM_DATA_OB;
typedef KM_DATA_OB *KM_DATA_OP;
typedef class design_parameter_ob DESIGN_PARAMETER_OB;
typedef DESIGN_PARAMETER_OB *DESIGN_PARAMETER_OP;
typedef class kontext_ob KONTEXT_OB;
typedef KONTEXT_OB *KONTEXT_OP;
typedef class class_rep_ob CLASS_REP_OB;
typedef CLASS_REP_OB *CLASS_REP_OP;
typedef class group_canonic_form_ob GROUP_CANONIC_FORM_OB;
typedef GROUP_CANONIC_FORM_OB *GROUP_CANONIC_FORM_OP;
typedef class group_info_ob GROUP_INFO_OB;
typedef GROUP_INFO_OB *GROUP_INFO_OP;
typedef class set_orbits_ob SET_ORBITS_OB;
typedef SET_ORBITS_OB *SET_ORBITS_OP;
typedef class km_matrix_ob KM_MATRIX_OB;
typedef KM_MATRIX_OB *KM_MATRIX_OP;
typedef class bitvec_ob BITVEC_OB;
typedef BITVEC_OB *BITVEC_OP;
typedef class sgo_info_ob SGO_INFO_OB;
typedef SGO_INFO_OB *SGO_INFO_OP;
typedef class group_selection_ob GROUP_SELECTION_OB;
typedef GROUP_SELECTION_OB *GROUP_SELECTION_OP;
typedef class unipoly_ob UNIPOLY_OB;
typedef UNIPOLY_OB *UNIPOLY_OP;
typedef class tree_ob TREE_OB;
typedef TREE_OB *TREE_OP;
typedef class geo_by_base_blocks_ob GEO_BY_BASE_BLOCKS_OB;
typedef GEO_BY_BASE_BLOCKS_OB *GEO_BY_BASE_BLOCKS_OP;
typedef class generators_ob GENERATORS_OB;
typedef GENERATORS_OB *GENERATORS_OP;
typedef class design_parameter_source_ob DESIGN_PARAMETER_SOURCE_OB;
typedef DESIGN_PARAMETER_SOURCE_OB *DESIGN_PARAMETER_SOURCE_OP;
typedef class solid_ob SOLID_OB;
typedef SOLID_OB *SOLID_OP;

typedef class chunk_memh CHUNK_MEMH;
/* cp.h */



//==============================================================================

//  M_V_KONTEXT (BW)
//
//==============================================================================

typedef class m_v_kontext_ob M_V_KONTEXT_OB;
typedef M_V_KONTEXT_OB *M_V_KONTEXT_OP;

typedef class scale_ob SCALE_OB;
typedef SCALE_OB    *SCALE_OP;

//==============================================================================
//
//  L_F_KONTEXT (BW)
//
//==============================================================================

typedef class l_f_kontext_ob L_F_KONTEXT_OB;
typedef L_F_KONTEXT_OB *L_F_KONTEXT_OP;


#endif /* DISCRETA_TYPES_INCLUDED */

