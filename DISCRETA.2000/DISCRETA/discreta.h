/* discreta.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */



#ifndef DISCRETA_INCLUDED
#define DISCRETA_INCLUDED

#include <DISCRETA/discreta_types.h>

extern INT base_kind[SYM_MAX_KIND]; // discreta.C
extern BYTE *discreta_home;
extern BYTE *discreta_arch;

#if TEXDOCU
#else
SYM\_OB is the basic class of DISCRETA objects. 
A pointer to an instance of type SYM\_OB is called SYM\_OP.
Each DISCRETA class is derived from this class. 
Each DISCRETA class has an object kind, which is 
simply an integer which stands for the type. 
The SYM\_OB instances have two members:
//CODE
	OBJECTKIND ob_kind;
	OBJECTSELF ob_self;
///CODE
ob\_kind is an integer, OBJECTKIND, describing the kind (or type) 
of the object. The objectkind values for all DISCRETA classes are 
defined in discreta\_types.h.
Here are a few of them:
//CODE
#define EMPTY 0         /* 290590 */
#define INTEGER 1
#define VECTOR 2
#define PARTITION 3
#define BRUCH 4
#define PERMUTATION 6
#define SKEWPARTITION 7  /* 020488 */
#define TABLEAUX 8       /* 020488 */
#define POLYNOM 9
///CODE

ob\_self is the representation of the object. 
This representation depends heavily on the objectkind. 
For instance, an object of kind INTEGER is just the integer value 
itself. Other objects, like vector, may have deeply nested 
data describing the object.
OBJECTSELF is a union defined in discreta\_types.h 
which has the following definition:
//CODE
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
///CODE
The objectself is just one of these entries, 
detremined by the kind. For instance, an INTEGER has 
its value in ob\_self.ob\_INT. 
However, one should never access the data through 
pb\_self.something itself. There are access function in 
the DISCRETA classes which does the elementary access 
like storing and retrieving of data, allocating or freeing 
the memory for the representation.
DISCRETA has a destructor which cleans up all the different 
DISCRETA objects. 


Here is the definition of the class SYM\_OB:
#endif

#if TEXDOCU
class sym_ob {
public:
	OBJECTKIND ob_kind; // kind of object as integer
	OBJECTSELF ob_self; // representation of object

	sym_ob() { ob_kind = EMPTY; };
		// default constructor
	~sym_ob() { freeself(); };
		// default destructor, frees the memory
	
	OBJECTKIND s_obj_k() { return(ob_kind); };
		// select object kind
	INT s_obj_bk() { return base_kind[s_obj_k()]; };
		// select object base kind
	OBJECTSELF s_obj_s() { return(ob_self); };
		// select object self
	void c_obj_k(OBJECTKIND kind) { ob_kind = kind; };
		// compute object kind
	void c_obj_s(OBJECTSELF self) { ob_self = self; };
		// compute object self
	void b_ks(OBJECTKIND kind, 
		OBJECTSELF self) {
		if (!emptyp()) freeself();
		c_obj_k(kind);
		c_obj_s(self); };
		// build kind self object

	/* io.C: */
	INT field_name(INT i, INT j, BYTE *str);
	INT calc_len(INT *len);
	SYM_OP get_ijth(INT i, INT j);
	INT sscan(BYTE *s);
	INT print();
		// print the object
	INT println();
		// print with a new-line in the end
	INT fprint(FILE *of);
		// print to FILE
	INT fprintln(FILE *of);
		// print to FILE with newline
	INT sprint(BYTE *s);
		// print (append) to a string
	INT fprint_GAP(FILE *fp);
		// print to FILE in GAP format
		// this can be used as input for GAP
	INT latex(FILE *fp);
		// print to FILE in tex format
		// the output can be used as input for tex
	INT sprint_latex(BYTE *s);
		// print (append) to string in tex format
	
	/* io2.C: */
	INT sprintobjectkind(BYTE *s);
		// print the object kind to string
	INT printobjectkind();
		// print the object kind
	
	/* iof.C: */
	INT pack(MEM_OP mem, INT f_verbose, INT debug_depth);
	INT unpack(MEM_OP mem, INT f_verbose, INT debug_depth);
	INT calc_size_on_file();
	INT write_mem(MEM_OP mem, INT debug_depth);
	INT read_mem(MEM_OP mem, INT debug_depth);
	INT save(BYTE *fname);
		// save the object to a binary file 
	INT save_quiet(BYTE *fname);
		// save without noise
	INT load(BYTE *fname);
		// loads a binary file previously stored with save
	INT load_quiet(BYTE *fname);
		// loads quietly
	INT save_ascii(FILE *fp);
		// save as ASCII file
	INT load_ascii(FILE *fp);
		// loads ASCII file
	
	/*
	 * discreta.C:
	 */
	INT emptyp() { return(ob_kind == EMPTY); };
	INT swap(SYM_OP b);
		// swaps this and b
	INT freeself();
		// destroys the object (frees the memory)
	INT freeself_debug(INT print_depth);
	INT copy(SYM_OP b);
		// copy this to b
	INT nullp();
		// tests for 0
	INT einsp();
		// tests for 1
	INT negeinsp();
		// tests for -1
	INT sym_comp(SYM_OP b);
		// compares this with b
		// returns -1 if this < b, 0 if this = b, 1 if this > b
	INT eq(SYM_OP b);
		// tests if equal
	INT neq(SYM_OP b);
		// tests if not equal
	INT ge(SYM_OP b);
		// tests if greater than or equal
	INT gt(SYM_OP b);
		// tests if greater than
	INT le(SYM_OP b);
		// tests if less than or equal
	INT lt(SYM_OP b);
		// tests if less than
	INT init(OBJECTKIND kind);
	INT listp();

	INT power_int(SYM_OP res, INT l);
		// res becomes this to the l-th power
	INT power_int_apply(INT l);
		// this raised to the l-th power

	/*
	 * nu.C:
	 */
	INT order(INT *order);
	INT order_if_prime(INT *ord);
	INT order_if_prime_power(INT *ord, INT *prime, INT *k);


	INT add(SYM_OP b, SYM_OP d);
		// d := this + b
	INT add_apply(SYM_OP b);
		// b := this + b 
	INT addinvers(SYM_OP res);
		// res := -this
	INT addinvers_apply();
		// this := - this
	
	INT operator+=(SYM_OB &a);
	INT operator*=(SYM_OB &a);
	INT operator/=(SYM_OB &a);
	
	INT mult(SYM_OP b, SYM_OP d);
		// d := this * b
	INT mult_apply(SYM_OP b);
		// b := this * b
	INT invers(SYM_OP b);
		// b := this^{-1}
	INT invers_apply();
		// this := this^{-1}

	INT inc();
		// this := this + 1
		// a vector is extended by one element
	INT dec();
		// this := this - 1
	INT zero();
		// this := 0
	INT one();
		// this := 1
	INT m_one();
		// this := -1
	INT homo_z(INT z);
		// this := z
	INT sym_div(SYM_OP b, SYM_OP d);
	INT quores(SYM_OP b, SYM_OP c, SYM_OP d);
	/* c = ganzdiv(a,b)  d = mod(a,b) */
	INT ganzdiv(SYM_OP b, SYM_OP d);
	INT ganzdiv_integral(SYM_OP b, SYM_OP d);
	INT fakul(SYM_OP d);
	INT fakul_int(int i);
	INT sub(SYM_OP b, SYM_OP c);
	/* c = a - b */

	INT m_i_i(INT i);
		// this := i (as an INTEGER object)
	INT s_i_i();
		// select integer i 
		// (assumes that the object is of kind INTEGER)
};
#endif

/* io.C: */
INT print_kind(INT kind);
INT sprint_kind(INT kind, BYTE *str);
INT printeingabe(char *text);


/* discreta.C: */

#define SPEICHERSIZE 10L

extern SYM_OP *speicher;
extern INT speicherposition;

INT discreta_init(void);
INT discreta_exit(void); 

INT copy_ab(SYM_OP a, SYM_OP b);
INT compare_ab(SYM_OP a, SYM_OP b);
INT mult_abc(SYM_OP a, SYM_OP b, SYM_OP c);
INT addinvers_ab(SYM_OP a, SYM_OP b);
INT addinvers_apply_a(SYM_OP a);
SYM_OP callocobject(BYTE *memory_label);
INT freeall(SYM_OP a);

/* nu.C: */
INT n_choose_k_ob(SYM_OP n, SYM_OP k, SYM_OP a);
INT bezout(SYM_OP m, SYM_OP n, SYM_OP u, SYM_OP v, SYM_OP g);
/* g = u * m + v * n */
INT inverse_mod(SYM_OP a, SYM_OP b, SYM_OP m, INT *f_notinvertible);
INT prime_power_parts(SYM_OP g, VECTOR_OP gpp, INT type, void *data);
INT vec_to_vec_pp(VECTOR_OP V, VECTOR_OP Vpp, INT type, void *data);

#define DO_TYPE_SYM 0
#define DO_TYPE_PERM 1
// #define DO_TYPE_ENUM 2
// #define DO_TYPE_IPERM 3
#define DO_TYPE_FG 4

INT do_mult(SYM_OP a, SYM_OP b, SYM_OP c, INT type, void *data);
INT do_invers(SYM_OP a, SYM_OP b, INT type, void *data);
INT do_einsp(SYM_OP a, INT type, void *data);
INT do_one(SYM_OP a, INT type, void *data);
INT do_order(SYM_OP a, INT *order, INT type, void *data);
INT do_order_if_prime(SYM_OP a, INT *ord, INT type, void *data);
INT do_order_if_prime_power(SYM_OP a, 
	INT *ord, INT *prime, INT *k, 
	INT type, void *data);
INT do_power(SYM_OP a, INT i, SYM_OP b, INT type, void *data);
INT do_print(SYM_OP a, INT type, void *data);
INT do_sprint(SYM_OP a, BYTE *str, INT type, void *data);
INT do_latex(SYM_OP a, FILE *fp, INT type, void *data);
INT do_image_of(SYM_OP a, INT i, INT type, void *data);
INT do_comp(SYM_OP a, SYM_OP b, INT type, void *data);
INT do_copy(SYM_OP a, SYM_OP b, INT type, void *data);

#ifdef SINGER_TRUE
/* singer.C: */
int singer(long chi, long deg, long *pp);
#endif

// iof.C:
INT write_vec_file(VECTOR_OP V, BYTE *file_name, INT f_verbose, INT f_use_compress);
INT read_vec_file(VECTOR_OP V, BYTE *file_name, INT f_verbose, INT f_use_compress);
INT write_op_file(SYM_OP V, BYTE *file_name, INT f_verbose, INT f_use_compress);
INT read_op_file(SYM_OP V, BYTE *file_name, INT f_verbose, INT f_use_compress);

BYTE *Eostr(BYTE *s);
void Srfs(BYTE *s1, BYTE *s2);
void Srff(BYTE *s1, BYTE *s2);
INT tstropn(TSTRING **tp, BYTE *p);
INT tstrcls(TSTRING **tp);
BYTE *tstr(TSTRING *tp);
INT tstrlen(TSTRING *tp);
INT tstrfill(TSTRING **tp, INT len, INT c);
INT tstrcat(TSTRING **tp, BYTE *p);
INT tstrcpy(TSTRING **dst, TSTRING *src);
BYTE *eostr(BYTE *s);
BYTE upperchar(INT c);
BYTE lowerchar(INT c);
void strfill(BYTE *s, INT len, INT c);
void fill_char(void *v, INT cnt, INT c);
void fill_int(void *v, INT cnt, INT i);
void fill_long(void *v, INT cnt, INT l);
void move_char(void *s, void *d, INT cnt);
void move_int(void *s, void *d, INT cnt);
void move_long(void *s, void *d, INT cnt);
INT ij2k(INT i, INT j, INT n);
INT k2ij(INT k, INT *i, INT *j, INT n);
INT graph_add_edge(UINT *g, INT i, INT j, INT n);
INT graph_del_edge(UINT *g, INT i, INT j, INT n);
INT graph_has_edge(UINT g, INT i, INT j, INT n);
INT graph_code(VECTOR_OP V, UINT *g, INT n);
INT graph_decode(VECTOR_OP V, UINT g, INT n);
INT my_ggt(INT m, INT n);
void my_memcpy(void *dst, void *src, INT size);
INT s_scan_int(BYTE **s, INT *i);
INT s_scan_token(BYTE **s, BYTE *str);
INT s_scan_str(BYTE **s, BYTE *str);

extern BYTE *gl_hot_param; /* the command parameters */

INT parse_args(VECTOR_OP args);
void test_lo_iof();

// os.C:
typedef struct memory_entry MEMORY_ENTRY;

struct memory_entry {
	BYTE *comment;
	BYTE *file;
	int line;
	int num, size;
};


void test_swap(void);
void block_swap_bytes(SCHAR *ptr, INT size, INT no);
void f_read(char *buffer, int size, unsigned long items, FILE *ptr);
void f_write(char *buffer, int size, unsigned long items, FILE *ptr);
void f_gets(char *buffer, FILE *ptr);
void b_read(char *buffer, int size, unsigned long items, char **ptr);
void b_gets(char *buffer, char **ptr);
INT file_size(BYTE *name);
INT f_seek_set(FILE *f, INT offset);
INT f_seek_cur(FILE *f, INT offset);
INT f_seek_end(FILE *f, INT offset);
INT f_tell(FILE *f);

INT os_ticks(void);
INT os_ticks_system(void);
INT os_ticks_per_second(void);
INT os_ticks_to_hms(INT ticks, INT *h, INT *m, INT *s);
INT os_ticks_to_hms_tps(INT ticks, INT tps, INT *h, INT *m, INT *s);
void print_delta_time_100(INT l, BYTE *str);
void print_delta_time(INT l, BYTE *str);
void print_delta_time_tps(INT l, INT tps, BYTE *str);
void print_delta_time_tps_f_short(INT l, INT tps, BYTE *str, INT f_short);
void print_delta_time_tps_short(INT l, INT tps, BYTE *str);
INT call_system(BYTE *s);

#define my_malloc(a, comment) my_malloc5((a), FALSE, comment, __FILE__, __LINE__)
#define my_calloc(a, comment) my_malloc5((a), TRUE, comment, __FILE__, __LINE__)
void *my_malloc5(INT size, INT f_calloc, BYTE *comment, BYTE *file, int line);
void my_free(void *p);
void my_ptr_free(void **p);
INT memory_chain_length();
void *memory_chain_next(void *last);
void *memory_chain_p_next(void *p);
void *memory_chain_p_last(void *p);
BYTE *memory_chain_file(void *p);
BYTE *memory_chain_comment(void *p);
int memory_chain_line(void *p);
int memory_chain_size(void *p);
INT memory_usage();
int memory_entry_find(MEMORY_ENTRY *ME, INT l, 
	BYTE *comment, BYTE *file, int line, int *f_found);
void memory_chain_dump();

void wait_sec_4();
void wait_sec_2();
void wait_sec_n(INT n);
void date_as_string(BYTE *s);
INT runtime(long *l);
INT error(char *fehlertext);
INT f_error(FILE * fp_txt, char *fehlertext);
INT no_memory(void);


#include <DISCRETA/in.h>
#include <DISCRETA/vec.h>



#endif /* DISCRETA_INCLUDED */



