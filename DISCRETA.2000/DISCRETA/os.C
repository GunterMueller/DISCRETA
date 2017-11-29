/* os.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>


/*
 * OS - library
 */

#include <math.h>
#include <fcntl.h>
#ifdef SYSTEMUNIX
#include <fcntl.h>
/* #include <malloc.h> */
#include <unistd.h>
	/* for sysconf */
#include <limits.h>
	/* for CLK_TCK */
#include <sys/types.h>
#include <sys/times.h>
	/* for times() */
#include <stdlib.h>
#include <time.h>
	/* for time() */
#endif
#ifdef SYSTEMWINDOWS
#include <io.h>
#include <process.h>
#endif
#ifdef SYSTEMMAC
#include <console.h>
#include <time.h> // for clock() 
#include <unix.h>
#endif
#ifdef MSDOS
#include <time.h> // for clock()
#endif


/* the following two routines 
 * are meant for binary input and output
 * in accordance with a 
 * byte swap on a machine 
 * We assume that the size of 
 * unsigned long is always 4!! 
 * RG: Roland Grund */

/* 
 * we know the swap behaviour of the following machines:
 * DEC alpha          has a swap
 * LINUX PC           has a swap
 * Silicon Graphics   no swap
 * IBM RS6000         no swap
 */

static INT is_swap = 0;
	/* indicates a byte swap */

#if TEXDOCU
void test_swap(void)
#endif
{
    unsigned long test_long = 0x11223344L;
    SCHAR *ptr;
    
    ptr = (char *) &test_long;
    is_swap = (ptr[0] == 0x44);
}

/* block_swap_bytes:
 * turns round the bytes within 
 * "no" intervalls of "size" in the 
 * buffer pointed to by "ptr" */

#if TEXDOCU
void block_swap_bytes(SCHAR *ptr, INT size, INT no)
#endif
{
	SCHAR *ptr_end, *ptr_start;
	SCHAR chr;
	INT i;
	
	if ((is_swap) && (size > 1)) {

		for(; no--; ) {
	
			ptr_start = ptr;
			ptr_end = ptr_start + (size - 1);
			for(i = size / 2; i--; ) {
				chr = *ptr_start;
				*ptr_start++ = *ptr_end;
				*ptr_end-- = chr;
				}
			ptr += size;
			}
		}
}

#if TEXDOCU
void f_read(char *buffer, int size, unsigned long items, FILE *ptr)
#endif
{
	fread(buffer, size, items, ptr);

	block_swap_bytes(buffer, size, items);
}

#if TEXDOCU
void f_write(char *buffer, int size, unsigned long items, FILE *ptr)
#endif
{
	char *tmp_buf;

	if (size > 1) {
		tmp_buf = (char *) my_malloc(size * items * sizeof(char), "fwrite()");
		memcpy(tmp_buf, buffer, size * items);

		block_swap_bytes(tmp_buf, size, items);
		
		fwrite(tmp_buf, size, items, ptr);
		my_free(tmp_buf);
		}
	else
		fwrite(buffer, size, items, ptr);
}

#if TEXDOCU
void f_gets(char *buffer, FILE *ptr)
#endif
{
	fscanf(ptr, "%s", buffer);
	fseek(ptr, SEEK_CUR, 1);
}

#if TEXDOCU
void b_read(char *buffer, int size, unsigned long items, char **ptr)
#endif
{
	memcpy(buffer, *ptr, size * items);
	block_swap_bytes(buffer, size, items);
	(*ptr) += size * items;
}

#if TEXDOCU
void b_gets(char *buffer, char **ptr)
#endif
{
	sscanf(*ptr, "%s", buffer);
	(*ptr) += strlen(*ptr) + 1;
}


#if TEXDOCU
INT file_size(BYTE *name)
#endif
{
#ifdef SYSTEMUNIX
	INT handle, size;
	
	handle = open(name, O_RDWR/*mode*/);
	size = lseek(handle, 0L, SEEK_END);
	close(handle);
	return(size);
#endif
#ifdef SYSTEMMAC
	INT handle, size;
	
	handle = open(name, O_RDONLY);
		/* THINK C Unix Lib */
	size = lseek(handle, 0L, SEEK_END);
		/* THINK C Unix Lib */
	close(handle);
	return(size);
#endif
#ifdef SYSTEMWINDOWS
	INT handle = _open (name,_O_RDONLY);
	INT size   = _lseek (handle,0,SEEK_END);
	close (handle);
	return (size);
#endif
}

#if TEXDOCU
INT f_seek_set(FILE *f, INT offset)
#endif
{
	return(fseek(f, offset, SEEK_SET));
}

#if TEXDOCU
INT f_seek_cur(FILE *f, INT offset)
#endif
{
	return(fseek(f, offset, SEEK_CUR));
}

#if TEXDOCU
INT f_seek_end(FILE *f, INT offset)
#endif
{
	return(fseek(f, offset, SEEK_END));
}

#if TEXDOCU
INT f_tell(FILE *f)
#endif
{
	return((INT)ftell(f));
}

#if TEXDOCU
INT os_ticks()
#endif
{
#ifdef SYSTEMMAC
	clock_t t;
	
	t = clock();
	return((INT)t);
#endif
#ifdef SYSTEMUNIX
	struct tms tms_buffer;

	if (-1 == times(&tms_buffer))
		return(-1);
	return(tms_buffer.tms_utime);
#endif
	return(0);
}

static INT system_time0 = 0;

#if TEXDOCU
INT os_ticks_system()
#endif
{
#ifdef SYSTEMMAC
	clock_t t;
	
	t = clock();
	return((INT)t);
#endif
#ifdef SYSTEMUNIX
#if 0
	struct tms tms_buffer;

	if (-1 == times(&tms_buffer))
		return(-1);
	return(tms_buffer.tms_stime);
#endif
	INT t;

	t = time(NULL);
	if (system_time0 == 0) {
		system_time0 = t;
		}
	t -= system_time0;
	t *= os_ticks_per_second();
	return t;
#endif
	return(0);
}

#if TEXDOCU
INT os_ticks_per_second()
#endif
{
	INT clk_tck = 1;
	
#ifdef SYSTEMUNIX
	clk_tck = sysconf(_SC_CLK_TCK);
	/* printf("clk_tck = %ld\n", clk_tck); */
#endif
#ifdef SYSTEMMAC
	clk_tck = CLOCKS_PER_SEC;
#endif
	return(clk_tck);
}

#if TEXDOCU
INT os_ticks_to_hms(INT ticks, INT *h, INT *m, INT *s)
#endif
{
	os_ticks_to_hms_tps(ticks, os_ticks_per_second(), h, m, s);
	return OK;
}

#if TEXDOCU
INT os_ticks_to_hms_tps(INT ticks, INT tps, INT *h, INT *m, INT *s)
#endif
{
	INT l1;

	l1 = ticks / tps;
	*s = l1 % 60;
	l1 -= *s;
	l1 /= 60;
	*m = l1 % 60;
	l1 -= *m;
	l1 /= 60;
	*h = l1;
	return(OK);
}

#if TEXDOCU
void print_delta_time_100(INT l, BYTE *str)
#endif
{
	INT tps = os_ticks_per_second();
	INT hs = l % tps;
	hs = (INT) ((100 * hs) / tps);
	
	print_delta_time_tps(l, tps, str);
	sprintf(str + strlen(str), ".%02ld", hs);
}

#if TEXDOCU
void print_delta_time(INT l, BYTE *str)
#endif
{
	print_delta_time_tps(l, os_ticks_per_second(), str);
}

#if TEXDOCU
void print_delta_time_tps(INT l, INT tps, BYTE *str)
#endif
{
	INT h, m, s;

	os_ticks_to_hms_tps(l, tps, &h, &m, &s);
	sprintf(Eostr(str), "%ld:%02ld:%02ld", h, m, s);
}

#if TEXDOCU
void print_delta_time_tps_f_short(INT l, INT tps, BYTE *str, INT f_short)
#endif
{
	if (f_short) {
		print_delta_time_tps_short(l, tps, str);
		}
	else {
		print_delta_time_tps(l, tps, str);
		}
}

#if TEXDOCU
void print_delta_time_tps_short(INT l, INT tps, BYTE *str)
#endif
{
	INT h, m, s;

	os_ticks_to_hms_tps(l, tps, &h, &m, &s);
	if (h == 0) {
		if (m <= 1)
			return;
		if (m == 0) {
			sprintf(Eostr(str), "%ld", s);
			}
		else {
			sprintf(Eostr(str), "%ld:%02ld", m, s);
			}
		}
	else {
		sprintf(Eostr(str), "%ld:%02ld:%02ld", h, m, s);
		}
}

#ifdef malloc
#undef malloc
#endif
#ifdef free
#undef free
#endif


#if TEXDOCU
INT call_system(BYTE *s)
#endif
{
	printf("system('%s')\n", s);
	fflush(stdout);
	system(s);
	return OK;
}

static void *memory_list = NIL;
// the top element in the chain of memory chunks

#if TEXDOCU
void *my_malloc5(INT size, INT f_calloc, BYTE *comment, BYTE *file, int line)
#endif
{
	void *p0, *p, **pv;
	BYTE **pp, *pc;
	int *pi, i;
	int size0 = (int) size;
	
	if (size < 0) {
		Srfs("my_malloc5", "size < 0");
		return NIL;
		}
	if (size == 0) {
		/* Srfs("my_malloc", 
			"warning: size == 0") */ ;
		}
	size += 2 * sizeof(void *);
	size += sizeof(BYTE *);
	size += sizeof(BYTE *);
	size += 2 * sizeof(int);
#ifdef USE_OLD_MALLOC
	p = malloc(size);
#else
	p = new char[size];
#endif
	p0 = p;
	pc = (BYTE *)p;
	if (p == NIL) {
		Srfs("my_malloc5", "no memory");
		return NIL;
		}
	if (f_calloc) {
		for (i = 0; i < size; i++)
			pc[i] = 0;
		}
	
	// set the next pointer:
	pv = (void **) pc;
	*pv = memory_list;
	pc += sizeof(void *);
	
	// set the last pointer:
	pv = (void **) pc;
	*pv = NIL;
	pc += sizeof(void *);

	if (memory_list) {
		// update the last pointer of the 
		// previously top element:
		pv = (void **) memory_list;
		pv[1] = p0;
		}
	
	// our new chunk will be the new top element
	memory_list = p0;
	
	pp = (BYTE **) pc;
	pp[0] = comment;
	pp[1] = file;
	pc += 2 * sizeof(BYTE *);
	
	pi = (int *) pc;
	pi[0] = line;
	pi[1] = size0;
	pc += 2 * sizeof(int);
	return pc;
}

#if TEXDOCU
void my_free(void *p)
#endif
{
	BYTE *pc;
	void *p_next, *p_last;
	void **pv;
	
	if (p == NIL) {
		Srfs("my_free", "p == NIL");
		return;
		}
	pc = (BYTE *) p;
	pc -= 2 * sizeof(int);
	pc -= 2 * sizeof(BYTE *);
	pc -= 2 * sizeof(void *);
	pv = (void **) pc;
	p_next = pv[0]; // possibly NIL
	p_last = pv[1]; // possibly NIL if it is the top chunk
	if (p_next) {
		pv = (void **) p_next;
		pv[1] = p_last;
		}
	if (p_last) {
		pv = (void **) p_last;
		pv[0] = p_next;
		}
	else {
		memory_list = p_next;
		}
	
#ifdef USE_OLD_MALLOC
	free(pc);
#else
	delete pc;
#endif
}

#if TEXDOCU
void my_ptr_free(void **p)
#endif
{
	if (p == NIL)
		error("my_ptr_free(): p == NIL");
	if (*p) {
		my_free(*p);
		*p = NIL;
		}
}

#if TEXDOCU
INT memory_chain_length()
#endif
{
	void *p;
	INT i = 0;
	
	p = memory_chain_next(NIL);
	while (p != NIL) {
		p = memory_chain_next(p);
		i++;
		}
	return i;
}

#if TEXDOCU
void *memory_chain_next(void *last)
#endif
{
	if (last == NIL)
		return memory_list;
	return memory_chain_p_next(last);
}

#if TEXDOCU
void *memory_chain_p_next(void *p)
#endif
{
	void **pv = (void **) p;
	
	return pv[0];
}

#if TEXDOCU
void *memory_chain_p_last(void *p)
#endif
{
	void **pv = (void **) p;
	
	return pv[1];
}

#if TEXDOCU
BYTE *memory_chain_comment(void *p)
#endif
{
	void **pv = (void **) p;
	BYTE **pc;
	
	pc = (BYTE **) &pv[2];
	return pc[0];
}

#if TEXDOCU
BYTE *memory_chain_file(void *p)
#endif
{
	void **pv = (void **) p;
	BYTE **pc;
	
	pc = (BYTE **) &pv[2];
	return pc[1];
}

#if TEXDOCU
int memory_chain_line(void *p)
#endif
{
	void **pv = (void **) p;
	BYTE **pc;
	int *pi;
	
	pc = (BYTE **) &pv[2];
	pi = (int *) &pc[2];
	return pi[0];
}

#if TEXDOCU
int memory_chain_size(void *p)
#endif
{
	void **pv = (void **) p;
	BYTE **pc;
	int *pi;
	
	pc = (BYTE **) &pv[2];
	pi = (int *) &pc[2];
	return pi[1];
}

#if TEXDOCU
INT memory_usage()
#endif
{
	void *p;
	INT Size = 0;
	void *p_next, *p_last;
	BYTE *comment, *file;
	int line, size;
	int i = 0, j, f_found;
	INT l, ll;
	MEMORY_ENTRY *ME, *ME1;
	int size_comment, size_file;
	BYTE str1[1024];
	BYTE str2[1024];
	
	p = memory_chain_next(NIL);
	while (p != NIL) {
		p_next = memory_chain_p_next(p);
		p_last = memory_chain_p_last(p);
		comment = memory_chain_comment(p);
		file = memory_chain_file(p);
		line = memory_chain_line(p);
		size = memory_chain_size(p);
		Size += size;
		// printf("memory chunk %d: comment%s file=%s line=%d size=%d %ld\n", 
		// 	i, comment, file, line, size, Size);


		p = memory_chain_next(p);
		i++;
		}
	l = i;
	printf("\noverall memory usage: %ld bytes in %ld chunks\n", Size, l);
	
	ll = 0;
	ME = (MEMORY_ENTRY *) malloc(sizeof(MEMORY_ENTRY) * (ll + 1));
	p = memory_chain_next(NIL);
	while (p != NIL) {
		comment = memory_chain_comment(p);
		file = memory_chain_file(p);
		line = memory_chain_line(p);
		size = memory_chain_size(p);
		i = memory_entry_find(ME, ll, comment, file, line, &f_found);
		if (f_found) {
			ME[i].num++;
			ME[i].size += size;
			}
		else {
			ME1 = (MEMORY_ENTRY *) malloc(sizeof(MEMORY_ENTRY) * (ll + 1));
			for (j = 0; j < i; j++) {
				ME1[j] = ME[j];
				}
			for (j = i; j < ll; j++) {
				ME1[j + 1] = ME[j];
				}
			free(ME);
			ME1[i].comment = comment;
			ME1[i].file = file;
			ME1[i].line = line;
			ME1[i].num = 1;
			ME1[i].size = size;
			ME = ME1;
			ll++;
			}
		p = memory_chain_next(p);
		}
	
	printf("used in:\n");
	size_comment = 0;
	size_file = 0;
	for (i = 0; i < ll; i++) {
		size_comment = MAX(size_comment, strlen(ME[i].comment));
		size_file = MAX(size_file, strlen(ME[i].file));
		}
	for (i = 0; i < ll; i++) {
		strcpy(str1, ME[i].comment);
		strfill(str1, size_comment, ' ');
		strcpy(str2, ME[i].file);
		strfill(str2, size_file, ' ');
		printf("%s / %s / %4d: %5d times  %10d bytes\n", 
			str1, str2, ME[i].line, 
			ME[i].num, ME[i].size);
		}
	free(ME);
	return Size;
}

#if TEXDOCU
int memory_entry_find(MEMORY_ENTRY *ME, INT l, 
	BYTE *comment, BYTE *file, int line, int *f_found)
#endif
{
	int i, r;
	
	for (i = 0; i < l; i++) {
		r = strcmp(ME[i].file, file);
		if (r > 0) {
			*f_found = FALSE;
			return i;
			}
		if (r == 0) {
			r = ME[i].line - line;
			if (r > 0) {
				*f_found = FALSE;
				return i;
				}
			if (r == 0) {
				r = strcmp(ME[i].comment, comment);
				if (r > 0) {
					*f_found = FALSE;
					return i;
					}
				if (r == 0) {
					*f_found = TRUE;
					return i;
					}
				}
			}
		}
	*f_found = FALSE;
	return l;
}

#if TEXDOCU
void memory_chain_dump()
#endif
{
	void *p;
	INT Size = 0;
	void *p_next, *p_last;
	BYTE *comment, *file;
	int line, size;
	int i = 0;
	
	p = memory_chain_next(NIL);
	while (p != NIL) {
		p_next = memory_chain_p_next(p);
		p_last = memory_chain_p_last(p);
		comment = memory_chain_comment(p);
		file = memory_chain_file(p);
		line = memory_chain_line(p);
		size = memory_chain_size(p);
		Size += size;
		printf("memory chunk %d: comment%s file=%s line=%d size=%d %ld\n", 
			i, comment, file, line, size, Size);


		p = memory_chain_next(p);
		i++;
		}
}

#if TEXDOCU
void wait_sec_4()
#endif
{
	INT t0, t1, clk_tck;

	t0 = os_ticks();
	clk_tck = os_ticks_per_second();
	clk_tck >>= 2;
	while (TRUE) {
		t1 = os_ticks();
		if (t1 - t0 > clk_tck)
			return;
		}
}

#if TEXDOCU
void wait_sec_2()
#endif
{
	INT t0, t1, clk_tck;

	t0 = os_ticks();
	clk_tck = os_ticks_per_second();
	clk_tck >>= 1;
	while (TRUE) {
		t1 = os_ticks();
		if (t1 - t0 > clk_tck)
			return;
		}
}

#if TEXDOCU
void wait_sec_n(INT n)
#endif
{
	INT t0, t1, clk_tck;

	t0 = os_ticks();
	clk_tck = os_ticks_per_second();
	clk_tck *= n;
	while (TRUE) {
		t1 = os_ticks();
		if (t1 - t0 > clk_tck)
			return;
		}
}

#define BUF_SIZE 10000

#if TEXDOCU
void date_as_string(BYTE *s)
#endif
{
	BYTE str[BUF_SIZE];
	FILE *fp;
	INT l;
	
	system("date >a");
	fp = fopen("a", "r");
	fgets(str, BUF_SIZE, fp);
	fclose(fp);
	l = strlen(str);
	if (l && str[l - 1] == '\n')
		str[l - 1] = 0;
	strcpy(s, str);
}

#if TEXDOCU
INT runtime(long *l)
#endif
{
#ifdef SYSTEMUNIX
	struct tms *buffer = (struct tms *) my_malloc(sizeof(struct tms), "runtime");
	times(buffer);
	*l = (long) buffer->tms_utime;
	my_free(buffer);
	return OK;
#endif
#ifdef SYSTEMMAC
	*l = 0;
	return ERROR;
#endif
#ifdef MSDOS
	*l = (long) clock();
	return OK;
#endif /* MSDOS */
}

#ifndef M_PI
#define M_PI 3.1415927
#endif

#ifdef SYSTEMMAC
#include <math.h>
#endif

#if TEXDOCU
INT error(char *fehlertext)
#endif
{
	return f_error(stdout, fehlertext);
}

#if TEXDOCU
INT f_error(FILE * fp_txt, char *fehlertext)
#endif
/* if answer == a ==> abort
   if answer == g ==> go on
   if answer == f ==> go on forever
   else               exit */
{
	char antwort[2];
	static int forever=0;
	
	fflush(fp_txt);
#ifdef SYSTEMMAC
	INT i, j;
	double x;
	fprintf(stderr,"ERROR: %s\n",fehlertext);
	for (j = 0L; j < 3L; j++) {
		for (i = 0L; i < 180L; i++) {
			x = sin((double)i * M_PI / 180.);
			}
		}
	return ERROR;
#endif
	if (forever==2) return ERROR;
	fflush(stdout);
	fflush(stderr);
	fprintf(fp_txt, "\nenter a to abort with core dump, g to go, f to supress, else stop\n");
	fprintf(fp_txt,"ERROR: %s?: ",fehlertext);
	fflush(fp_txt);
	if (forever==1) return ERROR;
	scanf("%s",antwort);
	if (antwort[0] == 'a') abort();
	if (antwort[0] == 'f') {forever = 1; return ERROR;}
	if (antwort[0] == 's') {forever = 2; return ERROR;}
	if (antwort[0] == 'g') return ERROR;
	/* exit(1); */
	return OK;
}


#if TEXDOCU
INT no_memory()
#endif
{
	return error("no memory left");
}


