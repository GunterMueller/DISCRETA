/* iof.C: (file IO) */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>
#ifdef DIVS_TRUE
#include <DISCRETA/divs.h>
#endif
#ifdef MATRIXTRUE
#include <DISCRETA/ma.h>
#endif
#ifdef PERMTRUE
#include <DISCRETA/perm.h>
#endif
#ifdef BRUCHTRUE
#ifndef BRUCH_INCLUDED
#include <DISCRETA/bruch.h>
#endif
#endif
#ifdef LONGINTTRUE
#ifndef LO_INCLUDED
#include <DISCRETA/lo.h>
#endif
#endif
#include <ctype.h> // for isalpha etc.

#define USE_COMPRESS


#define ONE_BYTE_INT(a) (((a) > -126) && ((a) < 127))

#if TEXDOCU
INT SYM_OB::save(BYTE *fname)
#else
saves the object this into a binary file with name fname.
This binary file is a so-called DISCRETA file 
and is only readable by the corresponding 
DISCRETA routines (namely, load, for instance).
It is recommended that DISCRETA files should end 
in \lq .dsc\rq in their file names. 
But this is not a must for these routines.
The size of the generated file is written on the screen.
#endif
{
		return write_op_file(this, fname, 
			TRUE /* f_verbose */, FALSE /* f_use_compress */);
}

#if TEXDOCU
INT SYM_OB::save_quiet(BYTE *fname)
#else
saves the object whithout noise, i.e. 
no size of file is written on the screen.
#endif
{
		return write_op_file(this, fname, 
			FALSE /* f_verbose */, FALSE /* f_use_compress */);
}

#if TEXDOCU
INT SYM_OB::load(BYTE *fname)
#else
loads a DISCRETA binary file into the this object.
Pendant to save()
#endif
{
		return read_op_file(this, fname, 
			TRUE /* f_verbose */, FALSE /* f_use_compress */);
}

#if TEXDOCU
INT SYM_OB::load_quiet(BYTE *fname)
#else
loads without noise
#endif
{
		return read_op_file(this, fname, 
			FALSE /* f_verbose */, FALSE /* f_use_compress */);
}

#if TEXDOCU
INT write_vec_file(VECTOR_OP V, BYTE *file_name, 
	INT f_verbose, INT f_use_compress)
#else
#endif
{
	return write_op_file(V, file_name, f_verbose, f_use_compress);
}

#if TEXDOCU
INT read_vec_file(VECTOR_OP V, BYTE *file_name, 
	INT f_verbose, INT f_use_compress)
#else
#endif
{
	INT erg = OK, len;

	erg += read_op_file(V, file_name, f_verbose, f_use_compress);
	if (f_verbose && erg == OK) {
		len = V->s_li();
		printf("read vector of %ld items.\n", len);
		}
	return erg;
}

#if TEXDOCU
INT write_op_file(SYM_OP V, BYTE *file_name, 
	INT f_verbose, INT f_use_compress)
#else
one level below \lq save\rq: 
here, we can specify if compression 
should be applied or not
#endif
{
	MEM_OP mem = (MEM_OP) callocobject("write_op_file()");
	INT f_very_verbose = FALSE, size, debug_depth;
	BYTE *pc;
	FILE *f;
	
	if (f_very_verbose)
		debug_depth = 1;
	else
		debug_depth = 0;
	V->pack(mem, f_verbose, debug_depth);
	if (f_use_compress)
		mem->compress(f_verbose);
	size = mem->s_used_length_i();
	pc = mem->ob_self.ob_charpointer;
	
	f = fopen(file_name, "wb");

	fwrite(pc, 1 /* size */, size /* items */, f);
	
	fclose(f);
	if (file_size(file_name) != size) {
		printf("write_op_file(): file %s\n", file_name);
		fflush(stdout);
		return error("write_op_file(): file_size(file_name) != size");
		}
	if (f_verbose) {
		printf("write_op_file()|wrote file %s of size %ld\n", 
			file_name, file_size(file_name));
		fflush(stdout);
		}
	
	freeall(mem);
	return OK;
}

#if TEXDOCU
INT read_op_file(SYM_OP V, BYTE *file_name, 
	INT f_verbose, INT f_use_compress)
#else
one level below \lq load\rq
#endif
{
	MEM_OP mem = (MEM_OP) callocobject("read_op_file");
	INT f_very_verbose = FALSE, size, debug_depth;
	BYTE *pc;
	FILE *f;
	
	size = file_size(file_name);
	if (f_verbose) {
		printf("read_op_file(): reading file %s of size %ld\n", 
			file_name, size);
		}
	mem->alloc(size);
	pc = mem->ob_self.ob_charpointer;

	f = fopen(file_name, "rb");

	fread(pc, 1 /* size */, size /* items */, f);
	
	fclose(f);
	if (f_verbose) {
		printf("file read.\n");
		}
	mem->c_used_length(size);
	if (f_use_compress)
		mem->decompress(TRUE /* f_verbose */);
	mem->c_cur_pointer(0);
	if (f_very_verbose)
		debug_depth = 1;
	else
		debug_depth = 0;
	V->unpack(mem, f_verbose, debug_depth);
	
	freeall(mem);
	return OK;
}

#if TEXDOCU
INT SYM_OB::save_ascii(FILE *fp)
#else
writes in ASCII text format (uuencoded like) 
into the stream fp. 
#endif
{
	MEM_OP mem = (MEM_OP) callocobject("SYM::save_ascii");
	INT f_verbose = FALSE, f_very_verbose = FALSE;
	INT size, debug_depth;
	INT i;
	UINT a, a1, a2;
	UBYTE *pc, c1, c2;

	if (f_verbose) {
		printf("SYM_OB::save_ascii(): calculating memory size\n");
		fflush(stdout);
		}
	if (f_very_verbose)
		debug_depth = 1;
	else
		debug_depth = 0;
	pack(mem, f_verbose, debug_depth);
#ifdef USE_COMPRESS
	mem->compress(f_verbose);
#endif
	size = mem->s_used_length_i();
	pc = (UBYTE *) mem->ob_self.ob_charpointer;
	
	fprintf(fp, "ASCII %ld\n", size);
	for (i = 0; i < size; i++) {
		a = (UINT) pc[i];
		a1 = a % (UINT) 16;
		a2 = a >> 4;
		c1 = '0' + a1;
		c2 = '0' + a2;
		fprintf(fp, "%c%c", c1, c2);
		if ((i + 1) % 40 == 0)
			fprintf(fp, "\n");
		}
	fprintf(fp, "\nASCIIEND\n");
	fflush(fp);

	freeall(mem);
	return OK;
}

#define BUFSIZE 10000

#if TEXDOCU
INT SYM_OB::load_ascii(FILE *fp)
#else
reads ASCII style objects written with save-ascii
#endif
{
	MEM_OP mem = (MEM_OP) callocobject("SYM::load_ascii");
	BYTE buf[BUFSIZE];
	BYTE str[1024], *p;
	INT f_verbose = TRUE;
	INT f_very_verbose = FALSE, size, i, debug_depth;
	UBYTE *pc;
	UBYTE c;
	INT a;
	UINT a1, a2;
		
	if (fgets(buf, BUFSIZE, fp) == NULL)
		return error("SYM_OB::load_ascii(): error reading header");
	p = buf;
	s_scan_token(&p, str);
	if (strcmp(str, "ASCII") != 0)
		return error("SYM_OB::load_ascii(): error reading header: ASCII keyword not found");
	s_scan_int(&p, &size);
	if (f_verbose) {
		printf("SYM_OB::load_ascii(): reading ASCII file of size %ld\n", size);
		}
	mem->alloc(size);
	pc = (UBYTE *) mem->ob_self.ob_charpointer;
	for (i = 0; i < size; i++) {
		while (TRUE) {
			a = fgetc(fp);
			if (a == '\n')
				continue;
			if (a == EOF) {
				printf("i = %ld\n", i);
				return error("SYM_OB::load_ascii(): EOF occured");
				}
			break;
			}
		a1 = (UINT) a;
		a2 = fgetc(fp);
		if (a2 == EOF)
			return error("SYM_OB::load_ascii(): EOF occured reading a2");
		a1 = a1 - '0';
		a2 = a2 - '0';
		a = a2 << 4;
		a += a1;
		c = (UBYTE) a;
		pc[i] = c; 
		}
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NULL)
			return error("SYM_OB::load_ascii(): error reading footer");
		if (buf[0] != '\n')
			break;
		}
	p = buf;
	s_scan_token(&p, str);
	if (strcmp(str, "ASCIIEND") != 0)
		return error("SYM_OB::load_ascii(): error reading footer: ASCIIEND keyword not found");

			
	
	if (f_verbose) {
		printf("file read.\n");
		}
	mem->c_used_length(size);
#ifdef USE_COMPRESS
	mem->decompress(TRUE /* f_verbose */);
#endif
	mem->c_cur_pointer(0);
	if (f_very_verbose)
		debug_depth = 1;
	else
		debug_depth = 0;
	unpack(mem, f_verbose, debug_depth);
	
	freeall(mem);
	return OK;
}

#if TEXDOCU
INT SYM_OB::pack(MEM_OP mem, INT f_verbose, INT debug_depth)
#else
used to pack (i.e. to linearize) objects into (binary) strings 
in MEM-OB objects.
#endif
{
	INT size, size0;
	
	if (f_verbose) {
		printf("SYM_OB::pack(): calculating memory size\n");
		fflush(stdout);
		}
	size0 = calc_size_on_file();
	mem->init(0, NIL);
	if (f_verbose) {
		printf("SYM::pack(): allocating memory of size %ld\n", size0);
		fflush(stdout);
		}
	if (mem->alloc(size0) != OK)
		return error("SYM::pack(): no memory");
	mem->c_used_length(0);
	write_mem(mem, debug_depth);
	size = mem->s_used_length_i();
	if (size != size0) {
		printf("SYM::pack(): size = %ld != size0 = %ld", size, size0);
		return ERROR;
		}
	return OK;
}

#if TEXDOCU
INT SYM_OB::unpack(MEM_OP mem, INT f_verbose, INT debug_depth)
#else
unpacks an object from a binary representation in a MEM-OB
#endif
{
	read_mem(mem, debug_depth);
	return OK;
}

#if TEXDOCU
INT SYM_OB::calc_size_on_file()
#else
Computes the size of the linearized version of the object.
If the size is known in adnvance, no realloc is needed in pack.
This speeds up pack a lot !
#endif
{
	INT size = 0;
	
	size += 1; /* obj_kind (now a char) */
	switch (s_obj_bk()) {
	case EMPTY:
		break;
	case INTEGER:
		size += 4;
		break;
#ifdef DIVS_TRUE
	case STRING:
		size += strlen(((STRING_OP)this)->s_str()) + 1;
		break;
#endif
#ifdef BRUCHTRUE
	case BRUCH:
		size += ((BRUCH_OP)this)->s_o()->calc_size_on_file();
		size += ((BRUCH_OP)this)->s_u()->calc_size_on_file();
		break;
#endif
#ifdef VECTORTRUE
	case VECTOR:
		size += ((VECTOR_OP)this)->calc_size_on_file();
		break;
#endif
	case UNIPOLY:
		size += ((VECTOR_OP)this)->calc_size_on_file();
		break;
#ifdef MATRIXTRUE
	case MATRIX:
		size += ((MATRIX_OP)this)->calc_size_on_file();
		break;
#endif
#ifdef PERMTRUE
	case PERMUTATION:
		size += ((PERMUTATION_OP)this)->calc_size_on_file();
		break;
#endif
#ifdef LONGINTTRUE
	case LONGINT:
		size += ((LONGINT_OP)this)->calc_size_on_file();
		break;
#endif
	case MEM:
		size += ((MEM_OP)this)->calc_size_on_file();
		break;
	default:
		printobjectkind();
		return error("SYM::calc_size_on_file(): nyi");
		
	}
	return size;
}

#if TEXDOCU
INT SYM_OB::write_mem(MEM_OP mem, INT debug_depth)
#endif
{
	INT erg = OK, k;
	
	k = s_obj_k();
	if (!ONE_BYTE_INT(k))
		return error("write_mem(): obj_kind not 1 byte");
	erg += mem->write_char((BYTE) k);
	switch (s_obj_bk()) {
	case EMPTY:
		break;
	case INTEGER:
		erg += mem->write_int(((INTEGER_OP)this)->s_i());
		break;
#ifdef DIVS_TRUE
	case STRING:
		erg += ((STRING_OP)this)->write_mem(mem);
		break;
#endif
#ifdef BRUCHTRUE
	case BRUCH:
		((BRUCH_OP)this)->s_o()->write_mem(mem, debug_depth);
		((BRUCH_OP)this)->s_u()->write_mem(mem, debug_depth);
		break;
#endif
#ifdef VECTORTRUE
	case VECTOR:
		((VECTOR_OP)this)->write_mem(mem, debug_depth);
		break;
#endif
	case UNIPOLY:
		((VECTOR_OP)this)->write_mem(mem, debug_depth);
		break;
#ifdef MATRIXTRUE
	case MATRIX:
		((MATRIX_OP)this)->write_mem(mem, debug_depth);
		break;
#endif
#ifdef PERMTRUE
	case PERMUTATION:
		((PERMUTATION_OP)this)->write_mem(mem, debug_depth);
		break;
#endif
#ifdef LONGINTTRUE
	case LONGINT:
		((LONGINT_OP)this)->write_mem(mem, debug_depth);
		break;
#endif
	case MEM:
		((MEM_OP)this)->write_mem(mem);
		break;
	default:
		printobjectkind();
		return error("SYM::write_mem(): nyi");
		
	}
	return erg;
}

#if TEXDOCU
INT SYM_OB::read_mem(MEM_OP mem, INT debug_depth)
#endif
{
	INT erg = OK;
	INT obj_kind;
	INT i;
	BYTE c;
	
	if (!emptyp())
		freeself();
	if (mem->read_char(&c) != OK) {
		Srff("SYM::read_mem", "read_char");
		return ERROR;
		}
	obj_kind = (INT) c;
	if (obj_kind < 0 || obj_kind >= SYM_MAX_KIND) {
		printf("illegal obj_kind: %ld\n", obj_kind);
		return ERROR;
		}
	// switch (ik[obj_kind].ik_kind) {
	// switch (obj_kind) {
	switch (base_kind[obj_kind]) {
	case EMPTY:
		break;
	case INTEGER:
		erg += mem->read_int(&i);
		((INTEGER_OP)this)->m_i(i);
		break;
#ifdef DIVS_TRUE
	case STRING:
		erg += ((STRING_OP)this)->read_mem(mem);
		break;
#endif
#ifdef BRUCHTRUE
	case BRUCH:
		((BRUCH_OP)this)->b_ou(callocobject("read_mem(): BRUCH"), callocobject("read_mem(): BRUCH"));
		((BRUCH_OP)this)->s_o()->read_mem(mem, debug_depth);
		((BRUCH_OP)this)->s_u()->read_mem(mem, debug_depth);
		break;
#endif
#ifdef VECTORTRUE
	case VECTOR:
		((VECTOR_OP)this)->read_mem(mem, debug_depth);
		break;
#endif
	case UNIPOLY:
		((VECTOR_OP)this)->read_mem(mem, debug_depth);
		c_obj_k(UNIPOLY);
		break;
#ifdef MATRIXTRUE
	case MATRIX:
		((MATRIX_OP)this)->read_mem(mem, debug_depth);
		break;
#endif
#ifdef PERMTRUE
	case PERMUTATION:
		((PERMUTATION_OP)this)->read_mem(mem, debug_depth);
		break;
#endif
#ifdef LONGINTTRUE
	case LONGINT:
		((LONGINT_OP)this)->read_mem(mem, debug_depth);
		break;
#endif
	case MEM:
		((MEM_OP)this)->read_mem(mem);
		break;
	default:
		printf("SYM::read_mem(): nyi obj_kind = ");
		print_kind(obj_kind);
		printf("\n");
		return ERROR;
		
	}
	c_obj_k(obj_kind);
	return erg;
}

/*
 * STRING
 */

#ifdef DIVS_TRUE

INT STRING_OB::write_mem(MEM_OP mem)
{
	INT l;
	INT erg = OK;
	BYTE str[5];
	
#ifdef DEBUG_WRITE
	printf("STRING::write_mem()|%s\n", s_str());
#endif
	if (ob_self.ob_charpointer == NIL) {
		str[0] = 0;
		erg += mem->append(1, str);
		}
	else {
		l = strlen(ob_self.ob_charpointer);
		erg += mem->append(l + 1, ob_self.ob_charpointer);
		}
	return erg;
}

INT STRING_OB::read_mem(MEM_OP mem)
{
	INT cur_p, used_length, l, l1;
	BYTE *cp;
	
	cur_p = mem->s_cur_pointer_i();
	used_length = mem->s_used_length_i();
	if (cur_p > used_length) {
		printf("STRING::read_mem()|cur_p > used_length\n");
		init(" ");
		return OK;
		}
	l1 = used_length - cur_p;
	cp = mem->ob_self.ob_charpointer;
	for (l = 0; l < l1; l++) {
		if (cp[cur_p + l] == 0) {
			if (l == 0) {
				/* printf("STRING::read_mem()|warning: "
				"string of length 0\n"); */
				init("");
				}
			else {
				init(cp + cur_p);
				}
			mem->c_cur_pointer(cur_p + l + 1);
#ifdef DEBUG_READ
			printf("STRING::read_mem():%s\n", s_str());
#endif
			return OK;
			}
		}
	printf("STRING::read_mem()|no null byte\n");
	init(" ");
	return OK;
}
#endif

/*
 * MEM
 */

INT MEM_OB::calc_size_on_file()
{
	INT size = 0;
	INT l;
	
	l = s_used_length_i();
	size += 4; /* l */
	size += l; /* data */
	return size;
}

INT MEM_OB::write_mem(MEM_OP mem)
{
	INT erg = OK;
	INT l;
	BYTE *pc;
	
	l = s_used_length_i();
	pc = ob_self.ob_charpointer;
	erg += mem->write_int(l);
	if (mem->append(l, pc) != OK)
		return error("MEM::write_mem(): error in mem->append()");
	return erg;
}

INT MEM_OB::read_mem(MEM_OP mem)
{
	INT erg = OK;
	INT l, cur_pointer;
	BYTE *pc, *data;
	
	erg += mem->read_int(&l);
	cur_pointer = mem->s_cur_pointer_i();
	pc = mem->ob_self.ob_charpointer;
	data = pc + cur_pointer;
	if (init(l, data) != OK)
		return error("MEM::read_mem(): error in init()");
	mem->c_cur_pointer(cur_pointer + l);
	return erg;
}

/*
 * VEC
 */

#ifdef VECTORTRUE

INT VECTOR_OB::hip()
/* homogeneous integer vector predicate */
{
	INT l, len;
	
	len = s_li();
	for (l = 0; l < len; l++) {
		if (s_i(l)->s_obj_k() != INTEGER)
			return FALSE;
		}
	return TRUE;
}

INT VECTOR_OB::hip1()
/* homogeneous integer vector predicate, 
 * test for 1 byte numbers; 
 * only to apply if hip TRUE. */
{
	INT l, len, k;
	
	len = s_li();
	for (l = 0; l < len; l++) {
		if (s_i(l)->s_obj_k() != INTEGER)
			return error("VECTOR::hip1(): not INTEGER");
		k = s_ii(l);
		if (!ONE_BYTE_INT(k))
			return FALSE;
		}
	return TRUE;
}

INT VECTOR_OB::calc_size_on_file()
{
	INT l, len;
	BYTE f_hip, f_hip1;
	INT size = 0;
	
	len = s_li();
	size += 4; /* len */
	f_hip = (BYTE) hip();
	size += 1; /* f_hip */
	if (f_hip) {
		f_hip1 = (BYTE) hip1();
		size += 1; /* f_hip1 */
		if (f_hip1)
			size += 1 * len;
		else
			size += 4 * len;
		}
	else {
		for (l = 0; l < len; l++)
			size += s_i(l)->calc_size_on_file();
		}
	return size;
}

INT VECTOR_OB::write_mem(MEM_OP mem, INT debug_depth)
{
	INT l, len, k;
	BYTE f_hip, f_hip1;
	INT erg = OK;
	
	len = s_li();
	erg += mem->write_int(len);
	f_hip = (BYTE) hip();
	if (f_hip)
		f_hip1 = (BYTE) hip1();
	if (debug_depth > 0) {
		printf("writing ");
		if (f_hip) {
			if (f_hip1)
				printf("hip1 ");
			else
				printf("hip ");
			}
		printf("vector of length %ld\n", len);
		}
	erg += mem->write_char(f_hip);
	if (f_hip) {
		erg += mem->write_char(f_hip1);
		if (f_hip1) {
			for (l = 0; l < len; l++) {
				k = s_ii(l);
				erg += mem->write_char((BYTE) k);
				}
			}
		else {
			for (l = 0; l < len; l++) {
				erg += mem->write_int(s_ii(l));
				}
			}
		}
	else {
		for (l = 0; l < len; l++) {
			if (debug_depth > 0) {
				printf("%ld ", l);
				if ((l % 20) == 0)
					printf("\n");
				fflush(stdout);
				}
			if (s_i(l)->write_mem(mem, debug_depth - 1) != OK)
				return error("VECTOR::write_mem()|"
					"error in SYM::write_mem()");
			}
		}
	return erg;
}

INT VECTOR_OB::read_mem(MEM_OP mem, INT debug_depth)
{
	INT l, len, k;
	BYTE c, f_hip, f_hip1;
	INT erg = OK;
	
	if (mem->read_int(&len) != OK) {
		Srff("VECTOR::read_mem", "mem->read_int");
		return ERROR;
		}
	erg += m_il(len);
	if (mem->read_char(&f_hip) != OK) {
		Srff("VECTOR::read_mem", "mem->read_char");
		return ERROR;
		}
	if (f_hip) {
		if (mem->read_char(&f_hip1) != OK) {
			Srff("VECTOR::read_mem", "mem->read_char");
			return ERROR;
			}
		}
	if (debug_depth > 0) {
		printf("reading ");
		if (f_hip) {
			if (f_hip1)
				printf("hip1 ");
			else
				printf("hip ");
			}
		printf("vector of length %ld\n", len);
		}
	if (f_hip) {
		if (f_hip1) {
			for (l = 0; l < len; l++) {
				erg += mem->read_char(&c);
				k = (INT) c;
				m_ii(l, k);
				}
			}
		else {
			for (l = 0; l < len; l++) {
				if (mem->read_int(&k) != OK) {
					Srff("VECTOR::read_mem", "mem->read_int");
					return ERROR;
					}
				m_ii(l, k);
				}
			}
		}
	else {
		for (l = 0; l < len; l++) {
			if (debug_depth > 0) {
				printf("%ld ", l);
				if ((l % 20) == 0)
					printf("\n");
				fflush(stdout);
				}
			if (s_i(l)->read_mem(mem, debug_depth - 1) != OK)
				return error("VECTOR::read(): error in SYM::read_mem()");
			}
		}
	return erg;
}
#endif

/*
 * MA
 */

#ifdef MATRIXTRUE

INT MATRIX_OB::hip()
/* homogeneous integer matrix predicate */
{
	INT i, j, l, h;
	
	l = s_li();
	h = s_hi();
	for (i = 0; i < h; i++)
		for (j = 0; j < l; j++)
			if (s_ij(i, j)->s_obj_k() != INTEGER)
				return FALSE;
	return TRUE;
}

INT MATRIX_OB::hip1()
/* homogeneous integer matrix predicate, 
 * test for 1 byte numbers; 
 * only to apply if hip TRUE. */
{
	INT i, j, l, h, k;
	
	l = s_li();
	h = s_hi();
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			if (s_ij(i, j)->s_obj_k() != INTEGER)
				return error("MATRIX::hip1(): not INTEGER");
			k = s_iji(i, j);
			if (!ONE_BYTE_INT(k))
				return FALSE;
			}
		}
	return TRUE;
}

INT MATRIX_OB::calc_size_on_file()
{
	INT size = 0;
	INT i, j, l, h;
	BYTE f_hip, f_hip1;
	
	l = s_li();
	h = s_hi();
	size += 8; /* l, h */
	f_hip = (BYTE) hip();
	size += 1; /* f_hip */
	if (f_hip) {
		f_hip1 = (BYTE) hip1();
		size += 1; /* f_hip1 */
		if (f_hip1)
			size += 1 * l * h;
		else
			size += 4 * l * h;
		}
	else {
		for (i = 0; i < h; i++) {
			for (j = 0; j < l; j++) {
				size += s_ij(i, j)->
					calc_size_on_file();
				}
			}
		}
	return size;
}

INT MATRIX_OB::write_mem(MEM_OP mem, INT debug_depth)
{
	INT erg = OK;
	INT i, j, l, h, k;
	BYTE f_hip, f_hip1;
	
	l = s_li();
	h = s_hi();
	erg += mem->write_int(l);
	erg += mem->write_int(h);
	f_hip = (BYTE) hip();
	if (f_hip)
		f_hip1 = (BYTE) hip1();
	if (debug_depth > 0) {
		printf("writing %ld x %ld ", h, l);
		if (f_hip) {
			if (f_hip1)
				printf("hip1 ");
			else
				printf("hip ");
			}
		printf("matrix\n");
		}
	erg += mem->write_char(f_hip);
	if (f_hip) {
		erg += mem->write_char(f_hip1);
		if (f_hip1) {
			for (i = 0; i < h; i++) {
				for (j = 0; j < l; j++) {
					k = s_iji(i, j);
					erg += mem->write_char((BYTE) k);
					}
				}
			}
		else {
			for (i = 0; i < h; i++) {
				for (j = 0; j < l; j++) {
					erg += mem->write_int(s_iji(i, j));
					}
				}
			}
		}
	else {
		for (i = 0; i < h; i++) {
			for (j = 0; j < l; j++) {
#if 0
				if (debug_depth > 0) {
					printf("(%ld,%ld) ", i, j);
					if ((j % 20) == 0)
						printf("\n");
					fflush(stdout);
					}
#endif
				if (s_ij(i, j)->write_mem(mem, debug_depth - 1) != OK)
					return error("MATRIX::write_mem()|"
						"error in SYM::write_mem()");
				}
			}
		}
	return erg;
}

INT MATRIX_OB::read_mem(MEM_OP mem, INT debug_depth)
{
	INT erg = OK;
	INT i, j, l, h, k;
	BYTE c, f_hip, f_hip1;
	
	erg += mem->read_int(&l);
	erg += mem->read_int(&h);
	erg += m_ilih(l, h);
	if (mem->read_char(&f_hip) != OK) {
		Srff("MATRIX::read_mem", "mem->read_char");
		return ERROR;
		}
	if (f_hip) {
		if (mem->read_char(&f_hip1) != OK) {
			Srff("MATRIX::read_mem", "mem->read_char");
			return ERROR;
			}
		}
	if (debug_depth > 0) {
		printf("reading %ld x %ld ", h, l);
		if (f_hip) {
			if (f_hip1)
				printf("hip1 ");
			else
				printf("hip ");
			}
		printf("matrix\n");
		}
	if (f_hip) {
		if (f_hip1) {
			for (i = 0; i < h; i++) {
				for (j = 0; j < l; j++) {
					if (mem->read_char(&c) != OK) {
						Srff("MATRIX::read_mem", "mem->read_char");
						return ERROR;
						}
					k = (INT) c;
					m_iji(i, j, k);
					}
				}
			}
		else {
			for (i = 0; i < h; i++) {
				for (j = 0; j < l; j++) {
					if (mem->read_int(&k) != OK) {
						Srff("MATRIX::read_mem", "mem->read_int");
						return ERROR;
						}
					m_iji(i, j, k);
					}
				}
			}
		}
	else {
		for (i = 0; i < h; i++) {
			for (j = 0; j < l; j++) {
				if (debug_depth > 0) {
					printf("(%ld,%ld) ", i, j);
					if ((j % 20) == 0)
						printf("\n");
					fflush(stdout);
					}
				if (s_ij(i, j)->read_mem(mem, debug_depth - 1) != OK)
					return error("MATRIX::read_mem(): error in SYM::read_mem()");
				}
			}
		}
	return erg;
}
#endif

/*
 * PERM 
 */

#ifdef PERMTRUE

INT PERMUTATION_OB::calc_size_on_file()
{
	INT len;
	INT size = 0;
	
	len = s_li();
	size += 4; /* len */
	if (ONE_BYTE_INT(len))
		size += len * 1;
	else
		size += len * 4;
	return size;
}

INT PERMUTATION_OB::write_mem(MEM_OP mem, INT debug_depth)
{
	INT l, len, i;
	INT erg = OK;
	
	len = s_li();
	erg += mem->write_int(len);
	if (ONE_BYTE_INT(len)) {
		for (l = 0; l < len; l++) {
			i = s_ii(l);
			erg += mem->write_char((BYTE) i);
			}
		}
	else {
		for (l = 0; l < len; l++) {
			i = s_ii(l);
			erg += mem->write_int(i);
			}
		}
	return erg;
}

INT PERMUTATION_OB::read_mem(MEM_OP mem, INT debug_depth)
{
	INT l, len, i;
	INT erg = OK;
	BYTE c;
	
	erg += mem->read_int(&len);
	m_il(len);
	if (ONE_BYTE_INT(len)) {
		for (l = 0; l < len; l++) {
			erg += mem->read_char(&c);
			i = (INT) c;
			m_ii(l, i);
			}
		}
	else {
		for (l = 0; l < len; l++) {
			erg += mem->read_int(&i);
			m_ii(l, i);
			}
		}
	return erg;
}
#endif

/*
 * LONGINT 
 */

#undef DEBUG_LONGINT_IOF

#ifdef LONGINTTRUE

INT LONGINT_OB::calc_size_on_file()
{
	INT size;
	STRING_OB s;
	
	lo2string(&s);
	size = s.calc_size_on_file();
#ifdef DEBUG_LONGINT_IOF
	printf("LO::calc_size_on_file() ");
	print();
	printf(" size = %ld\n", size);
	fflush(stdout);
#endif
	return size;
}

INT LONGINT_OB::write_mem(MEM_OP mem, INT debug_depth)
{
	STRING_OB s;
	
#ifdef DEBUG_LONGINT_IOF
	printf("LO::write_mem() ");
	print();
	fflush(stdout);
#endif
	lo2string(&s);
#ifdef DEBUG_LONGINT_IOF
	printf(" as string: ");
	s.println();
	fflush(stdout);
#endif
	s.write_mem(mem);
	return OK;
}

INT LONGINT_OB::read_mem(MEM_OP mem, INT debug_depth)
{
	STRING_OB s;
	
	s.read_mem(mem);
#ifdef DEBUG_LONGINT_IOF
	printf("LO::read_mem() ");
	s.print();
	fflush(stdout);
#endif
	string2lo(&s);
#ifdef DEBUG_LONGINT_IOF
	printf(" as longint: ");
	println();
	fflush(stdout);
#endif
	return OK;
}
#endif

/*
 * various utility functions
 */

#if TEXDOCU
#else
various utility functions:
#endif

#if TEXDOCU
BYTE *Eostr(BYTE *s)
#endif
{
	return(&s[strlen(s)]);
}

#if TEXDOCU
void Srfs(BYTE *s1, BYTE *s2)
#endif
{
	printf("%s()|%s\n", s1, s2);
}

#if TEXDOCU
void Srff(BYTE *s1, BYTE *s2)
#endif
{
	printf("%s()|error in %s()\n", s1, s2);
}

#if TEXDOCU
INT tstropn(TSTRING **tp, BYTE *p)
#endif
{
	BYTE *p1 = NIL;
	INT l;
	
	if (tp == NIL || p == NIL) {
		Srfs("tstropn", "args NIL");
		return(FALSE);
		}
#if FALSE
	if (*tp != NIL) {
		Srfs("tstropn", "Warning: Overwriting existing tstring");
		}
#endif
	tstrcls(tp);
		/* Neu: Vorher wurde *tp = NIL gesetzt !!! */
	l = strlen(p);
	if (l > 10000) {
		Srfs("tstropn", "Warning: strlen(p) > 10000|no tstropn");
		return(TRUE);
		}
	if (l == 0L) {
		*tp = NIL;
		return(TRUE);
		}
	p1 = (BYTE *)my_malloc(l + 1, "tstropn()");
	if (p1 == NIL) {
		Srfs("tstropn", "no memory");
		return(FALSE);
		}
	strcpy(p1, p);
	*tp = p1;
	return(TRUE);
}

#if TEXDOCU
INT tstrcls(TSTRING **tp)
#endif
{
	if (tp == NIL) {
		Srfs("tstrcls", "args NIL");
		return(FALSE);
		}
	if (*tp != NIL) {
		my_free(*tp);
		*tp = NIL;
		}
	return(TRUE);
}

#if TEXDOCU
BYTE *tstr(TSTRING *tp)
#endif
{
	if (tp == NIL)
		return("");
	else
		return(tp);
}

#if TEXDOCU
INT tstrlen(TSTRING *tp)
#endif
{
	INT l;
	
	if (tp == NIL)
		return(0L);
	l = strlen(tp);
	return(l);
}

#if TEXDOCU
INT tstrfill(TSTRING **tp, INT len, INT c)
#endif
{
	BYTE *p1 = NIL;
	
	if (tp == NIL) {
		Srfs("tstrfill", "args NIL");
		return(FALSE);
		}
	if (len == 0L) {
		tstrcls(tp);
		return(TRUE);
		}
	p1 = (BYTE *)my_malloc(len + 1, "tstrfill");
	if (p1 == NIL) {
		Srfs("tstrfill", "no memory");
		return(FALSE);
		}
	strcpy(p1, tstr(*tp));
	strfill(p1, len, c);
	tstrcls(tp);
	*tp = p1;
	return(TRUE);
}

#if TEXDOCU
INT tstrcat(TSTRING **tp, BYTE *p)
#endif
{
	BYTE *p1 = NIL;
	INT len, len1, len2;
	
	if (tp == NIL || p == NIL) {
		Srfs("tstrcat", "args NIL");
		return(FALSE);
		}
	len1 = tstrlen(*tp);
	if (len1 > 2000) {
		Srfs("tstrcat", 
			"Warning: len1 > 2000|no tstrcat");
		return(TRUE);
		}
	len2 = strlen(p);
	if (len2 > 2000) {
		Srfs("tstrcat", 
			"Warning: len2 > 2000|no tstrcat");
		return(TRUE);
		}
	len = len1 + len2;
	if (len == 0L) {
		tstrcls(tp);
		return(TRUE);
		}
	p1 = (BYTE *)my_malloc(len + 1L, "tstrcat");
	if (p1 == NIL) {
		Srfs("tstrcat", "no memory");
		return(FALSE);
		}
	strcpy(p1, tstr(*tp));
	strcat(p1, p);
	tstrcls(tp);
	*tp = p1;
	return(TRUE);
}

#if TEXDOCU
INT tstrcpy(TSTRING **dst, TSTRING *src)
#endif
/* Liefert NIL-Pointer 
 * als Leerstringdarstellung. */
{
#if FALSE
	if (*dst)
		tstrcls(dst);
	/* Unnoetig, da tstropn() 
	 * jetzt vorher freigibt ! */
#endif
	if (!tstropn(dst, tstr(src))) {
		Srfs("tstrcat", "tstropn");
		return(FALSE);
		}
	return(TRUE);
}

#if TEXDOCU
BYTE *eostr(BYTE *s)
#endif
{
	return(&s[strlen(s)]);
}

#if TEXDOCU
BYTE upperchar(INT c)
#endif
{
	BYTE c1 = (BYTE)c;
	
	if (c1 >= 'a' && c1 <= 'z')
		return(c1 - 'a' + 'A');
	return(c1);
}

#if TEXDOCU
BYTE lowerchar(INT c)
#endif
{
	BYTE c1 = (BYTE)c;
	
	if (c1 >= 'A' && c1 <= 'Z')
		return(c1 - 'A' + 'a');
	return(c1);
}

#if TEXDOCU
void strfill(BYTE *s, INT len, INT c)
#endif
{
	INT l;
	BYTE c1 = (BYTE)c;
	
	l = strlen(s);
	while (l < len)
		s[l++] = c1;
	s[l] = '\0';
}

#if TEXDOCU
void fill_char(void *v, INT cnt, INT c)
#endif
{
	BYTE ch = (BYTE)c;
	BYTE *s;
	
	s = (BYTE *)v;
	cnt++;
	while (--cnt)
		*s++ = ch;
}

#if TEXDOCU
void fill_int(void *v, INT cnt, INT i)
#endif
{
	INT *s;
	
	s = (INT *)v;
	cnt++;
	while (--cnt)
		*s++ = i;
}

#if TEXDOCU
void fill_long(void *v, INT cnt, INT l)
#endif
{
	INT *s;
	
	s = (INT *)v;
	cnt++;
	while (--cnt)
		*s++ = l;
}

#if TEXDOCU
void move_char(void *s, void *d, INT cnt)
#endif
{
	BYTE *ps;
	BYTE *pd;
	
	if (d < s) {
		ps = (BYTE *)s;
		pd = (BYTE *)d;
		cnt++;
		while (--cnt)
			*pd++ = *ps++;
		}
	else if (d > s){
		ps = (BYTE *)s + cnt;
		pd = (BYTE *)d + cnt;
		cnt++;
		while (--cnt)
			*(--pd) = *(--ps);
		}
	/* else d == s: tue nichts */
}

#if TEXDOCU
void move_int(void *s, void *d, INT cnt)
#endif
{
	INT *ps;
	INT *pd;
	
	if (d <= s) {
		ps = (INT *)s;
		pd = (INT *)d;
		cnt++;
		while (--cnt)
			*pd++ = *ps++;
		}
	else if (d > s){
		ps = (INT *)s + cnt;
		pd = (INT *)d + cnt;
		cnt++;
		while (--cnt)
			*(--pd) = *(--ps);
		}
	/* else d == s: tue nichts */
}

#if TEXDOCU
void move_long(void *s, void *d, INT cnt)
#endif
{
	INT *ps;
	INT *pd;
	
	if (d <= s) {
		ps = (INT *)s;
		pd = (INT *)d;
		cnt++;
		while (--cnt)
			*pd++ = *ps++;
		}
	else if (d > s){
		ps = (INT *)s + cnt;
		pd = (INT *)d + cnt;
		cnt++;
		while (--cnt)
			*(--pd) = *(--ps);
		}
	/* else d == s: tue nichts */
}

#if TEXDOCU
INT ij2k(INT i, INT j, INT n)
#endif
{
	if (i == j)
		return error("ij2k() i == j");
	if (i > j)
		return ij2k(j, i, n);
	return ((n - i) * i + 
		((i * (i - 1)) >> 1) + 
		j - i - 1);
}


#if TEXDOCU
INT k2ij(INT k, INT *i, INT *j, INT n)
#endif
{
	INT ii;
	
	for (ii = 0; ii < n; ii++) {
		if (k < n - ii - 1) {
			*i = ii;
			*j = k + ii + 1;
			return OK;
			}
		k -= (n - ii - 1);
		}
	return error("k too large");
}

#if TEXDOCU
INT graph_add_edge(UINT *g, INT i, INT j, INT n)
#endif
{
	INT k;
	
	k = ij2k(i, j, n);
	*g |= ((UINT) 1L << k);
	return OK;
}

#if TEXDOCU
INT graph_del_edge(UINT *g, INT i, INT j, INT n)
#endif
{
	INT k = ij2k(i, j, n);
	
	*g &= ~((UINT) 1L << k);
	return OK;
}

#if TEXDOCU
INT graph_has_edge(UINT g, INT i, INT j, INT n)
#endif
{
	INT k = ij2k(i, j, n);
	
	return ( (g & ((UINT) 1L << k)) ? TRUE : FALSE );
}

#if TEXDOCU
INT graph_code(VECTOR_OP V, UINT *g, INT n)
#endif
{
	INT l, len;
	VECTOR_OP V1;
	
	len = V->s_li();
	*g = 0;
	for (l = 0; l < len; l++) {
		V1 = (VECTOR_OP) V->s_i(l);
		graph_add_edge(g, V1->s_ii(0), V1->s_ii(1), n);
		}
	return OK;
}

#if TEXDOCU
INT graph_decode(VECTOR_OP V, UINT g, INT n)
#endif
{
	VECTOR_OB V1;
	UINT g1 = g;
	INT k = 0, i, j;
	
	V->m_il(0);
	while (g1) {
		if (g1 & (UINT) 1L) {
			V1.m_il(2);
			k2ij(k, &i, &j, n);
			V1.m_ii(0, i);
			V1.m_ii(1, j);
			V->inc();
			V1.swap(V->s_i(V->s_li() - 1));
			}
		g1 >>= 1;
		k++;
		}
	return OK;
}

#if TEXDOCU
INT my_ggt(INT m, INT n)
#endif
{
	INT r, s;
	
	if (n > m) {
		r = m;
		m = n;
		n = r;
		}
	while (TRUE) {
		s = m / n;
		r = m - (s * n);
		if (r == 0)
			return(n);
		m = n;
		n = r;
		}
}

#if TEXDOCU
void my_memcpy(void *dst, void *src, INT size)
#endif
{
	memcpy(dst, src, size);
}

#if TEXDOCU
INT s_scan_int(BYTE **s, INT *i)
#endif
{
	BYTE str1[512];
	
	if (!s_scan_token(s, str1))
		return FALSE;
	sscanf(str1, "%ld", i);
	return TRUE;
}

#if TEXDOCU
INT s_scan_token(BYTE **s, BYTE *str)
#endif
{
	BYTE c;
	INT len;
	
	while (TRUE) {
		c = **s;
		if (c == 0) {
			return(FALSE);
			}
		if (c == ' ' || c == '\t' || 
			c == '\r' || c == 10 || c == 13) {
			(*s)++;
			continue;
			}
		break;
		}
	len = 0;
	c = **s;
	if (isalpha(c)) {
		while (isalnum(c) || c == '_') {
			str[len] = c;
			len++;
			(*s)++;
			c = **s;
			}
		str[len] = 0;
		}
	else if (isdigit(c) || c == '-') {
		str[len++] = c;
		(*s)++;
		c = **s;
		while (isdigit(c)) {
			str[len] = c;
			len++;
			(*s)++;
			c = **s;
			}
		str[len] = 0;
		}
	else {
		str[0] = c;
		str[1] = 0;
		(*s)++;		
		}
	// printf("token = \"%s\"\n", str);
	return TRUE;
}

#if TEXDOCU
INT s_scan_str(BYTE **s, BYTE *str)
#endif
{
	BYTE c;
	INT len, f_break;
	
	while (TRUE) {
		c = **s;
		if (c == 0) {
			return(FALSE);
			}
		if (c == ' ' || c == '\t' || 
			c == '\r' || c == 10 || c == 13) {
			(*s)++;
			continue;
			}
		break;
		}
	if (c != '\"') {
		Srfs("s_scan_str", "c != '\"'");
		return(FALSE);
		}
	(*s)++;
	len = 0;
	f_break = FALSE;
	while (TRUE) {
		c = **s;
		if (c == 0) {
			break;
			}
		if (c == '\\') {
			(*s)++;
			c = **s;
			str[len] = c;
			len++;
			}
		else if (c == '\"') {
			f_break = TRUE;
			}
		else {
			str[len] = c;
			len++;
			}
		(*s)++;
		if (f_break)
			break;
		}
	str[len] = 0;
	return TRUE;
}

BYTE *gl_hot_param; /* the command parameters */

#if TEXDOCU
INT parse_args(VECTOR_OP args)
#endif
{
	BYTE str[10000];
	BYTE str1[10000];
	BYTE *p;
	INT i, ii, l, al = 0;
	BYTE c;
	STRING_OB s_ob;
	INTEGER_OB int_ob;

	p = gl_hot_param;
	args->m_il(0);
	if (!s_scan_token(&p, str))
		return OK;
	l = strlen(str);
	for (i = 0; i < l; i++) {
		c = str[i];
		if (c == 's') {
			if (!s_scan_str(&p, str1))
				return OK;
			printf("parsing string arg %s\n", str1);
			fflush(stdout);
			s_ob.init(str1);
			args->inc();
			s_ob.swap((STRING_OP) args->s_i(al));
			al++;
			}
		else if (c == 'i') {
			if (!s_scan_int(&p, &ii))
				return OK;
			printf("parsing arg %ld\n", ii);
			fflush(stdout);
			int_ob.m_i(ii);
			args->inc();
			int_ob.swap((INTEGER_OP) args->s_i(al));
			al++;
			}
		else
			printf("illegal format parameter %c, ignored\n", c);
		}
	return OK;
}

#if TEXDOCU
void test_lo_iof()
#endif
{
	printf("testing i/o for LONGINT:\n");
	fflush(stdout);
	{
		INT i, k = 20;
		SYM_OB a, b, c;
		/* LONGINT_OB a, b, c; */
		BYTE *fname = "lo_test.dsc";
		
		printf("%ld! = \n", k);
		fflush(stdout);
		/* a.m_i(1); */
		((INTEGER_OP) &a)->m_i(1);
		for (i = 2; i <= k; i++) {
			/* b.m_i(i); */
			((INTEGER_OP) &b)->m_i(i);
			b.mult_apply(&a);
			}
		a.println();
		fflush(stdout);
		write_op_file(&a, fname, TRUE, FALSE /* f_use_compress */);
		printf("written file %s of size %ld\n", fname, file_size(fname));
		fflush(stdout);
		read_op_file(&c, fname, TRUE, FALSE /* f_use_compress */);
		printf("read file %s of size %ld\n", fname, file_size(fname));
		fflush(stdout);
		c.println();
		fflush(stdout);
	}
}


