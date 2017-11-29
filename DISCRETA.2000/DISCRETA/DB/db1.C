/* db1.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef DB_TRUE

#include <DISCRETA/db.h>

INT i_custom_db(DATABASE_OP *db_fg)
{
	DATABASE_OP db;
	BAYERTREE_OP bt;
	
	db = (DATABASE_OP) callocobject("db1.C i_custom_db()");
	db->init("custom.db", VECTOR);

	db->add_btree("custom0.idx", 
		TRUE /* f_duplicatekeys */, 
		0 /* btree_idx */ );
	bt = db->s_btree_i(0);
	bt->add_key_STRING(20, 0 /* field1: */, 0);

	db->add_btree("custom1.idx", 
		TRUE /* f_duplicatekeys */, 
		0 /* btree_idx */ );
	bt = db->s_btree_i(1);
	bt->add_key_STRING(20, 1 /* field1: */, 0);

	db->add_btree("custom2.idx", 
		TRUE /* f_duplicatekeys */, 
		0 /* btree_idx */ );
	bt = db->s_btree_i(2);
	bt->add_key_STRING(20, 2 /* field1: */, 0);

	db->add_btree("custom3.idx", 
		TRUE /* f_duplicatekeys */, 
		0 /* btree_idx */ );
	bt = db->s_btree_i(3);
	bt->add_key_STRING(20, 3 /* field1: */, 0);

	*db_fg = db;
	return OK;
}

void e_custom_db(DATABASE_OP db)
{
	freeall(db);
}

INT c_db_custom_create_and_close()
{
	DATABASE_OP db;

	if (i_custom_db(&db) != OK) {
		Srfs("c_db_custom_create_and_close", 
			"i_custom_db");
		return ERROR;
		}
	if (db->create() != OK) {
		Srfs("c_db_custom_create_and_close", 
			"DB::create");
		return ERROR;
		}
	if (db->close() != OK) {
		Srfs("c_db_custom_create_and_close", 
			"DB::close");
		return ERROR;
		}
	e_custom_db(db);
	return OK;
}

INT c_db_custom_info()
{
	DATABASE_OP db;
	
	if (i_custom_db(&db) != OK) {
		Srfs("c_db_custom_info", "i_custom_db");
		return ERROR;
		}
	db->print(0);
	e_custom_db(db);
	return OK;
}

#define BUFSIZE 5000

#define CDB_DEBUG 1

/* file format:
 * ASCII text, 
 * lines beginning with 'comment_mark' 
 *   are ignored throughout
 * \'begin_mark' starts a dataset
 * \'end_mark' ends a dataset (appends to database)
 * \0key or \1key or \2key or \3key are keys
 */

INT db_custom_import(BYTE *fname, 
	BYTE begin_mark, BYTE end_mark, 
	BYTE eof_mark, BYTE comment_mark)
{
	BYTE buf[BUFSIZE];
	BYTE c;
	FILE *f;
	INT i, f_read_mode = FALSE;
	VECTOR_OB V;
	STRING_OP s;
	DATABASE_OP db;
	
	if (i_custom_db(&db) != OK) {
		Srfs("db_custom_import", "i_custom_db");
		return ERROR;
		}
	if (db->open() != OK) {
		Srfs("db_custom_import", "DB::open");
		return ERROR;
		}
	f = fopen(fname, "r");
	while (TRUE) {
		buf[0] = 0;
		if (fgets(buf, BUFSIZE, f) == NULL)
			break;
		/* discard trailing '\n': */
		if (strlen(buf) > 1)
			buf[strlen(buf) - 1] = 0;
#ifdef CDB_DEBUG
		printf("%s\n", buf);
#endif
		if (strlen(buf) > 1) {
			c = buf[0];
			if (c == '\\') {
				c = buf[1];
				if (c == eof_mark) {
					if (f_read_mode) {
						printf("warning: eof "
						"mark occured in read mode !\n");
						break;
						}
					break;
					}
				else if (c == begin_mark) {
					if (f_read_mode) {
						printf("begin mark "
						"occured in read mode !\n");
						break;
						}
					V.m_il(4);
					s = (STRING_OP) V.s_i(0);
					s->init("");
					s = (STRING_OP) V.s_i(1);
					s->init("");
					s = (STRING_OP) V.s_i(2);
					s->init("");
					s = (STRING_OP) V.s_i(3);
					s->init("");
					f_read_mode = TRUE;
					continue;
					}
				else if (c == end_mark) {
					if (!f_read_mode) {
						printf("end mark occured "
						"not in read mode !\n");
						break;
						}
					if (db->add_op(&V) != OK)
						break;
					f_read_mode = FALSE;
					continue;
					}
				else if (c == '0' || 
					c == '1' || 
					c == '2' || 
					c == '3') {
					if (!f_read_mode) {
						printf("key mark occured "
						"not in read mode !\n");
						break;
						}
					i = c - '0';
					s = (STRING_OP) V.s_i(i);
					s->init(buf + 2);
					continue;
					}
				}
			else if (c == comment_mark)
				continue;
			if (f_read_mode) {
				V.inc();
				s = (STRING_OP) V.s_i(V.s_li() - 1);
				s->init(buf);
				}
			}
		}
	fclose(f);
	if (db->close() != OK) {
		Srfs("db_custom_import", "DB::close");
		return ERROR;
		}
	
	e_custom_db(db);
	return OK;
}

static INT HTML_name(
	DATABASE_OP db, BAYERTREE_OP btree0);
static INT HTML_adress(
	DATABASE_OP db, BAYERTREE_OP btree1);
static INT HTML_country(
	DATABASE_OP db, BAYERTREE_OP btree2);
static INT HTML_seminaire(
	DATABASE_OP db, BAYERTREE_OP btree0);

INT db_custom_HTML()
{
	BAYERTREE_OP btree0; /* name */
	BAYERTREE_OP btree1; /* address */
	BAYERTREE_OP btree2; /* country */
	DATABASE_OP db;
	
	if (i_custom_db(&db) != OK)
		return error("db_custom_HTML(): "
		"error in i_custom_db()");
	btree0 = db->s_btree_i(0);
	btree1 = db->s_btree_i(1);
	btree2 = db->s_btree_i(2);
	/* btree->print_pages(); */

	if (db->open_DB() != OK)
		return error("db_custom_HTML(): "
		"error in db->open_DB()");
	if (btree0->open() != OK ||
		btree1->open() != OK ||
		btree2->open() != OK)
		return error("db_custom_HTML(): "
		"error in btree->open()");

	HTML_seminaire(db, btree0);
	HTML_name(db, btree0);
	HTML_adress(db, btree1);
	HTML_country(db, btree2);

	db->close_DB();
	btree0->close();
	btree1->close();
	btree2->close();
	e_custom_db(db);
	return OK;
}

static INT HTML_name(
	DATABASE_OP db, BAYERTREE_OP btree0)
{
	BYTE *fname = "name.html";
	BYTE label[256];
	FILE *fp;
	KEYTYPE key;
	DATATYPE data;
	BYTE str[1024];
	VECTOR_OB V;
	STRING_OP s0, s1, s2, s;
	INT i, len1;
	BYTE *p_name;

	fp = fopen(fname, "w");
	
	fputs("<Head>\n", fp);
	fputs("<Title>\n", fp);
	fputs("Seminaire, name\n", fp);
	fputs("</Title>\n", fp);
	fputs("</Head>\n", fp);
	
	fputs("<Body>\n", fp);
	fputs("<H1>\n", fp);
	fputs("Seminaire, adress\n", fp);
	fputs("</H1>\n", fp);
	
	btree0->len(&len1);
	fputs("<ul>\n", fp);
	for (i = 0; i < len1; i++) {
		btree0->ith(i, &key, &data);
		db->get_op(&data, &V);
		s0 = (STRING_OP) V.s_i(0);
		s1 = (STRING_OP) V.s_i(1);
		s2 = (STRING_OP) V.s_i(2);
		if (V.s_li() > 4) {
			s = (STRING_OP) V.s_i(4);
			p_name = s->s_str();
			}
		else
			p_name = NIL;
		sprintf(label, "%s_%ld", s0->s_str(), i);
		sprintf(str, "<li> <A HREF="
			"\"seminaire.html#%ld\"> %s, %s </A> ", 
			(INT) data.datref, s0->s_str(), p_name);
		fputs(str, fp);
		fputs("\n", fp);
		fflush(fp);
		}
	fputs("</ul>\n", fp);
	fputs("</Body>\n\n", fp);
	fflush(fp);
	fclose(fp);
	printf("written file %s of size %ld\n", 
		fname, file_size(fname));
	fflush(stdout);
	return OK;
}

static INT HTML_adress(
	DATABASE_OP db, BAYERTREE_OP btree1)
{
	BYTE *fname = "adress.html";
	BYTE label[256];
	FILE *fp;
	KEYTYPE key;
	DATATYPE data;
	BYTE str[1024];
	VECTOR_OB V;
	STRING_OP s0, s1, s2, s;
	INT i, len1;
	BYTE *p_name;

	fp = fopen(fname, "w");
	
	fputs("<Head>\n", fp);
	fputs("<Title>\n", fp);
	fputs("Seminaire, place\n", fp);
	fputs("</Title>\n", fp);
	fputs("</Head>\n", fp);
	
	fputs("<Body>\n", fp);
	fputs("<H1>\n", fp);
	fputs("Seminaire, place\n", fp);
	fputs("</H1>\n", fp);
	
	btree1->len(&len1);
	fputs("<ul>\n", fp);
	for (i = 0; i < len1; i++) {
		btree1->ith(i, &key, &data);
		db->get_op(&data, &V);
		s0 = (STRING_OP) V.s_i(0);
		s1 = (STRING_OP) V.s_i(1);
		s2 = (STRING_OP) V.s_i(2);
		if (V.s_li() > 4) {
			s = (STRING_OP) V.s_i(4);
			p_name = s->s_str();
			}
		else
			p_name = NIL;
		sprintf(label, "%s_%ld", s0->s_str(), i);
		sprintf(str, "<li> "
			"<A HREF=\"seminaire.html#%ld\"> %s, %s </A> ", 
			(INT) data.datref, s1->s_str(), p_name);
		fputs(str, fp);
		fputs("\n", fp);
		fflush(fp);
		}
	fputs("</ul>\n", fp);
	fputs("</Body>\n\n", fp);
	fflush(fp);
	fclose(fp);
	printf("written file %s of size %ld\n", 
		fname, file_size(fname));
	fflush(stdout);
	return OK;
}

static INT HTML_country(
	DATABASE_OP db, BAYERTREE_OP btree2)
{
	BYTE *fname = "country.html";
	BYTE label[256];
	FILE *fp;
	KEYTYPE key;
	DATATYPE data;
	BYTE str[1024];
	VECTOR_OB V;
	STRING_OP s0, s1, s2, s;
	INT i, len1;
	BYTE *p_name;

	fp = fopen(fname, "w");
	
	fputs("<Head>\n", fp);
	fputs("<Title>\n", fp);
	fputs("Seminaire, country\n", fp);
	fputs("</Title>\n", fp);
	fputs("</Head>\n", fp);
	
	fputs("<Body>\n", fp);
	fputs("<H1>\n", fp);
	fputs("Seminaire, place\n", fp);
	fputs("</H1>\n", fp);
	
	btree2->len(&len1);
	fputs("<ul>\n", fp);
	for (i = 0; i < len1; i++) {
		btree2->ith(i, &key, &data);
		db->get_op(&data, &V);
		s0 = (STRING_OP) V.s_i(0);
		s1 = (STRING_OP) V.s_i(1);
		s2 = (STRING_OP) V.s_i(2);
		if (V.s_li() > 4) {
			s = (STRING_OP) V.s_i(4);
			p_name = s->s_str();
			}
		else
			p_name = NIL;
		sprintf(label, "%s_%ld", s0->s_str(), i);
		sprintf(str, "<li> "
			"<A HREF=\"seminaire.html#%ld\"> %s, %s </A> ", 
			(INT) data.datref, s2->s_str(), p_name);
		fputs(str, fp);
		fputs("\n", fp);
		fflush(fp);
		}
	fputs("</ul>\n", fp);
	fputs("</Body>\n\n", fp);
	fflush(fp);
	fclose(fp);
	printf("written file %s of size %ld\n", 
		fname, file_size(fname));
	fflush(stdout);
	return OK;
}

static INT HTML_seminaire(
	DATABASE_OP db, BAYERTREE_OP btree0)
{
	BYTE *fname = "seminaire.html";
	BYTE label[256];
	FILE *fp;
	KEYTYPE key;
	DATATYPE data;
	BYTE str[1024];
	VECTOR_OB V;
	STRING_OP s0, s1, s2, s;
	INT i, len1, j;

	fp = fopen(fname, "w");
	
	fputs("<Head>\n", fp);
	fputs("<Title>\n", fp);
	fputs("Seminaire, adress\n", fp);
	fputs("</Title>\n", fp);
	fputs("</Head>\n", fp);
	
	fputs("<Body>\n", fp);
	fputs("<H1>\n", fp);
	fputs("Seminaire, adress\n", fp);
	fputs("</H1>\n", fp);
	
	btree0->len(&len1);
	for (i = 0; i < len1; i++) {
		btree0->ith(i, &key, &data);
		db->get_op(&data, &V);
		s0 = (STRING_OP) V.s_i(0);
		s1 = (STRING_OP) V.s_i(1);
		s2 = (STRING_OP) V.s_i(2);
		sprintf(label, "%ld", (INT) data.datref);
		sprintf(str, "<A NAME=\"%s\">*</A>", label);
		fputs(str, fp);
		fputs("\n", fp);
		fputs("<PRE>\n", fp);
		for (j = 4; j < V.s_li(); j++) {
			s = (STRING_OP) V.s_i(j);
			sprintf(str, "%s", s->s_str());
			fputs(str, fp);
			fputs("\n", fp);
			}
		fputs("</PRE>\n", fp);
		fputs("<P>\n", fp);
		fflush(fp);
		}
	fputs("</Body>\n\n", fp);
	fflush(fp);
	fclose(fp);
	printf("written file %s of size %ld\n", 
		fname, file_size(fname));
	fflush(stdout);
	return OK;
}

INT db_custom_exp_tex(INT btree_idx, BYTE *tex_file_name)
{
	BAYERTREE_OP btree;
	DATABASE_OP db;
	INT i, len_db;
	FILE *f;
	KEYTYPE key_type;
	DATATYPE data_type;
	VECTOR_OB V;
	STRING_OP s;
	BYTE *p;
	INT j, l, f_large;
	
	if (i_custom_db(&db) != OK)
		return error("db_custom_exp_tex(): "
		"error in i_custom_db()");
	btree = db->s_btree_i(btree_idx);
	/* btree->print_pages(); */

	if (db->open() != OK)
		return error("db_custom_exp_tex(): "
		"error in db->open()");
	f = fopen(tex_file_name, "w");
	fprintf(f, "%%format latex\n\n");
	fprintf(f, "\\documentstyle[twocolumn,"
		"art10]{article}\n");
	fprintf(f, "\\textheight=640pt\n");
	fprintf(f, "\\textwidth=440pt\n");
	fprintf(f, "\\topmargin=0pt\n");
	fprintf(f, "\\headsep=18pt\n");
	fprintf(f, "\\footskip=45pt\n");
	fprintf(f, "\\mathsurround=1pt\n");
	fprintf(f, "\\evensidemargin=15pt\n");
	fprintf(f, "\\oddsidemargin=15pt\n");
	fprintf(f, "\\parindent=0pt\n");
	fprintf(f, "\\columnseprule=0.1mm\n");
	fprintf(f, "\\newcommand{\\WIDTH}{6cm}\n");
	fprintf(f, "\\begin{document}\n");
	fprintf(f, "\\twocolumn\n");
	if (btree->len(&len_db) != OK) {
		Srfs("db_custom_exp_tex", "BT::len");
		return ERROR;
		}
	for (i = 0; i < len_db; i++) {
		if (btree->ith(i, 
			&key_type, &data_type) != OK) {
			Srff("db_custom_exp_tex", "BT::ith");
			return ERROR;
			}
		if (db->get_op(&data_type, &V) != OK) {
			Srff("db_custom_exp_tex", "DB::get_op");
			return ERROR;
			}
		fprintf(f, "\\parbox[b]{\\WIDTH}{%%\n");
		l = V.s_li();
		for (j = 4; j < l; j++) {
			s = (STRING_OP) V.s_i(j);
			p = s->s_str();
			if (*p == '#')
				continue;
			f_large = FALSE;
			
			if (!strncmp(p, "\\l", 2)) {
				f_large = TRUE;
				p += 2;
				}
			if (j == 4)
				f_large = TRUE;
			if (f_large)
				fprintf(f, "{\\large %s}\\\\\n", p);
			else
				fprintf(f, "%s\\\\\n", p);
			}
		fprintf(f, "}\n");
		fprintf(f, "\\par\n\n");
		}

	fprintf(f, "\\end{document}\n");
	fclose(f);
	printf("db_custom_exp_tex()|wrote file "
		"%s of size %ld\n", 
		tex_file_name, 
		file_size(tex_file_name));
	db->close();
	e_custom_db(db);
	return OK;
}

#endif /* DB_TRUE */




