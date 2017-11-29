
INTEGER_OP Len = NIL;
VECTOR_OP Key = NIL;
VECTOR_OP Vater = NIL;
VECTOR_OP Mutter = NIL;
VECTOR_OP Basis = NIL;
INTEGER_OP LenB = NIL;
static INT id_max;

static INT db_dp_tvt_closure1(DATABASE_OP db, BAYERTREE_OP btree);
INT lambda_test(INT t,INT v, INT k,INT lambda);
INT print_ausgabe_design(FILE *fp, INT i, char* text);
INT tvt_test(INT t, INT v, INT k, INT lambda,INT id);
INT tvt_ausgabe(FILE *fp1 /* new_design */, FILE *fp2 /* fikt_design */);

void db_dp_tvt_closure()
{
	DATABASE_OP db;
	BAYERTREE_OP btree0, btree1, btree2;
	
	i_db_dp(&db, DP_PATH);
	
	btree0 = db->s_btree_i(0);
	btree1 = db->s_btree_i(1);
	btree2 = db->s_btree_i(2);
	/* btree->print_pages(); */
	
	if (db->open() != OK) {
		error("db_dp_tvt_closure(): error in db->open()");
		return;
		}
	
	db_dp_tvt_closure1(db, btree0);
	
	db->close();
	e_db_dp(db);
	
}


static INT db_dp_tvt_closure1(DATABASE_OP db, BAYERTREE_OP btree)
{
	KEYTYPE key;
	DATATYPE data;
	DESIGN_PARAMETER_OB dp;
	INT i, len1;
	INT t, v, k, lambda, id;
	FILE *fp1, *fp2;
	BYTE *fnew = "new_design.txt";
	BYTE *ffikt = "fikt_design.txt";
	BYTE cmd[256];



	LenB = (INTEGER_OP) callocobject("LenB");
	Basis = (VECTOR_OP) callocobject("Basis");
	Len = (INTEGER_OP) callocobject("Len");
	Key = (VECTOR_OP) callocobject("Key");
	Vater = (VECTOR_OP) callocobject("Vater");
	Mutter = (VECTOR_OP) callocobject("Mutter");

	LenB->m_i(0);
	Basis->m_il(VECTOR_OVERSIZE);	
	Len->m_i(0);
	Key->m_il(VECTOR_OVERSIZE);
	Vater->m_il(VECTOR_OVERSIZE);
	Mutter->m_il(VECTOR_OVERSIZE);
	
	
	btree->len(&len1);
	
	btree->ith(len1 - 1, &key, &data);
	db->get_op(&data, &dp);
	id_max = dp.s_id_i();
	printf("highest id in database: %ld\n", id_max);
	
	for (i = 0; i < len1; i++) {
		btree->ith(i, &key, &data);
		db->get_op(&data, &dp);
		t = dp.s_t_i();
		v = dp.s_v_i();
		k = dp.s_k_i();
		lambda = dp.s_lambda_i();
		id = dp.s_id_i();
		printf("db_dp_tvt_closure() t=%ld v=%ld k=%ld lambda=%ld id=%ld\n", 
			t, v, k, lambda, id);
		
		if (!tvt_test(t, v, k, lambda,id))
			printf("Fehler beim Einfuegen");
		}
	
	fp1 = fopen(fnew, "w");
	fp2 = fopen(ffikt, "w");
	
	
	tvt_ausgabe(fp1, fp2);	

	fclose(fp1);
	fclose(fp2);

	freeall(LenB);
	freeall(Basis);
	freeall(Len);
	freeall(Key);
	freeall(Vater);
	freeall(Mutter);
	
	printf("written file %s of size %ld\n", 
		fnew, file_size(fnew));
	fflush(stdout);

	sprintf(cmd, "cat %s", fnew);
	call_system(cmd);
	
	return OK;
}

/* tvt.C
   (c) Anton Betten , Ulrich Mechtold 22.04.1997 
 */

// #include <DISCRETA/discreta.h>
#include <ctype.h>

#ifndef BUFSIZE
#define BUFSIZE 10000
#endif

//FILE *fp1,*fp2;


INT lambda_test(INT t,INT v, INT k,INT lambda)
{
        INT i, nom, denom, n, d, g, lambda_new;

	nom = 1;
	denom = 1;
	n = v - t;
	d = k - t;
	for (i = 0; i < (k - t); i++) {
		nom *= n;
		denom *= d;
		n--;
		d--;
		ggt_iipi(nom, denom, &g);
		if (g != 1 && g != -1) {
			nom /= g;
			denom /= g;
			}
		}
	if (denom != 1)
		return error("DP::complement() error: denom != 1");
	lambda_new = nom - lambda;
	if (lambda_new <= 0)
		return error("DP::complement() error: lambda_new <= 0");
	if (lambda_new <lambda) return lambda_new;
	else return lambda;
	
}



INT print_ausgabe_design(FILE *fp, INT i, char* text)
{
	VECTOR_OP key;
	
	key = (VECTOR_OP) Key->s_i(i);
	fprintf(fp, "%ld %ld %ld %ld fiktiv von  %s \n", 
		key->s_ii(0), 
		key->s_ii(1), 
		key->s_ii(2), 
		key->s_ii(3),
		text);
	return OK;
}

INT tvt_test(INT t, INT v, INT k, INT lambda,INT id)
{
	INTEGER_OB f_v, f_m;
	VECTOR_OB V, m, ur;
	INT idx, f_found, vlambda,mlambda, urlambda;
	INT hilf;
        
	urlambda=lambda_test(t, v, k, lambda);
	
	ur.m_il(4);
	ur.m_ii(0, t );
	ur.m_ii(1, v );
	ur.m_ii(2, k );
	ur.m_ii(3, urlambda);
	
	Basis->search(LenB->s_i(), TRUE, &ur, &idx,&f_found);
	if (!f_found) {
		Basis->insert_at(LenB, idx, &ur);
 		} 

        vlambda=lambda_test(t+1, v+1, k+1, lambda);
	
	V.m_il(4);
	V.m_ii(0, t+1 );
	V.m_ii(1, v+1 );
	V.m_ii(2, k+1 );
	V.m_ii(3, vlambda);
	
	Key->search(Len->s_i(), TRUE/* f_ascending */, &V,
        	&idx, &f_found);
	if (!f_found) {
		Key->insert_at(Len, idx, &V);
		f_v.m_i(id);
		f_m.m_i(-1);
		Len->dec();
		Vater->insert_at(Len, idx, &f_v);
		Len->dec();
		Mutter->insert_at(Len, idx, &f_m);
		}
	else {
		idx--;
		Vater->m_ii(idx, id);
		
		}	
        
	hilf= ((lambda *(k-t))%(v-k+1));
	if (hilf==0) {
	
		mlambda=((lambda *(k-t))/(v-k+1));
        	mlambda=lambda_test(t+1, v+1, k,mlambda );
       
		m.m_il(4);
		m.m_ii(0, t+1 );
		m.m_ii(1, v+1 );
		m.m_ii(2, k);
        	m.m_ii(3, mlambda);

	

		Key->search(Len->s_i(), TRUE/* f_ascending */, &m,
        	&idx, &f_found);
		if (!f_found) {
			Key->insert_at(Len, idx, &m);
			f_v.m_i(-1);
			f_m.m_i(id);
			Len->dec();
			Vater->insert_at(Len, idx, &f_v);
			Len->dec();
			Mutter->insert_at(Len, idx, &f_m);
			}
			else {
			idx--;
			Mutter->m_ii(idx, id);
			
			}
	}
	return (1);
}


#if 0

INT read_base(BYTE *fname)
{
	FILE *fp;
	INT i, l, t, v, k, lambda, line = 0, f_has_comment;
	DATABASE_OP db;
	BYTE buf[BUFSIZE], *p, *comment;
	INT f_v = TRUE;
	INT id=0;
	fp = fopen(fname, "r");

	
	while (TRUE) {
		buf[0] = 0;
		line++;
		if (fgets(buf, BUFSIZE, fp) == NULL) {
			break;
			}
		l = strlen(buf);
		if (l) {
			if (buf[l - 1] == '\n')
				buf[l - 1] = 0;
			}
		printf("%ld: %s\n", line, buf);
		fflush(stdout);
		p = buf;
		
		if (!s_scan_int(&p, &t)) {
			printf("can't read t !\n");
			break;
			}
		if (t <= 0)
			break;
		if (!s_scan_int(&p, &v)) {
			printf("can't read t !\n");
			break;
			}
		if (!s_scan_int(&p, &k)) {
			printf("can't read t !\n");
			break;
			}
		if (!s_scan_int(&p, &lambda)) {
			printf("can't read t !\n");
			break;
			}
		f_has_comment = FALSE;
		for (i = 0; i < strlen(p); i++) {
			if (!isspace(p[i])) {
				f_has_comment = TRUE;
				break;
				}
			}
		printf("%d *** %ld-(%ld,%ld,%ld)", id,t, v, k, lambda);

		if (!tvt_test(t, v, k, lambda,id))
			printf("Fehler beim Einfuegen");
                  
		id++;
	
		if (f_has_comment) {
			comment = p;
			printf(" %s:\n", comment);
			}
		else {
			comment = NIL;
			printf(":\n");
			}
		}
	fclose(fp);
	return OK;
}

#endif


INT tvt_ausgabe(FILE *fp1 /* new_design */, FILE *fp2 /* fikt_design */)
{
	VECTOR_OP key; 
	VECTOR_OB alt;
	INT i, l, fertig,z,iterat, neueanz;
	INT t, k, v, lambda;
	INT t_neu, v_neu, k_neu, lambda_neu,id,id1,id2,idx,alt_found;
	char buffer[120];	

	id=id_max + 1;
	z=0;
	neueanz=0;
	iterat=1;
	fertig=FALSE;

	while((!fertig) && (iterat < 20))
	{
		fertig=TRUE; 
		l = Len->s_i();
		/*fprintf(stdout, "Laenge %d \n",l);*/ 
		for (i = 0; i < l; i++) {
			if ((Vater->s_ii(i)!=-1) && (Mutter->s_ii(i)!=-1)) {
				key = (VECTOR_OP) Key->s_i(i);
		
				t=key->s_ii(0); 
				v=key->s_ii(1); 
				k=key->s_ii(2); 
				lambda=key->s_ii(3);
				id1=Vater->s_ii(i);
				id2=Mutter->s_ii(i);

				/*  red */

				t_neu=t-1 ;		
				v_neu=v;
				k_neu=k;
				lambda_neu=(lambda * (v - t + 1)) /(k - t +1);	

				lambda_neu=lambda_test(t_neu,v_neu,k_neu,lambda_neu);
				
				alt.m_il(4);
				alt.m_ii(0, t_neu );
				alt.m_ii(1, v_neu );
				alt.m_ii(2, k_neu );
				alt.m_ii(3, lambda_neu);
	
				Basis->search(LenB->s_i(), TRUE,
&alt,&idx,&alt_found);
				if (!alt_found) {
				
					/* Ausgabe */
				
					sprintf(buffer,"%d-(%d,%d,%d) $\\cup$  "
					"%d-(%d,%d,%d)",
					t-1,v-1,k-1,lambda,
					t-1,v-1,k,((lambda* (v - t +1))/(k-t+1) - lambda));


					fprintf(fp1, "%ld %ld %ld %ld TvT: %s \n", t_neu, v_neu, k_neu, lambda_neu, buffer); 
			
					sprintf(buffer," ID %d der: %d-(%d,%d,%d) %d, "
					"res: %d-(%d,%d,%d) %d red: %d-(%d,%d,%d)",
					id,t-1,v-1,k-1,lambda,id1,
					t-1,v-1,k,((lambda* (v - t +1))/(k-t+1) - lambda),id2,t_neu,v_neu,k_neu,lambda_neu );

					printf( "Neu: %ld %ld %ld %ld TvT %s \n", t_neu, v_neu, k_neu, lambda_neu, buffer); 
					/* fflush(stdout); */
	
					print_ausgabe_design(fp2,i,buffer);
					neueanz++;
				}
			
				Key->delete_ith(Len,i);
				Len->inc();
				Vater->delete_ith(Len,i);
				Len->inc();
				Mutter->delete_ith(Len,i); 
				
				
				/* printf(" %d. aus Nummer %d \n", z,i);*/
				
				
				if (!alt_found) {
							
					if(tvt_test( t_neu,v_neu,k_neu,lambda_neu,id)) {
					id++;
					}		
					else fprintf(stderr,"Fehler bei tvt_test");
				} 

				l = Len->s_i();
				fertig=FALSE;
				z++;

										
			}
			/*else fprintf(stdout, "keines gefunden \n");*/
		}
	iterat++;
	}
  printf("Iterationen: %d , gefunden insg: %d , davon neu: %d \n",iterat,z, neueanz);
	fflush(fp1);
	fflush(fp2);
  return OK;
}


#if 0
int main(int argc, char **argv)
{
	INT t0, t1, user_time;
	BYTE s[256];
	
	
	discreta_init();
	{
	
	t0 = os_ticks();

	LenB = (INTEGER_OP) callocobject("LenB");
	Basis = (VECTOR_OP) callocobject("Basis");
	Len = (INTEGER_OP) callocobject("Len");
	Key = (VECTOR_OP) callocobject("Key");
	Vater = (VECTOR_OP) callocobject("Vater");
	Mutter = (VECTOR_OP) callocobject("Mutter");

	LenB->m_i(0);
	Basis->m_il(VECTOR_OVERSIZE);	
	Len->m_i(0);
	Key->m_il(VECTOR_OVERSIZE);
	Vater->m_il(VECTOR_OVERSIZE);
	Mutter->m_il(VECTOR_OVERSIZE);
	
	read_base("design.txt");
	
	fp1 = fopen("new_design.txt", "w");
	fp2 = fopen("fikt_design.txt", "w");
	
	
	tvt_ausgabe();	

	fclose(fp1);
	fclose(fp2);


	t1 = os_ticks();
	user_time = t1 - t0;
	s[0] = 0;
	print_delta_time(user_time, s);
	printf("total computing time: %s\n", s);
	fflush(stdout);
	}
	discreta_exit();
	return 0;
}
#endif



