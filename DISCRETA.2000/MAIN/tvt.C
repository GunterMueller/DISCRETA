// tvt.C

#include <sys/stat.h>
#include <unistd.h> 



static INT flush_dp_data(VECTOR_OP dp_data);
static INT search_in_all_data(VECTOR_OP all_data, INT t, INT v, INT k, INT lambda);
static INT add_in_all_data(VECTOR_OP all_data, INT t, INT v, INT k, INT lambda);
static INT add_to_db(DATABASE_OP db, VECTOR_OP dp_data, 
	INT t, INT v, INT k, INT lambda, INT id, BYTE *text);
INT suche_rechten_partner(VECTOR_OP all_data, INT t, INT v, INT k, INT l, 
	INT *t_new, INT *v_new, INT *k_new, INT *l_new, BYTE *text);
INT suche_linken_partner(VECTOR_OP all_data, INT t, INT v, INT k, INT l, 
	INT *t_new, INT *v_new, INT *k_new, INT *l_new, BYTE *text);
INT design_da(INT t, INT v, INT k, INT l);


#if 0
INTEGER_OP Len = NIL;
VECTOR_OP Key = NIL;
VECTOR_OP Vater = NIL;
VECTOR_OP Mutter = NIL;
VECTOR_OP Basis = NIL;
INTEGER_OP LenB = NIL;
static INT id_max;
#endif



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
	INT t_new, v_new, k_new, lambda_new;
	FILE *fp1, *fp2;
	BYTE *fnew = "new_design.txt";
	BYTE *ffikt = "fikt_design.txt";
	BYTE cmd[256];
	INT id_max, id_cur;
	BYTE text[1024];
	
	VECTOR_OB all_data;
	VECTOR_OB d;
	
	
	VECTOR_OB dp_data;

	btree->len(&len1);
	
	btree->ith(len1 - 1, &key, &data);
	db->get_op(&data, &dp);
	id_max = dp.s_id_i();
	printf("highest id in database: %ld\n", id_max);
	id_cur = id_max;
	
	call_system("rm tvt.log");
	
	dp_data.m_il(0);

	printf("loading complete database\n");
	fflush(stdout);
	all_data.m_il(len1);
	for (i = 0; i < len1; i++) {
		btree->ith(i, &key, &data);
		db->get_op(&data, &dp);
		t = dp.s_t_i();
		v = dp.s_v_i();
		k = dp.s_k_i();
		lambda = dp.s_lambda_i();
		// id = dp.s_id_i();
		d.m_il(4);
		d.m_ii(0, t);
		d.m_ii(1, v);
		d.m_ii(2, k);
		d.m_ii(3, lambda);
		d.copy((VECTOR_OP) all_data.s_i(i));
		}
	printf("sorting vector of length %ld\n", len1); 
	fflush(stdout);
	all_data.quicksort(len1, TRUE /* ascending */);
	printf("finished\n");
	fflush(stdout);
	
	for (i = 0; i < len1; i++) {
		btree->ith(i, &key, &data);
		db->get_op(&data, &dp);
		t = dp.s_t_i();
		v = dp.s_v_i();
		k = dp.s_k_i();
		lambda = dp.s_lambda_i();
		id = dp.s_id_i();
		// printf("db_dp_tvt_closure() t=%ld v=%ld k=%ld lambda=%ld id=%ld\n", t, v, k, lambda, id);
		
		if ((i + 1) % 10 == 0) {
			printf(",");
			fflush(stdout);
			}
		else
			printf(".");
		if ((i + 1) % 50 == 0) {
			printf(" %ld\n", i + 1);
			fflush(stdout);
			}
		if (suche_rechten_partner(&all_data, (INT) t, (INT) v, (INT) k,(INT) lambda, 
			&t_new, &v_new, &k_new, &lambda_new, text)) {
			add_to_db(db, &dp_data, 
				t_new, v_new, k_new, lambda_new, ++id_cur, text);
			flush_dp_data(&dp_data);
			add_in_all_data(&all_data, t_new, v_new, k_new, lambda_new);
			}
		
		
		
		if (suche_linken_partner(&all_data, (INT) t, (INT) v, (INT) k, (INT) lambda, 
			&t_new, &v_new, &k_new, &lambda_new, text)) {
			add_to_db(db, &dp_data, 
				t_new, v_new, k_new, lambda_new, ++id_cur, text);
			flush_dp_data(&dp_data);
			add_in_all_data(&all_data, t_new, v_new, k_new, lambda_new);
			}
		
		}
	
	





	return OK;
}


static INT flush_dp_data(VECTOR_OP dp_data)
{
	DESIGN_PARAMETER_OP dp1;
	INT l = dp_data->s_li();
	INT t_, v_, k_, lambda_, id_;
	VECTOR_OB vec_data;
	MATRIX_OB mat_data;
	INT i;
	
	mat_data.m_ilih(5, l);
	for (i = 0; i < l; i++) {
		dp1 = (DESIGN_PARAMETER_OP) dp_data->s_i(i);
		v_ = dp1->s_v_i();
		t_ = dp1->s_t_i();
		k_ = dp1->s_k_i();
		lambda_ = dp1->s_lambda_i();
		id_ = dp1->s_id_i();
		mat_data.m_iji(i, 0, t_);
		mat_data.m_iji(i, 1, v_);
		mat_data.m_iji(i, 2, k_);
		mat_data.m_iji(i, 3, lambda_);
		mat_data.m_iji(i, 4, id_);
		}
	vec_data.m_il(1);
	mat_data.swap((MATRIX_OP) vec_data.s_i(0));
	
	db_add_vec_data(&vec_data);
	dp_data->m_il(0);
	return OK;
}

static INT search_in_all_data(VECTOR_OP all_data, INT t, INT v, INT k, INT lambda)
{
	VECTOR_OB d;
	INT idx, f_found;
	
	d.m_il(4);
	d.m_ii(0, t);
	d.m_ii(1, v);
	d.m_ii(2, k);
	d.m_ii(3, lambda);
	all_data->search(all_data->s_li(), TRUE /* f_ascending */, &d, 
		&idx, &f_found);
	return f_found;
}

static INT add_in_all_data(VECTOR_OP all_data, INT t, INT v, INT k, INT lambda)
{
	VECTOR_OB d;
	INT idx, f_found, l, i;
	
	d.m_il(4);
	d.m_ii(0, t);
	d.m_ii(1, v);
	d.m_ii(2, k);
	d.m_ii(3, lambda);
	all_data->search(all_data->s_li(), TRUE /* f_ascending */, &d, 
		&idx, &f_found);
	l = all_data->s_li();
	all_data->inc();
	for (i = l - 1; i >= idx; i--) {
		all_data->s_i(i)->swap(all_data->s_i(i + 1));
		}
	d.swap((VECTOR_OP) all_data->s_i(idx));
	return OK;
}

static INT add_to_db(DATABASE_OP db, VECTOR_OP dp_data, 
	INT t, INT v, INT k, INT lambda, INT id, BYTE *text)
{
	DESIGN_PARAMETER_OB p;
	DESIGN_PARAMETER_SOURCE_OB S;
	
	p.init();
	p.s_v()->m_i(v);
	p.s_t()->m_i(t);
	p.s_k()->m_i(k);
	p.s_lambda()->m_i(lambda);
	p.s_id()->m_i(id);


	S.init();
	S.s_prev()->m_i(-1);
	S.s_rule()->m_i(DP_RULE_TRUNG);
	S.s_comment()->init(text);
	p.s_source()->append_in_place(&S);


	if (db->add_op(&p) != OK) {
		db->close();
		return error("add_to_vec_data() error in DB::add_op");
		}
	dp_data->inc();
	p.copy((DESIGN_PARAMETER_OP) dp_data->s_i(dp_data->s_li() - 1));
	return OK;
}


INT suche_rechten_partner(VECTOR_OP all_data, INT t, INT v, INT k, INT l, 
	INT *t_new, INT *v_new, INT *k_new, INT *l_new, BYTE *text)
{
    // sucht einen rechten Partner fuer t-(v,k,l) in der DB fuer TvT
    FILE *datei;
    INT i, size, size_alt, size_design, l1;
    char c;
    struct stat buf;
    INT a, b;
    
    // das lambda fuer rechten Partner wird berechnet
    a = l*(v-k);
    b = (k+1-t);
    if (b==0) return FALSE;
    if ((a % b) != 0) {
    	// printf("warning: a % b != 0\n");
	return FALSE;
    	}
    l1=a/b;    
    
    if (!search_in_all_data(all_data, t, v, k + 1, l1))
    	return FALSE;
	
#if 0
    // Abfragedatei fuer die Designs
    datei=fopen("tvt_r.sdf","w+");
    fprintf(datei,"from despar\nto tvt_rp.out\nwhere t\n>= %i\n<= %i\n",t,t);
    fprintf(datei,"where v\n>= %i\n<= %i\nwhere k\n>= %i\n<= %i\n",v,v,k+1,k+1);
    fprintf(datei,"where lambda\n>= %i\n<= %i\n",l1,l1);
    fclose(datei);        
    // Ausfuehren der Abfrage
    system("discreta_dbselect.out tvt_r.sdf >tvtselect.log");   
    // In tvt_r.out steht das gefundene Design (oder keins)
    datei=fopen("tvt_rp.out","rb");
    stat("tvt_rp.out",&buf); 
    fclose(datei);
    size=buf.st_size;
    // Falls size==0 wurde kein Design gefunden
    if (size==0) return FALSE;
#endif
    if (!design_da(t,v+1,k+1,l+l1))     
    {
       datei=fopen("tvt.log","a");
       sprintf(text, "TvT: %i-(%i, %i, %i) with %i-(%i, %i, %i)",t,v,k,l,t,v,k+1,l1);
       fprintf(datei,"%d %d %d %d %s\n", t,v+1,k+1,l+l1, text);
       fclose(datei); 
       // daten_eintragen();
       *t_new = t;
       *v_new = v + 1;
       *k_new = k + 1;
       *l_new = l + l1;
       return TRUE;
    }   
    return FALSE;
}

INT suche_linken_partner(VECTOR_OP all_data, INT t, INT v, INT k, INT l, 
	INT *t_new, INT *v_new, INT *k_new, INT *l_new, BYTE *text)
{
    // sucht einen linken Partner fuer t-(v,k,l) in der DB fuer TvT 
    FILE *datei;
    INT i, size, size_alt, size_design, l1;
    char c;
    struct stat buf;
    // das lambda fuer linken Partner wird berechnet
    
    if ((v-k+1)==0) return FALSE;
    
    l1=l*(k-t)/(v-k+1);
    if (!search_in_all_data(all_data, t, v, k - 1, l1))
    	return FALSE;
	
#if 0
    // Abfragedatei fuer die Designs
    datei=fopen("tvt_l.sdf","w+");
    fprintf(datei,"from despar\nto tvt_lp.out\nwhere t\n>= %i\n<= %i\n",t,t);
    fprintf(datei,"where v\n>= %i\n<= %i\nwhere k\n>= %i\n<= %i\n",v,v,k-1,k-1);
    fprintf(datei,"where lambda\n>= %i\n<= %i\n",l1,l1);
    fclose(datei);    
    // Ausfuehren der Abfrage 
    system("discreta_dbselect.out tvt_l.sdf >tvtselect.log"); 
    // In tvt_r.out steht das gefundene Design (oder keins)
    datei=fopen("tvt_lp.out","rb");
    stat("tvt_lp.out",&buf); 
    fclose(datei);
    size=buf.st_size;
    // Falls size==0 wurde kein Design gefunden
    if (size==0) return FALSE;
#endif
    if (!design_da(t,v+1,k,l+l1)) 
    {
         datei=fopen("tvt.log","a");
	 sprintf(text, "TvT: %i-(%i, %i, %i) with %i-(%i, %i, %i)",t,v,k,l,t,v,k-1,l1);
         fprintf(datei,"%d %d %d %d %s\n", t,v+1,k,l+l1,text);
         fclose(datei);
         // daten_eintragen(t,v+1,k,l+l1);
       *t_new = t;
       *v_new = v + 1;
       *k_new = k;
       *l_new = l + l1;
       return TRUE;
    }	 
    return FALSE; 
}


INT design_da(INT t, INT v, INT k, INT l)
{  
   FILE *datei;
   INT size, wert=0;
   struct stat buf;
   datei=fopen("tvt_t.sdf","w+");
   fprintf(datei,"from despar\nto tvt_t.out\nwhere t\n>= %i\n<= %i\n",t,t);
   fprintf(datei,"where v\n>= %i\n<= %i\nwhere k\n>= %i\n<= %i\n",v,v,k,k);
   fprintf(datei,"where lambda\n>= %i\n<= %i\n",l,l);
   fclose(datei); 
   system("discreta_dbselect.out tvt_t.sdf >tvtselect.log"); 
   stat("tvt_t.out",&buf);     
   size=buf.st_size;
   // Falls size==0 wurde kein Design gefunden
   if (size!=0) wert=1;
   return wert;
}




#if 0
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



#endif

