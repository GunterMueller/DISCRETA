#define DISCRETA_OUTPUT 1 \

#define max_level 200
#define max_degree 500
#define max_cols 3000
#define max_nodes 200000
#define name_size 4 \
 \

#define root col_array[0] \

#define buf_size (name_size+1) *max_cols+name_size \

#define panic(m) {fprintf(stderr,"%s!\n%s",m,buf) ;exit(-1) ;} \

/*1:*/
#line 32 "dance.w"

#include <stdio.h> 
#include <stdlib.h> 
#include <ctype.h> 
#include <string.h> 
/*3:*/
#line 101 "dance.w"

typedef struct node_struct{
struct node_struct*left,*right;
struct node_struct*up,*down;
int rownumber;
struct col_struct*col;
}node;

/*:3*//*4:*/
#line 122 "dance.w"

typedef struct col_struct{
node head;
int len;
char name[name_size+1];
struct col_struct*prev,*next;
}column;

/*:4*/
#line 37 "dance.w"

/*2:*/
#line 64 "dance.w"

int verbose;
int SILENT;
FILE*txt;

int count= 0;
#if 0
unsigned int updates;
#endif
double updates;
int spacing= 1;
int profile[max_level][max_degree];
#if 0
unsigned int upd_prof[max_level];
#endif
double upd_prof[max_level];
int maxb= 0;
int maxl= 0;

/*:2*//*10:*/
#line 175 "dance.w"

column col_array[max_cols+2];
node node_array[max_nodes];
char buf[buf_size];

/*:10*//*14:*/
#line 250 "dance.w"

int number_of_rows;
int*solution_vector;
#if DISCRETA_OUTPUT
int*selectedcolumns;
int overallrows;
#endif

/*:14*//*18:*/
#line 321 "dance.w"

int level;
node*choice[max_level];

/*:18*/
#line 38 "dance.w"

/*6:*/
#line 140 "dance.w"

void print_row(p)
node*p;
{register node*q= p;
register int k;

do{
printf(" %s",q->col->name);
q= q->right;
}while(q!=p);
for(q= p->col->head.down,k= 1;q!=p;k++)
if(q==&(p->col->head)){
printf("\n");return;
}else q= q->down;
printf(" (%d of %d)\n",k,p->col->len);
}

/*:6*//*19:*/
#line 329 "dance.w"

void cover(c)
column*c;
{register column*l,*r;
register node*rr,*nn,*uu,*dd;
register int k= 1;
l= c->prev;r= c->next;
l->next= r;r->prev= l;
for(rr= c->head.down;rr!=&(c->head);rr= rr->down)
for(nn= rr->right;nn!=rr;nn= nn->right){
uu= nn->up;dd= nn->down;
uu->down= dd;dd->up= uu;
k++;
nn->col->len--;
}
updates+= k;
upd_prof[level]+= k;
}

/*:19*//*20:*/
#line 352 "dance.w"

void uncover(c)
column*c;
{register column*l,*r;
register node*rr,*nn,*uu,*dd;
for(rr= c->head.up;rr!=&(c->head);rr= rr->up)
for(nn= rr->left;nn!=rr;nn= nn->left){
uu= nn->up;dd= nn->down;
uu->down= dd->up= nn;
nn->col->len++;
}
l= c->prev;r= c->next;
l->next= r->prev= c;
}

/*:20*/
#line 39 "dance.w"
;

int main(argc,argv)
int argc;
char*argv[];
{
/*12:*/
#line 210 "dance.w"

register column*cur_col;
register char*p,*q;
register node*cur_node;
int primary;
int i;

/*:12*//*17:*/
#line 317 "dance.w"

register column*best_col;
register node*pp;

/*:17*//*24:*/
#line 398 "dance.w"

register int minlen;
register int j,k,x;

/*:24*/
#line 45 "dance.w"
;

SILENT= 0;
for(i= 1;i<argc;i++){
if(strcmp(argv[i],"silent")==0)SILENT= 1;
}
verbose= argc-1-SILENT;
if(verbose)sscanf(argv[1+SILENT],"%d",&spacing);
/*9:*/
#line 169 "dance.w"

/*11:*/
#line 182 "dance.w"

cur_col= col_array+1;
fgets(buf,buf_size,stdin);
if(buf[strlen(buf)-1]!='\n')panic("Input line too long");
for(p= buf,primary= 1;*p;p++){
while(isspace(*p))p++;
if(!*p)break;
if(*p=='|'){
primary= 0;
if(cur_col==col_array+1)panic("No primary columns");
(cur_col-1)->next= &root,root.prev= cur_col-1;
continue;
}
for(q= p+1;!isspace(*q);q++);
if(q> p+name_size)panic("Column name too long");
if(cur_col>=&col_array[max_cols])panic("Too many columns");
for(q= cur_col->name;!isspace(*p);q++,p++)*q= *p;
cur_col->head.up= cur_col->head.down= &cur_col->head;
cur_col->len= 0;
if(primary)cur_col->prev= cur_col-1,(cur_col-1)->next= cur_col;
else cur_col->prev= cur_col->next= cur_col;
cur_col++;
}
if(primary){
if(cur_col==col_array+1)panic("No primary columns");
(cur_col-1)->next= &root,root.prev= cur_col-1;
}

/*:11*/
#line 170 "dance.w"
;
/*13:*/
#line 217 "dance.w"

cur_node= node_array;
i= 0;
while(fgets(buf,buf_size,stdin)&&strstr(buf,"SELECTEDCOLUMNS")==NULL){
register column*ccol;
register node*row_start;
if(buf[strlen(buf)-1]!='\n')panic("Input line too long");
row_start= NULL;
for(p= buf;*p;p++){
while(isspace(*p))p++;
if(!*p)break;
for(q= p+1;!isspace(*q);q++);
if(q> p+name_size)panic("Column name too long");
for(q= cur_col->name;!isspace(*p);q++,p++)*q= *p;
*q= '\0';
for(ccol= col_array;strcmp(ccol->name,cur_col->name);ccol++);
if(ccol==cur_col)panic("Unknown column name");
if(cur_node==&node_array[max_nodes])panic("Too many nodes");
if(!row_start)row_start= cur_node;
else cur_node->left= cur_node-1,(cur_node-1)->right= cur_node;
cur_node->col= ccol;
cur_node->up= ccol->head.up,ccol->head.up->down= cur_node;
cur_node->rownumber= i;
ccol->head.up= cur_node,cur_node->down= &ccol->head;
ccol->len++;
cur_node++;
}
if(!row_start)panic("Empty row");
row_start->left= cur_node-1,(cur_node-1)->right= row_start;
i++;
}
/*15:*/
#line 260 "dance.w"

number_of_rows= i;
solution_vector= (int*)calloc(number_of_rows,sizeof(int));
#if DISCRETA_OUTPUT
if(strstr(buf,"SELECTEDCOLUMNS")!=NULL){
scanf("%d\n",&overallrows);
selectedcolumns= (int*)calloc(overallrows,sizeof(int));
for(i= 0;i<overallrows;i++)
scanf("%d\n",&(selectedcolumns[i]));
}
#endif

/*:15*/
#line 248 "dance.w"
;

/*:13*/
#line 171 "dance.w"
;

/*:9*/
#line 53 "dance.w"
;

/*7:*/
#line 157 "dance.w"

#if DISCRETA_OUTPUT
if(!SILENT)txt= fopen("solutions","w");
#endif
/*:7*/
#line 55 "dance.w"
;
/*16:*/
#line 289 "dance.w"

level= 0;
forward:/*23:*/
#line 378 "dance.w"

minlen= max_nodes;
if(verbose> 2)printf("Level %d:",level);
for(cur_col= root.next;cur_col!=&root;cur_col= cur_col->next){
if(verbose> 2)printf(" %s(%d)",cur_col->name,cur_col->len);
if(cur_col->len<minlen)best_col= cur_col,minlen= cur_col->len;
}
if(verbose){
if(level> maxl){
if(level>=max_level)panic("Too many levels");
maxl= level;
}
if(minlen> maxb){
if(minlen>=max_degree)panic("Too many branches");
maxb= minlen;
}
profile[level][minlen]++;
if(verbose> 2)printf(" branching on %s(%d)\n",best_col->name,minlen);
}

/*:23*/
#line 291 "dance.w"
;
cover(best_col);
cur_node= choice[level]= best_col->head.down;
advance:if(cur_node==&(best_col->head))goto backup;
if(verbose> 1){
printf("L%d:",level);
print_row(cur_node);
}
/*21:*/
#line 367 "dance.w"

for(pp= cur_node->right;pp!=cur_node;pp= pp->right)cover(pp->col);

/*:21*/
#line 299 "dance.w"
;
if(root.next==&root)
#if DISCRETA_OUTPUT
/*26:*/
#line 416 "dance.w"

{
count++;
if(verbose)profile[level+1][0]++;
if(!SILENT&&count%spacing==0){
for(k= 0;k<number_of_rows;k++)solution_vector[k]= 0;
for(k= 0;k<=level;k++){
solution_vector[choice[k]->rownumber]= 1;
}

k= 0;
for(i= 0;i<overallrows;i++){
if(selectedcolumns[i]==0){
printf("0");
fprintf(txt,"0");
}else{
printf("%d",solution_vector[k]);
fprintf(txt,"%d",solution_vector[k]);
k++;
}
}
printf(" : %d\n",count);
fprintf(txt,"\n");
}
goto recover;
}

/*:26*/
#line 302 "dance.w"
;
#else
/*25:*/
#line 402 "dance.w"

{
count++;
if(verbose){
profile[level+1][0]++;

if(!SILENT&&count%spacing==0){
printf("%d:\n",count);
for(k= 0;k<=level;k++)print_row(choice[k]);
}
}
goto recover;
}

/*:25*/
#line 304 "dance.w"
;
#endif
level++;
goto forward;
backup:uncover(best_col);
if(level==0)goto done;
level--;
cur_node= choice[level];best_col= cur_node->col;
recover:/*22:*/
#line 375 "dance.w"

for(pp= cur_node->left;pp!=cur_node;pp= pp->left)uncover(pp->col);

/*:22*/
#line 312 "dance.w"
;
cur_node= choice[level]= cur_node->down;goto advance;
done:if(verbose> 3)
/*27:*/
#line 443 "dance.w"

{
printf("Final column lengths");
for(cur_col= root.next;cur_col!=&root;cur_col= cur_col->next)
printf(" %s(%d)",cur_col->name,cur_col->len);
printf("\n");
}

/*:27*/
#line 315 "dance.w"
;

/*:16*/
#line 56 "dance.w"
;
/*8:*/
#line 161 "dance.w"

#if DISCRETA_OUTPUT
if(!SILENT)fclose(txt);
#endif

/*:8*/
#line 57 "dance.w"
;

printf("Altogether %d solutions, after %.0f updates.\n",count,updates);
if(verbose)/*28:*/
#line 451 "dance.w"

{
x= 1;
for(level= 1;level<=maxl+1;level++){
j= 0;
for(k= 0;k<=maxb;k++){
printf("%6d",profile[level][k]);
j+= profile[level][k];
}
printf("%8d nodes, %f updates\n",j,upd_prof[level-1]);
x+= j;
}
printf("Total %d nodes.\n",x);
}

/*:28*/
#line 60 "dance.w"
;
exit(0);
}

/*:1*/
