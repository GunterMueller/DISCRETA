/* perm2.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef PERMTRUE

#include <DISCRETA/ma.h> /* for PGL_2_p_as_perm() */
#include <DISCRETA/perm.h>
#include <DISCRETA/part.h>

#if TEXDOCU
INT PERMUTATION_OB::m_n_cycle(INT n)
#else
make an $n$-cycle. Computes the permutation $(1,2,3\ldots,n-1)$ into this.
#endif
{
	INT i;
	
	m_il(n);
	for (i = 0; i < n; i++) {
		if (i < n - 1)
			m_ii(i, i + 2);
		else
			m_ii(n - 1, 1);
		}
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::reordering_for_cyclic_permutation(PERMUTATION_OP q)
#endif
{
	PERMUTATION_OB p1;
	INT i, j, l;
	
	l = s_li();
	p1.m_il(l);
	j = 0;
	for (i = 0; i < l; i++) {
		p1.m_ii(i, j + 1);
		j = s_ii(j) - 1;
		}
	p1.invers(q);
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::reorder_to_cycles(PERMUTATION_OP q, VECTOR_OP v, INT f_v)
#endif
{
	INT i, j, l, o, ol, c, next, f;
	VECTOR_OB f_met;
	
	f = 0;
	o = 0;
	l = s_li();
	f_met.m_il_n(l);
	q->m_il(l);
	v->m_il(l);
	for (i = 0; i < l; i++) {
		if (f_met.s_ii(i))
			continue;
		ol = 0;
		c = i;
		while (TRUE) {
			q->m_ii(f, c + 1);
			f++;
			ol++;
			f_met.m_ii(c, 1);
			next = s_ii(c) - 1;
			if (f_met.s_ii(next) == 0) {
				c = next;
				}
			else {
				v->m_ii(o, ol);
				o++;
				break;
				}
			}
		}
	v->realloc_z(o);
	if (f_v) {
		printf("PERMUTATION_OB::reorder_to_cycles()\n");
		printf("this=");
		println();
		printf("q=");
		q->print_list();
		printf("v=");
		v->println();
		fflush(stdout);
		}
	q->invers_apply();
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::cyclic_multiplicator(INT n, INT a)
#else
computes the map $x \mapsto a \cdot x$ mod $n$.
#endif
{
	INT i, ai;
	
	m_il(n);
	printf("cyclic multiplicator n = %ld a = %ld\n", n, a);
	fflush(stdout);
	for (i = 0; i < n; i++) {
		ai = a * i % n;
		/* printf("i = %ld ai = %ld\n", i, ai);
		fflush(stdout); */
		m_ii(i, ai + 1);
		}
	/* println();
	fflush(stdout); */
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::number_of_fixpoints()
#else
Returns the 
number of fixpoints of the permutation this.
#endif
{
	INT i, l, k = 0;
	
	l = s_li();
	for (i = 0; i < l; i++) {
		if (s_ii(i) == i + 1)
			k++;
		}
	return k;
}

#ifdef PARTTRUE

#if TEXDOCU
INT PERMUTATION_OB::class_rep(PARTITION_OP type)
#else
The $S_n$ conjugacy classes are labelled by partitions. 
Given a partition (type), this routine computes a representative  
for the corresponding class, i.~e. 
an element in $S_n$ ($n$ is the weight of the partition) 
with cycle type $type$.
#endif
{
	INT w, erg = OK, i, j, k, t;
	
	if (type->s_k() == EXPONENT) {
		PARTITION_OB type1;
		
		type->t_EXPONENT_VECTOR(&type1);
		return class_rep(&type1);
		}
	
	w = type->weight_i();
	erg += m_il(w);
	k = 0;
	for (i = 0; i < type->s_li(); i++) {
		/* k ist naechste freie stelle */
		t = type->s_ii(i);
		for (j = 1; j < t; j++)
			m_ii(k + j - 1, k + j + 1);
		m_ii(k + t - 1, k + 1);
		k += t;
		}
	return OK;
}
#endif

#if TEXDOCU
INT PERMUTATION_OB::induce_action_on_columns(PERMUTATION_OP gg, MATRIX_OP I)
#endif
{
	VECTOR_OB B;
	PERMUTATION_OB q, qv, g, g1;
	INT v, b, f_v = FALSE;
	
	v = I->s_hi();
	b = I->s_li();
	I->calc_blocks(&B, &q, &qv, f_v);
	induce_action_on_blocks(&g, &B);
	q.mult(&g, &g1);
	g1.mult(&qv, gg);
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::induce_action_on_blocks(PERMUTATION_OP gg, VECTOR_OP B)
#else
Computes the induced action on the blocks of a design. 
this contains the point-permutation, B the sorted list 
of blocks of the simple (no repeated blocks) design. 
gg will contain the action of degree B$->$s\_li(); \\
Important: B is a sorted vector of sorted blocks !!!
This routine works fine only for \lq simple\rq designs, 
e.g. no block occurs twice 
exception: empty blocks are treated correctly, 
they lie in the first positions of B (due to the sorting).
#endif
{
	INT i, j, l, b, a, aa, ii;
	INT nb_empty, idx, f_found;
	VECTOR_OP bl;
	VECTOR_OB b1;
	
	b = B->s_li();
	gg->m_il(b);
	nb_empty = 0;
	for (i = 0; i < b; i++) {
		bl = (VECTOR_OP) B->s_i(i);
		l = bl->s_li();
		if (l == 0) {
			idx = nb_empty;
			nb_empty++;
			}
		else {
			b1.m_il(l);
			for (j = 0; j < l; j++) {
				a = bl->s_ii(j);
				aa = s_ii(a) - 1;
				b1.m_ii(j, aa);
				}
			b1.quicksort(l, TRUE /* f_ascending */);
			B->search(b, TRUE, &b1, &idx, &f_found);
			if (!f_found) {
				printf("perm: ");
				println();
				printf("block: ");
				printf("[");
				for (ii = 0; ii < l; ii++) {
					printf("%ld ", bl->s_ii(ii) + 1);
					}
				printf("]\n");
				// bl->println();
				printf("block image under perm: ");
				printf("[");
				for (ii = 0; ii < l; ii++) {
					printf("%ld ", b1.s_ii(ii) + 1);
					}
				printf("]\n");
				// b1.println();
				fflush(stdout);
				return error("PERM::induce_action_on_blocks(): not found");
				}
			idx--;
			}
		gg->m_ii(i, idx + 1);
		}
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::induce3(PERMUTATION_OP b)
#else
induction on three sets. Gives a permutation of degree $(n choose 3)$ 
where $n$ is the degree of the permutation in this.
#endif
{
	INT k, l, n, n3, m1, m2, m3, bm1, bm2, bm3, x;

	n = s_li();
	n3 = n * (n - 1) * (n - 2);
	n3 /= 2;
	n3 /= 3;
	b->m_il(n3);
	k = 0;
	for (m1 = 1; m1 <= n; m1++) {
		for (m2 = m1 + 1; m2 <= n; m2++) {
			for (m3 = m2 + 1; m3 <= n; m3++) {
				bm1 = s_ii(m1 - 1);
				bm2 = s_ii(m2 - 1);
				bm3 = s_ii(m3 - 1);
				if (bm1 > bm2) {
					x = bm1; bm1 = bm2; bm2 = x;
					}
				if (bm2 > bm3) {
					x = bm2; bm2 = bm3; bm3 = x;
					}
				if (bm1 > bm2) {
					x = bm1; bm1 = bm2; bm2 = x;
					}
				x = bm3 - bm2;
				for (l = bm1 + 1; l < bm2; l++)
					x += n - l;
				for (l = 1; l < bm1; l++) {
					if (n - l > 1)
						x += ((n - l) * (n - l - 1)) >> 1;
					}
				b->m_ii(k, x);
				k++;
				} // next m3
			} // next m2
		} // next m1
	if (k != n3)
		return error("PERM::induce3() k != n3");
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::induce2(PERMUTATION_OP b, INT n)
#else
a is in fact only a permutation of $1, .. n$. 
It computes the induced action of a on the pairs $(i,j)$, 
$1 \le i < j \le n$, which are enumerated in the following way:
$\{ \{1,2\}, \{1,3\}, \ldots,\{2,3\},\ldots,\{n-1,n\}\}$.
#endif
{
	INT i, j, k, i1, j1, k1, m;
	
	m = (n * (n - 1)) >> 1;
	b->m_il(m);
	for (i = 0; i < n - 1; i++) {
		i1 = s_ii(i) - 1;
		for (j = i + 1; j < n; j++) {
			j1 = s_ii(j) - 1;
			k = ij2k(i, j, n);
			k1 = ij2k(i1, j1, n);
			b->m_ii(k, k1 + 1);
			}
		}
	return OK;
}

#if TEXDOCU
INT support_on_2_sets(VECTOR_OP R, INT n, VECTOR_OP support, INT *supp)
#else
Some specialized routine needed for graphical designs. 
Given a set R of edges of a graph, this routine 
computes the support, i.e. the set of vertices 
v such that $(v,i)$ is an edge in the graph R.
#endif
// on n points
{
	INT i, k, l, a, b, su = 0;

	l = R->s_li();
	support->m_il_n(n);
	for (i = 0; i < l; i++) {
		k = R->s_ii(i);
		k2ij(k, &a, &b, n);
		support->m_ii(a, 1);
		support->m_ii(b, 1);
		}
	for (i = 0; i < n; i++) {
		if (support->s_ii(i))
			su++;
		}
	*supp = su;
	return OK;
}

#if TEXDOCU
INT tuple2_rank(INT rank, INT *i, INT *j, INT n, INT f_injective)
#else
enumeration of 2-tuples $(i,j)$ (f\_injective TRUE iff $i=j$ forbidden).
this routine produces the tuple with number ``rank'' into $i$ and $j$. 
$n$ is the number of points $1 \le i,j \le n$.
#endif
{
	INT a;
	
	if (f_injective) {
		a = rank % (n - 1);
		*i = (rank - a) / (n - 1);
		if (a <*i)
			*j = a;
		else
			*j = a + 1;
		}
	else {
		a = rank % n;
		*i = (rank - a) / n;
		*j = a;
		}
	return OK;
}

#if TEXDOCU
INT tuple2_unrank(INT i, INT j, INT *rank, INT n, INT f_injective)
#else
inverse function of tuple2\_rank(): 
returns the rank of a given tuple $(i,j)$.
#endif
{
	if (f_injective) {
		*rank = i * (n - 1);
		if (j < i)
			*rank += j;
		else if (j == i)
			return error("tuple2_unrank() not injective !");
		else
			*rank += j - 1;
		}
	else {
		*rank = i * n + j;
		}
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::induce_on_2tuples(PERMUTATION_OP p, INT f_injective)
#else
computes induction on two sets.
tuple2\_rank and tuple2\_unrank are used.
#endif
{
	INT n, m, i, j, rank, i1, j1, rank1;

	n = s_li();
	if (f_injective)
		m = n * (n - 1);
	else
		m = n * n;
	p->m_il(m);
	for (rank = 0; rank < m; rank++) {
		tuple2_rank(rank, &i, &j, n, f_injective);
		i1 = s_ii(i) - 1;
		j1 = s_ii(j) - 1;
		tuple2_unrank(i1, j1, &rank1, n, f_injective);
		p->m_ii(rank, rank1 + 1);
		}
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::induce2_01(PERMUTATION_OP b, INT n)
#else
Induced action on mappings $\{0,1\}^{{n \choose 2}}$.
#endif
{
	INT n2 = (n * (n - 1)) >> 1;
	INT m = (UINT) 1L << n2;
	INT i, j, k, k1;
	UINT g, g1;
	
	/* printf("induce2_01(): n = %ld s_li() = %ld "
	"n2 = %ld m = %ld\n", n, s_li(), n2, m);
	fflush(stdout); */
	b->m_il(m);
	for (g = 0; g < m; g++) {
		g1 = 0;
		for (k = 0; k < n2; k++) {
			if (g & ((UINT) 1L << k)) {
				k2ij(k, &i, &j, n);
				i = s_ii(i) - 1;
				j = s_ii(j) - 1;
				k1 = ij2k(i, j, n);
				g1 |= ((UINT) 1 << k1);
				}
			}
		b->m_ii(g, g1 + 1);
		}
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::join(PERMUTATION_OP b, PERMUTATION_OP c)
#else
Let this and $b$ act on disjoint sets, and put the resulting permutation into $c$.
#endif
{
	INT i, a, l1, l2, l3;

	l1 = s_li();
	l2 = b->s_li();
	l3 = l1 + l2;
	c->m_il(l3);
	for (i = 0; i < l1; i++) {
		a = s_ii(i);
		c->m_ii(i, a);
		}
	for (i = 0; i < l2; i++) {
		a = b->s_ii(i);
		c->m_ii(l1 + i, l1 + a);
		}
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::add_fixpoint_in_front(PERMUTATION_OP b)
#else
Adds a fixpoint as the new first point of the premutation. 
All other points are shifted up by one element.
#endif
{
	add_n_fixpoints_in_front(b, 1);
	return OK;
#if 0
	INT i, j, l;

	l = s_li();
	b->m_il(l + 1);
	b->m_ii(0, 1);
	for (i = 0; i < l; i++) {
		j = s_ii(i);
		b->m_ii(i + 1, j + 1);
		}
	return OK;
#endif
}

#if TEXDOCU
INT PERMUTATION_OB::embed_at(PERMUTATION_OP b, INT n, INT at)
#else
adds at fixpoints at the beginning, 
n - at - l fixpoints at the end. 
l is the length of the this permutation.
Result is b.
#endif
{
	PERMUTATION_OB q;
	INT l, m;

	l = s_li();
	if (n < l)
		return error("PERM::embed_at() n < l");
	if (at + l >= n)
		return error("PERM::embed_at() at + l >= n");
	m = n - at - l; // this is >= 0 !
	add_n_fixpoints_in_front(&q, at);
	q.add_n_fixpoints_at_end(b, m);
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::add_n_fixpoints_in_front(PERMUTATION_OP b, INT n)
#endif
{
	INT i, j, l;

	l = s_li();
	b->m_il(n + l);
	for (i = 0; i < n; i++)
		b->m_ii(i, i + 1);
	for (i = 0; i < l; i++) {
		j = s_ii(i);
		b->m_ii(n + i, n + j);
		}
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::add_n_fixpoints_at_end(PERMUTATION_OP b, INT n)
#endif
{
	INT i, j, l;

	l = s_li();
	b->m_il(n + l);
	for (i = 0; i < l; i++) {
		j = s_ii(i);
		b->m_ii(i, j);
		}
	for (i = 0; i < n; i++)
		b->m_ii(l + i, l + i + 1);
	return OK;
}

#if TEXDOCU
INT PERMUTATION_OB::remove_fixpoint(PERMUTATION_OP b, INT i)
#endif
{
	INT j, k, l;

	l = s_li();
	if (s_ii(i) - 1 != i)
		return error("PERM::remove_fixpoint(): i is not a fixpoint");
	b->m_il(l - 1);
	for (j = 0; j < l; j++) {
		if (j == i)
			continue;
		k = s_ii(j);
		if (k > i) // k > j + 1 also possible
			k--;
		if (j > i) // (j > i + 1 also possible because j == i + 1 never occures)
			b->m_ii(j - 1, k);
		else
			b->m_ii(j, k);
		}
	return OK;
}

#if TEXDOCU
INT bruhat_comp_perm(PERMUTATION_OP a, PERMUTATION_OP b)
#else
returns 1 if $a>b$,  $0$ if $a=b$,  $-1$ if $a<b$, NONCOMPARABLE otherwise.
#endif
{
        INT erg,erg2;
	
        erg = bru_comp(a,b);
        erg2 = bru_comp(b,a);
        if ((erg == TRUE) && (erg2 == TRUE)) return (INT) 0;
        if (erg == TRUE) return (INT ) 1;
        if ((erg == FALSE) && (erg2 == FALSE)) return NONCOMPARABLE;
        return (INT) -1;
}

#if TEXDOCU
INT bru_comp(PERMUTATION_OP a, PERMUTATION_OP c)
#else
returns TRUE  if $a\ge c$  in the Bruhat order 
AND condition when $c$ is not long enough.
#endif
{
        INT i, j, k, x, y1, y2;
	
	k = a->s_li();
	y1 = a->s_ii(0);
	y2 = a->s_ii(k-1);
        if (c->s_ii(0) > y1 ) return (FALSE);

        if ( k < c->s_li() ) {
                for (j=k; j < c->s_li(); j++)
                        if (j != c->s_ii(j) - 1) return FALSE;
                }
        if ( (c->s_li() == k) && (c->s_ii(k-1) < y2)) return (FALSE);


        if (c->s_li() < k) k = c->s_li();

        for (i=0; i<k; i++)  {
                x = 0;
                for (j=0; j<k; j++){
                        if (a->s_ii(j) >i ) x++;
                        if (c->s_ii(j) >i ) x--;
                        if (x<0)  return (FALSE);
                }
        }

        return (TRUE);
}

#if TEXDOCU
INT fastrectr(PERMUTATION_OP a, VECTOR_OP v)
#else
Alain Lascoux, 3/1997
input: permutation $a \in S_n$
output v a vector of rectrices, which are vectors of length 3 over integers.
the sum of the elements (plus an additional fourth element) 
is constant (equal to n).
#endif
{
	PERMUTATION_OB b;
	VECTOR_OB u;
	INT i, k, x, y, z, iv, i1;
	
	a->invers(&b);
	u.m_il(3);
	v->m_il(0);
	iv = 0;
	for (i = 0; i < a->s_li() - 1; i++) {
		if (a->s_ii(i) > a->s_ii(i + 1)) {
			z = a->s_ii(i);
			x = a->s_ii(i + 1);
			for (k = z; k >= x; k--) {
				if  (b.s_ii(k-1) >= i + 2 && b.s_ii(k) <= i + 1) {
					y = 0;
					for (i1 = 0; i1 <= i; i1++) { 
						if (a->s_ii(i1) < k) 
							y++;
						} 
					u.m_ii(0, y);
					u.m_ii(1, i + 1 - y);
					u.m_ii(2, k - y);
					v->inc();
					u.copy((VECTOR_OP) v->s_i(iv));
					iv++;
					} // if
				} // next k
			} // if 
		} // next i
	return OK;
}

#if 0
INT fastrectr(a,v)  OP a,v;
{  OP b,u; INT i,j,k,x,y,z,iv,i1;  b=callocobject();u=callocobject();

invers(a,b);  init(VECTOR,v);m_il_v(3L,u);iv=0L;
for(i=0L;i<S_P_LI(a)-1L;i++)
  { if( S_P_II(a,i)>S_P_II(a,i+1))
    {z= S_P_II(a,i); x=S_P_II(a,i+1);
      for (k=z;k>=x;k--)
       {if  ( S_P_II(b,k-1) >= i+2 && S_P_II(b,k) <=i+1)
         { y=0; for(i1=0;i1<=i;i1++) { if( S_P_II(a,i1) <k) y++;}
              M_I_I(y,S_V_I(u,0L));
               M_I_I(i+1-y,S_V_I(u,1L));
                   M_I_I(k-y,S_V_I(u,2L));
   inc(v);copy(u,S_V_I(v,iv)); iv++;     }}} }
                          freeall(b);freeall(u);
}
#
# Rectr ## rectrices pour le groupe sym
Rectr:=proc(perm)
  local
    i, i1, k, nvperm, n, pi, pii, y, lr;
  lr:=NULL;
  invperm:=`SG/InvPerm`(perm);
  n:=nops(perm);
  for i from 0 to n-2 do
    if (perm[i+1]>perm[i+2]) then
      pi:=perm[i+1];
      pii:=perm[i+2];
      for k from pi by -1 to pii do
        if ((invperm[k]>=i+2) and (invperm[k+1]<=i+1)) then
          y:=0;
          for i1 from 0 to i do
            if (perm[i1+1]<k) then
              y:=y+1
            fi
          od; 
          lr:=lr,[y, i+1-y, k-y]
        fi;
      od;  
    fi;  
  od;
  RETURN([lr])
end;
#endif

#endif /* PERMTRUE */


