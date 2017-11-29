// utils.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#include <iostream.h>

#include "packcover.H"
#include "utils.H"
#include "limits.H"

void errorExit(const string& msg)
{
  cout << msg << "\n";
  abort();
}

int binCoef(int a, int b)
{
  if(b > a) return 0;
  int res = 1, al = a - b, bl = 0;
  while(bl != b)
    res = res * ++al / ++bl;
  return res;
}

Checked<int> binCoef(Checked<int> a, Checked<int> b)
{
  if(a <= 0 || b < 0)
    errorExit("invalid parameters: binCoef");
  if(b > a) return Checked<int>(0);
  Checked<int> res = 1, al = a - b, bl = 0;
  while(bl != b)
    res = res * ++al / ++bl;
  return res;
}

BinCoefTable binCoefTable;

BinCoefTable::BinCoefTable()
{
  int v,k;

  for(v = 0; v <= maxv; v++) {
    binCoef[v][0] = binCoef[v][v] = 1;
    binCoef[v][v+1] = 0;
    for(k = 1; k <= v - 1; k++) {
      binCoef[v][k] = binCoef[v - 1][k - 1] + binCoef[v - 1][k];
      if(binCoef[v][k] < binCoef[v - 1][k - 1] ||
	 binCoef[v][k] < binCoef[v - 1][k] ||
	 binCoef[v - 1][k - 1] == 0 ||
	 (binCoef[v - 1][k] == 0 && k < v))
	binCoef[v][k] = 0; /* there was an overflow */
    }
  }
}

rankType rankSubset(varietyType *subset, int card)
{
  int i;
  rankType rank = 0;

  for(i = 0; i < card; i++)
    rank += binCoefTable.value(*subset++, i + 1);
  return rank;
}


/*
** `getFirstSubset' gets the subset with rank 0 and the given cardinality.
**
*/

void getFirstSubset(varietyType *subset, int card)
{
  int i;

  for(i = 0; i < card; i++)
    *subset++ = i;
  *subset = maxv + 1; /* sentinel for getNextSubset() */
}


/*
** `getNextSubset' gets the subset that has rank one bigger than the rank
** of the given subset. `v' is the number of varieties. Varieties range
** from 0 to v - 1. It returns 0, if this was the greatest possible rank
** for this type subset, otherwise the return value is 1.
*/

int getNextSubset(varietyType *subset, int card, int v)
{
  int i,j;

  if(subset[0] >= v - card)
    return 0;
  else {
    j = 0;
    while(subset[j + 1] <= subset[j] + 1)
      j++;
    subset[j]++;
    for(i = 0; i < j; i++)
      subset[i] = i;
    return 1;
  }
}


/*
** `unrankSubset' makes a subset out of a rank.
**
*/

void unrankSubset(rankType rank, varietyType *subset, int card)
{
  int p, m, i;

  m = rank;
  for(i = card - 1; i >= 0; i--) {
    p = i;
    while(binCoefTable.value(p + 1, i + 1) <= m)
      p++;
    m -= binCoefTable.value(p, i + 1);
    subset[i] = p;
  }
}


/*
** `makeComplement' forms the complement of a given subset `s' to
** space pointed by `c'. The number of varieties is `v'.
**
*/

void makeComplement(varietyType *s, varietyType *c, int v)
{
  int i;

  for(i = 0; i < v; i++)
    if(*s == (varietyType) i)
      s++;
    else
      *c++ = (varietyType) i;
  *c = maxv + 1; /* sentinel */
}

CharWeight::CharWeight()
{
  const unsigned charMax = numeric_limits<unsigned char>::max();
  weights = new int[charMax+1];
  fill(weights, weights + charMax + 1, 0);
  for(unsigned i = 0; i <= charMax; ++i) {
    unsigned test = i;
    while(test) {
      if(test & 1) weights[i]++;
      test >>= 1;
    }
  }
}

CharWeight charWeight;

