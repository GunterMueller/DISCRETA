// utils.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_utils_h
#define PACKCOVER_utils_h

#include <string>
#include <stack.h>
#include <list.h>
#include "checked.H"

// debug etc. routines

void errorExit(const string& msg);

int binCoef(int a, int b);
Checked<int> binCoef(Checked<int> a, Checked<int> b);

#define maxv 1001

class BinCoefTable {
private:
  int binCoef[maxv + 1][maxv + 2];
public:
  BinCoefTable();
  int value(int n, int m) { return binCoef[n][m]; }
};

extern BinCoefTable binCoefTable;

typedef unsigned short varietyType;
typedef int rankType;
#define solX 0

rankType rankSubset(varietyType *subset, int card);
void getFirstSubset(varietyType *subset, int card);
int getNextSubset(varietyType *subset, int card, int v);
void unrankSubset(rankType rank, varietyType *subset, int card);
void makeComplement(varietyType *s, varietyType *c, int v);

class CharWeight {
private:
  int *weights;
public:
  CharWeight();
  ~CharWeight() { delete[] weights; }
  weight(unsigned char c) { return weights[c]; }
};

extern CharWeight charWeight;

template <class T>class PopTopStack : public stack<T> {
public:
  T popTop() {
    if(empty()) errorExit("parse error: stack empty");
    T tmp = top(); pop(); return tmp; 
  }
};

template<class T> inline vector<T> makeVector(const T& a)
{
  vector<T> vec(1);
  vec[0] = a;
  return vec;
};

template<class T> inline vector<T> makeVector(const T& a,const T& b)
{
  vector<T> vec(2);
  vec[0] = a; vec[1] = b;
  return vec;
}

template<class T> inline vector<T> makeVector(const T& a,const T& b,const T& c)
{
  vector<T> vec(3);
  vec[0] = a; vec[1] = b; vec[2] = c;
  return vec;
}

template<class T> inline vector<T> makeVector(const T& a,const T& b,const T& c,
				       const T& d)
{
  vector<T> vec(4);
  vec[0] = a; vec[1] = b; vec[2] = c; vec[3] = d;
  return vec;
}

template<class T> inline vector<T> makeVector(const T& a,const T& b,const T& c,
				       const T& d,const T& e)
{
  vector<T> vec(5);
  vec[0] = a; vec[1] = b; vec[2] = c; vec[3] = d; vec[4] = e;
  return vec;
}

template<class T> 
inline ostream& operator<<(ostream& out, vector<T> v) {
  out << "vector: ";
  copy(v.begin(), v.end(), ostream_iterator<T>(out, " "));
  return out;
}

template<class T> 
inline ostream& operator<<(ostream& out, list<T> v) {
  out << "list: ";
  copy(v.begin(), v.end(), ostream_iterator<T>(out, " "));
  return out;
}

template<class T>
inline ostream& operator<<(ostream& out, deque<T> q) {
  out << "deque: ";
  copy(q.begin(), q.end(), ostream_iterator<T>(out, " "));
  return out;
}

#endif
