// tenure.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_tenure_h
#define PACKCOVER_tenure_h

#include <vector.h>
#include <deque.h>
#include <map.h>

#include "utils.H"

class Tenure {
public:
  Tenure() {};
  virtual int time() const = 0;
  virtual void newPenalty(penaltyType) {}
};

class ConstantTenure : public Tenure {
private:
  const int tenureTime;
public:
  ConstantTenure(int t) : tenureTime(t) {}
  virtual int time() const { return tenureTime; }
};

class CompareWindows {
  // do a Wilcoxon rank sum test for two successive sliding windows of
  // equal size. If the test cannot show difference, assume that
  // windows are ewual (not very good asumption always...)
private:
  const int windowLength;
  vector<deque<penaltyType> > sample;
  // should use a queue instead of deque, but it is not possible to
  // print out a queue since it does not provide iterators! This makes
  // debugging cumbersome.
  map<pair<penaltyType, int>, int> counts;
  typedef typename map<pair<penaltyType, int>, int>::iterator count_iterator;
public:
  CompareWindows(int wl);
  bool changed();
  void add(penaltyType pen);
};

class DynamicTenure : public Tenure {
private:
  const int windowLength;
  const vector<int> tenures;
  int currTenure, tenIx;
  int refreshWindow;
  CompareWindows *win;
public:
  DynamicTenure(vector<int> tens, int wl) : windowLength(wl), tenures(tens) {
    if(windowLength < 3) errorExit("window length must be at least 3");
    win = new CompareWindows(windowLength);
    tenIx = 0;
    currTenure = tenures[tenIx];
    refreshWindow = windowLength;
  }
  virtual ~DynamicTenure() { delete win; }
  virtual int time() const { return currTenure; }
  virtual void newPenalty(penaltyType penalty);
  bool stagnated() { return !win->changed(); }
  void registerPenalty(penaltyType penalty);
};

class ShiftTenure : public Tenure {
private:
  const int changeLength;
  const vector<int> tenures;
  int currTenure, tenIx;
  int count;
public:
  ShiftTenure(vector<int> tens, int wl) : changeLength(wl), tenures(tens) {
    tenIx = 0;
    currTenure = tenures[tenIx];
    count = 0;
  }
  virtual ~ShiftTenure() {}
  virtual int time() const { return currTenure; }
  virtual void newPenalty(penaltyType) {
    if(++count >= changeLength) {
      if(verbose)
	cout << "changing tenure from " << currTenure;
      count = 0;
      tenIx = ++tenIx % tenures.size();
      currTenure = tenures[tenIx];
      if(verbose)
	cout << " to " << currTenure << '\n';
    }
  }
};

#endif
