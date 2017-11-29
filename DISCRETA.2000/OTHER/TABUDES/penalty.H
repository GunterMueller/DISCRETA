// penalty.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_penalty_h
#define PACKCOVER_penalty_h

#include <vector.h>
#include <math.h>
#include "utils.H"

class CoverPenalty {
public:
  virtual bool isUniform() const { return false; }
  virtual penaltyType penalty(const pointType point, 
			      const lambdaType coverTimes,
			      const lambdaType coverGoal) const = 0;
  virtual penaltyType* penaltyTable(const int, const int) {
    errorExit("penaltyTable() called for non-uniform penalty");
    return 0;
  }
};

class UniformCoverPenalty : public CoverPenalty {
protected:
  penaltyType* table;
  int currMin, currMax;
public:
  UniformCoverPenalty() : table(0), currMin(0), currMax(0) {}
  virtual ~UniformCoverPenalty() { delete[] table; }
  virtual bool isUniform() const { return true; }
  virtual penaltyType penalty(const pointType point, 
			      const lambdaType coverTimes,
			      const lambdaType coverGoal) const {
    errorExit("penalty() called for uniform penalty class");
    return penaltyType(0);
  }
  virtual void checkBounds(const int min, const int max) const {
    if(min >= max || min >= 0 || max < 0)
      errorExit("impossible bounds for UniformCoverPenalty");
  }
  virtual void allocTable(const int min, const int max) {
    checkBounds(min, max);
    delete[] table;
    table = new penaltyType[max-min+1];
    for(int i = 0; i < max-min+1; i++)
      table[i] = penaltyType(0);
    currMin = min;
    currMax = max;
  }
  virtual penaltyType* newPenaltyTable(const int min, const int max) = 0;
  virtual penaltyType* penaltyTable(const int min, const int max) {
    checkBounds(min, max);
    if(min < currMin || max > currMax)
      return newPenaltyTable(min, max);
    else
      return table - currMin;
  }
};

class CoveringPenalty : public UniformCoverPenalty {
public:
  virtual ~CoveringPenalty() {}
  virtual penaltyType* newPenaltyTable(const int min, const int max);
};

class PackingPenalty : public UniformCoverPenalty {
public:
  virtual ~PackingPenalty() {}
  virtual penaltyType* newPenaltyTable(const int min, const int max);
};

class DesignPenalty : public UniformCoverPenalty {
public:
  virtual ~DesignPenalty() {}
  virtual penaltyType* newPenaltyTable(const int min, const int max);
};

class AsymHighDesignPenalty : public UniformCoverPenalty {
private:
  const penaltyType padd;
public:
  AsymHighDesignPenalty(penaltyType pa) : padd(pa) {}
  virtual ~AsymHighDesignPenalty() {}
  virtual penaltyType* newPenaltyTable(const int min, const int max);
};

class AsymLowDesignPenalty : public UniformCoverPenalty {
private:
  const penaltyType padd;
public:
  AsymLowDesignPenalty(penaltyType pa) : padd(pa) {}
  virtual ~AsymLowDesignPenalty() {}
  virtual penaltyType* newPenaltyTable(const int min, const int max);
};

#endif
