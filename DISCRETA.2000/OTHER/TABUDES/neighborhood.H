// neighborhood.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_neighborhood_h
#define PACKCOVER_neighborhood_h

#include <set.h>
#include <hash_set.h>
#include "limits.h"
#include "attr.H"
#include "warning.H"
#include "dp.H"
#include "solution.H"

// Neighborhood class calling sequence:

// 1: constructor
// 2: [setNoLimits + limitAdd/limitRemove]
// 3: inititer()/initRandom()/initNextRandom()
// 4: a) or b)
// 4 a): findNextMove() (repeat n times)
// 4 b): findRandomMove() (repeat n times)
// 4 c): findNextRandomMove() (repeat n times)
// go back to 2.

// makeAttr and makeTestAttr can be called any time after the first
// initIter, but they refer to the solution (if solution is needed) at
// the time of the last initIter

class Neighborhood {
private:
  const int attrId;
public:
  Neighborhood() : attrId(AttrIdProvider::getId()) {}
  int getAttrId() const { return attrId; }
  virtual ~Neighborhood() {}
  virtual void initIter() = 0;
  virtual void initRandom() { initIter(); }
  virtual bool findNextMove(Move& move) = 0;
  virtual bool findRandomMove(Move& move) {
    // if no random neighborhood defined, then use standard
    // neighborhood
    if(findNextMove(move))
      return true;
    else { // start over
      initIter();
      return findNextMove(move);
    }
  }

  virtual void initNextRandom() {
    errorExit("initNextRandom not implemented for this neighborhood");
  }
  virtual bool findNextRandomMove(Move& move) {
    errorExit("findNextRandom not implemented for this neighborhood");
    return false;
  }

  // makeAttr and testAttr give the attributes of the move specified
  // by the parameters (see Attribute class in tabulist.H). They
  // should include getAttrId() as the last element of each attribute
  // vector

  virtual list<vector<int> > makeAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty)
    const { return list <vector<int> > (1, makeVector(attrId)); }
  virtual list<vector<int> > makeTestAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty)
    const { return makeAttr(m, sol, newVolume, newPenalty); }

  virtual int neighborhoodSize() const { return -1; }
  virtual bool isSymmetric() const { return false; }

  // After call of limitRemove or limitAdd, only those moves that
  // remove (or add) the specified patches are allowed. Each call of
  // limitRemove or discards any previous removal limitations
  // (similarly for limitAdd). Not all neighborhoods support
  // limitations. setNoLimits() discards all limits. It is not allowed
  // to set an empty limit vector.

  virtual void limitRemove(const vector<int> alwRemIx) {
    errorExit("Neighborhood::limitRemove() not implemented");
  }
  virtual void limitAdd(const vector<patchType> alwAddPatches) {
    errorExit("Neighborhood::limitAdd() not implemented");
  }
  virtual void setNoLimits() {
    warning->print("Neighborhood::setNoLimits() called (no effect)");
  }
};

class EmptyNeighborhood : public Neighborhood {
public:
  EmptyNeighborhood() {}
  virtual void initIter() {}
  virtual bool findNextMove(Move&) { return false; }
  virtual bool findRandomMove(Move&) { return false; }
  virtual int neighborhoodSize() const { return 0; }
  virtual bool isSymmetric() const { return true; }
  virtual void limitRemove(const vector<int>) {}
  virtual void limitAdd(const vector<patchType>) {}
  virtual void setNoLimits() {}
};

class AddManyNeighborhood : public Neighborhood {
private:
  const DataProvider& dataProvider;
  int N;
  vector<int> changeIxs;
  int nextChange;
  vector<patchType> addVect;
  int patchCount() const { return addVect.empty() ? 
			     dataProvider.totalPatchCount() : addVect.size(); }
  vector<patchType> nextRandom;
  int nextRandomCount;
public:
  AddManyNeighborhood(const DataProvider& dp, int n = 1) : 
    dataProvider(dp), N(n) { initIter(); }
  virtual void initIter() { 
    changeIxs.resize(0); changeIxs.resize(N, patchType(0)); 
    changeIxs[0] = -1; nextChange = 0;
  }
  virtual bool findNextMove(Move& move);
  virtual bool findRandomMove(Move& move);
  virtual void setCount(int n) {
    if(n <= 0) errorExit("AddManyNeighborhood::setCount()");
    N = n; initIter();
  }
  virtual int neighborhoodSize() const {
    return binCoef(Checked<int>(patchCount()) +
		   Checked<int>(N) - Checked<int>(1),
		   Checked<int>(N))();
  }
  virtual void limitRemove(const vector<int> alwRemIx) {}
  virtual void limitAdd(const vector<int> alwAddPatches) {
    addVect = alwAddPatches;
  }
  virtual void setNoLimits() { addVect.clear(); }
  virtual void initNextRandom();
  virtual bool findNextRandomMove(Move& move);
};

class AddManyUniqueNeighborhood : public Neighborhood {
private:
  const DataProvider& dataProvider;
  int N;
  vector<patchType> patches;
  int nextChange;
  bool justInitialized;
  vector<patchType> addVect;
  int patchCount() const { return addVect.empty() ? 
			     dataProvider.totalPatchCount() : addVect.size(); }
public:
  AddManyUniqueNeighborhood(const DataProvider& dp, int n = 1) : 
    dataProvider(dp), N(n) { initIter(); }
  virtual void initIter() { 
    patches.resize(0); patches.resize(N);
    justInitialized = true;
  }
  virtual bool findNextMove(Move& move);
  virtual bool findRandomMove(Move& move);
  virtual void setCount(int n) {
    if(n <= 0) 
      errorExit("AddManyNeighborhood::setCount()");
    N = n; initIter();
  }
  virtual int neighborhoodSize() const {
    return binCoef(Checked<int>(patchCount()),
		   Checked<int>(N))();
  }
  virtual void limitRemove(const vector<int> alwRemIx) {}
  virtual void limitAdd(const vector<int> alwAddPatches) {
    addVect = alwAddPatches;
  }
  virtual void setNoLimits() { addVect.clear(); }
};

class RemoveManyNeighborhood : public Neighborhood {
private:
  const Solution& solution;
  int N;
  vector<int> changeIxs;
  bool justInitialized;
  vector<int> remVect;
  int ixCount() const { 
    return remVect.empty() ? solution.numPatches() : remVect.size();
  }
public:
  RemoveManyNeighborhood(const Solution& sol, int n = 1) : 
    solution(sol), N(n) { initIter(); }
  virtual void initIter() {
    justInitialized = true;
    changeIxs.resize(N);
  }
  virtual bool findNextMove(Move& move);
  virtual bool findRandomMove(Move& move);
  virtual void setCount(int n) {
    if(n <= 0) errorExit("RemoveManyNeighborhood::setCount()");
    N = n; initIter();
  }
  virtual int neighborhoodSize() const {
    return ixCount() ? binCoef(Checked<int>(ixCount()),
			       Checked<int>(N))() : 0;
  }
  virtual void limitAdd(const vector<patchType> alwAddPatches) {}
  virtual void limitRemove(const vector<int> alwRemIx) {
    remVect = alwRemIx;
  }
  virtual void setNoLimits() { remVect.clear(); }
};

class CompositeNeighborhood : public Neighborhood {
  // tries all combinations of single moves from each neighborhood of
  // nbs
private:
  vector<Neighborhood*> nbs;
  vector<Move> moves;
  bool firstTry;
public:
  CompositeNeighborhood() {}
  virtual void addNeighborhood(Neighborhood *n) {
    nbs.push_back(n);
    moves.resize(nbs.size());
  }
  virtual void initIter();
  virtual bool findNextMove(Move& move);
  virtual bool findRandomMove(Move& move);
  virtual int neighborhoodSize() const;
  virtual list<vector<int> > makeAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty)const{
    // this may or may not work, depending on the implementation of
    // makeAttr for each neighborhood of nbs
    list<vector<int> > tmp, tmp2;
    for(vector<Neighborhood*>::const_iterator I = nbs.begin(); 
	I != nbs.end(); ++I) {
      tmp2 = (*I)->makeAttr(m, sol, newVolume, newPenalty);
      tmp.splice(tmp.end(), tmp2);
    }
    return tmp;
  }
  virtual list<vector<int> > makeTestAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty)const{
    // this may or may not work, depending on the implementation of
    // makeTestAttr for each neighborhood of nbs
    list<vector<int> > tmp, tmp2;
    for(vector<Neighborhood*>::const_iterator I = nbs.begin(); 
	I != nbs.end(); ++I) {
      tmp2 = (*I)->makeTestAttr(m, sol, newVolume, newPenalty);
      tmp.splice(tmp.end(), tmp2);
    }
    return tmp;
  }    
  virtual bool isSymmetric() const {
    warning->print("don't know if CompositeNeighborhood is symmetric");
    return false; 
  }
  virtual void limitRemove(const vector<int> alwRemIx) {
    for(unsigned i = 0; i < nbs.size(); ++i)
      nbs[i]->limitRemove(alwRemIx);
  }
  virtual void limitAdd(const vector<patchType> alwAddPatches) {
    for(unsigned i = 0; i < nbs.size(); ++i)
      nbs[i]->limitAdd(alwAddPatches);
  }
  virtual void setNoLimits() {
    for(unsigned i = 0; i < nbs.size(); ++i)
      nbs[i]->setNoLimits();
  }
};

class ChangeManyNeighborhood : public CompositeNeighborhood {
private:
  AddManyNeighborhood *an;
  RemoveManyNeighborhood *rn;
  ChangeManyNeighborhood() {}
  ChangeManyNeighborhood(const ChangeManyNeighborhood&) {}
  ChangeManyNeighborhood& operator=(const ChangeManyNeighborhood&)
    { return *this; }
public:
  ChangeManyNeighborhood(const Solution& sol, const DataProvider& dp, 
			 int n = 1) {
    addNeighborhood(an = new AddManyNeighborhood(dp, n));
    addNeighborhood(rn = new RemoveManyNeighborhood(sol, n));
  }
  virtual ~ChangeManyNeighborhood() {
    delete an;
    delete rn;
  }
  virtual bool isSymmetric() const { return true; }
};

class ChangeManyUniqueNeighborhood : public CompositeNeighborhood {
private:
  AddManyUniqueNeighborhood *an;
  RemoveManyNeighborhood *rn;
  ChangeManyUniqueNeighborhood() {}
  ChangeManyUniqueNeighborhood(const ChangeManyUniqueNeighborhood&) {}
  ChangeManyUniqueNeighborhood& operator=(const ChangeManyUniqueNeighborhood&)
    { return *this; }
public:
  ChangeManyUniqueNeighborhood(const Solution& sol, const DataProvider& dp, 
			       int n = 1) {
    addNeighborhood(an = new AddManyUniqueNeighborhood(dp, n));
    addNeighborhood(rn = new RemoveManyNeighborhood(sol, n));
  }
  virtual ~ChangeManyUniqueNeighborhood() {
    delete an;
    delete rn;
  }
  virtual bool isSymmetric() const { return false; }
};

class VChangeOneNeighborhood : public Neighborhood {
private:
  const Solution& solution;
  const DataProvider& dataProvider;
  int changeIx, changeTo;
  unsigned iterIx;
  bool firstIter;
  const vector<patchType>* vp;
  vector<int> remVect;
  unsigned chix;
  hash_set<patchType> addSet;
public:
  VChangeOneNeighborhood(const Solution& sol, const DataProvider& dp) :
    solution(sol), dataProvider(dp) {}
  virtual void initIter() { firstIter = true; }
  virtual bool findNextMove(Move& move);
  virtual bool findRandomMove(Move& move);
  virtual int neighborhoodSize() const;
  virtual bool isSymmetric() const { return true; }
  virtual void limitRemove(const vector<int> alwRemIx) {
    remVect = alwRemIx;
  }
  virtual void limitAdd(const vector<patchType> alwAddPatches) {
    addSet.clear();
    for(vector<patchType>::const_iterator I = alwAddPatches.begin();
	I != alwAddPatches.end(); ++I)
      addSet.insert(*I);
  }
  virtual void setNoLimits() {
    remVect.clear();
    addSet.clear();
  }
};

class AdditiveNeighborhood : public Neighborhood {
  // tries the union of the moves of the neighborhoods in nbs
private:
  vector<Neighborhood*> nbs;
  unsigned currNb;
  int nbSize;
public:
  AdditiveNeighborhood() {}
  virtual void addNeighborhood(Neighborhood *n) {
    nbs.push_back(n);
  }
  virtual void initIter();
  virtual bool findNextMove(Move& move);
  virtual bool findRandomMove(Move& move);
  virtual int neighborhoodSize() const;
  virtual list<vector<int> > makeAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty)const{
    errorExit("AdditiveNeighborhood::makeAttr() called");
    return list<vector<int> >();
  }
  virtual list<vector<int> > makeTestAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty)const{
    errorExit("AdditiveNeighborhood::makeTestAttr() called");
    return list<vector<int> >();
  }
  virtual bool isSymmetric() const {
    bool symm = true;
    for(unsigned i = 0; i < nbs.size() && symm; ++i)
      symm = symm && nbs[i]->isSymmetric();
    return symm;
  }
  virtual void limitRemove(const vector<int> alwRemIx) {
    for(unsigned i = 0; i < nbs.size(); ++i)
      nbs[i]->limitRemove(alwRemIx);
  }
  virtual void limitAdd(const vector<patchType> alwAddPatches) {
    for(unsigned i = 0; i < nbs.size(); ++i)
      nbs[i]->limitAdd(alwAddPatches);
  }
  virtual void setNoLimits() {
    for(unsigned i = 0; i < nbs.size(); ++i)
      nbs[i]->setNoLimits();
  }
};

class CorrectCoveringNeighborhood : public Neighborhood {
protected:
  const CoverSolution& solution;
private:
  Neighborhood *incCoverNeigh, *decCoverNeigh;
  const CoverMatrix& matrix;
protected:
  int checkPoint;
private:
  Neighborhood* currNeighborhood;
protected:
  virtual void nextCheckPoint();
public:
  CorrectCoveringNeighborhood
  (const CoverSolution& sol, Neighborhood* incCN, Neighborhood* decCN,
   const CoverMatrix& mat) :
    solution(sol), incCoverNeigh(incCN), decCoverNeigh(decCN),
    matrix(mat), checkPoint(-1), currNeighborhood(0) {}
  virtual void initIter();
  virtual bool findNextMove(Move& move);
  virtual bool findRandomMove(Move& move);
  virtual int neighborhoodSize() const;
};

class CorrectRandomCoveringNeighborhood : public CorrectCoveringNeighborhood {
protected:
  virtual void nextCheckPoint();
public:
  CorrectRandomCoveringNeighborhood
  (const CoverSolution& sol, Neighborhood* incCN, Neighborhood* decCN,
   const CoverMatrix& mat) : 
    CorrectCoveringNeighborhood(sol, incCN, decCN, mat) {}
};

class LargeCorrectCoveringNeighborhood : public AdditiveNeighborhood {
private:
  const CoverSolution& solution;
  Neighborhood *incCoverNeigh, *decCoverNeigh;
  const CoverMatrix& matrix;
  bool useIncCover, useDecCover;
public:
  LargeCorrectCoveringNeighborhood
  (const CoverSolution& sol, Neighborhood* incCN, Neighborhood* decCN,
   const CoverMatrix& mat) :
    solution(sol), incCoverNeigh(incCN), decCoverNeigh(decCN),
    matrix(mat) {
    addNeighborhood(incCoverNeigh);
    addNeighborhood(decCoverNeigh);
  }
  virtual void initIter();
  virtual bool findNextMove(Move& move);
  virtual bool findRandomMove(Move& move);
  virtual int neighborhoodSize() const;
  virtual void initNextRandom();
  virtual bool findNextRandomMove(Move& move);
};

class BalanceCoveringNeighborhood : public CompositeNeighborhood {
protected:
  const CoverSolution& solution;
  Neighborhood *const incCoverNeigh, *const decCoverNeigh;
  const CoverMatrix& matrix;
  int addCheck, decCheck;
public:
  BalanceCoveringNeighborhood
  (const CoverSolution& sol, Neighborhood* incCN, Neighborhood* decCN, 
   const CoverMatrix& mat) :
    CompositeNeighborhood(),
    solution(sol), incCoverNeigh(incCN), decCoverNeigh(decCN), matrix(mat),
    addCheck(-1), decCheck(-1)
    {
      addNeighborhood(incCoverNeigh);
      addNeighborhood(decCoverNeigh);
    }
  virtual void initIter();
};

class LargeBalanceCoveringNeighborhood : public BalanceCoveringNeighborhood {
private:
  bool useIncCover, useDecCover;
public:
  LargeBalanceCoveringNeighborhood
  (const CoverSolution& sol, Neighborhood* incCN, Neighborhood* decCN, 
   const CoverMatrix& mat) : 
    BalanceCoveringNeighborhood(sol, incCN, decCN, mat) {}
  virtual void initIter();
  virtual bool findNextMove(Move& move);
  virtual bool findRandomMove(Move& move);
  virtual int neighborhoodSize() const;
};

class MultiNeighborhood : public Neighborhood {
private:
  vector<Neighborhood*> neighs;
  const int goalSize;
  Neighborhood *currNPtr;
public:
  MultiNeighborhood(int gs) : goalSize(gs), currNPtr(0) {}
  virtual void addNeighborhood(Neighborhood *const n) {
    neighs.push_back(n);
  }
  virtual void initIter();
  virtual bool findNextMove(Move& move) 
    { return currNPtr->findNextMove(move); }
  virtual bool findRandomMove(Move& move) 
    { return currNPtr->findRandomMove(move); }
  virtual int neighborhoodSize() const
   { return currNPtr->neighborhoodSize(); }
};

#endif
