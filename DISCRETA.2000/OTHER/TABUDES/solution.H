// solution.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_solution_h
#define PACKCOVER_solution_h

#include <iostream.h>
#include <multiset.h>
#include <iterator.h>
#include <algo.h>

#include "move.H"
#include "penalty.H"
#include "dp.H"
#include "matrix.H"

class Solution {
protected:
  vector<patchType> patches;
  volumeType currVolume;
  penaltyType currPenalty;

  void setCurrentVolume(volumeType vol) { currVolume = vol; }
  void setCurrentPenalty(penaltyType pen) { currPenalty = pen; }
  void addVolume(volumeType vol) { currVolume += vol; }
  void addPenalty(penaltyType pen) { currPenalty += pen; }

  void removePatch(int index) {
    patches[index] = patches[patches.size()-1];
    patches.resize(patches.size() - 1);
  }
  void removePatches(const vector<int>& indexes);
  void removePatches(const list<int>& indexes);
  void addPatch(const patchType patch) {
    if(patches.capacity() <= patches.size())
      patches.reserve(patches.size()*2);
    patches.resize(patches.size() + 1);
    patches[patches.size() - 1] = patch;
  }
  void changePatches(const Move& move);
  void addPatches(const list<patchType>& patches);
public:
  Solution(int initSize = 1) : currVolume(0), currPenalty(0)
    { patches.reserve(initSize); }
  //Solution(const Solution& orig);
  // memberwise copy constructor OK
  Solution& operator=(const Solution& orig) {
    patches = orig.patches;
    currVolume = orig.currVolume;
    currPenalty = orig.currPenalty;
    recalculateThings();
    return *this;
  }
  virtual ~Solution() {}
  virtual int numPatches() const { return patches.size(); }
  virtual bool empty() const { return patches.empty(); }
  virtual int getPatch(int index) const { return patches[index]; }

  virtual volumeType getCurrentVolume() const { return currVolume; }
  virtual penaltyType getCurrentPenalty() const { return currPenalty; }

  virtual void whatIfChange(Move& move, volumeType& volChange, 
			    penaltyType& penChange)
    { errorExit("Solution::whatIfChange called"); }
  virtual void recalculateThings() {} // nothing to recalculate here
  virtual void makeMove(Move& move) {
    volumeType volCh;
    penaltyType penCh;
    whatIfChange(move, volCh, penCh);
    makeMoveVP(move, volCh, penCh);
  }
  virtual void makeMoveVP(Move& move, volumeType volCh, penaltyType penCh) {
    addVolume(volCh);
    addPenalty(penCh);
    if(move.addCount() == move.removeCount())
      changePatches(move);
    else {
      removePatches(move.getRemoveIndexes());
      addPatches(move.getAdded());
    }
  }
  virtual void clear() { patches.clear(); setCurrentVolume(0); }
  virtual impactType impact(const Move& move) 
    { return move.removeCount() + move.addCount(); }

  virtual void initSolAtMostVol(const DataProvider& dp, 
				const volumeType maxvol);
  virtual void initSolAtLeastVol(const DataProvider& dp, 
				 const volumeType maxvol);
  virtual void initSolAtMostPen(const DataProvider& dp, 
				const penaltyType maxPen);
  virtual void initSolPartitions(const DataProvider& dp,
				 const int partSize);
  virtual void initSolVolumes(const DataProvider& dp, const string arg);
  // stream I/O (textual)

  virtual void print(ostream& out, int format = 0) const
  {
    copy(patches.begin(), patches.end(), 
	 ostream_iterator<patchType> (out, " "));
  }
  virtual void printDiscretaOutput(const string filename, 
				   const DataProvider *dp) const;
  virtual void printStatistics(ostream& out) const {}

  // file I/O (binary to conserve space)
  virtual void writeToFile(const string fileName);
  virtual void readFromFile(const string fileName);

  // check the volume and penalty of the (final) solution
  virtual void check();
};

ostream& operator<<(ostream&, const Solution&);

class CoverSolution : public Solution {
protected:
  const CoverMatrix& matrix;
  const CoveringDataProvider& dataProvider;
  CoverPenalty& penalty;
  vector<lambdaType> *covered;
  const vector<lambdaType>* coverGoal;
  int coverMax;
public:
  virtual void recalculateCoverings();
  virtual void recomputePenalty(int extraPatches = 0);
  virtual void recomputeVolume();
private:
  CoverSolution(const CoverSolution& orig) : 
    matrix(orig.matrix), dataProvider(orig.dataProvider), 
    penalty(orig.penalty), 
    covered(new vector<lambdaType>(*orig.covered)),
    coverGoal(orig.coverGoal), coverMax(orig.coverMax) {}
  virtual void CoverSolution::updateCovered
  (const vector<pointType>* points, const vector<lambdaType>* lambdas,
   const int addRem);
public:
  CoverSolution(const CoverMatrix& m, CoveringDataProvider& dp, 
		CoverPenalty& pen) : 
    matrix(m), dataProvider(dp), penalty(pen) { 
    covered = new vector<lambdaType>(dataProvider.totalPointCount(), 
				     lambdaType(0));
    coverGoal = dataProvider.coveringGoal();
    coverMax = *max_element((*coverGoal).begin(), (*coverGoal).end());
    recomputePenalty();
  }
  virtual ~CoverSolution() {
    delete covered;
  }

  virtual void whatIfChange(Move& move, volumeType& volChange, 
			    penaltyType& penChange);
  virtual void makeMoveVP(Move& move, volumeType volCh, penaltyType penCh);
  // slight speedup: write specialized makeMove(Move&)

  virtual void recalculateThings() {
    recomputeVolume();
    recalculateCoverings();
    recomputePenalty();
  }
  virtual void clear() {
    Solution::clear();
    fill(covered->begin(), covered->end(), lambdaType(0));
    recomputePenalty();
  }

  virtual void print(ostream& out, int format = 0) const {
    for(vector<patchType>::const_iterator I = patches.begin();
  	I != patches.end(); ++I != patches.end() ? out << ", " : out)
      dataProvider.printPatch(out, *I, format);
  }
  
  virtual void printStatistics(ostream& out) const {
    matrix.printStatistics(out);
  }

  virtual const lambdaType& coveringGoal(pointType point) const {
    return (*coverGoal)[point];
  }
  virtual const lambdaType& currentlyCovered(pointType point) const {
    return (*covered)[point];
  }
  virtual penaltyType pointPenalty(pointType point) const;
  virtual int numPoints() const { return (*covered).size(); }

  virtual void initSolAtLeastTotalCoverage(int reqCov);
};

#endif
