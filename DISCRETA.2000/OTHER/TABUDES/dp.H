// dp.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_dp_h
#define PACKCOVER_dp_h

#include "warning.H"
#include "utils.H"
#include <hash_map.h>

class DataProvider {
protected:
  vector<volumeType> volumes;
  hash_map<volumeType,vector<patchType>* > volVectors;
  virtual void computeVolVectors();

  // N.B. classes derived from DataProvider should compute volumes and
  // volVectors (possibly using computeVolVectors()), or provide new
  // definitions of volPatches(), getRandomPatch(vol), and getVolume()

public:
  virtual int totalPatchCount() const = 0;
  virtual void printPatch(ostream& out, const patchType patch, 
			  const int format = 0) const { out << patch; }
  virtual ~DataProvider() {
    for(hash_map<volumeType,vector<patchType>*>::iterator I=volVectors.begin();
	I != volVectors.end(); ++I)
      delete I->second;
  }
  const vector<patchType>* volPatches(volumeType vol) const {
    hash_map<volumeType,vector<patchType>* >::const_iterator I =
      volVectors.find(vol);
    if(I != volVectors.end())
      return I->second;
    else
      return 0;
  }
  virtual patchType getRandomPatch(const volumeType volume) const {
    hash_map<volumeType,vector<patchType>* >::const_iterator I =
      volVectors.find(volume);
    if(I == volVectors.end() || I->second->empty())
      errorExit("getRandomPatch(volume): no patches of required volume");
    const vector<patchType>* vp = I->second;
    return (*vp)[rnd(vp->size())];
  }
  virtual patchType getRandomPatch() const { return rnd(totalPatchCount()); }
  virtual volumeType getVolume(const patchType patch) const 
    { return volumes[patch]; }
  virtual void check(const Solution& sol) const {
    warning->print("no check performed by DataProvider");
  }
};

class CoveringDataProvider : public DataProvider {
public:
  virtual int totalPointCount() const = 0;
  virtual void printPoint(ostream& out, const pointType point, 
			  const int format = 0) const { out << point; }
  virtual void getCoverings
  (const patchType patch, const vector<pointType>*& cpoints,
   const vector<lambdaType>*& coverTimes, bool& sorted) const = 0;

  // If all points are covered the same number of times, coverTimes
  // contains only one element, the number of coverings. Otherwise,
  // the number each element is covered is in coverTimes with the same
  // indice as the point id in cpoints. If sorted == true, the
  // covering vector should contain a sentinel (a value larger than
  // any point number) following the last element of the vector (use
  // reserve and resize on the vector to achieve this). See also the
  // "merge" part of solution.C.

  // CoveringDataProvider allocates memory for the vectors, deallocation
  // is done when the cache size is too large in CoverMatrix

  virtual int maxCovering() const = 0;
  // returns the largest number of points that a single patch can
  // cover
  
  virtual long memNeeded(const patchType patch) const { return 0L; }
  // if easy, calculates the amount of memory (bytes) needed by the
  // covering data

  virtual const vector<lambdaType>* coveringGoal() const = 0;
  // there should be one element for each point

  virtual vector<patchType> coveringPatches(const pointType point) const {
    errorExit("coveringPatches() not implemented for this class");
    return vector<patchType>(0);
  }
  // some neighborhoods require that the patches that cover a given
  // point at least once are known (see also CoverMemoryMatrix)
};

#endif
