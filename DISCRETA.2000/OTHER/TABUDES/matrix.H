// matrix.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_matrix_h
#define PACKCOVER_matrix_h

#include <hash_map.h>
#include "limits.H"
#include "checked.H"
#include "dp.H"

class PatchMatrix {
private:
  const long maxMemory;
public:
  PatchMatrix(long mxmem) : maxMemory(mxmem) {}
  getMaxMemory() const { return maxMemory; }
  virtual void print(ostream& out) const { 
    errorExit("matrix print not implemented yet"); 
  }
};

inline ostream& operator<<(ostream& out, PatchMatrix& matrix) {
  matrix.print(out);
  return out;
}

class CoverMatrix : public PatchMatrix {
protected:
  const CoveringDataProvider& dataProvider;
  const patchType noPatch;
public:
  CoverMatrix(CoveringDataProvider& dp, long mxmem) : 
    PatchMatrix(mxmem), dataProvider(dp), 
    noPatch((Checked<patchType>(dp.totalPatchCount()))()) {};

  // N.B. the vectors returned by coveredBy() are only guaranteed to
  // exist until the next invocation of coveredBy()!!! (unless
  // lockCache() and releaseCache() are used)

  virtual void coveredBy
  (const patchType patch, const vector<pointType>*& points, 
   const vector<lambdaType>*& coverTimes, bool& sorted) const = 0;

  virtual void lockCache() const {}
  // coverings calculated after lockCache() cannot be erased until
  // releaseCache() is called
  virtual void releaseCache() const {}

  virtual void printDebug() {}
  virtual void printStatistics(ostream& out) const {}

  virtual const vector<patchType> coveringPatches
  (const pointType point) const {
    return dataProvider.coveringPatches(point);
  }
};

class CoverMemoryMatrix : public CoverMatrix {
private:
  typedef vector<pointType> pointVectorType;
  typedef vector<patchType> patchVectorType;
  typedef vector<lambdaType> lambdaVectorType;
  vector<pointVectorType*> pointVectors;
  vector<lambdaVectorType*> lambdaVectors;
  vector<bool> sortedVector;
  mutable vector<patchVectorType*> inverseVectors;
  mutable long memSum;
public:
  CoverMemoryMatrix(CoveringDataProvider& dp, long mxmem);
  virtual ~CoverMemoryMatrix() {
    for(unsigned i = 0; i < pointVectors.size(); ++i) delete pointVectors[i];
    for(unsigned i = 0; i < lambdaVectors.size(); ++i) delete lambdaVectors[i];
  }

  virtual void coveredBy
  (const patchType patch, const vector<pointType>*& points, 
   const vector<lambdaType>*& coverTimes, bool& sorted) const {
    points = pointVectors[patch];
    coverTimes = lambdaVectors[patch];
    sorted = sortedVector[patch];
  }
  virtual void calculateInverses() const;
  virtual const vector<patchType> coveringPatches
  (const pointType point) const {
    if(inverseVectors.empty())
      calculateInverses();
    return *inverseVectors[point];
  }
  virtual void print(ostream& out) const;
};

#endif
