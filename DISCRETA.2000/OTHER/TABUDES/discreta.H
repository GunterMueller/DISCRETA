// discreta.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_discreta_h
#define PACKCOVER_discreta_h

#include "dp.H"

class DiscretaDataProvider : public CoveringDataProvider {
private:
  void read(istream&, long);
  int rows, cols, maxCover;
  vector<lambdaType>* coverGoal;
  vector<vector<pointType> > covers;
  vector<vector<lambdaType> > coverTimes;
public:
  DiscretaDataProvider(const string infile, const long maxmem);
  virtual int totalPatchCount() const { return cols; }
  virtual int totalPointCount() const { return rows; }
  virtual void getCoverings
  (const patchType patch, const vector<pointType>*& cpoints,
   const vector<lambdaType>*& ctimes, bool& sorted) const {
    cpoints = &covers[patch];
    ctimes = &coverTimes[patch];
    sorted = true;
  }
  virtual int maxCovering() const { return maxCover; }
  virtual long memNeeded(const patchType patch) const { return 0L; }
  virtual const vector<lambdaType>* coveringGoal() const { return coverGoal; }
};

#endif
