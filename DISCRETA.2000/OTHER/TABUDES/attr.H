// attr.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_attr_h
#define PACKCOVER_attr_h

#include <vector.h>
#include <list.h>
#include <deque.h>
#include "utils.H"
#include "limits.H"

class AttrIdProvider {
private:
  static int attrIds;
public:
  static int getId() { return attrIds++; }
};

class Attribute {
private:
  int attrId;
protected:
  // binary copy ok
public:
  Attribute() { attrId = AttrIdProvider::getId(); }
  virtual ~Attribute() {}
  int getAttrId() const { return attrId; }

  // makeAttr(m1,...) and makeTestAttr(m2,...) should return the same
  // attribute when m2 is the inverse of m1 (this is more complicated
  // when sol, newVolume, and/or newPenalty are used in the
  // construction of the attribute)

  // N.B. both makeAttr and makeTestAttr should be called BEFORE the
  // move is actually made. Both functions should include the value of
  // attrId as the last entry of each vector (so that the different
  // attributes won't mix).

  virtual list<vector<int> > makeAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const
    { 
      list <vector<int> > tmp;
      tmp.push_back(makeVector(getAttrId())); 
      return tmp; 
    }
  virtual list<vector<int> > makeTestAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const
    { return makeAttr(m, sol, newVolume, newPenalty); }
};

class RemoveAttribute : public Attribute {
  // if a move adds a patch in the solution, then a move removing the
  // same patch is forbidden
public:
  virtual list<vector<int> > makeAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
  virtual list<vector<int> > makeTestAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
};

class AddAttribute : public Attribute {
  // if a move removes a patch from the solution, then the same patch
  // cannot be added in the solution
public:
  virtual list<vector<int> > makeAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
  virtual list<vector<int> > makeTestAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
};

class AddAllAttribute : public Attribute {
  // if a move removes a set of patches from the solution, then
  // excactly the same set of patches cannot be added in the solution
public:
  virtual list<vector<int> > makeAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
  virtual list<vector<int> > makeTestAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
};

class IndexAttribute : public Attribute {
  // if the patches removed from a solution are indexed i1, i2, ...,
  // then a move removing any of these _indexes_ is forbidden
public:
  virtual list<vector<int> > makeAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
  virtual list<vector<int> > makeTestAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const
    {
      return makeAttr(m, sol, newVolume, newPenalty);
    }
};

class ChangeAttribute : public Attribute {
  // if a move removes patch a and adds patch b, then another move
  // cannot remove b and add a
public:
  virtual list<vector<int> > makeAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
  virtual list<vector<int> > makeTestAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
};

class PenaltyVolumeLoopAttribute : public Attribute {
  // move is forbidden, if it would make a sequence of last (penalty,
  // volume) pairs of length checkLen identical to a sequence on
  // current tabu list.

  // N.B. it is assumed that makeAttr is called
  // _exactly once_ for each committed move.
private:
  int checkLen;
  mutable deque<int> lastValues;
public:
  PenaltyVolumeLoopAttribute(int cl) : Attribute(), checkLen(cl) {
      if(cl <= 0) 
      errorExit("PenaltyVolumeLoopAttribute: positive loop length required");
    for(int i = 0; i < 2*checkLen; ++i) 
      lastValues.push_back(numeric_limits<int>::max());
  }
  virtual list<vector<int> > makeAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
  virtual list<vector<int> > makeTestAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
};

class CompositeAttribute : public Attribute {
private:
  vector<Attribute*> attrs;
public:
  void addAttribute(Attribute *a) { attrs.push_back(a); }
  virtual list<vector<int> > makeAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
  virtual list<vector<int> > makeTestAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
};

class NeighborhoodAttribute : public Attribute {
private:
  const Neighborhood& neighborhood;
public:
  NeighborhoodAttribute(const Neighborhood& nb) : neighborhood(nb) {}
  virtual list<vector<int> > makeAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
  virtual list<vector<int> > makeTestAttr
  (Move& m, Solution& sol, volumeType newVolume, penaltyType newPenalty) const;
};

#endif
