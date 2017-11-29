// move.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_move_h
#define PACKCOVER_move_h

#include <list.h>
#include <vector.h>
#include <algo.h>

#define no_cObj -1

class Solution;

class Move {
private:
  list<int> removed;
  list<patchType> added;
  int remCnt, addCnt;
public:
  Move() : remCnt(0), addCnt(0) {}
  Move(const Move& orig) {
    removed = orig.removed;
    added = orig.added;
    remCnt = orig.remCnt;
    addCnt = orig.addCnt;
  }
  Move& operator=(const Move& orig) {
    removed = orig.removed;
    added = orig.added;
    remCnt = orig.remCnt;
    addCnt = orig.addCnt;
    return *this;
  }
  ~Move() {}

  typedef typename list<patchType>::iterator patchIterator;
  typedef typename list<int>::iterator indexIterator;
  typedef typename list<patchType>::const_iterator c_patchIterator;
  typedef typename list<int>::const_iterator c_indexIterator;

  bool pureRemove() const { return added.size() == 0; }
  bool pureAdd() const { return removed.size() == 0; }
  bool changes() const { return added.size() > 0 && removed.size() > 0; }
  bool empty() const { return added.empty() && removed.empty(); }

  void addAddPatch(const patchType& patch) { 
    added.push_back(patch);
    ++addCnt;
  }
  void addRemoveIndex(int index) {
    ++remCnt;
    removed.push_back(index); 
  }

  const list<int>& getRemoveIndexes() const { return removed; }
  const list<patchType>& getAdded() const { return added; }

  void clear() { removed.clear(); added.clear(); addCnt = remCnt = 0; }
  int removeCount() const { return remCnt; }
  int addCount() const { return addCnt; }
  void append(const Move& move);
  void checkConsistency();
  void printRemovePatches(ostream& out, Solution& solution);
  bool operator<(const Move& rhs) const {
    if(remCnt < rhs.remCnt) return true;
    else if(remCnt > rhs.remCnt) return false;
    else if(addCnt < rhs.addCnt) return true;
    else if(addCnt > rhs.addCnt) return false;
    else if(removed < rhs.removed) return true;
    else if(rhs.removed < removed) return false;
    else if(added < rhs.added) return true;
    else if(rhs.added < added) return false;
    else return false;
  }
  void print(Solution& solution, ostream& out);
};

ostream& operator<<(ostream& out, Move& move);

#endif
