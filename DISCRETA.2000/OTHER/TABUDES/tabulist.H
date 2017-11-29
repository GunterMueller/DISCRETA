// tabulist.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_tabulist_h
#define PACKCOVER_tabulist_h

#include <vector>
#include <multiset.h>
#include <deque.h>
#include <list.h>

class TabuEntry {
private:
  const vector<int> attr;
  const int removalTime;
public:
  TabuEntry(const vector<int>& a, int rt) : attr(a), removalTime(rt) {}

  const vector<int>& getAttr() const { return attr; }
  int getRemovalTime() const { return removalTime; }

  bool CmpAttrib(const TabuEntry& lhs, const TabuEntry& rhs) const
    { return lhs.attr < rhs.attr; }
  bool CmpTime(const TabuEntry& lhs, const TabuEntry& rhs) const
    { return lhs.removalTime < rhs.removalTime; }
  // for some strange reason cannot access private members of
  // TabuEntry directly from AttribCmp member functions
  // (compiler bug?)
  struct AttribCmp {
    bool operator() (const TabuEntry& lhs, const TabuEntry& rhs) const
      { return lhs.CmpAttrib(lhs, rhs); }
  };
  struct TimeCmp {
    bool operator() (const TabuEntry& lhs, const TabuEntry& rhs) const
      { return lhs.CmpTime(lhs, rhs); }
  };
};

ostream& operator<<(ostream&, const TabuEntry&);

// for some strange reason, giving prototypes here and definitions in
// tabulist.C doesn't seem to work (compiler bug?)

class TabuList {
private:
  // by default, do not allow copy construction or assignment
  TabuList(const TabuList&) {}
  TabuList& operator=(const TabuList&) { return *this; }

  int lastPurge;
  void purgeTime(int);

  list<pair<Attribute *const, Tenure *const> > attributes;

  multiset<TabuEntry, TabuEntry::AttribCmp> attrSet;
  multiset<TabuEntry, TabuEntry::TimeCmp> timeSet;

  void purgeOld(int now) {
    while(lastPurge < now)
      purgeTime(++lastPurge);
  }

  void addTabu(const TabuEntry& te)
    { attrSet.insert(te); timeSet.insert(te); }
  void addTabu(const vector<int>& attr, int removalTime)
    { TabuEntry tmp(attr, removalTime); attrSet.insert(tmp); 
    timeSet.insert(tmp); }
  bool isTabu(const vector<int>& attr) const
    { TabuEntry tmp(attr, 0); return attrSet.find(tmp) != attrSet.end(); }
  bool isTabu(const TabuEntry& te) const
    { return attrSet.find(te) != attrSet.end(); }

public:

  TabuList() : lastPurge(-1) {}
  ~TabuList() {}

  void makeMove(Move& m, Solution& sol, volumeType newVolume, 
		penaltyType newPenalty, int now);
  void skipMove(int now);
  bool isTabuMove(Move& m, Solution& sol, volumeType newVolume, 
		  penaltyType newPenalty) const;

  void addAttributeToList(Attribute *const a, Tenure *const t) {
    attributes.push_back(pair<Attribute *const,Tenure *const>(a, t));
  }

  friend ostream& operator<<(ostream&, TabuList&);
};

#endif
