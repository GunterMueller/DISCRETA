// neighborhood.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#include "packcover.H"
#include "neighborhood.H"
#include <cstdlib>

int VChangeOneNeighborhood::neighborhoodSize() const
{
  int size = 0;
  for(int i = 0; i < solution.numPatches(); ++i) {
    patchType p = solution.getPatch(i);
    size += dataProvider.volPatches(dataProvider.getVolume(p))->size() - 1;
  }
  return size;
}

bool VChangeOneNeighborhood::findNextMove(Move& move)
{
  if(!firstIter)
    goto endloop;
  firstIter = false;
  for(iterIx = 0; 
      iterIx < (remVect.empty() ? solution.numPatches() : remVect.size()); 
      ++iterIx) {
    changeIx = (remVect.empty() ? iterIx : remVect[iterIx]);
    vp = dataProvider.volPatches
      (dataProvider.getVolume(solution.getPatch(changeIx)));
    for(chix = 0; chix < vp->size(); ++chix)
      if((*vp)[chix] != solution.getPatch(changeIx)) {
	if(!addSet.empty() && addSet.find((*vp)[chix]) == addSet.end())
	  continue;
	move.clear();
	move.addAddPatch((*vp)[chix]);
	move.addRemoveIndex(changeIx);
	return true;
      endloop:
	;
      }
  }
  return false;
}

bool VChangeOneNeighborhood::findRandomMove(Move& move)
{
  patchType p;
  do {
    if(remVect.empty())
      changeIx = rnd(solution.numPatches());
    else
      changeIx = remVect[rnd(remVect.size())];
    p = solution.getPatch(changeIx);
    changeTo = dataProvider.getRandomPatch(dataProvider.getVolume(p));
  } while(p == changeTo);
  move.clear();
  move.addAddPatch(changeTo);
  move.addRemoveIndex(changeIx);
  return true;
}

bool AddManyNeighborhood::findNextMove(Move& move)
{
  if(changeIxs[nextChange] < patchCount()-1)
    changeIxs[nextChange]++;
  else {
    while(changeIxs[++nextChange] >= patchCount() - 1 &&
	  nextChange < N) ;
    if(nextChange == N) return false;
    changeIxs[nextChange]++;
    for(int i = nextChange-1; i >= 0; --i)
      changeIxs[i] = changeIxs[nextChange];
    nextChange = 0;
  }
  move.clear();
  for(vector<int>::const_iterator I = changeIxs.begin();
      I != changeIxs.end(); ++I)
    if(addVect.empty())
      move.addAddPatch(patchType(*I));
    else
      move.addAddPatch(addVect[*I]);
  return true;
}

bool AddManyNeighborhood::findRandomMove(Move& move)
{
  move.clear();
  if(!patchCount()) return false;
  for(int i = 0; i < N; i++)
    if(addVect.empty())
      move.addAddPatch(patchType(rnd(dataProvider.totalPatchCount())));
    else
      move.addAddPatch(addVect[rnd(addVect.size())]);
  return true;
}

void AddManyNeighborhood::initNextRandom()
{
  if(N > 1)
    errorExit("AddManyNeighborhood::initNextRandom(): N > 1");
  nextRandom.resize(dataProvider.totalPatchCount());
  if(addVect.empty())
    for(int i = 0; i < dataProvider.totalPatchCount(); ++i)
      nextRandom[i] = patchType(i);
  else
    copy(addVect.begin(), addVect.end(), nextRandom.begin());
  nextRandomCount = dataProvider.totalPatchCount();
}

bool AddManyNeighborhood::findNextRandomMove(Move& move)
{
  if(nextRandomCount == 0) return false;
  move.clear();
  int ix = rnd(nextRandomCount);
  patchType add = nextRandom[ix];
  nextRandom[ix] = nextRandom[--nextRandomCount];
  move.addAddPatch(add);
  return true;
}

bool RemoveManyNeighborhood::findNextMove(Move& move)
{
  if(solution.empty()) return false;
  if(justInitialized) {
    justInitialized = false;
    for(int i = 0; i < N; ++i)
      changeIxs[i] = i;
  } else {
    if(changeIxs[0] >= ixCount() - N) return false;
    int j = 0;
    while(changeIxs[j+1] <= changeIxs[j] + 1 && j < N-1) ++j;
    ++changeIxs[j];
    for(int i = 0; i < j; ++i)
      changeIxs[i] = i;
  }
  move.clear();
  for(int i = 0; i < N; ++i)
    if(remVect.empty())
      move.addRemoveIndex(changeIxs[i]);
    else
      move.addRemoveIndex(remVect[changeIxs[i]]);
  return true;
}

bool RemoveManyNeighborhood::findRandomMove(Move& move)
{
  // Inefficient if N is large and very near to
  // solution.numPatches(). This maybe should be implemented by
  // choosing a random rank and then by unranking that subset
  if(solution.empty()) return false;
  move.clear();
  if(ixCount() < N) return false;
  set<int> soFar;
  for(int i = 0; i < N; i++) {
    int cand;
    do {
      if(remVect.empty())
	cand = rnd(solution.numPatches());
      else
	cand = remVect[rnd(remVect.size())];
    } while(soFar.find(cand) != soFar.end());
    soFar.insert(cand);
    move.addRemoveIndex(cand);
  }
  return true;
}

bool AddManyUniqueNeighborhood::findNextMove(Move& move)
{
  if(justInitialized) {
    justInitialized = false;
    for(int i = 0; i < N; ++i)
      patches[i] = i;
  } else {
    if(patches[0] >= patchCount() - N) return false;
    int j = 0;
    while(patches[j+1] <= patches[j] + 1 && j < N-1) ++j;
    ++patches[j];
    for(int i = 0; i < j; ++i)
      patches[i] = i;
  }
  move.clear();
  for(int i = 0; i < N; ++i)
    if(addVect.empty())
      move.addAddPatch(patches[i]);
    else
      move.addAddPatch(addVect[patches[i]]);
  return true;
}

bool AddManyUniqueNeighborhood::findRandomMove(Move& move)
{
  // Inefficient if N is large and very near to
  // solution.numPatches(). This maybe should be implemented by
  // choosing a random rank and then by unranking that subset
  move.clear();
  if(patchCount() < N) return false;
  set<patchType> soFar;
  for(int i = 0; i < N; i++) {
    patchType cand;
    do {
      if(addVect.empty())
	cand = patchType(rnd(dataProvider.totalPatchCount()));
      else
	cand = addVect[rnd(addVect.size())];
    } while(soFar.find(cand) != soFar.end());
    soFar.insert(cand);
    move.addAddPatch(cand);
  }
  return true;
}

void CompositeNeighborhood::initIter()
{
  for(vector<Neighborhood*>::iterator I = nbs.begin();
      I != nbs.end(); ++I)
    (*I)->initIter();
  firstTry = true;
}

bool CompositeNeighborhood::findNextMove(Move& move)
{
  if(firstTry) {
    if(nbs.empty()) return false;
    firstTry = false;
    for(unsigned i = 0; i < nbs.size(); ++i)
      if(!nbs[i]->findNextMove(moves[i])) return false;
  } else {
    unsigned tryIx = 0;
    while(tryIx < nbs.size() && !nbs[tryIx]->findNextMove(moves[tryIx]))
      ++tryIx;
    if(tryIx >= nbs.size()) return false;
    for(unsigned i = 0; i < tryIx; ++i) {
      nbs[i]->initIter();
      if(!nbs[i]->findNextMove(moves[i])) return false;
    }
  }
  move.clear();
  for(vector<Move>::const_iterator I = moves.begin(); I != moves.end(); ++I)
    move.append(*I);
  return true;
}

bool CompositeNeighborhood::findRandomMove(Move& move)
{
  if(nbs.empty()) return false;
  move.clear();
  Move m;
  for(vector<Neighborhood*>::const_iterator I = nbs.begin();
      I != nbs.end(); ++I) {
    if(!(*I)->findRandomMove(m)) return false;
    move.append(m);
  }
  move.checkConsistency();
  return true;
}

int CompositeNeighborhood::neighborhoodSize() const
{
  if(nbs.empty()) return 0;
  Checked<int> res(1);
  for(vector<Neighborhood*>::const_iterator I = nbs.begin();
      I != nbs.end(); ++I)
    res *= Checked<int>((*I)->neighborhoodSize());
  return res();
}

bool AdditiveNeighborhood::findRandomMove(Move& move)
  // Assume that the nbs neighborhood sizes stay constant! This is a
  // O(n) algorithm. O(lg(n)) is achieved by modifying a binary search
  // routine.
{
  if(nbs.empty()) return false;
  int nselect = rnd(nbSize);
  unsigned curr = 0;
  while(true) {
    nselect -= nbs[curr]->neighborhoodSize();
    if(nselect < 0) break;
    if(++curr >= nbs.size()) 
      errorExit("variable size child neighborhood in AdditiveNeighborhood");
  }
  if(!nbs[curr]->findRandomMove(move))
    errorExit("variable size child neighborhood in AdditiveNeighborhood");
  return true;
}

void AdditiveNeighborhood::initIter()
{
  nbSize = 0;
  for(unsigned i = 0; i < nbs.size(); ++i) {
    nbs[i]->initIter();
    nbSize += nbs[i]->neighborhoodSize();
  }
  currNb = 0;
}

bool AdditiveNeighborhood::findNextMove(Move& move)
{
  if(nbs.empty()) return false;
  while(currNb < nbs.size() && !(nbs[currNb]->findNextMove(move))) ++currNb;
  if(currNb >= nbs.size()) 
    return false;
  else
    return true;
}

int AdditiveNeighborhood::neighborhoodSize() const
{
  int nb = 0;
  for(unsigned i = 0; i < nbs.size(); ++i)
    nb += nbs[i]->neighborhoodSize();
  return nb; 
}

void LargeBalanceCoveringNeighborhood::initIter()
{
  set<patchType> alwIncs;
  set<int> alwDecs;
  for(int point = 0; point < solution.numPoints(); ++point)
    if(solution.pointPenalty(pointType(point)) != 0)
      if(solution.currentlyCovered(pointType(point)) <
	 solution.coveringGoal(pointType(point))) {
	vector<patchType> coverPatches = 
	  matrix.coveringPatches(pointType(point));
	for(unsigned i = 0; i < coverPatches.size(); ++i)
	  alwIncs.insert(coverPatches[i]);
      } else if(solution.currentlyCovered(pointType(point)) >
		solution.coveringGoal(pointType(point))) {
	vector<patchType> coverPatches = 
	  matrix.coveringPatches(pointType(point));
	sort(coverPatches.begin(), coverPatches.end());
	for(int i = 0; i < solution.numPatches(); ++i)
	  if(binary_search(coverPatches.begin(), coverPatches.end(), 
			   solution.getPatch(i)))
	    alwDecs.insert(i);
      } else 
	errorExit("LargeBalanceCoveringNeighborhood::initIter():  error");
  vector<patchType> limitIncs(alwIncs.size());
  vector<int> limitDecs(alwDecs.size());
  copy(alwIncs.begin(), alwIncs.end(), limitIncs.begin());
  copy(alwDecs.begin(), alwDecs.end(), limitDecs.begin());
  useIncCover = !limitIncs.empty();
  useDecCover = !limitDecs.empty();
  if(useIncCover) {
    incCoverNeigh->setNoLimits();
    incCoverNeigh->limitAdd(limitIncs);
  }
  if(useDecCover) {
    decCoverNeigh->setNoLimits();
    decCoverNeigh->limitRemove(limitDecs);
  }
  if(useIncCover && useDecCover)
    CompositeNeighborhood::initIter();
  else {
    if(useIncCover)
      incCoverNeigh->initIter();
    if(useDecCover)
      decCoverNeigh->initIter();
  }
}

bool LargeBalanceCoveringNeighborhood::findNextMove(Move& move)
{
  if(useIncCover && useDecCover)
    return CompositeNeighborhood::findNextMove(move);
  else if(useIncCover)
    return incCoverNeigh->findNextMove(move);
  else if(useDecCover)
    return decCoverNeigh->findNextMove(move);
  else
    errorExit("LargeBalanceCoveringNeighborhood: internal error");
  return false;
}

bool LargeBalanceCoveringNeighborhood::findRandomMove(Move& move)
{
  if(useIncCover && useDecCover)
    return CompositeNeighborhood::findRandomMove(move);
  else if(useIncCover)
    return incCoverNeigh->findRandomMove(move);
  else if(useDecCover)
    return decCoverNeigh->findRandomMove(move);
  else
    errorExit("LargeBalanceCoveringNeighborhood: internal error");
  return false;
}

int LargeBalanceCoveringNeighborhood::neighborhoodSize() const
{
  if(useIncCover && useDecCover)
    return CompositeNeighborhood::neighborhoodSize();
  else if(useIncCover)
    return incCoverNeigh->neighborhoodSize();
  else if(useDecCover)
    return decCoverNeigh->neighborhoodSize();
  else
    errorExit("LargeBalanceCoveringNeighborhood: internal error");
  return 0;
}

void BalanceCoveringNeighborhood::initIter()
{
  int counted = 0;
  do {
    if(++counted > solution.numPoints()) {
      addCheck = rnd(solution.numPoints());
      break;
    }
    if(++addCheck >= solution.numPoints())
      addCheck = 0;
  } while(solution.pointPenalty(pointType(addCheck)) == 0 ||
	  solution.currentlyCovered(addCheck) >
	  solution.coveringGoal(addCheck));
  //cout << "point " << addCheck
  //     << ": too few covers\n";
  vector<patchType> coverPatches = 
     matrix.coveringPatches(pointType(addCheck));
  //cout << "coverPatches: " << coverPatches << '\n';
  if(coverPatches.empty())
    errorExit("internal error 1: BalanceCoveringNeighborhood::initIter()");
  incCoverNeigh->setNoLimits();
  incCoverNeigh->limitAdd(coverPatches);
  counted = 0;
  do {
    if(++counted > solution.numPoints()) {
      patchType p = solution.getPatch(rnd(solution.numPatches()));
      const vector<pointType>* pts;
      const vector<lambdaType>* coverTimes;
      bool sorted;
      matrix.coveredBy(p, pts, coverTimes, sorted);
      decCheck = (*pts)[rnd(pts->size())];
      break;
    }
    if(++decCheck >= solution.numPoints())
      decCheck = 0;
  } while(solution.pointPenalty(pointType(decCheck)) == 0 ||
	  solution.currentlyCovered(decCheck) <
	  solution.coveringGoal(decCheck));
  //cout << "point " << decCheck
  //     << ": too many covers\n";
  decCoverNeigh->setNoLimits();
  vector<int> limitVect;
  limitVect.reserve(solution.numPatches());
  coverPatches = matrix.coveringPatches(pointType(decCheck));
  sort(coverPatches.begin(), coverPatches.end());
  for(int i = 0; i < solution.numPatches(); ++i)
    if(binary_search(coverPatches.begin(), coverPatches.end(), 
		     solution.getPatch(i)))
      limitVect.push_back(patchType(i));
  //cout << "limitVect: " << limitVect << '\n';
  if(limitVect.empty()) 
    errorExit("internal error 2: BalanceCoveringNeighborhood::initIter()");
  decCoverNeigh->limitRemove(limitVect);
  CompositeNeighborhood::initIter();
}

void LargeCorrectCoveringNeighborhood::initIter()
{
  set<patchType> alwIncs;
  set<int> alwDecs;
  for(int point = 0; point < solution.numPoints(); ++point)
    if(solution.pointPenalty(pointType(point)) != 0)
      if(solution.currentlyCovered(pointType(point)) <
	 solution.coveringGoal(pointType(point))) {
	vector<patchType> coverPatches = 
	  matrix.coveringPatches(pointType(point));
	for(unsigned i = 0; i < coverPatches.size(); ++i)
	  alwIncs.insert(coverPatches[i]);
      } else if(solution.currentlyCovered(pointType(point)) >
		solution.coveringGoal(pointType(point))) {
	vector<patchType> coverPatches = 
	  matrix.coveringPatches(pointType(point));
	sort(coverPatches.begin(), coverPatches.end());
	for(int i = 0; i < solution.numPatches(); ++i)
	  if(binary_search(coverPatches.begin(), coverPatches.end(), 
			   solution.getPatch(i)))
	    alwDecs.insert(i);
      } else 
	errorExit("LargeCorrectCoveringNeighborhood::initIter():  error");
  vector<patchType> limitIncs(alwIncs.size());
  vector<int> limitDecs(alwDecs.size());
  copy(alwIncs.begin(), alwIncs.end(), limitIncs.begin());
  copy(alwDecs.begin(), alwDecs.end(), limitDecs.begin());
  useIncCover = !limitIncs.empty();
  useDecCover = !limitDecs.empty();
  if(useIncCover) {
    incCoverNeigh->setNoLimits();
    incCoverNeigh->limitAdd(limitIncs);
  }
  if(useDecCover) {
    decCoverNeigh->setNoLimits();
    decCoverNeigh->limitRemove(limitDecs);
  }
  if(useIncCover && useDecCover)
    AdditiveNeighborhood::initIter();
  else {
    if(useIncCover)
      incCoverNeigh->initIter();
    else
      decCoverNeigh->initIter();
  }
}

void LargeCorrectCoveringNeighborhood::initNextRandom()
{
  initIter();
  if(useIncCover && useDecCover)
    AdditiveNeighborhood::initNextRandom();
  else if(useIncCover)
    incCoverNeigh->initNextRandom();
  else if(useDecCover)
    decCoverNeigh->initNextRandom();
  else
    errorExit("LargeCorrectCoveringNeighborhood: internal error");
}

bool LargeCorrectCoveringNeighborhood::findNextRandomMove(Move& move)
{
  if(useIncCover && useDecCover)
    return AdditiveNeighborhood::findNextRandomMove(move);
  else if(useIncCover)
    return incCoverNeigh->findNextRandomMove(move);
  else if(useDecCover)
    return decCoverNeigh->findNextRandomMove(move);
  else
    errorExit("LargeCorrectCoveringNeighborhood: internal error");
  return false;
}

bool LargeCorrectCoveringNeighborhood::findNextMove(Move& move)
{
  if(useIncCover && useDecCover)
    return AdditiveNeighborhood::findNextMove(move);
  else if(useIncCover)
    return incCoverNeigh->findNextMove(move);
  else if(useDecCover)
    return decCoverNeigh->findNextMove(move);
  else
    errorExit("LargeCorrectCoveringNeighborhood: internal error");
  return false;
}

bool LargeCorrectCoveringNeighborhood::findRandomMove(Move& move)
{
  if(useIncCover && useDecCover)
    return AdditiveNeighborhood::findRandomMove(move);
  else if(useIncCover)
    return incCoverNeigh->findRandomMove(move);
  else if(useDecCover)
    return decCoverNeigh->findRandomMove(move);
  else
    errorExit("LargeCorrectCoveringNeighborhood: internal error");
  return false;
}

int LargeCorrectCoveringNeighborhood::neighborhoodSize() const
{
  if(useIncCover && useDecCover)
    return AdditiveNeighborhood::neighborhoodSize();
  else if(useIncCover)
    return incCoverNeigh->neighborhoodSize();
  else if(useDecCover)
    return decCoverNeigh->neighborhoodSize();
  else
    errorExit("LargeCorrectCoveringNeighborhood: internal error");
  return 0;
}

void CorrectCoveringNeighborhood::nextCheckPoint()
  // caution: infloops if penalty is zero!
{
  do {
    if(++checkPoint >= solution.numPoints())
      checkPoint = 0;
  } while(solution.pointPenalty(pointType(checkPoint)) == 0);
}

void CorrectCoveringNeighborhood::initIter()
{
  nextCheckPoint();
  if(solution.currentlyCovered(checkPoint) < 
     solution.coveringGoal(checkPoint)) {
    //cout << "point " << checkPoint
    // << ": too few covers\n";
    currNeighborhood = incCoverNeigh;
    incCoverNeigh->setNoLimits();
    vector<patchType> coverPatches = 
      matrix.coveringPatches(pointType(checkPoint));
    //cout << "coverPatches: " << coverPatches << '\n';
    if(coverPatches.empty())
      errorExit("internal error: CorrectCoveringNeighborhood::initIter()");
    incCoverNeigh->limitAdd(coverPatches);
    incCoverNeigh->initIter();
  }
  else if(solution.currentlyCovered(checkPoint) > 
	  solution.coveringGoal(checkPoint)) {
    //cout << "point " << checkPoint
    // << ": too many covers\n";
    currNeighborhood = decCoverNeigh;
    decCoverNeigh->setNoLimits();
    vector<int> limitVect;
    limitVect.reserve(solution.numPatches());
    vector<patchType> coverPatches = 
      matrix.coveringPatches(pointType(checkPoint));
    sort(coverPatches.begin(), coverPatches.end());
    for(int i = 0; i < solution.numPatches(); ++i)
      if(binary_search(coverPatches.begin(), coverPatches.end(), 
		       solution.getPatch(i)))
	limitVect.push_back(i);
    //cout << "limitVect: " << limitVect << '\n';
    if(limitVect.empty()) 
      errorExit("internal error: CorrectCoveringNeighborhood::initIter()");
    decCoverNeigh->limitRemove(limitVect);
    decCoverNeigh->initIter();
  }
  else
    errorExit("internal error: CorrectCoveringNeighborhood::initIter()");
}

bool CorrectCoveringNeighborhood::findNextMove(Move& move)
{
  return currNeighborhood->findNextMove(move);
}

bool CorrectCoveringNeighborhood::findRandomMove(Move& move)
{
  return currNeighborhood->findRandomMove(move);
}

int CorrectCoveringNeighborhood::neighborhoodSize() const
{
  return currNeighborhood->neighborhoodSize();
}

void CorrectRandomCoveringNeighborhood::nextCheckPoint()
  // caution: infloops if penalty is zero or the PRNG is bad!
{
  do {
    checkPoint = rnd(solution.numPoints());
    // this could be done more efficiently
  } while(solution.pointPenalty(pointType(checkPoint)) == 0);
}

void MultiNeighborhood::initIter()
{
  int minDiff = numeric_limits<int>::max(), diff;
  currNPtr = 0;
  for(vector<Neighborhood *const>::const_iterator I = neighs.begin();
      I != neighs.end(); ++I) {
    diff = abs((*I)->neighborhoodSize() - goalSize);
    if(diff < minDiff) {
      minDiff = diff;
      currNPtr = *I;
    }
  }
  if(!currNPtr) errorExit("MultiNeighborhood: could not select neighborhood");
  currNPtr->initIter();
}
