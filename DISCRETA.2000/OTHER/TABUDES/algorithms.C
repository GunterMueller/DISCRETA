// algorithms.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#include <iostream.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "packcover.H"
#include "algorithms.H"
#include "args.H"
#include "utils.H"
#include "move.H"
#include "solution.H"
#include "neighborhood.H"
#include "tabulist.H"

void randomWalk(Solution& solution, vector<Neighborhood*>& nvect,
		int maxiter)
{
  if(nvect.size() != 1) errorExit("single neighborhood expected");
  Neighborhood& neighborhood = *nvect[0];
  int iterCount = 0;
  volumeType volCh;
  penaltyType penCh;
  Move move;
  if(verbose)
    cout << "initial penalty = " << solution.getCurrentPenalty()
	 << ", initial volume = " << solution.getCurrentVolume()
	 << '\n';
  while(iterCount++ < maxiter) {
    neighborhood.initIter();
    if(!neighborhood.findRandomMove(move))
      errorExit("no moves!");
    solution.whatIfChange(move, volCh, penCh);
    solution.makeMoveVP(move, volCh, penCh);
    if(verbose)
      cout << "penalty = " << solution.getCurrentPenalty()
	   << ", volume = " << solution.getCurrentVolume() << '\n';
  }
}

void localOptimize(Solution& solution, vector<Neighborhood*>& nvect)
{
  if(nvect.size() != 1) errorExit("single neighborhood expected");
  Neighborhood& neighborhood = *nvect[0];
  Move move;
  volumeType volCh;
  penaltyType penCh;
  if(verbose) 
    cout << "initial penalty: " << solution.getCurrentPenalty() << '\n';
  while(true) {
    neighborhood.initIter();
    if(verbose) 
      cout << "neighbors: " << neighborhood.neighborhoodSize() << '\n';
    bool found = false;
    int tryCount = 0;
    for(int i = neighborhood.neighborhoodSize(); i > 0; --i) {
      neighborhood.findRandomMove(move);
      ++tryCount;
      solution.whatIfChange(move, volCh, penCh);
      if(penCh < penaltyType(0)) {
	found = true;
	break;
      }
    }
    if(found) {
      solution.makeMoveVP(move, volCh, penCh);
      if(verbose)
	cout << "penalty: " << solution.getCurrentPenalty()
	     << " (" << tryCount << " tries)\n";
    } else {
      if(verbose) cout << "searching complete neighborhood\n";
      while(neighborhood.findNextMove(move)) {
	solution.whatIfChange(move, volCh, penCh);
	++tryCount;
	if(penCh < penaltyType(0)) {
	  found = true;
	  break;
	}
      }
      if(found) {
	solution.makeMoveVP(move, volCh, penCh);
	cout << "penalty: " << solution.getCurrentPenalty()
	     << " (" << tryCount << " tries)\n";
      } else {
	if(verbose)
	  cout << "no improving moves\n";
	return;
      }
    }
  }
}

void twoWaySearch
(Solution& solution, vector<Neighborhood*>& nvect, penaltyType endLimit)
{
  if(nvect.size() != 2) 
    errorExit("twoWaySearch requires two neighborhoods");
  Neighborhood& incNeigh = *(nvect[0]);
  Neighborhood& decNeigh = *(nvect[1]);
  int iterCount = 0, maxIter = options->get("maxiter", 1000);
  patchType lastRemoval = numeric_limits<patchType>::max();
  volumeType volCh;
  penaltyType penCh;
  penaltyType bestPenalty = numeric_limits<penaltyType>::max();
  if(verbose)
    cout << "initial penalty = " << solution.getCurrentPenalty()
	 << ", volume = " << solution.getCurrentVolume() << '\n';
  while(solution.getCurrentPenalty() > endLimit && iterCount < maxIter) {
    incNeigh.initNextRandom();
    Move move;
    // try to find an improving move in incNeigh
    while(incNeigh.findNextRandomMove(move)) {
      solution.whatIfChange(move, volCh, penCh);
      if(/*move.getAdded().front() != lastRemoval && */penCh <= 0) {
	++iterCount;
	solution.makeMoveVP(move, volCh, penCh);
	if(verbose)
	  cout << "add patch " << move.getAdded().front() << '\n';
	goto nextIter;
      }
    }
    // no improving move found, make one move in decNeigh
    decNeigh.initIter();
    if(!decNeigh.findRandomMove(move))
      errorExit("twoWaySearch: no moves in decNeigh");
    lastRemoval = solution.getPatch(move.getRemoveIndexes().front());
    if(verbose)
      cout << "remove patch " << lastRemoval << '\n';
    solution.makeMove(move);
    ++iterCount;
  nextIter:
    if(verbose)
      cout << "penalty " << solution.getCurrentPenalty()
	   << ", volume " << solution.getCurrentVolume() << '\n';
    if(solution.getCurrentPenalty() < bestPenalty)
      bestPenalty = solution.getCurrentPenalty();
  }
  cout << "best penalty: " << bestPenalty << '\n'
       << iterCount << " iterations\n";
}

void basicTabuSearch(Solution& solution, vector<Neighborhood*>& nvect,
		     TabuList& tabuList, penaltyType endLimit)
  // if multiple neighborhoods are given, the first ones in nvect have
  // preference over the last ones in case of equal costs
{
  int iterCount = 0, maxIter = options->get("maxiter", 1000);
  Solution bestSol = solution;
  Move move, selectMove;
  vector<Move> bestMoves;
  vector<volumeType> bestVols;
  penaltyType penCh, propPen, bestPen, bestVal;
  volumeType volCh, propVol, selectVol;
  int tabuCount, neighCount, imprCount, zeroCount;
  if(verbose) {
    cout << "starting tabu search...\n";
    cout << "curr neighs volu tabu impr zero best alwd bcnt\n"
	 << "----------------------------------------------\n";
  }
  // cout << "initial solution:\n" << solution << '\n';
  penaltyType smallPenalty = options->get
    ("smallpenalty", numeric_limits<penaltyType>::max());
  int largeRemaining = 0;
  const int optLargeRemaining = options->get("largeremaining", int());;
  while(solution.getCurrentPenalty() > endLimit && iterCount < maxIter) {
    bestPen = bestVal = numeric_limits<penaltyType>::max();
    neighCount = tabuCount = imprCount = zeroCount = 0;
    bestMoves.clear();
    bestVols.clear();

    for(unsigned nvectIx = 0; nvectIx < nvect.size(); ++nvectIx) {
      if(nvectIx > 0 && 
	 (solution.getCurrentPenalty() > smallPenalty ||
	  solution.getCurrentPenalty() != 
	  min(bestSol.getCurrentPenalty(),bestPen) ||
	  !largeRemaining))
	continue;
      //cout << largeRemaining << '\n';
      Neighborhood& neighborhood = *nvect[nvectIx];
      neighborhood.initIter();
      //cout << "neighborhood size: " << neighborhood.neighborhoodSize()
      //	   << '\n';
      while(neighborhood.findNextMove(move)) {
	//cout << "move: ";
	//move.print(solution, cout);
	//cout << '\n';
	if(move.empty())
	  warning->print("empty move encountered");
	++neighCount;
	solution.whatIfChange(move, volCh, penCh);
	propVol = solution.getCurrentVolume() + volCh;
	propPen = solution.getCurrentPenalty() + penCh;
	if(propPen <= bestPen && 
	   !tabuList.isTabuMove(move, solution, propVol, propPen)) {
	  if(propPen < bestPen) {
	    bestMoves.clear();
	    bestVols.clear();
	  }
	  bestPen = propPen;
	  bestMoves.push_back(move);
	  bestVols.push_back(propVol);
	}
	if(verbose > 1) {
	  if(tabuList.isTabuMove(move, solution, propVol, propPen)) {
	    ++tabuCount;
	    //cout << "tabu move: ";
	    //move.print(solution, cout);
	    //cout << '\n';
	  }
	  if(propPen < solution.getCurrentPenalty())
	    ++imprCount;
	  if(propPen == solution.getCurrentPenalty())
	    ++zeroCount;
	  if(propPen < bestVal)
	    bestVal = propPen;
	}
      }
    }
    if(bestPen == numeric_limits<penaltyType>::max()) {
      warning->print("all moves tabu");
      tabuList.skipMove(iterCount++);
    } else {
      if(verbose > 1) {
	printf("%4d %6d %4d %4d %4d %4d %4d %4d %4u ",
	       solution.getCurrentPenalty(),
	       neighCount, solution.getCurrentVolume(), tabuCount,
	       imprCount, zeroCount, bestVal, bestPen, bestMoves.size());
	// solution.printStatistics(cout);
	cout << '\n';
      } else if(verbose) {
	cout << solution.getCurrentPenalty() << " ";
	cout.flush();
      }
      int select = rnd(bestMoves.size());
      //cout << bestMoves.size() << " moves: selecting number " << select 
      //     << '\n';
      selectMove = bestMoves[select];
      selectVol = bestVols[select];
      //cout << "made move: ";
      //selectMove.print(solution, cout);
      //cout << '\n';
      tabuList.makeMove(selectMove, solution, selectVol, bestPen, iterCount++);
      if(solution.getCurrentPenalty() > bestPen)
	largeRemaining = optLargeRemaining;
      else if(largeRemaining) --largeRemaining;
      solution.makeMoveVP(selectMove, selectVol - solution.getCurrentVolume(),
			  bestPen - solution.getCurrentPenalty());
      if(solution.getCurrentPenalty() < bestSol.getCurrentPenalty())
	bestSol = solution;
    }
    //cout << "Tabu list:\n" << tabuList << '\n';
  }
  if(verbose)
    cout << "\ntabu search completed (" 
	 << (solution.getCurrentPenalty() <= endLimit ? "endLimit" : 
	     "iterations")
	 << ")\n";
  solution = bestSol;
  cout << "penalty " << solution.getCurrentPenalty() << ", "
       << iterCount << " iterations\n";
}

void iterateWithinPenaltyLimit
(Solution& solution, vector<Neighborhood*>& nvect, penaltyType minPen,
 penaltyType maxPen, int maxiter)
  // Finds all moves that keep the penalty within the penalty limits
  // and then selects one at random. Repeats this until no such moves
  // can be found.
{
  if(nvect.size() != 1) errorExit("single neighborhood expected");
  Neighborhood& neighborhood = *nvect[0];
  vector<Move> moves;
  Move move;
  penaltyType penCh;
  volumeType volCh;
  int iterCount = 0;
  while(iterCount++ < maxiter) {
    neighborhood.initIter();
    if(verbose)
      cout << "checking " << neighborhood.neighborhoodSize() << " moves...\n";
    penaltyType currPen = solution.getCurrentPenalty();
    while(neighborhood.findNextMove(move)) {
      solution.whatIfChange(move, volCh, penCh);
      if(currPen + penCh <= maxPen && currPen + penCh >= minPen)
	moves.push_back(move);
    }
    if(moves.empty())
      return;
    if(verbose)
      cout << "found " << moves.size() << " feasible moves, taking one\n";
    solution.makeMove(moves[rnd(moves.size())]);
    moves.clear();
  }
}

void randomWalkWithinPenaltyLimit
(Solution& solution, vector<Neighborhood*>& nvect, penaltyType minPen,
 penaltyType maxPen, int maxiter)
  // Walks randomly within the neighborhood and the penalty limits
  // until maxiter moves have been made or no neighbors with
  // acceptable penalty exist.
{
  if(nvect.size() != 1) errorExit("single neighborhood expected");
  Neighborhood& neighborhood = *nvect[0];
  set<Move> checkedMoves;
  Move move;
  penaltyType penCh;
  volumeType volCh;
  int iterCount = 0;
  while(iterCount++ < maxiter) {
    bool found = false;
    neighborhood.initIter();
    int L = neighborhood.neighborhoodSize(), tryCount = 0;
    penaltyType currPen = solution.getCurrentPenalty();
    for(int i = 0; i < L; ++i) {
      ++tryCount;
      neighborhood.findRandomMove(move);
      if(checkedMoves.find(move) == checkedMoves.end()) {
	solution.whatIfChange(move, volCh, penCh);
	if(currPen + penCh <= maxPen && currPen + penCh >= minPen) {
	  found = true;
	  solution.makeMoveVP(move, volCh, penCh);
	  break;
	}
	checkedMoves.insert(move);
      }
    }
    if(verbose)
      if(found) cout << tryCount << " tries, penalty " 
		     << solution.getCurrentPenalty() << ", volume "
		     << solution.getCurrentVolume() << '\n';
      else cout << tryCount << " tries, not found, making complete check\n";
    if(!found) {
      neighborhood.initIter();
      while(neighborhood.findNextMove(move))
	if(checkedMoves.find(move) == checkedMoves.end()) {
	  solution.whatIfChange(move, volCh, penCh);
	  if(currPen + penCh <= maxPen && currPen + penCh >= minPen) {
	    found = true;
	    solution.makeMoveVP(move, volCh, penCh);
	    break;
	  }
	}
      if(verbose)
	if(found)
	  cout << "move found after " << checkedMoves.size() 
	       << " penalty change calculations, penalty "
	       << solution.getCurrentPenalty() << ", volume "
	       << solution.getCurrentVolume() << '\n';
	else
	  cout << "no acceptable moves left, terminating search\n";
    }
    checkedMoves.clear();
    if(!found) break;
  }
}

#define T_LIFESAVER 1.0
#define T_ITER 300
#define T_PRINT_ITER 300

void basicSimulatedAnnealing(Solution& solution, vector<Neighborhood*>& nvect,
			     penaltyType endLimit)
  // constant-volume simulated annealing algorithm, tries to get the
  // penalty to or below acceptPenalty
{
  if(nvect.size() != 1) errorExit("single neighborhood expected");
  Neighborhood& neighborhood = *nvect[0];
  neighborhood.initIter();
  double coolFact = options->get("CF", 0.99);
  int iterLength;
  if(!options->getNext("L", iterLength))
    iterLength = int(options->getMust("LF", double()) * 
		     neighborhood.neighborhoodSize() + .5);
  int frozen = options->get("frozen", 10);
  float deltaF_;
  float r, D;
  penaltyType costDelta, currCost, lastCost, bestSeen;
  int notChanged = 0, i, m1, m2, m3, m0;
  float T;
  Solution bestSol = solution;
  volumeType volCh;
  Move move;
  int iterCounter = 0;

  if(verbose)
    cout << "Starting annealing...\n\n";
  currCost = solution.getCurrentPenalty();
  bestSeen = lastCost = currCost;
  
  if(!options->getNext("IT", T)) {
    int i, m2;
    T = 0.0;
    m2 = 0;
    for(i = 0; i < T_ITER; i++) {
      neighborhood.findRandomMove(move);
      solution.whatIfChange(move,  volCh, costDelta);
      if(costDelta > 0) {
	m2++;
	T += -costDelta;
      }
    }
    if(m2 == 0) {
      T = T_LIFESAVER;
      warning->print("approximate temperature not found");
    }
    else
      T = T / m2 / log(options->getMust("IP", double()));
  }

  /* tests the probability for cost increasing moves */
  if(verbose) {
    m1 = m2 = 0;
    deltaF_ = 0.0;
    for(i = 0; i < T_PRINT_ITER; i++) {
      neighborhood.findRandomMove(move);
      solution.whatIfChange(move,  volCh, costDelta);
      if(costDelta > 0) {
	m1++;
	if(random01() < exp(-costDelta / T))
	  m2++;
      }
    }
    if(m1 == 0) {
      m1 = 1;
      m2 = 0;
    }
    printf("initial inc%%  = %.2f\n\n", (double) m2 / (double) m1);
  }

  if(verbose)
    printf("      T      cost   best   inc%%   tot%%   0-m%%\n"
	   "    ------------------------------------------\n");
  while(notChanged < frozen) {
    m1 = m2 = m3 = m0 = 0;
    for(i = 0; i < iterLength; i++) {
      neighborhood.findRandomMove(move);
      solution.whatIfChange(move,  volCh, costDelta);
      iterCounter++;
      if(costDelta <= 0) {
	m3++;
	solution.makeMoveVP(move, volCh, costDelta);
	neighborhood.initIter();
	currCost += costDelta;
	if(currCost <= endLimit) {
	  if(verbose)
	    printf("\n...annealing accomplished.\n\n");
	  return;       /* a good enough final solution was found */
	}
	if(costDelta < 0) {
	  notChanged = 0;
	  if(currCost < bestSeen) {
	    bestSeen = currCost;
	    bestSol = solution;
	  }
	}
	else
	  m0++;
      }
      else {
	r = random01();
	D = costDelta / T;
	if(r < exp(-D)) {
	  solution.makeMoveVP(move, volCh, costDelta);
	  neighborhood.initIter();
	  m1++;
	  currCost += costDelta;
	}
	else
	  m2++;
      }
    }
    if(lastCost <= currCost)
      notChanged++;
    lastCost = currCost;
    if(m2 == 0)
      m2 = 1; /* prevent division by zero */
    if(verbose)
      printf("    %5.2f   %4d   %4d    %4.1f   %4.1f   %4.3f\n", 
	     T, currCost, bestSeen, (double) m1 / (double) (m1 + m2) * 100.0,
	     (double) (m3 + m1) / (double) (m1 + m2 + m3) * 100.0,
	     (double) (m0) / (double) (m1 + m2 + m3));
    T *= coolFact;
  }
  if(verbose)
    printf("\n...annealing accomplished.\n\n");
  solution = bestSol;
}

void greatDeluge(Solution& solution, vector<Neighborhood*>& nvect, 
		 penaltyType endLimit)
  // slightly modified version from Dueck's original algorithm: does not
  // get immediately stuck on a flat surface
{
  if(nvect.size() != 1) errorExit("single neighborhood expected");
  Neighborhood& neighborhood = *nvect[0];
  neighborhood.initIter();
  int frozen = options->get("frozen", 1000);
  double down = options->get("down", 0.001);
  int notChanged = 0;
  penaltyType curr = solution.getCurrentPenalty(), penCh;
  volumeType volCh;
  Move move;
  double waterLevel = options->get("initlevel", curr);

  while(notChanged < frozen && curr > endLimit) {
    if(verbose)
      cout << waterLevel << ": " << curr << '\n';
    if(!neighborhood.findRandomMove(move)) errorExit("no moves!");
    solution.whatIfChange(move, volCh, penCh);
    if(curr + penCh <= waterLevel) {
      solution.makeMoveVP(move, volCh, penCh);
      neighborhood.initIter();
      curr += penCh;
    }
    if(curr > waterLevel - down)
      ++notChanged;
    else {
      waterLevel = waterLevel - down;
      notChanged = 0;
    }
  }
  if(solution.getCurrentPenalty() != curr)
    errorExit("penalty calculation error");
}

void thresholdAccepting(Solution& solution, vector<Neighborhood*>& nvect,
			penaltyType endLimit)
{
  if(nvect.size() != 1) errorExit("single neighborhood expected");
  Neighborhood& neighborhood = *nvect[0];
  neighborhood.initIter();
  double threshold = options->getMust("threshold", double());
  double down = options->get("down", 0.01);
  int iterLength;
  if(!options->getNext("L", iterLength))
    iterLength = int(options->getMust("LF", double()) * 
		     neighborhood.neighborhoodSize() + .5);
  int frozen = options->get("frozen", 10);
  penaltyType costDelta, currCost, lastCost, bestSeen;
  int notChanged = 0, i, m1, m2, m3, m0;
  Solution bestSol = solution;
  volumeType volCh;
  Move move;
  int iterCounter = 0;

  if(verbose)
    cout << "Starting threshold accepting...\n\n";

  currCost = solution.getCurrentPenalty();
  bestSeen = lastCost = currCost;

  if(verbose)
    printf("      T      cost   best   inc%%   tot%%   0-m%%\n"
	   "    ------------------------------------------\n");
  while(notChanged < frozen) {
    m1 = m2 = m3 = m0 = 0;
    for(i = 0; i < iterLength; i++) {
      neighborhood.findRandomMove(move);
      solution.whatIfChange(move,  volCh, costDelta);
      iterCounter++;
      if(costDelta <= threshold) {
	if(costDelta <= 0)
	  m3++;
	else if(costDelta > 0)
	  m1++;
	else
	  m0++;
	solution.makeMoveVP(move, volCh, costDelta);
	neighborhood.initIter();
	currCost += costDelta;
	if(currCost <= endLimit) {
	  if(verbose)
	    printf("\n...threshold accepting accomplished.\n\n");
	  return;       /* a good enough final solution was found */
	}
	if(costDelta < 0) {
	  notChanged = 0;
	  if(currCost < bestSeen) {
	    bestSeen = currCost;
	    bestSol = solution;
	  }
	}
      } else
	++m2;
    }
    if(lastCost <= currCost)
      notChanged++;
    lastCost = currCost;
    if(m2 == 0)
      m2 = 1; /* prevent division by zero */
    if(verbose)
      printf("    %5.2f   %4d   %4d    %4.1f   %4.1f   %4.3f\n", 
	     threshold, currCost, bestSeen, 
	     (double) m1 / (double) (m1 + m2) * 100.0,
	     (double) (m3 + m1) / (double) (m1 + m2 + m3) * 100.0,
	     (double) (m0) / (double) (m1 + m2 + m3));
    threshold = max(0., threshold - down);
  }
  if(verbose)
    printf("\n...threshold accepting accomplished.\n\n");
  solution = bestSol;
}

void recordToRecordTravel(Solution& solution, vector<Neighborhood*>& nvect,
			  penaltyType endLimit)
{
  if(nvect.size() != 1) errorExit("single neighborhood expected");
  Neighborhood& neighborhood = *nvect[0];
  neighborhood.initIter();
  penaltyType curr, bestPenaltySeen = solution.getCurrentPenalty();
  curr = bestPenaltySeen;
  int notChanged = 0;
  int frozen = options->get("frozen", 1000);
  penaltyType deviation = options->get("deviation", penaltyType(0)), penCh;
  volumeType volCh;
  Move move;
  int specialRecord = options->get("specialrecord", 0);
  int checkCount = 0;

  while(notChanged++ < frozen && curr > endLimit) {
    ++checkCount;
    if(verbose)
      cout << bestPenaltySeen << ": " << curr 
	   << " (vol " << solution.getCurrentVolume() << ")\n";
    if(!neighborhood.findRandomMove(move)) errorExit("no moves!");
    solution.whatIfChange(move, volCh, penCh);
    if(curr + penCh <= bestPenaltySeen + deviation) {
      solution.makeMoveVP(move, volCh, penCh);
      if(!specialRecord)
	neighborhood.initIter();
      curr += penCh;
    }
    if(curr < bestPenaltySeen) {
      bestPenaltySeen = curr;
      notChanged = 0;
    }
    if(specialRecord && checkCount >= specialRecord) {
      checkCount = 0;
      neighborhood.initIter();
    }
  }
  if(solution.getCurrentPenalty() != curr)
    errorExit("penalty calculation error");
}
