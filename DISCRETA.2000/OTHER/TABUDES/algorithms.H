// algorithms.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#ifndef PACKCOVER_algorithms_h
#define PACKCOVER_algorithms_h

#include <vector.h>

void basicSimulatedAnnealing(Solution&, vector<Neighborhood*>&, penaltyType);

void thresholdAccepting(Solution&, vector<Neighborhood*>&, penaltyType);

void recordToRecordTravel(Solution&, vector<Neighborhood*>&, penaltyType);

void greatDeluge(Solution&, vector<Neighborhood*>&, penaltyType);

void iterateWithinPenaltyLimit
(Solution& solution, vector<Neighborhood*>& nvect, penaltyType minPen,
 penaltyType maxPen, int maxiter);

void randomWalk(Solution& solution, vector<Neighborhood*>& nvect,
		int maxiter);

void randomWalkWithinPenaltyLimit
(Solution& solution, vector<Neighborhood*>& nvect, penaltyType minPen,
 penaltyType maxPen, int maxiter);

void localOptimize(Solution& solution, vector<Neighborhood*>& nvect);

void basicTabuSearch(Solution& solution, vector<Neighborhood*>& nvect,
		     TabuList& tabuList, penaltyType endLimit);

void twoWaySearch(Solution& solution, vector<Neighborhood*>& nvect,
		  penaltyType endLimit);

#endif
