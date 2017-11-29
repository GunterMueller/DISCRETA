// penalty.C

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

#include "packcover.H"
#include "penalty.H"

template<class T>
T abs(const T arg)
{
  if(arg < 0)
    return -arg;
  else
    return arg;
}

penaltyType* CoveringPenalty::newPenaltyTable(const int min, const int max)
{
  allocTable(min, max);
  penaltyType* ptr = table;
  for(int i = min; i < 0; i++)
    *ptr++ = abs(penaltyType(i));
  return table - min;
}

penaltyType* PackingPenalty::newPenaltyTable(const int min, const int max)
{
  allocTable(min, max);
  penaltyType* ptr = table - min + 1;
  for(int i = 1; i <= max; i++)
    *ptr++ = abs(penaltyType(i));
  return table - min;
}

penaltyType* DesignPenalty::newPenaltyTable(const int min, const int max)
{
  allocTable(min, max);
  penaltyType* ptr = table;
  for(int i = min; i <= max; i++)
    *ptr++ = abs(penaltyType(i));
  return table - min;
}

penaltyType* AsymHighDesignPenalty::newPenaltyTable
(const int min, const int max)
{
  allocTable(min, max);
  penaltyType* ptr = table;
  for(int i = min; i <= max; i++)
    *ptr++ = (i <= 0 ? abs(penaltyType(i)) : padd+penaltyType(i));
  return table - min;
}

penaltyType* AsymLowDesignPenalty::newPenaltyTable
(const int min, const int max)
{
  allocTable(min, max);
  penaltyType* ptr = table;
  for(int i = min; i <= max; i++)
    *ptr++ = (i >= 0 ? penaltyType(i) : padd+abs(penaltyType(i)));
  return table - min;
}
