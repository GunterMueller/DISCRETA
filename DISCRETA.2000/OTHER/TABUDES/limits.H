// limits.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

// this file is needed if <limits> is not implemented in the system,
// otherwise, #include <limits>

// #include <limits>

#ifndef PACKCOVER_limits_h
#define PACKCOVER_limits_h

#include <nan.h>

template<class T>
class numeric_limits {
// this works only for integral types in most common architectures
public:
  inline static T min() { 
    T tmp(1), mn(0);
    while(tmp) { tmp <<= 1; if(tmp < mn) mn = tmp; }
    return mn;
  }
  inline static T max() {
    T tmp(1), mx(0);
    while(true) { 
      tmp = tmp << 1 | T(1);
      if(tmp <= mx) break;
      mx = tmp; 
    }
    return mx;
  }
  inline static int digits() {
    T tmp(1); 
    int ds = 0;
    while(tmp) { tmp <<= 1; ++ds; }
    return ds;
  }
  inline static int compactDigits() {
    return digits();
  }
};

class numeric_limits<bool> {
public:
  inline static bool min() { return false; }
  inline static bool max() { return true; }
  inline static int digits() { 
    return sizeof(bool) * numeric_limits<char>::digits();
  }
  inline static int compactDigits() { return 1; }
};

class numeric_limits<double> {
public:
  inline static double signaling_NaN() { return NAN; }
};

#endif
