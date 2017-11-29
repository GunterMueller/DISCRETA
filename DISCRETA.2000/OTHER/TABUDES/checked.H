// checked.H

// Local search routines for solving discrete covering and packing
// problems.
//
// (c) 1999 Kari J. Nurmela

// integral types for overflow-checked, not time-critical arithmetic
// operations

// this implementation is inefficient, but should be very portable

#ifndef PACKCOVER_checked_h
#define PACKCOVER_checked_h

template<class T>
class Checked {
private:
  T val;
  Checked check(const T testT, const double testDouble, 
		const char* str = "overFlow: Checked") const {
    if(fabs(double(testT)-testDouble) >= 0.5)
      errorExit(str);
    return Checked(testT);
  }
public:
  Checked(const T& p_val) { val = p_val; }
  Checked(const Checked& orig) {val = orig.val;}
  Checked(const double f) { check(T(f), f); val = T(f); }
  Checked& operator=(const Checked& orig) { val = orig.val; return *this; }
  T operator()() { return val; }
  Checked operator+(const Checked& rhs) const { return check(val + rhs.val, 
					double(val)+double(rhs.val)); }
  Checked operator-(const Checked& rhs) const { return check(val - rhs.val, 
					double(val)-double(rhs.val)); }
  Checked operator*(const Checked& rhs) const { return check(val * rhs.val, 
					double(val)*double(rhs.val)); }
  Checked operator/(const Checked& rhs) const { return check(val / rhs.val, 
					double(val)/double(rhs.val),
					"nonintegral division: Checked"); }
  Checked operator+=(const Checked& rhs) { *this = *this + rhs;
                                            return *this; }
  Checked operator-=(const Checked& rhs) { *this = *this - rhs;
                                            return *this; }
  Checked operator*=(const Checked& rhs) { *this = *this * rhs;
                                            return *this; }
  Checked operator/=(const Checked& rhs) { *this = *this / rhs;
                                            return *this; }

  bool operator<(const Checked& rhs) const { return val < rhs.val; }
  bool operator<=(const Checked& rhs) const { return val <= rhs.val; }
  bool operator==(const Checked& rhs) const { return val == rhs.val; }
  bool operator>(const Checked& rhs) const { return val > rhs.val; }
  bool operator>=(const Checked& rhs) const { return val >= rhs.val; }
  bool operator!=(const Checked& rhs) const { return val != rhs.val; }

  Checked& operator++() { *this = *this + T(1); return *this; }
  Checked operator++(int) { Checked tmp = *this; ++*this; return tmp; }
  Checked& operator--() { *this = *this - T(1); return *this; }
  Checked operator--(int) { Checked tmp = *this; --*this; return tmp; }
};

#endif
