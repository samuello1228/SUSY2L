#ifndef ANAHELPER_H
#define ANAHELPER_H
// #include <iostream>

class anaHelper{
 public:
  static double get_mT2(
    const double mVis1, const double pxVis1, const double pyVis1,
    const double mVis2, const double pxVis2, const double pyVis2,
    const double pxMiss, const double pyMiss,
    const double mInvis1, const double mInvis2,
    const double desiredPrecisionOnMT2=0,
    const bool useDeciSectionsInitially=true);
};
#endif // ANAHELPER_H
