#include <multiLepSearch/anaHelper.h>
#include <multiLepSearch/MT2.h>

// static double anaHelper::get_MT2(const double mVis1, const double pxVis1, const double pyVis1, const double mVis2, const double pxVis2, const double pyVis2, const double pxMiss, const double pyMiss, const double mInvis1, const double mInvis2, const double desiredPrecisionOnMT2, const bool useDeciSectionsInitially){
double anaHelper::get_mT2(const double mVis1, const double pxVis1, const double pyVis1, const double mVis2, const double pxVis2, const double pyVis2, const double pxMiss, const double pyMiss, const double mInvis1, const double mInvis2, const double desiredPrecisionOnMT2, const bool useDeciSectionsInitially){
  return asymm_mt2_lester_bisect::get_mT2(mVis1, pxVis1, pyVis1, mVis2, pxVis2, pyVis2, pxMiss, pyMiss, mInvis1, mInvis2, desiredPrecisionOnMT2, useDeciSectionsInitially);
}
