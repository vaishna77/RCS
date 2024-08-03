#ifndef __RRQR_HPP__
#define __RRQR_HPP__
#include "includes.hpp"

class RRQR {
public:
  Mat A;
  int fillInRank;
  int partition;
  double inputPivot;
  double maxPivot;
  double maxPivotPartition1;
  double maxPivotPartition2;
  double tolerance1;
  double tolerance2;
  Mat Q_trunc;
  Mat R_trunc;
  Mat P;
  int flag;
  RRQR(Mat& A, int TOL_POW1, int TOL_POW2, int partition, double inputPivot, int flag, int fillInRank);
  void colPivHouseholderQRWithPartition();
  void printVec(std::vector<int>& vecToBePrinted);
  void printVec(std::vector<double>& vecToBePrinted);
  void gethouselholderVec(Vec x, int i, int j, Vec& v);
  int mySign(double x);
  double errorInQR();
};

#endif
