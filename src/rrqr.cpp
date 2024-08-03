#include "rrqr.hpp"

RRQR::RRQR(Mat &A, int TOL_POW1, int TOL_POW2, int partition, double inputPivot, int flag, int fillInRank) {
  this->flag = flag;
  this->A = A;
  this->tolerance1 = pow(10.0,-TOL_POW1);
  this->tolerance2 = pow(10.0,-TOL_POW2);
  this->inputPivot = inputPivot;
  this->fillInRank = fillInRank;
  // std::cout << "tolerance: " << tolerance << std::endl;
  this->partition = partition;
}

// Vec& RRQR::getNormSquareColumns(Mat& X) {
//   Vec normSquareColumns(X.cols());
//   for (size_t i = 0; i < X.cols(); i++) {
//     normSquareColumns(i) = X.col(i).norm
//   }
// }

void RRQR::printVec(std::vector<double>& vecToBePrinted) {
  for (size_t i = 0; i < vecToBePrinted.size(); i++) {
    std::cout << vecToBePrinted[i] << ", " << std::endl;
  }
  std::cout << std::endl;
  return;
}

void RRQR::printVec(std::vector<int>& vecToBePrinted) {
  for (size_t i = 0; i < vecToBePrinted.size(); i++) {
    std::cout << vecToBePrinted[i] << ", " << std::endl;
  }
  std::cout << std::endl;
  return;
}

int RRQR::mySign(double x) {
  if (x >= 0) {
    return 1;
  }
  else {
    return -1;
  }
}

void RRQR::gethouselholderVec(Vec x, int i, int j, Vec& v) {
  int n = x.size();
  v = Vec::Zero(n);
  // std::cout << "x: " << std::endl << x << std::endl;
  v.segment(i,j-i) = x.segment(i,j-i);
  // std::cout << "v: " << std::endl << v << std::endl;
  Vec shortVec = x.segment(i,j-i);
  // std::cout << "shortVec.norm(): " << shortVec.norm() << std::endl;
  // std::cout << "v(i): " << v(i) << std::endl;
  // std::cout << "mySign(x(i)): " <<mySign(x(i)) << std::endl;
  #ifdef USE_DOUBLE
    v(i) = v(i) + mySign(x(i)) * shortVec.norm();
  #elif USE_COMPLEX64
    v(i) = v(i) + std::exp(I*std::arg(x(i))) * shortVec.norm();
  #endif

  // std::cout << "v: " << std::endl << v << std::endl;
  if (v.norm() > 0) {
    v = v * std::sqrt(2)/v.norm();
  }
  // std::cout << "v: " << std::endl << v << std::endl;
}

void RRQR::colPivHouseholderQRWithPartition() {
  int n_rows = A.rows();
  int n_cols = A.cols();
  int steps;
  if (n_cols < n_rows) {
    steps = n_cols;
  }
  else {
    steps = n_rows;
  }
  Mat R = A;
  Mat Q = Mat::Identity(n_rows, n_rows);
  P = Mat::Identity(n_cols, n_cols);
  Eigen::VectorXd EigenNormSquareColumns = R.colwise().squaredNorm();
  //convert Eigen Vec to std vec
  std::vector<double> normSquareColumns(EigenNormSquareColumns.size());
  for (size_t i = 0; i < EigenNormSquareColumns.size(); i++) {
    normSquareColumns[i] = EigenNormSquareColumns(i);
  }
  int rank_partition1 = std::min(partition, n_rows);
  // std::cout << "1" << std::endl;

  // double maxPivot, maxPivotPartition1, maxPivotPartition2;
  for (size_t j = 0; j < rank_partition1; j++) {

    std::vector<int> indicesAfterSorting(normSquareColumns.size());
    std::size_t n(0);
    std::generate(std::begin(indicesAfterSorting), std::end(indicesAfterSorting), [&]{ return n++; });

    std::sort(std::begin(indicesAfterSorting)+j, std::begin(indicesAfterSorting)+partition, [&](int i1, int i2) { return normSquareColumns[i1] > normSquareColumns[i2]; } );
    std::sort(std::begin(normSquareColumns)+j, std::begin(normSquareColumns)+partition, std::greater<double>());

    Eigen::PermutationMatrix<Dynamic,Dynamic> perm(n_cols);
    VectorXi EigenIndicesAfterSorting = Eigen::Map<VectorXi>(indicesAfterSorting.data(), n_cols);
    perm.indices() = EigenIndicesAfterSorting;
    R = R*perm;
    P = P*perm;
    Vec householderVec;
    gethouselholderVec(R.col(j), j, n_rows, householderVec);
    Mat Rtemp = R - householderVec * (householderVec.adjoint()*R);
    if (j==0) {
      maxPivotPartition1 = std::abs(Rtemp(j,j));
    }
    // maxPivotPartition1 = std::max(maxPivotPartition1,inputPivot);
    maxPivot = maxPivotPartition1;
    double pivot = std::abs(Rtemp(j,j));
    if (pivot/maxPivotPartition1 < tolerance1) {
      rank_partition1 = j-1;
      break;
    }
    R = Rtemp;
    for (size_t i = j+1; i < n_rows; i++) {
      R(i,j) = 0.0;
    }
    Q = Q - (Q*householderVec) * householderVec.adjoint();
    for (size_t i = j+1; i < n_cols; i++) {
      // normSquareColumns[i] = normSquareColumns[i] - R(j,i)*R(j,i);
      normSquareColumns[i] = normSquareColumns[i] - std::abs(R(j,i))*std::abs(R(j,i));
    }
  }
  //////////////////////////////////////////////////////////////////////
  for (size_t i = rank_partition1; i < partition; i++) {
    normSquareColumns[i] = 0.0;
  }
  steps = steps-(partition-rank_partition1);
  int first = rank_partition1;
  // int rank_partition2 = std::min(steps-partition, n_rows);
  int rank_partition2 = std::min(steps-first, n_rows);

  for (size_t j = first; j < steps; j++) {
    std::vector<int> indicesAfterSorting(normSquareColumns.size());
    std::size_t n(0);
    std::generate(std::begin(indicesAfterSorting), std::end(indicesAfterSorting), [&]{ return n++; });

    std::sort(std::begin(indicesAfterSorting)+j, std::begin(indicesAfterSorting)+n_cols, [&](int i1, int i2) { return normSquareColumns[i1] > normSquareColumns[i2]; } );
    std::sort(std::begin(normSquareColumns)+j, std::begin(normSquareColumns)+n_cols, std::greater<double>());

    Eigen::PermutationMatrix<Dynamic,Dynamic> perm(n_cols);
    VectorXi EigenIndicesAfterSorting = Eigen::Map<VectorXi>(indicesAfterSorting.data(), n_cols);
    perm.indices() = EigenIndicesAfterSorting;
    R = R*perm;
    P = P*perm;
    Vec householderVec;
    gethouselholderVec(R.col(j), j, n_rows, householderVec);
    Mat Rtemp = R - householderVec * (householderVec.adjoint()*R);
    if (j==first) {
      maxPivotPartition2 = std::abs(Rtemp(first,first));
      // maxPivot = std::max(maxPivot,maxPivotPartition2);
      // maxPivot = std::max(maxPivot,inputPivot);
      maxPivot = std::max(inputPivot,maxPivotPartition2);
    }
    double pivot = std::abs(Rtemp(j,j));
    // if (j-first >= fillInRank) {//exits if the tolerance is met or the specified fillInRank is reached
    if (pivot/maxPivot < tolerance2 || j-first >= fillInRank) {//exits if the tolerance is met or the specified fillInRank is reached
    // if (pivot/maxPivot < tolerance2) { // with only tolerance as the deciding factor
        rank_partition2 = j-first;
        break;
    }
    R = Rtemp;
    for (size_t i = j+1; i < n_rows; i++) {
      R(i,j) = 0.0;
    }
    Q = Q - (Q*householderVec) * householderVec.adjoint();
    for (size_t i = j+1; i < n_cols; i++) {
      // normSquareColumns[i] = normSquareColumns[i] - R(j,i)*R(j,i);
      normSquareColumns[i] = normSquareColumns[i] - std::abs(R(j,i))*std::abs(R(j,i));
    }
  }
  int rank_A = rank_partition1 + rank_partition2;

  Q_trunc = Q.block(0,0,Q.rows(),rank_A);
  R_trunc = R.block(0,0,rank_A,R.cols());

}

double RRQR::errorInQR() {
  Mat Err = A*P - Q_trunc*R_trunc;
  return Err.norm()/A.norm();
}
