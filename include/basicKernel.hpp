//
//  kernel.hpp
//
//
//  Created by Vaishnavi Gujjula on 1/4/21.
//
//
#ifndef __BASICKERNEL_HPP__
#define __BASICKERNEL_HPP__

#include "kernel.hpp"

class basicKernel: public kernel {
private:
public:
  Eigen::MatrixXd centroids;
  basicKernel(Eigen::MatrixXd& centroids);
  dtype getMatrixEntry(const unsigned i, const unsigned j);
  ~basicKernel() {};
};

#endif
