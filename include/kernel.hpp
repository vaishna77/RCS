//
//  kernel.hpp
//
//
//  Created by Vaishnavi Gujjula on 1/4/21.
//
//
#ifndef __KERNEL_HPP__
#define __KERNEL_HPP__

#include "includes.hpp"

class kernel {
public:
  bool isTrans;		//	Checks if the kernel is translation invariant, i.e., the kernel is K(r).
	bool isHomog;		//	Checks if the kernel is homogeneous, i.e., K(r) = r^{alpha}.
	bool isLogHomog;	//	Checks if the kernel is log-homogeneous, i.e., K(r) = log(r^{alpha}).
	double alpha;		//	Degree of homogeneity of the kernel.
  double a;

  std::vector<pts2D> particles_X;
	std::vector<pts2D> particles_Y;

  kernel() {
	}

	virtual dtype getMatrixEntry(const unsigned i, const unsigned j) {
		std::cout << "virtual getInteraction" << std::endl;
		return 0.0;
	}

  Vec getRow(const int n, const int k);

  Vec getRow(const int j, std::vector<int> col_indices);

  Vec getCol(const int k, std::vector<int> row_indices);

  Vec getCol(const int n, const int k);

  Mat getMatrix(std::vector<int> row_indices, std::vector<int> col_indices);

  Mat getMatrix(int row_start_index, int col_start_index, int row_end_index, int col_end_index);

  void getMatrix(Mat& mat, int N);

  ~kernel() {};
};

#endif
