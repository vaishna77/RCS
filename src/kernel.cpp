//
//  kernel.hpp
//
//
//  Created by Vaishnavi Gujjula on 1/4/21.
//
//
#include "kernel.hpp"


Vec kernel::getRow(const int j, std::vector<int> col_indices) {
	int n_cols = col_indices.size();
	Vec row(n_cols);
  #pragma omp parallel for
  for(int k = 0; k < n_cols; k++) {
      row(k) = this->getMatrixEntry(j, col_indices[k]);
  }
  return row;
}

Vec kernel::getCol(const int k, std::vector<int> row_indices) {
	int n_rows = row_indices.size();
  Vec col(n_rows);
  #pragma omp parallel for
  for (int j=0; j<n_rows; ++j) {
		col(j) = this->getMatrixEntry(row_indices[j], k);
  }
  return col;
}

Vec kernel::getCol(const int n, const int k) {
  Vec col(n);
  // #pragma omp parallel for
  for (int j=0; j<n; ++j) {
		col(j) = this->getMatrixEntry(j, k);
  }
  return col;
}

Vec kernel::getRow(const int n, const int k) {
  Vec row(n);
  // #pragma omp parallel for
  for (int j=0; j<n; ++j) {
		row(j) = this->getMatrixEntry(k, j);
  }
  return row;
}

Mat kernel::getMatrix(std::vector<int> row_indices, std::vector<int> col_indices) {
	int n_rows = row_indices.size();
	int n_cols = col_indices.size();
  Mat mat(n_rows, n_cols);
  #pragma omp parallel for
  for (int j=0; j < n_rows; ++j) {
      #pragma omp parallel for
      for (int k=0; k < n_cols; ++k) {
          mat(j,k) = this->getMatrixEntry(row_indices[j], col_indices[k]);
      }
  }
  return mat;
}

Mat kernel::getMatrix(int row_start_index, int col_start_index, int row_end_index, int col_end_index) {
	Mat mat(row_end_index-row_start_index, col_end_index-col_start_index);
	// #pragma omp parallel for
	for (int j=row_start_index; j < row_end_index; ++j) {
			// #pragma omp parallel for
			for (int k=col_start_index; k < col_end_index; ++k) {
					mat(j,k) = this->getMatrixEntry(j, k);
			}
	}
	return mat;
}

void kernel::getMatrix(Mat& mat, int Nsize) {
  mat = Mat::Zero(Nsize,Nsize);
  for (size_t i = 0; i < Nsize; i++) {
    for (size_t j = 0; j < Nsize; j++) {
      mat(i,j) = this->getMatrixEntry(i, j);
    }
  }
  return;
}
