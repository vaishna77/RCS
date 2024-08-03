#ifndef __ACA_HPP__
#define __ACA_HPP__

#include "kernel.hpp"

template <typename kerneltype>
class LowRank {
public:
	kerneltype* K;
	double tol_ACA;
	double tol_robust;
	double tolerance_or_rank;
	std::vector<int> row_indices;
	std::vector<int> col_indices;

	LowRank(kerneltype* K, double tolerance_or_rank, std::vector<int>& row_indices, std::vector<int>& col_indices);

	void maxAbsVector(const Vec& v, const std::set<int>& allowed_indices,
																dtype& max, int& index
															);

	void rookPiv(std::vector<int>& row_ind, std::vector<int>& col_ind, int& computed_rank, Mat& L, Mat& R);

	void ACA_only_nodes(std::vector<int>& row_bases, std::vector<int>& col_bases, int &computed_rank, Mat &Ac, Mat &Ar, Mat& L, Mat& R);

	void ACA_only_nodes2(std::vector<int>& row_bases, std::vector<int>& col_bases, int &computed_rank, Mat &Ac, Mat &Ar);

	void ACA_only_nodes1(std::vector<int>& row_bases, std::vector<int>& col_bases, int &computed_rank, Mat &Ac, Mat &Ar);
};

#endif
