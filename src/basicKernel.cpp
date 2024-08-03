//
//  kernel.hpp
//
//
//  Created by Vaishnavi Gujjula on 1/4/21.
//
//
#include "basicKernel.hpp"


basicKernel::basicKernel(Eigen::MatrixXd& centroids) {
  this->centroids = centroids;
}

dtype basicKernel::getMatrixEntry(const unsigned i, const unsigned j) {
  if (i==j) {
    // return 100.0;
		int N = centroids.size();
		// return pow(1000*N, 1.0/2.0);
    return pow(10*N, 1.0/2.0);
  }
  else {
    // std::cout << "in getMatrixEntry" << std::endl;
    Node nodei = centroids.row(i);
    Node nodej = centroids.row(j);
    double distanceBetweenCentroids = getDistance(nodei, nodej);
    // std::cout << "1.0/distanceBetweenCentroids: " << 1.0/distanceBetweenCentroids << std::endl;
    // return 1.0;
    return exp(I*distanceBetweenCentroids)/distanceBetweenCentroids;
  }
}

// dtype basicKernel::getMatrixEntry(const unsigned i, const unsigned j) {
//   if (i==j) {
//     // return 100.0;
// 		int N = centroids.size();
// 		// return pow(1000*N, 1.0/2.0);
//     return pow(10*N, 1.0/2.0);
//   }
//   else {
//     // std::cout << "in getMatrixEntry" << std::endl;
//     Node nodei = centroids.row(i);
//     Node nodej = centroids.row(j);
//     double distanceBetweenCentroids = getDistance(nodei, nodej);
//     // std::cout << "1.0/distanceBetweenCentroids: " << 1.0/distanceBetweenCentroids << std::endl;
//     // return 1.0;
//     return 1.0/distanceBetweenCentroids;
//   }
// }
