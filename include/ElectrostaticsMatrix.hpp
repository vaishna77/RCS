#ifndef __ELECTROSTATICSMATRIX_HPP__
#define __ELECTROSTATICSMATRIX_HPP__

#include "includes.hpp"
#include "Patch.hpp"
#include "kernel.hpp"

class ElectrostaticsMatrix: public kernel {
public:
  std::vector<Node> nodeList;
  std::vector<Patch> patchList;
  int N;
  std::set<int> conductorIndices;
  int numConductors;
  ElectrostaticsMatrix(double Lx, double Ly, double Lz, int Nx, int Ny);
  ElectrostaticsMatrix(double Lx, double Ly, double Lz, int Nx, int Ny, int Nz);
  ElectrostaticsMatrix(inputsToKernelClass inputsToKernel);
  void generateNodeList(std::string fileName);
  void generateNodeList(double Lx, double Ly, double Lz, int Nx, int Ny);
  void generateNodeList(double Lx, double Ly, double Lz, int Nx, int Ny, int Nz);
  void generatePatchList(int Nx, int Ny);
  void generatePatchList(std::string fileName);
  void getSizeOfMatrix();
  Node getCentroid(Node& vertex1, Node& vertex2, Node& vertex3);
  void getCentroids(Eigen::MatrixX3d& centroids);
  double getAreaOfTriangle(Node& vertex1, Node& vertex2, Node& vertex3);
  double getLargestEdge(Node& vertex1, Node& vertex2, Node& vertex3);
  double onePointGaussQuadratureRule(Node centroidTrianglei, Patch& patchj);
  double sevenPointGaussQuadratureRule(Node centroidTrianglei, double areaOfTrianglej, Node& vertex1Trianglej, Node& vertex2Trianglej, Node& vertex3Trianglej);
  double analyticIntegration(Node& centroidTrianglei, Node& vertex1Trianglej, Node& vertex2Trianglej, Node& vertex3Trianglej);
  double getMatrixEntry(const unsigned i, const unsigned j);
  void viewNodeList();
  void viewPatchList();
  // void getMatrix(Eigen::MatrixXd& mat);
  void getchargeDensityUsingLU(Eigen::VectorXd& rhs, Eigen::VectorXd& chargeDensity);
  void getchargeDensityUsingLU(Eigen::MatrixXd& rhs, Eigen::MatrixXd& chargeDensity);
  double getCapacitance(Eigen::VectorXd& chargeDensity);
  void getRhs(int conductorIndex, Eigen::VectorXd& rhs);
  void getCapacitanceMatrix(Eigen::MatrixXd chargeDensityMatrix, Eigen::MatrixXd& capacitanceMatrix);
};


#endif
