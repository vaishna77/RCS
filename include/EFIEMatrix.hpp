#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

#include "includes.hpp"
#include "Patch.hpp"
#include "kernel.hpp"

class Edge {
public:
  std::pair<int, int> indicesOfNodes;
  double length;
  std::vector<int> indicesOfAssociatedPatches;
};

struct compare {
    bool operator() (const Edge &edge1, const Edge &edge2) const{
      std::pair<int, int> p1,p2;
      p1 = edge1.indicesOfNodes;
      p2 = edge2.indicesOfNodes;
      return (p1.first < p2.first);
    }
};

struct edgeProperties {
  double length;
  std::vector<int> indicesOfAssociatedPatches;
};

class EFIEMatrix: public kernel {
public:
  std::vector<Node> nodeList;
  std::vector<Patch> patchList;
  std::map<int, Eigen::MatrixXd> M2L;
  std::set<int> conductorIndices;
  int numConductors;
  // std::set<std::pair<int, int> > edgeList;
  std::set<Edge, compare> edgeList;

  std::map<std::pair<int, int>, edgeProperties > edge2PatchMap;
  // std::set<std::pair<int, int> > edgeList;
  // std::vector<std::vector<int> > edge2PatchList;

  // std::unordered_set<
  //     std::pair<int, int>,
  //     boost::hash< std::pair<int, int> >
  //     > edgeList;
  int n;  // only for square plate // number of dscretizations made, i.e. h=electricalSize*lambda/n;
  double L;  // only for square plate // length of side of square plate
  // std::map<std::pair<int, int>, std::complex<double> > lookUpTable[9];  // only for square plate // looup table to store and reuse some matrix entries, from which the entire matrix can be constructed
  std::multimap<std::pair<int, int>, std::complex<double> > lookUpTable[9];  // only for square plate // looup table to store and reuse some matrix entries, from which the entire matrix can be constructed
  // std::unordered_map<std::pair<int, int>, std::complex<double>, hash_pair > lookUpTable[9];

  double getMatrixEntryFromScratchTime;// for monitoring time only
  double getMatrixEntryTime;// for monitoring time only
  int ngetMatrixEntryFromScratch;
  int ngetMatrixEntry;
  int etaAnalInteg;
  int N;
  EFIEMatrix(double Lx, double Ly, double Lz, int Nx, int Ny);
  EFIEMatrix(inputsToKernelClass inputsToKernel);
  void generateNodeList(std::string fileName);
  void generateNodeList(double Lx, double Ly, double Lz, int Nx, int Ny);
  void generatePatchList(int Nx, int Ny);
  void generatePatchList(std::string fileName);
  void generateEdgeList();
  void deleteNonRWGEdges();
  // void generateEdge2PatchMap();
  void getSizeOfMatrix();
  Node getCentroid(Node& vertex1, Node& vertex2, Node& vertex3);
  double getAreaOfTriangle(Node& vertex1, Node& vertex2, Node& vertex3);
  double getLargestEdge(Node& vertex1, Node& vertex2, Node& vertex3);

  template<scalarFunctionPointer scalarFPtr>
  cdouble onePointGaussQuadratureRule(Node& centroidTrianglei, int indexOfSourceTriangle);

  template<vectorRealFunctionPointer vectorFPtr>
  Node onePointGaussQuadratureRule(Node& centroidTrianglei, int indexOfSourceTriangle, int indexOfNonRWGVertex);

  template<scalarFunctionPointer scalarFPtr>
  cdouble sevenPointGaussQuadratureRuleScalar(Node& centroidTrianglei, int indexOfSourceTriangle);

  template<vectorComplexFunctionPointer vectorFPtr>
  NodeC sevenPointGaussQuadratureRuleComplex(Node& centroidTrianglei, int indexOfSourceTriangle, int indexOfNonRWGVertex);

  template<vectorRealFunctionPointer vectorFPtr>
  Node sevenPointGaussQuadratureRuleReal(Node& centroidTrianglei, int indexOfSourceTriangle, int indexOfNonRWGVertex);

  cdouble analyticIntegration(Node& centroidTrianglei, Node& vertex1Trianglej, Node& vertex2Trianglej, Node& vertex3Trianglej);
  cdouble getMatrixEntryFromScratch(const unsigned i, const unsigned j);
  void viewNodeList();
  void viewPatchList();
  void viewEdgeList();
  int getIndexOfThirdNode(int indexOfTestPositiveTriangle, int indexOfTestEdgeNode1, int indexOfTestEdgeNode2);
  cdouble getPhiIntegration(int indexOfTestTriangle, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex);
  cdouble getAIntegration(int indexOfTestTriangle, int indexOfSourceTriangle, int indexOfTestNonRWGVertex, int indexOfSourceNonRWGVertex);
  cdouble getScalarIntegrationGreensFunction(Node& centroid, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex);
  NodeC getVectorIntegrationGreensFunction(Node& testCentroid, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex);
  Node getVectorIntegrationSingularPart(Node& testCentroid, int indexOfSourceTriangle, int indexOfNonRWGVertex);
  NodeC getVectorIntegrationNonSingularPart(Node& testCentroid, int indexOfSourceTriangle, int indexOfNonRWGVertex);
  double analyticIntegrationLaplaceFunction(Node& centroid, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex);
  Node analyticIntegrationRhoLaplaceFunction(Node& centroid, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex);
  // void getMatrix(Eigen::MatrixXcd& mat);
  void getRhs(Node& Einc, Eigen::VectorXcd& rhs);
  cdouble getRhs(Node& Einc, int i);
  void getCoefficients(Eigen::VectorXcd& rhs, Eigen::VectorXcd& coefficients);
  void getJAtSevenPointQudraturePoints(Eigen::VectorXcd& coefficients, std::vector<Eigen::MatrixXcd>& JAtSevenPointQudraturePoints);
  void getJAtSevenPointQudraturePoints_One(Eigen::VectorXcd& coefficients, std::vector<Eigen::MatrixXcd>& JAtSevenPointQudraturePoints);
  double getRCS_SevenPointQuadrature_One(double theta, double phi, Eigen::VectorXcd& coefficients);
  double getRCS_SevenPointQuadrature(double theta, double phi, Eigen::VectorXcd& coefficients);
  double getRCS_OnePointQuadrature(double theta, double phi, Eigen::MatrixXcd& JAtCentroids);
  double getRCS(double theta, double phi, Eigen::MatrixXcd& JAtCentroids);
  double getRCSMy(double theta, double phi, Eigen::VectorXcd& coefficients);
  void createLookUpTable();
  cdouble getMatrixEntry(const unsigned i, const unsigned j);
  void getJAtCentroids(Eigen::VectorXcd& coefficients, Eigen::MatrixXcd& JAtCentroids);
  Node findProjectionOnTriangle(Node& centroidTrianglei, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex);
  void check();
};


#endif
