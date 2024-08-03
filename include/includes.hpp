#ifndef __INCLUDES__
#define __INCLUDES__

#include <cfloat>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <complex>
#include <Eigen/Dense>
#include <set>
#include <iterator>
#include <map>
#include <unordered_map>
#include <utility>
using namespace std;
using namespace Eigen;

typedef Eigen::Vector3d Node;
typedef Eigen::Vector3cd NodeC; //complex valued vector of 3 dimensions
typedef std::complex<double> cdouble;
#ifdef USE_DOUBLE
	typedef double dtype;
	using Mat=Eigen::MatrixXd;
	using Vec=Eigen::VectorXd;
#elif USE_COMPLEX64
	typedef std::complex<double> dtype;
	using Mat=Eigen::MatrixXcd;
	using Vec=Eigen::VectorXcd;
#endif

struct pts2D {
	double x,y;
};

struct pts3D {
	double x,y,z;
};

#ifdef USE_TWO
	const int nDims = 2;
	typedef pts2D ptsdD;
#elif USE_THREE
	const int nDims = 3;
	typedef pts3D ptsdD;
#endif

typedef cdouble (*scalarFunctionPointer)(Node&, Node&);
typedef Node (*vectorRealFunctionPointer)(Node&, Node&, Node&);
typedef NodeC (*vectorComplexFunctionPointer)(Node&, Node&, Node&);

const double epsilon = 8.854187817e-12;
const std::complex<double> I = 1i;
const double PI = 4.0*std::atan(1.0);
const double c = 3.0e8;
// const double freq = 3.0e8;
const double freq = 1.0e7;
const double lambda = c/freq;
const double kappa = 2*PI/lambda;
const double C0 = 1.0/(4*PI*epsilon);
const double mu = 4*PI*1e-7;
const double eta_0 = std::sqrt(mu/epsilon);
double ElectrostaticGreensFunction(Node& ri, Node& rj);
cdouble EFIEGreensFunction(Node& ri, Node& rj);
cdouble nonSingularGreensFunction(Node& ri, Node& rj);
double getDistance(Node& p1, Node& p2);
NodeC rhoGreensFunction(Node& ri, Node& rj, Node& r0);
NodeC rhoNonSingularGreensFunction(Node& ri, Node& rj, Node& r0);
Node rhoLaplaceFunction(Node& ri, Node& rj, Node& r3);
NodeC farFieldIntegrand(Node& ri, Node& rj, Node& r0);

struct inputsToSolverClass {
  double xmin,xmax,ymin,ymax,zmin,zmax;
  int leafSize, tol_pow, nLevelsInput, fillin_rank, fillin_tol;
  Eigen::MatrixXd locations;
};

struct inputsToKernelClass {
  std::string nodeListFilename;
  std::string patchListFilename;
	int n; // only for square plate // number of dscretizations made, i.e. h=electricalSize*lambda/n;
	double L;
};

// struct hash_pair {
//   template <class T1, class T2>
//   size_t operator()(const pair<T1, T2>& p) const {
//     auto hash1 = hash<T1>()(p.first);
//     auto hash2 = hash<T2>()(p.second);
//     return hash1 ^ hash2;
//   }
// };

#endif
