//
//  FMM2DTree.hpp
//
//
//  Created by Vaishnavi Gujjula on 1/4/21.
//
//
#ifndef _AIFMMBox_HPP__
#define _AIFMMBox_HPP__

#include "ACA.hpp"

class FMM2DBox {
public:
	Vec L2P_weights;
	Vec P2M_weights_adjoint;
	double maxPivot_L2P;
	double maxPivot_P2M;
	bool ILActive;
	#ifdef USE_TWO
		double xmin, xmax, ymin, ymax, diameter;
	#elif USE_THREE
		double xmin, xmax, ymin, ymax, zmin, zmax, diameter;
	#endif

	int numPoints,level;
	int boxNumber;
	int parentNumber;
	int childrenNumbers[8];
	std::vector<int> interactionList;
	std::vector<int> neighborList; // considering self in neighbors
	bool Eliminated;
	int NumParticles;
	int NumMultipoles; // = NumLocals
	int NumLocals; // = NumLocals
	// int NumLocals;
	ptsdD center;
	// Vec charges;
	Vec particles;
	Vec multipoles;
	Vec locals;
	Vec particle_rhs;
	Vec multipole_rhs; //intially value will be 0
	Vec local_rhs; //intially value will be 0

	std::map<int, Mat> M2L;
	std::map<int, Mat> P2P;
	std::map<int, Mat> M2P;
	std::map<int, Mat> P2L;
  std::map<int, Mat> L2P_f;
	Mat L2M_f; //fill-ins that occur in U.transpose() equation due to elimination of x
	std::map<int, Mat> P2M_f; //fill-ins that occur in U.transpose() equation due to elimination of x
	std::map<int, Mat> M2M_f;
	Mat L2P;					//	Transfer from locals of parent to locals of children.
	Mat P2M;					//	Transfer from multipoles of 4 children to multipoles of parent.
	//	The following will be stored only at the leaf nodes
	Eigen::ColPivHouseholderQR<Mat> U_qr;
	Eigen::ColPivHouseholderQR<Mat> V_qr;
  std::vector<ptsdD> particle_loc;
  std::vector<int> chargeLocations;

	std::vector<int> incoming_chargePoints;//equivalent points {y_{k}^{B,i}}
	std::vector<int> incoming_checkPoints;//check points {x_{k}^{B,i}}
	std::vector<int> outgoing_chargePoints;
	std::vector<int> outgoing_checkPoints;
	Mat Ac;
	Mat Ar;
	int neighborListSizeAccounted;
	int indexInMatrix; // indices of starting row and column indices when the charges of bxes are in a sequential order
	FMM2DBox ();
};

#endif
