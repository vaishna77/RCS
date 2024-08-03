//
//  FMM2DTree.hpp
//
//
//  Created by Vaishnavi Gujjula on 1/4/21.
//
//
#ifndef _AIFMMTree_HPP__
#define _AIFMMTree_HPP__

#include "ACA.hpp"
#include "AIFMMBox.hpp"
#include "rrqr.hpp"

template <typename kerneltype>
class FMM2DTree {
public:
	kerneltype* K;
	int nLevels;			//	Number of levels in the tree.
	int N;					//	Number of particles.
	double L;				//	Semi-length of the simulation box.
	std::vector<double> alpha;			//	Number of boxes at each level in the tree.
	std::vector<int> nBoxesPerLevel;			//	Number of boxes at each level in the tree.
	std::vector<double> boxSideLength;				//	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	std::vector<std::vector<FMM2DBox> >* tree;	//	The tree storing all the information.

	int nParticlesInLeafAlong1D;
  int nParticlesInLeaf;
	std::vector<double> Nodes1D;
	int TOL_POW;
	double RRQR_threshold, LU_threshold, fillin_RRQR_threshold;
	Vec b_error_check;
	double* locations;
	std::vector<std::pair<int, int> > P2L_M2P;
	std::vector<std::pair<int, int> > P2P;
	int dims;
	int startLevel;
	double eta;
	int childLevel;
	int parentLevel;
	int grandParentLevel;
	int greatGrandParentLevel;
	int nLevelsInput;
	int nChild;
	double machinePrecision, fillInTolerance, roundOffError;
	Eigen::MatrixXd particles;
	int fillin_rank, fill_in_TOL_POW;
	double timeInCompressFillin, timeQR, timeFillInQR;
	double domain_xmin, domain_xmax, domain_ymin, domain_ymax, domain_zmin, domain_zmax;
// public:
	FMM2DTree(inputsToSolverClass inputs, kerneltype* K) {
		timeFillInQR = 0.0;
		timeQR = 0.0;
		timeInCompressFillin = 0.0;
		particles = inputs.locations;
		tree = new std::vector<std::vector<FMM2DBox> >();
		machinePrecision = 1e-16;
		fillin_rank = inputs.fillin_rank;
		roundOffError = machinePrecision;
		nChild = 8;// for 3d ni. of children is 8
		nLevelsInput = inputs.nLevelsInput;
		childLevel = 0;
		parentLevel = 0;
		this->K					=	K;
		this->TOL_POW = inputs.tol_pow;
		this->fill_in_TOL_POW = inputs.fillin_tol;
		this->RRQR_threshold = pow(10,-1.0*TOL_POW);
		this->fillin_RRQR_threshold = pow(10,-1.0*(TOL_POW));
		fillInTolerance = pow(10,-1.0*TOL_POW);
		this->LU_threshold = pow(10,-1.0*TOL_POW);

		domain_xmin = inputs.xmin;
		domain_xmax = inputs.xmax;
		domain_ymin = inputs.ymin;
		domain_ymax = inputs.ymax;
		domain_zmin = inputs.zmin;
		domain_zmax = inputs.zmax;

		FMM2DBox root;
		dims = inputs.locations.cols();
		root.xmin = inputs.xmin;
		root.xmax = inputs.xmax;
		root.ymin = inputs.ymin;
		root.ymax = inputs.ymax;
		root.zmin = inputs.zmin;
		root.zmax = inputs.zmax;
		root.level = 0;
		root.numPoints = inputs.locations.rows();

		// append indices as the 4th column
		Eigen::MatrixXd locationsAndIndices(root.numPoints, 4);
		locationsAndIndices.block(0,0,root.numPoints,3) = inputs.locations;
		Eigen::VectorXd indices(root.numPoints);
		for (size_t i = 0; i < root.numPoints; i++) {
			indices(i) = i;
		}
		locationsAndIndices.block(0,3,root.numPoints,1) = indices; //{0,1,...numPoints}
		for (size_t i = 0; i < root.numPoints; i++) {
			root.chargeLocations.push_back(indices(i));
		}
		this->nLevels		=	nLevelsInput;
		std::cout << "nLevels: " << nLevels << std::endl;
		double Lx = root.xmax - root.xmin;
		double Ly = root.ymax - root.ymin;
		double Lz = root.zmax - root.zmin;
		L = std::max({Lx, Ly, Lz});

		this->TOL_POW = inputs.tol_pow;
		nBoxesPerLevel.push_back(1);
		boxSideLength.push_back(L);
		for (int k=1; k<=nLevels; ++k) {
			nBoxesPerLevel.push_back(nChild*nBoxesPerLevel[k-1]);
			boxSideLength.push_back(0.5*boxSideLength[k-1]);
		}
		this->N				=	inputs.locations.rows();
		eta = std::sqrt(3);
	}

	void createTree() {
		//	First create root and add to tree
		FMM2DBox root;
		root.boxNumber		=	0;
		root.parentNumber	=	-1;
		for (int l=0; l<nChild; ++l) {
			root.childrenNumbers[l]	=	l;
		}
		std::vector<FMM2DBox> rootLevel;
		rootLevel.push_back(root);
		tree->push_back(rootLevel);

		for (int j=1; j<=nLevels; ++j) {
			std::vector<FMM2DBox> level;
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				FMM2DBox box;
				box.boxNumber		=	k;
				box.parentNumber	=	k/nChild;
				for (int l=0; l<nChild; ++l) {
					box.childrenNumbers[l]	=	nChild*k+l;
				}
				level.push_back(box);
			}
			tree->push_back(level);
		}
	}

	bool istouch(int j, int k, int b) {
		if (k==b) {
			return true;
		}
		double roundOffError = 1e-13;
		int surfaceFlag = 0;
		if (std::abs(tree->at(j)[k].xmin - tree->at(j)[b].xmax) < roundOffError ||
				std::abs(tree->at(j)[k].xmax - tree->at(j)[b].xmin) < roundOffError) {
				surfaceFlag = 1;
		}
		if (std::abs(tree->at(j)[k].ymin - tree->at(j)[b].ymax) < roundOffError ||
				std::abs(tree->at(j)[k].ymax - tree->at(j)[b].ymin) < roundOffError) {
				surfaceFlag = 2;
		}
		if (std::abs(tree->at(j)[k].zmin - tree->at(j)[b].zmax) < roundOffError ||
				std::abs(tree->at(j)[k].zmax - tree->at(j)[b].zmin) < roundOffError) {
				surfaceFlag = 3;
		}
		if (surfaceFlag == 0) {
			return false;
		}
		else {
			if (surfaceFlag == 3) {
				if (!((tree->at(j)[b].xmax < tree->at(j)[k].xmin || tree->at(j)[k].xmax < tree->at(j)[b].xmin) || (tree->at(j)[b].ymax < tree->at(j)[k].ymin || tree->at(j)[k].ymax < tree->at(j)[b].ymin))) {
					return true;
				}
				else {
					return false;
				}
			}

			else if (surfaceFlag == 2) {
				if (!((tree->at(j)[b].xmax < tree->at(j)[k].xmin || tree->at(j)[k].xmax < tree->at(j)[b].xmin) || (tree->at(j)[b].zmax < tree->at(j)[k].zmin || tree->at(j)[k].zmax < tree->at(j)[b].zmin))) {
					return true;
				}
				else {
					return false;
				}
			}

			else {
				if (!((tree->at(j)[b].ymax < tree->at(j)[k].ymin || tree->at(j)[k].ymax < tree->at(j)[b].ymin) || (tree->at(j)[b].zmax < tree->at(j)[k].zmin || tree->at(j)[k].zmax < tree->at(j)[b].zmin))) {
					return true;
				}
				else {
					return false;
				}
			}
		}
	}

	void check3() {
		bool touch = istouch(5,4,20);
		std::cout << "istouch(5,4,20): " << touch << std::endl;
	}

	double distanceBetweenCenters(int jk, int k, int jb, int b) {
		double distCentersX = tree->at(jk)[k].center.x - tree->at(jb)[b].center.x;
		double distCentersY = tree->at(jk)[k].center.y - tree->at(jb)[b].center.y;
		double distCentersZ = tree->at(jk)[k].center.z - tree->at(jb)[b].center.z;
		double distCenters = distCentersX*distCentersX + distCentersY*distCentersY + distCentersZ*distCentersZ;
		return std::sqrt(distCenters);
	}


	void assign_Child_Interaction(int j, int k, int whichChild) {
		int child = nChild*k+whichChild; //tree->at(j+1)[2*k]
		for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) {
			int nn = tree->at(j)[k].neighborList[n]; //neighbor number
			for (size_t c = 0; c < nChild; c++) { //iterate through children of parent's neighbors
				int other = nn*nChild+c;
				admissibilityCondition(j+1, child, other); // test.cpp works fine = cube with volume data
			}
		}
	}

	bool admissibilityCondition(int j, int child, int other) {
		// j is parent level of child and other
		// finds if other belongs to the IL or N of child
		double distCenters = distanceBetweenCenters(j, child, j, other);
		if (distCenters > eta*boxSideLength[j] + boxSideLength[nLevels]/2.0) {
			tree->at(j)[child].interactionList.push_back(other);
			return true;
		}
		else {
			tree->at(j)[child].neighborList.push_back(other);
			return false;
		}
	}

	//	Assigns the interactions for the children of a box
	void assign_Box_Interactions(int j, int k) {
		for (size_t i = 0; i < nChild; i++) {
			assign_Child_Interaction(j,k,i);
		}
	}

	//	Assigns the interactions for the children all boxes at a given level
	void assign_Level_Interactions(int j) {
		// #pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			assign_Box_Interactions(j,k);
		}
	}

	//	Assigns the interactions for the children all boxes in the tree
	void assign_Tree_Interactions() {
		tree->at(0)[0].neighborList.push_back(0);
		for (int j=0; j<nLevels; ++j) {
			assign_Level_Interactions(j);
		}
		getStartLevel();
	}


	void getStartLevel() {
		int flag = -1;
		for (size_t j = 0; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				if (flag != 0 && tree->at(j)[k].interactionList.size() != 0) {
					startLevel = j;
					std::cout << "Interaction list starts to exist at level: " << j << std::endl;
					flag = 0;
				}
			}
		}
	}

	void checkILandNSize() {
		int *avgILSize = new int[nLevels+1];
		int *avgNSize = new int[nLevels+1];
		for (size_t j = 0; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				avgNSize[j] += tree->at(j)[k].neighborList.size();
				avgILSize[j] += tree->at(j)[k].interactionList.size();
			}
			avgNSize[j] = avgNSize[j]/nBoxesPerLevel[j];
			avgILSize[j] = avgILSize[j]/nBoxesPerLevel[j];
			std::cout << "j: " << j << "	avgILSize: " << avgILSize[j] << "	avgNSize: " << avgNSize[j] << std::endl;
		}
		delete avgILSize;
		delete avgNSize;
	}

	void checkSize() {
		for (size_t i = 0; i <= nLevels; i++) {
			std::cout << "boxSideLength[i]: " << boxSideLength[i] << std::endl;
		}
		for (size_t j = 0; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				 	std::cout << "j: " << j << "	k: " << k << "	chargeSize: " << tree->at(j)[k].chargeLocations.size() << "	incoming_checkPoints: " << tree->at(j)[k].incoming_checkPoints.size() << "	outgoing_chargePoints: " << tree->at(j)[k].outgoing_chargePoints.size() << std::endl;
					std::cout << "center: " << tree->at(j)[k].center.x << "," << tree->at(j)[k].center.y << "," << tree->at(j)[k].center.z << std::endl;
					std::cout << std::endl << "----------------------------" << std::endl;
			}
		}
	}

	void check() {
		for (size_t i = 0; i <= nLevels; i++) {
			std::cout << "boxSideLength[i]: " << boxSideLength[i] << std::endl;
		}
		for (size_t j = 0; j <= 1; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				 	std::cout << "j: " << j << "	k: " << k << "	chargeSize: " << tree->at(j)[k].chargeLocations.size() << "	numPoints: " << tree->at(j)[k].numPoints << std::endl;
					std::cout << "center: " << tree->at(j)[k].center.x << "," << tree->at(j)[k].center.y << "," << tree->at(j)[k].center.z << std::endl;
					std::cout << std::endl << "neighborList.size(): " << tree->at(j)[k].neighborList.size() << std::endl;
					std::cout << "interactionList.size(): " << tree->at(j)[k].interactionList.size() << std::endl;
				 	std::cout << "neighborList: ";
				 	for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) {
						std::cout << tree->at(j)[k].neighborList[n] << ",";
					}
					std::cout << std::endl << "interactionList: ";
					for (size_t n = 0; n < tree->at(j)[k].interactionList.size(); n++) {
						std::cout << tree->at(j)[k].interactionList[n] << ",";
					}
					std::cout << std::endl << "----------------------------" << std::endl;
			}
		}
	}

	void check22() {
		int j=7;
		int k=10;
		std::cout << std::endl << "neighborlist of " << j << "," << k << "----------------" << std::endl;
		for (size_t i = 0; i < tree->at(j)[k].neighborList.size(); i++) {
			int n = tree->at(j)[k].neighborList[i];
			std::cout << n << ",";
		}
		std::cout << std::endl;
		std::cout << std::endl << "IL of " << j << "," << k << "----------------" << std::endl;
		for (size_t i = 0; i < tree->at(j)[k].interactionList.size(); i++) {
			int n = tree->at(j)[k].interactionList[i];
			std::cout << n << ",";
		}
		std::cout << std::endl;
	}

	void assign_Center_Location() {
		int J;
		tree->at(0)[0].center.x	=	(domain_xmin+domain_xmax)/2.0;
		tree->at(0)[0].center.y	=	(domain_ymin+domain_ymax)/2.0;
		tree->at(0)[0].center.z	=	(domain_zmin+domain_zmax)/2.0;
		for (int j=0; j<nLevels; ++j) {
			J	=	j+1;
			double shift	=	0.5*boxSideLength[J];
			// #pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				tree->at(J)[8*k].center.x		=	tree->at(j)[k].center.x-shift;
				tree->at(J)[8*k+1].center.x	=	tree->at(j)[k].center.x+shift;
				tree->at(J)[8*k+2].center.x	=	tree->at(j)[k].center.x+shift;
				tree->at(J)[8*k+3].center.x	=	tree->at(j)[k].center.x-shift;
				tree->at(J)[8*k+4].center.x	=	tree->at(j)[k].center.x-shift;
				tree->at(J)[8*k+5].center.x	=	tree->at(j)[k].center.x+shift;
				tree->at(J)[8*k+6].center.x	=	tree->at(j)[k].center.x+shift;
				tree->at(J)[8*k+7].center.x	=	tree->at(j)[k].center.x-shift;

				tree->at(J)[8*k].center.y		=	tree->at(j)[k].center.y-shift;
				tree->at(J)[8*k+1].center.y	=	tree->at(j)[k].center.y-shift;
				tree->at(J)[8*k+2].center.y	=	tree->at(j)[k].center.y+shift;
				tree->at(J)[8*k+3].center.y	=	tree->at(j)[k].center.y+shift;
				tree->at(J)[8*k+4].center.y	=	tree->at(j)[k].center.y-shift;
				tree->at(J)[8*k+5].center.y	=	tree->at(j)[k].center.y-shift;
				tree->at(J)[8*k+6].center.y	=	tree->at(j)[k].center.y+shift;
				tree->at(J)[8*k+7].center.y	=	tree->at(j)[k].center.y+shift;

				tree->at(J)[8*k].center.z		=	tree->at(j)[k].center.z-shift;
				tree->at(J)[8*k+1].center.z	=	tree->at(j)[k].center.z-shift;
				tree->at(J)[8*k+2].center.z	=	tree->at(j)[k].center.z-shift;
				tree->at(J)[8*k+3].center.z	=	tree->at(j)[k].center.z-shift;
				tree->at(J)[8*k+4].center.z	=	tree->at(j)[k].center.z+shift;
				tree->at(J)[8*k+5].center.z	=	tree->at(j)[k].center.z+shift;
				tree->at(J)[8*k+6].center.z	=	tree->at(j)[k].center.z+shift;
				tree->at(J)[8*k+7].center.z	=	tree->at(j)[k].center.z+shift;
			}
		}
	}

void reorder(Vec &potential) {
	Vec potentialTemp = potential;
	int start = 0;
	for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
		for (size_t i = 0; i < tree->at(nLevels)[k].chargeLocations.size(); i++) {
			int index = tree->at(nLevels)[k].chargeLocations[i];
			potential(index) = potentialTemp(start);
			start++;
		}
	}
}

	int getMaxRank() {
		int max = 0;
		for (int j = nLevels; j >= startLevel; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				if (max < tree->at(j)[k].multipoles.size()) {
					max = tree->at(j)[k].multipoles.size();
				}
			}
		}
		return max;
	}

	int getAvgLeafSize() {
		int avg = 0;
		int nBox = 0;
		int j = nLevels;
		for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
			if (tree->at(j)[k].chargeLocations.size() > 0) {
				nBox++;
				avg += tree->at(j)[k].chargeLocations.size();
			}
		}
		return avg/nBox;
	}

	int getAvgRank() {
		int avg = 0;
		int nBox = 0;
		for (int j = nLevels; j >= startLevel; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				if (tree->at(j)[k].multipoles.size() > 0) {
					nBox++;
					avg += tree->at(j)[k].multipoles.size();
				}
			}
		}
		return avg/nBox;
	}

	int getAvgLocals() {
		int avg = 0;
		int nBox = 0;
		for (int j = nLevels; j >= startLevel; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				if (tree->at(j)[k].L2P.cols() > 0) {
					nBox++;
					avg += tree->at(j)[k].L2P.cols();
				}
			}
		}
		return avg/nBox;
	}

	int getAvgMultipoles() {
		int avg = 0;
		int nBox = 0;
		for (int j = nLevels; j >= startLevel; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				if (tree->at(j)[k].P2M.rows() > 0) {
					nBox++;
					avg += tree->at(j)[k].P2M.rows();
				}
			}
		}
		return avg/nBox;
	}

	void check55() {
		for (int j = nLevels; j >= 2; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				std::cout << "j: " << j << "	k: " << k << "	cL: " << tree->at(j)[k].chargeLocations.size() << std::endl;
				for (size_t i = 0; i < tree->at(j)[k].chargeLocations.size(); i++) {
					std::cout << tree->at(j)[k].chargeLocations[i] << ",";
				}
				std::cout << std::endl;
				// tree->at(j)[k].chargeLocations.clear();
			}
		}
		for (int j = nLevels; j >= 2; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				std::cout << "j: " << j << "	k: " << k << "	pr: " << tree->at(j)[k].particle_rhs.size() << std::endl;
				for (size_t i = 0; i < tree->at(j)[k].particle_rhs.size(); i++) {
					std::cout << tree->at(j)[k].particle_rhs[i] << ",";
				}
				std::cout << std::endl;
			}
		}
	}

	void assignChargeLocations() {
		for (size_t i = 0; i < N; i++) {
			tree->at(0)[0].chargeLocations.push_back(i);
		}
		for (size_t j = 0; j < nLevels; j++) { //assign particles to its children
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				int J = j+1;
				int Kp = 8*k;
				for (size_t i = 0; i < tree->at(j)[k].chargeLocations.size(); i++) {
					int index = tree->at(j)[k].chargeLocations[i];
					if (particles(index,2) <= tree->at(j)[k].center.z) { //children 0,1,2,3
						if (particles(index,0) <= tree->at(j)[k].center.x) { //children 0,3
							if (particles(index,1) <= tree->at(j)[k].center.y) { //child 0
								tree->at(J)[Kp].chargeLocations.push_back(index);
							}
							else { //child 3
								tree->at(J)[Kp+3].chargeLocations.push_back(index);
							}
						}
						else { //children 1,2
							if (particles(index,1) <= tree->at(j)[k].center.y) { //child 1
								tree->at(J)[Kp+1].chargeLocations.push_back(index);
							}
							else { //child 2
								tree->at(J)[Kp+2].chargeLocations.push_back(index);
							}
						}
					}
					else {//children 4,5,6,7
						if (particles(index,0) <= tree->at(j)[k].center.x) { //children 4,7
							if (particles(index,1) <= tree->at(j)[k].center.y) { //child 4
								tree->at(J)[Kp+4].chargeLocations.push_back(index);
							}
							else { //child 7
								tree->at(J)[Kp+7].chargeLocations.push_back(index);
							}
						}
						else { //children 5,6
							if (particles(index,1) <= tree->at(j)[k].center.y) { //child 5
								tree->at(J)[Kp+5].chargeLocations.push_back(index);
							}
							else { //child 6
								tree->at(J)[Kp+6].chargeLocations.push_back(index);
							}
						}
					}
				}
			}
		}
	}

	void assignNonLeafChargeLocations() {
		for (int j = nLevels-1; j >= 1; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				tree->at(j)[k].chargeLocations.clear();
				for (size_t i = 0; i < nChild; i++) {
					tree->at(j)[k].chargeLocations.insert(tree->at(j)[k].chargeLocations.end(), tree->at(j+1)[nChild*k+i].chargeLocations.begin(), tree->at(j+1)[nChild*k+i].chargeLocations.end());
				}
			}
		}
	}

	void assign_Leaf_rhs(Vec &charges) {
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			tree->at(nLevels)[k].particle_rhs	=	Vec::Zero(tree->at(nLevels)[k].chargeLocations.size());
			for (size_t i = 0; i < tree->at(nLevels)[k].chargeLocations.size(); i++) {
				int index = tree->at(nLevels)[k].chargeLocations[i];
				tree->at(nLevels)[k].particle_rhs[i]	=	charges[index];
			}
		}
	}


	void get_L2P_P2M_box(int j, int k) {
		if (tree->at(j)[k].incoming_checkPoints.size() < tree->at(j)[k].outgoing_chargePoints.size()) {
			tree->at(j)[k].outgoing_checkPoints.erase(tree->at(j)[k].outgoing_checkPoints.begin()+tree->at(j)[k].incoming_checkPoints.size(), tree->at(j)[k].outgoing_checkPoints.end());
			tree->at(j)[k].outgoing_chargePoints.erase(tree->at(j)[k].outgoing_chargePoints.begin()+tree->at(j)[k].incoming_checkPoints.size(), tree->at(j)[k].outgoing_chargePoints.end());
		}
		else {
			tree->at(j)[k].incoming_checkPoints.erase(tree->at(j)[k].incoming_checkPoints.begin()+tree->at(j)[k].outgoing_chargePoints.size(), tree->at(j)[k].incoming_checkPoints.end());
			tree->at(j)[k].incoming_chargePoints.erase(tree->at(j)[k].incoming_chargePoints.begin()+tree->at(j)[k].outgoing_chargePoints.size(), tree->at(j)[k].incoming_chargePoints.end());
		}
		if (tree->at(j)[k].ILActive == true) {
			Mat temp1 = tree->at(j)[k].Ac.block(0,0,tree->at(j)[k].Ac.rows(), tree->at(j)[k].incoming_checkPoints.size());
			tree->at(j)[k].L2P = Mat::Zero(tree->at(j)[k].Ac.rows(),tree->at(j)[k].incoming_checkPoints.size());
			if (tree->at(j)[k].incoming_checkPoints.size() > 0) {
				Mat D1 = K->getMatrix(tree->at(j)[k].incoming_checkPoints, tree->at(j)[k].incoming_chargePoints);
				Eigen::PartialPivLU<Mat> D1_lu = Eigen::PartialPivLU<Mat>(D1.adjoint());
				tree->at(j)[k].L2P = D1_lu.solve(temp1.adjoint()).adjoint();
			}
			Mat temp2 = tree->at(j)[k].Ar.block(0,0,tree->at(j)[k].outgoing_checkPoints.size(),tree->at(j)[k].Ar.cols());
			tree->at(j)[k].P2M = Mat(tree->at(j)[k].outgoing_checkPoints.size(),tree->at(j)[k].Ar.cols());
			if (tree->at(j)[k].outgoing_checkPoints.size() > 0) {
				Mat D2 = K->getMatrix(tree->at(j)[k].outgoing_checkPoints, tree->at(j)[k].outgoing_chargePoints);
				Eigen::PartialPivLU<Mat> D2_lu = Eigen::PartialPivLU<Mat>(D2);
				tree->at(j)[k].P2M = D2_lu.solve(temp2);
			}
		}
		else {
			tree->at(j)[k].L2P = Mat::Identity(tree->at(j)[k].incoming_checkPoints.size(), tree->at(j)[k].incoming_checkPoints.size());
			tree->at(j)[k].P2M = Mat::Identity(tree->at(j)[k].outgoing_chargePoints.size(), tree->at(j)[k].outgoing_chargePoints.size());
		}
		if (tree->at(j)[k].Ac.size() > 0) {
			Mat truncatedAc = tree->at(j)[k].Ac.block(0,0,tree->at(j)[k].Ac.rows(), tree->at(j)[k].incoming_checkPoints.size());
			Eigen::ColPivHouseholderQR<Mat> qr_Ac(truncatedAc.rows(), truncatedAc.cols());
			qr_Ac.compute(truncatedAc);
			Vec R_Ac_diagonal = qr_Ac.matrixQR().diagonal();
			tree->at(j)[k].maxPivot_L2P = std::abs(R_Ac_diagonal(0));
		}
		if (tree->at(j)[k].Ar.size() > 0) {
			Mat truncatedAr = tree->at(j)[k].Ar.block(0,0,tree->at(j)[k].outgoing_checkPoints.size(),tree->at(j)[k].Ar.cols());
			Eigen::ColPivHouseholderQR<Mat> qr_Ar_adjoint(truncatedAr.cols(), truncatedAr.rows());
			qr_Ar_adjoint.compute(truncatedAr.adjoint());
			Vec R_Ar_diagonal = qr_Ar_adjoint.matrixQR().diagonal();
			tree->at(j)[k].maxPivot_P2M = std::abs(R_Ar_diagonal(0));//weights of P2M_adjoint
		}
	}

	void write_L2P_P2M_box(int j, int k) {
		std::string filename = "L2P_" + std::to_string(k);
		std::ofstream myfile;
		myfile.open(filename.c_str());
		myfile << tree->at(j)[k].L2P << std::endl;

		std::string filename2 = "P2M_" + std::to_string(k);
		std::ofstream myfile2;
		myfile2.open(filename2.c_str());
		myfile2 << tree->at(j)[k].P2M << std::endl;
	}


	void check1() {
		for (int j=nLevels; j>=2; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				std::cout << "j: " << j << "	k: " << k << "	L2P: " << tree->at(j)[k].L2P.rows() << ", " << tree->at(j)[k].L2P.cols() << "	P2M: " << tree->at(j)[k].P2M.rows() << ", " << tree->at(j)[k].P2M.cols() << "	D: " << tree->at(j)[k].L2P.rows() - tree->at(j)[k].P2M.cols() << ", " << tree->at(j)[k].L2P.cols() - tree->at(j)[k].P2M.rows() << std::endl;
			}
		}
	}

	void get_L2P_P2M_level(int j) {
		// #pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			get_L2P_P2M_box(j, k);
		}
	}

	void write_extended_sparse_matrix() {
		int sizeA = 0;
		int j = nLevels;
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			sizeA += tree->at(j)[k].chargeLocations.size() + tree->at(j)[k].L2P.cols() + tree->at(j)[k].P2M.rows();
		}
		Mat A = Mat::Zero(sizeA, sizeA);
		std::vector<int> rowIndex(nBoxesPerLevel[j]);
		std::vector<int> colIndex(nBoxesPerLevel[j]);
		rowIndex[0] = 0;
		colIndex[0] = 0;
		for (int k=1; k<nBoxesPerLevel[j]; ++k) {
			rowIndex[k] = rowIndex[k-1] + tree->at(j)[k-1].chargeLocations.size() + tree->at(j)[k-1].P2M.rows();
			colIndex[k] = colIndex[k-1] + tree->at(j)[k-1].chargeLocations.size() + tree->at(j)[k-1].L2P.cols();
		}
		// P2Ps
		{
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			int row = rowIndex[k];
			for (size_t i = 0; i < tree->at(j)[k].neighborList.size(); i++) {
				int n = tree->at(j)[k].neighborList[i];
				int col = colIndex[n];
				A.block(row, col, tree->at(j)[k].P2P[n].rows(), tree->at(j)[k].P2P[n].cols()) = tree->at(j)[k].P2P[n];
			}
		}
	}
		// M2Ls
		{
		int col_offset = 0;
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			col_offset += tree->at(j)[k].chargeLocations.size() + tree->at(j)[k].L2P.cols();
		}
		int row_offset = 0;
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			row_offset += tree->at(j)[k].chargeLocations.size() + tree->at(j)[k].P2M.rows();
		}
		std::vector<int> rowIndexM2L(nBoxesPerLevel[j]);
		std::vector<int> colIndexM2L(nBoxesPerLevel[j]);
		rowIndexM2L[0] = 0;
		colIndexM2L[0] = 0;
		for (int k=1; k<nBoxesPerLevel[j]; ++k) {
			rowIndexM2L[k] = rowIndexM2L[k-1] + tree->at(j)[k-1].L2P.cols();
			colIndexM2L[k] = colIndexM2L[k-1] + tree->at(j)[k-1].P2M.rows();
		}

		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			int row = row_offset;
			for (size_t i = 0; i < tree->at(j)[k].interactionList.size(); i++) {
				int IL = tree->at(j)[k].interactionList[i];
				int col = col_offset + colIndexM2L[IL];
				A.block(row, col, tree->at(j)[k].M2L[IL].rows(), tree->at(j)[k].M2L[IL].cols()) = tree->at(j)[k].M2L[IL];
			}
			row_offset += tree->at(j)[k].L2P.cols();
		}
	}
		// L2Ps
		{
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			int col = colIndex[k] + tree->at(j)[k].chargeLocations.size();
			int row = rowIndex[k];
			A.block(row, col, tree->at(j)[k].L2P.rows(), tree->at(j)[k].L2P.cols()) = tree->at(j)[k].L2P;
		}
	}
		// P2Ms
		{
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			int col = colIndex[k];
			int row = rowIndex[k] + tree->at(j)[k].chargeLocations.size();
			A.block(row, col, tree->at(j)[k].P2M.rows(), tree->at(j)[k].P2M.cols()) = tree->at(j)[k].P2M;
		}
	}
		// Is
		{
		int col_offset = 0;
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			col_offset += tree->at(j)[k].chargeLocations.size() + tree->at(j)[k].L2P.cols();
		}
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			int col = col_offset;
			int row = rowIndex[k] + tree->at(j)[k].chargeLocations.size();
			A.block(row, col, tree->at(j)[k].P2M.rows(), tree->at(j)[k].P2M.rows()) = -Mat::Identity(tree->at(j)[k].P2M.rows(), tree->at(j)[k].P2M.rows());
			col_offset += tree->at(j)[k].P2M.rows();
		}
	}
		// Is
	{
		int row_offset = 0;
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			row_offset += tree->at(j)[k].chargeLocations.size() + tree->at(j)[k].P2M.rows();
		}
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			int row = row_offset;
			int col = colIndex[k] + tree->at(j)[k].chargeLocations.size();
			A.block(row, col, tree->at(j)[k].L2P.cols(), tree->at(j)[k].L2P.cols()) = -Mat::Identity(tree->at(j)[k].L2P.cols(), tree->at(j)[k].L2P.cols());
			row_offset += tree->at(j)[k].L2P.cols();
		}
	}

	std::string filename = "extendedSparseMatrix_" + std::to_string(N);
	std::ofstream myfile;
	myfile.open(filename.c_str());
	myfile << A << std::endl;

}

void write_extended_rhs() {
	int sizeA = 0;
	int j = nLevels;
	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
		sizeA += tree->at(j)[k].chargeLocations.size() + tree->at(j)[k].L2P.cols() + tree->at(j)[k].P2M.rows();
	}
	Vec extendedRhs = Vec::Zero(sizeA);
	int index = 0;
	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
		extendedRhs.segment(index, tree->at(nLevels)[k].particle_rhs.size()) = tree->at(nLevels)[k].particle_rhs;
		index += tree->at(nLevels)[k].particle_rhs.size() + tree->at(nLevels)[k].P2M.rows();
	}
	std::string filename2 = "extendedRhs_" + std::to_string(N);
	std::ofstream myfile2;
	myfile2.open(filename2.c_str());
	myfile2 << extendedRhs << std::endl;
}


void write_NumParticles() {
	int j = nLevels;
	Vec numParticles(nBoxesPerLevel[j]);
	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
		numParticles(k) = tree->at(nLevels)[k].chargeLocations.size();
	}
	std::string filename2 = "numParticles_" + std::to_string(N) + ".txt";
	std::ofstream myfile2;
	myfile2.open(filename2.c_str());
	myfile2 << numParticles << std::endl;
}

void write_NumMultipoles() {
	int j = nLevels;
	Vec numMultipoles(nBoxesPerLevel[j]);
	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
		numMultipoles(k) = tree->at(nLevels)[k].P2M.rows();
	}
	std::string filename2 = "numMultipoles_" + std::to_string(N) + ".txt";
	std::ofstream myfile2;
	myfile2.open(filename2.c_str());
	myfile2 << numMultipoles << std::endl;
}

	void write_L2P_P2M_level() {
		int j = nLevels;
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			write_L2P_P2M_box(j, k);
		}
	}

	void getNodes() {
		for (int j=nLevels; j>=startLevel; j--) {
			getNodes_outgoing_level(j);
			getNodes_incoming_level(j);
			get_L2P_P2M_level(j);
		}
	}

	void getNodes_outgoing_level(int j) { //LFR; box interactions
    // #pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			getNodes_outgoing_box(j, k);
		}
	}

	void getNodes_incoming_level(int j) { //LFR; box interactions
    // #pragma omp parallel for
    for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			getNodes_incoming_box(j, k);
		}
	}

	void getParticlesFromChildren_incoming_row(int j, int k, std::vector<int>& searchNodes) {
		if (j==nLevels) {
			searchNodes.insert(searchNodes.end(), tree->at(j)[k].chargeLocations.begin(), tree->at(j)[k].chargeLocations.end());
		}
		else {
			int J = j+1;
			for (int c = 0; c < nChild; c++) {
				searchNodes.insert(searchNodes.end(), tree->at(J)[nChild*k+c].incoming_checkPoints.begin(), tree->at(J)[nChild*k+c].incoming_checkPoints.end());
			}
		}
	}

	void getParticlesFromChildren_incoming_col(int j, int k, std::vector<int>& searchNodes) {
		if (j==nLevels) {
			searchNodes.insert(searchNodes.end(), tree->at(j)[k].chargeLocations.begin(), tree->at(j)[k].chargeLocations.end());
		}
		else {
			int J = j+1;
			for (int c = 0; c < nChild; c++) {
				searchNodes.insert(searchNodes.end(), tree->at(J)[nChild*k+c].outgoing_chargePoints.begin(), tree->at(J)[nChild*k+c].outgoing_chargePoints.end());
			}
		}
	}

	void getParticlesFromChildren_outgoing_row(int j, int k, std::vector<int>& searchNodes) {
		if (j==nLevels) {
			searchNodes.insert(searchNodes.end(), tree->at(j)[k].chargeLocations.begin(), tree->at(j)[k].chargeLocations.end());
		}
		else {
			int J = j+1;
			for (int c = 0; c < nChild; c++) {
				searchNodes.insert(searchNodes.end(), tree->at(J)[nChild*k+c].incoming_checkPoints.begin(), tree->at(J)[nChild*k+c].incoming_checkPoints.end());
			}
		}
	}

	void getParticlesFromChildren_outgoing_col(int j, int k, std::vector<int>& searchNodes) {
		if (j==nLevels) {
			searchNodes.insert(searchNodes.end(), tree->at(j)[k].chargeLocations.begin(), tree->at(j)[k].chargeLocations.end());
		}
		else {
			int J = j+1;
			for (int c = 0; c < nChild; c++) {
				searchNodes.insert(searchNodes.end(), tree->at(J)[nChild*k+c].outgoing_chargePoints.begin(), tree->at(J)[nChild*k+c].outgoing_chargePoints.end());
			}
		}
	}

	void getNodes_outgoing_box(int j, int k) {
		int n_rows, n_cols, ComputedRank;
    std::vector<int> boxA_Nodes;
		getParticlesFromChildren_outgoing_col(j, k, boxA_Nodes);
		std::vector<int> IL_Nodes;//indices
		for(int in=0; in<tree->at(j)[k].interactionList.size(); ++in) {
			int kIL = tree->at(j)[k].interactionList[in];
			std::vector<int> chargeLocations;
			getParticlesFromChildren_outgoing_row(j, kIL, chargeLocations);
			IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
		}
		n_rows = IL_Nodes.size();
		n_cols = boxA_Nodes.size();
		if (n_cols > 0 && n_rows > 0) {
			std::vector<int> row_indices, col_indices;
			row_indices = IL_Nodes;
			col_indices = boxA_Nodes;//object of base class G_LowRank
			std::vector<int> row_bases, col_bases;
			Mat Ac;
			LowRank<kerneltype>* LR		=	new LowRank<kerneltype>(K, pow(10,-TOL_POW), row_indices, col_indices);
			// LR->rookPiv(row_bases, col_bases, ComputedRank, Ac, tree->at(j)[k].Ar);
				LR->ACA_only_nodes2(row_bases, col_bases, ComputedRank, Ac, tree->at(j)[k].Ar);

			int minN = n_rows;
			if (n_rows > n_cols) {
				minN = n_cols;
			}
			if(ComputedRank > 0) {
				for (int r = 0; r < row_bases.size(); r++) {
					tree->at(j)[k].outgoing_checkPoints.push_back(IL_Nodes[row_bases[r]]);
				}
				for (int c = 0; c < col_bases.size(); c++) {
					tree->at(j)[k].outgoing_chargePoints.push_back(boxA_Nodes[col_bases[c]]);
				}


				// Mat D = K->getMatrix(tree->at(j)[k].outgoing_checkPoints, tree->at(j)[k].outgoing_chargePoints);
				// Eigen::PartialPivLU<Mat> D_lu = Eigen::PartialPivLU<Mat>(D);
				// Mat mat = K->getMatrix(row_indices,col_indices);
				// Mat Err = mat - Ac*D_lu.solve(tree->at(j)[k].Ar);
				// double err = Err.norm()/mat.norm();
				// std::cout << "O j: " << j << "	k: " << k << "	rank: " << ComputedRank << "	err: " << err << std::endl;

			}
			// else {
			// 	tree->at(j)[k].P2M = Mat(tree->at(j)[k].outgoing_checkPoints.size(),tree->at(j)[k].Ar.cols());
			// }
			// std::cout << "j: " << j << "	k: "<< k << "	ComputedRank: " << ComputedRank << std::endl;
		}
		else if (n_rows == 0) {
			tree->at(j)[k].ILActive = false;
			tree->at(j)[k].outgoing_chargePoints = boxA_Nodes;
		}
	}

	void writeHmatrix(Mat& mat) {
		// assuming nlevels = 2
		int j = nLevels;
		mat = Mat(N,N);
		std::vector<int> rowIndex(nBoxesPerLevel[j]);
		std::vector<int> colIndex(nBoxesPerLevel[j]);
		for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
			rowIndex[0] = 0;
			colIndex[0] = 0;
			for (int k=1; k<nBoxesPerLevel[j]; ++k) {
				rowIndex[k] = rowIndex[k-1] + tree->at(j)[k-1].chargeLocations.size();
				colIndex[k] = colIndex[k-1] + tree->at(j)[k-1].chargeLocations.size();
			}
		}

		// Neighbors
		{
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				int row = rowIndex[k];
				for (size_t i = 0; i < tree->at(j)[k].neighborList.size(); i++) {
					int n = tree->at(j)[k].neighborList[i];
					int col = colIndex[n];
					mat.block(row, col, tree->at(j)[k].chargeLocations.size(), tree->at(j)[n].chargeLocations.size()) = tree->at(j)[k].P2P[n];
				}
			}
		}

		// IL
		{
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				int row = rowIndex[k];
				for (size_t i = 0; i < tree->at(j)[k].interactionList.size(); i++) {
					int n = tree->at(j)[k].interactionList[i];
					int col = colIndex[n];
					mat.block(row, col, tree->at(j)[k].chargeLocations.size(), tree->at(j)[n].chargeLocations.size()) = tree->at(j)[k].L2P*tree->at(j)[k].M2L[n]*tree->at(j)[n].P2M;
				}
			}
		}
	}

	void writeActualMatrix(Mat& mat) {
		// assuming nlevels = 2
		int j = nLevels;
		mat = Mat(N,N);
		std::vector<int> rowIndex(nBoxesPerLevel[j]);
		std::vector<int> colIndex(nBoxesPerLevel[j]);
		for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
			rowIndex[0] = 0;
			colIndex[0] = 0;
			for (int k=1; k<nBoxesPerLevel[j]; ++k) {
				rowIndex[k] = rowIndex[k-1] + tree->at(j)[k-1].chargeLocations.size();
				colIndex[k] = colIndex[k-1] + tree->at(j)[k-1].chargeLocations.size();
			}
		}

		// Neighbors
		{
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				int row = rowIndex[k];
				for (size_t n = 0; n < nBoxesPerLevel[j]; n++) {
					int col = colIndex[n];
					mat.block(row, col, tree->at(j)[k].chargeLocations.size(), tree->at(j)[n].chargeLocations.size()) = K->getMatrix(tree->at(j)[k].chargeLocations, tree->at(j)[n].chargeLocations);
				}
			}
		}
	}

	void findErrorInHMatrix() {
		Mat actualMatrix, Hmatrix;
		writeActualMatrix(actualMatrix);
		writeHmatrix(Hmatrix);
		double err = (Hmatrix-actualMatrix).norm()/actualMatrix.norm();
		std::cout << "Err in H matrix: " << err << std::endl;
		return;
	}

	void getNodes_incoming_box(int j, int k) {
  // void getNodes_incoming_box(int j, int k, int& n_rows, int& n_cols, int& ComputedRank) {
		int n_rows, n_cols, ComputedRank;
    std::vector<int> boxA_Nodes;
		getParticlesFromChildren_incoming_row(j, k, boxA_Nodes);
		std::vector<int> IL_Nodes;//indices
		for(int in=0; in<tree->at(j)[k].interactionList.size(); ++in) {
				int kIL = tree->at(j)[k].interactionList[in];
				std::vector<int> chargeLocations;
				getParticlesFromChildren_incoming_col(j, kIL, chargeLocations);
				IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
		}
		n_rows = boxA_Nodes.size();
		n_cols = IL_Nodes.size();

		if (n_rows > 0 && n_cols > 0) {
			std::vector<int> row_indices, col_indices;
			row_indices = boxA_Nodes;//object of base class G_LowRank
			col_indices = IL_Nodes;
			std::vector<int> row_bases, col_bases;
			Mat Ar;
			LowRank<kerneltype>* LR		=	new LowRank<kerneltype>(K, pow(10,-TOL_POW), row_indices, col_indices);
			LR->ACA_only_nodes2(row_bases, col_bases, ComputedRank, tree->at(j)[k].Ac, Ar);
			// LR->rookPiv(row_bases, col_bases, ComputedRank, tree->at(j)[k].Ac, Ar);

			int minN = n_rows;
			if (n_rows > n_cols) {
				minN = n_cols;
			}
			if(ComputedRank > 0) {
				for (int r = 0; r < row_bases.size(); r++) {
					tree->at(j)[k].incoming_checkPoints.push_back(boxA_Nodes[row_bases[r]]);
				}
				for (int c = 0; c < col_bases.size(); c++) {
					tree->at(j)[k].incoming_chargePoints.push_back(IL_Nodes[col_bases[c]]);
				}
				// Mat D = K->getMatrix(tree->at(j)[k].incoming_checkPoints, tree->at(j)[k].incoming_chargePoints);
				// Eigen::PartialPivLU<Mat> D_T = Eigen::PartialPivLU<Mat>(D.transpose());
				// tree->at(j)[k].L2P = D_T.solve(tree->at(j)[k].Ac.transpose()).transpose();

				// Mat D = K->getMatrix(tree->at(j)[k].incoming_checkPoints, tree->at(j)[k].incoming_chargePoints);
				// Eigen::PartialPivLU<Mat> D_lu = Eigen::PartialPivLU<Mat>(D);
				// Mat mat = K->getMatrix(row_indices,col_indices);
				// Mat Err = mat - tree->at(j)[k].Ac*D_lu.solve(Ar);
				// double err = Err.norm()/mat.norm();
				// std::cout << "I j: " << j << "	k: " << k << "	rank: " << ComputedRank << "	err: " << err << std::endl;

			}
			// else {
			// 	tree->at(j)[k].L2P = Mat(tree->at(j)[k].Ac.rows(), tree->at(j)[k].incoming_checkPoints.size());
			// }
	    // std::cout << "j: " << j << "	k: "<< k << "	ComputedRank: " << ComputedRank << std::endl;
		}
		else if (n_cols == 0) {
			tree->at(j)[k].ILActive = false;
			tree->at(j)[k].incoming_checkPoints = boxA_Nodes;
			// if (j==8 && k==23) {
			// 	std::cout << "tree->at(j)[k].incoming_checkPoints.size(): " << tree->at(j)[k].incoming_checkPoints.size() << std::endl;
			// 	std::cout << "boxA_Nodes.size(): " << boxA_Nodes.size() << ", IL_Nodes.size(): " << IL_Nodes.size() << std::endl;
			// }
		}
		// if (j==5 && k==8) {
		// 	std::cout << "tree->at(j)[k].ILActive: " << tree->at(j)[k].ILActive << std::endl;
		// 	std::cout << "j: " << j << "	k: "<< k << "	ComputedRank: " << ComputedRank << "	n_rows: " << n_rows << "	n_cols: " << n_cols << std::endl;
			// 	std::cout << "boxA_Nodes.size(): " << boxA_Nodes.size() << std::endl;
		// 	std::cout << "tree->at(j+1)[2*k].incoming_checkPoints.size(): " << tree->at(j+1)[2*k].incoming_checkPoints.size() << std::endl;
		// 	std::cout << "tree->at(j+1)[2*k+1].incoming_checkPoints.size(): " << tree->at(j+1)[2*k+1].incoming_checkPoints.size() << std::endl;
		// }
	}

	void assemble_M2L() {
		// int numComputations = 0;
		// int M2LSize = 0;
		for (size_t j = startLevel; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				for(int in=0; in<tree->at(j)[k].interactionList.size(); ++in) {
						// numComputations++;
						int kIL = tree->at(j)[k].interactionList[in];
						tree->at(j)[k].M2L[kIL] = K->getMatrix(tree->at(j)[k].incoming_checkPoints, tree->at(j)[kIL].outgoing_chargePoints);
						// M2LSize += tree->at(j)[k].M2L[kIL].size();
				}
			}
		}
		// std::cout << "numComputations: " << numComputations << " M2LSize: " << M2LSize << std::endl;
	}

	void make_basis_unitary() {
		for (size_t j = startLevel; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				Eigen::FullPivHouseholderQR<Mat> qr(tree->at(j)[k].L2P.rows(), tree->at(j)[k].L2P.cols());
				tree->at(j)[k].U_qr = qr.compute(tree->at(j)[k].L2P);
				// Mat Q_u = qr.householderQ(); // Use it with ColPivHouseholderQR
				Mat Q_u = qr.matrixQ(); // Use it with FullPivHouseholderQR
				if(tree->at(j)[k].L2P.rows() > tree->at(j)[k].L2P.cols()) {
					Q_u = Q_u.block(0,0,tree->at(j)[k].L2P.rows(), tree->at(j)[k].L2P.cols());
				}
				tree->at(j)[k].L2P = Q_u;
			}
		}
	}

	void update_M2L_to_unitary_basis() {
		// #pragma omp parallel for
		for (size_t j = startLevel; j <= nLevels; j++) {
			// #pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				// #pragma omp parallel for
				for(int in=0; in<tree->at(j)[k].interactionList.size(); ++in) {
						int kIL = tree->at(j)[k].interactionList[in];
						// if (tree->at(j)[k].M2L[kIL].size() == 0) {
							Mat R_u = tree->at(j)[k].U_qr.matrixQR().triangularView<Eigen::Upper>();
							Mat R_v = tree->at(j)[kIL].U_qr.matrixQR().triangularView<Eigen::Upper>();
							if(tree->at(j)[k].L2P.rows() > tree->at(j)[k].L2P.cols()) {
								R_u = R_u.block(0,0,tree->at(j)[k].L2P.cols(), tree->at(j)[k].L2P.cols());
							}
							if(tree->at(j)[kIL].L2P.rows() > tree->at(j)[kIL].L2P.cols()) {
								R_v = R_v.block(0,0,tree->at(j)[kIL].L2P.cols(), tree->at(j)[kIL].L2P.cols());
							}
							tree->at(j)[k].M2L[kIL] = R_u * tree->at(j)[k].M2L[kIL] * R_v.transpose();
				}
			}
		}
	}

	void getUnitaryBasis(Mat &R, Mat &Lq, Mat &Lr) {
		// finds R = Lq * Lr
		Eigen::FullPivHouseholderQR<Mat> L(R.rows(), R.cols());
		L.setThreshold(RRQR_threshold);
		L.compute(R);
		// Mat LqOld = L.householderQ(); //new tree[j][n].L2P // Use it with ColPivHouseholderQR
		Mat LqOld = L.matrixQ(); //new tree[j][n].L2P // Use it with FullPivHouseholderQR
		// Mat Lr = L.matrixQR().triangularView<Eigen::Upper>();
		Lq = LqOld.block(0,0,LqOld.rows(),L.rank());
		Lr = Lq.adjoint()*R;
	}

	void findM2L(int j, int k, int b, Mat M2L_old, Mat &M2L_k_b, Mat Rk, Mat Rn, Mat Rn_prime) {
		if (M2L_old.size() > 0) {
			M2L_k_b = Rk*M2L_old*Rn.adjoint() + Rn_prime.adjoint();
		}
	}

	// ColPivHouseholderQR with weights; P2P row basis is found using P2P instead of Rk_prime
	void compress_P2P_with_weights(int j, int k, int b, Mat old_L2P, Mat old_P2M, Mat P2PToBeCompressed, Mat &new_L2P, Mat &new_P2M_transpose, Mat &Rk, Mat &Rn, Mat &Rn_prime) {
		// comment the below line
		// P2PToBeCompressed = 0*P2PToBeCompressed;

		Mat A_col(old_L2P.rows(), old_L2P.cols()+P2PToBeCompressed.cols());
		if (k==3 || b==3) {
			std::cout << "j: " << j << "	k: " << k << "	b: " << b << std::endl;
			std::cout << "old_L2P: " << old_L2P.rows() << "," << old_L2P.cols() << std::endl;
			std::cout << "P2PToBeCompressed: " << P2PToBeCompressed.rows() << "," << P2PToBeCompressed.cols() << std::endl;
		}
		// std::cout << "tree->at(j)[k].L2P_weights: " << std::endl << tree->at(j)[k].L2P_weights << std::endl;
		Mat old_L2P_temp = old_L2P;
		Vec weights1 = tree->at(j)[k].L2P_weights.head(old_L2P.cols());
		old_L2P = old_L2P*weights1.asDiagonal();
		A_col << old_L2P, P2PToBeCompressed;

		Eigen::ColPivHouseholderQR<Mat> qr_A_col(A_col.rows(), A_col.cols());
		// Eigen::HouseholderQR<Mat> qr_A_col(A_col.rows(), A_col.cols());
		// qr_A_col.setThreshold(RRQR_threshold);
		qr_A_col.compute(A_col);
		Mat new_L2P_old = qr_A_col.householderQ() ; //new tree->at(j)[k].L2P // Use it with ColPivHouseholderQR
		// Mat new_L2P_old = qr_A_col.matrixQ() ; //new tree->at(j)[k].L2P // Use it with FullPivHouseholderQR
		int rank = qr_A_col.rank();
		// Mat R_col = qr_A_col.matrixQR().triangularView<Eigen::Upper>();
		Mat R_col = qr_A_col.matrixR().topLeftCorner(qr_A_col.rank(), A_col.cols()).triangularView<Eigen::Upper>();
		new_L2P = new_L2P_old.block(0,0,new_L2P_old.rows(),rank);
		// new_L2P = new_L2P_old;
		tree->at(j)[k].L2P_weights = R_col.diagonal();
		// for (size_t i = 0; i < tree->at(j)[k].L2P_weights.size(); i++) {
		// 	if (tree->at(j)[k].L2P_weights(i) < 0) {
		// 		new_L2P.col(i) = -new_L2P.col(i);
		// 		tree->at(j)[k].L2P_weights(i) = -tree->at(j)[k].L2P_weights(i);
		// 	}
		// }

		// if (rank-tree->at(j)[k].L2P_weights.size() != 0) {
		// 	std::cout << "rank: " << rank << "	tree->at(j)[k].L2P_weights.size(): " << tree->at(j)[k].L2P_weights.size() << std::endl;
		// }
		// std::cout << "1 done" << std::endl;
		// if (k==3) {
		// 	Eigen::ColPivHouseholderQR<Mat>::PermutationType Pmat(qr_A_col.colsPermutation());
		// 	ArrayXi Pmat_indices = Pmat.indices();
		// 	ArrayXi Pmat_keep = Pmat_indices.head(rank);
		// 	std::cout << "Pmat_keep: " << Pmat_keep << std::endl;
		// }
		// std::cout << "Pmat_keep[0]: " << Pmat_keep[0] << std::endl;
		// std::cout << "Pmat_keep.size(): " << Pmat_keep.size() << std::endl;
		// Vec weightsTemp = tree->at(j)[k].L2P_weights;
		// for (size_t i = 0; i < Pmat_keep.size(); i++) {
		// 	int index = Pmat_keep[i];
		// 	if (index < old_L2P.cols()) { // column belongs to old_L2P
		// 		tree->at(j)[k].L2P_weights[i] = weightsTemp[index];
		// 	}
		// 	else { // column belongs to fill-in
		// 		tree->at(j)[k].L2P_weights[i] = weightsTemp[index];
		// 	}
		// }
		// std::cout << "R_col: " << std::endl << R_col << std::endl;
		// for (int i = 0; i < std::min(R_col.rows(),R_col.cols()); i++) {
		// 	// std::cout << "R_col(i,i): " << R_col(i,i) << std::endl;
		// 	if (std::abs(R_col(i,i)) > RRQR_threshold) {
		// 	}
		// 	else {
		// 		// std::cout << "i: " << i << std::endl;
		// 		rank = i+1;
		// 		break;
		// 	}
		// }

		// int P2P_rank1 = P2PToBeCompressed.cols()*0.00;
		// int rank1 = std::min(old_L2P.rows(), old_L2P.cols()+P2P_rank1);
		// new_L2P = new_L2P_old.block(0,0,new_L2P_old.rows(),rank1);

		if (k==3 || b==3) {
			std::cout << "new_L2P: " << new_L2P.rows() << ", " << new_L2P.cols() << std::endl;
		}

		Rk = new_L2P.adjoint()*old_L2P;
		Mat Rk_prime = new_L2P.adjoint()*P2PToBeCompressed;

		Mat Lij(Rk.rows(), Rk.cols()+Rk_prime.cols());
		Lij << Rk, Rk_prime;
		// std::cout << "A_col: " << std::endl << A_col << std::endl;
		// std::cout << "new_L2P: " << std::endl << new_L2P << std::endl;
		// std::cout << "Lij: " << std::endl << Lij << std::endl;
		// std::cout << "weights1: " << std::endl << weights1 << std::endl;
		Mat errMat1 = A_col - new_L2P*Lij;
		double err1 = errMat1.norm()/A_col.norm();
		if (err1 > 1e-10) {
			std::cout << "Err in QR1: " << err1 << std::endl;
			// exit(0);
		}

		Mat A_row(old_P2M.cols(), old_P2M.rows()+P2PToBeCompressed.rows());
		Vec weights2 = tree->at(j)[b].P2M_weights_adjoint.adjoint().head(old_P2M.rows());
		Mat old_P2M_temp = old_P2M;
		old_P2M = weights2.asDiagonal() * old_P2M;

		A_row << old_P2M.adjoint(), P2PToBeCompressed.adjoint();

		Eigen::HouseholderQR<Mat> qr_A_row(A_row.rows(), A_row.cols());
		qr_A_row.compute(A_row);
		Mat new_P2M_transpose_old = qr_A_row.householderQ(); //new tree[j][n].L2P // Use it with ColPivHouseholderQR
		Mat R_row = qr_A_row.matrixQR();
		new_P2M_transpose = new_P2M_transpose_old;
		int rank2 = std::min(A_row.rows(), A_row.cols());

		// Eigen::ColPivHouseholderQR<Mat> qr_A_row(A_row.rows(), A_row.cols());
		// qr_A_row.setThreshold(RRQR_threshold);
		// qr_A_row.compute(A_row);
		// Mat new_P2M_transpose_old = qr_A_row.householderQ(); //new tree[j][n].L2P // Use it with ColPivHouseholderQR
		// int rank2 = qr_A_row.rank();
		// Mat R_row = qr_A_row.matrixR().topLeftCorner(qr_A_row.rank(), A_row.cols()).triangularView<Eigen::Upper>();

		new_P2M_transpose = new_P2M_transpose_old.block(0,0,new_P2M_transpose_old.rows(),rank2);
		tree->at(j)[b].P2M_weights_adjoint = R_row.diagonal();
		// for (size_t i = 0; i < tree->at(j)[b].P2M_weights_adjoint.size(); i++) {
		// 	if (tree->at(j)[b].P2M_weights_adjoint(i) < 0) {
		// 		new_P2M_transpose.col(i) = -new_P2M_transpose.col(i);
		// 		tree->at(j)[b].P2M_weights_adjoint(i) = -tree->at(j)[b].P2M_weights_adjoint(i);
		// 	}
		// }
		// if (rank2-tree->at(j)[b].P2M_weights_adjoint.size() != 0) {
		// 	std::cout << "rank2: " << rank2 << "	tree->at(j)[b].P2M_weights_adjoint.size(): " << tree->at(j)[b].P2M_weights_adjoint.size() << std::endl;
		// }
		// std::cout << "2 done" << std::endl;

		// for (int i = 0; i < std::min(R_row.rows(),R_row.cols()); i++) {
		// 	if (std::abs(R_row(i,i)) > RRQR_threshold) {
		// 	}
		// 	else {
		// 		rank2 = i+1;
		// 		break;
		// 	}
		// }

		// int P2P_rank2 = Rk_prime.rows()*0.00;
		// int rank2 = std::min(old_P2M.cols(), old_P2M.rows() + P2P_rank2);
		// new_P2M_transpose = new_P2M_transpose_old.block(0,0,new_P2M_transpose_old.rows(),rank2);


		Rn = new_P2M_transpose.adjoint()*old_P2M.adjoint();
		Rn_prime = new_P2M_transpose.adjoint()*P2PToBeCompressed.adjoint();

		Mat Rij(Rn.rows(), Rn.cols()+Rn_prime.cols());
		Rij << Rn, Rn_prime;
		Mat errMat2 = A_row - new_P2M_transpose*Rij;
		double err2 = errMat2.norm()/A_row.norm();
		if (err2 > 1e-10) {
			std::cout << "Err in QR2: " << err2 << std::endl;
		}
		Rn_prime = Rn_prime * new_L2P;

		Rk = Rk * weights1.cwiseInverse().asDiagonal();
		Rn = Rn * weights2.adjoint().cwiseInverse().asDiagonal();
		// std::cout << "new_L2P: " << new_L2P.rows() << "," << new_L2P.cols() << std::endl;

	}

	// Q of P2P is considered and a fixed rank of it is considered for combined HouseholderQR; with QR for row basis
	void compress_P2P_orthogonal_fillin(int j, int k, int b, Mat old_L2P, Mat old_P2M, Mat P2PToBeCompressed, Mat &new_L2P, Mat &new_P2M_transpose, Mat &Rk, Mat &Rn, Mat &Rn_prime, Mat& newM2L) {
		Eigen::ColPivHouseholderQR<Mat> qr_P2P(P2PToBeCompressed.rows(), P2PToBeCompressed.cols());
		// Eigen::HouseholderQR<Mat> qr_P2P(P2PToBeCompressed.rows(), P2PToBeCompressed.cols());
		qr_P2P.setThreshold(fillin_RRQR_threshold);
		qr_P2P.compute(P2PToBeCompressed);
		Mat Q_P2P_old = qr_P2P.householderQ();
		// int rank = std::min(P2PToBeCompressed.rows(), P2PToBeCompressed.cols());
		// int rank = std::max(qr_P2P.rank(),1);
		int rank = std::min(fillin_rank,int(qr_P2P.rank()));
		// std::cout << "j: " << j << "	qr_P2P.rank(): " << qr_P2P.rank() << std::endl;
		// int rank = qr_P2P.rank();
		Mat Q_P2P = Q_P2P_old.block(0,0,Q_P2P_old.rows(),rank);
		Mat R_P2P = Q_P2P.adjoint()*P2PToBeCompressed; //new tree[j][k].L2P
		// std::cout << "1" << std::endl;
		Mat A_col(old_L2P.rows(), old_L2P.cols()+Q_P2P.cols());
		A_col << old_L2P, Q_P2P;

		// Eigen::ColPivHouseholderQR<Mat> qr_A_col(A_col.rows(), A_col.cols());
		Eigen::HouseholderQR<Mat> qr_A_col(A_col.rows(), A_col.cols());
		// qr_A_col.setThreshold(RRQR_threshold);
		qr_A_col.compute(A_col);
		Mat new_L2P_old = qr_A_col.householderQ() ; //new tree->at(j)[k].L2P // Use it with ColPivHouseholderQR
		Mat R_col = qr_A_col.matrixQR().triangularView<Eigen::Upper>();
		Vec R_col_diagonal = R_col.diagonal();
		// std::cout << "R_col_diagonal: " << std::endl << R_col_diagonal << std::endl;
		// exit(0);
		// int rank_c1 = qr_A_col.rank();
		int rank_c1 = std::min(A_col.rows(), A_col.cols());
		new_L2P = new_L2P_old.block(0,0,new_L2P_old.rows(),rank_c1);
		Rk = new_L2P.adjoint()*old_L2P;
		Mat Rk_prime = new_L2P.adjoint()*Q_P2P;
		// std::cout << "3" << std::endl;

		Mat Rk_new = Rk_prime*R_P2P;
		Eigen::ColPivHouseholderQR<Mat> qr_P2P2(Rk_new.adjoint().rows(), Rk_new.adjoint().cols());
		// Eigen::HouseholderQR<Mat> qr_P2P2(Rk_new.adjoint().rows(), Rk_new.adjoint().cols());
		qr_P2P2.setThreshold(fillin_RRQR_threshold);
		qr_P2P2.compute(Rk_new.adjoint());
		Mat Q_Rk_old = qr_P2P2.householderQ();
		// int rank2 = qr_P2P2.rank();
		int rank2 = std::min(fillin_rank,int(qr_P2P2.rank()));
		// std::cout << "j: " << j << "	qr_P2P2.rank(): " << qr_P2P2.rank() << std::endl;
		// int rank2 = std::min(Rk_new.adjoint().rows(), Rk_new.adjoint().cols());
		Mat Q_Rk = Q_Rk_old.block(0,0,Q_Rk_old.rows(),rank2);
		Mat R_Rk = Q_Rk.adjoint()*Rk_new.adjoint(); //new tree[j][k].L2P
		// std::cout << "4" << std::endl;

		Mat A_row(old_P2M.cols(), old_P2M.rows()+Q_Rk.cols());
		A_row << old_P2M.adjoint(), Q_Rk;
		// Eigen::ColPivHouseholderQR<Mat> qr_A_row(A_row.rows(), A_row.cols());
		Eigen::HouseholderQR<Mat> qr_A_row(A_row.rows(), A_row.cols());
		// qr_A_row.setThreshold(RRQR_threshold);
		qr_A_row.compute(A_row);
		// std::cout << "5" << std::endl;
		// int rank_c2 = qr_A_row.rank();
		int rank_c2 = std::min(A_row.rows(), A_row.cols());

		Mat new_P2M_transpose_old = qr_A_row.householderQ(); //new tree[j][n].L2P // Use it with ColPivHouseholderQR
		new_P2M_transpose = new_P2M_transpose_old.block(0,0,new_P2M_transpose_old.rows(),rank_c2);
		Rn = new_P2M_transpose.adjoint()*old_P2M.adjoint();
		Rn_prime = new_P2M_transpose.adjoint()*Q_Rk;
		if (tree->at(j)[k].M2L[b].size() > 0) {
			newM2L = Rk*tree->at(j)[k].M2L[b]*Rn.adjoint() + R_Rk.adjoint()*Rn_prime.adjoint();
		}
		else {
			newM2L = R_Rk.adjoint()*Rn_prime.adjoint();
		}
	}

	// custom RRQR with no Q from P2P fill-in
	void compress_P2P_myRRQR(int j, int k, int b, Mat old_L2P, Mat old_P2M, Mat P2PToBeCompressed, Mat &new_L2P, Mat &new_P2M_transpose, Mat &Rk, Mat &Rn, Mat &Rn_prime, Mat& newM2L) {
		Mat A_col(old_L2P.rows(), old_L2P.cols()+P2PToBeCompressed.cols());
		A_col << old_L2P, P2PToBeCompressed;

		// if (j==3 && k==231 && b==239) {
			// std::cout << "1" << std::endl;
		// }
		RRQR *rrqr = new RRQR(A_col, TOL_POW, fill_in_TOL_POW, old_L2P.cols(), tree->at(j)[k].maxPivot_L2P,0);
		// std::cout << "start" << std::endl;
		// if (j==3 && k==231 && b==239) {
			// std::cout << "2" << std::endl;
		// }
		rrqr->colPivHouseholderQRWithPartition();
		double err1 = rrqr->errorInQR();
		std::cout << "err1: " << err1 << std::endl;
		// if (j==3 && k==231 && b==239) {
			// std::cout << "3" << std::endl;
		// }
		// std::cout << "done" << std::endl;
	  // Mat Q = rrqr->Q_trunc;
	  // Mat P = rrqr->P;
	  // Mat R = rrqr->R_trunc;
	  // std::cout << "Q: "  << std::endl << Q << std::endl;
	  // std::cout << "R: "  << std::endl << R << std::endl;
	  // std::cout << "P: "  << std::endl << P << std::endl;
	  // Mat Err1 = A_col*P - Q*R;
	  // double errQR = Err1.norm()/A_col.norm();
	  // Mat Err2 = Q.transpose()*Q - Mat::Identity(Q.cols(), Q.cols());
	  // double errQ = Err2.norm();
	  // std::cout << "errQR: " << errQR << std::endl;
	  // std::cout << "errQ: " << errQ << std::endl;
		// exit(0);
		tree->at(j)[k].maxPivot_L2P = std::max(tree->at(j)[k].maxPivot_L2P, rrqr->maxPivot);
		new_L2P = rrqr->Q_trunc;
		Rk = new_L2P.adjoint()*old_L2P;
		Mat Rk_prime = new_L2P.adjoint()*P2PToBeCompressed;


		// Eigen::ColPivHouseholderQR<Mat> qr_P2P2(Rk_prime.adjoint().rows(), Rk_prime.adjoint().cols());
		//
		// qr_P2P2.setThreshold(fillin_RRQR_threshold);
		// qr_P2P2.compute(Rk_prime.adjoint());
		// Mat Q_Rk = qr_P2P2.householderQ();
		//
		// Mat R_Rk = Q_Rk.adjoint()*Rk_prime.adjoint(); //new tree[j][k].L2P

		Mat A_row(old_P2M.cols(), old_P2M.rows()+Rk_prime.adjoint().cols());
		A_row << old_P2M.adjoint(), Rk_prime.adjoint();

		// if (j==3 && k==231 && b==239) {
			// std::cout << "4" << std::endl;
		// }
		// RRQR *rrqr2;
		// if (j==3 && k==231 && b==239) {
		// 	rrqr2 = new RRQR(A_row, TOL_POW, fill_in_TOL_POW, old_P2M.rows(), tree->at(j)[b].maxPivot_P2M,1);
		// }
		// else {
		// 	rrqr2 = new RRQR(A_row, TOL_POW, fill_in_TOL_POW, old_P2M.rows(), tree->at(j)[b].maxPivot_P2M,0);
		// }
		RRQR *rrqr2 = new RRQR(A_row, TOL_POW, fill_in_TOL_POW, old_P2M.rows(), tree->at(j)[b].maxPivot_P2M,0);
		// if (j==3 && k==231 && b==239) {
			// std::cout << "5" << std::endl;
		// }
		rrqr2->colPivHouseholderQRWithPartition();
		double err2 = rrqr2->errorInQR();
		std::cout << "err2: " << err2 << std::endl;
		// if (j==3 && k==231 && b==239) {
			// std::cout << "6" << std::endl;
		// }
		// std::cout << "done" << std::endl;

		tree->at(j)[b].maxPivot_P2M = std::max(tree->at(j)[b].maxPivot_P2M, rrqr2->maxPivot);
		new_P2M_transpose = rrqr2->Q_trunc;
		Rn = new_P2M_transpose.adjoint()*old_P2M.adjoint();
		Rn_prime = new_P2M_transpose.adjoint()*Rk_prime.adjoint();
		if (tree->at(j)[k].M2L[b].size() > 0) {
			newM2L = Rk*tree->at(j)[k].M2L[b]*Rn.adjoint() + Rn_prime.adjoint();
		}
		else {
			newM2L = Rn_prime.adjoint();
		}
	}

	void make_L2P_unitary(int j, int k) {
		Eigen::ColPivHouseholderQR<Mat> qr_L2P(tree->at(j)[k].L2P.rows(), tree->at(j)[k].L2P.cols());
		// Eigen::HouseholderQR<Mat> qr_L2P(tree->at(j)[k].L2P.rows(), tree->at(j)[k].L2P.cols());
		// qr_L2P.setThreshold(RRQR_threshold);
		qr_L2P.compute(tree->at(j)[k].L2P);
		Mat Q_L2P_old = qr_L2P.householderQ();
		// int rank = qr_L2P.rank();
		int rank = tree->at(j)[k].L2P.cols();//qr_L2P.rank();
		// Mat Q_L2P = Q_L2P_old.block(0,0,Q_L2P_old.rows(),qr_L2P.rank());
		Mat Q_L2P = Q_L2P_old.block(0,0,Q_L2P_old.rows(),rank);
		// std::cout << "QR done" << std::endl;
		// std::cout << "old rank: " << tree->at(j)[k].L2P.cols() << "	qr_L2P.rank(): " << qr_L2P.rank() << std::endl;
		// Mat R_L2P = Q_L2P.adjoint()*tree->at(j)[k].L2P;
		update_other_M2Ls_due_L2P_make_basis_unitary(j, k, Q_L2P);
		// std::cout << "update_other_M2Ls_due_L2P done" << std::endl;
		if(j>startLevel) {
			update_parent_L2P_old(j, k, Q_L2P);
		}
		// std::cout << "update_parent_L2P_old done" << std::endl;
		tree->at(j)[k].L2P = Q_L2P;
	}

	void make_P2M_unitary(int j, int k) {
		Eigen::ColPivHouseholderQR<Mat> qr_P2M(tree->at(j)[k].P2M.adjoint().rows(), tree->at(j)[k].P2M.adjoint().cols());
		// Eigen::HouseholderQR<Mat> qr_P2M(tree->at(j)[k].P2M.adjoint().rows(), tree->at(j)[k].P2M.adjoint().cols());
		// qr_P2M.setThreshold(RRQR_threshold);
		qr_P2M.compute(tree->at(j)[k].P2M.adjoint());
		Mat Q_P2M_old = qr_P2M.householderQ();
		// int rank = qr_P2M.rank();
		int rank = tree->at(j)[k].P2M.adjoint().cols();//qr_P2M.rank();
		Mat Q_P2M = Q_P2M_old.block(0,0,Q_P2M_old.rows(),rank);
		// Mat Q_P2M = Q_P2M_old;
		// Mat R_P2M = Q_L2P.adjoint()*tree->at(j)[k].L2P;
		// std::cout << "QR done" << std::endl;
		// std::cout << "old rank: " << tree->at(j)[k].P2M.adjoint().cols() << "	qr_P2M.rank(): " << qr_P2M.rank() << std::endl;
		update_other_M2Ls_due_P2M_make_basis_unitary(j, k, Q_P2M);
		// std::cout << "update_other_M2Ls_due_P2M done" << std::endl;
		if(j>startLevel) {
			update_parent_P2M_old(j, k, Q_P2M);
		}
		// std::cout << "update_parent_P2M_old done" << std::endl;
		tree->at(j)[k].P2M = Q_P2M.adjoint();
		// tree->at(j)[k].P2M_weights_adjoint = qr_P2M.matrixR().diagonal();
	}

	void make_basis_unitary2() {
		for (size_t j = 2; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				// std::cout << "j: " << j << "	k: " << k << std::endl;
				if (tree->at(j)[k].L2P.size() > 0) {
					make_L2P_unitary(j,k);
					// std::cout << "make_L2P_unitary done" << std::endl;
					make_P2M_unitary(j,k);
					// std::cout << "make_P2M_unitary done" << std::endl;
				}
			}
		}
	}

	// only old basis, fill-ins are projected with orthogonal old basis
	void compress_P2P_projection(int j, int k, int b, Mat old_L2P, Mat old_P2M, Mat P2PToBeCompressed, Mat &new_L2P, Mat &new_P2M_transpose, Mat &Rk, Mat &Rn, Mat &Rn_prime) {
		Mat Rk_prime = old_L2P.adjoint()*P2PToBeCompressed;
		Rk = Mat::Identity(old_L2P.cols(),old_L2P.cols());

		Rn = Mat::Identity(old_P2M.rows(),old_P2M.rows());
		Rn_prime = old_P2M*Rk_prime.adjoint();
	}

	// only projections of fill-ins found
	void compress_P2P_Two_Ways_projection(int j, int k, int b) {
		if (tree->at(j)[b].M2L[k].size() == 0 || tree->at(j)[k].M2L[b].size() == 0) {
			return;
		}
		Mat new_L2P_k;
		Mat new_P2M_transpose_b;
		Mat new_L2P_k_modified;
		Mat new_P2M_transpose_b_modified;
		Mat M2L_k_b;
		Mat Rk;
		Mat Rn;
		Mat Rn_prime;

		compress_P2P_projection(j, k, b, tree->at(j)[k].L2P, tree->at(j)[b].P2M, tree->at(j)[k].P2P[b], new_L2P_k, new_P2M_transpose_b, Rk, Rn, Rn_prime);
		findM2L(j, k, b, tree->at(j)[k].M2L[b], tree->at(j)[k].M2L[b], Rk, Rn, Rn_prime);

		Mat new_L2P_b;
		Mat new_P2M_transpose_k;
		Mat new_L2P_b_modified;
		Mat new_P2M_transpose_k_modified;
		Mat M2L_b_k;
		compress_P2P_projection(j, b, k, tree->at(j)[b].L2P, tree->at(j)[k].P2M, tree->at(j)[b].P2P[k], new_L2P_b, new_P2M_transpose_k, Rk, Rn, Rn_prime);
		findM2L(j, b, k, tree->at(j)[b].M2L[k], tree->at(j)[b].M2L[k], Rk, Rn, Rn_prime);

		// tree->at(j)[k].M2L[b] = M2L_k_b.block(0,0,new_L2P_k_modified.cols(),new_P2M_transpose_b_modified.cols());
		// tree->at(j)[b].M2L[k] = M2L_b_k.block(0,0,new_L2P_b_modified.cols(),new_P2M_transpose_k_modified.cols());
	}

	// new basis for fill-ins found
	void compress_P2P_Two_Ways_orthogonal_fillin(int j, int k, int b) {
		// if (k==13 || k==48 || k==50 ||k==53 || b==13 || b==48 || b==50 || b==53) {
		// 	return;
		// }
		if (tree->at(j)[b].M2L[k].size() == 0 || tree->at(j)[k].M2L[b].size() == 0) {
			return;
		}
		Mat new_L2P_k;
		Mat new_P2M_transpose_b;
		Mat new_L2P_k_modified;
		Mat new_P2M_transpose_b_modified;
		Mat M2L_k_b;
		// old_L2P : tree->at(j)[k].L2P
		// old_P2M : tree->at(j)[b].P2M
		// old_M2L : tree->at(j)[k].M2L[b]
		// P2PToBeCompressed : tree->at(j)[k].P2P[b]
		Mat Rk;
		Mat Rn;
		Mat Rn_prime;
		// checks
		// Mat old_full_mat_kb = tree->at(j)[k].L2P*tree->at(j)[k].M2L[b]*tree->at(j)[b].P2M + tree->at(j)[k].P2P[b];
		// Mat old_full_mat_bk = tree->at(j)[b].L2P*tree->at(j)[b].M2L[k]*tree->at(j)[k].P2M + tree->at(j)[b].P2P[k];

		// std::cout << "start 1" << std::endl;
		compress_P2P_orthogonal_fillin(j, k, b, tree->at(j)[k].L2P, tree->at(j)[b].P2M, tree->at(j)[k].P2P[b], new_L2P_k, new_P2M_transpose_b, Rk, Rn, Rn_prime, M2L_k_b);
		// compress_P2P(j, k, b, tree->at(j)[k].L2P, tree->at(j)[b].P2M, tree->at(j)[k].P2P[b], new_L2P_k, new_P2M_transpose_b, Rk, Rn, Rn_prime);
		// std::cout << "compress_P2P done" << std::endl;
		// findM2L(j, k, b, tree->at(j)[k].M2L[b], M2L_k_b, Rk, Rn, Rn_prime);
		// compress_P2P(j, k, b, new_L2P_k, new_P2M_transpose_b, M2L_k_b);
		// std::cout << "findM2L done" << std::endl;

		Mat new_L2P_b;
		Mat new_P2M_transpose_k;
		Mat new_L2P_b_modified;
		Mat new_P2M_transpose_k_modified;
		Mat M2L_b_k;
		// std::cout << "start 2" << std::endl;
		compress_P2P_orthogonal_fillin(j, b, k, tree->at(j)[b].L2P, tree->at(j)[k].P2M, tree->at(j)[b].P2P[k], new_L2P_b, new_P2M_transpose_k, Rk, Rn, Rn_prime, M2L_b_k);
		// compress_P2P(j, b, k, tree->at(j)[b].L2P, tree->at(j)[k].P2M, tree->at(j)[b].P2P[k], new_L2P_b, new_P2M_transpose_k, Rk, Rn, Rn_prime);
		// std::cout << "compress_P2P done" << std::endl;
		// findM2L(j, b, k, tree->at(j)[b].M2L[k], M2L_b_k, Rk, Rn, Rn_prime);
		// std::cout << "findM2L done" << std::endl;
		// compress_P2P(j, b, k, new_L2P_b, new_P2M_transpose_k, M2L_b_k);

		if (new_L2P_k.cols() < new_P2M_transpose_k.cols()) {
			new_P2M_transpose_k_modified = new_P2M_transpose_k.block(0,0,new_P2M_transpose_k.rows(),new_L2P_k.cols());
			new_L2P_k_modified = new_L2P_k;
		}
		else {
			new_L2P_k_modified = new_L2P_k.block(0,0,new_L2P_k.rows(),new_P2M_transpose_k.cols());
			new_P2M_transpose_k_modified = new_P2M_transpose_k;
		}
		if (new_L2P_b.cols() < new_P2M_transpose_b.cols()) {
			new_P2M_transpose_b_modified = new_P2M_transpose_b.block(0,0,new_P2M_transpose_b.rows(),new_L2P_b.cols());
			new_L2P_b_modified = new_L2P_b;
		}
		else {
			new_L2P_b_modified = new_L2P_b.block(0,0,new_L2P_b.rows(),new_P2M_transpose_b.cols());
			new_P2M_transpose_b_modified = new_P2M_transpose_b;
		}
		////////////////////

		if (M2L_k_b.size() > 0) {
			tree->at(j)[k].M2L[b] = M2L_k_b.block(0,0,new_L2P_k_modified.cols(),new_P2M_transpose_b_modified.cols());
		}

		// std::vector<Mat> mats_kb_a1, mats_kb_b1;
		// update_other_M2Ls_at_k_old_check(j, k, b, mats_kb_a1);
		// update_other_M2Ls_at_b_old_check(j, k, b, mats_kb_b1);


		// std::cout << "M2L crop done" << std::endl;
		update_other_M2Ls_at_k_old(j, k, b, new_L2P_k_modified);
		// std::cout << "update_other_M2Ls_at_k_old done" << std::endl;
		update_other_M2Ls_at_b_old(j, k, b, new_P2M_transpose_b_modified);
		// std::cout << "update_other_M2Ls_at_b_old done" << std::endl;
		if(j>startLevel) {
			update_parent_L2P_old(j, k, new_L2P_k_modified);
			update_parent_P2M_old(j, b, new_P2M_transpose_b_modified);
		}
		// std::cout << "update_parent_P2M_old done" << std::endl;
		tree->at(j)[k].L2P = new_L2P_k_modified;
		tree->at(j)[b].P2M = new_P2M_transpose_b_modified.adjoint();

		tree->at(j)[k].NumLocals = tree->at(j)[k].L2P.cols();
		tree->at(j)[b].NumMultipoles = tree->at(j)[b].P2M.rows();
		//////////////////////////////////////////////////////////////////////
		if (M2L_b_k.size() > 0) {
			tree->at(j)[b].M2L[k] = M2L_b_k.block(0,0,new_L2P_b_modified.cols(),new_P2M_transpose_k_modified.cols());
		}


		// std::vector<Mat> mats_bk_a1, mats_bk_b1;
		// update_other_M2Ls_at_k_old_check(j, b, k, mats_bk_a1);
		// update_other_M2Ls_at_b_old_check(j, b, k, mats_bk_b1);

		update_other_M2Ls_at_k_old(j, b, k, new_L2P_b_modified);
		update_other_M2Ls_at_b_old(j, b, k, new_P2M_transpose_k_modified);
		// std::cout << "update_other_M2Ls_at_b_old done" << std::endl;
		if(j>startLevel) {
			update_parent_L2P_old(j, b, new_L2P_b_modified);
			update_parent_P2M_old(j, k, new_P2M_transpose_k_modified);
		}
		// std::cout << "update_parent_P2M_old done" << std::endl;
		tree->at(j)[b].L2P = new_L2P_b_modified;
		tree->at(j)[k].P2M = new_P2M_transpose_k_modified.adjoint();

		tree->at(j)[b].NumLocals = tree->at(j)[b].L2P.cols();
		tree->at(j)[k].NumMultipoles = tree->at(j)[k].P2M.rows();
		// tree->at(j)[k].P2P[b] = Mat::Zero(tree->at(j)[k].P2P[b].rows(), tree->at(j)[k].P2P[b].cols());
		// tree->at(j)[b].P2P[k] = Mat::Zero(tree->at(j)[b].P2P[k].rows(), tree->at(j)[b].P2P[k].cols());
		// std::cout << "DONE" << std::endl;


		// checks
		// Mat new_full_mat_kb = tree->at(j)[k].L2P*tree->at(j)[k].M2L[b]*tree->at(j)[b].P2M;
		// Mat ErrMat_kb = new_full_mat_kb - old_full_mat_kb;
		// double err_kb = ErrMat_kb.norm()/old_full_mat_kb.norm();
		// if (err_kb > 1e-10) {
		// 	std::cout << "err in P2P_kb: " << err_kb << std::endl;
		// }
		//
		// Mat new_full_mat_bk = tree->at(j)[b].L2P*tree->at(j)[b].M2L[k]*tree->at(j)[k].P2M;
		// Mat ErrMat_bk = new_full_mat_bk - old_full_mat_bk;
		// double err_bk = ErrMat_bk.norm()/old_full_mat_bk.norm();
		// if (err_bk > 1e-10) {
		// 	std::cout << "err in P2P_bk: " << err_bk << std::endl;
		// }
		//
		// std::vector<Mat> mats_kb_a2, mats_kb_b2;
		// update_other_M2Ls_at_k_old_check(j, k, b, mats_kb_a2);
		// update_other_M2Ls_at_b_old_check(j, k, b, mats_kb_b2);
		// std::vector<double> errVec_kb_a, errVec_kb_b;
		// std::vector<Mat> err_kb_a, err_kb_b;
		// for (size_t i = 0; i < mats_kb_a2.size(); i++) {
		// 	err_kb_a.push_back(mats_kb_a2[i] - mats_kb_a1[i]);
		// 	errVec_kb_a.push_back(err_kb_a[i].norm()/mats_kb_a1[i].norm());
		// 	if (err_kb_a[i].norm() > 0 && errVec_kb_a[i] > 1e-10) {
		// 		std::cout << "i: " << i << "	err in other M2L_kb_a: " << errVec_kb_a[i] << std::endl;
		// 	}
		// }
		//
		// for (size_t i = 0; i < mats_kb_b2.size(); i++) {
		// 	err_kb_b.push_back(mats_kb_b2[i] - mats_kb_b1[i]);
		// 	errVec_kb_b.push_back(err_kb_b[i].norm()/mats_kb_b1[i].norm());
		// 	if (err_kb_b[i].norm() > 0 && errVec_kb_b[i] > 1e-10) {
		// 		std::cout << "i: " << i << "	err in other M2L_kb_b: " << errVec_kb_b[i] << std::endl;
		// 	}
		// }
		//
		//
		// std::vector<Mat> mats_bk_a2, mats_bk_b2;
		// update_other_M2Ls_at_k_old_check(j, b, k, mats_bk_a2);
		// update_other_M2Ls_at_b_old_check(j, b, k, mats_bk_b2);
		// std::vector<double> errVec_bk_a, errVec_bk_b;
		// std::vector<Mat> err_bk_a, err_bk_b;
		// for (size_t i = 0; i < mats_bk_a2.size(); i++) {
		// 	err_bk_a.push_back(mats_bk_a2[i] - mats_bk_a1[i]);
		// 	errVec_bk_a.push_back(err_bk_a[i].norm()/mats_bk_a1[i].norm());
		// 	if (err_bk_a[i].norm() > 0 && errVec_bk_a[i] > 1e-10) {
		// 		std::cout << "i: " << i << "	err in other M2L_bk_a: " << errVec_bk_a[i] << std::endl;
		// 	}
		// }
		//
		// for (size_t i = 0; i < mats_bk_b2.size(); i++) {
		// 	err_bk_b.push_back(mats_bk_b2[i] - mats_bk_b1[i]);
		// 	errVec_bk_b.push_back(err_bk_b[i].norm()/mats_bk_b1[i].norm());
		// 	if (err_bk_b[i].norm() > 0 && errVec_bk_b[i] > 1e-10) {
		// 		std::cout << "i: " << i << "	err in other M2L_bk_b: " << errVec_bk_b[i] << std::endl;
		// 	}
		// }

	}

	//with weights
	void compress_P2P_Two_Ways_weights(int j, int k, int b) {
		// if (k==13 || k==48 || k==50 ||k==53 || b==13 || b==48 || b==50 || b==53) {
		// 	return;
		// }
		if (tree->at(j)[b].M2L[k].size() == 0 || tree->at(j)[k].M2L[b].size() == 0) {
			return;
		}
		Mat new_L2P_k;
		Mat new_P2M_transpose_b;
		Mat new_L2P_k_modified;
		Mat new_P2M_transpose_b_modified;
		Mat M2L_k_b;
		// old_L2P : tree->at(j)[k].L2P
		// old_P2M : tree->at(j)[b].P2M
		// old_M2L : tree->at(j)[k].M2L[b]
		// P2PToBeCompressed : tree->at(j)[k].P2P[b]
		Mat Rk;
		Mat Rn;
		Mat Rn_prime;
		std::cout << "j: " << j << "	k: " << k << "	b: " << b << std::endl;
		// checks
		// tree->at(j)[k].P2P[b] = tree->at(j)[k].P2P[b]*0;
		// tree->at(j)[b].P2P[k] = tree->at(j)[b].P2P[k]*0;
		Mat old_full_mat_kb = tree->at(j)[k].L2P*tree->at(j)[k].M2L[b]*tree->at(j)[b].P2M + tree->at(j)[k].P2P[b];
		Mat old_full_mat_bk = tree->at(j)[b].L2P*tree->at(j)[b].M2L[k]*tree->at(j)[k].P2M + tree->at(j)[b].P2P[k];

		// std::cout << "start 1" << std::endl;
		compress_P2P_with_weights(j, k, b, tree->at(j)[k].L2P, tree->at(j)[b].P2M, tree->at(j)[k].P2P[b], new_L2P_k, new_P2M_transpose_b, Rk, Rn, Rn_prime);
		// compress_P2P(j, k, b, tree->at(j)[k].L2P, tree->at(j)[b].P2M, tree->at(j)[k].P2P[b], new_L2P_k, new_P2M_transpose_b, Rk, Rn, Rn_prime);
		// std::cout << "compress_P2P done" << std::endl;
		findM2L(j, k, b, tree->at(j)[k].M2L[b], M2L_k_b, Rk, Rn, Rn_prime);
		// compress_P2P(j, k, b, new_L2P_k, new_P2M_transpose_b, M2L_k_b);
		// std::cout << "findM2L done" << std::endl;

		Mat new_L2P_b;
		Mat new_P2M_transpose_k;
		Mat new_L2P_b_modified;
		Mat new_P2M_transpose_k_modified;
		Mat M2L_b_k;
		// std::cout << "start 2" << std::endl;
		compress_P2P_with_weights(j, b, k, tree->at(j)[b].L2P, tree->at(j)[k].P2M, tree->at(j)[b].P2P[k], new_L2P_b, new_P2M_transpose_k, Rk, Rn, Rn_prime);
		// compress_P2P(j, b, k, tree->at(j)[b].L2P, tree->at(j)[k].P2M, tree->at(j)[b].P2P[k], new_L2P_b, new_P2M_transpose_k, Rk, Rn, Rn_prime);
		// std::cout << "compress_P2P done" << std::endl;
		findM2L(j, b, k, tree->at(j)[b].M2L[k], M2L_b_k, Rk, Rn, Rn_prime);
		// std::cout << "findM2L done" << std::endl;
		// compress_P2P(j, b, k, new_L2P_b, new_P2M_transpose_k, M2L_b_k);

		if (new_L2P_k.cols() < new_P2M_transpose_k.cols()) {
			new_P2M_transpose_k_modified = new_P2M_transpose_k.block(0,0,new_P2M_transpose_k.rows(),new_L2P_k.cols());
			new_L2P_k_modified = new_L2P_k;
		}
		else {
			new_L2P_k_modified = new_L2P_k.block(0,0,new_L2P_k.rows(),new_P2M_transpose_k.cols());
			new_P2M_transpose_k_modified = new_P2M_transpose_k;
		}
		if (new_L2P_b.cols() < new_P2M_transpose_b.cols()) {
			new_P2M_transpose_b_modified = new_P2M_transpose_b.block(0,0,new_P2M_transpose_b.rows(),new_L2P_b.cols());
			new_L2P_b_modified = new_L2P_b;
		}
		else {
			new_L2P_b_modified = new_L2P_b.block(0,0,new_L2P_b.rows(),new_P2M_transpose_b.cols());
			new_P2M_transpose_b_modified = new_P2M_transpose_b;
		}
		////////////////////

		if (M2L_k_b.size() > 0) {
			tree->at(j)[k].M2L[b] = M2L_k_b.block(0,0,new_L2P_k_modified.cols(),new_P2M_transpose_b_modified.cols());
		}

		std::vector<Mat> mats_kb_a1, mats_kb_b1;
		update_other_M2Ls_at_k_old_check(j, k, b, mats_kb_a1);
		update_other_M2Ls_at_b_old_check(j, k, b, mats_kb_b1);


		// std::cout << "M2L crop done" << std::endl;
		update_other_M2Ls_at_k_old(j, k, b, new_L2P_k_modified);
		// std::cout << "update_other_M2Ls_at_k_old done" << std::endl;
		update_other_M2Ls_at_b_old(j, k, b, new_P2M_transpose_b_modified);
		// std::cout << "update_other_M2Ls_at_b_old done" << std::endl;
		if(j>startLevel) {
			update_parent_L2P_old(j, k, new_L2P_k_modified);
			update_parent_P2M_old(j, b, new_P2M_transpose_b_modified);
		}
		// std::cout << "update_parent_P2M_old done" << std::endl;
		tree->at(j)[k].L2P = new_L2P_k_modified;
		tree->at(j)[b].P2M = new_P2M_transpose_b_modified.adjoint();

		tree->at(j)[k].NumLocals = tree->at(j)[k].L2P.cols();
		tree->at(j)[b].NumMultipoles = tree->at(j)[b].P2M.rows();
		//////////////////////////////////////////////////////////////////////
		if (M2L_b_k.size() > 0) {
			tree->at(j)[b].M2L[k] = M2L_b_k.block(0,0,new_L2P_b_modified.cols(),new_P2M_transpose_k_modified.cols());
		}


		std::vector<Mat> mats_bk_a1, mats_bk_b1;
		update_other_M2Ls_at_k_old_check(j, b, k, mats_bk_a1);
		update_other_M2Ls_at_b_old_check(j, b, k, mats_bk_b1);

		update_other_M2Ls_at_k_old(j, b, k, new_L2P_b_modified);
		update_other_M2Ls_at_b_old(j, b, k, new_P2M_transpose_k_modified);
		// std::cout << "update_other_M2Ls_at_b_old done" << std::endl;
		if(j>startLevel) {
			update_parent_L2P_old(j, b, new_L2P_b_modified);
			update_parent_P2M_old(j, k, new_P2M_transpose_k_modified);
		}
		// std::cout << "update_parent_P2M_old done" << std::endl;
		if (k==3) {
			std::cout << "new tree->at(j)[k].L2P after compression because of P2M: " << tree->at(j)[k].L2P.rows() << ", " << tree->at(j)[k].L2P.cols() << std::endl;
			std::cout << "new tree->at(j)[b].L2P after compression because of P2M: " << tree->at(j)[b].L2P.rows() << ", " << tree->at(j)[b].L2P.cols() << std::endl;
		}
		tree->at(j)[b].L2P = new_L2P_b_modified;
		tree->at(j)[k].P2M = new_P2M_transpose_k_modified.adjoint();

		tree->at(j)[b].NumLocals = tree->at(j)[b].L2P.cols();
		tree->at(j)[k].NumMultipoles = tree->at(j)[k].P2M.rows();
		// tree->at(j)[k].P2P[b] = Mat::Zero(tree->at(j)[k].P2P[b].rows(), tree->at(j)[k].P2P[b].cols());
		// tree->at(j)[b].P2P[k] = Mat::Zero(tree->at(j)[b].P2P[k].rows(), tree->at(j)[b].P2P[k].cols());
		// std::cout << "DONE" << std::endl;


		// checks
		Mat new_full_mat_kb = tree->at(j)[k].L2P*tree->at(j)[k].M2L[b]*tree->at(j)[b].P2M;
		Mat ErrMat_kb = new_full_mat_kb - old_full_mat_kb;
		double err_kb = ErrMat_kb.norm()/old_full_mat_kb.norm();
		// std::cout << "old_full_mat_kb: " << std::endl << old_full_mat_kb << std::endl;
		// std::cout << "new_full_mat_kb: " << std::endl << new_full_mat_kb << std::endl;
		// if (err_kb > 1e-10) {
			std::cout << "err in P2P_kb: " << err_kb << std::endl;
			// exit(0);
		// }

		Mat new_full_mat_bk = tree->at(j)[b].L2P*tree->at(j)[b].M2L[k]*tree->at(j)[k].P2M;
		Mat ErrMat_bk = new_full_mat_bk - old_full_mat_bk;
		double err_bk = ErrMat_bk.norm()/old_full_mat_bk.norm();
		// if (err_bk > 1e-10) {
			std::cout << "err in P2P_bk: " << err_bk << std::endl;
			// exit(0);
		// }
		// std::vector<Mat> mats_kb_a2, mats_kb_b2;
		// update_other_M2Ls_at_k_old_check(j, k, b, mats_kb_a2);
		// update_other_M2Ls_at_b_old_check(j, k, b, mats_kb_b2);
		// std::vector<double> errVec_kb_a, errVec_kb_b;
		// std::vector<Mat> err_kb_a, err_kb_b;
		// for (size_t i = 0; i < mats_kb_a2.size(); i++) {
		// 	err_kb_a.push_back(mats_kb_a2[i] - mats_kb_a1[i]);
		// 	errVec_kb_a.push_back(err_kb_a[i].norm()/mats_kb_a1[i].norm());
		// 	if (err_kb_a[i].norm() > 0 && errVec_kb_a[i] > 1e-10) {
		// 		std::cout << "i: " << i << "	err in other M2L_kb_a: " << errVec_kb_a[i] << std::endl;
		// 	}
		// }
		//
		// for (size_t i = 0; i < mats_kb_b2.size(); i++) {
		// 	err_kb_b.push_back(mats_kb_b2[i] - mats_kb_b1[i]);
		// 	errVec_kb_b.push_back(err_kb_b[i].norm()/mats_kb_b1[i].norm());
		// 	if (err_kb_b[i].norm() > 0 && errVec_kb_b[i] > 1e-10) {
		// 		std::cout << "i: " << i << "	err in other M2L_kb_b: " << errVec_kb_b[i] << std::endl;
		// 	}
		// }
		//
		//
		// std::vector<Mat> mats_bk_a2, mats_bk_b2;
		// update_other_M2Ls_at_k_old_check(j, b, k, mats_bk_a2);
		// update_other_M2Ls_at_b_old_check(j, b, k, mats_bk_b2);
		// std::vector<double> errVec_bk_a, errVec_bk_b;
		// std::vector<Mat> err_bk_a, err_bk_b;
		// for (size_t i = 0; i < mats_bk_a2.size(); i++) {
		// 	err_bk_a.push_back(mats_bk_a2[i] - mats_bk_a1[i]);
		// 	errVec_bk_a.push_back(err_bk_a[i].norm()/mats_bk_a1[i].norm());
		// 	if (err_bk_a[i].norm() > 0 && errVec_bk_a[i] > 1e-10) {
		// 		std::cout << "i: " << i << "	err in other M2L_bk_a: " << errVec_bk_a[i] << std::endl;
		// 	}
		// }
		//
		// for (size_t i = 0; i < mats_bk_b2.size(); i++) {
		// 	err_bk_b.push_back(mats_bk_b2[i] - mats_bk_b1[i]);
		// 	errVec_bk_b.push_back(err_bk_b[i].norm()/mats_bk_b1[i].norm());
		// 	if (err_bk_b[i].norm() > 0 && errVec_bk_b[i] > 1e-10) {
		// 		std::cout << "i: " << i << "	err in other M2L_bk_b: " << errVec_bk_b[i] << std::endl;
		// 	}
		// }

	}

	//custom RRQR
	void compress_P2P_Two_Ways_myRRQR(int j, int k, int b) {
		// if (k==13 || k==48 || k==50 ||k==53 || b==13 || b==48 || b==50 || b==53) {
		// 	return;
		// }
		if (tree->at(j)[b].M2L[k].size() == 0 || tree->at(j)[k].M2L[b].size() == 0) {
			return;
		}
		Mat new_L2P_k;
		Mat new_P2M_transpose_b;
		Mat new_L2P_k_modified;
		Mat new_P2M_transpose_b_modified;
		Mat M2L_k_b;
		// old_L2P : tree->at(j)[k].L2P
		// old_P2M : tree->at(j)[b].P2M
		// old_M2L : tree->at(j)[k].M2L[b]
		// P2PToBeCompressed : tree->at(j)[k].P2P[b]
		Mat Rk;
		Mat Rn;
		Mat Rn_prime;
		// checks
		// Mat old_full_mat_kb = tree->at(j)[k].L2P*tree->at(j)[k].M2L[b]*tree->at(j)[b].P2M + tree->at(j)[k].P2P[b];
		// Mat old_full_mat_bk = tree->at(j)[b].L2P*tree->at(j)[b].M2L[k]*tree->at(j)[k].P2M + tree->at(j)[b].P2P[k];

		// if (j==3 && k==231 && b==239) {
			// std::cout << "start 1" << std::endl;
		// }
		compress_P2P_myRRQR(j, k, b, tree->at(j)[k].L2P, tree->at(j)[b].P2M, tree->at(j)[k].P2P[b], new_L2P_k, new_P2M_transpose_b, Rk, Rn, Rn_prime, M2L_k_b);
		// compress_P2P(j, k, b, tree->at(j)[k].L2P, tree->at(j)[b].P2M, tree->at(j)[k].P2P[b], new_L2P_k, new_P2M_transpose_b, Rk, Rn, Rn_prime);
		// if (j==3 && k==231 && b==239) {
			// std::cout << "compress_P2P done" << std::endl;
		// }
		// findM2L(j, k, b, tree->at(j)[k].M2L[b], M2L_k_b, Rk, Rn, Rn_prime);
		// compress_P2P(j, k, b, new_L2P_k, new_P2M_transpose_b, M2L_k_b);
		// std::cout << "findM2L done" << std::endl;

		Mat new_L2P_b;
		Mat new_P2M_transpose_k;
		Mat new_L2P_b_modified;
		Mat new_P2M_transpose_k_modified;
		Mat M2L_b_k;
		// if (j==3 && k==231 && b==239) {
			// std::cout << "start 2" << std::endl;
		// }
		compress_P2P_myRRQR(j, b, k, tree->at(j)[b].L2P, tree->at(j)[k].P2M, tree->at(j)[b].P2P[k], new_L2P_b, new_P2M_transpose_k, Rk, Rn, Rn_prime, M2L_b_k);
		// compress_P2P(j, b, k, tree->at(j)[b].L2P, tree->at(j)[k].P2M, tree->at(j)[b].P2P[k], new_L2P_b, new_P2M_transpose_k, Rk, Rn, Rn_prime);
		// if (j==3 && k==231 && b==239) {
			// std::cout << "compress_P2P done" << std::endl;
		// }
		// findM2L(j, b, k, tree->at(j)[b].M2L[k], M2L_b_k, Rk, Rn, Rn_prime);
		// std::cout << "findM2L done" << std::endl;
		// compress_P2P(j, b, k, new_L2P_b, new_P2M_transpose_k, M2L_b_k);

		if (new_L2P_k.cols() < new_P2M_transpose_k.cols()) {
			new_P2M_transpose_k_modified = new_P2M_transpose_k.block(0,0,new_P2M_transpose_k.rows(),new_L2P_k.cols());
			new_L2P_k_modified = new_L2P_k;
		}
		else {
			new_L2P_k_modified = new_L2P_k.block(0,0,new_L2P_k.rows(),new_P2M_transpose_k.cols());
			new_P2M_transpose_k_modified = new_P2M_transpose_k;
		}
		if (new_L2P_b.cols() < new_P2M_transpose_b.cols()) {
			new_P2M_transpose_b_modified = new_P2M_transpose_b.block(0,0,new_P2M_transpose_b.rows(),new_L2P_b.cols());
			new_L2P_b_modified = new_L2P_b;
		}
		else {
			new_L2P_b_modified = new_L2P_b.block(0,0,new_L2P_b.rows(),new_P2M_transpose_b.cols());
			new_P2M_transpose_b_modified = new_P2M_transpose_b;
		}
		////////////////////

		// if (M2L_k_b.size() > 0) {
			tree->at(j)[k].M2L[b] = M2L_k_b.block(0,0,new_L2P_k_modified.cols(),new_P2M_transpose_b_modified.cols());
		// }

		// std::vector<Mat> mats_kb_a1, mats_kb_b1;
		// update_other_M2Ls_at_k_old_check(j, k, b, mats_kb_a1);
		// update_other_M2Ls_at_b_old_check(j, k, b, mats_kb_b1);


		// std::cout << "M2L crop done" << std::endl;
		update_other_M2Ls_at_k_old(j, k, b, new_L2P_k_modified);
		// if (j==3 && k==231 && b==239) {
		// 	std::cout << "update_other_M2Ls_at_k_old done" << std::endl;
		// }
		update_other_M2Ls_at_b_old(j, k, b, new_P2M_transpose_b_modified);
		// if (j==3 && k==231 && b==239) {
		// 	std::cout << "update_other_M2Ls_at_b_old done" << std::endl;
		// }
		if(j>startLevel) {
			update_parent_L2P_old(j, k, new_L2P_k_modified);
			update_parent_P2M_old(j, b, new_P2M_transpose_b_modified);
		}
		// if (j==3 && k==231 && b==239) {
		// 	std::cout << "update_parent_P2M_old done" << std::endl;
		// }
		tree->at(j)[k].L2P = new_L2P_k_modified;
		tree->at(j)[b].P2M = new_P2M_transpose_b_modified.adjoint();

		tree->at(j)[k].NumLocals = tree->at(j)[k].L2P.cols();
		tree->at(j)[b].NumMultipoles = tree->at(j)[b].P2M.rows();
		//////////////////////////////////////////////////////////////////////
		// if (M2L_b_k.size() > 0) {
			tree->at(j)[b].M2L[k] = M2L_b_k.block(0,0,new_L2P_b_modified.cols(),new_P2M_transpose_k_modified.cols());
		// }


		// std::vector<Mat> mats_bk_a1, mats_bk_b1;
		// update_other_M2Ls_at_k_old_check(j, b, k, mats_bk_a1);
		// update_other_M2Ls_at_b_old_check(j, b, k, mats_bk_b1);

		update_other_M2Ls_at_k_old(j, b, k, new_L2P_b_modified);
		update_other_M2Ls_at_b_old(j, b, k, new_P2M_transpose_k_modified);
		// if (j==3 && k==231 && b==239) {
		// 	std::cout << "update_other_M2Ls_at_b_old done" << std::endl;
		// }
		if(j>startLevel) {
			update_parent_L2P_old(j, b, new_L2P_b_modified);
			update_parent_P2M_old(j, k, new_P2M_transpose_k_modified);
		}
		// if (j==3 && k==231 && b==239) {
		// 	std::cout << "update_parent_P2M_old done" << std::endl;
		// }
		tree->at(j)[b].L2P = new_L2P_b_modified;
		tree->at(j)[k].P2M = new_P2M_transpose_k_modified.adjoint();

		tree->at(j)[b].NumLocals = tree->at(j)[b].L2P.cols();
		tree->at(j)[k].NumMultipoles = tree->at(j)[k].P2M.rows();
		// tree->at(j)[k].P2P[b] = Mat::Zero(tree->at(j)[k].P2P[b].rows(), tree->at(j)[k].P2P[b].cols());
		// tree->at(j)[b].P2P[k] = Mat::Zero(tree->at(j)[b].P2P[k].rows(), tree->at(j)[b].P2P[k].cols());
		// std::cout << "DONE" << std::endl;


		// checks
		// Mat new_full_mat_kb = tree->at(j)[k].L2P*tree->at(j)[k].M2L[b]*tree->at(j)[b].P2M;
		// Mat ErrMat_kb = new_full_mat_kb - old_full_mat_kb;
		// double err_kb = ErrMat_kb.norm()/old_full_mat_kb.norm();
		// if (err_kb > 1e-10) {
		// 	std::cout << "err in P2P_kb: " << err_kb << std::endl;
		// }
		//
		// Mat new_full_mat_bk = tree->at(j)[b].L2P*tree->at(j)[b].M2L[k]*tree->at(j)[k].P2M;
		// Mat ErrMat_bk = new_full_mat_bk - old_full_mat_bk;
		// double err_bk = ErrMat_bk.norm()/old_full_mat_bk.norm();
		// if (err_bk > 1e-10) {
		// 	std::cout << "err in P2P_bk: " << err_bk << std::endl;
		// }
		//
		// std::vector<Mat> mats_kb_a2, mats_kb_b2;
		// update_other_M2Ls_at_k_old_check(j, k, b, mats_kb_a2);
		// update_other_M2Ls_at_b_old_check(j, k, b, mats_kb_b2);
		// std::vector<double> errVec_kb_a, errVec_kb_b;
		// std::vector<Mat> err_kb_a, err_kb_b;
		// for (size_t i = 0; i < mats_kb_a2.size(); i++) {
		// 	err_kb_a.push_back(mats_kb_a2[i] - mats_kb_a1[i]);
		// 	errVec_kb_a.push_back(err_kb_a[i].norm()/mats_kb_a1[i].norm());
		// 	if (err_kb_a[i].norm() > 0 && errVec_kb_a[i] > 1e-10) {
		// 		std::cout << "i: " << i << "	err in other M2L_kb_a: " << errVec_kb_a[i] << std::endl;
		// 	}
		// }
		//
		// for (size_t i = 0; i < mats_kb_b2.size(); i++) {
		// 	err_kb_b.push_back(mats_kb_b2[i] - mats_kb_b1[i]);
		// 	errVec_kb_b.push_back(err_kb_b[i].norm()/mats_kb_b1[i].norm());
		// 	if (err_kb_b[i].norm() > 0 && errVec_kb_b[i] > 1e-10) {
		// 		std::cout << "i: " << i << "	err in other M2L_kb_b: " << errVec_kb_b[i] << std::endl;
		// 	}
		// }
		//
		//
		// std::vector<Mat> mats_bk_a2, mats_bk_b2;
		// update_other_M2Ls_at_k_old_check(j, b, k, mats_bk_a2);
		// update_other_M2Ls_at_b_old_check(j, b, k, mats_bk_b2);
		// std::vector<double> errVec_bk_a, errVec_bk_b;
		// std::vector<Mat> err_bk_a, err_bk_b;
		// for (size_t i = 0; i < mats_bk_a2.size(); i++) {
		// 	err_bk_a.push_back(mats_bk_a2[i] - mats_bk_a1[i]);
		// 	errVec_bk_a.push_back(err_bk_a[i].norm()/mats_bk_a1[i].norm());
		// 	if (err_bk_a[i].norm() > 0 && errVec_bk_a[i] > 1e-10) {
		// 		std::cout << "i: " << i << "	err in other M2L_bk_a: " << errVec_bk_a[i] << std::endl;
		// 	}
		// }
		//
		// for (size_t i = 0; i < mats_bk_b2.size(); i++) {
		// 	err_bk_b.push_back(mats_bk_b2[i] - mats_bk_b1[i]);
		// 	errVec_bk_b.push_back(err_bk_b[i].norm()/mats_bk_b1[i].norm());
		// 	if (err_bk_b[i].norm() > 0 && errVec_bk_b[i] > 1e-10) {
		// 		std::cout << "i: " << i << "	err in other M2L_bk_b: " << errVec_bk_b[i] << std::endl;
		// 	}
		// }

	}

	// only projection is considered where the old basis are orthogonal
	void compress_M2P_P2L_projection(int j, int k, int b) {
		if (tree->at(j)[b].M2L[k].size() == 0 || tree->at(j)[k].M2L[b].size() == 0) {
			return;
		}
		Mat new_L2P_k;
		Mat new_L2P_k_modified;
		Mat M2L_k_b;
		compress_M2P_projection(j, k, b, new_L2P_k, tree->at(j)[k].M2L[b]);
		Mat new_P2M_transpose_k;
		Mat new_P2M_transpose_k_modified;
		Mat M2L_b_k;
		compress_P2L_projection(j, b, k, new_P2M_transpose_k, tree->at(j)[b].M2L[k]);
	}

	// fill-ins are orthogonalised and then compressed
	void compress_M2P_P2L_orthogonal_fillin(int j, int k, int b) {
		// std::cout << "0" << std::endl;
		// Mat old_full_mat_kb = tree->at(j)[k].L2P*tree->at(j)[k].M2L[b] + tree->at(j)[k].M2P[b];
		// Mat old_full_mat_bk = tree->at(j)[b].M2L[k]*tree->at(j)[k].P2M + tree->at(j)[b].P2L[k];
		// std::cout << "old L2P: " << std::endl << tree->at(j)[k].L2P << std::endl;
		// std::cout << "old M2L: " << std::endl << tree->at(j)[k].M2L[b] << std::endl;
		// std::cout << "M2P: " << std::endl << tree->at(j)[k].M2P[b] << std::endl;
		if (tree->at(j)[b].M2L[k].size() == 0 || tree->at(j)[k].M2L[b].size() == 0) {
			return;
		}
		// std::cout << "tree->at(j)[k].M2L[b]: " << tree->at(j)[k].M2L[b].rows() << "," << tree->at(j)[k].M2L[b].cols() << std::endl;
		Mat new_L2P_k;
		Mat new_L2P_k_modified;
		Mat M2L_k_b;
		// std::cout << "0a" << std::endl;
		compress_M2P_orthogonal_fillin(j, k, b, new_L2P_k, M2L_k_b);
		// std::cout << "1" << std::endl;
		Mat new_P2M_transpose_k;
		Mat new_P2M_transpose_k_modified;
		Mat M2L_b_k;

		// std::cout << "tree->at(j)[k].M2L[b]: " << tree->at(j)[k].M2L[b].rows() << "," << tree->at(j)[k].M2L[b].cols() << std::endl;

		compress_P2L_orthogonal_fillin(j, b, k, new_P2M_transpose_k, M2L_b_k);
		// std::cout << "tree->at(j)[k].M2L[b]: " << tree->at(j)[k].M2L[b].rows() << "," << tree->at(j)[k].M2L[b].cols() << std::endl;
		// std::cout << "2" << std::endl;
		// if (k==25 && b==12) {
		// 	new_P2M_transpose_k_modified = new_P2M_transpose_k;
		// 	new_L2P_k_modified = new_L2P_k;
		// 	int min_size = std::min(new_L2P_k_modified.cols(), new_P2M_transpose_k_modified.cols());
		// 	std::cout << "min_size: " << min_size << std::endl;
		// 	std::cout << "size: " << new_L2P_k_modified.cols() << ", " << new_P2M_transpose_k_modified.cols() << std::endl;
		// }
		// else {
		// std::cout << "3" << std::endl;

			if (new_L2P_k.cols() < new_P2M_transpose_k.cols()) {
				new_P2M_transpose_k_modified = new_P2M_transpose_k.block(0,0,new_P2M_transpose_k.rows(),new_L2P_k.cols());
				new_L2P_k_modified = new_L2P_k;
			}
			else {
				new_L2P_k_modified = new_L2P_k.block(0,0,new_L2P_k.rows(),new_P2M_transpose_k.cols());
				new_P2M_transpose_k_modified = new_P2M_transpose_k;
			}
			// std::cout << "4" << std::endl;

		// }

		// if (k==25 && b==12) {
		// 	int min_size = std::min(new_L2P_k.cols(), new_P2M_transpose_k.cols());
		// 	std::cout << "min_size: " << min_size << std::endl;
		// if (new_L2P_k.cols() != new_P2M_transpose_k.cols()) {
		// 	std::cout << "size: " << new_L2P_k.cols() << ", " << new_P2M_transpose_k.cols() << std::endl;
		// }
		// }
		// std::vector<Mat> mats_kb_a1, mats_bk_b1;
		// update_other_M2Ls_at_k_old_check(j, k, b, mats_kb_a1);
		// update_other_M2Ls_at_b_old_check(j, b, k, mats_bk_b1);
		// std::cout << "M2L_k_b: " << M2L_k_b.rows() << "," << M2L_k_b.cols() << std::endl;
		// std::cout << "new_L2P_k_modified: " << new_L2P_k_modified.rows() << "," << new_L2P_k_modified.cols() << std::endl;
		// std::cout << "tree->at(j)[k].M2L[b]: " << tree->at(j)[k].M2L[b].rows() << "," << tree->at(j)[k].M2L[b].cols() << std::endl;

		tree->at(j)[k].M2L[b] = M2L_k_b.block(0,0,new_L2P_k_modified.cols(),M2L_k_b.cols());
		// std::cout << "M2L update done" << std::endl;
		update_other_M2Ls_at_k_old(j, k, b, new_L2P_k_modified);
		// std::cout << "update_other_M2Ls_at_k_old done" << std::endl;
		if(j>startLevel) {
			update_parent_L2P_old(j, k, new_L2P_k_modified);
		}
		// std::cout << "update_parent_L2P_old done" << std::endl;
		tree->at(j)[k].L2P = new_L2P_k_modified;
		tree->at(j)[k].NumLocals = tree->at(j)[k].L2P.cols();
		///////////////////////////////////////////////////////
		tree->at(j)[b].M2L[k] = M2L_b_k.block(0,0,M2L_b_k.rows(),new_P2M_transpose_k_modified.cols());
		update_other_M2Ls_at_b_old(j, b, k, new_P2M_transpose_k_modified);
		// std::cout << "update_other_M2Ls_at_b_old done" << std::endl;
		if(j>startLevel) {
			update_parent_P2M_old(j, k, new_P2M_transpose_k_modified);
		}
		// std::cout << "update_parent_P2M_old done" << std::endl;
		tree->at(j)[k].P2M = new_P2M_transpose_k_modified.adjoint();
		tree->at(j)[k].NumMultipoles = tree->at(j)[k].P2M.rows();
		// tree->at(j)[k].M2P[b] = Mat::Zero(tree->at(j)[k].M2P[b].rows(), tree->at(j)[k].M2P[b].cols());
		// tree->at(j)[b].P2L[k] = Mat::Zero(tree->at(j)[b].P2L[k].rows(), tree->at(j)[b].P2L[k].cols());
		// // checks
		/*
		Mat new_full_mat_kb = tree->at(j)[k].L2P*tree->at(j)[k].M2L[b];
		Mat ErrMat_kb = new_full_mat_kb - old_full_mat_kb;
		double err_kb = ErrMat_kb.norm()/old_full_mat_kb.norm();
		// if (err_kb > 1e-5) {
			// std::cout << "k: " << k << "	b: " << b << "	err in self M2L kb: " << err_kb << std::endl;
		// }

		Mat new_full_mat_bk = tree->at(j)[b].M2L[k]*tree->at(j)[k].P2M;
		Mat ErrMat_bk = new_full_mat_bk - old_full_mat_bk;
		double err_bk = ErrMat_bk.norm()/old_full_mat_bk.norm();
		// if (err_bk > 1e-5) {
			// std::cout << "k: " << k << "	b: " << b << "	err in self M2L bk: " << err_bk << std::endl;
		// }

		std::vector<Mat> mats_kb_a2, mats_kb_b2;
		update_other_M2Ls_at_k_old_check(j, k, b, mats_kb_a2);
		// update_other_M2Ls_at_b_old_check(j, k, b, mats_kb_b2);
		std::vector<double> errVec_kb_a, errVec_kb_b;
		std::vector<Mat> err_kb_a, err_kb_b;
		for (size_t i = 0; i < mats_kb_a2.size(); i++) {
			err_kb_a.push_back(mats_kb_a2[i] - mats_kb_a1[i]);
			errVec_kb_a.push_back(err_kb_a[i].norm()/mats_kb_a1[i].norm());
			if (err_kb_a[i].norm() > 0 && errVec_kb_a[i] > 1e-5) {
				// std::cout << "i: " << i << "	err in other M2L_kb_a: " << errVec_kb_a[i] << std::endl;
			}
		}

		std::vector<Mat> mats_bk_a2, mats_bk_b2;
		// update_other_M2Ls_at_k_old_check(j, b, k, mats_bk_a2);
		update_other_M2Ls_at_b_old_check(j, b, k, mats_bk_b2);
		std::vector<double> errVec_bk_a, errVec_bk_b;
		std::vector<Mat> err_bk_a, err_bk_b;
		for (size_t i = 0; i < mats_bk_b2.size(); i++) {
			err_bk_b.push_back(mats_bk_b2[i] - mats_bk_b1[i]);
			errVec_bk_b.push_back(err_bk_b[i].norm()/mats_bk_b1[i].norm());
			if (err_bk_b[i].norm() > 0 && errVec_bk_b[i] > 1e-5) {
				// std::cout << "i: " << i << "	err in other M2L_bk_b: " << errVec_bk_b[i] << std::endl;
			}
		}
		// std::cout << "3" << std::endl;
		*/
	}

	// custom RRQR
	void compress_M2P_P2L_myRRQR(int j, int k, int b) {
		// std::cout << "0" << std::endl;
		// Mat old_full_mat_kb = tree->at(j)[k].L2P*tree->at(j)[k].M2L[b] + tree->at(j)[k].M2P[b];
		// Mat old_full_mat_bk = tree->at(j)[b].M2L[k]*tree->at(j)[k].P2M + tree->at(j)[b].P2L[k];
		// std::cout << "old L2P: " << std::endl << tree->at(j)[k].L2P << std::endl;
		// std::cout << "old M2L: " << std::endl << tree->at(j)[k].M2L[b] << std::endl;
		// std::cout << "M2P: " << std::endl << tree->at(j)[k].M2P[b] << std::endl;
		if (tree->at(j)[b].M2L[k].size() == 0 || tree->at(j)[k].M2L[b].size() == 0) {
			return;
		}
		// std::cout << "tree->at(j)[k].M2L[b]: " << tree->at(j)[k].M2L[b].rows() << "," << tree->at(j)[k].M2L[b].cols() << std::endl;
		Mat new_L2P_k;
		Mat new_L2P_k_modified;
		Mat M2L_k_b;
		// if (j==6 && k==151527 &&	b==151359) {
		// 	std::cout << "0a" << std::endl;
		// }
		compress_M2P_myRRQR(j, k, b, new_L2P_k, M2L_k_b);
		// if (j==6 && k==151527 &&	b==151359) {
		// 	std::cout << "1" << std::endl;
		// }
		Mat new_P2M_transpose_k;
		Mat new_P2M_transpose_k_modified;
		Mat M2L_b_k;

		// std::cout << "tree->at(j)[k].M2L[b]: " << tree->at(j)[k].M2L[b].rows() << "," << tree->at(j)[k].M2L[b].cols() << std::endl;
		compress_P2L_myRRQR(j, b, k, new_P2M_transpose_k, M2L_b_k);
		// std::cout << "tree->at(j)[k].M2L[b]: " << tree->at(j)[k].M2L[b].rows() << "," << tree->at(j)[k].M2L[b].cols() << std::endl;
		// if (j==6 && k==151527 &&	b==151359) {
		// 	std::cout << "2" << std::endl;
		// }
		// if (k==25 && b==12) {
		// 	new_P2M_transpose_k_modified = new_P2M_transpose_k;
		// 	new_L2P_k_modified = new_L2P_k;
		// 	int min_size = std::min(new_L2P_k_modified.cols(), new_P2M_transpose_k_modified.cols());
		// 	std::cout << "min_size: " << min_size << std::endl;
		// 	std::cout << "size: " << new_L2P_k_modified.cols() << ", " << new_P2M_transpose_k_modified.cols() << std::endl;
		// }
		// else {
		// std::cout << "3" << std::endl;

			if (new_L2P_k.cols() < new_P2M_transpose_k.cols()) {
				new_P2M_transpose_k_modified = new_P2M_transpose_k.block(0,0,new_P2M_transpose_k.rows(),new_L2P_k.cols());
				new_L2P_k_modified = new_L2P_k;
			}
			else {
				new_L2P_k_modified = new_L2P_k.block(0,0,new_L2P_k.rows(),new_P2M_transpose_k.cols());
				new_P2M_transpose_k_modified = new_P2M_transpose_k;
			}
			// std::cout << "4" << std::endl;

		// }

		// if (k==25 && b==12) {
		// 	int min_size = std::min(new_L2P_k.cols(), new_P2M_transpose_k.cols());
		// 	std::cout << "min_size: " << min_size << std::endl;
		// if (new_L2P_k.cols() != new_P2M_transpose_k.cols()) {
		// 	std::cout << "size: " << new_L2P_k.cols() << ", " << new_P2M_transpose_k.cols() << std::endl;
		// }
		// }
		// std::vector<Mat> mats_kb_a1, mats_bk_b1;
		// update_other_M2Ls_at_k_old_check(j, k, b, mats_kb_a1);
		// update_other_M2Ls_at_b_old_check(j, b, k, mats_bk_b1);
		// std::cout << "M2L_k_b: " << M2L_k_b.rows() << "," << M2L_k_b.cols() << std::endl;
		// std::cout << "new_L2P_k_modified: " << new_L2P_k_modified.rows() << "," << new_L2P_k_modified.cols() << std::endl;
		// std::cout << "tree->at(j)[k].M2L[b]: " << tree->at(j)[k].M2L[b].rows() << "," << tree->at(j)[k].M2L[b].cols() << std::endl;

		tree->at(j)[k].M2L[b] = M2L_k_b.block(0,0,new_L2P_k_modified.cols(),M2L_k_b.cols());
		// if (j==6 && k==151527 &&	b==151359) {
		// 	std::cout << "M2L update done" << std::endl;
		// }
		update_other_M2Ls_at_k_old(j, k, b, new_L2P_k_modified);
		// if (j==6 && k==151527 &&	b==151359) {
		// 	std::cout << "update_other_M2Ls_at_k_old done" << std::endl;
		// }
		if(j>startLevel) {
			update_parent_L2P_old(j, k, new_L2P_k_modified);
		}
		// if (j==6 && k==151527 &&	b==151359) {
		// 	std::cout << "update_parent_L2P_old done" << std::endl;
		// }
		tree->at(j)[k].L2P = new_L2P_k_modified;
		tree->at(j)[k].NumLocals = tree->at(j)[k].L2P.cols();
		///////////////////////////////////////////////////////
		tree->at(j)[b].M2L[k] = M2L_b_k.block(0,0,M2L_b_k.rows(),new_P2M_transpose_k_modified.cols());
		update_other_M2Ls_at_b_old(j, b, k, new_P2M_transpose_k_modified);
		// if (j==6 && k==151527 &&	b==151359) {
		// 	std::cout << "update_other_M2Ls_at_b_old done" << std::endl;
		// }
		if(j>startLevel) {
			update_parent_P2M_old(j, k, new_P2M_transpose_k_modified);
		}
		// if (j==6 && k==151527 &&	b==151359) {
		// 	std::cout << "update_parent_P2M_old done" << std::endl;
		// }
		tree->at(j)[k].P2M = new_P2M_transpose_k_modified.adjoint();
		tree->at(j)[k].NumMultipoles = tree->at(j)[k].P2M.rows();
		// tree->at(j)[k].M2P[b] = Mat::Zero(tree->at(j)[k].M2P[b].rows(), tree->at(j)[k].M2P[b].cols());
		// tree->at(j)[b].P2L[k] = Mat::Zero(tree->at(j)[b].P2L[k].rows(), tree->at(j)[b].P2L[k].cols());
		// // checks
		/*
		Mat new_full_mat_kb = tree->at(j)[k].L2P*tree->at(j)[k].M2L[b];
		Mat ErrMat_kb = new_full_mat_kb - old_full_mat_kb;
		double err_kb = ErrMat_kb.norm()/old_full_mat_kb.norm();
		// if (err_kb > 1e-5) {
			// std::cout << "k: " << k << "	b: " << b << "	err in self M2L kb: " << err_kb << std::endl;
		// }

		Mat new_full_mat_bk = tree->at(j)[b].M2L[k]*tree->at(j)[k].P2M;
		Mat ErrMat_bk = new_full_mat_bk - old_full_mat_bk;
		double err_bk = ErrMat_bk.norm()/old_full_mat_bk.norm();
		// if (err_bk > 1e-5) {
			// std::cout << "k: " << k << "	b: " << b << "	err in self M2L bk: " << err_bk << std::endl;
		// }

		std::vector<Mat> mats_kb_a2, mats_kb_b2;
		update_other_M2Ls_at_k_old_check(j, k, b, mats_kb_a2);
		// update_other_M2Ls_at_b_old_check(j, k, b, mats_kb_b2);
		std::vector<double> errVec_kb_a, errVec_kb_b;
		std::vector<Mat> err_kb_a, err_kb_b;
		for (size_t i = 0; i < mats_kb_a2.size(); i++) {
			err_kb_a.push_back(mats_kb_a2[i] - mats_kb_a1[i]);
			errVec_kb_a.push_back(err_kb_a[i].norm()/mats_kb_a1[i].norm());
			if (err_kb_a[i].norm() > 0 && errVec_kb_a[i] > 1e-5) {
				// std::cout << "i: " << i << "	err in other M2L_kb_a: " << errVec_kb_a[i] << std::endl;
			}
		}

		std::vector<Mat> mats_bk_a2, mats_bk_b2;
		// update_other_M2Ls_at_k_old_check(j, b, k, mats_bk_a2);
		update_other_M2Ls_at_b_old_check(j, b, k, mats_bk_b2);
		std::vector<double> errVec_bk_a, errVec_bk_b;
		std::vector<Mat> err_bk_a, err_bk_b;
		for (size_t i = 0; i < mats_bk_b2.size(); i++) {
			err_bk_b.push_back(mats_bk_b2[i] - mats_bk_b1[i]);
			errVec_bk_b.push_back(err_bk_b[i].norm()/mats_bk_b1[i].norm());
			if (err_bk_b[i].norm() > 0 && errVec_bk_b[i] > 1e-5) {
				// std::cout << "i: " << i << "	err in other M2L_bk_b: " << errVec_bk_b[i] << std::endl;
			}
		}
		// std::cout << "3" << std::endl;
		*/
	}

	// with weights
	void compress_M2P_P2L_weights(int j, int k, int b) {
		// std::cout << "0" << std::endl;
		// Mat old_full_mat_kb = tree->at(j)[k].L2P*tree->at(j)[k].M2L[b] + tree->at(j)[k].M2P[b];
		// Mat old_full_mat_bk = tree->at(j)[b].M2L[k]*tree->at(j)[k].P2M + tree->at(j)[b].P2L[k];
		// std::cout << "old L2P: " << std::endl << tree->at(j)[k].L2P << std::endl;
		// std::cout << "old M2L: " << std::endl << tree->at(j)[k].M2L[b] << std::endl;
		// std::cout << "M2P: " << std::endl << tree->at(j)[k].M2P[b] << std::endl;
		if (tree->at(j)[b].M2L[k].size() == 0 || tree->at(j)[k].M2L[b].size() == 0) {
			return;
		}
		// std::cout << "j: " << j << "	k: " << k << "	b: " << b << std::endl;
		// std::cout << "tree->at(j)[k].M2P[b]: " << tree->at(j)[k].M2P[b].rows() << "," << tree->at(j)[k].M2P[b].cols() << std::endl;
		// std::cout << "tree->at(j)[b].P2L[k]: " << tree->at(j)[b].P2L[k].rows() << "," << tree->at(j)[b].P2L[k].cols() << std::endl;
		// std::cout << "tree->at(j)[k].L2P: " << tree->at(j)[k].L2P.rows() << "," << tree->at(j)[k].L2P.cols() << std::endl;
		Mat new_L2P_k;
		Mat new_L2P_k_modified;
		Mat M2L_k_b;
		// std::cout << "0a" << std::endl;
		compress_M2P_weights(j, k, b, new_L2P_k, M2L_k_b);
		// std::cout << "1" << std::endl;
		Mat new_P2M_transpose_k;
		Mat new_P2M_transpose_k_modified;
		Mat M2L_b_k;

		// std::cout << "tree->at(j)[k].M2L[b]: " << tree->at(j)[k].M2L[b].rows() << "," << tree->at(j)[k].M2L[b].cols() << std::endl;

		compress_P2L_weights(j, b, k, new_P2M_transpose_k, M2L_b_k);
		// std::cout << "tree->at(j)[k].M2L[b]: " << tree->at(j)[k].M2L[b].rows() << "," << tree->at(j)[k].M2L[b].cols() << std::endl;
		// std::cout << "2" << std::endl;
		// if (k==25 && b==12) {
		// 	new_P2M_transpose_k_modified = new_P2M_transpose_k;
		// 	new_L2P_k_modified = new_L2P_k;
		// 	int min_size = std::min(new_L2P_k_modified.cols(), new_P2M_transpose_k_modified.cols());
		// 	std::cout << "min_size: " << min_size << std::endl;
		// 	std::cout << "size: " << new_L2P_k_modified.cols() << ", " << new_P2M_transpose_k_modified.cols() << std::endl;
		// }
		// else {
		// std::cout << "3" << std::endl;

			if (new_L2P_k.cols() < new_P2M_transpose_k.cols()) {
				new_P2M_transpose_k_modified = new_P2M_transpose_k.block(0,0,new_P2M_transpose_k.rows(),new_L2P_k.cols());
				new_L2P_k_modified = new_L2P_k;
			}
			else {
				new_L2P_k_modified = new_L2P_k.block(0,0,new_L2P_k.rows(),new_P2M_transpose_k.cols());
				new_P2M_transpose_k_modified = new_P2M_transpose_k;
			}
			// std::cout << "4" << std::endl;

		// }

		// if (k==25 && b==12) {
		// 	int min_size = std::min(new_L2P_k.cols(), new_P2M_transpose_k.cols());
		// 	std::cout << "min_size: " << min_size << std::endl;
		// if (new_L2P_k.cols() != new_P2M_transpose_k.cols()) {
		// 	std::cout << "size: " << new_L2P_k.cols() << ", " << new_P2M_transpose_k.cols() << std::endl;
		// }
		// }
		// std::vector<Mat> mats_kb_a1, mats_bk_b1;
		// update_other_M2Ls_at_k_old_check(j, k, b, mats_kb_a1);
		// update_other_M2Ls_at_b_old_check(j, b, k, mats_bk_b1);
		// std::cout << "M2L_k_b: " << M2L_k_b.rows() << "," << M2L_k_b.cols() << std::endl;
		// std::cout << "new_L2P_k_modified: " << new_L2P_k_modified.rows() << "," << new_L2P_k_modified.cols() << std::endl;
		// std::cout << "tree->at(j)[k].M2L[b]: " << tree->at(j)[k].M2L[b].rows() << "," << tree->at(j)[k].M2L[b].cols() << std::endl;

		tree->at(j)[k].M2L[b] = M2L_k_b.block(0,0,new_L2P_k_modified.cols(),M2L_k_b.cols());
		// std::cout << "M2L update done" << std::endl;
		update_other_M2Ls_at_k_old(j, k, b, new_L2P_k_modified);
		// std::cout << "update_other_M2Ls_at_k_old done" << std::endl;
		if(j>startLevel) {
			update_parent_L2P_old(j, k, new_L2P_k_modified);
		}
		// std::cout << "update_parent_L2P_old done" << std::endl;
		tree->at(j)[k].L2P = new_L2P_k_modified;
		tree->at(j)[k].NumLocals = tree->at(j)[k].L2P.cols();
		///////////////////////////////////////////////////////
		tree->at(j)[b].M2L[k] = M2L_b_k.block(0,0,M2L_b_k.rows(),new_P2M_transpose_k_modified.cols());
		update_other_M2Ls_at_b_old(j, b, k, new_P2M_transpose_k_modified);
		// std::cout << "update_other_M2Ls_at_b_old done" << std::endl;
		if(j>startLevel) {
			update_parent_P2M_old(j, k, new_P2M_transpose_k_modified);
		}
		// std::cout << "update_parent_P2M_old done" << std::endl;
		tree->at(j)[k].P2M = new_P2M_transpose_k_modified.adjoint();
		tree->at(j)[k].NumMultipoles = tree->at(j)[k].P2M.rows();
		// tree->at(j)[k].M2P[b] = Mat::Zero(tree->at(j)[k].M2P[b].rows(), tree->at(j)[k].M2P[b].cols());
		// tree->at(j)[b].P2L[k] = Mat::Zero(tree->at(j)[b].P2L[k].rows(), tree->at(j)[b].P2L[k].cols());
		// // checks
		/*
		Mat new_full_mat_kb = tree->at(j)[k].L2P*tree->at(j)[k].M2L[b];
		Mat ErrMat_kb = new_full_mat_kb - old_full_mat_kb;
		double err_kb = ErrMat_kb.norm()/old_full_mat_kb.norm();
		// if (err_kb > 1e-5) {
			// std::cout << "k: " << k << "	b: " << b << "	err in self M2L kb: " << err_kb << std::endl;
		// }

		Mat new_full_mat_bk = tree->at(j)[b].M2L[k]*tree->at(j)[k].P2M;
		Mat ErrMat_bk = new_full_mat_bk - old_full_mat_bk;
		double err_bk = ErrMat_bk.norm()/old_full_mat_bk.norm();
		// if (err_bk > 1e-5) {
			// std::cout << "k: " << k << "	b: " << b << "	err in self M2L bk: " << err_bk << std::endl;
		// }

		std::vector<Mat> mats_kb_a2, mats_kb_b2;
		update_other_M2Ls_at_k_old_check(j, k, b, mats_kb_a2);
		// update_other_M2Ls_at_b_old_check(j, k, b, mats_kb_b2);
		std::vector<double> errVec_kb_a, errVec_kb_b;
		std::vector<Mat> err_kb_a, err_kb_b;
		for (size_t i = 0; i < mats_kb_a2.size(); i++) {
			err_kb_a.push_back(mats_kb_a2[i] - mats_kb_a1[i]);
			errVec_kb_a.push_back(err_kb_a[i].norm()/mats_kb_a1[i].norm());
			if (err_kb_a[i].norm() > 0 && errVec_kb_a[i] > 1e-5) {
				// std::cout << "i: " << i << "	err in other M2L_kb_a: " << errVec_kb_a[i] << std::endl;
			}
		}

		std::vector<Mat> mats_bk_a2, mats_bk_b2;
		// update_other_M2Ls_at_k_old_check(j, b, k, mats_bk_a2);
		update_other_M2Ls_at_b_old_check(j, b, k, mats_bk_b2);
		std::vector<double> errVec_bk_a, errVec_bk_b;
		std::vector<Mat> err_bk_a, err_bk_b;
		for (size_t i = 0; i < mats_bk_b2.size(); i++) {
			err_bk_b.push_back(mats_bk_b2[i] - mats_bk_b1[i]);
			errVec_bk_b.push_back(err_bk_b[i].norm()/mats_bk_b1[i].norm());
			if (err_bk_b[i].norm() > 0 && errVec_bk_b[i] > 1e-5) {
				// std::cout << "i: " << i << "	err in other M2L_bk_b: " << errVec_bk_b[i] << std::endl;
			}
		}
		// std::cout << "3" << std::endl;
		*/
	}

	// fill-ins are projected and the old basis are orthogonal
	void compress_M2P_projection(int j, int k, int b, Mat &new_L2P, Mat &M2L) {
		Mat Rk_prime = tree->at(j)[k].L2P.adjoint()*tree->at(j)[k].M2P[b];
		M2L = tree->at(j)[k].M2L[b] + Rk_prime;
	}

	// fill-ins are projected and the old basis are orthogonal
	void compress_P2L_projection(int j, int k, int b, Mat &new_P2M_transpose, Mat &M2L) {
		Mat Rn_prime = tree->at(j)[b].P2M*tree->at(j)[k].P2L[b].adjoint();
		M2L = tree->at(j)[k].M2L[b] + Rn_prime.adjoint();
	}

	// HouseholderQR
	void compress_M2P_weights(int j, int k, int b, Mat &new_L2P, Mat &M2L) {
		Vec weights1 = tree->at(j)[k].L2P_weights.head(tree->at(j)[k].L2P.cols());
		// std::cout << "weights1: " << weights1.size() << std::endl;

		Mat A_col(tree->at(j)[k].L2P.rows(), tree->at(j)[k].L2P.cols()+tree->at(j)[k].M2P[b].cols());
		A_col << tree->at(j)[k].L2P*weights1.asDiagonal(), tree->at(j)[k].M2P[b];

		// Eigen::HouseholderQR<Mat> qr_A_col(A_col.rows(), A_col.cols());
		// qr_A_col.compute(A_col);
		// Mat new_L2P_old = qr_A_col.householderQ() ; //new tree->at(j)[k].L2P // Use it with ColPivHouseholderQR
		// int rank = std::min(A_col.rows(), A_col.cols());
		// Mat R_col = qr_A_col.matrixQR().triangularView<Eigen::Upper>();

		Eigen::ColPivHouseholderQR<Mat> qr_A_col(A_col.rows(), A_col.cols());
		qr_A_col.setThreshold(RRQR_threshold);
		qr_A_col.compute(A_col);
		int rank = qr_A_col.rank();
		Mat new_L2P_old = qr_A_col.householderQ() ; //new tree->at(j)[k].L2P // Use it with ColPivHouseholderQR
		Mat R_col = qr_A_col.matrixR().topLeftCorner(qr_A_col.rank(), A_col.cols()).triangularView<Eigen::Upper>();

		new_L2P = new_L2P_old.block(0,0,new_L2P_old.rows(),rank);
		// std::cout << "R_col: " << R_col.rows() << "," << R_col.cols() << std::endl;
		// std::cout << "new_L2P: " << new_L2P.rows() << "," << new_L2P.cols() << std::endl;
		// std::cout << "R_col.diagonal().size(): " << R_col.diagonal().size() << std::endl;
		tree->at(j)[k].L2P_weights = R_col.diagonal().head(new_L2P.cols());
		// std::cout << "ewwww" << std::endl;

		Mat Rk_prime = new_L2P.adjoint()*tree->at(j)[k].M2P[b];
		Mat Rk = new_L2P.adjoint()*tree->at(j)[k].L2P;

		M2L = Rk*tree->at(j)[k].M2L[b] + Rk_prime;
	}

	// HouseholderQR
	void compress_P2L_weights(int j, int k, int b, Mat &new_P2M_transpose, Mat &M2L) {
		Mat A_row(tree->at(j)[b].P2M.cols(), tree->at(j)[b].P2M.rows()+tree->at(j)[k].P2L[b].rows());
		Vec weights2 = tree->at(j)[b].P2M_weights_adjoint.adjoint().head(tree->at(j)[b].P2M.rows());
		A_row << tree->at(j)[b].P2M.adjoint()*weights2.adjoint().asDiagonal(), tree->at(j)[k].P2L[b].adjoint();

		// Eigen::HouseholderQR<Mat> qr_A_row(A_row.rows(), A_row.cols());
		// qr_A_row.compute(A_row);
		// Mat new_P2M_transpose_old = qr_A_row.householderQ(); //new tree->at(j)[b].L2P // Use it with ColPivHouseholderQR
		// Mat R_row = qr_A_row.matrixQR().triangularView<Eigen::Upper>();
		// int rank = std::min(A_row.rows(), A_row.cols());

		Eigen::ColPivHouseholderQR<Mat> qr_A_row(A_row.rows(), A_row.cols());
		qr_A_row.setThreshold(RRQR_threshold);
		qr_A_row.compute(A_row);
		Mat new_P2M_transpose_old = qr_A_row.householderQ(); //new tree->at(j)[b].L2P // Use it with ColPivHouseholderQR
		Mat R_row = qr_A_row.matrixQR().topLeftCorner(qr_A_row.rank(), A_row.cols()).triangularView<Eigen::Upper>();
		int rank = qr_A_row.rank();
		new_P2M_transpose = new_P2M_transpose_old.block(0,0,new_P2M_transpose_old.rows(),rank);

		tree->at(j)[b].P2M_weights_adjoint = R_row.diagonal();

		Mat Rn = new_P2M_transpose.adjoint()*tree->at(j)[b].P2M.adjoint();
		Mat Rn_prime = new_P2M_transpose.adjoint()*tree->at(j)[k].P2L[b].adjoint();
		M2L = tree->at(j)[k].M2L[b]*Rn.adjoint() + Rn_prime.adjoint();
	}

	// a predefined rank for fill-in is used; ColPivHouseholderQR where Q of fill-in is considered
	void compress_M2P_orthogonal_fillin(int j, int k, int b, Mat &new_L2P, Mat &M2L) {
		double start, end;
		start = omp_get_wtime();
		Eigen::ColPivHouseholderQR<Mat> qr_M2P(tree->at(j)[k].M2P[b].rows(), tree->at(j)[k].M2P[b].cols());
		// Eigen::HouseholderQR<Mat> qr_M2P(tree->at(j)[k].M2P[b].rows(), tree->at(j)[k].M2P[b].cols());
		qr_M2P.setThreshold(fillin_RRQR_threshold);
		qr_M2P.compute(tree->at(j)[k].M2P[b]);
		Mat Q_M2P_old = qr_M2P.householderQ();
		// std::cout << "compress_M2P 1" << std::endl;

		// int rank = qr_M2P.rank();
		// int rank = std::min(Q_M2P_old.rows(),tree->at(j)[k].M2P[b].cols());
		int rank = std::min(fillin_rank, int(qr_M2P.rank()));
		// std::cout << "j: " << j << "	qr_M2P.rank(): " << qr_M2P.rank() << std::endl;
		Mat Q_M2P = Q_M2P_old.block(0,0,Q_M2P_old.rows(),rank);
		Mat R_M2P = Q_M2P.adjoint()*tree->at(j)[k].M2P[b]; //new tree[j][k].L2P
		end = omp_get_wtime();
		timeFillInQR += end - start;
		// std::cout << "compress_M2P 2" << std::endl;

		Mat A_col(tree->at(j)[k].L2P.rows(), tree->at(j)[k].L2P.cols()+Q_M2P.cols());
		A_col << tree->at(j)[k].L2P, Q_M2P;
		// Eigen::ColPivHouseholderQR<Mat> qr_A_col(A_col.rows(), A_col.cols());
		Eigen::HouseholderQR<Mat> qr_A_col(A_col.rows(), A_col.cols());
		// qr_A_col.setThreshold(RRQR_threshold);
		qr_A_col.compute(A_col);
		Mat new_L2P_old = qr_A_col.householderQ() ; //new tree[j][k].L2P
		// int rank_c = qr_A_col.rank();
		int rank_c = std::min(A_col.rows(), A_col.cols());
		new_L2P = new_L2P_old.block(0,0,new_L2P_old.rows(),rank_c);
		// std::cout << "compress_M2P 3" << std::endl;

		Mat Rk = new_L2P.adjoint()*tree->at(j)[k].L2P;
		Mat Rk_prime = new_L2P.adjoint()*Q_M2P;
			M2L = Rk*tree->at(j)[k].M2L[b] + Rk_prime*R_M2P;
	}

	// a predefined rank for fill-in is used; ColPivHouseholderQR where Q of fill-in is considered
	void compress_P2L_orthogonal_fillin(int j, int k, int b, Mat &new_P2M_transpose, Mat &M2L) {
		double start, end;
		start = omp_get_wtime();
		// Eigen::HouseholderQR<Mat> qr_P2L(tree->at(j)[k].P2L[b].adjoint().rows(), tree->at(j)[k].P2L[b].adjoint().cols());
		Eigen::ColPivHouseholderQR<Mat> qr_P2L(tree->at(j)[k].P2L[b].adjoint().rows(), tree->at(j)[k].P2L[b].adjoint().cols());
		qr_P2L.setThreshold(fillin_RRQR_threshold);
		qr_P2L.compute(tree->at(j)[k].P2L[b].adjoint());
		Mat Q_P2L_old = qr_P2L.householderQ();
		// std::cout << "Q_P2L_old: " << Q_P2L_old.rows() << "," << Q_P2L_old.cols() << ", " << qr_P2L.rank() << std::endl;
		// int rank = qr_P2L.rank();
		int rank = std::min(fillin_rank, int(qr_P2L.rank()));
		// std::cout << "j: " << j << "	qr_P2L.rank(): " << qr_P2L.rank() << std::endl;
		Mat Q_P2L = Q_P2L_old.block(0,0,Q_P2L_old.rows(),rank);
		// Mat Q_P2L = Q_P2L_old;
		Mat R_P2L = Q_P2L.adjoint()*tree->at(j)[k].P2L[b].adjoint(); //new tree[j][k].L2P
		end = omp_get_wtime();
		timeFillInQR += end - start;

		Mat A_row(tree->at(j)[b].P2M.cols(), tree->at(j)[b].P2M.rows()+Q_P2L.cols());
		A_row << tree->at(j)[b].P2M.adjoint(), Q_P2L;
		// Eigen::ColPivHouseholderQR<Mat> qr_A_row(A_row.rows(), A_row.cols());
		Eigen::HouseholderQR<Mat> qr_A_row(A_row.rows(), A_row.cols());
		// qr_A_row.setThreshold(RRQR_threshold);
		qr_A_row.compute(A_row);
		Mat new_P2M_transpose_old = qr_A_row.householderQ(); //new tree[j][b].L2P
		int rank_c = std::min(A_row.rows(), A_row.cols());
		// int rank_c = qr_A_row.rank();
		new_P2M_transpose = new_P2M_transpose_old.block(0,0,new_P2M_transpose_old.rows(),rank_c);

		Mat Rn = new_P2M_transpose.adjoint()*tree->at(j)[b].P2M.adjoint();
		Mat Rn_prime = new_P2M_transpose.adjoint()*Q_P2L;
			M2L = tree->at(j)[k].M2L[b]*Rn.adjoint() + R_P2L.adjoint()*Rn_prime.adjoint();
	}

	// custom rrqr
	void compress_M2P_myRRQR(int j, int k, int b, Mat &new_L2P, Mat &M2L) {
		Mat A_col(tree->at(j)[k].L2P.rows(), tree->at(j)[k].L2P.cols()+tree->at(j)[k].M2P[b].cols());
		A_col << tree->at(j)[k].L2P, tree->at(j)[k].M2P[b];

		RRQR *rrqr = new RRQR(A_col, TOL_POW, fill_in_TOL_POW, tree->at(j)[k].L2P.cols(), tree->at(j)[k].maxPivot_L2P,0);
		rrqr->colPivHouseholderQRWithPartition();
		double err1 = rrqr->errorInQR();
		std::cout << "err1: " << err1 << std::endl;

		new_L2P = rrqr->Q_trunc;
		tree->at(j)[k].maxPivot_L2P = std::max(tree->at(j)[k].maxPivot_L2P, rrqr->maxPivot);
		// Mat Q = rrqr->Q;
		// Mat P = rrqr->P;
		// Mat R = rrqr->R_trunc;
		//
		// Mat Err1 = A_col*P - Q*R;
		// double errQR = Err1.norm()/A_col.norm();
		// Mat Err2 = Q.transpose()*Q - Mat::Identity(Q.cols(), Q.cols());
		// double errQ = Err2.norm();

		Mat Rk = new_L2P.adjoint()*tree->at(j)[k].L2P;
		Mat Rk_prime = new_L2P.adjoint()*tree->at(j)[k].M2P[b];
		M2L = Rk*tree->at(j)[k].M2L[b] + Rk_prime;
	}

	// custom rrqr
	void compress_P2L_myRRQR(int j, int k, int b, Mat &new_P2M_transpose, Mat &M2L) {
		Mat A_row(tree->at(j)[b].P2M.cols(), tree->at(j)[b].P2M.rows()+tree->at(j)[k].P2L[b].adjoint().cols());
		A_row << tree->at(j)[b].P2M.adjoint(), tree->at(j)[k].P2L[b].adjoint();
		// if (j==6 && b==151527 && k==151359) {
		// 	std::cout << "A" << std::endl;
		// }
		RRQR *rrqr = new RRQR(A_row, TOL_POW, fill_in_TOL_POW, tree->at(j)[b].P2M.adjoint().cols(), tree->at(j)[b].maxPivot_P2M,0);

		// RRQR *rrqr;
		// if (j==6 && b==151527 && k==151359) {
		// 	rrqr = new RRQR(A_row, TOL_POW, fill_in_TOL_POW, tree->at(j)[b].P2M.adjoint().cols(), tree->at(j)[b].maxPivot_P2M,1);
		// }
		// else {
		// 	rrqr = new RRQR(A_row, TOL_POW, fill_in_TOL_POW, tree->at(j)[b].P2M.adjoint().cols(), tree->at(j)[b].maxPivot_P2M,0);
		// }
		// if (j==6 && b==151527 && k==151359) {
		// 	std::cout << "A_row: " << A_row.rows() << "," << A_row.cols() << std::endl;
		// 	std::cout << "tree->at(j)[b].P2M.adjoint(): " << tree->at(j)[b].P2M.adjoint().rows() << "," << tree->at(j)[b].P2M.adjoint().cols() << std::endl;
		// 	std::cout << "tree->at(j)[b].maxPivot_P2M: " <<  tree->at(j)[b].maxPivot_P2M << std::endl;
		// 	std::cout << "B" << std::endl;
		// }
		rrqr->colPivHouseholderQRWithPartition();
		double err1 = rrqr->errorInQR();
		std::cout << "err1: " << err1 << std::endl;
		// if (j==6 && b==151527 && k==151359) {
		// 	std::cout << "C" << std::endl;
		// }
		new_P2M_transpose = rrqr->Q_trunc;
		// if (j==6 && b==151527 && k==151359) {
		// 	std::cout << "D" << std::endl;
		// }
		tree->at(j)[b].maxPivot_P2M = std::max(tree->at(j)[b].maxPivot_P2M, rrqr->maxPivot);
		// Mat Q = rrqr->Q;
		// Mat P = rrqr->P;
		// Mat R = rrqr->R_trunc;
		//
		// Mat Err1 = A_col*P - Q*R;
		// double errQR = Err1.norm()/A_col.norm();
		// Mat Err2 = Q.transpose()*Q - Mat::Identity(Q.cols(), Q.cols());
		// double errQ = Err2.norm();

		Mat Rn = new_P2M_transpose.adjoint()*tree->at(j)[b].P2M.adjoint();
		// if (j==6 && b==151527 && k==151359) {
		// 	std::cout << "E" << std::endl;
		// }
		Mat Rn_prime = new_P2M_transpose.adjoint()*tree->at(j)[k].P2L[b].adjoint();
		// if (j==6 && b==151527 && k==151359) {
		// 	std::cout << "F" << std::endl;
		// }
		M2L = tree->at(j)[k].M2L[b]*Rn.adjoint() + Rn_prime.adjoint();
		// if (j==6 && b==151527 && k==151359) {
		// 	std::cout << "G" << std::endl;
		// }
	}

	void update_other_M2Ls_at_k(int j, int k, int b, Mat &L) {
		// #pragma omp parallel for
		// std::cout << "L: " << L.rows() << "," << L.cols() << std::endl;
		for(int in=0; in<tree->at(j)[k].interactionList.size(); ++in) {
				int kIL = tree->at(j)[k].interactionList[in];
				if(kIL == b) { //worked out in update_M2L_and_basis method
					continue;
				}
				// std::cout << "L: " << L.rows() << "," << L.cols() << std::endl;
				// std::cout << "tree->at(j)[k].M2L[kIL]: " << tree->at(j)[k].M2L[kIL].rows() << "," << tree->at(j)[k].M2L[kIL].cols() << std::endl;
				tree->at(j)[k].M2L[kIL] = L * tree->at(j)[k].M2L[kIL];
		}
		//k's Neighbors
		// #pragma omp parallel for
		for(int n=0; n<tree->at(j)[k].neighborList.size(); ++n) {
			int kN = tree->at(j)[k].neighborList[n];
			if((kN != k) && tree->at(j)[kN].Eliminated && tree->at(j)[k].Eliminated) {
				// if (k==0 && kN==10) {
				// 	std::cout << "update_other_M2Ls_at_k 4: tree->at(j)[k].M2L[kN]: " << tree->at(j)[k].M2L[kN].rows() << "," << tree->at(j)[k].M2L[kN].cols() << std::endl;
				// }
				// std::cout << "L: " << L.rows() << "," << L.cols() << std::endl;
				// std::cout << "tree->at(j)[k].M2L[kN]: " << tree->at(j)[k].M2L[kN].rows() << "," << tree->at(j)[k].M2L[kN].cols() << std::endl;
				tree->at(j)[k].M2L[kN] = L * tree->at(j)[k].M2L[kN];
				// tree->at(j)[kN].M2L[k] = tree->at(j)[k].M2L[kN].transpose(); // compressing P2P_{nn,k}
			}
		}
		//k's self
		if(tree->at(j)[k].Eliminated) {
			// std::cout << "tree->at(j)[k].M2L[k]: " << tree->at(j)[k].M2L[k].rows() << "," << tree->at(j)[k].M2L[k].cols() << std::endl;
			tree->at(j)[k].M2L[k] = L * tree->at(j)[k].M2L[k];
		}
	}

	void update_other_M2Ls_at_b(int j, int k, int b, Mat &R) {
		// #pragma omp parallel for
		for(int in=0; in<tree->at(j)[b].interactionList.size(); ++in) {
			int kIL = tree->at(j)[b].interactionList[in];
			// std::cout << "kIL: " << kIL << " b: " << b << std::endl;
			// std::cout << "tree->at(j)[kIL].M2L[b]: " << tree->at(j)[kIL].M2L[b].rows() << "," << tree->at(j)[kIL].M2L[b].cols() << std::endl;
			// std::cout << "R.adjoint(): " << R.adjoint().rows() << "," << R.adjoint().cols() << std::endl;
				if(kIL == k) { //worked out in update_M2L_and_basis method
					continue;
				}
				tree->at(j)[kIL].M2L[b] = tree->at(j)[kIL].M2L[b] * R.adjoint();
		}
		// std::cout << "IL done" << std::endl;
		//k's Neighbors
		// #pragma omp parallel for
		for(int n=0; n<tree->at(j)[b].neighborList.size(); ++n) {
			int kN = tree->at(j)[b].neighborList[n];
			// std::cout << "kN: " << kN << std::endl;
			if(kN != b && tree->at(j)[kN].Eliminated && tree->at(j)[b].Eliminated) {
				// if (kN==0 && b==10) {
				// 	std::cout << "update_other_M2Ls_at_b 4: tree->at(j)[kN].M2L[b]: " << tree->at(j)[kN].M2L[b].rows() << "," << tree->at(j)[kN].M2L[b].cols() << std::endl;
				// }
				tree->at(j)[kN].M2L[b] = tree->at(j)[kN].M2L[b] * R.adjoint();
			}
		}
		// std::cout << "neighbor done" << std::endl;
		//k's self
		if(tree->at(j)[b].Eliminated) {
			tree->at(j)[b].M2L[b] = tree->at(j)[b].M2L[b] * R.adjoint();
		}
		// std::cout << "self done" << std::endl;
	}

	void update_parent_L2P(int j, int k, Mat &L) {
		int kP = k/nChild;//k's parent
		int c_k = k%nChild; // k is c_k^{th} child of k_parent

		Mat U[nChild];
		int row_index = 0;
		for (size_t c = 0; c < nChild; c++) {
			// std::cout << "c: " << c << std::endl;
			// std::cout << "tree->at(j-1)[kP].L2P.size(): " << tree->at(j-1)[kP].L2P.rows() << ", " << tree->at(j-1)[kP].L2P.cols() << std::endl;
			// std::cout << "tree->at(j)[kP*2+c].NumLocals: " << tree->at(j)[kP*2+c].NumLocals << std::endl;
			U[c] = tree->at(j-1)[kP].L2P.block(row_index,0,tree->at(j)[kP*nChild+c].NumLocals,tree->at(j-1)[kP].L2P.cols());
			row_index += tree->at(j)[kP*nChild+c].NumLocals;
		}
		U[c_k] = L * U[c_k];
		int n_rows = 0;
		for (size_t c = 0; c < nChild; c++) {
			n_rows += U[c].rows();
		}
		tree->at(j-1)[kP].L2P = Mat(n_rows, U[0].cols());
		tree->at(j-1)[kP].L2P << U[0], U[1], U[2], U[3], U[4], U[5], U[6], U[7];
	}

	void update_parent_P2M(int j, int b, Mat &R) {
		int bP = b/nChild;//k's parent
		int c_b = b%nChild; // k is c_b^{th} child of k_parent

		Mat V[nChild];
		int col_index = 0;
		for (size_t c = 0; c < nChild; c++) {
			// std::cout << "c: " << c << std::endl;
			// std::cout << "tree->at(j-1)[bP].P2M.size(): " << tree->at(j-1)[bP].P2M.rows() << ", " << tree->at(j-1)[bP].P2M.cols() << std::endl;
			// std::cout << "tree->at(j)[bP*2+c].NumMultipoles: " << tree->at(j)[bP*2+c].NumMultipoles << std::endl;
			V[c] = tree->at(j-1)[bP].P2M.block(0,col_index,tree->at(j-1)[bP].P2M.rows(),tree->at(j)[bP*nChild+c].NumMultipoles);
			col_index += tree->at(j)[bP*nChild+c].NumMultipoles;
		}
		V[c_b] = V[c_b] * R.adjoint();
		int n_cols = 0;
		for (size_t c = 0; c < nChild; c++) {
			n_cols += V[c].cols();
		}
		tree->at(j-1)[bP].P2M = Mat(V[0].rows(), n_cols);
		tree->at(j-1)[bP].P2M << V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7];
	}

	void update_other_M2Ls_at_k(int j, int k, Mat &new_L2P) {
		// #pragma omp parallel for
		for (auto const& [key, val] : tree->at(j)[k].P2L) {
			tree->at(j)[k].P2L[key] = new_L2P.adjoint() * tree->at(j)[k].L2P * tree->at(j)[k].P2L[key];
		}
		for(int in=0; in<tree->at(j)[k].interactionList.size(); ++in) {
				int kIL = tree->at(j)[k].interactionList[in];
				// std::cout << "new_L2P.adjoint(): " << new_L2P.adjoint().rows() << "," << new_L2P.adjoint().cols() << std::endl;
				// std::cout << "tree->at(j)[k].L2P: " << tree->at(j)[k].L2P.rows() << "," << tree->at(j)[k].L2P.cols() << std::endl;
				// std::cout << "tree->at(j)[k].M2L[kIL]: " << tree->at(j)[k].M2L[kIL].rows() << "," << tree->at(j)[k].M2L[kIL].cols() << std::endl;
				tree->at(j)[k].M2L[kIL] = new_L2P.adjoint() * tree->at(j)[k].L2P * tree->at(j)[k].M2L[kIL];
		}
		// std::cout << "IL done" << std::endl;
		//k's Neighbors
		// #pragma omp parallel for
		for(int n=0; n<tree->at(j)[k].neighborList.size(); ++n) {
			int kN = tree->at(j)[k].neighborList[n];
			if((kN != k) && tree->at(j)[kN].Eliminated && tree->at(j)[k].Eliminated) {
				// std::cout << "new_L2P.adjoint(): " << new_L2P.adjoint().rows() << "," << new_L2P.adjoint().cols() << std::endl;
				// std::cout << "tree->at(j)[k].L2P: " << tree->at(j)[k].L2P.rows() << "," << tree->at(j)[k].L2P.cols() << std::endl;
				// std::cout << "tree->at(j)[k].M2L[kN]: " << tree->at(j)[k].M2L[kN].rows() << "," << tree->at(j)[k].M2L[kN].cols() << std::endl;
				tree->at(j)[k].M2L[kN] = new_L2P.adjoint() * tree->at(j)[k].L2P * tree->at(j)[k].M2L[kN];
				// tree->at(j)[kN].M2L[k] = tree->at(j)[k].M2L[kN].transpose(); // compressing P2P_{nn,k}
			}
		}
		// std::cout << "N done" << std::endl;
		//k's self
		if(tree->at(j)[k].Eliminated) {
			// std::cout << "tree->at(j)[k].M2L[k]: " << tree->at(j)[k].M2L[k].rows() << "," << tree->at(j)[k].M2L[k].cols() << std::endl;
			tree->at(j)[k].M2L[k] = new_L2P.adjoint() * tree->at(j)[k].L2P * tree->at(j)[k].M2L[k];
		}
		// std::cout << "self done" << std::endl;
	}

	void update_other_M2Ls_at_k_old_check(int j, int k, int b, std::vector<Mat>& mats) {
		for (auto const& [key, val] : tree->at(j)[k].P2L) {
			mats.push_back(tree->at(j)[k].L2P*tree->at(j)[k].P2L[key]);
		}

		// #pragma omp parallel for
		for(int in=0; in<tree->at(j)[k].interactionList.size(); ++in) {
				int kIL = tree->at(j)[k].interactionList[in];
				if(kIL == b) { //worked out in update_M2L_and_basis method
					continue;
				}
				// std::cout << "update_other_M2Ls_at_k_old_check k: " << k << "	kIL: " << kIL << std::endl;
				mats.push_back(tree->at(j)[k].L2P*tree->at(j)[k].M2L[kIL]);
		}
		//k's Neighbors
		// #pragma omp parallel for
		for(int n=0; n<tree->at(j)[k].neighborList.size(); ++n) {
			int kN = tree->at(j)[k].neighborList[n];
			if((kN != k) && tree->at(j)[kN].Eliminated && tree->at(j)[k].Eliminated) {
				// std::cout << "update_other_M2Ls_at_k_old_check k: " << k << "	kN: " << kN << std::endl;
				mats.push_back(tree->at(j)[k].L2P*tree->at(j)[k].M2L[kN]);
			}
		}
		//k's self
		if(tree->at(j)[k].Eliminated) {
			// std::cout << "self M2L k: " << tree->at(j)[k].M2L[k].norm() << std::endl;
			// std::cout << "update_other_M2Ls_at_k_old_check k: " << k << std::endl;
			mats.push_back(tree->at(j)[k].L2P*tree->at(j)[k].M2L[k]);
		}
	}

	void update_other_M2Ls_due_L2P_make_basis_unitary(int j, int k, Mat &new_L2P) {
		for(int in=0; in<tree->at(j)[k].interactionList.size(); ++in) {
			int kIL = tree->at(j)[k].interactionList[in];
			tree->at(j)[k].M2L[kIL] = new_L2P.adjoint() * tree->at(j)[k].L2P * tree->at(j)[k].M2L[kIL];
		}
	}


	void update_other_M2Ls_due_P2M_make_basis_unitary(int j, int b, Mat &new_P2M_transpose) {
		for(int in=0; in<tree->at(j)[b].interactionList.size(); ++in) {
			int kIL = tree->at(j)[b].interactionList[in];
			tree->at(j)[kIL].M2L[b] = tree->at(j)[kIL].M2L[b] * tree->at(j)[b].P2M * new_P2M_transpose;
		}
	}

	void update_other_M2Ls_at_k_old(int j, int k, int b, Mat &new_L2P) {
		for (auto const& [key, val] : tree->at(j)[k].P2L) {
			tree->at(j)[k].P2L[key] = new_L2P.adjoint() * tree->at(j)[k].L2P * tree->at(j)[k].P2L[key];
		}

		// #pragma omp parallel for
		for(int in=0; in<tree->at(j)[k].interactionList.size(); ++in) {
				int kIL = tree->at(j)[k].interactionList[in];
				if(kIL == b) { //worked out in update_M2L_and_basis method
					continue;
				}
				tree->at(j)[k].M2L[kIL] = new_L2P.adjoint() * tree->at(j)[k].L2P * tree->at(j)[k].M2L[kIL];
		}
		//k's Neighbors
		// #pragma omp parallel for
		for(int n=0; n<tree->at(j)[k].neighborList.size(); ++n) {
			int kN = tree->at(j)[k].neighborList[n];
			if((kN != k) && tree->at(j)[kN].Eliminated && tree->at(j)[k].Eliminated) {
				// if (k==0 && kN==10) {
				// 	std::cout << "update_other_M2Ls_at_k_old 4: tree->at(j)[k].M2L[kN]: " << tree->at(j)[k].M2L[kN].rows() << "," << tree->at(j)[k].M2L[kN].cols() << std::endl;
				// }
				tree->at(j)[k].M2L[kN] = new_L2P.adjoint() * tree->at(j)[k].L2P * tree->at(j)[k].M2L[kN];
				// tree->at(j)[kN].M2L[k] = tree->at(j)[k].M2L[kN].transpose(); // compressing P2P_{nn,k}
			}
		}
		//k's self
		if(tree->at(j)[k].Eliminated) {
			tree->at(j)[k].M2L[k] = new_L2P.adjoint() * tree->at(j)[k].L2P * tree->at(j)[k].M2L[k];
		}
	}

	void update_other_M2Ls_at_b(int j, int b, Mat &new_P2M_transpose) {
		// #pragma omp parallel for
		for (auto const& [key, val] : tree->at(j)[b].P2L) {
			tree->at(j)[key].M2P[b] = tree->at(j)[key].M2P[b] * tree->at(j)[b].P2M * new_P2M_transpose;
		}

		for(int in=0; in<tree->at(j)[b].interactionList.size(); ++in) {
				int kIL = tree->at(j)[b].interactionList[in];
				// std::cout << "tree->at(j)[kIL].M2L[b]: " << tree->at(j)[kIL].M2L[b].rows() << "," << tree->at(j)[kIL].M2L[b].cols() << std::endl;
				// std::cout << "tree->at(j)[b].P2M: " << tree->at(j)[b].P2M.rows() << "," << tree->at(j)[b].P2M.cols() << std::endl;
				// std::cout << "new_P2M_transpose: " << new_P2M_transpose.rows() << "," << new_P2M_transpose.cols() << std::endl;
				tree->at(j)[kIL].M2L[b] = tree->at(j)[kIL].M2L[b] * tree->at(j)[b].P2M * new_P2M_transpose;
		}
		// std::cout << "IL done" << std::endl;
		//k's Neighbors
		// #pragma omp parallel for
		for(int n=0; n<tree->at(j)[b].neighborList.size(); ++n) {
			int kN = tree->at(j)[b].neighborList[n];
			if(kN != b && tree->at(j)[kN].Eliminated && tree->at(j)[b].Eliminated) {
				// if (kN==0 && b==10) {
				// 	std::cout << "update_other_M2Ls_at_b 3: tree->at(j)[kN].M2L[b]: " << tree->at(j)[kN].M2L[b].rows() << "," << tree->at(j)[kN].M2L[b].cols() << std::endl;
				// }
				tree->at(j)[kN].M2L[b] = tree->at(j)[kN].M2L[b] * tree->at(j)[b].P2M * new_P2M_transpose;
				// tree->at(j)[kN].M2L[k] = tree->at(j)[k].M2L[kN].transpose(); // compressing P2P_{nn,k}
			}
		}
		// std::cout << "N done" << std::endl;
		//k's self
		if(tree->at(j)[b].Eliminated) {
			// std::cout << "tree->at(j)[b].M2L[b]: " << tree->at(j)[b].M2L[b].rows() << "," << tree->at(j)[b].M2L[b].cols() << std::endl;
			tree->at(j)[b].M2L[b] = tree->at(j)[b].M2L[b] * tree->at(j)[b].P2M * new_P2M_transpose;
		}
		// std::cout << "self done" << std::endl;
	}

	void update_other_M2Ls_at_b_old_check(int j, int k, int b, std::vector<Mat>& mats) {
		for (auto const& [key, val] : tree->at(j)[b].P2L) {
			mats.push_back(tree->at(j)[key].M2P[b]*tree->at(j)[b].P2M);
		}

		for(int in=0; in<tree->at(j)[b].interactionList.size(); ++in) {
				int kIL = tree->at(j)[b].interactionList[in];
				if(kIL == k) { //worked out in update_M2L_and_basis method
					continue;
				}
				// std::cout << "update_other_M2Ls_at_b_old_check kIL: " << kIL << "	b: " << b << std::endl;
				mats.push_back(tree->at(j)[kIL].M2L[b]*tree->at(j)[b].P2M);
		}
		for(int n=0; n<tree->at(j)[b].neighborList.size(); ++n) {
			int kN = tree->at(j)[b].neighborList[n];
			if(kN != b && tree->at(j)[kN].Eliminated && tree->at(j)[b].Eliminated) {
				// std::cout << "update_other_M2Ls_at_b_old_check kN: " << kN << "	b: " << b << std::endl;
				mats.push_back(tree->at(j)[kN].M2L[b]*tree->at(j)[b].P2M);
			}
		}
		if(tree->at(j)[b].Eliminated) {
			// std::cout << "self M2L b: " << tree->at(j)[b].M2L[b].norm() << std::endl;
			// std::cout << "update_other_M2Ls_at_b_old_check b: " << b << std::endl;
			mats.push_back(tree->at(j)[b].M2L[b]*tree->at(j)[b].P2M);
		}
	}

	void update_other_M2Ls_at_b_old(int j, int k, int b, Mat &new_P2M_transpose) {
		for (auto const& [key, val] : tree->at(j)[b].P2L) {
			tree->at(j)[key].M2P[b] = tree->at(j)[key].M2P[b] * tree->at(j)[b].P2M * new_P2M_transpose;
		}
		// #pragma omp parallel for
		for(int in=0; in<tree->at(j)[b].interactionList.size(); ++in) {
				int kIL = tree->at(j)[b].interactionList[in];
				if(kIL == k) { //worked out in update_M2L_and_basis method
					continue;
				}
				tree->at(j)[kIL].M2L[b] = tree->at(j)[kIL].M2L[b] * tree->at(j)[b].P2M * new_P2M_transpose;
				// if (j==3 && kIL==6 && b==10) {
				// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree->at(j)[k].M2L[b]: " << std::endl << tree->at(j)[kIL].M2L[b] << std::endl << std::endl;
				// }
				// tree->at(j)[kIL].M2L[k] = tree->at(j)[k].M2L[kIL].transpose(); // compressing P2P_{nn,k}
		}
		// std::cout << "inner done" << std::endl;
		// std::cout << "outer done" << std::endl;
		//k's Neighbors
		// #pragma omp parallel for
		for(int n=0; n<tree->at(j)[b].neighborList.size(); ++n) {
			int kN = tree->at(j)[b].neighborList[n];
			if(kN != b && tree->at(j)[kN].Eliminated && tree->at(j)[b].Eliminated) {
				// if (kN==0 && b==10) {
				// 	std::cout << "update_other_M2Ls_at_b_old 4: tree->at(j)[kN].M2L[b]: " << tree->at(j)[kN].M2L[b].rows() << "," << tree->at(j)[kN].M2L[b].cols() << std::endl;
				// }
				tree->at(j)[kN].M2L[b] = tree->at(j)[kN].M2L[b] * tree->at(j)[b].P2M * new_P2M_transpose;
				// tree->at(j)[kN].M2L[k] = tree->at(j)[k].M2L[kN].transpose(); // compressing P2P_{nn,k}
			}
		}
		// std::cout << "neighbor done" << std::endl;
		//k's self
		if(tree->at(j)[b].Eliminated) {
			tree->at(j)[b].M2L[b] = tree->at(j)[b].M2L[b] * tree->at(j)[b].P2M * new_P2M_transpose;
		}
		// std::cout << "self done" << std::endl;
	}

	void update_parent_L2P_old(int j, int k, Mat &new_L2P) {
		int kP = k/nChild;//k's parent
		int c_k = k%nChild; // k is c_k^{th} child of k_parent

		Mat U[nChild];
		int row_index = 0;
		for (size_t c = 0; c < nChild; c++) {
			// std::cout << "c: " << c << std::endl;
			// std::cout << "tree->at(j-1)[kP].L2P.size(): " << tree->at(j-1)[kP].L2P.rows() << ", " << tree->at(j-1)[kP].L2P.cols() << std::endl;
			// std::cout << "tree->at(j)[kP*2+c].NumLocals: " << tree->at(j)[kP*2+c].NumLocals << std::endl;
			U[c] = tree->at(j-1)[kP].L2P.block(row_index,0,tree->at(j)[kP*nChild+c].NumLocals,tree->at(j-1)[kP].L2P.cols());
			row_index += tree->at(j)[kP*nChild+c].NumLocals;
		}
		U[c_k] = new_L2P.adjoint() * tree->at(j)[k].L2P * U[c_k];
		int n_rows = 0;
		for (size_t c = 0; c < nChild; c++) {
			n_rows += U[c].rows();
		}
		tree->at(j-1)[kP].L2P = Mat(n_rows, U[0].cols());
		tree->at(j-1)[kP].L2P << U[0], U[1], U[2], U[3], U[4], U[5], U[6], U[7];
	}

	void update_parent_P2M_old(int j, int b, Mat &new_P2M_transpose) {
		int bP = b/nChild;//k's parent
		int c_b = b%nChild; // k is c_b^{th} child of k_parent

		Mat V[nChild];
		int col_index = 0;
		for (size_t c = 0; c < nChild; c++) {
			// std::cout << "c: " << c << std::endl;
			// std::cout << "tree->at(j-1)[bP].P2M.size(): " << tree->at(j-1)[bP].P2M.rows() << ", " << tree->at(j-1)[bP].P2M.cols() << std::endl;
			// std::cout << "tree->at(j)[bP*2+c].NumMultipoles: " << tree->at(j)[bP*2+c].NumMultipoles << std::endl;
			V[c] = tree->at(j-1)[bP].P2M.block(0,col_index,tree->at(j-1)[bP].P2M.rows(),tree->at(j)[bP*nChild+c].NumMultipoles);
			col_index += tree->at(j)[bP*nChild+c].NumMultipoles;
		}
		V[c_b] = V[c_b] * tree->at(j)[b].P2M * new_P2M_transpose;
		int n_cols = 0;
		for (size_t c = 0; c < nChild; c++) {
			n_cols += V[c].cols();
		}
		tree->at(j-1)[bP].P2M = Mat(V[0].rows(), n_cols);
		tree->at(j-1)[bP].P2M << V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7];
	}


	void initialise_rhs() {
		// #pragma omp parallel for
		for (size_t j = startLevel; j <= nLevels; j++) {
			// #pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				tree->at(j)[k].local_rhs = Vec::Zero(tree->at(j)[k].NumMultipoles);
				tree->at(j)[k].multipole_rhs = Vec::Zero(tree->at(j)[k].NumLocals);
			}
		}
	}

	void initialise_phase() {
		// #pragma omp parallel for
		for (size_t j = startLevel; j <= nLevels; j++) {
			// #pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				tree->at(j)[k].NumMultipoles = tree->at(j)[k].P2M.rows();
				tree->at(j)[k].NumLocals = tree->at(j)[k].L2P.cols();
        if (tree->at(j)[k].NumMultipoles != tree->at(j)[k].NumLocals) {
					// std::cout << "j: " << j << "	k: " << k << "	M: " << tree->at(j)[k].NumMultipoles << "	L: " << tree->at(j)[k].NumLocals << std::endl;
				}
				// std::cout << "j: " << j << "	k: " << k << "	M: " << tree->at(j)[k].NumMultipoles << "	L: " << tree->at(j)[k].NumLocals << std::endl;
			}
		}
		initialise_rhs();//initialise local_rhs and multipole_rhs
		// std::cout << "tree->at(3)[385].NumMultipoles: " << tree->at(3)[385].NumMultipoles << std::endl;
		// std::cout << "tree->at(3)[345].M2L[385]: " << tree->at(3)[345].M2L[385].rows() << "," << tree->at(3)[345].M2L[385].cols() << std::endl;
	}

	void initialise_P2P_Leaf_Level() {
		// #pragma omp parallel for
		// std::cout << "tree->at(3)[385].NumMultipoles: " << tree->at(3)[385].NumMultipoles << std::endl;
		// std::cout << "tree->at(3)[345].M2L[385]: " << tree->at(3)[345].M2L[385].rows() << "," << tree->at(3)[345].M2L[385].cols() << std::endl;
		for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
			// #pragma omp parallel for
			for(int n=0; n<tree->at(nLevels)[k].neighborList.size(); ++n) {
				int nn = tree->at(nLevels)[k].neighborList[n];//P2P_{k,nn}
				if(nn != k) {
					tree->at(nLevels)[k].P2P[nn] = K->getMatrix(tree->at(nLevels)[k].chargeLocations, tree->at(nLevels)[nn].chargeLocations);
				}
			}
			tree->at(nLevels)[k].P2P[k] = K->getMatrix(tree->at(nLevels)[k].chargeLocations, tree->at(nLevels)[k].chargeLocations);
		}
	}

	void write_P2P_Leaf_Level() {
		// #pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
			// #pragma omp parallel for
			for(int n=0; n<tree->at(nLevels)[k].neighborList.size(); ++n) {
				int nn = tree->at(nLevels)[k].neighborList[n];//P2P_{k,nn}
				if(nn != k) {
					tree->at(nLevels)[k].P2P[nn] = K->getMatrix(tree->at(nLevels)[k].chargeLocations, tree->at(nLevels)[nn].chargeLocations);

					std::string filename = "P2P_" + std::to_string(k) + "_" + std::to_string(nn);
					std::ofstream myfile;
					myfile.open(filename.c_str());
					myfile << tree->at(nLevels)[k].P2P[nn] << std::endl;

				}
			}
			tree->at(nLevels)[k].P2P[k] = K->getMatrix(tree->at(nLevels)[k].chargeLocations, tree->at(nLevels)[k].chargeLocations);
			std::string filename = "P2P_" + std::to_string(k) + "_" + std::to_string(k);
			std::ofstream myfile;
			myfile.open(filename.c_str());
			myfile << tree->at(nLevels)[k].P2P[k] << std::endl;

		}
	}

	void initialise_P2P_NonLeafLevel(int j) {
		// #pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			// #pragma omp parallel for
			for(int n=0; n<tree->at(j)[k].neighborList.size(); ++n) {
				int nn = tree->at(j)[k].neighborList[n];//P2P_{k,nn}
				// if (j==2 && k==43 && nn==48) {
				// 	std::cout << "k: " << k << "	nn: " << nn << std::endl;
				// }
				if(nn != k) {
					int n_rows = 0;
					int n_cols = 0;

					// int n_rows_temp2 = 0;
					// int n_cols_temp2 = 0;


				for(int a=0; a<nChild; ++a) {
					int n_cols_temp = 0;
					// int n_cols_temp2= 0;
					for(int b=0; b<nChild; ++b) {
						// if (tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].size() > 0) {
						// if (j==2 && k==43 && nn==48 && a==1 && b==1) {
						// 	std::cout << "nChild*k+a: " << nChild*k+a << std::endl;
						// 	std::cout << "nChild*nn+b: " << nChild*nn+b << std::endl;
						// 	std::cout << "tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b]: " << tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b] << std::endl;
						// 	std::cout << "tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b]: " << tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].rows() << "," << tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].cols() << std::endl;
						// 	std::cout << "tree->at(j+1)[nChild*nn+b].NumMultipoles: " << tree->at(j+1)[nChild*nn+b].NumMultipoles << std::endl;
						// }
							n_cols_temp += tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].cols();
							// if (j==2 && k==43 && nn==48) {
							// 	n_cols_temp2 += tree->at(j+1)[nChild*nn+b].NumMultipoles;
							// }
						// }
						// if (j==2 && k==43 && nn==48) {
						// 	std::cout << "a: " << a << "	b: " << b << "	diff: " << tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].cols()-tree->at(j+1)[nChild*nn+b].NumMultipoles << std::endl;
						// 	std::cout << "a: " << a << "	b: " << b << "	n_cols_temp: " << n_cols_temp << std::endl;
						// }
					}
					// if (j==2 && k==43 && nn==48) {
					// 	std::cout << "a: " << a << "	n_cols_temp: " << n_cols_temp << std::endl;
					// }
					// if (n_cols_temp > 0) {
						n_cols = std::max(n_cols, n_cols_temp);
						// break;
					// }
				}
				for(int b=0; b<nChild; ++b) {
					int n_rows_temp = 0;
					// int n_rows_temp2 = 0;
					for(int a=0; a<nChild; ++a) {
						// if (tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].size() > 0) {
							n_rows_temp += tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].rows();
							// if (j==2 && k==43 && nn==48) {
							// 	n_rows_temp2 += tree->at(j+1)[nChild*k+a].NumLocals;
							// 	// std::cout << << std::endl;
							// }
							// if (j==2 && k==43 && nn==48) {
							// 	std::cout << "a: " << a << "	b: " << b << "	diff: " << tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].rows()-tree->at(j+1)[nChild*k+a].NumLocals << std::endl;
							// 	std::cout << "a: " << a << "	b: " << b  << "	n_rows_temp: " << n_rows_temp << std::endl;
							// }
						// }
					}
					// if (j==2 && k==43 && nn==48) {
					// 	std::cout << "b: " << b << "	n_rows_temp: " << n_rows_temp << std::endl;
					// }
					// if (n_rows_temp > 0) {
						// n_rows = n_rows_temp;
						n_rows = std::max(n_rows, n_rows_temp);
					// 	break;
					// }
				}

				// if (j==2 && k==43 && nn==48) {
				// 		std::cout << "n_rows_temp2, n_cols_temp2: " << n_rows << "," << n_cols << std::endl;
				// }


					// for(int a=0; a<nChild; ++a) {
					// 	int b = 0;
					// 	// n_rows += tree->at(j+1)[nChild*k+a].NumLocals;//tree->at(j+1)[4*k+a].M2L[4*nn].rows();
					// 	// if (tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].size() > 0) {
					// 		n_rows += tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].rows();
					// 	// }
					// }
					// for(int b=0; b<nChild; ++b) {
					// 	int a = 0;
					// 	// n_cols += tree->at(j+1)[nChild*nn+b].NumMultipoles;//tree[j+1][4*k].M2L[4*nn+b].cols();
					// 	// if (tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].size() > 0) {
					// 		n_cols += tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].cols();
					// 	// }
					// }
					tree->at(j)[k].P2P[nn] = Mat::Zero(n_rows, n_cols);
					// if (j==2 && k==43 && nn==48) {
					// 	std::cout << "n_rows, n_cols: " << n_rows << "," << n_cols << std::endl;
					// }
					int row_index = 0;
					int col_index = 0;
					for(int a=0; a<nChild; ++a) {
						for(int b=0; b<nChild; ++b) {
							if (tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].size() > 0) {
								// if (j==2 && k==43 && nn==48) {
								// 	std::cout << "a: " << a << "	b: " << b << "	block: " << row_index << "," << col_index << "," << tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].rows() << "," << tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].cols() << std::endl;
								// }
								tree->at(j)[k].P2P[nn].block(row_index,col_index,tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].rows(),tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].cols()) = tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b];
								// tree->at(j)[k].P2P[nn].block(row_index,col_index,tree->at(j+1)[nChild*k+a].NumLocals,tree->at(j+1)[nChild*nn+b].NumMultipoles) = tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b];
								// col_index += tree->at(j+1)[nChild*nn+b].NumMultipoles;//tree->at(j+1)[4*k+a].M2L[4*nn+b].cols();
								col_index += tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b].cols();
							}
						}
						// row_index += tree->at(j+1)[nChild*k+a].NumLocals;//tree->at(j+1)[4*k+a].M2L[4*nn].rows();
						row_index += tree->at(j+1)[nChild*k+a].M2L[nChild*nn+0].rows();//tree->at(j+1)[4*k+a].M2L[4*nn].rows();
						col_index = 0;
					}
				}
			}
			// std::cout << "-----------self-------------------------" << std::endl;
			// P2P_self
			{
				int nn = k;
				int n_rows = 0;
				int n_cols = 0;
				for(int a=0; a<nChild; ++a) {
					n_rows += tree->at(j+1)[nChild*k+a].NumLocals;//tree->at(j+1)[4*k+a].M2L[4*nn].rows();
				}
				for(int b=0; b<nChild; ++b) {
					n_cols += tree->at(j+1)[nChild*nn+b].NumMultipoles;//tree[j+1][4*k].M2L[2*nn+b].cols();
				}
				tree->at(j)[k].P2P[nn] = Mat::Zero(n_rows, n_cols);
				int row_index = 0;
				int col_index = 0;
				for(int a=0; a<nChild; ++a) {
					for(int b=0; b<nChild; ++b) {
						tree->at(j)[k].P2P[nn].block(row_index,col_index,tree->at(j+1)[nChild*k+a].NumLocals,tree->at(j+1)[nChild*nn+b].NumMultipoles) = tree->at(j+1)[nChild*k+a].M2L[nChild*nn+b];
						col_index += tree->at(j+1)[nChild*nn+b].NumMultipoles;//tree->at(j+1)[4*k+a].M2L[4*nn+b].cols();
					}
					row_index += tree->at(j+1)[nChild*k+a].NumLocals;//tree->at(j+1)[4*k+a].M2L[4*nn].rows();
					col_index = 0;
				}
			}
		}
    // std::cout << "-----------DONE-------------------------" << std::endl;
		// assign_particle_rhs_NonLeafLevel(j);
    // std::cout << "-----------DONE-------------------------" << std::endl;
	}

	void assign_particle_rhs_NonLeafLevel(int j) {
    // #pragma omp parallel for
    for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
			int particle_rhs_length = 0;
			for (size_t c = 0; c < nChild; c++) {
				// particle_rhs_length += tree->at(j+1)[4*k+c].NumLocals;
        particle_rhs_length += tree->at(j+1)[nChild*k+c].multipole_rhs.size();
			}
			tree->at(j)[k].particle_rhs = Vec(particle_rhs_length);
			tree->at(j)[k].particle_rhs << tree->at(j+1)[nChild*k+0].multipole_rhs, tree->at(j+1)[nChild*k+1].multipole_rhs, tree->at(j+1)[nChild*k+2].multipole_rhs, tree->at(j+1)[nChild*k+3].multipole_rhs, tree->at(j+1)[nChild*k+4].multipole_rhs, tree->at(j+1)[nChild*k+5].multipole_rhs, tree->at(j+1)[nChild*k+6].multipole_rhs, tree->at(j+1)[nChild*k+7].multipole_rhs;
		}
	}

  void rhs_eliminate_cluster(int j, int k) {
		rhs_eliminate_x(j, k);
		// std::cout << "eliminate_x done" << std::endl;
		// if (tree->at(j)[k].P2M.rows() > 0) {
			rhs_eliminate_z(j, k);
			// std::cout << "eliminate_z done" << std::endl;
		// }
		// else {
		// 	std::cout << "eliminate_z NOT done" << std::endl;
		// }
	}

	void eliminate_cluster(int j, int k) {
		// double start1 = omp_get_wtime();
		eliminate_x(j, k);
		// double end = omp_get_wtime();
		// time_In_eliminate_x += end-start1;
		// if (j==6 && k==14) {
		// 	std::cout << "eliminate_x done" << std::endl;
		// }
		// if (tree->at(j)[k].P2M.rows() > 0) {
		// double start2 = omp_get_wtime();
		eliminate_z(j, k);
		// end = omp_get_wtime();
		// time_In_eliminate_z += end-start2;
		// time_In_eliminate_cluster += end - start1;

		// if (j==6 && k==14) {
		// 	std::cout << "eliminate_z done" << std::endl;
		// }
		// }
		// else {
		// 	std::cout << "eliminate_z NOT done" << std::endl;
		// }
	}

	void eliminate_x(int j, int k) {
		// in P2P(j,k) we have terms: P2P_SELF, U, P2P_neighbor;
		// self is eliminated, so we have fill-ins/updates with respect to U, P2P_neighbor
		// in equation U.transpose()(j,k)
		// if (j==2 && k==14) {
		// 	std::cout << "tree->at(j)[k].P2P[k].S: " << tree->at(j)[k].P2P[k].rows() << ", " << tree->at(j)[k].P2P[k].cols() << std::endl;
		// }
		if (tree->at(j)[k].P2P[k].size() != 0) {
			Eigen::PartialPivLU<Mat> P2P_self_QR = tree->at(j)[k].P2P[k].lu();
			// double start = omp_get_wtime();
			filings_in_equation_P2M_due_to_x(j,k,P2P_self_QR);//equation z
			// double end = omp_get_wtime();
			// time_In_eliminate_x1 += end-start;
			// if (j==2 && k==14) {
			// 	std::cout << "filings_in_equation_P2M_due_to_x done" << std::endl;
			// }
			// std::cout << "here" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) { //loop over P2P(j,N(k))
				int nn = tree->at(j)[k].neighborList[n];//P2P(j,nn)
				if(nn != k) {
					// if (j==2 && k==14) {
					// 	std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_x START" << std::endl;
					// }
					// double start = omp_get_wtime();
					filings_in_equation_P2P_due_to_x(j,k,nn,P2P_self_QR);
					// if (j==2 && k==14) {
					// 	std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_x DONE" << std::endl;
					// }
					// double end = omp_get_wtime();
					// time_In_eliminate_x2 += end-start;
				}
			}
		}
		else {
			// std::cout << "here" << std::endl;
			// Eigen::FullPivHouseholderQR<Mat> P2P_self_QR = tree->at(j)[k].P2P[k].colPivHouseholderQr();
			// std::cout << "here" << std::endl;
			filings_in_equation_P2M_due_to_x(j,k);//equation z
			// std::cout << "filings_in_equation_P2M_due_to_x done" << std::endl;
			// std::cout << "here" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) { //loop over P2P(j,N(k))
				int nn = tree->at(j)[k].neighborList[n];//P2P(j,nn)
				if(nn != k) {
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_x START" << std::endl;
					filings_in_equation_P2P_due_to_x(j,k,nn);
				}
			}
		}
	}

  void rhs_eliminate_x(int j, int k) {
		// in P2P(j,k) we have terms: P2P_SELF, U, P2P_neighbor;
		// self is eliminated, so we have fill-ins/updates with respect to U, P2P_neighbor
		// in equation U.transpose()(j,k)
		// if (j==3 && k==39) {
			// std::cout << "tree->at(j)[k].P2P[k].S: " << tree->at(j)[k].P2P[k].rows() << ", " << tree->at(j)[k].P2P[k].cols() << std::endl;
		// }
		if (tree->at(j)[k].P2P[k].size() != 0) {
      // std::cout << "x here !=0" << std::endl;
			Eigen::PartialPivLU<Mat> P2P_self_QR = tree->at(j)[k].P2P[k].lu();
			rhs_filings_in_equation_P2M_due_to_x(j,k,P2P_self_QR);//equation z
      // std::cout << "x here" << std::endl;
      // std::cout << "filings_in_equation_P2M_due_to_x done" << std::endl;
			// std::cout << "here" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) { //loop over P2P(j,N(k))
				int nn = tree->at(j)[k].neighborList[n];//P2P(j,nn)
				if(nn != k) {
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_x START" << std::endl;
					rhs_filings_in_equation_P2P_due_to_x(j,k,nn,P2P_self_QR);
          // std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_x DONE" << std::endl;
				}
			}
		}
		else {
			// std::cout << "here 0" << std::endl;
			// Eigen::FullPivHouseholderQR<Mat> P2P_self_QR = tree->at(j)[k].P2P[k].colPivHouseholderQr();
			// std::cout << "here" << std::endl;
			rhs_filings_in_equation_P2M_due_to_x(j,k);//equation z
			// std::cout << "filings_in_equation_P2M_due_to_x done" << std::endl;
			// std::cout << "here" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) { //loop over P2P(j,N(k))
				int nn = tree->at(j)[k].neighborList[n];//P2P(j,nn)
				if(nn != k) {
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_x START" << std::endl;
					rhs_filings_in_equation_P2P_due_to_x(j,k,nn);
				}
			}
		}
	}

////////////////////////////////////////////////////////////////////
void filings_in_equation_P2M_due_to_x(int j, int k) {
	filing_in_equation_P2M_due_to_L2P(j,k);//U fill-in/update; in P2P(j,n) equation
	// std::cout << "filing_in_equation_P2M_due_to_L2P done" << std::endl;
	filing_in_equation_P2M_due_to_P2P_neighbors(j,k);//P2P_neighbor fill-in/update; in P2P(j,n) equation
	// std::cout << "filing_in_equation_P2M_due_to_P2P_neighbors done" << std::endl;
	// rhs_update_in_equation_P2M_due_to_x(j,k);
	// std::cout << "rhs_update_in_equation_P2M_due_to_x done" << std::endl;
}

void rhs_filings_in_equation_P2M_due_to_x(int j, int k) {
	rhs_update_in_equation_P2M_due_to_x(j,k);
	// std::cout << "rhs_update_in_equation_P2M_due_to_x done" << std::endl;
}

void rhs_update_in_equation_P2M_due_to_x(int j, int k) {
	tree->at(j)[k].local_rhs = Vec::Zero(tree->at(j)[k].P2M.rows());//- tree->at(j)[k].P2M * P2P_self_QR.solve(tree->at(j)[k].particle_rhs);
}

void filing_in_equation_P2M_due_to_L2P(int j, int k) {
	tree->at(j)[k].L2M_f = Mat::Zero(tree->at(j)[k].P2M.rows(),tree->at(j)[k].L2P.cols());//-tree->at(j)[k].P2M * P2P_self_QR.solve(tree->at(j)[k].L2P);
}

void filing_in_equation_P2M_due_to_P2P_neighbors(int j, int k) {
	for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) {//neighbors of k
		int nn = tree->at(j)[k].neighborList[n];//P2P(j,nn)
		if (nn != k) {
			if (!tree->at(j)[nn].Eliminated) {
				tree->at(j)[k].P2M_f[nn] = Mat::Zero(tree->at(j)[k].P2M.rows(),tree->at(j)[k].P2P[nn].cols());//-tree->at(j)[k].P2M * P2P_self_QR.solve(tree->at(j)[k].P2P[nn]);
			}
			else {
				tree->at(j)[k].M2M_f[nn] = Mat::Zero(tree->at(j)[k].P2M.rows(),tree->at(j)[k].M2P[nn].cols());//-tree->at(j)[k].P2M * P2P_self_QR.solve(tree->at(j)[k].M2P[nn]);
			}
		}
	}
}

void filings_in_equation_P2P_due_to_x(int j, int k, int nn) {
	filing_due_to_L2P(j,k,nn);//U fill-in/update; in P2P(j,n) equation
	// std::cout << "filing_due_to_L2P done" << std::endl;
	filing_due_to_P2P_neighbors(j,k,nn);//P2P_neighbor fill-in/update; in P2P(j,n) equation
	// std::cout << "filing_due_to_P2P_neighbors done" << std::endl;
	// rhs_update_in_equation_P2P_due_to_x(j,k,nn);
	// std::cout << "rhs_update_in_equation_P2P_due_to_x done" << std::endl;
}

void rhs_filings_in_equation_P2P_due_to_x(int j, int k, int nn) {
	rhs_update_in_equation_P2P_due_to_x(j,k,nn);
	// std::cout << "rhs_update_in_equation_P2P_due_to_x done" << std::endl;
}

void rhs_update_in_equation_P2P_due_to_x(int j, int k, int nn) {
	if (!tree->at(j)[nn].Eliminated) {
		tree->at(j)[nn].particle_rhs = tree->at(j)[nn].particle_rhs;// - tree->at(j)[nn].P2P[k] * P2P_self_QR.solve(tree->at(j)[k].particle_rhs);
	}
	else {
		tree->at(j)[nn].multipole_rhs = tree->at(j)[nn].multipole_rhs;// - tree->at(j)[nn].P2L[k] * P2P_self_QR.solve(tree->at(j)[k].particle_rhs);
	}
}

void filing_due_to_L2P(int j, int k, int nn) {//in P2P(j,nn) equation; nn:eq number
	if (!tree->at(j)[nn].Eliminated) {
		tree->at(j)[nn].L2P_f[k] = Mat::Zero(tree->at(j)[nn].P2P[k].rows(),tree->at(j)[k].L2P.cols());//-tree->at(j)[nn].P2P[k] * P2P_self_QR.solve(tree->at(j)[k].L2P); //temporary; this gets eliminated when z is eliminated
	}
	else {
		// if (k==1 && nn==0) {
		// 	std::cout << "tree->at(j)[nn].P2L[k].S: " << tree->at(j)[nn].P2L[k].rows() << ", " << tree->at(j)[nn].P2L[k].cols() << std::endl;
		// 	std::cout << "tree->at(j)[k].P2P[k].S: " << tree->at(j)[k].P2P[k].rows() << ", " << tree->at(j)[k].P2P[k].cols() << std::endl;
		// 	std::cout << "tree->at(j)[k].L2P.S: " << tree->at(j)[nn].L2P.rows() << ", " << tree->at(j)[nn].L2P.cols() << std::endl;
		// }
		tree->at(j)[nn].L2P_f[k] = Mat::Zero(tree->at(j)[nn].P2L[k].rows(),tree->at(j)[k].L2P.cols());//-tree->at(j)[nn].P2L[k] * P2P_self_QR.solve(tree->at(j)[k].L2P);
	}// 0-K21 K11^{-1} U1
}

void filing_due_to_P2P_neighbors(int j, int k, int nn) {//in P2P(j,nn) equation
	// Eq: P2P_n
	// workout for fill-ins due to all neighbors of (j,k)
	for (size_t p = 0; p < tree->at(j)[k].neighborList.size(); p++) { //loop over all neighbors of (j,k)
		int pp = tree->at(j)[k].neighborList[p];
		if (!tree->at(j)[nn].Eliminated) {
			if (pp != k) {
				if (!tree->at(j)[pp].Eliminated) {
					if (tree->at(j)[nn].P2P[pp].size() == 0) {
						tree->at(j)[nn].P2P[pp] = Mat::Zero(tree->at(j)[nn].P2P[k].rows(),tree->at(j)[k].P2P[pp].cols());//-tree->at(j)[nn].P2P[k] * P2P_self_QR.solve(tree->at(j)[k].P2P[pp]);
					}
					else {
						tree->at(j)[nn].P2P[pp] = tree->at(j)[nn].P2P[pp];// - tree->at(j)[nn].P2P[k] * P2P_self_QR.solve(tree->at(j)[k].P2P[pp]);
					}
					// if(is_well_separated(j,nn,pp)) {
					// 	if (nn < pp) { // using symmetry
					// 		if(tree->at(j)[nn].P2P[pp].norm() > fillInTolerance) {
					// 			std::pair<int, int> g(nn,pp);
					// 			P2P.push_back(g);
					// 		}
					// 	}
					// }
				}
				else {
					if (tree->at(j)[nn].M2P[pp].size() == 0) {
						tree->at(j)[nn].M2P[pp] = Mat::Zero(tree->at(j)[nn].P2P[k].rows(),tree->at(j)[k].M2P[pp].cols());//-tree->at(j)[nn].P2P[k] * P2P_self_QR.solve(tree->at(j)[k].M2P[pp]);
					}
					else {
						tree->at(j)[nn].M2P[pp] = tree->at(j)[nn].M2P[pp];// - tree->at(j)[nn].P2P[k] * P2P_self_QR.solve(tree->at(j)[k].M2P[pp]);
					}
				}
			}
		}
		else {
			if (pp != k) {
				if (!tree->at(j)[pp].Eliminated) {
					if (tree->at(j)[nn].P2L[pp].size() == 0) {
						tree->at(j)[nn].P2L[pp] = Mat::Zero(tree->at(j)[nn].P2L[k].rows(),tree->at(j)[k].P2P[pp].cols());//-tree->at(j)[nn].P2L[k] * P2P_self_QR.solve(tree->at(j)[k].P2P[pp]);
					}
					else {
						tree->at(j)[nn].P2L[pp] = tree->at(j)[nn].P2L[pp];// - tree->at(j)[nn].P2L[k] * P2P_self_QR.solve(tree->at(j)[k].P2P[pp]);
					}
					// if(is_well_separated(j,nn,pp)) {
					// 	if (tree->at(j)[nn].P2L[pp].norm() > fillInTolerance) {
					// 		std::pair<int, int> g(nn,pp);
					// 		P2L_M2P.push_back(g);
					// 	}
					// }
				}
				else {
					tree->at(j)[nn].M2L[pp] = tree->at(j)[nn].M2L[pp];// - tree->at(j)[nn].P2L[k] * P2P_self_QR.solve(tree->at(j)[k].M2P[pp]);
				}
			}
		}
		// p : 2,4; nn: 2
		// Eq 2:(n=2)
		// K22-K21 K11^{-1}K12 :p=2; K22 update
		// -K21 K11^{-1}K14 :p=4; K24 update
		// Eq 4:(n=4)
		// -K41 K11^{-1}K12 :p=2		 -Kn1 K11^{-1}K1p
		// K44-K41 K11^{-1}K14 :p=4		Knp-Kn1 K11^{-1}K1p
	}
}


///////////////////////////////////////////////////////////////////
void rhs_filings_in_equation_P2M_due_to_x(int j, int k, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
  rhs_update_in_equation_P2M_due_to_x(j,k,P2P_self_QR);
  // std::cout << "rhs_update_in_equation_P2M_due_to_x done" << std::endl;
}

	void filings_in_equation_P2M_due_to_x(int j, int k, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
		filing_in_equation_P2M_due_to_L2P(j,k,P2P_self_QR);//U fill-in/update; in P2P(j,n) equation
		// std::cout << "filing_in_equation_P2M_due_to_L2P done" << std::endl;
		filing_in_equation_P2M_due_to_P2P_neighbors(j,k,P2P_self_QR);//P2P_neighbor fill-in/update; in P2P(j,n) equation
		// std::cout << "filing_in_equation_P2M_due_to_P2P_neighbors done" << std::endl;
		// rhs_update_in_equation_P2M_due_to_x(j,k,P2P_self_QR);
		// std::cout << "rhs_update_in_equation_P2M_due_to_x done" << std::endl;
	}

	void rhs_update_in_equation_P2M_due_to_x(int j, int k, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
		tree->at(j)[k].local_rhs = - tree->at(j)[k].P2M * P2P_self_QR.solve(tree->at(j)[k].particle_rhs);
	}

	void filing_in_equation_P2M_due_to_L2P(int j, int k, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
		// std::cout << "tree->at(j)[k].P2M: " << tree->at(j)[k].P2M.rows() << ", " << tree->at(j)[k].P2M.cols() << std::endl;
		// std::cout << "tree->at(j)[k].L2P: " << tree->at(j)[k].L2P.rows() << ", " << tree->at(j)[k].L2P.cols() << std::endl;
		tree->at(j)[k].L2M_f = -tree->at(j)[k].P2M * P2P_self_QR.solve(tree->at(j)[k].L2P);
		// if (j==7 && k==10) {
		// 	std::cout << "filing_in_equation_P2M_due_to_L2P tree->at(j)[k].L2M_f: " << tree->at(j)[k].L2M_f.rows() << "," << tree->at(j)[k].L2M_f.cols() << std::endl;
		// }

	}

	void filing_in_equation_P2M_due_to_P2P_neighbors(int j, int k, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
		// if (j==3 && k==27) {
		// 	std::cout << "tttttt" << std::endl;
		// }
		for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) {//neighbors of k
			int nn = tree->at(j)[k].neighborList[n];//P2P(j,nn)
			// if (j==3 && k==27) {
			// 	std::cout << "k: " << k << "	nn: " << nn << std::endl;
			// 	std::cout << "is_well_separated: " << is_well_separated(j,k,nn) << std::endl;
			// }
			if (nn != k) {
				if (!tree->at(j)[nn].Eliminated) {
					// if (j==3 && k==27) {
						// std::cout << "tree->at(j)[k].P2M.S: " << tree->at(j)[k].P2M.rows() << ", " << tree->at(j)[k].P2M.cols() << std::endl;
						// std::cout << "tree->at(j)[k].P2P[k].S: " << tree->at(j)[k].P2P[k].rows() << ", " << tree->at(j)[k].P2P[k].cols() << std::endl;
						// std::cout << "tree->at(j)[k].P2P[nn].S: " << tree->at(j)[k].P2P[nn].rows() << ", " << tree->at(j)[k].P2P[nn].cols() << std::endl;
					// }
					tree->at(j)[k].P2M_f[nn] = -tree->at(j)[k].P2M * P2P_self_QR.solve(tree->at(j)[k].P2P[nn]);
				}
				else {
					// if (j==3 && k==27) {
					// 	std::cout << "tree->at(j)[k].P2M.S: " << tree->at(j)[k].P2M.rows() << ", " << tree->at(j)[k].P2M.cols() << std::endl;
					// 	std::cout << "tree->at(j)[k].P2P[k].S: " << tree->at(j)[k].P2P[k].rows() << ", " << tree->at(j)[k].P2P[k].cols() << std::endl;
					// 	std::cout << "tree->at(j)[k].M2P[nn].S: " << tree->at(j)[k].M2P[nn].rows() << ", " << tree->at(j)[k].M2P[nn].cols() << std::endl;
					// }
					tree->at(j)[k].M2M_f[nn] = -tree->at(j)[k].P2M * P2P_self_QR.solve(tree->at(j)[k].M2P[nn]);
				}
			}
		}
	}

	void filings_in_equation_P2P_due_to_x(int j, int k, int nn, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
		filing_due_to_L2P(j,k,nn,P2P_self_QR);//U fill-in/update; in P2P(j,n) equation
		// if (j==2 && k==14 && nn==43) {
		// 	std::cout << "filing_due_to_L2P done" << std::endl;
		// }
		filing_due_to_P2P_neighbors(j,k,nn,P2P_self_QR);//P2P_neighbor fill-in/update; in P2P(j,n) equation
		// if (j==2 && k==14 && nn==43) {
		// 	std::cout << "filing_due_to_P2P_neighbors done" << std::endl;
		// }
		// rhs_update_in_equation_P2P_due_to_x(j,k,nn,P2P_self_QR);
		// std::cout << "rhs_update_in_equation_P2P_due_to_x done" << std::endl;
	}

  void rhs_filings_in_equation_P2P_due_to_x(int j, int k, int nn, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
		rhs_update_in_equation_P2P_due_to_x(j,k,nn,P2P_self_QR);
		// std::cout << "rhs_update_in_equation_P2P_due_to_x done" << std::endl;
	}

	void rhs_update_in_equation_P2P_due_to_x(int j, int k, int nn, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
    // if (j==2 && k==8 && nn==10) {
    //   std::cout << "tree->at(j)[nn].Eliminated: " << tree->at(j)[nn].Eliminated << std::endl;
    //   std::cout << "tree->at(j)[nn].P2P[k]: " << tree->at(j)[nn].P2P[k].rows() << "," << tree->at(j)[nn].P2P[k].cols() << std::endl;
    //   std::cout << "tree->at(j)[nn].P2L[k]: " << tree->at(j)[nn].P2L[k].rows() << "," << tree->at(j)[nn].P2L[k].cols() << std::endl;
    //   std::cout << "tree->at(j)[nn].particle_rhs: " << tree->at(j)[nn].particle_rhs.size() << std::endl;
    //   std::cout << "tree->at(j)[nn].multipole_rhs: " << tree->at(j)[nn].multipole_rhs.size() << std::endl;
    // }
		if (!tree->at(j)[nn].Eliminated) {
			tree->at(j)[nn].particle_rhs = tree->at(j)[nn].particle_rhs - tree->at(j)[nn].P2P[k] * P2P_self_QR.solve(tree->at(j)[k].particle_rhs);
		}
		else {
			tree->at(j)[nn].multipole_rhs = tree->at(j)[nn].multipole_rhs - tree->at(j)[nn].P2L[k] * P2P_self_QR.solve(tree->at(j)[k].particle_rhs);
		}
	}

	void filing_due_to_L2P(int j, int k, int nn, Eigen::PartialPivLU<Mat>& P2P_self_QR) {//in P2P(j,nn) equation; nn:eq number
		if (!tree->at(j)[nn].Eliminated) {
			tree->at(j)[nn].L2P_f[k] = -tree->at(j)[nn].P2P[k] * P2P_self_QR.solve(tree->at(j)[k].L2P); //temporary; this gets eliminated when z is eliminated
		}
		else {
			// if (k==1 && nn==0) {
			// 	std::cout << "tree->at(j)[nn].P2L[k].S: " << tree->at(j)[nn].P2L[k].rows() << ", " << tree->at(j)[nn].P2L[k].cols() << std::endl;
			// 	std::cout << "tree->at(j)[k].P2P[k].S: " << tree->at(j)[k].P2P[k].rows() << ", " << tree->at(j)[k].P2P[k].cols() << std::endl;
			// 	std::cout << "tree->at(j)[k].L2P.S: " << tree->at(j)[nn].L2P.rows() << ", " << tree->at(j)[nn].L2P.cols() << std::endl;
			// }
			tree->at(j)[nn].L2P_f[k] = -tree->at(j)[nn].P2L[k] * P2P_self_QR.solve(tree->at(j)[k].L2P);
		}// 0-K21 K11^{-1} U1
	}

	void filing_due_to_P2P_neighbors(int j, int k, int nn, Eigen::PartialPivLU<Mat>& P2P_self_QR) {//in P2P(j,nn) equation
		// Eq: P2P_n
		// workout for fill-ins due to all neighbors of (j,k)
		// if (j==7 && (k==10 || k==11)) {
		// 	std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ tree->at(j)[11].M2P[10]: " << tree->at(j)[11].M2P[10].rows() << ", " << tree->at(j)[11].M2P[10].cols() << std::endl;
		// }
		for (size_t p = 0; p < tree->at(j)[k].neighborList.size(); p++) { //loop over all neighbors of (j,k)
			int pp = tree->at(j)[k].neighborList[p];
			// if (j==7 && k==11) {
				// std::cout << "pp: " << pp << "	E: " << tree->at(j)[pp].Eliminated << std::endl;
				// std::cout << "nn: " << nn << "	E: " << tree->at(j)[nn].Eliminated << std::endl;
			// 	std::cout << "is_well_separated(j,k,nn): " << is_well_separated(j,k,nn) << std::endl;
			// 	std::cout << "is_well_separated(j,nn,pp): " << is_well_separated(j,nn,pp) << std::endl;
			// 	std::cout << "is_well_separated(j,k,pp): " << is_well_separated(j,k,pp) << std::endl;
			// }
			// if (j==2 && k==14 && nn==43) {
			// 	std::cout << "pp: " << pp << "	E: " << tree->at(j)[pp].Eliminated << std::endl;
			// 	std::cout << "nn: " << nn << "	E: " << tree->at(j)[nn].Eliminated << std::endl;
			// 	std::cout << "tree->at(j)[nn].P2P[pp].size(): " << tree->at(j)[nn].P2P[pp].rows() << "," << tree->at(j)[nn].P2P[pp].cols() << std::endl;
			// 	std::cout << "tree->at(j)[k].P2P[k].size(): " << tree->at(j)[k].P2P[k].rows() << "," << tree->at(j)[k].P2P[k].cols() << std::endl;
			// 	std::cout << "tree->at(j)[k].P2P[pp].size(): " << tree->at(j)[k].P2P[pp].rows() << "," << tree->at(j)[k].P2P[pp].cols() << std::endl;
			// }

			if (!tree->at(j)[nn].Eliminated) {
				if (pp != k) {
					if (!tree->at(j)[pp].Eliminated) {
						if (tree->at(j)[nn].P2P[pp].size() == 0) {
							// double start = omp_get_wtime();
							tree->at(j)[nn].P2P[pp] = -tree->at(j)[nn].P2P[k] * P2P_self_QR.solve(tree->at(j)[k].P2P[pp]);
							// double end = omp_get_wtime();
							// timeToSolve += end - start;
							// nSolve++;
						}
						else {
							// double start = omp_get_wtime();
							tree->at(j)[nn].P2P[pp] = tree->at(j)[nn].P2P[pp] - tree->at(j)[nn].P2P[k] * P2P_self_QR.solve(tree->at(j)[k].P2P[pp]);
							// double end = omp_get_wtime();
							// timeToSolve += end - start;
							// nSolve++;
						}
						// if(is_well_separated(j,nn,pp)) {
						// 	if (nn < pp) { // using symmetry
						// 		if(tree->at(j)[nn].P2P[pp].norm() > fillInTolerance) {
						// 			std::pair<int, int> g(nn,pp);
						// 			P2P.push_back(g);
						// 		}
						// 	}
						// }
					}
					else {
						// std::cout << "tree->at(j)[nn].P2P[k].S: " << tree->at(j)[nn].P2P[k].rows() << ", " << tree->at(j)[nn].P2P[k].cols() << std::endl;
						// std::cout << "tree->at(j)[k].P2P[k].S: " << tree->at(j)[k].P2P[k].rows() << ", " << tree->at(j)[k].P2P[k] << std::endl;
						// std::cout << "tree->at(j)[k].M2P[pp].S: " << tree->at(j)[k].M2P[pp].rows() << ", " << tree->at(j)[k].M2P[pp].cols() << std::endl;
						if (tree->at(j)[nn].M2P[pp].size() == 0) {
							// double start = omp_get_wtime();
							tree->at(j)[nn].M2P[pp] = -tree->at(j)[nn].P2P[k] * P2P_self_QR.solve(tree->at(j)[k].M2P[pp]);
							// double end = omp_get_wtime();
							// timeToSolve += end - start;
							// nSolve++;
						}
						else {
							// double start = omp_get_wtime();
							tree->at(j)[nn].M2P[pp] = tree->at(j)[nn].M2P[pp] - tree->at(j)[nn].P2P[k] * P2P_self_QR.solve(tree->at(j)[k].M2P[pp]);
							// double end = omp_get_wtime();
							// timeToSolve += end - start;
							// nSolve++;
						}
						// if (nn==11 && pp==10) {
						// 	std::cout << "1111111111111111111111111111 tree->at(j)[nn].M2P[pp]: " << tree->at(j)[nn].M2P[pp].rows() << ", " << tree->at(j)[nn].M2P[pp].cols() << std::endl;
						// }
					}
				}
			}
			else {
				if (pp != k) {
					if (!tree->at(j)[pp].Eliminated) {
						if (tree->at(j)[nn].P2L[pp].size() == 0) {
							// double start = omp_get_wtime();
							tree->at(j)[nn].P2L[pp] = -tree->at(j)[nn].P2L[k] * P2P_self_QR.solve(tree->at(j)[k].P2P[pp]);
							// double end = omp_get_wtime();
							// timeToSolve += end - start;
							// nSolve++;
						}
						else {
							// double start = omp_get_wtime();
							tree->at(j)[nn].P2L[pp] = tree->at(j)[nn].P2L[pp] - tree->at(j)[nn].P2L[k] * P2P_self_QR.solve(tree->at(j)[k].P2P[pp]);
							// double end = omp_get_wtime();
							// timeToSolve += end - start;
							// nSolve++;
						}
						// if(is_well_separated(j,nn,pp)) {
						// 	if (tree->at(j)[nn].P2L[pp].norm() > fillInTolerance) {
						// 		std::pair<int, int> g(nn,pp);
						// 		P2L_M2P.push_back(g);
						// 	}
						// }
					}
					else {
						// if (j==7 && k==11 && pp==10) {
						// 	std::cout << "tree->at(j)[nn].M2L[pp].S: " << tree->at(j)[nn].M2L[pp].rows() << ", " << tree->at(j)[nn].M2L[pp].cols() << std::endl;
						// 	std::cout << "tree->at(j)[nn].P2L[k].S: " << tree->at(j)[nn].P2L[k].rows() << ", " << tree->at(j)[nn].P2L[k].cols() << std::endl;
						// 	std::cout << "tree->at(j)[k].P2P[k].S: " << tree->at(j)[k].P2P[k].rows() << ", " << tree->at(j)[nn].P2P[k].cols() << std::endl;
						// 	std::cout << "tree->at(j)[k].M2P[pp].S: " << tree->at(j)[k].M2P[pp].rows() << ", " << tree->at(j)[k].M2P[pp].cols() << std::endl;
						// }
						if (tree->at(j)[nn].M2L[pp].size() > 0) {
							// double start = omp_get_wtime();
							tree->at(j)[nn].M2L[pp] = tree->at(j)[nn].M2L[pp] - tree->at(j)[nn].P2L[k] * P2P_self_QR.solve(tree->at(j)[k].M2P[pp]);
							// double end = omp_get_wtime();
							// timeToSolve += end - start;
							// nSolve++;
						}
						else {
							// double start = omp_get_wtime();
							tree->at(j)[nn].M2L[pp] = - tree->at(j)[nn].P2L[k] * P2P_self_QR.solve(tree->at(j)[k].M2P[pp]);
							// double end = omp_get_wtime();
							// timeToSolve += end - start;
							// nSolve++;
						}
						// if (j==3 && nn==6 && pp==10) {
						// 	std::cout << "eliminate_x: k: " << k << std::endl;
						// 	std::cout << "tree->at(j)[nn].P2L[k]: " << std::endl << tree->at(j)[nn].P2L[k] << std::endl << std::endl;
						// 	std::cout << "tree->at(j)[k].P2P[k]: " << std::endl << tree->at(j)[k].P2P[k] << std::endl << std::endl;
						// 	std::cout << "tree->at(j)[k].M2P[pp]: " << std::endl << tree->at(j)[k].M2P[pp] << std::endl << std::endl;
						// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree->at(j)[nn].M2L[pp]: " << std::endl << tree->at(j)[nn].M2L[pp] << std::endl << std::endl;
						// }
					}
				}
			}
			// p : 2,4; nn: 2
			// Eq 2:(n=2)
			// K22-K21 K11^{-1}K12 :p=2; K22 update
			// -K21 K11^{-1}K14 :p=4; K24 update
			// Eq 4:(n=4)
			// -K41 K11^{-1}K12 :p=2		 -Kn1 K11^{-1}K1p
			// K44-K41 K11^{-1}K14 :p=4		Knp-Kn1 K11^{-1}K1p
		}
	}

	bool is_well_separated(int j, int k, int l) {
		if(k==l) { //self
			return false;
		}
		for (size_t i = 0; i < tree->at(j)[k].neighborList.size(); i++) {
			int nn = tree->at(j)[k].neighborList[i];
			if(l==nn) {
				return false;
			}
		}
		return true;
	}

	// bool is_well_separated(int j, int k, int l) {
	// 	bool ILFlag = false;
	// 	bool NFlag = false;
	// 	for (size_t i = 0; i < tree->at(j)[k].interactionList.size(); i++) {
	// 		int nn = tree->at(j)[k].interactionList[i];
	// 		if(l==nn) {
	// 			// ILFlag = true;
	// 			return true;
	// 		}
	// 	}
	// 	// COMMENT below
	// 	for (size_t i = 0; i < tree->at(j)[k].neighborList.size(); i++) {
	// 		int nn = tree->at(j)[k].neighborList[i];
	// 		if(l==nn) {
	// 			NFlag = true;
	// 		}
	// 	}
	// 	if (!ILFlag && !NFlag) {
	// 		std::cout << "############################  j,l,k: " << j << "," << l << "," << k << "   ###############################" << std::endl;
	// 	}
	// 	//
	// 	return false;
	// }

  void eliminate_z(int j, int k) {
		// if(j==2 && k==0) {
		// 	std::cout << "tree->at(j)[k].L2M_f.S: " << tree->at(j)[k].L2M_f.rows() << ", " << tree->at(j)[k].L2M_f.cols() << std::endl;
		// }
		if (tree->at(j)[k].L2M_f.size() > 0) {
			Eigen::PartialPivLU<Mat> L2M_f_QR = tree->at(j)[k].L2M_f.lu();
			// std::cout << "there" << std::endl;
			filings_in_equation_M2L_due_to_z(j,k,L2M_f_QR);
			// std::cout << "filings_in_equation_M2L_due_to_z done" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) { //loop over P2P(j,N(k))
				int nn = tree->at(j)[k].neighborList[n];//P2P(j,nn)
				if(nn != k) {
					filings_in_equation_P2P_due_to_z(j,k,nn,L2M_f_QR);
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_z done" << std::endl;
				}
			}
		}
		else {
			// Eigen::FullPivHouseholderQR<Mat> L2M_f_QR = tree->at(j)[k].L2M_f.colPivHouseholderQr();
			// std::cout << "there" << std::endl;
			filings_in_equation_M2L_due_to_z(j,k);
			// std::cout << "filings_in_equation_M2L_due_to_z done" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) { //loop over P2P(j,N(k))
				int nn = tree->at(j)[k].neighborList[n];//P2P(j,nn)
				if(nn != k) {
					filings_in_equation_P2P_due_to_z(j,k,nn);
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_z done" << std::endl;
				}
			}
		}
	}

	void rhs_eliminate_z(int j, int k) {
		// if(j==2 && k==0) {
			// std::cout << "tree->at(j)[k].L2M_f.S: " << tree->at(j)[k].L2M_f.rows() << ", " << tree->at(j)[k].L2M_f.cols() << std::endl;
		// }
		if (tree->at(j)[k].L2M_f.size() > 0) {
			Eigen::PartialPivLU<Mat> L2M_f_QR = tree->at(j)[k].L2M_f.lu();
			// std::cout << "THERE " << k << std::endl;
			rhs_filings_in_equation_M2L_due_to_z(j,k,L2M_f_QR);
			// std::cout << "filings_in_equation_M2L_due_to_z done" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) { //loop over P2P(j,N(k))
				int nn = tree->at(j)[k].neighborList[n];//P2P(j,nn)
				if(nn != k) {
          // std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_z START" << std::endl;
					rhs_filings_in_equation_P2P_due_to_z(j,k,nn,L2M_f_QR);
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_z done" << std::endl;
				}
			}
		}
		else {
			// Eigen::FullPivHouseholderQR<Mat> L2M_f_QR = tree->at(j)[k].L2M_f.colPivHouseholderQr();
			// std::cout << "there" << std::endl;
			rhs_filings_in_equation_M2L_due_to_z(j,k);
			// std::cout << "filings_in_equation_M2L_due_to_z done" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) { //loop over P2P(j,N(k))
				int nn = tree->at(j)[k].neighborList[n];//P2P(j,nn)
				if(nn != k) {
					rhs_filings_in_equation_P2P_due_to_z(j,k,nn);
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_z done" << std::endl;
				}
			}
		}
	}

	///////////////////////////////////////////////
  void rhs_filings_in_equation_M2L_due_to_z(int j, int k) {
		rhs_update_in_equation_M2L(j,k);
		// std::cout << "rhs_update_in_equation_M2L done" << std::endl;
	}

	void filings_in_equation_M2L_due_to_z(int j, int k) {
		filings_in_equation_M2L_due_to_y(j,k);
		// std::cout << "filings_in_equation_M2L_due_to_y done" << std::endl;
		filings_in_equation_M2L_due_to_neighbors(j,k);
		// std::cout << "filings_in_equation_M2L_due_to_neighbors done" << std::endl;
		// rhs_update_in_equation_M2L(j,k);
		// std::cout << "rhs_update_in_equation_M2L done" << std::endl;
	}

	void rhs_update_in_equation_M2L(int j, int k) {
		tree->at(j)[k].multipole_rhs = Vec(0);//L2M_f_QR.solve(tree->at(j)[k].local_rhs);
	}

	void filings_in_equation_M2L_due_to_y(int j, int k) {
		tree->at(j)[k].M2L[k] = Mat::Zero(0,tree->at(j)[k].NumMultipoles);//L2M_f_QR.solve(-Mat::Identity(tree->at(j)[k].NumMultipoles, tree->at(j)[k].NumMultipoles));
	}

	void filings_in_equation_M2L_due_to_neighbors(int j, int k) {
		for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) {
			int nn = tree->at(j)[k].neighborList[n];//P2P(j,nn)
			if(nn != k) {
				if (!tree->at(j)[nn].Eliminated) {
					tree->at(j)[k].P2L[nn] = Mat::Zero(0,tree->at(j)[k].P2M_f[nn].cols());//L2M_f_QR.solve(tree->at(j)[k].P2M_f[nn]);
					// if(is_well_separated(j,k,nn)) {
					// 	if (tree->at(j)[k].P2L[nn].norm() > fillInTolerance) {
					// 		std::pair<int, int> g(k,nn);
					// 		P2L_M2P.push_back(g);
					// 	}
					// 	// tree->at(j)[k].P2L[nn] = Mat(0,0);
					// }
				}
				else {
					tree->at(j)[k].M2L[nn] = Mat::Zero(0,tree->at(j)[k].M2M_f[nn].cols());//L2M_f_QR.solve(tree->at(j)[k].M2M_f[nn]);
				}
			}
		}
	}

  void rhs_filings_in_equation_P2P_due_to_z(int j, int k, int nn) {
		rhs_update_in_equation_P2P_due_to_z(j,k,nn);
		// std::cout << "rhs_update_in_equation_P2P_due_to_z done" << std::endl;
	}

	void filings_in_equation_P2P_due_to_z(int j, int k, int nn) {
		filings_in_equation_P2P_due_to_y(j,k,nn);
		// std::cout << "filings_in_equation_P2P_due_to_y done" << std::endl;
		filings_in_equation_P2P_due_to_neighbors(j,k,nn);
		// std::cout << "filings_in_equation_P2P_due_to_neighbors done" << std::endl;
		// rhs_update_in_equation_P2P_due_to_z(j,k,nn);
		// std::cout << "rhs_update_in_equation_P2P_due_to_z done" << std::endl;
	}

	void rhs_update_in_equation_P2P_due_to_z(int j, int k, int nn) {
		if (!tree->at(j)[nn].Eliminated) {
			tree->at(j)[nn].particle_rhs = tree->at(j)[nn].particle_rhs;// - tree->at(j)[nn].L2P_f * L2M_f_QR.solve(tree->at(j)[k].local_rhs);
		}
		else {
			tree->at(j)[nn].multipole_rhs = tree->at(j)[nn].multipole_rhs;// - tree->at(j)[nn].L2P_f * L2M_f_QR.solve(tree->at(j)[k].local_rhs);
		}
	}

	void filings_in_equation_P2P_due_to_neighbors(int j, int k, int nn) {
		// nn: equation number
		// fill-ins between neighbors of k
		for (size_t p = 0; p < tree->at(j)[k].neighborList.size(); p++) {
			int pp = tree->at(j)[k].neighborList[p];//P2P(j,nn)
			if (!tree->at(j)[nn].Eliminated) {
				if(pp != k) {
					if (!tree->at(j)[pp].Eliminated) { //case 4
						tree->at(j)[nn].P2P[pp] = tree->at(j)[nn].P2P[pp];// - tree->at(j)[nn].L2P_f * L2M_f_QR.solve(tree->at(j)[k].P2M_f[pp]);
						// if(is_well_separated(j,nn,pp)) {
						// 	if (nn < pp) { // using symmetry
						// 		if(tree->at(j)[nn].P2P[pp].norm() > fillInTolerance) {
						// 			// std::cout << "compress_P2P here" << std::endl;
						// 			// compress_P2P(j,nn,pp);
						// 			std::pair<int, int> g(nn,pp);
						// 			P2P.push_back(g);
						// 			// std::cout << "there" << std::endl;
						// 			// std::cout << "j: " << j << "	nn: " << nn << "	pp: " << pp << "	compress_P2P" << std::endl;
						// 		}
						// 	}
						// 	// tree->at(j)[nn].P2P[pp] = Mat(0,0);
						// }
					}
					else { //case 2/3
						tree->at(j)[nn].M2P[pp] = tree->at(j)[nn].M2P[pp];// - tree->at(j)[nn].L2P_f * L2M_f_QR.solve(tree->at(j)[k].M2M_f[pp]);
						if(is_well_separated(j,nn,pp)) {
							// if (tree->at(j)[nn].M2P[pp].norm() > 1e-10) {
							// 	// std::cout << "there" << std::endl;
							// }
							// tree->at(j)[nn].M2P[pp] = Mat(0,0);
						}
					}
				}
			}
			else {
				if(pp != k) {
					if (!tree->at(j)[pp].Eliminated) { //case 2/3
						tree->at(j)[nn].P2L[pp] = tree->at(j)[nn].P2L[pp];// - tree->at(j)[nn].L2P_f * L2M_f_QR.solve(tree->at(j)[k].P2M_f[pp]);
						// if(is_well_separated(j,nn,pp)) {
						// 	if (tree->at(j)[nn].P2L[pp].norm() > fillInTolerance) {
						// 		std::pair<int, int> g(nn,pp);
						// 		P2L_M2P.push_back(g);
						// 		// std::cout << "there" << std::endl;
						// 	}
						// 	// tree->at(j)[nn].P2L[pp] = Mat(0,0);
						// }
					}
					else { // case 1
						tree->at(j)[nn].M2L[pp] = tree->at(j)[nn].M2L[pp];// - tree->at(j)[nn].L2P_f * L2M_f_QR.solve(tree->at(j)[k].M2M_f[pp]);
					}
				}
			}
		}
	}

	void filings_in_equation_P2P_due_to_y(int j, int k, int nn) {
		// fill-ins between k and its neighbors
		if (!tree->at(j)[nn].Eliminated) { //case 2/3
			if(tree->at(j)[nn].M2P[k].size() == 0) {
				tree->at(j)[nn].M2P[k] = Mat::Zero(tree->at(j)[nn].L2P_f[k].rows(),tree->at(j)[k].NumMultipoles);//- tree->at(j)[nn].L2P_f * L2M_f_QR.solve(-Mat::Identity(tree->at(j)[k].NumMultipoles, tree->at(j)[k].NumMultipoles));
			}
			else {
				tree->at(j)[nn].M2P[k] = tree->at(j)[nn].M2P[k];// - tree->at(j)[nn].L2P_f * L2M_f_QR.solve(-Mat::Identity(tree->at(j)[k].NumMultipoles, tree->at(j)[k].NumMultipoles));
			}
			if(is_well_separated(j,nn,k)) {
				// if (tree->at(j)[nn].M2P[k].norm() > 1e-10) {
				// 	// std::cout << "there" << std::endl;
				// }
				// tree->at(j)[nn].M2P[k] = Mat(0,0);
			}
		}
		else { //case 1
			if (tree->at(j)[nn].M2L[k].size() == 0) {
				tree->at(j)[nn].M2L[k] = Mat::Zero(tree->at(j)[nn].L2P_f[k].rows(), tree->at(j)[k].NumMultipoles);//- tree->at(j)[nn].L2P_f * L2M_f_QR.solve(-Mat::Identity(tree->at(j)[k].NumMultipoles, tree->at(j)[k].NumMultipoles));
			}
			else {
				tree->at(j)[nn].M2L[k] = tree->at(j)[nn].M2L[k];// - tree->at(j)[nn].L2P_f * L2M_f_QR.solve(-Mat::Identity(tree->at(j)[k].NumMultipoles, tree->at(j)[k].NumMultipoles));
			}
		}
	}

	///////////////////////////////////////////////
  void rhs_filings_in_equation_M2L_due_to_z(int j, int k, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		rhs_update_in_equation_M2L(j,k,L2M_f_QR);
		// std::cout << "rhs_update_in_equation_M2L done" << std::endl;
	}

	void filings_in_equation_M2L_due_to_z(int j, int k, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		filings_in_equation_M2L_due_to_y(j,k,L2M_f_QR);
		// std::cout << "filings_in_equation_M2L_due_to_y done" << std::endl;
		filings_in_equation_M2L_due_to_neighbors(j,k,L2M_f_QR);
		// std::cout << "filings_in_equation_M2L_due_to_neighbors done" << std::endl;
		// rhs_update_in_equation_M2L(j,k,L2M_f_QR);
		// std::cout << "rhs_update_in_equation_M2L done" << std::endl;
	}

	void rhs_update_in_equation_M2L(int j, int k, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		// if (j==7 && k==10) {
		// 	std::cout << "tree->at(j)[k].L2M_f.S: " << tree->at(j)[k].L2M_f.rows() << ", " << tree->at(j)[k].L2M_f.cols() << std::endl;
		// 	std::cout << "tree->at(j)[k].local_rhs: " << tree->at(j)[k].local_rhs.size() << std::endl;
		// 	std::cout << "tree->at(j)[k].NumMultipoles: " << tree->at(j)[k].NumMultipoles << std::endl;
		// 	std::cout << "tree->at(j)[k].NumLocals: " << tree->at(j)[k].NumLocals << std::endl;
		// 	std::cout << "tree->at(j)[k].L2P: " << tree->at(j)[k].L2P.rows() << "," << tree->at(j)[k].L2P.cols() << std::endl;
		// 	std::cout << "tree->at(j)[k].P2M: " << tree->at(j)[k].P2M.rows() << "," << tree->at(j)[k].P2M.cols() << std::endl;
		// 	std::cout << "tree->at(j)[k].P2P[k]: " << tree->at(j)[k].P2P[k].rows() << "," << tree->at(j)[k].P2P[k].cols() << std::endl;
		// }
		tree->at(j)[k].multipole_rhs = L2M_f_QR.solve(tree->at(j)[k].local_rhs);
	}

	void filings_in_equation_M2L_due_to_y(int j, int k, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		// std::cout << "tree->at(j)[k].L2M_f.S: " << tree->at(j)[k].L2M_f.rows() << ", " << tree->at(j)[k].L2M_f.cols() << std::endl;
		// std::cout << "tree->at(j)[k].NumMultipoles: " << tree->at(j)[k].NumMultipoles << std::endl;
		tree->at(j)[k].M2L[k] = L2M_f_QR.solve(-Mat::Identity(tree->at(j)[k].NumMultipoles, tree->at(j)[k].NumMultipoles));
		// std::cout << "done" << std::endl;
	}

	void filings_in_equation_M2L_due_to_neighbors(int j, int k, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) {
			int nn = tree->at(j)[k].neighborList[n];//P2P(j,nn)
			if(nn != k) {
				if (!tree->at(j)[nn].Eliminated) {
					tree->at(j)[k].P2L[nn] = L2M_f_QR.solve(tree->at(j)[k].P2M_f[nn]);
					if(is_well_separated(j,k,nn)) {
						if (tree->at(j)[k].M2L[nn].norm() > 0 && tree->at(j)[k].P2L[nn].norm()/tree->at(j)[k].M2L[nn].norm() > fillInTolerance) {
							// std::cout << "k: " << k <<	"	nn: " << nn << "	compress P2L_M2P tree->at(j)[k].P2L[nn].norm(): " << tree->at(j)[k].P2L[nn].norm() << std::endl;
							// std::cout << "P2L_M2P.push_back j: " << j << "	k: " << k << "	nn: " << nn << std::endl;
							std::pair<int, int> g(k,nn);
							P2L_M2P.push_back(g);
						}
						// tree->at(j)[k].P2L[nn] = Mat(0,0);
					}
				}
				else {
					tree->at(j)[k].M2L[nn] = L2M_f_QR.solve(tree->at(j)[k].M2M_f[nn]);
					// if (j==3 && k==6 && nn==10) {
					// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree->at(j)[k].M2L[nn]: " << std::endl << tree->at(j)[k].M2L[nn] << std::endl << std::endl;
					// }
				}
			}
		}
	}

  void rhs_filings_in_equation_P2P_due_to_z(int j, int k, int nn, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
    // std::cout << "rhs_update_in_equation_P2P_due_to_z START" << std::endl;
		rhs_update_in_equation_P2P_due_to_z(j,k,nn,L2M_f_QR);
		// std::cout << "rhs_update_in_equation_P2P_due_to_z done" << std::endl;
	}

	void filings_in_equation_P2P_due_to_z(int j, int k, int nn, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		filings_in_equation_P2P_due_to_y(j,k,nn,L2M_f_QR);
		// std::cout << "filings_in_equation_P2P_due_to_y done" << std::endl;
		filings_in_equation_P2P_due_to_neighbors(j,k,nn,L2M_f_QR);
    // if (k==1 && nn==4) {
    //   std::cout << "tree->at(j)[nn].Eliminated: " << tree->at(j)[nn].Eliminated << std::endl;
    //   std::cout << "tree->at(j)[k].L2M_f: " << tree->at(j)[k].L2M_f.rows() << "," << tree->at(j)[k].L2M_f.cols() << std::endl;
    //   std::cout << "tree->at(j)[nn].L2P_f[k]: " << tree->at(j)[nn].L2P_f[k].rows() << "," << tree->at(j)[nn].L2P_f[k].cols() << std::endl;
    //   std::cout << "tree->at(j)[nn].local_rhs: " << tree->at(j)[k].local_rhs.size() << std::endl;
    //   std::cout << "tree->at(j)[nn].multipole_rhs: " << tree->at(j)[k].multipole_rhs.size() << std::endl;
    //   std::cout << "tree->at(j)[nn].particle_rhs: " << tree->at(j)[k].particle_rhs.size() << std::endl;
    // }
		// std::cout << "filings_in_equation_P2P_due_to_neighbors done" << std::endl;
		// rhs_update_in_equation_P2P_due_to_z(j,k,nn,L2M_f_QR);
		// std::cout << "rhs_update_in_equation_P2P_due_to_z done" << std::endl;
	}

	void rhs_update_in_equation_P2P_due_to_z(int j, int k, int nn, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
    // std::cout << "tree->at(j)[nn].Eliminated: " << tree->at(j)[nn].Eliminated << std::endl;
    // std::cout << "tree->at(j)[k].L2M_f: " << tree->at(j)[k].L2M_f.rows() << "," << tree->at(j)[k].L2M_f.cols() << std::endl;
    // std::cout << "tree->at(j)[nn].L2P_f[k]: " << tree->at(j)[nn].L2P_f[k].rows() << "," << tree->at(j)[nn].L2P_f[k].cols() << std::endl;
    // std::cout << "tree->at(j)[nn].local_rhs: " << tree->at(j)[k].local_rhs.size() << std::endl;
    // std::cout << "tree->at(j)[nn].multipole_rhs: " << tree->at(j)[k].multipole_rhs.size() << std::endl;
    // std::cout << "tree->at(j)[nn].particle_rhs: " << tree->at(j)[k].particle_rhs.size() << std::endl;
		if (!tree->at(j)[nn].Eliminated) {
			tree->at(j)[nn].particle_rhs = tree->at(j)[nn].particle_rhs - tree->at(j)[nn].L2P_f[k] * L2M_f_QR.solve(tree->at(j)[k].local_rhs);
		}
		else {
			tree->at(j)[nn].multipole_rhs = tree->at(j)[nn].multipole_rhs - tree->at(j)[nn].L2P_f[k] * L2M_f_QR.solve(tree->at(j)[k].local_rhs);
		}
	}

	void filings_in_equation_P2P_due_to_neighbors(int j, int k, int nn, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		// nn: equation number
		// fill-ins between neighbors of k
		for (size_t p = 0; p < tree->at(j)[k].neighborList.size(); p++) {
			int pp = tree->at(j)[k].neighborList[p];//P2P(j,nn)
			if (!tree->at(j)[nn].Eliminated) {
				if(pp != k) {
					if (!tree->at(j)[pp].Eliminated) { //case 4
						tree->at(j)[nn].P2P[pp] = tree->at(j)[nn].P2P[pp] - tree->at(j)[nn].L2P_f[k] * L2M_f_QR.solve(tree->at(j)[k].P2M_f[pp]);
						if(is_well_separated(j,nn,pp)) {
							if (nn < pp) { // using symmetry
								// std::cout << "compress_P2P here tree->at(j)[nn].P2P[pp].norm(): " << tree->at(j)[nn].P2P[pp].norm() << std::endl;
								if(tree->at(j)[nn].M2L[pp].norm() > 0 && tree->at(j)[nn].P2P[pp].norm()/tree->at(j)[nn].M2L[pp].norm() > fillInTolerance) {
									// std::cout << "P2P.push_back j: " << j << "	k: " << k << "	nn: " << nn << "	pp: " << pp << " tree->at(j)[nn].P2P[pp].size(): " << tree->at(j)[nn].P2P[pp].size() << std::endl;
									// exit(0);
									std::pair<int, int> g(nn,pp);
									P2P.push_back(g);
								}
							}
							// tree->at(j)[nn].P2P[pp] = Mat(0,0);
						}
					}
					else { //case 2/3
						// if (nn==52 && pp==12) {
						// 	std::cout << "----------------------------------- nn: " << nn <<	"	pp: " << pp << "	tree->at(j)[nn].M2P[pp].size(): " << tree->at(j)[nn].M2P[pp].rows() << "," << tree->at(j)[nn].M2P[pp].cols() << ", " << tree->at(j)[nn].M2P[pp].norm() << std::endl;
						// 	std::cout << "----------------------------------- k: " << k << "	pp: " << pp << "	tree->at(j)[k].M2M_f[pp].size(): " << tree->at(j)[k].M2M_f[pp].rows() << "," << tree->at(j)[k].M2M_f[pp].cols() << ", " << tree->at(j)[k].M2M_f[pp].norm() << std::endl;
						// 	std::cout << "----------------------------------- k: " << k << "	tree->at(j)[k].L2M_f.size(): " << tree->at(j)[k].L2M_f.rows() << "," << tree->at(j)[k].L2M_f.cols() << ", " << tree->at(j)[k].L2M_f.norm() << std::endl;
						// 	// exit(0);
						// }
						tree->at(j)[nn].M2P[pp] = tree->at(j)[nn].M2P[pp] - tree->at(j)[nn].L2P_f[k] * L2M_f_QR.solve(tree->at(j)[k].M2M_f[pp]);
						// if (nn==11 && pp==10) {
						// 	std::cout << "2222222222222222222222222222 tree->at(j)[nn].M2P[pp]: " << tree->at(j)[nn].M2P[pp].rows() << ", " << tree->at(j)[nn].M2P[pp].cols() << std::endl;
						// }
						if(is_well_separated(j,nn,pp)) {
							if (tree->at(j)[nn].M2P[pp].norm() > fillInTolerance) {
								// std::cout << "nn: " << nn <<	"	pp: " << pp << "	compress P2L_M2P tree->at(j)[nn].M2P[pp].norm(): " << tree->at(j)[nn].M2P[pp].norm() << std::endl;
							}
							// tree->at(j)[nn].M2P[pp] = Mat(0,0);
						}
					}
				}
			}
			else {
				if(pp != k) {
					if (!tree->at(j)[pp].Eliminated) { //case 2/3
						tree->at(j)[nn].P2L[pp] = tree->at(j)[nn].P2L[pp] - tree->at(j)[nn].L2P_f[k] * L2M_f_QR.solve(tree->at(j)[k].P2M_f[pp]);
						if(is_well_separated(j,nn,pp)) {
							if (tree->at(j)[nn].M2L[pp].norm() > 0 && tree->at(j)[nn].P2L[pp].norm()/tree->at(j)[nn].M2L[pp].norm() > fillInTolerance) {
								// std::cout << "nn: " << nn <<	"	pp: " << pp << "	compress P2L_M2P tree->at(j)[nn].P2L[pp].norm(): " << tree->at(j)[nn].P2L[pp].norm() << std::endl;
								// std::cout << "P2L_M2P.push_back j: " << j << "	nn: " << nn << "	pp: " << pp << std::endl;
								std::pair<int, int> g(nn,pp);
								P2L_M2P.push_back(g);
							}
							// tree->at(j)[nn].P2L[pp] = Mat(0,0);
						}
					}
					else { // case 1
						tree->at(j)[nn].M2L[pp] = tree->at(j)[nn].M2L[pp] - tree->at(j)[nn].L2P_f[k] * L2M_f_QR.solve(tree->at(j)[k].M2M_f[pp]);
						// if (j==3 && nn==6 && pp==10) {
						// 	std::cout << "eliminate_z: k: " << k << std::endl;
						// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree->at(j)[nn].M2L[pp]: " << std::endl << tree->at(j)[nn].M2L[pp] << std::endl << std::endl;
						// }
					}
				}
			}
		}
	}

	void filings_in_equation_P2P_due_to_y(int j, int k, int nn, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		// fill-ins between k and its neighbors
		if (!tree->at(j)[nn].Eliminated) { //case 2/3
		// 	if (nn==11 && k==10) {
		// 		std::cout << "before 3333333333333333333333333333 tree->at(j)[nn].M2P[k]: " << tree->at(j)[nn].M2P[k].rows() << ", " << tree->at(j)[nn].M2P[k].cols() << std::endl;
		// 		std::cout << "tree->at(j)[nn].L2P_f[k]: " << tree->at(j)[nn].L2P_f[k].rows() << ", " << tree->at(j)[nn].L2P_f[k].cols() << std::endl;
		// 		std::cout << "tree->at(j)[k].L2M_f: " << tree->at(j)[k].L2M_f.rows() << ", " << tree->at(j)[k].L2M_f.cols() << std::endl;
		// }
		if (nn==52 && k==12) {
			std::cout << "----------------------------------- nn: " << nn <<	"	k: " << k << "	tree->at(j)[nn].M2P[k].size(): " << tree->at(j)[nn].M2P[k].rows() << "," << tree->at(j)[nn].M2P[k].cols() << ", " << tree->at(j)[nn].M2P[k].norm() << std::endl;
			std::cout << "----------------------------------- nn: " << nn << "	k: " << k << "	tree->at(j)[nn].L2P_f[k].size(): " << tree->at(j)[nn].L2P_f[k].rows() << "," << tree->at(j)[nn].L2P_f[k].cols() << ", " << tree->at(j)[nn].L2P_f[k].norm() << std::endl;
			std::cout << "----------------------------------- k: " << k << "	tree->at(j)[k].L2M_f.size(): " << tree->at(j)[k].L2M_f.rows() << "," << tree->at(j)[k].L2M_f.cols() << ", " << tree->at(j)[k].L2M_f.norm() << std::endl;
		}
			if(tree->at(j)[nn].M2P[k].size() == 0) {
				tree->at(j)[nn].M2P[k] = - tree->at(j)[nn].L2P_f[k] * L2M_f_QR.solve(-Mat::Identity(tree->at(j)[k].NumMultipoles, tree->at(j)[k].NumMultipoles));
			}
			else {
				tree->at(j)[nn].M2P[k] = tree->at(j)[nn].M2P[k] - tree->at(j)[nn].L2P_f[k] * L2M_f_QR.solve(-Mat::Identity(tree->at(j)[k].NumMultipoles, tree->at(j)[k].NumMultipoles));
			}

			if(is_well_separated(j,nn,k)) {
				if (tree->at(j)[nn].M2P[k].norm() > fillInTolerance) {
					// std::cout << "nn: " << nn <<	"	k: " << k << "	compress P2L_M2P tree->at(j)[nn].M2P[k].norm(): " << tree->at(j)[nn].M2P[k].norm() << std::endl;
				}
				// tree->at(j)[nn].M2P[k] = Mat(0,0);
			}
		}
		else { //case 1
			if (tree->at(j)[nn].M2L[k].size() == 0) {
				tree->at(j)[nn].M2L[k] = - tree->at(j)[nn].L2P_f[k] * L2M_f_QR.solve(-Mat::Identity(tree->at(j)[k].NumMultipoles, tree->at(j)[k].NumMultipoles));
			}
			else {
				tree->at(j)[nn].M2L[k] = tree->at(j)[nn].M2L[k] - tree->at(j)[nn].L2P_f[k] * L2M_f_QR.solve(-Mat::Identity(tree->at(j)[k].NumMultipoles, tree->at(j)[k].NumMultipoles));
			}
			// if (j==3 && nn==6 && k==10) {
			// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree->at(j)[nn].M2L[k]: " << std::endl << tree->at(j)[nn].M2L[k] << std::endl << std::endl;
			// }
		}
	}

	double findMemory() {
		double sum = 0.0; // number of double sized data to be stored
		for(int j=nLevels; j>=startLevel; --j) {
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				sum += tree->at(j)[k].L2P.size();
				sum += tree->at(j)[k].P2M.size();
			}
		}
		// IL
		for(int j=nLevels; j>=startLevel; --j) {
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				for (size_t il = 0; il < tree->at(j)[k].interactionList.size(); il++) {
					int n = tree->at(j)[k].interactionList[il];
					sum += tree->at(j)[k].M2L[n].size();
				}
			}
		}
		// Neighbors
		int j=nLevels;
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			for (size_t il = 0; il < tree->at(j)[k].neighborList.size(); il++) {
				int n = tree->at(j)[k].neighborList[il];
				sum += tree->at(j)[k].P2P[n].size();
			}
		}
		sum = sum/8; // Number of bytes
		sum = sum*1e-9; // memory in GB
		return sum;
	}

	double findMemoryAfterFactorize() {
		double sum = 0.0; // number of double sized data to be stored
		for(int j=nLevels; j>=startLevel; --j) {
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				sum += tree->at(j)[k].L2P.size();
				sum += tree->at(j)[k].P2M.size();
			}
		}
		// IL
		// for(int j=nLevels; j>=startLevel; --j) {
		// 	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
		// 		for (size_t il = 0; il < tree->at(j)[k].interactionList.size(); il++) {
		// 			int n = tree->at(j)[k].interactionList[il];
		// 			sum += tree->at(j)[k].M2L[n].size();
		// 		}
		// 	}
		// }
		// Neighbors
		for(int j=nLevels; j>=startLevel; --j) {
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				for (size_t il = 0; il < tree->at(j)[k].neighborList.size(); il++) {
					int n = tree->at(j)[k].neighborList[il];
					sum += tree->at(j)[k].P2P[n].size();
				}
			}
		}
		sum = sum/8; // Number of bytes
		sum = sum*1e-9; // memory in GB
		return sum;
	}

	void deleteUnnecessary(int j) {
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			// std::cout << "j: " << j << "	k: " << k << std::endl;
			for (auto itr = tree->at(j)[k].M2L.begin(); itr != tree->at(j)[k].M2L.end(); ++itr) {
				tree->at(j)[k].M2L[itr->first] = Mat(0,0);
			}
			// for (auto itr = tree->at(j)[k].P2L.begin(); itr != tree->at(j)[k].P2L.end(); ++itr) {
			// 	tree->at(j)[k].P2L[itr->first] = Mat(0,0);
			// }
			// for (auto itr = tree->at(j)[k].L2P_f.begin(); itr != tree->at(j)[k].L2P_f.end(); ++itr) {
			// 	tree->at(j)[k].L2P_f[itr->first] = Mat(0,0);
			// }
			for (auto itr = tree->at(j)[k].P2M_f.begin(); itr != tree->at(j)[k].P2M_f.end(); ++itr) {
				tree->at(j)[k].P2M_f[itr->first] = Mat(0,0);
			}
			for (auto itr = tree->at(j)[k].M2M_f.begin(); itr != tree->at(j)[k].M2M_f.end(); ++itr) {
				tree->at(j)[k].M2M_f[itr->first] = Mat(0,0);
			}

			// std::cout << "delete done" << std::endl;
			// tree->at(j)[k].L2P_f.clear();
			// tree->at(j)[k].P2M_f.clear();
			// tree->at(j)[k].M2M_f.clear();
		}
	}

	void eliminate_phase_efficient() {
		// std::string filename = "M2P_P2L_check" + std::to_string(N);
		// std::ofstream myfile;
		// myfile.open(filename.c_str());

		// std::cout << "startLevel: " << startLevel << std::endl;
		for(int j=nLevels; j>=startLevel; --j) {
			// std::cout << "j: " << j << "	start" << std::endl;
			// std::cout << "tree->at(3)[385].NumMultipoles: " << tree->at(3)[385].NumMultipoles << std::endl;
			// std::cout << "tree->at(3)[385].P2M.rows(): " << tree->at(3)[385].P2M.rows() << std::endl;
			// std::cout << "tree->at(3)[345].M2L[385]: " << tree->at(3)[345].M2L[385].rows() << "," << tree->at(3)[345].M2L[385].cols() << std::endl;

			if(j != nLevels) {
				initialise_P2P_NonLeafLevel(j);
				////// deleting unnecessary
				// deleteUnnecessary(j+1);
				////// deleting unnecessary
			}
			// std::cout << "initialise_P2P_NonLeafLevel done" << std::endl;
			// std::cout << "tree->at(3)[385].NumMultipoles: " << tree->at(3)[385].NumMultipoles << std::endl;
			// std::cout << "tree->at(3)[345].M2L[385]: " << tree->at(3)[345].M2L[385].rows() << "," << tree->at(3)[345].M2L[385].cols() << std::endl;

			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				// if (k==0) {
					std::cout << "j: " << j << "	k: " << k << std::endl;
				// 	std::cout << "tree->at(3)[385].NumMultipoles: " << tree->at(3)[385].NumMultipoles << std::endl;
				// 	std::cout << "tree->at(3)[385].P2M.rows(): " << tree->at(3)[385].P2M.rows() << std::endl;
				// 	std::cout << "tree->at(3)[345].M2L[385]: " << tree->at(3)[345].M2L[385].rows() << "," << tree->at(3)[345].M2L[385].cols() << std::endl;
				// }
				// 	check22();
				// 	std::cout << "before compression: NumMulitpoles, NumLocals: " << tree->at(j)[k].NumMultipoles << ", " << tree->at(j)[k].NumLocals << std::endl;
				// }
				// if (j==7 && k==11) {
				// 	std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ tree->at(j)[11].M2P[10]: " << tree->at(j)[11].M2P[10].rows() << ", " << tree->at(j)[11].M2P[10].cols() << std::endl;
				// }
				// if (j==6 && k==151527) {
				// if (j==2) {
					// std::cout << "j: " << j << "	k: " << k << "	----------------------------------" << std::endl;
				// }
				// std::cout << "j: " << j << "	k: " << k << "	P2P.size(): " << P2P.size() << "	P2L_M2P.size(): " << P2L_M2P.size() << std::endl;
				// std::cout << "compression START" << std::endl;
				// std::cout << "P2P.size(): " << P2P.size() << std::endl;
				// std::cout << "tree->at(j)[1049].M2L[1256]: " << tree->at(j)[1049].M2L[1256].rows() << "," << tree->at(j)[1049].M2L[1256].cols() << std::endl;
				double start  = omp_get_wtime();
				for (size_t i = 0; i < P2P.size(); i++) {
					std::cout << "P2P compression" << std::endl;
					// std::cout << "P2P size(): " << P2P.size() << "	i: " << i << std::endl;
					int nn = P2P[i].first;
					int pp = P2P[i].second;
					// int cc = 13; // 13, 48, 50, 53
					if (k == nn || k == pp) {
							P2P.erase(P2P.begin()+i);
							i--;
							if (tree->at(j)[nn].P2P[pp].size() != 0) {
								// std::cout << "begin" << std::endl;
								// compress_P2P_Two_Ways_weights(j,nn,pp);
								// if (j==3 && nn==231 && pp==239) {
								// 	std::cout << "compress_P2P j: " << j << " k: " << k << "	nn: " << nn << "	pp: " << pp << ", " << tree->at(j)[nn].P2P[pp].norm() << ", " << tree->at(j)[pp].P2P[nn].norm() << std::endl;
								// }
								// if (j==nLevels) {
								// 	compress_P2P_Two_Ways_projection(j,nn,pp);
								// }
								// else {
									// compress_P2P_Two_Ways_orthogonal_fillin(j,nn,pp);
									compress_P2P_Two_Ways_myRRQR(j,nn,pp); // this one is the one in conference
								// }
								// if (j==2) {
								// 	std::cout << "DONE" << std::endl;
								// }
								tree->at(j)[nn].P2P[pp] = Mat(0,0);
								tree->at(j)[pp].P2P[nn] = Mat(0,0);
							}
					}
				}
				// if (j==2) {
					// std::cout << "P2P compression DONE" << std::endl;
				// }
				for (size_t i = 0; i < P2L_M2P.size(); i++) {
					std::cout << "P2L_M2P compression " << std::endl;
					int nn = P2L_M2P[i].first;
					int pp = P2L_M2P[i].second;
					if (k == nn || k == pp) {
						P2L_M2P.erase(P2L_M2P.begin()+i);
						i--;
						if (tree->at(j)[nn].P2L[pp].size() != 0) {
							// std::cout << "compress_M2P_P2L j: " << j << " k: " << k << "	pp: " << pp << "	nn: " << nn << std::endl;
							// std::cout << "tree->at(j)[pp].M2P[nn]: " << tree->at(j)[pp].M2P[nn].rows() << "," << tree->at(j)[pp].M2P[nn].cols() << std::endl;
							// std::cout << "tree->at(j)[nn].P2L[pp]: " << tree->at(j)[nn].P2L[pp].rows() << "," << tree->at(j)[nn].P2L[pp].cols() << std::endl;
							// if (j==nLevels) {
							// 	compress_M2P_P2L_projection(j, pp, nn);
							// }
							// else {
								// compress_M2P_P2L_orthogonal_fillin(j, pp, nn);
								compress_M2P_P2L_myRRQR(j, pp, nn); // this one is the one in conference
								// compress_M2P_P2L_weights(j, pp, nn);
								// std::cout << "done" << std::endl;
							// }
							tree->at(j)[nn].P2L[pp] = Mat(0,0);
							tree->at(j)[pp].M2P[nn] = Mat(0,0);
						}
					}
				}
				double end  = omp_get_wtime();
				timeInCompressFillin +=  end-start;
				// if (j==6 && k==151527) {
					// std::cout << "P2L_M2P compression DONE" << std::endl;
				// }
				// if (j==7 && k==10) {
				// 	std::cout << "after compression: NumMulitpoles, NumLocals: " << tree->at(j)[k].NumMultipoles << ", " << tree->at(j)[k].NumLocals << std::endl;
				// }
				// if (j==2) {
				// 	std::cout << "compression DONE" << std::endl;
				// }
				// std::cout << "P2P.size(): " << P2P.size() << std::endl;
				// std::cout << "P2L_M2P.size(): " << P2L_M2P.size() << std::endl;
				// P2P.clear();
				// P2L_M2P.clear();
				tree->at(j)[k].Eliminated = true;
				eliminate_cluster(j,k);
				// if (j==2) {
	   			// std::cout << "elimination DONE" << std::endl;
				// }
			}
			// std::cout << "j: " << j << "	done" << std::endl;
		}
		// std::cout << "childLevel: " << childLevel << "	parentLevel: " << parentLevel << " grandParentLevel: " << grandParentLevel << "	greatGrandParentLevel: " << greatGrandParentLevel << std::endl;
	}

  void rhs_eliminate_phase_efficient() {
    for(int j=nLevels; j>=startLevel; --j) {
      // std::cout << "here" << std::endl;
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
        tree->at(j)[k].Eliminated = false;
      }
      // std::cout << "here" << std::endl;
    }
		for(int j=nLevels; j>=startLevel; --j) {
      if(j != nLevels) {
        assign_particle_rhs_NonLeafLevel(j);
        // std::cout << "-----------DONE-------------------------" << std::endl;
			}
      for (int k=0; k<nBoxesPerLevel[j]; ++k) {
        // std::cout << "k: " << k << std::endl;
        tree->at(j)[k].Eliminated = true;
        rhs_eliminate_cluster(j,k);
			}
		}
	}

	void solve_particles_at_base_level(int j) {
		int NumMultipoles = 0;
		int NumLocals = 0;
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			NumMultipoles += tree->at(j)[k].NumMultipoles;
			NumLocals += tree->at(j)[k].NumLocals;
		}
    if (NumLocals == 0 || NumMultipoles == 0) {
      return;
    }
    else {
  		Mat A = Mat::Zero(NumLocals, NumMultipoles);
  		Vec b = Vec::Zero(NumLocals);
  		int r_index = 0;
  		for (int r=0; r<nBoxesPerLevel[j]; ++r) {
  			int s_index = 0;
  			for (int s=0; s<nBoxesPerLevel[j]; ++s) {
  				A.block(r_index, s_index, tree->at(j)[r].NumLocals, tree->at(j)[s].NumMultipoles) = tree->at(j)[r].M2L[s];
  				s_index += tree->at(j)[s].NumMultipoles;
  			}
  			b.segment(r_index, tree->at(j)[r].NumLocals) = tree->at(j)[r].multipole_rhs;
  			r_index += tree->at(j)[r].NumLocals;
  		}
  		// Eigen::PartialPivLU<Mat> A_lu = A.lu();
  		// std::cout << "A.s: " << A.rows() << ", " << A.cols() << std::endl;
  		Eigen::PartialPivLU<Mat> A_qr = A.lu();
  		Vec multipoles = A_qr.solve(b);
			// std::cout << "multipoles: " << std::endl << multipoles << std::endl;
  		// std::cout << "A: " << std::endl << A << std::endl;
  		// std::cout << "b: " << std::endl << b << std::endl;
  		// std::cout << "multipoles.N: " << multipoles.norm() << std::endl;
  		// std::cout << "symmetry err: " << (A-A.transpose()).norm() << std::endl;
  		int index = 0;
  		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
  			tree->at(j)[k].multipoles = multipoles.segment(index, tree->at(j)[k].NumMultipoles);
  			index += tree->at(j)[k].NumMultipoles;
  		}
    }
	}

	void back_substitution_phase() {
    //solve for particles at level 0; which are nothing but the multipoles at level 1
    int j = startLevel;//base level
    solve_particles_at_base_level(j);//now multipoles at level 1 are available
		for (size_t j = startLevel; j <= nLevels; j++) {
			for (int k=nBoxesPerLevel[j]-1; k>=0; k--) {
				// obtain x and z using P2P and U.transpose() equations
				// Eigen::PartialPivLU<Mat> P2P_lu = tree->at(j)[k].P2P[k].lu();
				// std::cout << "j: " << j << "	k: " << k << "	------------------------" << std::endl;
				// std::cout << "tree->at(j)[k].P2P[k].S: " << tree->at(j)[k].P2P[k].rows() << ", " << tree->at(j)[k].P2P[k].cols() << std::endl;
				if (tree->at(j)[k].P2P[k].size() != 0) {
					Eigen::PartialPivLU<Mat> P2P_qr = tree->at(j)[k].P2P[k].lu();
					Vec rhs1;
					rhs1 = tree->at(j)[k].particle_rhs;
					for (size_t n = 0; n < tree->at(j)[k].neighborList.size(); n++) {
						int nn = tree->at(j)[k].neighborList[n];
						if (nn != k) {
							if (!tree->at(j)[nn].Eliminated) {
								rhs1 = rhs1 - tree->at(j)[k].P2P[nn]*tree->at(j)[nn].particles;
							}
							else {
								rhs1 = rhs1 - tree->at(j)[k].M2P[nn]*tree->at(j)[nn].multipoles;
							}
						}
					}
					Vec rhs2 = tree->at(j)[k].multipoles - tree->at(j)[k].P2M * (P2P_qr .solve(rhs1));
					Mat Atemp = -tree->at(j)[k].P2M * (P2P_qr .solve(tree->at(j)[k].L2P));
					// std::cout << "Atemp.S: " << Atemp.rows() << ", " << Atemp.cols() << std::endl;
					if (Atemp.size() != 0) {
						Eigen::PartialPivLU<Mat> Atemp_qr = Atemp.lu();
						tree->at(j)[k].locals = Atemp_qr .solve(rhs2);
					}
					else {
						tree->at(j)[k].locals = Vec(0);
					}
					// tree->at(j)[k].locals = Atemp.lu() .solve(rhs2);
					tree->at(j)[k].particles = P2P_qr .solve(rhs1-tree->at(j)[k].L2P * tree->at(j)[k].locals);
				}
				else {
					tree->at(j)[k].particles = Vec(0);
				}
				tree->at(j)[k].Eliminated = false;
				if (nLevels != j) {
					int index = 0;
					for (size_t c = 0; c < nChild; c++) { // get y^{k+1} from x^{k}
						tree->at(j+1)[nChild*k+c].multipoles = tree->at(j)[k].particles.segment(index, tree->at(j+1)[nChild*k+c].NumMultipoles);
						index += tree->at(j+1)[nChild*k+c].NumMultipoles;
					}
				}
			}
		}
	}

	void check2() {
		std::cout << "multipoles and locals of leaf level boxes" << std::endl << std::endl;
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			std::cout << "k: " << k << ";	" << tree->at(nLevels)[k].NumMultipoles << ", " << tree->at(nLevels)[k].NumLocals << ", " << tree->at(nLevels)[k].NumMultipoles-tree->at(nLevels)[k].NumLocals << std::endl;
		}
	}

  void getx(Vec &x) {
		// std::cout << "getx entered" << std::endl;
    int *indexVec = new int[nBoxesPerLevel[nLevels]];
		// std::cout << "indexVec start" << std::endl;
		int sumCL = 0;
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			// std::cout << "k: " << k << std::endl;
			sumCL += tree->at(nLevels)[k].chargeLocations.size();
    	if(k==0) {
        indexVec[k] = 0;
      }
      else {
        indexVec[k] = indexVec[k-1] + tree->at(nLevels)[k-1].chargeLocations.size();
      }
    }
		// std::cout << "sum of chargeLocations: " << sumCL << std::endl;
		// std::cout << "indexVec computed" << std::endl;

    x = Vec::Zero(N); //all particles
    // #pragma omp parallel for
		int sum = 0;
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			sum += tree->at(nLevels)[k].particles.size();
		}
		// std::cout << "sum of particles: " << sum << std::endl;
    for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			if (indexVec[k] + tree->at(nLevels)[k].particles.size()-1 > N) {
				// std::cout << "k: " << k << "	err occured" << std::endl;
			}
      x.segment(indexVec[k], tree->at(nLevels)[k].particles.size()) = tree->at(nLevels)[k].particles;
      // index += tree->at(nLevels)[k].particles.size();
    }
		delete indexVec;
  }

	// void write_charge_indices() {
	// 	int sizeA = 0;
	// 	int j = nLevels;
	// 	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
	// 		sizeA += tree->at(j)[k].chargeLocations.size() + tree->at(j)[k].L2P.cols() + tree->at(j)[k].P2M.rows();
	// 	}
	// 	Vec x(sizeA);
	// 	int index = 0;
	// 	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
	// 		x.segment(index, tree->at(nLevels)[k].particles.size()) = Vec::Ones(tree->at(nLevels)[k].particles.size());
	// 		index += tree->at(nLevels)[k].particles.size();
	// 		x.segment(index, tree->at(nLevels)[k].locals.size()) = 0*tree->at(nLevels)[k].locals;
	// 		index += tree->at(nLevels)[k].locals.size();
	// 	}
	// 	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
	// 		x.segment(index, tree->at(nLevels)[k].multipoles.size()) = 0*tree->at(nLevels)[k].multipoles;
	// 		index += tree->at(nLevels)[k].multipoles.size();
	// 	}
	// 	std::string filename = "charge_indices_" + std::to_string(N) + ".txt";
	// 	std::ofstream myfile;
	// 	myfile.open(filename.c_str());
	// 	myfile << x << std::endl;
	// }
	//
	// void write_multipole_indices() {
	// 	int sizeA = 0;
	// 	int j = nLevels;
	// 	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
	// 		sizeA += tree->at(j)[k].chargeLocations.size() + tree->at(j)[k].L2P.cols() + tree->at(j)[k].P2M.rows();
	// 	}
	// 	Vec x(sizeA);
	// 	int index = 0;
	// 	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
	// 		x.segment(index, tree->at(nLevels)[k].particles.size()) = Vec::Zero(tree->at(nLevels)[k].particles.size());
	// 		index += tree->at(nLevels)[k].particles.size();
	// 		x.segment(index, tree->at(nLevels)[k].locals.size()) = Vec::Zero(tree->at(nLevels)[k].locals.size());
	// 		index += tree->at(nLevels)[k].locals.size();
	// 	}
	// 	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
	// 		x.segment(index, tree->at(nLevels)[k].multipoles.size()) = Vec::Ones(tree->at(nLevels)[k].multipoles.size());
	// 		index += tree->at(nLevels)[k].multipoles.size();
	// 	}
	// 	std::string filename = "multipole_indices_" + std::to_string(N) + ".txt";
	// 	std::ofstream myfile;
	// 	myfile.open(filename.c_str());
	// 	myfile << x << std::endl;
	// }
	//
	// void write_local_indices() {
	// 	int sizeA = 0;
	// 	int j = nLevels;
	// 	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
	// 		sizeA += tree->at(j)[k].chargeLocations.size() + tree->at(j)[k].L2P.cols() + tree->at(j)[k].P2M.rows();
	// 	}
	// 	Vec x(sizeA);
	// 	int index = 0;
	// 	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
	// 		x.segment(index, tree->at(nLevels)[k].particles.size()) = Vec::Zero(tree->at(nLevels)[k].particles.size());
	// 		index += tree->at(nLevels)[k].particles.size();
	// 		x.segment(index, tree->at(nLevels)[k].locals.size()) = Vec::Ones(tree->at(nLevels)[k].locals.size());
	// 		index += tree->at(nLevels)[k].locals.size();
	// 	}
	// 	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
	// 		x.segment(index, tree->at(nLevels)[k].multipoles.size()) = Vec::Zero(tree->at(nLevels)[k].multipoles.size());
	// 		index += tree->at(nLevels)[k].multipoles.size();
	// 	}
	// 	std::string filename = "local_indices_" + std::to_string(N) + ".txt";
	// 	std::ofstream myfile;
	// 	myfile.open(filename.c_str());
	// 	myfile << x << std::endl;
	// }
	//
	// void write_x() {
	// 	int sizeA = 0;
	// 	int j = nLevels;
	// 	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
	// 		sizeA += tree->at(j)[k].chargeLocations.size() + tree->at(j)[k].L2P.cols() + tree->at(j)[k].P2M.rows();
	// 	}
	// 	Vec x(sizeA);
	// 	int index = 0;
	// 	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
	// 		x.segment(index, tree->at(nLevels)[k].particles.size()) = tree->at(nLevels)[k].particles;
	// 		index += tree->at(nLevels)[k].particles.size();
	// 		x.segment(index, tree->at(nLevels)[k].locals.size()) = tree->at(nLevels)[k].locals;
	// 		index += tree->at(nLevels)[k].locals.size();
	// 	}
	// 	for (int k=0; k<nBoxesPerLevel[j]; ++k) {
	// 		x.segment(index, tree->at(nLevels)[k].multipoles.size()) = tree->at(nLevels)[k].multipoles;
	// 		index += tree->at(nLevels)[k].multipoles.size();
	// 	}
	// 	std::string filename = "x_" + std::to_string(N) + ".txt";
	// 	std::ofstream myfile;
	// 	myfile.open(filename.c_str());
	// 	myfile << x << std::endl;
	// }
	//
	void getLargestLeafSize(std::vector<int>& largestLeafSize) {
		for (int j = nLevels; j >= startLevel; j--) {
			int max = 0;
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				if (max < tree->at(j)[k].chargeLocations.size()) {
					max = tree->at(j)[k].chargeLocations.size();
				}
			}
			largestLeafSize.push_back(max);
		}
	}

	int getAverageLeafSize() {
		int j = nLevels;
		int averageLeafSize = 0;
		int nNonZero = 0;
		for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
			if (tree->at(j)[k].chargeLocations.size() > 0) {
				nNonZero++;
				averageLeafSize += tree->at(j)[k].chargeLocations.size();
			}
		}
		return averageLeafSize/nNonZero;
	}

};
#endif
