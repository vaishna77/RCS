//
//  testFMM2D.cpp
//
//
//  Created by Vaishnavi Gujjula on 1/4/21.
//
//
#include "ACA.hpp"
#include "AIFMM.hpp"
#include <filesystem>

int main(int argc, char* argv[]) {
	int TOL_POW = atoi(argv[1]);
	int fillin_tol = atoi(argv[2]);
	int fillin_rank = atoi(argv[3]); //atoi(argv[2]); // not used currently
	// int nLevelsInput = atoi(argv[3]);
	int leafSize = atoi(argv[4]);
	int numPatchesInSphere = atoi(argv[5]);
	// int numRWGEdges = atoi(argv[5]);
	double electricalSize = atof(argv[6]);
	int n = atoi(argv[7]);  // only for square plate // number of dscretizations made, i.e. h=electricalSize*lambda/n;
	double L = electricalSize*lambda;  // length of side of square plate

	int nLevelsInput = 0; // unassigned
	// int nLevelsInput = log(N/leafSize)/log(4);

	// int leafSize = 0; //atoi(argv[4]); // not used currently
	// int w = atoi(argv[5]); // not used currently

	double start, end;
	// std::cout << sizeof(std::complex<double>) << std::endl;

	omp_set_num_threads(10);
	// omp_set_place(omp_place_threads);
	// Create a parallel region
  // #pragma omp parallel
  // {
  //   // Print the thread ID
  //   printf("Thread ID: %d\n", omp_get_thread_num());
  // }
	// exit(0);
	// omp_set_dynamic(0);     // Explicitly disable dynamic teams
	// omp_set_num_threads(1);

	//////////// ElectrostaticsMatrix ////////////
	inputsToKernelClass inputsToKernel;

	// #ifdef USE_TWO
	// 	// plate - 2D
	// 	inputsToKernel.nodeListFilename = "../../RectangularPlateFromMATLAB/alpha_0_15/nodeList"+ std::to_string(numRWGEdges) +".txt"; // sphere
	//   inputsToKernel.patchListFilename = "../../RectangularPlateFromMATLAB/alpha_0_15/patchList"+ std::to_string(numRWGEdges) +".txt";
	// #elif USE_THREE
	// 	// sphere - 3D
	// 	inputsToKernel.nodeListFilename = "../../SphereDataFromMATLAB/sphereRadiusPoint2/nodeList"+ std::to_string(numPatchesInSphere) +".txt"; // sphere
	// 	inputsToKernel.patchListFilename = "../../SphereDataFromMATLAB/sphereRadiusPoint2/patchList"+ std::to_string(numPatchesInSphere) +".txt";
	// #endif
	inputsToKernel.n = n;
	inputsToKernel.L = L;


	// inputsToKernel.nodeListFilename = "../../CubePatches/nodeListCube"+ std::to_string(numPatchesInSphere) +".txt"; // sphere
	// inputsToKernel.patchListFilename = "../../CubePatches/patchListCube"+ std::to_string(numPatchesInSphere) +".txt";

	inputsToKernel.nodeListFilename = "../inputFiles/radius015/nodeList"+ std::to_string(numPatchesInSphere) +".txt"; // sphere
	inputsToKernel.patchListFilename = "../inputFiles/radius015/patchList"+ std::to_string(numPatchesInSphere) +".txt";

	// inputsToKernel.nodeListFilename = "../../SphereDataFromGmsh/nodeListBusStructure"+ std::to_string(numPatchesInSphere) + ".txt"; // sphere
	// inputsToKernel.patchListFilename = "../../SphereDataFromGmsh/patchListBusStructure"+ std::to_string(numPatchesInSphere) + ".txt";

	// inputsToKernel.nodeListFilename = "../../PCB_geometry/PCB_" + std::to_string(numPatchesInSphere) + "st_conductor_nodeList_36142" +".txt"; // sphere
	// inputsToKernel.patchListFilename = "../../PCB_geometry/PCB_" + std::to_string(numPatchesInSphere) + "st_conductor_patchList_36142" +".txt";

	// inputsToKernel.nodeListFilename = "../../PCB_geometry/PCB_nodeList_"+ std::to_string(numPatchesInSphere) +".txt"; // sphere
	// inputsToKernel.patchListFilename = "../../PCB_geometry/PCB_patchList_"+ std::to_string(numPatchesInSphere) +".txt";

	// inputsToKernel.nodeListFilename = "../../PCB_geometry/PCB2_nodeList.txt"; // sphere
	// inputsToKernel.patchListFilename = "../../PCB_geometry/PCB2_patchList.txt";
	// std::cout << "file accessed" << std::endl;

	// ElectrostaticsMatrix* mykernel		=	new ElectrostaticsMatrix(inputsToKernel);

	EFIEMatrix* mykernel		=	new EFIEMatrix(inputsToKernel);
	int N = mykernel->N;
	std::cout << "N: " << N << std::endl;

	// Eigen::MatrixXcd mat;
	// mat = mykernel->getMatrix(0,0,N,N);
	// std::cout << "mat: " << std::endl << std::fixed << mat << std::endl;
	// Eigen::IOFormat fmt(Eigen::FullPrecision, 0, ", ", "\n", "", "");
	//
	// string filename_real = "mat_real_" + std::to_string(N) + ".txt";
	// std::ofstream myfile_mat_real;
	// myfile_mat_real.open(filename_real.c_str());
	// myfile_mat_real << std::fixed << mat.real().format(fmt) << endl;
	// myfile_mat_real.close();
	//
	// string filename_complex = "mat_complex_" + std::to_string(N) + ".txt";
	// std::ofstream myfile_mat_complex;
	// myfile_mat_complex.open(filename_complex.c_str());
	// myfile_mat_complex << std::fixed << mat.imag().format(fmt) << endl;
	// myfile_mat_complex.close();

	// exit(0);
	////////////////////// look up tabel performance ////////////////////////////////////
	/*int sizeOfLookUpTable = 0;
	for (size_t i = 0; i < 9; i++) {
		sizeOfLookUpTable += mykernel->lookUpTable[i].size();
	}
	std::cout << "sizeOfLookUpTable: " << sizeOfLookUpTable << std::endl;
	std::cout << "N*N: " << N*N << std::endl;

	start = omp_get_wtime();
	Mat originalMatrix(N,N);
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			originalMatrix(i,j) = mykernel->getMatrixEntryFromScratch(i,j);
		}
	}
	end = omp_get_wtime();
	std::cout << "Time to construct matrix from scratch: " << end-start << std::endl;

	start = omp_get_wtime();
	Mat matrixFromLookUpTable(N,N);
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			matrixFromLookUpTable(i,j) = mykernel->getMatrixEntry(i,j);
		}
	}
	end = omp_get_wtime();
	std::cout << "Time to construct matrix from lookUpTable: " << end-start << std::endl;

	Mat Err = originalMatrix - matrixFromLookUpTable;
	std::cout << "Err in lookUpTable: " << Err.norm()/originalMatrix.norm() << std::endl;
	*/
	////////////////////////////////////////////////////////////////////////////////////////

	// mykernel->viewNodeList();
	// mykernel->viewPatchList();
	// mykernel->viewEdgeList();


	int numConductors = mykernel->conductorIndices.size();
	Eigen::MatrixXd midPointOfRWGEdgeMatrix(N,3);
	for (size_t i = 0; i < N; i++) {
		auto itrEdge = next(mykernel->edge2PatchMap.begin(), i);
		std::pair<int, int> edge = itrEdge->first;
		int indexOfEdgeNode1 = edge.first;
		int indexOfEdgeNode2 = edge.second;
		Node node1 = mykernel->nodeList[indexOfEdgeNode1];
	  Node node2 = mykernel->nodeList[indexOfEdgeNode2];
		midPointOfRWGEdgeMatrix.row(i) = (node1+node2)/2.0;
	}
	std::cout << "N: " << N << std::endl;
	std::cout << "No. of conductors: " << numConductors << std::endl;
	//////////// ElectrostaticsMatrix ////////////
	////////////////////////////////////////////////////////////


	Eigen::VectorXd maxCoeff = midPointOfRWGEdgeMatrix.colwise().maxCoeff();
	Eigen::VectorXd minCoeff = midPointOfRWGEdgeMatrix.colwise().minCoeff();

	inputsToSolverClass inputsToSolver;
	inputsToSolver.xmin = minCoeff(0);
	inputsToSolver.xmax = maxCoeff(0);
	inputsToSolver.ymin = minCoeff(1);
	inputsToSolver.ymax = maxCoeff(1);
	inputsToSolver.zmin = minCoeff(2);
	inputsToSolver.zmax = maxCoeff(2);
	inputsToSolver.leafSize = leafSize;
	inputsToSolver.tol_pow = TOL_POW;
	inputsToSolver.locations = midPointOfRWGEdgeMatrix;
	inputsToSolver.nLevelsInput = nLevelsInput;
	inputsToSolver.fillin_rank = fillin_rank;
	inputsToSolver.fillin_tol = fillin_tol;

	// std::cout << "numConductors: " << numConductors << std::endl;

	Node Einc; // Einc is in x direction (x polarization); plane wave is propogating in z direction --> Einc = exp(i\kappa z)\hat{x}
	Einc(0) = 1.0;
	Einc(1) = 0.0;
	Einc(2) = 0.0;
	Eigen::VectorXcd rhs;
	mykernel->getRhs(Einc, rhs);

	// std::cout << "rhs: " << std::endl << std::fixed << rhs << std::endl;
	//
	// string filename_rhs = "rhs_Cpp_" + std::to_string(N) + ".txt";
	// std::ofstream myfile_rhs;
	// myfile_rhs.open(filename_rhs.c_str());
	// myfile_rhs << std::fixed << rhs.real().format(fmt) << endl;
	// myfile_rhs.close();

	Eigen::VectorXcd coefficientsLU;
  mykernel->getCoefficients(rhs, coefficientsLU);
	// std::cout << "coefficients LU: " << std::endl;
	// for (size_t i = 0; i < coefficientsLU.size(); i++) {
	// 	std::cout << std::abs(coefficientsLU[i]) << std::endl;
	// }
	// std::cin.get();


	start		=	omp_get_wtime();

	AIFMM<EFIEMatrix> *aifmm = new AIFMM<EFIEMatrix>(inputsToSolver, mykernel);

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(1);

	// std::cin.get();

	// AIFMM<ElectrostaticsMatrix> *aifmm = new AIFMM<ElectrostaticsMatrix>(inputsToSolver, mykernel, rhs);
	end		=	omp_get_wtime();
	double timeAssemble =	(end-start);

	// std::string rankStructureFilename = "rankStructure" + std::to_string(numRWGEdges) + ".tex";
	// aifmm->viewRankStructure(rankStructureFilename);
	int avgLeafSize = aifmm->A->getAvgLeafSize();
	std::cout << std::endl << "Avg Leaf Size: " << avgLeafSize << std::endl << std::endl;
	std::cout << "Time taken to assemble is: " << timeAssemble << std::endl << std::endl;
	// std::cout << "aifmm->A->nLevels(): " << aifmm->A->nLevels	 << std::endl;
	// std::cout << "initial avg rank: " << aifmm->A->getAvgMultipoles() << std::endl << std::endl;
	std::cout << "initial avg rank: " << std::endl << std::endl;
	int avgRank = aifmm->A->getAvgMultipoles();
	double mem1 = 2*aifmm->A->findMemory(); // 2 for complex double
	std::cout << "Memory for the H2 matrix: " << mem1 << std::endl << std::endl;
	double cr1 = mem1/N/N/16.0/1e-9;
	std::cout << "CR 1: " << cr1 << std::endl << std::endl;


	// checking if the matrix is symmetric ---> no its not symmetric
	// Eigen::MatrixXcd mat;
	// aifmm->A->K->getMatrix(mat, N);
	// Eigen::MatrixXcd ErrS = mat - mat.transpose();
  // double errS = ErrS.norm()/mat.norm();
	// std::cout << "errS: " << errS << std::endl;
	// std::cout << "getMatrixEntry(0,1): " << aifmm->A->K->getMatrixEntry(0,1) << std::endl;
  // std::cout << "getMatrixEntry(1,0): " << aifmm->A->K->getMatrixEntry(1,0) << std::endl;

	start	=	omp_get_wtime();
	aifmm->factorize();
	end		=	omp_get_wtime();
	double timeFactorize =	(end-start);
	std::cout << std::endl << "Time taken to factorize is: " << timeFactorize << std::endl << std::endl;
	double mem2 = 2*aifmm->A->findMemoryAfterFactorize(); // 2 for complex double
	// std::cout << "final avg rank: " << aifmm->A->getAvgMultipoles() << std::endl << std::endl;
	std::cout << "final avg rank: " << std::endl << std::endl;
	int finalAvgRank = aifmm->A->getAvgMultipoles();
	std::cout << "Memory after factorization: " << mem2 << std::endl << std::endl;
	double cr2 = mem2/N/N/16.0/1e-9;
	std::cout << "CR 2: " << cr2 << std::endl << std::endl;

	double timeSolve;
	aifmm->backSubstitute1(rhs); //CONVERGENCE
	start	=	omp_get_wtime();
	aifmm->backSubstitute2(); //CONVERGENCE
	aifmm->backSubstitute3(); //CONVERGENCE
	end		=	omp_get_wtime();
	timeSolve =	(end-start);
	std::cout << std::endl << "Time taken to solve (back substitution phase): " << timeSolve << std::endl << std::endl;

	std::cout << N << " " << TOL_POW << " " << fillin_tol << " " << leafSize << " " << avgLeafSize << " " << avgRank << " " << mem1 << " " << cr1 << " " << finalAvgRank << " " << mem2 << " " << cr2 << " " << std::endl;
	std::cout << timeAssemble << " " << timeFactorize << " " << timeSolve << " " << std::endl;
	Eigen::VectorXcd coefficients;
	aifmm->getPhi(coefficients);



	// std::cout << "coefficients AIFMM: " << std::endl;
	// for (size_t i = 0; i < coefficients.size(); i++) {
	// 	std::cout << std::abs(coefficients[i]) << std::endl;
	// }

	omp_set_dynamic(1);     // Explicitly disable dynamic teams
	omp_set_num_threads(10);


	Eigen::MatrixXcd JAtCentroids;
	// mykernel->getJAtCentroids(coefficients, JAtCentroids);
	mykernel->getJAtCentroids(coefficientsLU, JAtCentroids);
	// std::cout << "JAtCentroids: " << std::endl << JAtCentroids.rowwise().norm() << std::endl;
	// std::cout << "JAtCentroids: " << std::endl << JAtCentroids << std::endl << std::endl  << std::endl;

	{ //constant theta
		double theta = PI/2.0;
		Eigen::VectorXd RCS = Eigen::VectorXd(360);
		Eigen::VectorXd phi = Eigen::VectorXd(360);
		int i = 0;
		for (int phiDegrees = 0; phiDegrees < 360; phiDegrees++) { // azimuthal
			phi(i) = phiDegrees/180.0*PI;
			// RCS(i) = mykernel->getRCS_SevenPointQuadrature(theta, phi(i), coefficients);
			// RCS(i) = mykernel->getRCS_SevenPointQuadrature_One(theta, phi(i), coefficientsLU);
			RCS(i) = mykernel->getRCS_OnePointQuadrature(theta, phi(i), JAtCentroids);
			++i;
		}
		Eigen::MatrixXd saveRCS(RCS.size(),2);
		saveRCS.col(0) = phi;
		saveRCS.col(1) = RCS;
		string filename = "RCS_sphere_theta_90_" + std::to_string(N) + ".txt";
		std::ofstream myfile;
		myfile.open(filename.c_str());
		myfile << saveRCS << endl;
		myfile.close();
	}

	{//constant phi
		double phi = 0.0;
		Eigen::VectorXd RCS = Eigen::VectorXd(181);
		Eigen::VectorXd theta = Eigen::VectorXd(181);
		int i = 0;
		for (int thetaDegrees = 0; thetaDegrees <= 180; thetaDegrees++) { // azimuthal
			theta(i) = thetaDegrees/180.0*PI;
			// RCS(i) = mykernel->getRCS_SevenPointQuadrature_One(theta(i), phi, coefficientsLU);
			// RCS(i) = mykernel->getRCS_SevenPointQuadrature(theta(i), phi, coefficients);
			RCS(i) = mykernel->getRCS_OnePointQuadrature(theta(i), phi, JAtCentroids);
			++i;
		}
		Eigen::MatrixXd saveRCS(RCS.size(),2);
		saveRCS.col(0) = theta;
		saveRCS.col(1) = RCS;
		string filename = "RCS_sphere_phi_0_" + std::to_string(N) + ".txt";
		std::ofstream myfile;
		myfile.open(filename.c_str());
		myfile << saveRCS << endl;
		myfile.close();
	}
	// std::cout << "RCS: " << std::endl << RCS << std::endl;


	// Backward Error computation
	// if (N < 1000) {
		Vec Ax = Vec::Zero(N);
		#pragma opm parallel for
		for (size_t j = 0; j < N; j++) {
			Ax = Ax + coefficients(j) * aifmm->A->K->getCol(N, j);
		}
		Vec backwardErrVec = Ax - rhs;
		double backwardErr = backwardErrVec.norm()/rhs.norm();
		std::cout << "Backward Error: " << backwardErr << std::endl << std::endl;
	// }
	// */



	// std::cout << "coefficients: " << std::endl << coefficients << std::endl;
	// std::cout << "RCS: " << std::endl << RCS << std::endl;
	// exit(0);

	// // if (N > 75000) {
	// if (N > 30000) {
	// 	std::cout << "Capacitance Matrix: " << std::endl << std::endl << capacitanceMatrix << std::endl;
	// 	std::cout << N << " " << aifmm->A->getAvgRank() << " " << timeAssemble << " " << timeFactorize << " " << timeSolve << std::endl;
	// 	exit(0);
	// }
	//
	// Eigen::MatrixXd chargeDensityMatrixUsingLU(N, numConductors);
	// start	=	omp_get_wtime();
	// mykernel->getchargeDensityUsingLU(rhsMatrix, chargeDensityMatrixUsingLU);
	// end		=	omp_get_wtime();
	// Eigen::MatrixXd capacitanceMatrixUsingLU;
	//
	// mykernel->getCapacitanceMatrix(chargeDensityMatrixUsingLU, capacitanceMatrixUsingLU);
	//
	// std::cout << std::endl << "Capacitance Matrix UsingLU: " << std::endl << std::endl << capacitanceMatrixUsingLU << std::endl << std::endl;
  // std::cout << "Capacitance Matrix: " << std::endl << std::endl << capacitanceMatrix << std::endl;
	// std::cout << std::endl << "Error in Capacitance Matrix: " << (capacitanceMatrix-capacitanceMatrixUsingLU).norm() << std::endl;
	// std::cout << std::endl << "Relative error in Capacitance Matrix: " << (capacitanceMatrix-capacitanceMatrixUsingLU).norm()/capacitanceMatrixUsingLU.norm() << std::endl;
	//
	// std::cout << std::endl;
	// for (size_t i = 0; i < numConductors; i++) {
	// 	std::cout << "Forward Error: " << (chargeDensityMatrix.col(i)-chargeDensityMatrixUsingLU.col(i)).norm()/chargeDensityMatrixUsingLU.col(i).norm() << std::endl;
	// }
	// std::cout << N << " " << aifmm->A->getAvgRank() << " " << timeAssemble << " " << timeFactorize << " " << timeSolve << " " << err << std::endl;
	//
	// delete aifmm;
	// delete mykernel;
}
