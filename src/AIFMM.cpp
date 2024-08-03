#include "AIFMM.hpp"
// #include "basicKernel.hpp"
#include "EFIEMatrix.hpp"
// #include "ElectrostaticsMatrix.hpp"

template <typename kerneltype>
  AIFMM<kerneltype>::AIFMM(inputsToSolverClass inputsToSolver, kerneltype* mykernel) {
    double start, end;
    // std::cout << "here" << std::endl;
    A	=	new FMM2DTree<kerneltype>(inputsToSolver, mykernel);
    // std::cout << "constructor done" << std::endl;
    A->createTree();
    // std::cout << "createTree done" << std::endl;
    A->assign_Center_Location();
    // std::cout << "assign_Center_Location done" << std::endl;
    A->assign_Tree_Interactions();
    // std::cout << "assign_Tree_Interactions done" << std::endl;

    A->assignChargeLocations();
    // int averageLeafSize = A->getAverageLeafSize();
    // std::cout << "averageLeafSize: " << averageLeafSize << std::endl;
    // exit(0);

    A->assignNonLeafChargeLocations();//actually it doesnt assign; clears it
    // A->rankGrowthBetweenBoxAndIL();
    // A->rankGrowthBetweenTwoBoxes();
    // exit(0);
    // A->check();
    // std::cout << "here" << std::endl;
    // A->checkILandNSize();
    // std::cout << "checkILandNSize done" << std::endl;
    // std::cout << "---------------- checking if tree interactions are rightly assigned ----------------" << std::endl << std::endl;
    // A->check5();
    // std::cout << "---------------- done ----------------" << std::endl << std::endl;
    // A->check();
  	// A->assignChargeLocations();
    // A->assignNonLeafChargeLocations();//actually it doesnt assign; clears it
    // std::cout << "assignNonLeafChargeLocations done" << std::endl;
    // start	=	omp_get_wtime();
    // std::cin.get();
    int avgLeafSize = A->getAvgLeafSize();
    std::cout << std::endl << "Avg Leaf Size: " << avgLeafSize << std::endl << std::endl;

    A->getNodes();
    // A->checkSize();
    // end		=	omp_get_wtime();
    // double time_getNodes =	(end-start);
    // std::cout << std::endl << "time_getNodes: " << time_getNodes << std::endl;

    // std::cout << "getNodes done" << std::endl;
    // std::cin.get();

    // start	=	omp_get_wtime();
    A->assemble_M2L();
    // std::cin.get();

    A->free_charge_and_checkPoints();
    // std::cin.get();

    // end		=	omp_get_wtime();
    // double time_assemble_M2L =	(end-start);
    // std::cout << std::endl << "time_assemble_M2L: " << time_assemble_M2L << std::endl;

    std::cout << "assemble_M2L done" << std::endl;
    // start	=	omp_get_wtime();
    A->initialise_phase();
    // std::cin.get();

    // end		=	omp_get_wtime();
    // double time_initialise_phase =	(end-start);
    // std::cout << std::endl << "time_initialise_phase: " << time_initialise_phase << std::endl;

    std::cout << "initialise_phase done" << std::endl;
    // start	=	omp_get_wtime();
    A->initialise_P2P_Leaf_Level();
    // A->make_basis_unitary2();
    // A->assign_Leaf_rhs(rhs);

    // A->write_extended_rhs();
    // A->write_extended_sparse_matrix();
    // A->write_NumParticles();
    // A->write_NumMultipoles();

    // A->findErrorInHMatrix();
    // exit(0);
    // A->write_P2P_Leaf_Level();
    // A->write_L2P_P2M_level();
    // end		=	omp_get_wtime();
  	// double time_initialise_P2P_Leaf_Level =	(end-start);
  	// std::cout << std::endl << "time_initialise_P2P_Leaf_Level: " << time_initialise_P2P_Leaf_Level << std::endl;

    std::cout << "initialise_P2P_Leaf_Level done" << std::endl;
  }

template <typename kerneltype>
void AIFMM<kerneltype>::factorize() {
  // std::cout << "here" << std::endl;
	// A->eliminate_phase();
  A->eliminate_phase_efficient();
}

template <typename kerneltype>
Vec AIFMM<kerneltype>::solve(Vec &rhs) {
  A->assign_Leaf_rhs(rhs);
  A->rhs_eliminate_phase_efficient();
  A->back_substitution_phase();
  Vec phi;
  A->getx(phi);
  return phi;
}

template <typename kerneltype>
void AIFMM<kerneltype>::backSubstitute(Vec &rhs) {
  A->assign_Leaf_rhs(rhs);
  A->rhs_eliminate_phase_efficient();
  A->back_substitution_phase();
}

template <typename kerneltype>
void AIFMM<kerneltype>::backSubstitute1(Vec &rhs) {
  A->assign_Leaf_rhs(rhs);
}

template <typename kerneltype>
void AIFMM<kerneltype>::backSubstitute2() {
  A->rhs_eliminate_phase_efficient();
}

template <typename kerneltype>
void AIFMM<kerneltype>::backSubstitute3() {
  A->back_substitution_phase();
}

template <typename kerneltype>
void AIFMM<kerneltype>::getPhi(Vec &phi) {
  // std::cout << "before getx" << std::endl;
  A->getx(phi);
  // std::cout << "after getx" << std::endl;

  // A->write_x();
  // A->write_charge_indices();
  // A->write_multipole_indices();
  // A->write_local_indices();
  A->reorder(phi);
  // std::cout << "reorder done" << std::endl;

}

template <typename kerneltype>
void AIFMM<kerneltype>::viewRankStructure(std::string filename) {
  A->viewRankStructure(filename);
};

// template <typename kerneltype>
// double AIFMM<kerneltype>::getError(Vec &rhs) {
//   A->assign_Leaf_rhs(rhs);
//   return A->error_check();
// }

template <typename kerneltype>
AIFMM<kerneltype>::~AIFMM() {
  delete A;
};

// template class AIFMM<ElectrostaticsMatrix>;
// template class AIFMM<basicKernel>;
template class AIFMM<EFIEMatrix>;

// #ifdef USE_COMPLEX64
// 	template class LowRank<EFIEMatrix>;
// #endif
// #ifdef USE_DOUBLE
// 	template class LowRank<ElectrostaticsMatrix>;
// #endif
