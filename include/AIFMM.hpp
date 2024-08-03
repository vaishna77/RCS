#ifndef __AIFMM_HPP__
#define __AIFMM_HPP__

// #include "basicKernel.hpp"
#include "EFIEMatrix.hpp"
// #include "ElectrostaticsMatrix.hpp"
// #include "AIFMMTree.hpp"
// #include "AIFMMTree2D_3D.hpp"
#include "AIFMMTree2D_3D_memory_optimized.hpp"

template <typename kerneltype>
class AIFMM {
public:
  FMM2DTree<kerneltype> *A;
  AIFMM(inputsToSolverClass inputsToSolver, kerneltype* mykernel);

  void factorize();

  Vec solve(Vec &rhs);

  void backSubstitute(Vec &rhs);

  void backSubstitute1(Vec &rhs);

  void backSubstitute2();

  void backSubstitute3();

  void getPhi(Vec &phi);

  double getError(Vec &rhs);

  void viewRankStructure(std::string filename);

  ~AIFMM();
};

#endif
