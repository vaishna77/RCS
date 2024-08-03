#include <iostream>
#include <algorithm>
#include "ElectrostaticsMatrix.hpp"

ElectrostaticsMatrix::ElectrostaticsMatrix(inputsToKernelClass inputsToKernel) {
  generateNodeList(inputsToKernel.nodeListFilename);
  generatePatchList(inputsToKernel.patchListFilename);
  getSizeOfMatrix();
}

ElectrostaticsMatrix::ElectrostaticsMatrix(double Lx, double Ly, double Lz, int Nx, int Ny) {
  generateNodeList(Lx, Ly, Lz, Nx, Ny);
  generatePatchList(Nx, Ny);
  getSizeOfMatrix();
}

void ElectrostaticsMatrix::generateNodeList(std::string fileName) {
  Node dummy = Node::Zero(3); // indexing in cpp starts at 0, though node indexes start from 1 in the patchlist
  nodeList.push_back(dummy); // so a dummy node is introduced at 0th position
  std::ifstream file(fileName);
  std::string line;
  while(std::getline(file, line)) {
    std::stringstream linestream(line);
    Node node;
    linestream >> node(0) >> node(1) >> node(2);
    nodeList.push_back(node);
  }
  file.close();
}

void ElectrostaticsMatrix::generateNodeList(double Lx, double Ly, double Lz, int Nx, int Ny) {
  std::vector<double> gridPointsx(Nx);
  std::vector<double> gridPointsy(Ny);
  double hx = Lx/(Nx-1);// grid size;
  double hy = Ly/(Ny-1);// grid size;
  gridPointsx[0] = 0.0;
  gridPointsx[0] = 0.0;
  for (size_t i = 1; i < Nx; i++) { // consider the grid points to be uniformly placed between 0 and Lx including 0 and Nx
    gridPointsx[i] = gridPointsx[i-1] + hx;
  }
  gridPointsy[0] = 0.0;
  gridPointsy[0] = 0.0;
  for (size_t i = 1; i < Ny; i++) { // consider the grid points to be uniformly placed between 0 and Lx including 0 and Nx
    gridPointsy[i] = gridPointsy[i-1] + hy;
  }
  Node node;
  for (size_t j = 0; j < Ny; j++) {
    for (size_t i = 0; i < Nx; i++) {
      node[0] = gridPointsx[i];
      node[1] = gridPointsy[j];
      node[2] = 0.0;
      nodeList.push_back(node);
    }
  }
}

void ElectrostaticsMatrix::generatePatchList(std::string fileName) {
  std::ifstream file(fileName);
  std::string line;
  int i = 0;
  while(std::getline(file, line)) {
    std::stringstream linestream(line);
    int indexOfVertex1, indexOfVertex2, indexOfVertex3, conductorIndex;
    linestream >> indexOfVertex1 >> indexOfVertex2 >> indexOfVertex3 >> conductorIndex;
    conductorIndices.insert(conductorIndex);

    Patch p(indexOfVertex1, indexOfVertex2, indexOfVertex3);
    Node vertex1 = nodeList[indexOfVertex1];
    Node vertex2 = nodeList[indexOfVertex2];
    Node vertex3 = nodeList[indexOfVertex3];

    p.area = getAreaOfTriangle(vertex1, vertex2, vertex3);
    p.centroid = getCentroid(vertex1, vertex2, vertex3);
    p.largestEdge = getLargestEdge(vertex1, vertex2, vertex3);
    p.conductorIndex = conductorIndex;
    patchList.push_back(p);
  }
  file.close();
  numConductors = conductorIndices.size();
}

void ElectrostaticsMatrix::getCentroids(Eigen::MatrixX3d& centroids) {
  centroids = Eigen::MatrixXd::Zero(patchList.size(),3);
  for (size_t i = 0; i < patchList.size(); i++) {
    Patch p = patchList[i];
    centroids.row(i) = p.centroid.transpose();
  }
}

void ElectrostaticsMatrix::getSizeOfMatrix() {
  N = patchList.size(); //assigning the size of the matrix
}

void ElectrostaticsMatrix::generatePatchList(int Nx, int Ny) {
  for (size_t j = 0; j < Ny-1; j++) { // to iterate over gridpoints in y direction to form a triangle
    for (size_t i = 0; i < Nx-1; i++) { // to iterate over gridpoints in x direction to form a triangle
      int indexOfVertex1 = j*Nx + i;
      int indexOfVertex2 = indexOfVertex1 + 1;
      int indexOfVertex3 = indexOfVertex1 + Nx;
      Patch p(indexOfVertex1, indexOfVertex2, indexOfVertex3);
      Node vertex1 = nodeList[indexOfVertex1];
      Node vertex2 = nodeList[indexOfVertex2];
      Node vertex3 = nodeList[indexOfVertex3];

      p.area = getAreaOfTriangle(vertex1, vertex2, vertex3);
      p.centroid = getCentroid(vertex1, vertex2, vertex3);
      p.largestEdge = getLargestEdge(vertex1, vertex2, vertex3);
      patchList.push_back(p);

      p.indexOfVertex1 = indexOfVertex3 + 1;
      vertex1 = nodeList[p.indexOfVertex1];

      p.area = getAreaOfTriangle(vertex1, vertex2, vertex3);
      p.centroid = getCentroid(vertex1, vertex2, vertex3);
      p.largestEdge = getLargestEdge(vertex1, vertex2, vertex3);
      patchList.push_back(p);
    }
  }
}

void ElectrostaticsMatrix::viewNodeList() {
  std::cout << "nodeList-------------------------" << std::endl;
  for (size_t i = 0; i < nodeList.size(); i++) {
    Node node = nodeList[i];
    std::cout << node[0] << " " << node[1] << " " << node[2] << std::endl;
  }
}

void ElectrostaticsMatrix::viewPatchList() {
  std::cout << "patchList-------------------------" << std::endl;
  for (size_t i = 0; i < patchList.size(); i++) {
    Patch p = patchList[i];
    int index1 = p.indexOfVertex1;
    int index2 = p.indexOfVertex2;
    int index3 = p.indexOfVertex3;
    std::cout << index1 << " " << index2 << " "<< index3 << " centroid: " << patchList[i].centroid << std::endl;
  }
}

double ElectrostaticsMatrix::getLargestEdge(Node& vertex1, Node& vertex2, Node& vertex3) {
  double lengthOfEdge1 = getDistance(vertex1, vertex2);
  double lengthOfEdge2 = getDistance(vertex2, vertex3);
  double lengthOfEdge3 = getDistance(vertex3, vertex1);
  double lengthOfEdges[3] = {lengthOfEdge1, lengthOfEdge2, lengthOfEdge3};
  double* pointerOflargestEdge = std::max_element(lengthOfEdges, lengthOfEdges+2);
  return *pointerOflargestEdge;
}


Node ElectrostaticsMatrix::getCentroid(Node& vertex1, Node& vertex2, Node& vertex3) {
  Node centroid = (vertex1 + vertex2 + vertex3)/3.0;
  return centroid;
}

double ElectrostaticsMatrix::getAreaOfTriangle(Node& vertex1, Node& vertex2, Node& vertex3) {
  Node edge1 = vertex2 - vertex1;
  Node edge2 = vertex3 - vertex1;
  Node crossProduct = edge1.cross(edge2);
  double area = 0.5*crossProduct.norm();
  return area;
}

double ElectrostaticsMatrix::onePointGaussQuadratureRule(Node centroidTrianglei, Patch& patchj) {
  return patchj.area*ElectrostaticGreensFunction(centroidTrianglei, patchj.centroid);
}

double ElectrostaticsMatrix::sevenPointGaussQuadratureRule(Node centroidTrianglei, double areaOfTrianglej, Node& vertex1Trianglej, Node& vertex2Trianglej, Node& vertex3Trianglej) {
  double alpha[7] = {0.3333333, 0.0597158, 0.470142, 0.470142, 0.7974269, 0.1012865, 0.1012865};
  double beta[7] = {0.3333333, 0.470142, 0.0597158, 0.470142, 0.1012865, 0.7974269, 0.1012865};
  double weights[7] = {0.225, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
  double sum = 0.0;
  for (size_t i = 0; i < 7; i++) {
    Node pointj = vertex1Trianglej + (vertex2Trianglej - vertex1Trianglej)*alpha[i] +(vertex3Trianglej-vertex1Trianglej)*beta[i];
    sum += areaOfTrianglej * weights[i] * ElectrostaticGreensFunction(centroidTrianglei, pointj);
  }
  return sum;
}

double ElectrostaticsMatrix::analyticIntegration(Node& centroidTrianglei, Node& vertex1Trianglej, Node& vertex2Trianglej, Node& vertex3Trianglej) {
  double areaOfTrianglej = getAreaOfTriangle(vertex1Trianglej, vertex2Trianglej, vertex3Trianglej);
  Node edge2Trianglej = vertex2Trianglej - vertex1Trianglej;
  Node edge3Trianglej = vertex3Trianglej - vertex1Trianglej;
  Node edge1Trianglej = centroidTrianglei - vertex1Trianglej;
  Node unitVector = edge2Trianglej.cross(edge3Trianglej);
  unitVector = -1.0*unitVector/unitVector.norm();
  // std::cout << "unitVector: " << std::endl << unitVector <<  std::endl;
  double dotProduct = edge1Trianglej.dot(unitVector);
  Node centroid_on_plane = centroidTrianglei-dotProduct*unitVector;
  // std::cout << "centroid_on_plane: " << std::endl << centroid_on_plane <<  std::endl;
  dotProduct = std::abs(dotProduct);
  // std::cout << "dotProduct: " << dotProduct <<  std::endl;
  Node origin(0,0,0);
  double MoM_entry = 0.0;

  // double **process_nodes = new double*[2];
  // process_nodes[0] = new double[3];
  // process_nodes[1] = new double[3];
  Eigen::Vector3d processNodes1, processNodes2;

  // Calculate the contribution of each edge
  for (int i = 0; i < 3; i++) {
    int flag=1;
    if (i == 0) {
      processNodes1 = vertex1Trianglej;
      processNodes2 = vertex2Trianglej;
    }
    else if (i == 1) {
      processNodes1 = vertex2Trianglej;
      processNodes2 = vertex3Trianglej;
    }
    else if (i == 2) {
      processNodes1 = vertex3Trianglej;
      processNodes2 = vertex1Trianglej;
    }
    double length_side = getDistance(processNodes1, processNodes2);
    // std::cout << "length_side: " << length_side <<  std::endl;
    // Find R_plus
    Eigen::Vector3d control = processNodes1 - centroidTrianglei;
    double R_plus = getDistance(control, origin);
    // std::cout << "R_plus: " << R_plus <<  std::endl;
    // Find R_minus
    control = processNodes2 - centroidTrianglei;
    double R_minus = getDistance(control, origin);
    // std::cout << "R_minus: " << R_minus <<  std::endl;
    // Find P_plus and P_minus
    Eigen::Vector3d P_plus = processNodes1 - centroid_on_plane;
    Eigen::Vector3d P_minus = processNodes2 - centroid_on_plane;
    // std::cout << "P_plus: " << std::endl << P_plus <<  std::endl;
    // std::cout << "P_minus: " << std::endl << P_minus <<  std::endl;
    // Find unit_l
    double distance=getDistance(processNodes1,processNodes2);
    Eigen::Vector3d unit_l = (processNodes1-processNodes2)/distance;
    // std::cout << "unit_l: " << std::endl << unit_l <<  std::endl;

    // Find L_plus and L_minus
    double L_plus = P_plus.dot(unit_l);
    double L_minus = P_minus.dot(unit_l);
    // std::cout << "L_plus: " << L_plus <<  std::endl;
    // std::cout << "L_minus: " << L_minus <<  std::endl;

    // Find unit_u
    Eigen::Vector3d unit_u = unit_l.cross(unitVector);
    // std::cout << "unit_u: " << std::endl << unit_u <<  std::endl;
    // Find P0 and R0
    double P0 = P_plus.dot(unit_u);
    if (P0 < 0) {
      P0 = P0*-1;
    }
    if (P0 < 1.0e-10*length_side) {
          flag = 0;
    }
    // std::cout << "P0: " << P0 <<  std::endl;
    // std::cout << "flag: " << flag <<  std::endl;
    // Find unit_p
    if(flag) {
      double R0 = std::sqrt((P0*P0)+(dotProduct*dotProduct));
      // std::cout << "R0: " << R0 <<  std::endl;
      Eigen::Vector3d unit_p = P_plus-L_plus*unit_l;
  		distance = getDistance(unit_p, origin);
      // std::cout << "distance: " << distance <<  std::endl;
      unit_p = unit_p/distance;
      double T1 = P0*log((R_plus+L_plus)/(R_minus+L_minus));
      // std::cout << "T1: " << T1 <<  std::endl;
      if (P0 < 1e-6*length_side) {
         if(R0 < 1e-4*length_side) {
            T1=0;
         }
      }
      // Calculate the part terms used in integration
      double T2 = std::atan(P0*L_plus/(R0*R0+dotProduct*R_plus));
      // std::cout << "T2: " << T2 <<  std::endl;
      double T3 = std::atan(P0*L_minus/(R0*R0+dotProduct*R_minus));
      // std::cout << "T3: " << T3 <<  std::endl;
      double T4 = unit_p.dot(unit_u);
      // std::cout << "T4: " << T4 <<  std::endl;
      MoM_entry = MoM_entry+(T4*(T1-dotProduct*(T2-T3)));
      // std::cout << "MoM_entry: " << MoM_entry <<  std::endl;
    }
  }
  return MoM_entry;
}

// void ElectrostaticsMatrix::getMatrix(Eigen::MatrixXd& mat) {
//   mat = Eigen::MatrixXd::Zero(N,N);
//   for (size_t i = 0; i < N; i++) {
//     for (size_t j = 0; j < N; j++) {
//       mat(i,j) = getMatrixEntry(i, j);
//     }
//   }
//   return;
// }

double ElectrostaticsMatrix::getMatrixEntry(const unsigned i, const unsigned j) {
  Patch patchi = patchList[i];// jth patch
  Patch patchj = patchList[j];// jth patch
  int index1Trianglej = patchj.indexOfVertex1;
  int index2Trianglej = patchj.indexOfVertex2;
  int index3Trianglej = patchj.indexOfVertex3;
  Node vertex1Trianglej = nodeList[index1Trianglej];
  Node vertex2Trianglej = nodeList[index2Trianglej];
  Node vertex3Trianglej = nodeList[index3Trianglej];
  double distanceBetweenCentroids = getDistance(patchi.centroid, patchj.centroid);
  double matrixEntry;
  // /*
  if(distanceBetweenCentroids < 4.0*patchj.largestEdge) {
    // matrixEntry = C0*analyticIntegration(patchi.centroid, vertex1Trianglej, vertex2Trianglej, vertex3Trianglej);

    // int N = patchList.size();
		// if (i==j) {
    //   matrixEntry = C0*analyticIntegration(patchi.centroid, vertex1Trianglej, vertex2Trianglej, vertex3Trianglej);
    //   // matrixEntry = sqrt(1000*N);
    // }
    // else {
      matrixEntry = C0*analyticIntegration(patchi.centroid, vertex1Trianglej, vertex2Trianglej, vertex3Trianglej);
      // matrixEntry = 1.0/distanceBetweenCentroids;
    // }
  }
  else if ((distanceBetweenCentroids >= 4.0*patchj.largestEdge) && (distanceBetweenCentroids < 10.0*patchj.largestEdge)) {
    matrixEntry = C0*sevenPointGaussQuadratureRule(patchi.centroid, patchj.area, vertex1Trianglej, vertex2Trianglej, vertex3Trianglej);
    // matrixEntry = 1.0/distanceBetweenCentroids;
  }
  else {
    matrixEntry = C0*onePointGaussQuadratureRule(patchi.centroid, patchj);
    // matrixEntry = 1.0/distanceBetweenCentroids;
  }
  return matrixEntry;
  // */
  /*
  if (i==j) {
    // return 100.0;
		int N = patchList.size();
		return pow(1000*N, 1.0/2.0);
    // return 1;
  }
  else {
    // return 0;
    return 1.0/distanceBetweenCentroids;
  }
   */
}

void ElectrostaticsMatrix::getchargeDensityUsingLU(Eigen::VectorXd& rhs, Eigen::VectorXd& chargeDensity) {
  Eigen::MatrixXd mat = getMatrix(0,0,N,N);
  Eigen::PartialPivLU<Eigen::MatrixXd> matLU = mat.lu();
  chargeDensity = matLU.solve(rhs);
}

void ElectrostaticsMatrix::getchargeDensityUsingLU(Eigen::MatrixXd& rhs, Eigen::MatrixXd& chargeDensity) {
  Eigen::MatrixXd mat = getMatrix(0,0,N,N);

  // std::string filename = "matrixSphere_" + std::to_string(N);
  // std::ofstream myfile;
  // myfile.open(filename.c_str());
  // myfile << mat << std::endl;

  Eigen::PartialPivLU<Eigen::MatrixXd> matLU = mat.lu();
  chargeDensity = matLU.solve(rhs);
  // std::cout << "Err in LU: " << (mat*chargeDensity.col(0) - rhs.col(0)).norm()/rhs.col(0).norm() << std::endl;
  // std::cout << "mat: " << std::endl << mat.block(0,0,10,10) << std::endl;
  // std::cout << "mat*: " << std::endl << mat.block(0,0,10,10)*chargeDensity.block(0,0,10,1) << std::endl;
  // std::cout << "rhs: " << std::endl << rhs.block(0,0,10,1) << std::endl;
  // std::cout << "chargeDensity: " << std::endl << chargeDensity.block(0,0,10,1) << std::endl;

}

double ElectrostaticsMatrix::getCapacitance(Eigen::VectorXd& chargeDensity) {
  Eigen::VectorXd areaVector(N);
  for (size_t i = 0; i < N; i++) {
    areaVector(i) = patchList[i].area;
  }
  return chargeDensity.dot(areaVector);
}

void ElectrostaticsMatrix::getRhs(int conductorIndex, Eigen::VectorXd& rhs) {
  rhs = Eigen::VectorXd::Zero(N);
  for (size_t i = 0; i < N; i++) {
    if (patchList[i].conductorIndex == conductorIndex) {
      rhs(i) = 1.0;
    }
  }
}

void ElectrostaticsMatrix::getCapacitanceMatrix(Eigen::MatrixXd chargeDensityMatrix, Eigen::MatrixXd& capacitanceMatrix) {
  Eigen::VectorXd areaVector(N);
  for (size_t i = 0; i < N; i++) {
    areaVector(i) = patchList[i].area;
  }
  for (size_t i = 0; i < numConductors; i++) { // index for chargedensity
    chargeDensityMatrix.col(i) = chargeDensityMatrix.col(i).cwiseProduct(areaVector);
  }
  capacitanceMatrix = Eigen::MatrixXd(numConductors, numConductors);
  // C_{ij}: charge on conductor i due to voltage on conductor j
  for (size_t i = 0; i < numConductors; i++) { // index for chargedensity
    for (size_t j = 0; j < numConductors; j++) { // index for rhs
      for (size_t k = 0; k < N; k++) {
        if (patchList[k].conductorIndex == i+1) {
          capacitanceMatrix(i,j) += chargeDensityMatrix(k,j);
        }
      }
    }
  }
}
