#include <iostream>
#include <algorithm>
#include "EFIEMatrix.hpp"

EFIEMatrix::EFIEMatrix(inputsToKernelClass inputsToKernel) {
  this->getMatrixEntryFromScratchTime = 0;
  this->getMatrixEntryTime = 0;
  this->ngetMatrixEntryFromScratch = 0;
  this->ngetMatrixEntry = 0;
  this->etaAnalInteg = 4;
  std::string nodeListFilename = inputsToKernel.nodeListFilename;
  std::string patchListFilename = inputsToKernel.patchListFilename;
  this->n = inputsToKernel.n;
  this->L = inputsToKernel.L;
  generateNodeList(nodeListFilename);
  generatePatchList(patchListFilename);
  generateEdgeList();
  deleteNonRWGEdges();
  getSizeOfMatrix();
}

EFIEMatrix::EFIEMatrix(double Lx, double Ly, double Lz, int Nx, int Ny) {
  generateNodeList(Lx, Ly, Lz, Nx, Ny);
  generatePatchList(Nx, Ny);
  generateEdgeList();
  deleteNonRWGEdges();
  getSizeOfMatrix();
  // createLookUpTable();
}

void EFIEMatrix::generateNodeList(std::string fileName) {
  Node dummy = Node::Zero(3); // indexing in cpp starts at 0, though node indexes start from 1 in the patchlist
  nodeList.push_back(dummy); // so a dummy node is introduced at 0th position

  std::ifstream file(fileName);
  std::string line;
  while(std::getline(file, line)) {
    std::stringstream linestream(line);
    Node node;
    linestream >> node[0] >> node[1] >> node[2];
    nodeList.push_back(node);
  }
  file.close();
}

void EFIEMatrix::generateNodeList(double Lx, double Ly, double Lz, int Nx, int Ny) {
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

void EFIEMatrix::generatePatchList(std::string fileName) {
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

void EFIEMatrix::getSizeOfMatrix() {
  N = edge2PatchMap.size(); //assigning the size of the matrix
}

void EFIEMatrix::generatePatchList(int Nx, int Ny) {
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
      vertex1 = nodeList[indexOfVertex1];


      p.area = getAreaOfTriangle(vertex1, vertex2, vertex3);
      p.centroid = getCentroid(vertex1, vertex2, vertex3);
      p.largestEdge = getLargestEdge(vertex1, vertex2, vertex3);
      patchList.push_back(p);
    }
  }
}

void EFIEMatrix::generateEdgeList() {
  // create edgeList
  for (size_t i = 0; i < patchList.size(); i++) {
    Patch p = patchList[i];
    int index1 = p.indexOfVertex1;
    int index2 = p.indexOfVertex2;
    int index3 = p.indexOfVertex3;
    // std::pair<int, int> edge1, edge2, edge3;
    Edge edge1, edge2, edge3;
    ////////////////// edge 1 of the patch //////////////////
    if (index1 < index2) {
      edge1.indicesOfNodes = std::pair<int, int> (index1, index2);
    }
    else {
      edge1.indicesOfNodes = std::pair<int, int> (index2, index1);
    }
    edge1.length = getDistance(nodeList[index1], nodeList[index2]);
    edgeList.insert(edge1);
    if (edge2PatchMap.find(edge1.indicesOfNodes) == edge2PatchMap.end()) {
      edgeProperties eProperties1;
      eProperties1.length = getDistance(nodeList[index1], nodeList[index2]);
      eProperties1.indicesOfAssociatedPatches.push_back(i);
      edge2PatchMap[edge1.indicesOfNodes] = eProperties1;
    }
    else {
      edge2PatchMap[edge1.indicesOfNodes].indicesOfAssociatedPatches.push_back(i);
    }


    ////////////////// edge 2 of the patch //////////////////
    if (index2 < index3) {
      edge2.indicesOfNodes = std::pair<int, int> (index2, index3);
    }
    else {
      edge2.indicesOfNodes = std::pair<int, int> (index3, index2);
    }
    edge2.length = getDistance(nodeList[index2], nodeList[index3]);
    edgeList.insert(edge2);
    if (edge2PatchMap.find(edge2.indicesOfNodes) == edge2PatchMap.end()) {
      edgeProperties eProperties2;
      eProperties2.length = getDistance(nodeList[index2], nodeList[index3]);
      eProperties2.indicesOfAssociatedPatches.push_back(i);
      edge2PatchMap[edge2.indicesOfNodes] = eProperties2;
    }
    else {
      edge2PatchMap[edge2.indicesOfNodes].indicesOfAssociatedPatches.push_back(i);
    }

    ////////////////// edge 3 of the patch //////////////////
    if (index3 < index1) {
      edge3.indicesOfNodes = std::pair<int, int> (index3, index1);
    }
    else {
      edge3.indicesOfNodes = std::pair<int, int> (index1, index3);
    }
    edge3.length = getDistance(nodeList[index3], nodeList[index1]);
    edgeList.insert(edge3);
    if (edge2PatchMap.find(edge3.indicesOfNodes) == edge2PatchMap.end()) {
      edgeProperties eProperties3;
      eProperties3.length = getDistance(nodeList[index3], nodeList[index1]);
      eProperties3.indicesOfAssociatedPatches.push_back(i);
      edge2PatchMap[edge3.indicesOfNodes] = eProperties3;
    }
    else {
      edge2PatchMap[edge3.indicesOfNodes].indicesOfAssociatedPatches.push_back(i);
    }
    // edge2PatchMap[edge3.indicesOfNodes].push_back(i);
  }
}

void EFIEMatrix::deleteNonRWGEdges() {
  for (auto itr = edge2PatchMap.begin(); itr != edge2PatchMap.end();) {
    if(itr->second.indicesOfAssociatedPatches.size() == 1) {
      itr = edge2PatchMap.erase(itr);
      // itr = edgeList.erase(itr);
    }
    else {
      ++itr;
    }
  }
}

void EFIEMatrix::viewEdgeList() {
  std::cout << "No. of patches: " << patchList.size() << std::endl;
  std::cout << "No. of edges: " << edge2PatchMap.size() << std::endl;
  std::cout << "Size of map: " << edge2PatchMap.size() << std::endl;
  std::cout << "-------EDGES--------------------------- " << std::endl;
  int rwg = 0;
  for (auto itr : edge2PatchMap) {
    std::cout << itr.first.first << " " << itr.first.second << " ";
    if (itr.second.indicesOfAssociatedPatches.size() == 2) {
      rwg++;
    }
    std::cout << std::endl;
  }
  std::cout << "---------------------------------- " << std::endl;
  std::cout << "rwg: " << rwg << std::endl;
}

void EFIEMatrix::viewNodeList() {
  std::cout << "nodeList-------------------------" << std::endl;
  for (size_t i = 0; i < nodeList.size(); i++) {
    Node node = nodeList[i];
    std::cout << node[0] << " " << node[1] << " " << node[2] << std::endl;
  }
}

void EFIEMatrix::viewPatchList() {
  std::cout << "patchList-------------------------" << std::endl;
  for (size_t i = 0; i < patchList.size(); i++) {
    Patch p = patchList[i];
    int index1 = p.indexOfVertex1;
    int index2 = p.indexOfVertex2;
    int index3 = p.indexOfVertex3;
    std::cout << index1 << " " << index2 << " "<< index3 << std::endl;
  }
}

double EFIEMatrix::getLargestEdge(Node& vertex1, Node& vertex2, Node& vertex3) {
  double lengthOfEdge1 = getDistance(vertex1, vertex2);
  double lengthOfEdge2 = getDistance(vertex2, vertex3);
  double lengthOfEdge3 = getDistance(vertex3, vertex1);
  double lengthOfEdges[3] = {lengthOfEdge1, lengthOfEdge2, lengthOfEdge3};
  double* pointerOflargestEdge = std::max_element(lengthOfEdges, lengthOfEdges+2);
  return *pointerOflargestEdge;
}


Node EFIEMatrix::getCentroid(Node& vertex1, Node& vertex2, Node& vertex3) {
  Node centroid = (vertex1 + vertex2 + vertex3)/3.0;
  return centroid;
}

double EFIEMatrix::getAreaOfTriangle(Node& vertex1, Node& vertex2, Node& vertex3) {
  Node edge1 = vertex2 - vertex1;
  Node edge2 = vertex3 - vertex1;
  Node crossProduct = edge1.cross(edge2);
  double area = 0.5*crossProduct.norm();
  return area;
}

template<scalarFunctionPointer scalarFPtr>
cdouble EFIEMatrix::onePointGaussQuadratureRule(Node& centroidTrianglei, int indexOfSourceTriangle) {
  Patch patchSource = patchList[indexOfSourceTriangle];// jth patch
  return patchSource.area*scalarFPtr(centroidTrianglei, patchList[indexOfSourceTriangle].centroid);
}

template<vectorRealFunctionPointer vectorFPtr>
Node EFIEMatrix::onePointGaussQuadratureRule(Node& centroidTrianglei, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex) {
  Patch patchSource = patchList[indexOfSourceTriangle];// jth patch
    return patchSource.area*vectorFPtr(centroidTrianglei, patchList[indexOfSourceTriangle].centroid, nodeList[indexOfSourceNonRWGVertex]);
}

template<scalarFunctionPointer scalarFPtr>
cdouble EFIEMatrix::sevenPointGaussQuadratureRuleScalar(Node& centroidTrianglei, int indexOfSourceTriangle) {
  double alpha[7] = {0.3333333, 0.0597158, 0.470142, 0.470142, 0.7974269, 0.1012865, 0.1012865};
  double beta[7] = {0.3333333, 0.470142, 0.0597158, 0.470142, 0.1012865, 0.7974269, 0.1012865};
  double weights[7] = {0.225, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
  cdouble sum = 0.0;
  Patch patchSource = patchList[indexOfSourceTriangle];// jth patch
  int index1SourceTriangle = patchSource.indexOfVertex1;
  int index2SourceTriangle = patchSource.indexOfVertex2;
  int index3SourceTriangle = patchSource.indexOfVertex3;
  Node vertex1SourceTriangle = nodeList[index1SourceTriangle];
  Node vertex2SourceTriangle = nodeList[index2SourceTriangle];
  Node vertex3SourceTriangle = nodeList[index3SourceTriangle];
  for (size_t i = 0; i < 7; i++) {
    Node point = vertex1SourceTriangle + (vertex2SourceTriangle - vertex1SourceTriangle)*alpha[i] +(vertex3SourceTriangle-vertex1SourceTriangle)*beta[i];
    if (scalarFPtr == nonSingularGreensFunction) {
      if (getDistance(point, centroidTrianglei) < 1e-12) {
        sum += weights[i] * (I*kappa)/4.0/PI;
      }
      else {
        sum += weights[i] * scalarFPtr(centroidTrianglei, point);
      }
    }
    else {
      sum += weights[i] * scalarFPtr(centroidTrianglei, point);
    }
  }
  return patchSource.area * sum;
}

template<vectorComplexFunctionPointer vectorFPtr>
NodeC EFIEMatrix::sevenPointGaussQuadratureRuleComplex(Node& centroidTrianglei, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex) {
  double alpha[7] = {0.3333333, 0.0597158, 0.470142, 0.470142, 0.7974269, 0.1012865, 0.1012865};
  double beta[7] = {0.3333333, 0.470142, 0.0597158, 0.470142, 0.1012865, 0.7974269, 0.1012865};
  double weights[7] = {0.225, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
  NodeC sum = Node::Zero(3);
  Patch patchSource = patchList[indexOfSourceTriangle];// jth patch
  int index1SourceTriangle = patchSource.indexOfVertex2;
  int index2SourceTriangle = patchSource.indexOfVertex3;
  int index3SourceTriangle = patchSource.indexOfVertex1;
  Node vertex1SourceTriangle = nodeList[index1SourceTriangle];
  Node vertex2SourceTriangle = nodeList[index2SourceTriangle];
  Node vertex3SourceTriangle = nodeList[index3SourceTriangle];
  for (size_t i = 0; i < 7; i++) {
    Node point = vertex1SourceTriangle + (vertex2SourceTriangle - vertex1SourceTriangle)*alpha[i] +(vertex3SourceTriangle-vertex1SourceTriangle)*beta[i];
    sum += weights[i] * vectorFPtr(centroidTrianglei, point, nodeList[indexOfSourceNonRWGVertex]);
  }
  return patchSource.area * sum;
}

template<vectorRealFunctionPointer vectorFPtr>
Node EFIEMatrix::sevenPointGaussQuadratureRuleReal(Node& centroidTrianglei, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex) {
  double alpha[7] = {0.3333333, 0.0597158, 0.470142, 0.470142, 0.7974269, 0.1012865, 0.1012865};
  double beta[7] = {0.3333333, 0.470142, 0.0597158, 0.470142, 0.1012865, 0.7974269, 0.1012865};
  double weights[7] = {0.225, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
  Node sum = Node::Zero(3);
  Patch patchSource = patchList[indexOfSourceTriangle];// jth patch
  int index1SourceTriangle = patchSource.indexOfVertex2;
  int index2SourceTriangle = patchSource.indexOfVertex3;
  int index3SourceTriangle = patchSource.indexOfVertex1;
  Node vertex1SourceTriangle = nodeList[index1SourceTriangle];
  Node vertex2SourceTriangle = nodeList[index2SourceTriangle];
  Node vertex3SourceTriangle = nodeList[index3SourceTriangle];
  for (size_t i = 0; i < 7; i++) {
    Node point = vertex1SourceTriangle + (vertex2SourceTriangle - vertex1SourceTriangle)*alpha[i] +(vertex3SourceTriangle-vertex1SourceTriangle)*beta[i];
    sum += weights[i] * vectorFPtr(centroidTrianglei, point, nodeList[indexOfSourceNonRWGVertex]);
  }
  return patchSource.area * sum;
}


int EFIEMatrix::getIndexOfThirdNode(int indexOfTriangle, int indexOfTestEdgeNode1, int indexOfTestEdgeNode2) {
  // int *arrayOfNodesOfTriangle = new int[3];
  int arrayOfNodesOfTriangle[3];
  arrayOfNodesOfTriangle[0] = patchList[indexOfTriangle].indexOfVertex1;
  arrayOfNodesOfTriangle[1] = patchList[indexOfTriangle].indexOfVertex2;
  arrayOfNodesOfTriangle[2] = patchList[indexOfTriangle].indexOfVertex3;

  // int *arrayOfNodesOfEdge = new int[2];
  int arrayOfNodesOfEdge[2];
  arrayOfNodesOfEdge[0] = indexOfTestEdgeNode1;
  arrayOfNodesOfEdge[1] = indexOfTestEdgeNode2;

  std::sort(arrayOfNodesOfTriangle,arrayOfNodesOfTriangle+3);
  std::sort(arrayOfNodesOfEdge,arrayOfNodesOfEdge+2);

  std::vector<int> diff(1);
  std::set_difference(arrayOfNodesOfTriangle, arrayOfNodesOfTriangle+3, arrayOfNodesOfEdge, arrayOfNodesOfEdge+2, diff.begin());
  // delete arrayOfNodesOfTriangle;
  // delete arrayOfNodesOfEdge;
  return diff[0];
}

cdouble EFIEMatrix::getAIntegration(int indexOfTestTriangle, int indexOfSourceTriangle, int indexOfTestNonRWGVertex, int indexOfSourceNonRWGVertex) {
  // (     (
  // | p(r)| G(r,r')p(r')dr' dr
  // )     )
  // T     S
  // using 1 point rule for outer integration
  Node rhoCentroid = patchList[indexOfTestTriangle].centroid - nodeList[indexOfTestNonRWGVertex];
  // std::cout << "rhoCentroid: " << rhoCentroid << std::endl;
  return patchList[indexOfTestTriangle].area * rhoCentroid.transpose()*getVectorIntegrationGreensFunction(patchList[indexOfTestTriangle].centroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex);
}

Node EFIEMatrix::getVectorIntegrationSingularPart(Node& testCentroid, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex) {
  // (
  // | 1/r p(r,r') dr'
  // )
  // S

  double distanceBetweenCentroids = getDistance(testCentroid, patchList[indexOfSourceTriangle].centroid);
  if(distanceBetweenCentroids < etaAnalInteg*patchList[indexOfSourceTriangle].largestEdge) {
    Node testCentroidOnPlaneOfSourceTriangle = findProjectionOnTriangle(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex);
    // return (-1.0/4.0/PI)*(analyticIntegrationRhoLaplaceFunction(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex) + analyticIntegrationLaplaceFunction(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex) * (testCentroid - nodeList[indexOfSourceNonRWGVertex]));
    return (-1.0/4.0/PI)*(analyticIntegrationRhoLaplaceFunction(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex) + analyticIntegrationLaplaceFunction(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex) * (testCentroidOnPlaneOfSourceTriangle - nodeList[indexOfSourceNonRWGVertex]));
    // return (-1.0/4.0/PI)*(analyticIntegrationRhoLaplaceFunction(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex) + analyticIntegrationLaplaceFunction(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex) * (testCentroidOnPlaneOfSourceTriangle));
  }
  else if ((distanceBetweenCentroids >= etaAnalInteg*patchList[indexOfSourceTriangle].largestEdge) && (distanceBetweenCentroids < 10.0*patchList[indexOfSourceTriangle].largestEdge)) {
    return sevenPointGaussQuadratureRuleReal<rhoLaplaceFunction>(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex);
  }
  else {
    return onePointGaussQuadratureRule<rhoLaplaceFunction>(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex);
  }
}

NodeC EFIEMatrix::getVectorIntegrationNonSingularPart(Node& testCentroid, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex) {
  // (
  // | (e()-1)/r p(r')dr' dr'
  // )
  // S

  double distanceBetweenCentroids = getDistance(testCentroid, patchList[indexOfSourceTriangle].centroid);
  NodeC rho = patchList[indexOfSourceTriangle].centroid - nodeList[indexOfSourceNonRWGVertex];
  // if(distanceBetweenCentroids < 1e-12) {
  //   return (1.0/4.0/PI)*(I*kappa)*rho;
  // }
  // else if ((distanceBetweenCentroids >= 4.0*patchList[indexOfSourceTriangle].largestEdge) && (distanceBetweenCentroids < 10.0*patchList[indexOfSourceTriangle].largestEdge)) {
    return sevenPointGaussQuadratureRuleComplex<rhoNonSingularGreensFunction>(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex);
  // }
  // else {
  //   return onePointGaussQuadratureRule<rhoNonSingularGreensFunction>(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex);
  // }
}

// NodeC EFIEMatrix::getVectorIntegrationGreensFunction(Node& testCentroid, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex) {
//   // (                                  (                  (
//   // | G(r_{c},r')p(r')dr' =  -p_{1}(n3)| G(r_{c},r')dr' + |G(r_{c},r')p_{2}(r')dr'
//   // )                                  )                  )
//   // S                                  S                  S
//   // p(r') = -p_{1}(n3) + p_{2}(r')
//   // where n3 is the vertex corresponding to the nonRWG edges
//   return -nodeList[indexOfSourceNonRWGVertex]*getScalarIntegrationGreensFunction(testCentroid, indexOfSourceTriangle) + getVectorIntegration2GreensFunction(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex);
// }

NodeC EFIEMatrix::getVectorIntegrationGreensFunction(Node& testCentroid, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex) {
  // (                                  (                      (
  // | G(r_{c},r')p(r')dr' =            | (e()-1)/r p(r')dr' + |1/r p(r')dr'
  // )                                  )                      )
  // S                                  S                  S
  // p(r') = -p_{1}(n3) + p_{2}(r')
  // where n3 is the vertex corresponding to the nonRWG edges
  // return NodeC(getVectorIntegrationSingularPart(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex));//+NodeC::Zero();
  // return getVectorIntegrationNonSingularPart(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex);
  return getVectorIntegrationNonSingularPart(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex) + getVectorIntegrationSingularPart(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex);
}

cdouble EFIEMatrix::getPhiIntegration(int indexOfTestTriangle, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex) {
  // ( (
  // | | G(r,r')dr' dr
  // ) )
  // T S
  // using 1 point rule for outer integration
  return patchList[indexOfTestTriangle].area * getScalarIntegrationGreensFunction(patchList[indexOfTestTriangle].centroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex);
}

Node EFIEMatrix::findProjectionOnTriangle(Node& centroidTrianglei, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex) {
  // finds projection of observation point on to the source triangle
  int indexOfVertex1 = patchList[indexOfSourceTriangle].indexOfVertex1;
  int indexOfVertex2 = patchList[indexOfSourceTriangle].indexOfVertex2;
  int indexOfVertex3 = patchList[indexOfSourceTriangle].indexOfVertex3;
  Node vertex1Trianglej = nodeList[indexOfVertex1];
  Node vertex2Trianglej = nodeList[indexOfVertex2];
  Node vertex3Trianglej = nodeList[indexOfVertex3];

  Node n3 = nodeList[indexOfSourceNonRWGVertex];
  Node edge2Trianglej = vertex2Trianglej - vertex1Trianglej;
  Node edge3Trianglej = vertex3Trianglej - vertex1Trianglej;
  Node edge1Trianglej = centroidTrianglei;
  // Node edge1Trianglej = centroidTrianglei - n3;
  Node unitVector = edge2Trianglej.cross(edge3Trianglej);
  unitVector = -1.0*unitVector/unitVector.norm();
  double dotProduct = edge1Trianglej.dot(unitVector);
  Node centroid_on_plane = edge1Trianglej-dotProduct*unitVector;
  return centroid_on_plane;
}

double EFIEMatrix::analyticIntegrationLaplaceFunction(Node& centroidTrianglei, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex) {
  // making indexOfVertex1 free node (n3)
  int indexOfVertex1 = indexOfSourceNonRWGVertex;
  int indexOfVertex2;
  int indexOfVertex3;
  if (patchList[indexOfSourceTriangle].indexOfVertex1 == indexOfVertex1) {
    indexOfVertex2 = patchList[indexOfSourceTriangle].indexOfVertex2;
    indexOfVertex3 = patchList[indexOfSourceTriangle].indexOfVertex3;
  }
  else if (patchList[indexOfSourceTriangle].indexOfVertex2 == indexOfVertex1) {
    indexOfVertex2 = patchList[indexOfSourceTriangle].indexOfVertex1;
    indexOfVertex3 = patchList[indexOfSourceTriangle].indexOfVertex3;
  }
  else {
    indexOfVertex2 = patchList[indexOfSourceTriangle].indexOfVertex1;
    indexOfVertex3 = patchList[indexOfSourceTriangle].indexOfVertex2;
  }

  Node vertex1Trianglej = nodeList[indexOfVertex1];
  Node vertex2Trianglej = nodeList[indexOfVertex2];
  Node vertex3Trianglej = nodeList[indexOfVertex3];

  double areaOfTrianglej = getAreaOfTriangle(vertex1Trianglej, vertex2Trianglej, vertex3Trianglej);
  Node edge2Trianglej = vertex2Trianglej - vertex1Trianglej;
  Node edge3Trianglej = vertex3Trianglej - vertex1Trianglej;
  Node edge1Trianglej = centroidTrianglei - vertex1Trianglej;
  Node unitVector = edge2Trianglej.cross(edge3Trianglej);
  unitVector = -1.0*unitVector/unitVector.norm();
  double dotProduct = edge1Trianglej.dot(unitVector);
  Node centroid_on_plane = centroidTrianglei-dotProduct*unitVector;
  dotProduct = std::abs(dotProduct);
  Node origin(0,0,0);
  double MoM_entry = 0.0;

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
    // Find R_plus
    Eigen::Vector3d control = processNodes1 - centroidTrianglei;
    double R_plus = getDistance(control, origin);
    // Find R_minus
    control = processNodes2 - centroidTrianglei;
    double R_minus = getDistance(control, origin);
    // Find P_plus and P_minus
    Eigen::Vector3d P_plus = processNodes1 - centroid_on_plane;
    Eigen::Vector3d P_minus = processNodes2 - centroid_on_plane;
    // Find unit_l
    double distance=getDistance(processNodes1,processNodes2);
    Eigen::Vector3d unit_l = (processNodes1-processNodes2)/distance;
    // std::cout << "unit_l: " << std::endl << unit_l <<  std::endl;

    // Find L_plus and L_minus
    double L_plus = P_plus.dot(unit_l);
    double L_minus = P_minus.dot(unit_l);
    // Find unit_u
    Eigen::Vector3d unit_u = unit_l.cross(unitVector);
    // Find P0 and R0
    double P0 = P_plus.dot(unit_u);
    if (P0 < 0) {
      P0 = P0*-1;
    }
    if (P0 < 1.0e-10*length_side) {
          flag = 0;
    }
    // Find unit_p
    if(flag) {
      double R0 = std::sqrt((P0*P0)+(dotProduct*dotProduct));
      Eigen::Vector3d unit_p = P_plus-L_plus*unit_l;
  		distance = getDistance(unit_p, origin);
      unit_p = unit_p/distance;
      double T1 = P0*log((R_plus+L_plus)/(R_minus+L_minus));
      if (P0 < 1e-6*length_side) {
         if(R0 < 1e-4*length_side) {
            T1=0;
         }
      }
      // Calculate the part terms used in integration
      double T2 = std::atan(P0*L_plus/(R0*R0+dotProduct*R_plus));
      double T3 = std::atan(P0*L_minus/(R0*R0+dotProduct*R_minus));
      double T4 = unit_p.dot(unit_u);
      MoM_entry = MoM_entry+(T4*(T1-dotProduct*(T2-T3)));
    }
  }
  return MoM_entry;
}

Node EFIEMatrix::analyticIntegrationRhoLaplaceFunction(Node& centroidTrianglei, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex) {
  // making indexOfVertex1 free node (n3)
  int indexOfVertex1 = indexOfSourceNonRWGVertex;
  int indexOfVertex2;
  int indexOfVertex3;
  if (patchList[indexOfSourceTriangle].indexOfVertex1 == indexOfVertex1) {
    indexOfVertex2 = patchList[indexOfSourceTriangle].indexOfVertex2;
    indexOfVertex3 = patchList[indexOfSourceTriangle].indexOfVertex3;
  }
  else if (patchList[indexOfSourceTriangle].indexOfVertex2 == indexOfVertex1) {
    indexOfVertex2 = patchList[indexOfSourceTriangle].indexOfVertex1;
    indexOfVertex3 = patchList[indexOfSourceTriangle].indexOfVertex3;
  }
  else {
    indexOfVertex2 = patchList[indexOfSourceTriangle].indexOfVertex1;
    indexOfVertex3 = patchList[indexOfSourceTriangle].indexOfVertex2;
  }

  Node vertex1Trianglej = nodeList[indexOfVertex1];
  Node vertex2Trianglej = nodeList[indexOfVertex2];
  Node vertex3Trianglej = nodeList[indexOfVertex3];

  double areaOfTrianglej = getAreaOfTriangle(vertex1Trianglej, vertex2Trianglej, vertex3Trianglej);
  Node edge2Trianglej = vertex2Trianglej - vertex1Trianglej;
  Node edge3Trianglej = vertex3Trianglej - vertex1Trianglej;
  Node edge1Trianglej = centroidTrianglei - vertex1Trianglej;
  Node unitVector = edge2Trianglej.cross(edge3Trianglej);
  unitVector = -1.0*unitVector/unitVector.norm();
  double dotProduct = edge1Trianglej.dot(unitVector);
  Node centroid_on_plane = centroidTrianglei-dotProduct*unitVector;
  dotProduct = std::abs(dotProduct);
  Node origin(0,0,0);
  Node MoM_entry(0,0,0);

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
    // double length_side = getDistance(processNodes1, processNodes2);
    // Find R_plus
    Eigen::Vector3d control = processNodes1 - centroidTrianglei;
    double R_plus = getDistance(control, origin);
    // Find R_minus
    control = processNodes2 - centroidTrianglei;
    double R_minus = getDistance(control, origin);
    // Find P_plus and P_minus
    Eigen::Vector3d P_plus = processNodes1 - centroid_on_plane;
    Eigen::Vector3d P_minus = processNodes2 - centroid_on_plane;
    // Find unit_l
    double distance=getDistance(processNodes1,processNodes2);
    Eigen::Vector3d unit_l = (processNodes1-processNodes2)/distance;

    // Find L_plus and L_minus
    double L_plus = P_plus.dot(unit_l);
    double L_minus = P_minus.dot(unit_l);
    // Find unit_u
    Eigen::Vector3d unit_u = unit_l.cross(unitVector);
    // Find P0 and R0
    double P0 = P_plus.dot(unit_u);

    if (P0 < 0) {
      P0 = P0*-1;
    }
    double sq_R0 = (P0*P0)+(dotProduct*dotProduct);
    double T1;
    if((R_minus+L_minus)>=1.0e-10) {
      T1=sq_R0*log((R_plus+L_plus)/(R_minus+L_minus));
    }
    else {
      T1=0;
    }

     double T2=(L_plus*R_plus)-(L_minus*R_minus);
     MoM_entry=MoM_entry+(unit_u*0.5*(T1+T2));
  }
  return MoM_entry;
}

cdouble EFIEMatrix::getScalarIntegrationGreensFunction(Node& testCentroid, int indexOfSourceTriangle, int indexOfSourceNonRWGVertex) {
  // (
  // | G(r,r')dr' dr
  // )
  // S
  double distanceBetweenCentroids = getDistance(testCentroid, patchList[indexOfSourceTriangle].centroid);
  if(distanceBetweenCentroids < etaAnalInteg*patchList[indexOfSourceTriangle].largestEdge) {
    return (-1.0/4.0/PI)*analyticIntegrationLaplaceFunction(testCentroid, indexOfSourceTriangle, indexOfSourceNonRWGVertex) + sevenPointGaussQuadratureRuleScalar<nonSingularGreensFunction>(testCentroid, indexOfSourceTriangle);
  }
  else if ((distanceBetweenCentroids >= etaAnalInteg*patchList[indexOfSourceTriangle].largestEdge) && (distanceBetweenCentroids < 10.0*patchList[indexOfSourceTriangle].largestEdge)) {
    return sevenPointGaussQuadratureRuleScalar<EFIEGreensFunction>(testCentroid, indexOfSourceTriangle);
  }
  else {
    return onePointGaussQuadratureRule<EFIEGreensFunction>(testCentroid, indexOfSourceTriangle);
  }
}


// rename the below function as getMatrixEntry for non square plate simulations
// the below function generates matrix entries using the anlytical/numerical integrations
cdouble EFIEMatrix::getMatrixEntryFromScratch(const unsigned i, const unsigned j) {
  auto itrTestEdge = next(edge2PatchMap.begin(), i);
  auto itrSourceEdge = next(edge2PatchMap.begin(), j);

  // test and source edges
  std::pair<int, int> testEdge = itrTestEdge->first;
  std::pair<int, int> sourceEdge = itrSourceEdge->first;

  // indices of test positive and negative triangles
  int indexOfTestPositiveTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[0];
  int indexOfTestNegativeTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[1];

  // indices of source positive and negative triangles
  int indexOfSourcePositiveTriangle = edge2PatchMap[sourceEdge].indicesOfAssociatedPatches[0];
  int indexOfSourceNegativeTriangle = edge2PatchMap[sourceEdge].indicesOfAssociatedPatches[1];

  // indices of nodes that form test edge
  int indexOfTestEdgeNode1 = testEdge.first;
  int indexOfTestEdgeNode2 = testEdge.second;

  // indices of nodes that form source edge
  int indexOfSourceEdgeNode1 = sourceEdge.first;
  int indexOfSourceEdgeNode2 = sourceEdge.second;

  double Lt = itrTestEdge->second.length; // length of test edge
  double Ls = itrSourceEdge->second.length; // length of source edge

  // nodes of test positive triangle
  Node testPositiveTriangleNode1 = nodeList[indexOfTestEdgeNode1];
  Node testPositiveTriangleNode2 = nodeList[indexOfTestEdgeNode2];
  int testPositiveTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfTestPositiveTriangle, indexOfTestEdgeNode1, indexOfTestEdgeNode2);
  Node testPositiveTriangleNode3 = nodeList[testPositiveTriangleIndexOfNonRWGVertex];

  // nodes of test negative triangle
  Node testNegativeTriangleNode1 = nodeList[indexOfTestEdgeNode1];
  Node testNegativeTriangleNode2 = nodeList[indexOfTestEdgeNode2];
  int testNegativeTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfTestNegativeTriangle, indexOfTestEdgeNode1, indexOfTestEdgeNode2);
  Node testNegativeTriangleNode3 = nodeList[testNegativeTriangleIndexOfNonRWGVertex];

  // nodes of source positive triangle
  Node sourcePositiveTriangleNode1 = nodeList[indexOfSourceEdgeNode1];
  Node sourcePositiveTriangleNode2 = nodeList[indexOfSourceEdgeNode2];
  int sourcePositiveTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfSourcePositiveTriangle, indexOfSourceEdgeNode1, indexOfSourceEdgeNode2);
  Node sourcePositiveTriangleNode3 = nodeList[getIndexOfThirdNode(indexOfSourcePositiveTriangle, indexOfSourceEdgeNode1, indexOfSourceEdgeNode2)];

  // nodes of source negative triangle
  Node sourceNegativeTriangleNode1 = nodeList[indexOfSourceEdgeNode1];
  Node sourceNegativeTriangleNode2 = nodeList[indexOfSourceEdgeNode2];
  int sourceNegativeTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfSourceNegativeTriangle, indexOfSourceEdgeNode1, indexOfSourceEdgeNode2);
  Node sourceNegativeTriangleNode3 = nodeList[getIndexOfThirdNode(indexOfSourceNegativeTriangle, indexOfSourceEdgeNode1, indexOfSourceEdgeNode2)];

  // rhoCentroid of test positive and negative triangles; needed to perform one-point rule for outer integral
  Node rhoCentroidTestPositiveTriangle = patchList[indexOfTestPositiveTriangle].centroid - testPositiveTriangleNode3;
  Node rhoCentroidTestNegativeTriangle = patchList[indexOfTestNegativeTriangle].centroid - testNegativeTriangleNode3;

  // constants related to phi
  cdouble phiConstantPlusPlus = -Lt*Ls/(I*2.0*PI*freq*epsilon*patchList[indexOfTestPositiveTriangle].area*patchList[indexOfSourcePositiveTriangle].area);
  cdouble phiConstantPlusMinus = Lt*Ls/(I*2.0*PI*freq*epsilon*patchList[indexOfTestPositiveTriangle].area*patchList[indexOfSourceNegativeTriangle].area);
  cdouble phiConstantMinusPlus = Lt*Ls/(I*2.0*PI*freq*epsilon*patchList[indexOfTestNegativeTriangle].area*patchList[indexOfSourcePositiveTriangle].area);
  cdouble phiConstantMinusMinus = -Lt*Ls/(I*2.0*PI*freq*epsilon*patchList[indexOfTestNegativeTriangle].area*patchList[indexOfSourceNegativeTriangle].area);

  // matrix entries related to phi
  cdouble matrixEntryPhiPlusPlus = phiConstantPlusPlus * getPhiIntegration(indexOfTestPositiveTriangle, indexOfSourcePositiveTriangle, sourcePositiveTriangleIndexOfNonRWGVertex);
  cdouble matrixEntryPhiPlusMinus = phiConstantPlusMinus * getPhiIntegration(indexOfTestPositiveTriangle, indexOfSourceNegativeTriangle, sourceNegativeTriangleIndexOfNonRWGVertex);
  cdouble matrixEntryPhiMinusPlus = phiConstantMinusPlus * getPhiIntegration(indexOfTestNegativeTriangle, indexOfSourcePositiveTriangle, sourcePositiveTriangleIndexOfNonRWGVertex);
  cdouble matrixEntryPhiMinusMinus = phiConstantMinusMinus * getPhiIntegration(indexOfTestNegativeTriangle, indexOfSourceNegativeTriangle, sourceNegativeTriangleIndexOfNonRWGVertex);

  // constants related to A
  cdouble AConstantPlusPlus = -(I*2.0*PI*freq*mu*Lt*Ls)/(4.0*patchList[indexOfTestPositiveTriangle].area*patchList[indexOfSourcePositiveTriangle].area);
  cdouble AConstantPlusMinus = (I*2.0*PI*freq*mu*Lt*Ls)/(4.0*patchList[indexOfTestPositiveTriangle].area*patchList[indexOfSourceNegativeTriangle].area);
  cdouble AConstantMinusPlus = (I*2.0*PI*freq*mu*Lt*Ls)/(4.0*patchList[indexOfTestNegativeTriangle].area*patchList[indexOfSourcePositiveTriangle].area);
  cdouble AConstantMinusMinus = -(I*2.0*PI*freq*mu*Lt*Ls)/(4.0*patchList[indexOfTestNegativeTriangle].area*patchList[indexOfSourceNegativeTriangle].area);

  // matrix entries related to A
  cdouble matrixEntryAPlusPlus = AConstantPlusPlus * getAIntegration(indexOfTestPositiveTriangle, indexOfSourcePositiveTriangle, testPositiveTriangleIndexOfNonRWGVertex, sourcePositiveTriangleIndexOfNonRWGVertex);
  cdouble matrixEntryAPlusMinus = AConstantPlusMinus * getAIntegration(indexOfTestPositiveTriangle, indexOfSourceNegativeTriangle, testPositiveTriangleIndexOfNonRWGVertex, sourceNegativeTriangleIndexOfNonRWGVertex);
  cdouble matrixEntryAMinusPlus = AConstantMinusPlus * getAIntegration(indexOfTestNegativeTriangle, indexOfSourcePositiveTriangle, testNegativeTriangleIndexOfNonRWGVertex, sourcePositiveTriangleIndexOfNonRWGVertex);
  cdouble matrixEntryAMinusMinus = AConstantMinusMinus * getAIntegration(indexOfTestNegativeTriangle, indexOfSourceNegativeTriangle, testNegativeTriangleIndexOfNonRWGVertex, sourceNegativeTriangleIndexOfNonRWGVertex);

  //adding all of them
  cdouble matrixEntry = matrixEntryPhiPlusPlus + matrixEntryPhiPlusMinus + matrixEntryPhiMinusPlus + matrixEntryPhiMinusMinus;
  matrixEntry += matrixEntryAPlusPlus + matrixEntryAPlusMinus + matrixEntryAMinusPlus + matrixEntryAMinusMinus;

  return matrixEntry;
}

void EFIEMatrix::check() {
  std::cout << "edge2PatchMap" << std::endl << std::endl;
  for (size_t i = 0; i < N; i++) {
    auto itrTestEdge = next(edge2PatchMap.begin(), i);
    std::pair<int, int> testEdge = itrTestEdge->first;
    int indexOfTestPositiveTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[0];
    int indexOfTestNegativeTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[1];
    std::cout << indexOfTestPositiveTriangle << ", " << indexOfTestNegativeTriangle << std::endl;
  }
}

cdouble EFIEMatrix::getMatrixEntry(const unsigned i, const unsigned j) {
  // iterators of edges i and j in the edgeList
  auto itrTestEdge = next(edge2PatchMap.begin(), i);
  auto itrSourceEdge = next(edge2PatchMap.begin(), j);

  // test and source edges
  std::pair<int, int> testEdge = itrTestEdge->first;
  std::pair<int, int> sourceEdge = itrSourceEdge->first;

  // indices of test positive and negative triangles
  int indexOfTestPositiveTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[0];
  int indexOfTestNegativeTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[1];

  // indices of source positive and negative triangles
  int indexOfSourcePositiveTriangle = edge2PatchMap[sourceEdge].indicesOfAssociatedPatches[0];
  int indexOfSourceNegativeTriangle = edge2PatchMap[sourceEdge].indicesOfAssociatedPatches[1];

  // indices of nodes that form test edge
  int indexOfTestEdgeNode1 = testEdge.first;
  int indexOfTestEdgeNode2 = testEdge.second;

  // indices of nodes that form source edge
  int indexOfSourceEdgeNode1 = sourceEdge.first;
  int indexOfSourceEdgeNode2 = sourceEdge.second;

  double Lt = itrTestEdge->second.length; // length of test edge
  double Ls = itrSourceEdge->second.length; // length of source edge

  // nodes of test positive triangle
  Node testPositiveTriangleNode1 = nodeList[indexOfTestEdgeNode1];
  Node testPositiveTriangleNode2 = nodeList[indexOfTestEdgeNode2];
  int testPositiveTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfTestPositiveTriangle, indexOfTestEdgeNode1, indexOfTestEdgeNode2);
  Node testPositiveTriangleNode3 = nodeList[testPositiveTriangleIndexOfNonRWGVertex];

  // nodes of test negative triangle
  Node testNegativeTriangleNode1 = nodeList[indexOfTestEdgeNode1];
  Node testNegativeTriangleNode2 = nodeList[indexOfTestEdgeNode2];
  int testNegativeTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfTestNegativeTriangle, indexOfTestEdgeNode1, indexOfTestEdgeNode2);
  Node testNegativeTriangleNode3 = nodeList[testNegativeTriangleIndexOfNonRWGVertex];

  // nodes of source positive triangle
  Node sourcePositiveTriangleNode1 = nodeList[indexOfSourceEdgeNode1];
  Node sourcePositiveTriangleNode2 = nodeList[indexOfSourceEdgeNode2];
  int sourcePositiveTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfSourcePositiveTriangle, indexOfSourceEdgeNode1, indexOfSourceEdgeNode2);
  Node sourcePositiveTriangleNode3 = nodeList[getIndexOfThirdNode(indexOfSourcePositiveTriangle, indexOfSourceEdgeNode1, indexOfSourceEdgeNode2)];

  // nodes of source negative triangle
  Node sourceNegativeTriangleNode1 = nodeList[indexOfSourceEdgeNode1];
  Node sourceNegativeTriangleNode2 = nodeList[indexOfSourceEdgeNode2];
  int sourceNegativeTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfSourceNegativeTriangle, indexOfSourceEdgeNode1, indexOfSourceEdgeNode2);
  Node sourceNegativeTriangleNode3 = nodeList[getIndexOfThirdNode(indexOfSourceNegativeTriangle, indexOfSourceEdgeNode1, indexOfSourceEdgeNode2)];

  // rhoCentroid of test positive and negative triangles; needed to perform one-point rule for outer integral
  Node rhoCentroidTestPositiveTriangle = patchList[indexOfTestPositiveTriangle].centroid - testPositiveTriangleNode3;
  Node rhoCentroidTestNegativeTriangle = patchList[indexOfTestNegativeTriangle].centroid - testNegativeTriangleNode3;

  // constants related to phi
  cdouble phiConstantPlusPlus = -Lt*Ls/(I*2.0*PI*freq*epsilon*patchList[indexOfTestPositiveTriangle].area*patchList[indexOfSourcePositiveTriangle].area);
  cdouble phiConstantPlusMinus = Lt*Ls/(I*2.0*PI*freq*epsilon*patchList[indexOfTestPositiveTriangle].area*patchList[indexOfSourceNegativeTriangle].area);
  cdouble phiConstantMinusPlus = Lt*Ls/(I*2.0*PI*freq*epsilon*patchList[indexOfTestNegativeTriangle].area*patchList[indexOfSourcePositiveTriangle].area);
  cdouble phiConstantMinusMinus = -Lt*Ls/(I*2.0*PI*freq*epsilon*patchList[indexOfTestNegativeTriangle].area*patchList[indexOfSourceNegativeTriangle].area);

  // matrix entries related to phi
  cdouble matrixEntryPhiPlusPlus = phiConstantPlusPlus * getPhiIntegration(indexOfTestPositiveTriangle, indexOfSourcePositiveTriangle, sourcePositiveTriangleIndexOfNonRWGVertex);
  cdouble matrixEntryPhiPlusMinus = phiConstantPlusMinus * getPhiIntegration(indexOfTestPositiveTriangle, indexOfSourceNegativeTriangle, sourceNegativeTriangleIndexOfNonRWGVertex);
  cdouble matrixEntryPhiMinusPlus = phiConstantMinusPlus * getPhiIntegration(indexOfTestNegativeTriangle, indexOfSourcePositiveTriangle, sourcePositiveTriangleIndexOfNonRWGVertex);
  cdouble matrixEntryPhiMinusMinus = phiConstantMinusMinus * getPhiIntegration(indexOfTestNegativeTriangle, indexOfSourceNegativeTriangle, sourceNegativeTriangleIndexOfNonRWGVertex);

  // constants related to A
  cdouble AConstantPlusPlus = -(I*2.0*PI*freq*mu*Lt*Ls)/(4.0*patchList[indexOfTestPositiveTriangle].area*patchList[indexOfSourcePositiveTriangle].area);
  cdouble AConstantPlusMinus = (I*2.0*PI*freq*mu*Lt*Ls)/(4.0*patchList[indexOfTestPositiveTriangle].area*patchList[indexOfSourceNegativeTriangle].area);
  cdouble AConstantMinusPlus = (I*2.0*PI*freq*mu*Lt*Ls)/(4.0*patchList[indexOfTestNegativeTriangle].area*patchList[indexOfSourcePositiveTriangle].area);
  cdouble AConstantMinusMinus = -(I*2.0*PI*freq*mu*Lt*Ls)/(4.0*patchList[indexOfTestNegativeTriangle].area*patchList[indexOfSourceNegativeTriangle].area);

  // matrix entries related to A
  cdouble matrixEntryAPlusPlus = AConstantPlusPlus * getAIntegration(indexOfTestPositiveTriangle, indexOfSourcePositiveTriangle, testPositiveTriangleIndexOfNonRWGVertex, sourcePositiveTriangleIndexOfNonRWGVertex);
  cdouble matrixEntryAPlusMinus = AConstantPlusMinus * getAIntegration(indexOfTestPositiveTriangle, indexOfSourceNegativeTriangle, testPositiveTriangleIndexOfNonRWGVertex, sourceNegativeTriangleIndexOfNonRWGVertex);
  cdouble matrixEntryAMinusPlus = AConstantMinusPlus * getAIntegration(indexOfTestNegativeTriangle, indexOfSourcePositiveTriangle, testNegativeTriangleIndexOfNonRWGVertex, sourcePositiveTriangleIndexOfNonRWGVertex);
  cdouble matrixEntryAMinusMinus = AConstantMinusMinus * getAIntegration(indexOfTestNegativeTriangle, indexOfSourceNegativeTriangle, testNegativeTriangleIndexOfNonRWGVertex, sourceNegativeTriangleIndexOfNonRWGVertex);

  //adding all of them
  cdouble matrixEntry = matrixEntryPhiPlusPlus + matrixEntryPhiPlusMinus + matrixEntryPhiMinusPlus + matrixEntryPhiMinusMinus;
  matrixEntry += matrixEntryAPlusPlus + matrixEntryAPlusMinus + matrixEntryAMinusPlus + matrixEntryAMinusMinus;

  return matrixEntry;
}

void EFIEMatrix::getRhs(Node& Einc, Eigen::VectorXcd& rhs) {
  rhs = Eigen::VectorXcd(N);
  for (size_t i = 0; i < N; i++) {
    rhs(i) = getRhs(Einc, i);
  }
}

cdouble EFIEMatrix::getRhs(Node& Einc, int i) {
  auto itrTestEdge = next(edge2PatchMap.begin(), i);
  double Lt = itrTestEdge->second.length; // length of test edge
  std::pair<int, int> testEdge = itrTestEdge->first;
  int indexOfTestPositiveTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[0];
  int indexOfTestNegativeTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[1];
  int indexOfTestEdgeNode1 = testEdge.first;
  int indexOfTestEdgeNode2 = testEdge.second;

  // nodes of test positive triangle
  Node testPositiveTriangleNode1 = nodeList[indexOfTestEdgeNode1];
  Node testPositiveTriangleNode2 = nodeList[indexOfTestEdgeNode2];
  int testPositiveTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfTestPositiveTriangle, indexOfTestEdgeNode1, indexOfTestEdgeNode2);
  Node testPositiveTriangleNode3 = nodeList[testPositiveTriangleIndexOfNonRWGVertex];

  // nodes of test negative triangle
  Node testNegativeTriangleNode1 = nodeList[indexOfTestEdgeNode1];
  Node testNegativeTriangleNode2 = nodeList[indexOfTestEdgeNode2];
  int testNegativeTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfTestNegativeTriangle, indexOfTestEdgeNode1, indexOfTestEdgeNode2);
  Node testNegativeTriangleNode3 = nodeList[testNegativeTriangleIndexOfNonRWGVertex];

  // rhoCentroid of test positive and negative triangles; needed to perform one-point rule for outer integral
  Node rhoCentroidTestPositiveTriangle = patchList[indexOfTestPositiveTriangle].centroid - testPositiveTriangleNode3;
  Node rhoCentroidTestNegativeTriangle = patchList[indexOfTestNegativeTriangle].centroid - testNegativeTriangleNode3;

  double theta = 0.0;
  double phi = 0.0;
  Node k_v;
  k_v(0) = kappa*sin(theta)*cos(phi);
  k_v(1) = kappa*sin(theta)*sin(phi);
  k_v(2) = cos(theta); // Stores the Phase information
  double K_V_dot_pos_centroid = k_v.transpose() * patchList[indexOfTestPositiveTriangle].centroid;
  double K_V_dot_neg_centroid = k_v.transpose() * patchList[indexOfTestNegativeTriangle].centroid;
  NodeC Einc_pos = Einc*exp(I*K_V_dot_pos_centroid);
  NodeC Einc_neg = Einc*exp(I*K_V_dot_neg_centroid);
  // return Lt/2.0*(Einc_pos.transpose()*rhoCentroidTestPositiveTriangle - Einc_neg.transpose()*rhoCentroidTestNegativeTriangle);
  return Eigen::VectorXcd(Lt/2.0*(Einc_pos.transpose()*rhoCentroidTestPositiveTriangle - Einc_neg.transpose()*rhoCentroidTestNegativeTriangle))(0);
}

void EFIEMatrix::getCoefficients(Eigen::VectorXcd& rhs, Eigen::VectorXcd& coefficients) {
  Eigen::MatrixXcd mat;
  mat = getMatrix(0,0,N,N);

  Eigen::PartialPivLU<Eigen::MatrixXcd> matLU = mat.lu();
  coefficients = matLU.solve(rhs);
}

double EFIEMatrix::getRCS_SevenPointQuadrature_One(double theta, double phi, Eigen::VectorXcd& coefficients) {
  int neta = 377;
  std::complex<double> sigma_theta_s = 0;
  std::complex<double> sigma_phi_s = 0;
  Eigen::Vector3d k_v1;
  k_v1 << kappa*sin(theta)*cos(phi), kappa*sin(theta)*sin(phi), kappa*cos(theta);
  std::vector<Eigen::MatrixXcd> JAtSevenPointQudraturePoints;
  getJAtSevenPointQudraturePoints_One(coefficients, JAtSevenPointQudraturePoints);
  double alpha[1] = {0.3333333};
  double beta[1] = {0.3333333};
  double weights[1] = {1};
  for (size_t j = 0; j < patchList.size(); j++) {
    Eigen::MatrixXcd JAtSevenPointsOfPatch = JAtSevenPointQudraturePoints[j];
    std::complex<double> sigma_1 = 0;
    std::complex<double> sigma_2 = 0;
    Patch patchSource = patchList[j];// jth patch
    int index1SourceTriangle = patchSource.indexOfVertex1;
    int index2SourceTriangle = patchSource.indexOfVertex2;
    int index3SourceTriangle = patchSource.indexOfVertex3;
    Node vertex1SourceTriangle = nodeList[index1SourceTriangle];
    Node vertex2SourceTriangle = nodeList[index2SourceTriangle];
    Node vertex3SourceTriangle = nodeList[index3SourceTriangle];

    for (size_t i = 0; i < 1; i++) {
      Eigen::VectorXcd JAtPoint = JAtSevenPointsOfPatch.row(i);
      Node point = vertex1SourceTriangle + (vertex2SourceTriangle - vertex1SourceTriangle)*alpha[i] +(vertex3SourceTriangle-vertex1SourceTriangle)*beta[i];
      sigma_1 += weights[i] * (JAtPoint(0)*cos(theta)*cos(phi)+JAtPoint(1)*cos(theta)*sin(phi)-JAtPoint(2)*sin(theta))*(exp(I*k_v1.dot(point)));
      sigma_2 += weights[i] * (-JAtPoint(0)*sin(phi) + JAtPoint(1)*cos(phi))*exp(I*k_v1.dot(point));
    }
    sigma_theta_s += patchSource.area * sigma_1;
    sigma_phi_s += patchSource.area * sigma_2;
  }
  double sigma_theta_s_abs_square = double(kappa*kappa/(4.0*PI)*neta*std::abs(sigma_theta_s)*neta*std::abs(sigma_theta_s));
  double sigma_phi_s_abs_square = double(kappa*kappa/(4.0*PI)*neta*std::abs(sigma_phi_s)*neta*std::abs(sigma_phi_s));

  return sigma_theta_s_abs_square + sigma_phi_s_abs_square;
}


double EFIEMatrix::getRCS_SevenPointQuadrature(double theta, double phi, Eigen::VectorXcd& coefficients) {
  int neta = 377;
  std::complex<double> sigma_theta_s = 0;
  std::complex<double> sigma_phi_s = 0;
  Eigen::Vector3d k_v1;
  k_v1 << kappa*sin(theta)*cos(phi), kappa*sin(theta)*sin(phi), kappa*cos(theta);
  std::vector<Eigen::MatrixXcd> JAtSevenPointQudraturePoints;
  getJAtSevenPointQudraturePoints(coefficients, JAtSevenPointQudraturePoints);
  double alpha[7] = {0.3333333, 0.0597158, 0.470142, 0.470142, 0.7974269, 0.1012865, 0.1012865};
  double beta[7] = {0.3333333, 0.470142, 0.0597158, 0.470142, 0.1012865, 0.7974269, 0.1012865};
  double weights[7] = {0.225, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
  for (size_t j = 0; j < patchList.size(); j++) {
    Eigen::MatrixXcd JAtSevenPointsOfPatch = JAtSevenPointQudraturePoints[j];
    std::complex<double> sigma_1 = 0;
    std::complex<double> sigma_2 = 0;
    Patch patchSource = patchList[j];// jth patch
    int index1SourceTriangle = patchSource.indexOfVertex1;
    int index2SourceTriangle = patchSource.indexOfVertex2;
    int index3SourceTriangle = patchSource.indexOfVertex3;
    Node vertex1SourceTriangle = nodeList[index1SourceTriangle];
    Node vertex2SourceTriangle = nodeList[index2SourceTriangle];
    Node vertex3SourceTriangle = nodeList[index3SourceTriangle];

    for (size_t i = 0; i < 7; i++) {
      Eigen::VectorXcd JAtPoint = JAtSevenPointsOfPatch.row(i);
      Node point = vertex1SourceTriangle + (vertex2SourceTriangle - vertex1SourceTriangle)*alpha[i] +(vertex3SourceTriangle-vertex1SourceTriangle)*beta[i];
      sigma_1 += weights[i] * (JAtPoint(0)*cos(theta)*cos(phi)+JAtPoint(1)*cos(theta)*sin(phi)-JAtPoint(2)*sin(theta))*(exp(I*k_v1.dot(point)));
      sigma_2 += weights[i] * (-JAtPoint(0)*sin(phi) + JAtPoint(1)*cos(phi))*exp(I*k_v1.dot(point));
    }
    sigma_theta_s += patchSource.area * sigma_1;
    sigma_phi_s += patchSource.area * sigma_2;
  }
  double sigma_theta_s_abs_square = double(kappa*kappa/(4.0*PI)*neta*std::abs(sigma_theta_s)*neta*std::abs(sigma_theta_s));
  double sigma_phi_s_abs_square = double(kappa*kappa/(4.0*PI)*neta*std::abs(sigma_phi_s)*neta*std::abs(sigma_phi_s));

  return sigma_theta_s_abs_square + sigma_phi_s_abs_square;
}

double EFIEMatrix::getRCS_OnePointQuadrature(double theta, double phi, Eigen::MatrixXcd& JAtCentroids) {
  int neta = 377;
  std::complex<double> sigma_theta_s = 0;
  std::complex<double> sigma_phi_s = 0;
  Eigen::Vector3d k_v1;
  k_v1 << kappa*sin(theta)*cos(phi), kappa*sin(theta)*sin(phi), kappa*cos(theta);
  for (size_t j = 0; j < patchList.size(); j++) {
    sigma_theta_s += (JAtCentroids(j,0)*cos(theta)*cos(phi)+JAtCentroids(j,1)*cos(theta)*sin(phi)-JAtCentroids(j,2)*sin(theta))*(exp(I*k_v1.dot(patchList[j].centroid)))*patchList[j].area;
    sigma_phi_s += (-JAtCentroids(j,0)*sin(phi) + JAtCentroids(j,1)*cos(phi))*exp(I*k_v1.dot(patchList[j].centroid))*patchList[j].area;
  }
  double sigma_theta_s_abs_square = double(kappa*kappa/(4.0*PI)*neta*std::abs(sigma_theta_s)*neta*std::abs(sigma_theta_s));
  double sigma_phi_s_abs_square = double(kappa*kappa/(4.0*PI)*neta*std::abs(sigma_phi_s)*neta*std::abs(sigma_phi_s));

  return sigma_theta_s_abs_square + sigma_phi_s_abs_square;
}

double EFIEMatrix::getRCSMy(double theta, double phi, Eigen::VectorXcd& coefficients) {
  NodeC sum = NodeC::Zero(3);
  for (size_t i = 0; i < coefficients.size(); i++) {
    auto itrTestEdge = next(edge2PatchMap.begin(), i);
    std::pair<int, int> testEdge = itrTestEdge->first;
    int indexOfTestPositiveTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[0];
    int indexOfTestNegativeTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[1];
    int indexOfTestEdgeNode1 = testEdge.first;
    int indexOfTestEdgeNode2 = testEdge.second;
    double Lt = itrTestEdge->second.length; // length of test edge

    Node testPositiveTriangleNode1 = nodeList[indexOfTestEdgeNode1];
    Node testPositiveTriangleNode2 = nodeList[indexOfTestEdgeNode2];
    int testPositiveTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfTestPositiveTriangle, indexOfTestEdgeNode1, indexOfTestEdgeNode2);
    Node testPositiveTriangleNode3 = nodeList[testPositiveTriangleIndexOfNonRWGVertex];

    Node testNegativeTriangleNode1 = nodeList[indexOfTestEdgeNode1];
    Node testNegativeTriangleNode2 = nodeList[indexOfTestEdgeNode2];
    int testNegativeTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfTestNegativeTriangle, indexOfTestEdgeNode1, indexOfTestEdgeNode2);
    Node testNegativeTriangleNode3 = nodeList[testNegativeTriangleIndexOfNonRWGVertex];

    Node testPoint;
    testPoint << sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta);
    sum = sum + coefficients[i]*Lt/2.0/patchList[indexOfTestPositiveTriangle].area*sevenPointGaussQuadratureRuleComplex<farFieldIntegrand>(testPoint, indexOfTestPositiveTriangle, testPositiveTriangleIndexOfNonRWGVertex);
    sum = sum + coefficients[i]*Lt/2.0/patchList[indexOfTestNegativeTriangle].area*sevenPointGaussQuadratureRuleComplex<farFieldIntegrand>(testPoint, indexOfTestNegativeTriangle, testNegativeTriangleIndexOfNonRWGVertex);
  }
  return PI*freq*freq*mu*mu*sum.squaredNorm();
}

void EFIEMatrix::getJAtSevenPointQudraturePoints_One(Eigen::VectorXcd& coefficients, std::vector<Eigen::MatrixXcd>& JAtSevenPointQudraturePoints) {
  JAtSevenPointQudraturePoints.resize(patchList.size());
  for (size_t i = 0; i < patchList.size(); i++) {
    JAtSevenPointQudraturePoints[i] = Eigen::MatrixXcd::Zero(1,3);
  }
  // JAtCentroids = Eigen::MatrixXcd::Zero(patchList.size(), 3);
  for (size_t i = 0; i < coefficients.size(); i++) {
    auto itrTestEdge = next(edge2PatchMap.begin(), i);
    std::pair<int, int> testEdge = itrTestEdge->first;
    int indexOfTestPositiveTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[0];
    int indexOfTestNegativeTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[1];
    int indexOfTestEdgeNode1 = testEdge.first;
    int indexOfTestEdgeNode2 = testEdge.second;
    double Lt = itrTestEdge->second.length; // length of test edge

    Node testPositiveTriangleNode1 = nodeList[indexOfTestEdgeNode1];
    Node testPositiveTriangleNode2 = nodeList[indexOfTestEdgeNode2];
    int testPositiveTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfTestPositiveTriangle, indexOfTestEdgeNode1, indexOfTestEdgeNode2);
    Node testPositiveTriangleNode3 = nodeList[testPositiveTriangleIndexOfNonRWGVertex];

    Node testNegativeTriangleNode1 = nodeList[indexOfTestEdgeNode1];
    Node testNegativeTriangleNode2 = nodeList[indexOfTestEdgeNode2];
    int testNegativeTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfTestNegativeTriangle, indexOfTestEdgeNode1, indexOfTestEdgeNode2);
    Node testNegativeTriangleNode3 = nodeList[testNegativeTriangleIndexOfNonRWGVertex];

    double Ap = patchList[indexOfTestPositiveTriangle].area;
    double An = patchList[indexOfTestNegativeTriangle].area;
    // NodeC sum = NodeC::Zero(3);
    //////////////////////////////////////////////////////////////////
    double alpha[1] = {0.3333333};
    double beta[1] = {0.3333333};
    double weights[1] = {1.0};
    { // indexOfTestPositiveTriangle
      Patch patchSource = patchList[indexOfTestPositiveTriangle];// jth patch
      int index1SourceTriangle = patchSource.indexOfVertex1;
      int index2SourceTriangle = patchSource.indexOfVertex2;
      int index3SourceTriangle = patchSource.indexOfVertex3;
      Node vertex1SourceTriangle = nodeList[index1SourceTriangle];
      Node vertex2SourceTriangle = nodeList[index2SourceTriangle];
      Node vertex3SourceTriangle = nodeList[index3SourceTriangle];
      Eigen::MatrixXcd JAtSevenPointOfPatch = Eigen::MatrixXcd::Zero(1,3);

      for (size_t j = 0; j < 1; j++) {
        Node point = vertex1SourceTriangle + (vertex2SourceTriangle - vertex1SourceTriangle)*alpha[j] +(vertex3SourceTriangle-vertex1SourceTriangle)*beta[j];
        NodeC rhop = point - testPositiveTriangleNode3;
        JAtSevenPointOfPatch.row(j) = coefficients(i)*rhop*Lt/2/Ap;
      }
      JAtSevenPointQudraturePoints[indexOfTestPositiveTriangle] += JAtSevenPointOfPatch;
    }

    { // indexOfTestNegativeTriangle
      Patch patchSource = patchList[indexOfTestNegativeTriangle];// jth patch
      int index1SourceTriangle = patchSource.indexOfVertex1;
      int index2SourceTriangle = patchSource.indexOfVertex2;
      int index3SourceTriangle = patchSource.indexOfVertex3;
      Node vertex1SourceTriangle = nodeList[index1SourceTriangle];
      Node vertex2SourceTriangle = nodeList[index2SourceTriangle];
      Node vertex3SourceTriangle = nodeList[index3SourceTriangle];
      Eigen::MatrixXcd JAtSevenPointOfPatch = Eigen::MatrixXcd::Zero(1,3);
      for (size_t j = 0; j < 1; j++) {
        Node point = vertex1SourceTriangle + (vertex2SourceTriangle - vertex1SourceTriangle)*alpha[j] +(vertex3SourceTriangle-vertex1SourceTriangle)*beta[j];
        NodeC rhon = point - testNegativeTriangleNode3;
        JAtSevenPointOfPatch.row(j) = -coefficients(i)*rhon*Lt/2/An;
      }
      JAtSevenPointQudraturePoints[indexOfTestNegativeTriangle] += JAtSevenPointOfPatch;
    }
  }
}

void EFIEMatrix::getJAtSevenPointQudraturePoints(Eigen::VectorXcd& coefficients, std::vector<Eigen::MatrixXcd>& JAtSevenPointQudraturePoints) {
  JAtSevenPointQudraturePoints.resize(patchList.size());
  for (size_t i = 0; i < patchList.size(); i++) {
    JAtSevenPointQudraturePoints[i] = Eigen::MatrixXcd::Zero(7,3);
  }
  // JAtCentroids = Eigen::MatrixXcd::Zero(patchList.size(), 3);
  for (size_t i = 0; i < coefficients.size(); i++) {
    auto itrTestEdge = next(edge2PatchMap.begin(), i);
    std::pair<int, int> testEdge = itrTestEdge->first;
    int indexOfTestPositiveTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[0];
    int indexOfTestNegativeTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[1];
    int indexOfTestEdgeNode1 = testEdge.first;
    int indexOfTestEdgeNode2 = testEdge.second;
    double Lt = itrTestEdge->second.length; // length of test edge

    Node testPositiveTriangleNode1 = nodeList[indexOfTestEdgeNode1];
    Node testPositiveTriangleNode2 = nodeList[indexOfTestEdgeNode2];
    int testPositiveTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfTestPositiveTriangle, indexOfTestEdgeNode1, indexOfTestEdgeNode2);
    Node testPositiveTriangleNode3 = nodeList[testPositiveTriangleIndexOfNonRWGVertex];

    Node testNegativeTriangleNode1 = nodeList[indexOfTestEdgeNode1];
    Node testNegativeTriangleNode2 = nodeList[indexOfTestEdgeNode2];
    int testNegativeTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfTestNegativeTriangle, indexOfTestEdgeNode1, indexOfTestEdgeNode2);
    Node testNegativeTriangleNode3 = nodeList[testNegativeTriangleIndexOfNonRWGVertex];

    double Ap = patchList[indexOfTestPositiveTriangle].area;
    double An = patchList[indexOfTestNegativeTriangle].area;
    // NodeC sum = NodeC::Zero(3);
    //////////////////////////////////////////////////////////////////
    double alpha[7] = {0.3333333, 0.0597158, 0.470142, 0.470142, 0.7974269, 0.1012865, 0.1012865};
    double beta[7] = {0.3333333, 0.470142, 0.0597158, 0.470142, 0.1012865, 0.7974269, 0.1012865};
    double weights[7] = {0.225, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
    { // indexOfTestPositiveTriangle
      Patch patchSource = patchList[indexOfTestPositiveTriangle];// jth patch
      int index1SourceTriangle = patchSource.indexOfVertex1;
      int index2SourceTriangle = patchSource.indexOfVertex2;
      int index3SourceTriangle = patchSource.indexOfVertex3;
      Node vertex1SourceTriangle = nodeList[index1SourceTriangle];
      Node vertex2SourceTriangle = nodeList[index2SourceTriangle];
      Node vertex3SourceTriangle = nodeList[index3SourceTriangle];
      Eigen::MatrixXcd JAtSevenPointOfPatch = Eigen::MatrixXcd::Zero(7,3);

      for (size_t j = 0; j < 7; j++) {
        Node point = vertex1SourceTriangle + (vertex2SourceTriangle - vertex1SourceTriangle)*alpha[j] +(vertex3SourceTriangle-vertex1SourceTriangle)*beta[j];
        NodeC rhop = point - testPositiveTriangleNode3;
        JAtSevenPointOfPatch.row(j) = coefficients(i)*rhop*Lt/2/Ap;
      }
      JAtSevenPointQudraturePoints[indexOfTestPositiveTriangle] += JAtSevenPointOfPatch;
    }

    { // indexOfTestNegativeTriangle
      Patch patchSource = patchList[indexOfTestNegativeTriangle];// jth patch
      int index1SourceTriangle = patchSource.indexOfVertex1;
      int index2SourceTriangle = patchSource.indexOfVertex2;
      int index3SourceTriangle = patchSource.indexOfVertex3;
      Node vertex1SourceTriangle = nodeList[index1SourceTriangle];
      Node vertex2SourceTriangle = nodeList[index2SourceTriangle];
      Node vertex3SourceTriangle = nodeList[index3SourceTriangle];
      Eigen::MatrixXcd JAtSevenPointOfPatch = Eigen::MatrixXcd::Zero(7,3);
      for (size_t j = 0; j < 7; j++) {
        // std::cout << "j: " << j << std::endl;
        Node point = vertex1SourceTriangle + (vertex2SourceTriangle - vertex1SourceTriangle)*alpha[j] +(vertex3SourceTriangle-vertex1SourceTriangle)*beta[j];
        NodeC rhon = point - testNegativeTriangleNode3;
        JAtSevenPointOfPatch.row(j) = -coefficients(i)*rhon*Lt/2/An;
      }
      JAtSevenPointQudraturePoints[indexOfTestNegativeTriangle] += JAtSevenPointOfPatch;
    }
  }
}

void EFIEMatrix::getJAtCentroids(Eigen::VectorXcd& coefficients, Eigen::MatrixXcd& JAtCentroids) {
  JAtCentroids = Eigen::MatrixXcd::Zero(patchList.size(), 3);
  for (size_t i = 0; i < coefficients.size(); i++) {
    auto itrTestEdge = next(edge2PatchMap.begin(), i);
    std::pair<int, int> testEdge = itrTestEdge->first;
    int indexOfTestPositiveTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[0];
    int indexOfTestNegativeTriangle = edge2PatchMap[testEdge].indicesOfAssociatedPatches[1];
    int indexOfTestEdgeNode1 = testEdge.first;
    int indexOfTestEdgeNode2 = testEdge.second;
    double Lt = itrTestEdge->second.length; // length of test edge

    Node testPositiveTriangleNode1 = nodeList[indexOfTestEdgeNode1];
    Node testPositiveTriangleNode2 = nodeList[indexOfTestEdgeNode2];
    int testPositiveTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfTestPositiveTriangle, indexOfTestEdgeNode1, indexOfTestEdgeNode2);
    Node testPositiveTriangleNode3 = nodeList[testPositiveTriangleIndexOfNonRWGVertex];

    Node testNegativeTriangleNode1 = nodeList[indexOfTestEdgeNode1];
    Node testNegativeTriangleNode2 = nodeList[indexOfTestEdgeNode2];
    int testNegativeTriangleIndexOfNonRWGVertex = getIndexOfThirdNode(indexOfTestNegativeTriangle, indexOfTestEdgeNode1, indexOfTestEdgeNode2);
    Node testNegativeTriangleNode3 = nodeList[testNegativeTriangleIndexOfNonRWGVertex];

    double Ap = patchList[indexOfTestPositiveTriangle].area;
    double An = patchList[indexOfTestNegativeTriangle].area;

    NodeC rhop = patchList[indexOfTestPositiveTriangle].centroid - testPositiveTriangleNode3;
    NodeC rhon = patchList[indexOfTestNegativeTriangle].centroid - testNegativeTriangleNode3;

    JAtCentroids.row(indexOfTestPositiveTriangle) += coefficients(i)*rhop*Lt/2/Ap;
    JAtCentroids.row(indexOfTestNegativeTriangle) += -coefficients(i)*rhon*Lt/2/An;
  }
}
