#ifndef __PATCH_HPP__
#define __PATCH_HPP__

#include "includes.hpp"

class Patch {
public:
  int indexOfVertex1, indexOfVertex2, indexOfVertex3, conductorIndex;
  Node centroid;
  double area;
  double largestEdge;
  Patch(int indexOfVertex1, int indexOfVertex2, int indexOfVertex3);
};
#endif
