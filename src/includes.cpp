#include <iostream>
#include "includes.hpp"

double getDistance(Node& p1, Node& p2) {
  Node edge = p1-p2;
  return edge.norm();
}

double ElectrostaticGreensFunction(Node& ri, Node& rj) {
  return 1.0/getDistance(ri, rj);
}

cdouble EFIEGreensFunction(Node& ri, Node& rj) {
  double dist = getDistance(ri, rj);
  return (-1.0/4.0/PI)*std::exp(-I*kappa*dist)/dist;
}

cdouble nonSingularGreensFunction(Node& ri, Node& rj) {
  double dist = getDistance(ri, rj);
  return (-1.0/4.0/PI)*(std::exp(-I*kappa*dist)-1.0)/dist;
}

NodeC rhoGreensFunction(Node& ri, Node& rj, Node& r0) {
  double dist = getDistance(ri, rj);
  NodeC rho = rj;
  return rho*(-1.0/4.0/PI)*std::exp(-I*kappa*dist)/dist;
}

Node rhoLaplaceFunction(Node& ri, Node& rj, Node& r3) {
  double dist = getDistance(ri, rj);
  Node rho = rj - r3;
  // std::cout << "value1: " << rho/dist << std::endl;
  return rho*(-1.0/4.0/PI)/dist;
}

NodeC rhoNonSingularGreensFunction(Node& ri, Node& rj, Node& r3) {
  double dist = getDistance(ri, rj);
  Node rho = rj - r3;
  NodeC value;
  // std::cout << "rho: " << rho << std::endl;
  // std::cout << "rj: " << rj << std::endl;
  if (dist < 1e-12) {
    value = rho*(-1.0/4.0/PI)*(-I*kappa);
    // std::cout << "value1: " << -value*4.0*PI << std::endl;
  }
  else {
    value = rho*(-1.0/4.0/PI)*(std::exp(-I*kappa*dist)-1.0)/dist;
    // std::cout << "value2: " << -value*4.0*PI << std::endl;
  }
  return value;
}

NodeC farFieldIntegrand(Node& ri, Node& rj, Node& r3) {
  // double dist = getDistance(ri, rj);
  Node rho = rj - r3;
  NodeC value;
  value = rho*std::exp(I*kappa*rj.dot(ri));
  return value;
}
