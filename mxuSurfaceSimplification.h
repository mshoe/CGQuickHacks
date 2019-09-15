#pragma once
#include <Eigen/Core>
#include <vector>
#include <memory>
#include <unordered_map>
// From a list of triangle indices into some vertex set V, determine the
// set of unique undirected edges.
//
// Inputs:
//   F  #F by 3 list of triangles indices into some vertex set V
// Returns E  #E by 2 list of undirected edge indices into V
//Eigen::MatrixXi edges(const Eigen::MatrixXi& F);

/*
Will try to implement the algorithm in "Simplifying Surfaces with Color and Texture using Quadric Error Metrics"

For now, my goal is just to simplify the surface meshes, and ignore the color and textures.

** Compute the "cost" of each edge and store each edge in a data structure that can efficiently return the lowest-cost edge.
*** Edges are stored by the two vertices it connects (v1, v2)

** Each time we contract an edge, our face count should go down by 2

** User is given current # of faces
** User specifies target # of faces

while (currentFaceCount > targetFaceCount) {
  contractLowestCostEdge();
  RecomputeEdgeCosts();
  currentFaceCount -= 2;
}

*/



struct Quadric {
  Quadric() {
    A = Eigen::Matrix3d::Zero();
    b = Eigen::Vector3d::Zero();
    c = 0.0;
  }
  Quadric(Eigen::Matrix3d _A, Eigen::Vector3d _b, double _c) : A(_A), b(_b), c(_c) {}

  double Evaluate(const Eigen::Vector3d& vec) {
    return vec.dot(A * vec) + 2.0 * b.dot(vec) + c;
  }

  Eigen::Matrix3d A;
  Eigen::Vector3d b;
  double c;
};

struct Vertex {
  Quadric Q;

  // vertex position
  Eigen::Vector3d V;

  std::vector<int> f; // vertex can have multiple faces
  std::vector<int> e; // vertex can have multiple edges
};

struct Face {
  int v1, v2, v3;
};

struct Edge {
  int v1, v2; // edge connects two vertices

  // when an edge is contracted:
    // - v1 must replace v2 in every face/edge
    // - v1 must get updated to the chosen position (subset or optimal placement)

  double quadricError;
};

void removeRow(Eigen::MatrixXi& matrix, unsigned int rowToRemove);

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);


int IterativeEdgeContraction(Eigen::MatrixXd& V, Eigen::MatrixXi& F, int targetNumFaces);