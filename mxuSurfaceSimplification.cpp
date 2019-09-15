#include "mxuSurfaceSimplification.h"
#include <igl/edges.h>
#include <igl/per_face_normals.h>
#include <iostream>



//typedef std::unordered_set<int> edge_set;
//
//Eigen::MatrixXi edges(const Eigen::MatrixXi& F)
//{
//
//
//  // Idea: use a map.
//  // - Assuming F is an |F| by 3 matrix, we know that there are 3 edges in a face
//
//  std::unordered_map<int, edge_set> edges;
//  int num_edges = 0;
//
//  for (int i = 0; i < F.rows(); i++) {
//
//    Eigen::RowVectorXi F_row = F.row(i);
//    assert(F_row.size() == 3);
//
//
//
//    // Traverse through the edges on the face, and add them into the edgeSet
//    for (int j = 0; j < 3; j++) {
//
//      // Need to check both directed edges, since we only want undirected edges
//      int v1 = F_row(j);
//      int v2 = F_row((j + 1) % 3);
//      if (edges.count(v1) > 0 || edges.count(v2) > 0) {
//        if (edges.count(v1) > 0 && edges[v1].count(v2) > 0) {
//          continue;
//        }
//        else if (edges.count(v2) > 0 && edges[v2].count(v1) > 0) {
//          continue;
//        }
//        else {
//          edges[v1].insert(v2);
//          num_edges++;
//        }
//      }
//      else {
//        edges[v1].insert(v2);
//        num_edges++;
//      }
//    }
//  }
//
//  Eigen::MatrixXi E(num_edges, 2);
//  int i = 0;
//  std::unordered_map<int, edge_set>::iterator edgeSet;
//  for (edgeSet = edges.begin(); edgeSet != edges.end(); edgeSet++) {
//    int v1 = (*edgeSet).first;
//    for (edge_set::iterator v2 = (*edgeSet).second.begin(); v2 != (*edgeSet).second.end(); v2++, i++) {
//
//      E(i, 0) = v1;
//      E(i, 1) = (*v2);
//    }
//  }
//
//  return E;
//}

void PrintEdges(const Eigen::MatrixXi& E) {
  for (int i = 0; i < E.rows(); i++) {
    std::cout << "Edge (" << std::to_string(i) << "): ";
    for (int j = 0; j < E.cols(); j++) {
      std::cout << std::to_string(E(i, j)) << ", ";
    }
    std::cout << std::endl;
  }
}

void removeRow(Eigen::MatrixXi& matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows() - 1;
  unsigned int numCols = matrix.cols();

  if (rowToRemove < numRows)
    matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

  matrix.conservativeResize(numRows, numCols);
}

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows() - 1;
  unsigned int numCols = matrix.cols();

  if (rowToRemove < numRows)
    matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

  matrix.conservativeResize(numRows, numCols);
}

int IterativeEdgeContraction(Eigen::MatrixXd& V, Eigen::MatrixXi& F, int targetNumFaces)
{
  Eigen::MatrixXi E;

  igl::edges(F, E);

  std::cout << "Initial edges: " << E.rows() << std::endl;

  // Get the normals of each face
  Eigen::MatrixXd N;
  igl::per_face_normals(V, F, N);

  // Get the constant d for the plane of each face
  Eigen::VectorXd d;
  d.resize(F.rows());
  for (int i = 0; i < F.rows(); i++) {
    Eigen::Vector3d vertexOnFace = V.row(F(i, 0));
    d(i) = -(N(i, 0) * vertexOnFace(0) + N(i, 1) * vertexOnFace(1) + N(i, 2) * vertexOnFace(2));
  }

  std::cout << "Done calculating the " << F.rows() << " values of d" << std::endl;

  // Now we have enough to calculate the quadric coefficients for each face => quadric coefficients for each vertex
  std::vector<Quadric> v_quadrics;
  v_quadrics.resize(V.rows());
  for (int i = 0; i < F.rows(); i++) {
    // first store the top right of symmetric matrix, then the diagonal
    Eigen::Matrix3d A = N.row(i).transpose() * N.row(i);
    Eigen::Vector3d b = N.row(i) * d(i);
    double c = d(i) * d(i);

    // assuming triangle mesh
    for (int j = 0; j < 3; j++) {
      int v = F(i, j);
      v_quadrics[v].A += A;
      v_quadrics[v].b += b;
      v_quadrics[v].c += c;
    }
    /*if (i == 1) {
      std::cout << "normal: " << N.row(i) << std::endl;
      std::cout << A << std::endl;
    }*/
  }
  std::cout << "Done calculating " << V.rows() << " quadric coefficients (1 quadric Q per vertex v)" << std::endl;


  // store 2 quadrics per edge, 1 for each choice of new v.
  // This is subset placement strategy.
  //std::vector<std::pair<Quadric, Quadric>> e_quadrics;

  double minError = DBL_MAX;
  int minErrorEdge = 0;
  Eigen::Vector3d V_min;

  for (int i = 0; i < E.rows(); i++) {
    //std::pair<Quadric, Quadric> qPair;
    Quadric q0, q1;
    int v0 = E(i, 0);
    int v1 = E(i, 1);

    q0 = v_quadrics[v0];
    q1 = v_quadrics[v1];

    Eigen::Vector3d V0 = V.row(v0);
    Eigen::Vector3d V1 = V.row(v1);

    Quadric q_total;
    q_total.A = q0.A + q1.A;
    q_total.b = q0.b + q1.b;
    q_total.c = q0.c + q1.c;

    // optimal placement
    if (/*q_total.A.determinant() != 0.0*/ false) {
      Eigen::Vector3d V_OPT = -q_total.A.inverse() * q_total.b;

      double error = -q_total.b.dot(V_OPT) + q_total.c;

      if (error < minError) {
        minError = error;
        minErrorEdge = i;
        V_min = V_OPT;
      }
    }
    else { // subset placement
      double error0 = q_total.Evaluate(V0);
      double error1 = q_total.Evaluate(V1);

      if (error0 < error1) {
        if (error0 < minError) {
          minError = error0;
          minErrorEdge = i;
          V_min = V.row(v0);
        }
      }
      else {
        if (error1 < minError) {
          minError = error1;
          minErrorEdge = i;
          V_min = V.row(v1);
        }
      }
    }
  }

  int v0 = E(minErrorEdge, 0);
  int v1 = E(minErrorEdge, 1);

  std::cout << "Found minErrorEdge: (" << v0 << ", " << v1 << ")";
  std::cout << " with error: " << minError << std::endl;

  // 1. move vertices v0 and v1 to V_min
  V.row(v0) = V_min;


  // 2. replace all occurences of v1 with v0

  // first replace all instances of v1 with v0.
  // then since v1 was removed from V, reduce the vertex of all vertices > v1 by 1
  // also remove all faces that have the same vertex (v0) twice (should be two)
  int f0 = 0;
  int f1 = 0;
  int sameVCounter = 0;
  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      if (F(i, j) == v1) {
        F(i, j) = v0;
      }
      if (F(i, j) > v1) {
        F(i, j) -= 1;
      }
    }
    // check if any repeated vertices (two keepVertex vertices)
    if (F(i, 0) == F(i, 1) || F(i, 1) == F(i, 2) || F(i, 2) == F(i, 0)) {
      sameVCounter++;
      if (sameVCounter == 1)
        f0 = i;
      else if (sameVCounter == 2)
        f1 = i;
    }
  }

  std::cout << "Found " << sameVCounter << " faces with repeated vertices." << std::endl;
  std::cout << "Faces: " << f0 << ", " << f1 << std::endl;

  // 3. remove v1 and any degenerate faces
  removeRow(F, f0);
  // f0 should be less than f1, so we will need to decrement f1
  f1--;
  removeRow(F, f1);
  std::cout << "Finished resizing F from " << F.rows() + 2 << " rows to " << F.rows() << " rows." << std::endl;

  removeRow(V, v1);
  std::cout << "Finished resizing V from " << V.rows() + 1 << " rows to " << V.rows() << " rows." << std::endl;


  return F.rows();
}
