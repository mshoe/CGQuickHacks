#pragma once
#include <Eigen/Core>
#include <vector>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <iterator>

#define MXU_DEBUG_OUTPUT 1
//#define MXU_FAST_SS_CHECKS 1

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
  currentFaceCount -= 2;//
}

*/

struct Quadric;
struct Vertex;
struct Face;
struct Edge;
struct VertexStructure;
struct FaceStructure;
struct EdgeStructure;


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

	int index; // need to set its index to nullptr in Vertex structure

	Quadric Q;

	// vertex position
	Eigen::Vector3d pos;

	std::unordered_set<std::shared_ptr<Edge>> edgeSet;
	std::unordered_set<std::shared_ptr<Face>> faceSet;

	void PrintEdgeSet();

};

struct Face {
	int index; // need to set its index to nullptr in Face structure

	std::shared_ptr<Vertex> v0 = nullptr, v1 = nullptr, v2 = nullptr;

	// plane equation of face: dot(n, v) + d = 0
	Eigen::Vector3d n; // normal n in plane equation
	double d; // d in plane equation

};

struct Edge {
	std::shared_ptr<Vertex> v0 = nullptr, v1 = nullptr;

	Eigen::Vector3d optimal_placement; // will be subset placement if optimal can't be computed
	double quadric_error;

	void ComputeErrorAndOptimalPlacement();
	bool edgeAlive = true;
	bool quadricErrorChanged = false;

	bool IsEquivalent(std::shared_ptr<Edge> otherEdge);
	void Print();

	//std::unordered_set<std::shared_ptr<Edge>> edges_to_be_removed; // excluding this edge from set, since this edge needs to be removed last
	//std::unordered_set<std::shared_ptr<Edge>> edges_to_be_updated;
	//std::unordered_set<std::shared_ptr<Face>> faces_to_be_removed;
	//std::unordered_set<std::shared_ptr<Face>> faces_to_be_updated;
};

void removeRow(Eigen::MatrixXi& matrix, unsigned int rowToRemove);

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);


int IterativeEdgeContraction(Eigen::MatrixXd& V, Eigen::MatrixXi& F, int targetNumFaces);

int MxuIterativeEdgeContraction(Eigen::MatrixXd& V, Eigen::MatrixXi& F, int targetNumFaces);

bool operator<(std::shared_ptr<Edge> a, std::shared_ptr<Edge> b);