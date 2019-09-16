#include "mxuSurfaceSimplification.h"
#include <igl/edges.h>
#include <igl/per_face_normals.h>
#include <iostream>
#include <igl/unique_edge_map.h>

#define MXU_DEBUG_OUTPUT 1

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
	// From profiling, recomputing all the values is what is most expensive.
	// Need to find a way to recompute values (normals, quadrics, etc.) without doing it for all elements
	Eigen::MatrixXi E;
	

	while (F.rows() > targetNumFaces) {

		igl::edges(F, E);
		

#ifdef MXU_DEBUG_OUTPUT
		std::cout << "Initial edges: " << E.rows() << std::endl;
#endif

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

#ifdef MXU_DEBUG_OUTPUT
		std::cout << "Done calculating the " << F.rows() << " values of d" << std::endl;
#endif

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
#ifdef MXU_DEBUG_OUTPUT
		std::cout << "Done calculating " << V.rows() << " quadric coefficients (1 quadric Q per vertex v)" << std::endl;
#endif

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

#ifdef MXU_DEBUG_OUTPUT
		std::cout << "Found minErrorEdge: (" << v0 << ", " << v1 << ")";
		std::cout << " with error: " << minError << std::endl;
#endif

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

#ifdef MXU_DEBUG_OUTPUT
		std::cout << "Found " << sameVCounter << " faces with repeated vertices." << std::endl;
		std::cout << "Faces: " << f0 << ", " << f1 << std::endl;
#endif

		// 3. remove v1 and any degenerate faces
		// Removal operations don't take that much time.
		removeRow(F, f0);
		// f0 should be less than f1, so we will need to decrement f1
		f1--;
		removeRow(F, f1);

#ifdef MXU_DEBUG_OUTPUT
		std::cout << "Finished resizing F from " << F.rows() + 2 << " rows to " << F.rows() << " rows." << std::endl;
#endif

		removeRow(V, v1);
#ifdef MXU_DEBUG_OUTPUT
		std::cout << "Finished resizing V from " << V.rows() + 1 << " rows to " << V.rows() << " rows." << std::endl;
#endif

		// Improvements:
		// 1. Can we remove the contracted edges from E now?
		// 2. Can we remove the normals corresponding to the removed faces in N now?
		// 3. Can we recalculate (only) the modified normals in N now?
		// 4. Can we recalculate (only) the modified quadrics in v_quadrics now?

		// Solutions:
		// 1. Replace v1 with v0 in all the edges
		//for (int i = 0; i < E.rows(); i++) {
		//	if (E(i, 0) == v1) {
		//		E(i, 0) = v0;
		//	}
		//	else if (E(i, 1) == v1) {
		//		E(i, 1) = v0;
		//	}
		//}

		//// 2.
		//removeRow(N, f0);
		//removeRow(N, f1);
	}

	return F.rows();
}

int MxuIterativeEdgeContraction(Eigen::MatrixXd& V, Eigen::MatrixXi& F, int targetNumFaces)
{
	// Function: move all vertex and face data to my own data structures to do this more efficiently


	

	// We want an vertex/face/edge data structure, where each edge will tell us how to modify F and E after an edge e: (v0, v1) is contracted.
		// If this edge is contracted:
			// v0.vec = v_OPT, v1.vec = v_OPT
		// which edges must be removed?
			// the edge (v0, v1)->(v0, v0) is this edge, and will be contracted ofc.
			// there will be two edge pairs like: (a, v0)->(a, v0), (a, v1)->(a, v0)
				// (a, v1) should be removed
		// which edges must be modified?
			// all edges of the form (a, v1)->(a, v0) or (v1, a)->(v0, a) that are unique (not in a pair after v1->v0)
		// which faces must be removed?
			// all (2) faces which contain a ptr to the edge e will be removed
		// which faces must be modified?
			// all faces that share an edge with the above 2 faces must update:
				// ptr v1->v0
				// recompute normals
				// recompute d in their plane equation
	
	// So when we constract an edge, we want to efficiently know:
		// set of ptrs to edges to be removed
		// set of ptrs to edges to be modified
		// set of ptrs to faces to be removed
		// set of ptrs to faces to be modified

#ifdef MXU_DEBUG_OUTPUT
	std::cout << "Initializing MxuIterativeEdgeContraction..." << std::endl;
#endif

	std::vector<std::shared_ptr<Vertex>> vertices;
	vertices.resize(V.rows());
	for (int i = 0; i < vertices.size(); i++) {
		vertices[i] = std::make_shared<Vertex>();
		vertices[i]->index = i;
		vertices[i]->pos = V.row(i);
	}

	std::vector<std::shared_ptr<Face>> faces;
	{
		// Get the normals of each face
		Eigen::MatrixXd N;
		igl::per_face_normals(V, F, N);

		faces.resize(F.rows());
		for (int i = 0; i < faces.size(); i++) {
			faces[i] = std::make_unique<Face>();
			faces[i]->index = i;
			faces[i]->v0 = vertices[F(i, 0)];
			faces[i]->v1 = vertices[F(i, 1)];
			faces[i]->v2 = vertices[F(i, 2)];

			faces[i]->v0->faceSet.insert(faces[i]);
			faces[i]->v1->faceSet.insert(faces[i]);
			faces[i]->v2->faceSet.insert(faces[i]);

			// Get the normal for the plane of the face
			faces[i]->n = N.row(i);

			// Get the constant d for the plane of the face
			Eigen::Vector3d vertexPosOnFace = faces[i]->v0->pos;
			faces[i]->d = -N.row(i).dot(vertexPosOnFace);

			// Calculate the quadric coefficients of the face
			Eigen::Matrix3d A = faces[i]->n * faces[i]->n.transpose();
			Eigen::Vector3d b = faces[i]->n * faces[i]->d;
			double c = faces[i]->d * faces[i]->d;

			// then add this faces quadrics to its vertices
			faces[i]->v0->Q.A += A;
			faces[i]->v0->Q.b += b;
			faces[i]->v0->Q.c += c;
			faces[i]->v1->Q.A += A;
			faces[i]->v1->Q.b += b;
			faces[i]->v1->Q.c += c;
			faces[i]->v2->Q.A += A;
			faces[i]->v2->Q.b += b;
			faces[i]->v2->Q.c += c;
		}
	}


	// First create all edges in a vector, then transfer them to a priority queue
	std::priority_queue<std::shared_ptr<Edge>> edgeQueue;
	{
		Eigen::MatrixXi E;
		igl::edges(F, E);

		// making sure it is working correctly or checking if the mesh is degenerate
		/*std::unordered_map<int, int> checkEdgeSet;
		int dups = 0;
		for (int i = 0; i < E.rows(); i++) {

			if (checkEdgeSet.count(E(i, 0)) && checkEdgeSet[E(i, 0)] == E(i, 1) ||
				checkEdgeSet.count(E(i, 1)) && checkEdgeSet[E(i, 1)] == E(i, 0)) {

				dups++;
				continue;
			}
			checkEdgeSet[E(i, 0)] = E(i, 1);
			checkEdgeSet[E(i, 1)] = E(i, 0);
		}
		std::cout << "DUPS: " << dups << std::endl;*/
		/*int badEdgeCounter = 0;
		for (int i = 0; i < E.rows(); i++) {
			if (E(i, 0) == E(i, 1))
				badEdgeCounter++;
		}
		std::cout << "BAD EDGES: " << badEdgeCounter << std::endl;*/

		std::vector<std::shared_ptr<Edge>> edges;
		edges.resize(E.rows());
		// first create all the edges
		for (int i = 0; i < edges.size(); i++) {
			edges[i] = std::make_shared<Edge>();
			edges[i]->v0 = vertices[E(i, 0)];
			edges[i]->v1 = vertices[E(i, 1)];

			edges[i]->v0->edgeSet.insert(edges[i]);
			edges[i]->v1->edgeSet.insert(edges[i]);
		}

		// Then evaluate all the edge quadrics
		for (int i = 0; i < edges.size(); i++) {
			edges[i]->ComputeErrorAndOptimalPlacement();
		}

		// now put all edges in the edge queue
		for (int i = 0; i < edges.size(); i++) {
			edgeQueue.push(edges[i]);
		}

		// CHECK IF EDGE QUEUE IS CORRECT
		//for (int i = 0; i < 100; i++) {
		//	std::shared_ptr<Edge> edge = edgeQueue.top();
		//	std::cout << "Edge: (" << edge->v0->index << ", " << edge->v1->index << ")" << std::endl;
		//	std::cout << "Quadric Error: " << edge->quadric_error << std::endl;
		//	edgeQueue.pop();
		//}
	}

#ifdef MXU_DEBUG_OUTPUT
	std::cout << "Initialization finished." << std::endl;
#endif

	int curNumVertices = vertices.size();
	int curNumFaces = faces.size();
	while (curNumFaces > targetNumFaces) {
		
		std::shared_ptr<Edge> minErrorEdge = edgeQueue.top();
		edgeQueue.pop();
		
		// Assigning a value of -2.0 to duplicate edges
		if (minErrorEdge->quadric_error < -1.0) {
			continue;
		}

		std::shared_ptr<Vertex> minv0, minv1;
		minv0 = minErrorEdge->v0;
		minv1 = minErrorEdge->v1;

#ifdef MXU_DEBUG_OUTPUT
		std::cout << std::endl << "Contracting edge: (" << minv0->index << ", " << minv1->index << ")" << std::endl;
		std::cout << "Quadric error: " << minErrorEdge->quadric_error << std::endl;
		std::cout << "Optimal placement: " << minErrorEdge->optimal_placement << std::endl;
#endif

		if (minv0 == minv1) {
			std::cout << "something went wrong" << std::endl;
		}

		// 1. v0 <- optimal placement
		minv0->pos = minErrorEdge->optimal_placement;

		// 2. replace v1 with v0 in all edges/faces

		// 2.1: faces
		// The only faces in v0's face set that have v1 are in v1's face set, so we can just iterate through v1's face set.
		for (std::shared_ptr<Face> v1face : minv1->faceSet) {
#ifdef MXU_DEBUG_OUTPUT
			if (v1face->v0 == v1face->v1 || v1face->v1 == v1face->v2 || v1face->v2 == v1face->v0) {
				std::cout << "error, degenerate face was found before replacing v1 with v0" << std::endl;
			}
#endif
			if (v1face->v0 == minv1) {
				v1face->v0 = minv0;
			}
			else if (v1face->v1 == minv1) {
				v1face->v1 = minv0;
			}
			else if (v1face->v2 == minv1) {
				v1face->v2 = minv0;
			}
		}

		
		// 2.2: edges
		// remove this edge from the edge sets
		minv0->edgeSet.erase(minErrorEdge);
		minv1->edgeSet.erase(minErrorEdge);

#ifdef MXU_DEBUG_OUTPUT
		std::cout << "Removed current edge from v0 and v1's edge sets." << std::endl;
#endif

		// The only edge in v0's edge set that has v1 is this current edge, so we can skip iterating through v0's edge set.
		// Iterate through v1's edge set and replace v1 with v0.
#ifdef MXU_DEBUG_OUTPUT
		for (std::shared_ptr<Edge> v0edge : minv0->edgeSet) {
			if (minv1->edgeSet.count(v0edge)) {
				std::cout << "error" << std::endl;
			}
		}
		int badEdgeCounter = 0;
		int changedEdgeCounter = 0;
#endif
		for (std::shared_ptr<Edge> searchEdge : minv1->edgeSet) {

#ifdef MXU_DEBUG_OUTPUT
			if (searchEdge->v0 == minv0 && searchEdge->v1 == minv1 ||
				searchEdge->v1 == minv0 && searchEdge->v0 == minv1) {
				std::cout << "how?" << std::endl;
			}
#endif

			if (searchEdge->v0 == minv1) {
				searchEdge->v0 = minv0;
#ifdef MXU_DEBUG_OUTPUT
				changedEdgeCounter++;
#endif
			}
			else if (searchEdge->v1 == minv1) {
				searchEdge->v1 = minv0;
#ifdef MXU_DEBUG_OUTPUT
				changedEdgeCounter++;
#endif
			}
#ifdef MXU_DEBUG_OUTPUT
			if (searchEdge->v0 == searchEdge->v1)
				badEdgeCounter++;
#endif
		}
#ifdef MXU_DEBUG_OUTPUT
		std::cout << "v1->edgeSet.size() == " << minv1->edgeSet.size() << std::endl;
		std::cout << "changed counter == " << changedEdgeCounter << std::endl;
		if (badEdgeCounter > 0) {
			std::cout << "Found " << badEdgeCounter << " bad edges in v1->edgeSet." << std::endl;
		}
#endif

		// 3. delete v1 and any degenerate faces/edges and any duplicate edges

		// 3.1: move edges from v1's edge set to v0
		// before deleting v1, move all of v1's edges to v0, except for duplicates.
		// in fact, mark duplicates to be skipped in the priority queue
		int duplicateEdgeCounter = 0;
		for (std::shared_ptr<Edge> v1Edge : minv1->edgeSet) {
			
			bool duplicate = false;
			for (std::shared_ptr<Edge> v0Edge : minv0->edgeSet) {

				bool duplicateCase1 = (v0Edge->v0 == v1Edge->v0 && v0Edge->v1 == v1Edge->v1);
				bool duplicateCase2 = (v0Edge->v1 == v1Edge->v0 && v0Edge->v0 == v1Edge->v1);

				if (duplicateCase1 || duplicateCase2) {
					duplicate = true;
					duplicateEdgeCounter++;
					break;
				}
			}

			if (!duplicate) {
				minv0->edgeSet.insert(v1Edge);
			}
			else {
				v1Edge->quadric_error = -2.0;
			}
		}
		
#ifdef MXU_DEBUG_OUTPUT
		std::cout << "Found " << duplicateEdgeCounter << " duplicate edges from from v1's edgeSet in v0's edgeSet." << std::endl;
#endif

		// 3.2.1: move faces from v1's face set to v0's face set.
		for (std::shared_ptr<Face> v1face : minv1->faceSet) {

			// ignore duplicates
			if (!minv0->faceSet.count(v1face)) {
				minv0->faceSet.insert(v1face);
			}
		}

		// 3.2.2: remove degenerate faces from v0's face set (should be 2)
		int degenFaces = 0;
		for (int i = 0; i < 3; i++) {
			for (std::shared_ptr<Face> v0face : minv0->faceSet) {
				if (v0face->v0 == v0face->v1 || v0face->v1 == v0face->v2 || v0face->v2 == v0face->v0) {
					// if the face is degenerate, remove all copies of its shared_ptr
					v0face->v0->faceSet.erase(v0face);
					v0face->v1->faceSet.erase(v0face);
					v0face->v2->faceSet.erase(v0face);
					faces[v0face->index] = nullptr;
					v0face = nullptr;
					degenFaces++;
					break;
				}
			}
		}
		if (degenFaces == 3) {
			std::cout << "is that possible?" << std::endl;
		}

#ifdef MXU_DEBUG_OUTPUT
		std::cout << "Removed " << degenFaces << " degenerate faces." << std::endl;
#endif

		curNumFaces -= degenFaces;

		// 3.3: finally remove v1 (its already been removed/replaced from all remaining faces/edges)
		vertices[minv1->index] = nullptr;
		minv1 = nullptr;
		curNumVertices -= 1;

		// 3.4: recompute plane equation on changed faces,
		// and recompute quadrics affected vertices
		for (std::shared_ptr<Face> v0face : minv0->faceSet) {

			// first remove this face's previous influence
			Eigen::Matrix3d A = v0face->n * v0face->n.transpose();
			Eigen::Vector3d b = v0face->n * v0face->d;
			double c = v0face->d * v0face->d;
			v0face->v0->Q.A -= A;
			v0face->v0->Q.b -= b;
			v0face->v0->Q.c -= c;
			v0face->v1->Q.A -= A;
			v0face->v1->Q.b -= b;
			v0face->v1->Q.c -= c;
			v0face->v2->Q.A -= A;
			v0face->v2->Q.b -= b;
			v0face->v2->Q.c -= c;
			
			// then add its new influence
			Eigen::Vector3d normal = (v0face->v1->pos - v0face->v0->pos).cross(v0face->v2->pos - v0face->v0->pos);
			double r = normal.norm();
			if (r == 0.0) {
				normal = Eigen::Vector3d::Zero();
			}
			else {
				normal /= r;
			}
			v0face->n = normal;
			v0face->d = -normal.dot(v0face->v0->pos);

			A = v0face->n * v0face->n.transpose();
			b = v0face->n * v0face->d;
			c = v0face->d * v0face->d;

			v0face->v0->Q.A += A;
			v0face->v0->Q.b += b;
			v0face->v0->Q.c += c;
			v0face->v1->Q.A += A;
			v0face->v1->Q.b += b;
			v0face->v1->Q.c += c;
			v0face->v2->Q.A += A;
			v0face->v2->Q.b += b;
			v0face->v2->Q.c += c;
		}

		// 3.6: recompute error and optimal placement on affected edges
		for (std::shared_ptr<Edge> v0edge : minv0->edgeSet) {
			v0edge->ComputeErrorAndOptimalPlacement();
		}

#ifdef MXU_DEBUG_OUTPUT
		std::cout << minv0->faceSet.size() << " Planes, ? quadrics, " << minv0->edgeSet.size() << " errors and optimal placements recomputed." << std::endl;
#endif

		
	}

	// Finally, reconstruct V and F.
#ifdef MXU_DEBUG_OUTPUT
	std::cout << "Initally " << V.rows() << " vertices and " << F.rows() << " faces." << std::endl;
	std::cout << "Now      " << curNumVertices << " vertices and " << curNumFaces << " faces." << std::endl;
#endif

	V = Eigen::MatrixXd();
	V.resize(curNumVertices, 3);

	F = Eigen::MatrixXi();
	F.resize(curNumFaces, 3);
	

	int j = 0;
	for (int i = 0; i < vertices.size(); i++) {
		if (vertices[i] == nullptr) {
			continue;
		}
		if (j >= curNumVertices) {
			std::cout << "something went wrong" << std::endl;
		}
		vertices[i]->index = j;
		V.row(j) = vertices[i]->pos;
		j++;
	}

	j = 0;
	for (int i = 0; i < faces.size(); i++) {
		if (faces[i] == nullptr) {
			continue;
		}
		F(j, 0) = faces[i]->v0->index;
		F(j, 1) = faces[i]->v1->index;
		F(j, 2) = faces[i]->v2->index;
		j++;
	}

	return curNumFaces;
}

// Using this for a min-heap/min-priority-queue of Edges
bool operator<(std::shared_ptr<Edge> a, std::shared_ptr<Edge> b)
{
	// shuffle nullptrs to the top of the queue
	if (b == nullptr) {
		return true;
	}
	if (a == nullptr) {
		return false;
	}

	return a->quadric_error > b->quadric_error;
}

void Edge::ComputeErrorAndOptimalPlacement()
{
	Quadric q_total;
	q_total.A = v0->Q.A + v1->Q.A;
	q_total.b = v0->Q.b + v1->Q.b;
	q_total.c = v0->Q.c + v1->Q.c;

	// optimal placement
	if (/*q_total.A.determinant() != 0.0*/ false) {
		optimal_placement = -q_total.A.inverse() * q_total.b;

		quadric_error = -q_total.b.dot(optimal_placement) + q_total.c;
	}
	else { // subset placement
		double error0 = q_total.Evaluate(v0->pos);
		double error1 = q_total.Evaluate(v1->pos);

		if (error0 < error1) {
			optimal_placement = v0->pos;
			quadric_error = error0;
		}
		else {
			optimal_placement = v1->pos;
			quadric_error = error1;
		}
	}
}

