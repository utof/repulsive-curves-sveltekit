Directory Structure:

└── ./
    ├── include
    │   ├── build_A_bar.h
    │   ├── build_weights.h
    │   ├── constraint_derivative.h
    │   ├── init_curve.h
    │   ├── loss_derivative.h
    │   ├── read_curve_file.h
    │   └── repulsion_iteration.h
    ├── src
    │   ├── build_A_bar.cpp
    │   ├── build_weights.cpp
    │   ├── constraint_derivative.cpp
    │   ├── init_curve.cpp
    │   ├── loss_derivative.cpp
    │   ├── read_curve_file.cpp
    │   └── repulsion_iteration.cpp
    └── main.cpp



---
File: /include/build_A_bar.h
---

#ifndef BUILD_A_BAR_H
#define BUILD_A_BAR_H
#include <Eigen/Core>
#include <vector>
// Build the Sobolev inner product matrix A_bar
//
// A_bar satisfies the equation:
// A_bar*g = (the derivative of the tangent point energy with respect to vertex positions)
// where g is the discrete fractional Sobolev gradient
//
// Inputs:
//   pt_num  the number of vertices in the curve (equal to #V)
//   E  #E by 2 list of edge indices for the curve
//   Ac  a vector of vectors containing a vector for each edge "I". This vector consists of
//	   the edge indices for all edges which do not intersect edge "I"
//	 T  #E by 3 list of unit tangent vectors for each edge of the curve. 
//	   Namely, T.row(e) = (V.row(E(e, 1)) - v.row(E(e, 0)))/L(e) where L is the next input
//	 L  #E by 1 list of the length of each edge in the curve
//   W  #E by #E list of weights used to construct the B matrix
//   W0  #E by #E list of weights used to construct the B0 matrix
// Outputs:
//   A_bar  (3*#V) by (3*#V) matrix corresponding to the fractional Sobolev inner product
//	  for a given curve, and given energy.
//   
void build_A_bar(
	const int& pt_num,
	const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& Ac,
	const Eigen::MatrixX3d& T,
	const Eigen::VectorXd& L,
	const Eigen::MatrixXd& W,
	const Eigen::MatrixXd& W0,
	Eigen::MatrixXd& A_bar);
#endif




---
File: /include/build_weights.h
---

#ifndef BUILD_WEIGHTS_H
#define BUILD_WEIGHTS_H
#include <Eigen/Core>
#include <vector>
// Build two weight matrices as described in the paper. It will be used in "build_A_bar" to build
// a matrix A_bar which is used for computing Sobolev inner products.
//
// W corresponds to the discrete high-order fractional Sobolev inner product
// W0 corresponds to the discrete low-order term (which will build a matrix B0 which will act
// like a regularizer for the more important matrix B). B and B0 will be used to construct A_bar.
// A_bar is described in "build_A_bar.h"
//
// Inputs:
//   alpha  double which changes the nature of the tangent point energy
//   beta  double which changes the nature of the tangent point energy like alpha
//   E  #E by 2 list of edge indices for the curve
//   Ac  a vector of vectors containing a vector for each edge "I". This vector consists of
//	   the edge indices for all edges which do not intersect edge "I"
//   V  #V by 3 list of curve vertex positions
//	 T  #E by 3 list of unit tangent vectors for each edge of the curve. 
//	   Namely, T.row(e) = (V.row(E(e, 1)) - v.row(E(e, 0)))/L(e) where L is the next input
//	 L  #E by 1 list of the length of each edge in the curve
// Outputs:
//   W  #E by #E list of weights used to construct the B matrix in build_A_bar.cpp
//   W0  #E by #E list of weights used to construct the B0 matrix in build_A_bar.cpp
void build_weights(
	const double& alpha,
	const double& beta,
	const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& Ac,
	const Eigen::MatrixX3d& V,
	const Eigen::MatrixX3d& T,
	const Eigen::VectorXd& L,
	Eigen::MatrixXd& W,
	Eigen::MatrixXd& W0);
#endif



---
File: /include/constraint_derivative.h
---

#ifndef CONSTRAINT_DERIVATIVE_H
#define CONSTRAINT_DERIVATIVE_H
#include <Eigen/Core>
#include <vector>

// Compute the derivative of each constraint with respect to vert positions.
// The constraints used are the total length constraint and the 
// individual edge length constraint. The formulas for these constraints are on
// page 13 of the paper, but I computed the derivative formulas by hand and then
// confirmed the result with Mathematica).
//
// Inputs:
//   E  #E by 2 list of edge indices for the curve
//   V  #V by 3 list of curve vertex positions
//   L  #E by 1 list of the length of each edge in the curve
//   E_adj  a vector of vectors where there is a vector for each point "p" consisting of 
//	   the edge indices for all edges containing the point "p"
//   constr_edges  list of the indices of the edges you want to constrain the lengths of.
//	   Note that this is different from constraining the total length of all edges. In the
//	   code, it is automatically set so that the edges which are in this list are the ones
//	   connected to a verticy of index 3 or more.
// Outputs:
//   C  (1 + #constr_edges) by 3 * #V list of the derivative of two constraints: the total length
//	   constraint, and the individual edge length constraints. The first row is the derivative
//	   of the total length constraint, and each following row corresponds to the derivative of
//	   the fixed edge length constraint for a different edge in "constr_edges"
void constraint_derivative(
	const Eigen::MatrixXi& E,
	const Eigen::MatrixX3d& V,
	const Eigen::VectorXd& L,
	const std::vector<std::vector<int>>& E_adj,
	const Eigen::VectorXi& constr_edges,
	Eigen::MatrixXd& C);
#endif



---
File: /include/init_curve.h
---

#ifndef INIT_CURVE_H
#define INIT_CURVE_H

#include <Eigen/Core>
#include <vector>
// Compute two vector of vectors (Ac and E_adj) which only depend
// on the curve's topology and thus won't have to be recomputed
// for each iteration of descent. 
// 
// Having these will save time on every iteration.
//
// Also, compute the length of every edge in the curve's initial state ("lengths").
// This will be used later when evaluating constraints.
//
// Inputs:
//   E  #E by 2 list of edge indices for the curve
//   V  #V by 3 list of curve vertex positions
// Outputs:
//   Ac  a vector of vectors containing a vector for each edge "I". This vector consists of
//	   the edge indices for all edges which do not intersect edge "I"
//   E_adj  a vector of vectors where there is a vector for each point "p" consisting of 
//	   the edge indices for all edges containing the point "p"
//   lengths  #E by 1 list of the length of each edge in the curve before any iterations were done
void init_curve(
	const Eigen::MatrixXi& E,
	const Eigen::MatrixX3d& V,
	std::vector<std::vector<int>>& Ac,
	std::vector<std::vector<int>>& E_adj,
	Eigen::VectorXd& lengths);
#endif



---
File: /include/loss_derivative.h
---

#ifndef LOSS_DERIVATIVE_H
#define LOSS_DERIVATIVE_H
#include <Eigen/Core>
#include <vector>
// Compute the derivative of the tangent point energy with respect to each vertex
// position, evaluated at the current vertex positions
//
// Inputs:
//   alpha  double which changes the nature of the tangent point energy
//   beta  double which changes the nature of the tangent point energy like alpha
//   E  #E by 2 list of edge indices for the curve
//   Ac  a vector of vectors containing a vector for each edge "I". This vector consists of
//	   the edge indices for all edges which do not intersect edge "I"
//   E_adj  a vector of vectors where there is a vector for each point "p" consisting of 
//	   the edge indices for all edges containing the point "p"
//   V  #V by 3 list of curve vertex positions
//	 T  #E by 3 list of unit tangent vectors for each edge of the curve. 
//	   Namely, T.row(e) = (V.row(E(e, 1)) - v.row(E(e, 0)))/L(e) where L is the next input
//	 L  #E by 1 list of the length of each edge in the curve
// Outputs:
//   Deriv  1 by #V * 3 list of the derivative of the tangent point energy with respect
//	   to each vertex position, evaluated at the current vertex positions
void loss_derivative(
	const double& alpha,
	const double& beta,
	const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& Ac,
	const std::vector<std::vector<int>>& E_adj,
	const Eigen::MatrixX3d& V,
	const Eigen::MatrixX3d& T,
	const Eigen::VectorXd& L,
	Eigen::RowVectorXd& Deriv);

#endif




---
File: /include/read_curve_file.h
---

#ifndef READ_CURVE_FILE_H
#define READ_CURVE_FILE_H
#include <Eigen/Core>
#include <string>

// Given a file path to an OBJ file containing a curve,
// find the vertex position and edges of this curve.
//
// THE FILE MUST BE AN OBJ FILE --- not another file format
//
// Inputs:
//   input  the file path to an OBJ file containing the curve to be used
// Outputs:
//   V  #V by 3 list of 3D positions of each vertex
//   E  #E by 2 list of indices for the virtices of each edge
//   
void read_curve_file(
    const std::string& input,
    Eigen::MatrixX3d& V,
    Eigen::MatrixXi& E);
#endif





---
File: /include/repulsion_iteration.h
---

#ifndef REPULSION_ITERATION_H
#define REPULSION_ITERATION_H

#include <Eigen/Core>
#include <vector>

// Conduct a single iteration of fractional Sobolev descent on a curve
// to minimize the tangent-point energy
//
// The first 6 inputs to this function are explained in much more
// detail in the main.cpp file, where they are meant to be tuned
//
// Inputs:
//   alpha  double which changes the nature of the tangent point energy
//   beta  double which changes the nature of the tangent point energy like alpha
//   a_const  parameter for my implementation of backtracking line search (must be between 0 and 0.5)
//	   used as the constant for the Armijo condition 
//   b_const  parameter for my implementation of backtracking line search (must be between 0 and 1)
//	   where in each step of line search, the step size gets multiplied by b_const
//   threshold  enforces that after the descent step, the norm of the constraints at the new vertex
//	   positions is less than this threshold
//   max_iters  max number of iterations to use when iteratively projecting 
//   E  #E by 2 list of edge indices for the curve
//   Ac  a vector of vectors containing a vector for each edge "I". This vector consists of
//	   the edge indices for all edges which do not intersect edge "I"
//   E_adj  a vector of vectors where there is a vector for each point "p" consisting of 
//	   the edge indices for all edges containing the point "p"
//   lengths  #E by 1 list of the length of each edge in the curve before any iterations were done
// Outputs:
//   V  #V by 3 list of curve vertex positions after the descent step 
void repulsion_iteration(
	const double& alpha,
	const double& beta,
	const double& a_const,
	const double& b_const,
	const double& threshold,
	const int& max_iters,
	const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& Ac,
	const std::vector<std::vector<int>>& E_adj,
	const Eigen::VectorXd& lengths,
	Eigen::MatrixX3d& V);

#endif



---
File: /src/build_A_bar.cpp
---

#include "build_A_bar.h"

#include <Eigen/Core>
#include <cmath>
#include <vector>

void build_A_bar(
	const int& pt_num,
	const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& Ac,
	const Eigen::MatrixX3d& T,
	const Eigen::VectorXd& L,
	const Eigen::MatrixXd& W,
	const Eigen::MatrixXd& W0,
	Eigen::MatrixXd& A_bar)
{

	Eigen::MatrixXd B;
	Eigen::MatrixXd B0;
	// Now we will build B and B0

	// Start with B and B0 as zero matrices
	B.setZero(pt_num, pt_num);
	B0.setZero(pt_num, pt_num);

	// B is the inner product matrix for the sobolev inner product
	// B0 acts like a regularizer (it's a lower order derivative) which makes the sum A = B + B0 nicer to work with

	// to build B and B0, we first iterate through all edges I, and then all edges J which don't intersect I
	for (size_t I = 0; I < Ac.size(); I++)
	{
		for (size_t J_ind = 0; J_ind < Ac[I].size(); J_ind++)
		{
			int J = Ac[I][J_ind];

			// then we just apply a formula from the paper to construct B and B0 based on W and W0
			for (int a = 0; a < 2; a++)
			{
				for (int b = 0; b < 2; b++)
				{
					B(E(I, a), E(I, b)) += std::pow(-1, a + b) * W(I, J) / std::pow(L(I), 2);
					B(E(J, a), E(J, b)) += std::pow(-1, a + b) * W(I, J) / std::pow(L(J), 2);

					B(E(I, a), E(J, b)) -= std::pow(-1, a + b) * W(I, J) * T.row(I).dot(T.row(J)) / (L(I) * L(J));
					B(E(J, a), E(I, b)) -= std::pow(-1, a + b) * W(I, J) * T.row(J).dot(T.row(I)) / (L(J) * L(I));


					B0(E(I, a), E(I, b)) += 0.25 * W0(I, J);
					B0(E(J, a), E(J, b)) += 0.25 * W0(I, J);

					B0(E(I, a), E(J, b)) -= 0.25 * W0(I, J);
					B0(E(J, a), E(I, b)) -= 0.25 * W0(I, J);
				}
			}
		}
	}

	// Now we construct matrix A by adding B and B0
	Eigen::MatrixXd A = B + B0;

	// Construct matrix A_bar
	A_bar.resize(A.cols() * 3, A.rows() * 3);
	A_bar << A, 0 * A, 0 * A,
		0 * A, A, 0 * A,
		0 * A, 0 * A, A;


}




---
File: /src/build_weights.cpp
---

#include "build_weights.h"
#include <Eigen/Dense>

#include <iostream>
#include <vector>

void build_weights(
	const double& alpha,
	const double& beta,
	const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& Ac,
	const Eigen::MatrixX3d& V,
	const Eigen::MatrixX3d& T,
	const Eigen::VectorXd& L,
	Eigen::MatrixXd& W,
	Eigen::MatrixXd& W0)
{
	int edge_num = E.rows();

	// sigma is defined in the paper
	// it is used to compute the value of the repulsion energy
	double sigma = (beta - 1) / alpha;

	// Now, when constructing W (which will happen on each iteration), we iterate through Ac
	W.setZero(edge_num, edge_num);
	// We will construct W0 at the same time (the weights used to make B0)
	W0.setZero(edge_num, edge_num);
	
	// loop through all edges
	for (size_t I = 0; I < Ac.size(); I++)
	{
		// loop through all edges which don't intersect edge I
		for (size_t J_ind = 0; J_ind < Ac[I].size(); J_ind++)
		{
			int J = Ac[I][J_ind]; // the index of the edge which doesn't intersect I

			double elt1 = 0;
			double elt2 = 0;

			// iterate through all combinations of endpoints of these 2 edges
			// use the coordinates of each combination of endpoints for formula from the paper for W and W0

			for (size_t a = 0; a < 2; a++)
			{
				for (size_t b = 0; b < 2; b++)
				{
					int i = E(I, a);
					int j = E(J, b);
					Eigen::Vector3d p = V.row(i);
					Eigen::Vector3d q = V.row(j);
					double diff_norm = (p - q).norm();

					// sometimes when diff_norm gets super small (by 2 vertices overlapping by chance), we divide by 0 and things mess up
					// this is a quick fix
					if (abs(diff_norm) < 0.00000001) {
						elt1 += 1000;
						elt2 += 1000;
					}

					elt1 += 1 / std::pow(diff_norm, 2 * sigma + 1);


					// use alpha = 2 and beta = 4, as specified for B0
					double alph = 2;
					double bet = 4;

					double k = std::pow(((p - q).cross(T.row(I))).norm(), alph) / std::pow(diff_norm, bet);
					elt2 += k / std::pow(diff_norm, 2 * sigma + 1);
				}
			}

			W(I, J) = 0.25 * L(I) * L(J) * elt1; // update the weight matrix W0
			W0(I, J) = 0.25 * L(I) * L(J) * elt2; // update the weight matrix W
		}
	}

}



---
File: /src/constraint_derivative.cpp
---

#include "constraint_derivative.h"
#include <Eigen/Dense>

#include <vector>
#include <iostream>

void constraint_derivative(
	const Eigen::MatrixXi& E,
	const Eigen::MatrixX3d& V,
	const Eigen::VectorXd& L,
	const std::vector<std::vector<int>>& E_adj,
	const Eigen::VectorXi& constr_edges,
	Eigen::MatrixXd& C)
{
	int pt_num = V.rows();

	// First we will find the derivative of the "total length constraint"
	// Phi_l = l0 - Sum over all edges I of {L(I)} 

	Eigen::RowVectorXd C_l;
	C_l.resize(3 * pt_num);

	// loop through all points
	for (size_t p = 0; p < pt_num; p++)
	{
		Eigen::RowVector3d deriv_p(0, 0, 0); // will store the derivative wrt phi_p

		// an edge will contribute to the derivative with respect to phi_p
		// if the edge is adjacent to p

		// loop through all edges adjacent to p
		for (size_t I_ind = 0; I_ind < E_adj[p].size(); I_ind++)
		{
			int I = E_adj[p][I_ind];

			// add or subtract from deriv_p based on the derivative formula

			if (E(I, 0) == p)
			{
				deriv_p -= (V.row(E(I, 0)) - V.row(E(I, 1))) / L(I);
			}
			else if (E(I, 1) == p)
			{
				deriv_p -= (V.row(E(I, 1)) - V.row(E(I, 0))) / L(I);
			}
		}

		// update the derivative matrix
		C_l(3 * p) = deriv_p(0);
		C_l(3 * p + 1) = deriv_p(1);
		C_l(3 * p + 2) = deriv_p(2);
	}


	// We do the same thing as above but now for constraining
	// the length of a specific set of edges "constr_edges"

	Eigen::MatrixXd C_e;
	C_e.setZero(constr_edges.rows(), 3 * pt_num);

	for (size_t i = 0; i < C_e.rows(); i++)
	{
		int I = constr_edges(i);

		int i1 = E(I, 0);
		int i2 = E(I, 1);

		Eigen::Vector3d diff = V.row(i1) - V.row(i2);

		// the formula for the derivative was computed by hand
		// the following code implements it

		C_e(i, i1 * 3 + 0) -= diff(0) / L(I);
		C_e(i, i1 * 3 + 1) -= diff(1) / L(I);
		C_e(i, i1 * 3 + 2) -= diff(2) / L(I);

		C_e(i, i2 * 3 + 0) += diff(0) / L(I);
		C_e(i, i2 * 3 + 1) += diff(1) / L(I);
		C_e(i, i2 * 3 + 2) += diff(2) / L(I);
	}

	// put each constraint in a row of C (the overall constraint derivative matrix which we are returning)
	C.resize(C_l.rows() + C_e.rows(), pt_num * 3);
	C << C_l, C_e;

	
	// the paper recommends that you add in a constraint for the barycenter
	// but I found it was faster to just adjust the barycenter manually after each iteration
	// for the applications I am using repulsive curves for

}



---
File: /src/init_curve.cpp
---

#include "init_curve.h"

#include <iostream>
#include <cmath>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>


void init_curve(
	const Eigen::MatrixXi& E,
	const Eigen::MatrixX3d& V,
	std::vector<std::vector<int>>& Ac,
	std::vector<std::vector<int>>& E_adj,
	Eigen::VectorXd& lengths)
{
	int pt_num = V.rows();
	int edge_num = E.rows();

	// Later, in each iteration of repulsion we will have to construct a weight matrix w_{I,J}
	// To do this we will need to iterate through all pairs of edges which are disjoing

	// In order to save time later, we will just store all edges J which are disjoint from
	// an edge I in a vector Ac[I]
	// So Ac is a vector of vectors of the form Ac[I]

	// Ac stands for "adjacency complement" because it's all the edges not adjacent to other edges
	// This will only depend on the topology of the mesh
	// So we won't need to recompute it

	// We use the "vector" class for this since different rows will have different lengths

	// build Ac
	for (size_t i = 0; i < edge_num; i++)
	{
		std::vector<int> row;
		for (size_t j = 0; j < edge_num; j++)
		{
			// if edge i and j don't share any vertices, add edge j to the row for edge i
			if (!(E(i, 0) == E(j, 0) || E(i, 0) == E(j, 1) || E(i, 1) == E(j, 0) || E(i, 1) == E(j, 1))) {
				row.push_back(j);
			}
		}
		Ac.push_back(row);
	}

	// Note that some of the rows of Ac will be empty, and this is fine


	// Define a vector of vectors called E_Adj
	// One row for each vertex
	// row i contains an element for each edge containing vertex i

	// Similar to "Ac", we only need to build this once since it only depends on the topology
	// of the curve, and it will save us time during each iteration when we need to find the
	// set of edges adjacent to a point.

	// iterate thorugh all points
	for (size_t i = 0; i < pt_num; i++)
	{
		std::vector<int> row;
		// iterate through all edges
		for (size_t j = 0; j < edge_num; j++)
		{
			// if edge i and j don't share any vertices, add edge j to the row for edge i
			if (E(j, 0) == i || E(j, 1) == i) {
				row.push_back(j);
			}

		}

		E_adj.push_back(row);
	}


	// We will get all the edge lengths now
	// This is important because during each iteartion we will want to compare
	// new edge lenghts to the edge lengths of the mesh when it was initialized

	Eigen::VectorXd L;
	L.resize(edge_num);
	for (size_t i = 0; i < edge_num; i++)
	{
		L(i) = (V.row(E(i, 0)) - V.row(E(i, 1))).norm();
	}

	// "lengths" is what gets outputted by the function
	lengths = L;
}



---
File: /src/loss_derivative.cpp
---

#include "loss_derivative.h"

#include <iostream>
#include <cmath>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>


void loss_derivative(
	const double& alpha,
	const double& beta,
	const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& Ac,
	const std::vector<std::vector<int>>& E_adj,
	const Eigen::MatrixX3d& V,
	const Eigen::MatrixX3d& T,
	const Eigen::VectorXd& L,
	Eigen::RowVectorXd& Deriv)
{
	// this file only works if we assume that there are no edges (u,v) such that u=v

	// I derived the derivative by hand so that the code would be faster
	// and we wouldn't have to waste computation power on computing the derivative
	// It is very messy because it had to be broken up into 4 cases
	// But, I checked with mathematica afterwards and it is correct.


	// we want to build the derivative of the repulsion energy with respect to the point gamma_i
	// the energy is a sum of many terms, and the derivative will be 0 if gamma_i isn't in a term
	// Thus, we loop through all disjoint pairs of edges where one of the edges contains gamma_i


	int pt_num = V.rows();


	// first loop through all the points we will be differentiating with respect to
	for (size_t p = 0; p < pt_num; p++)
	{
		Eigen::RowVector3d deriv_p(0, 0, 0);

		// now, loop through all the edges "I" containing vertex "i"
		for (size_t I_ind = 0; I_ind < E_adj[p].size(); I_ind++)
		{
			int I = E_adj[p][I_ind];

			// loop through all edges "J" not intersecting "I"
			for (size_t J_ind = 0; J_ind < Ac[I].size(); J_ind++)
			{
				int J = Ac[I][J_ind];

				// iterate through combinations of verts from the 2 edges
				for (size_t i = 0; i < 2; i++)
				{
					for (size_t j = 0; j < 2; j++)
					{

						// note that by construction, only I will have point p in it
						// the following if statement is just so that I can know which index of i has p

						if (E(I, i) == p) {

							// relabel the indices so "i1" is the vertex index equal to p
							// "i2" is the one that isn't p
							// this is just so it matches the math I wrote on paper
							int i1 = E(I, i);
							int i2 = E(I, (i + 1) % 2);

							// to fit with notation, we will also set
							int j1 = E(J, j);

							// what I'm about to write won't depend on j, so we don't need to break it up into cases for j

							// Case 1,1

							// the messy cross product in the numerator
							Eigen::Vector3d cross_term = (V.row(i2) - V.row(j1)).cross(V.row(i1)) - V.row(i2).cross(V.row(j1));

							// the difference on the denominator
							Eigen::Vector3d denom_diff = V.row(i1) - V.row(j1);

							Eigen::Vector3d e1(1, 0, 0);
							Eigen::Vector3d e2(0, 1, 0);
							Eigen::Vector3d e3(0, 0, 1);

							Eigen::Vector3d matrix_cross = (Eigen::Vector3d)(V.row(i2)-V.row(j1));

							Eigen::Matrix3d TMat = (matrix_cross.cross(e1)) * e1.transpose() + (matrix_cross.cross(e2)) * e2.transpose() + (matrix_cross.cross(e3)) * e3.transpose();

							Eigen::RowVector3d term1 = (1 - alpha)* std::pow(L(I), -1 * alpha - 1)* (V.row(i1) - V.row(i2))* std::pow(cross_term.norm(), alpha)* std::pow(denom_diff.norm(), -1 * beta);

							Eigen::RowVector3d term2 = alpha * std::pow(L(I), 1 - alpha) * std::pow(denom_diff.norm(), -1 * beta) * std::pow(cross_term.norm(), alpha - 2) * cross_term.transpose() * TMat;

							Eigen::RowVector3d term3 = -1 * beta * std::pow(L(I), 1 - alpha) * std::pow(denom_diff.norm(), -1 * beta - 2) * std::pow(cross_term.norm(), alpha) * denom_diff.transpose();

							deriv_p += 0.25 * L(J) * (term1 + term2 + term3);

							
							// Case 2,1

							denom_diff = V.row(i2) - V.row(j1);

							term1 = (1 - alpha) * std::pow(L(I), -1 * alpha - 1) * (V.row(i1) - V.row(i2)) * std::pow(cross_term.norm(), alpha) * std::pow(denom_diff.norm(), -1 * beta);

							term2 = alpha * std::pow(L(I), 1 - alpha) * std::pow(denom_diff.norm(), -1 * beta) * std::pow(cross_term.norm(), alpha - 2) * cross_term.transpose() * TMat;

							deriv_p += 0.25 * L(J) * (term1 + term2);


							// Now cases with J
							// Case J : 1, 1

							Eigen::Vector3d TJ = T.row(J);
							TMat = (TJ.cross(e1)) * e1.transpose() + (TJ.cross(e2)) * e2.transpose() + (TJ.cross(e3)) * e3.transpose();

							denom_diff = V.row(i1) - V.row(j1);

							cross_term = TJ.cross(denom_diff);

							term1 = std::pow(cross_term.norm(), alpha) * std::pow(denom_diff.norm(), -1 * beta) * denom_diff.transpose() / L(I);

							term2 = cross_term.transpose() * TMat;
							term2 = term2 * alpha * L(I) * std::pow(denom_diff.norm(), -1 * beta) * std::pow(cross_term.norm(), alpha - 2);

							term3 = -1 * beta * L(I) * std::pow(denom_diff.norm(), -1 * beta - 2) * std::pow(cross_term.norm(), alpha) * denom_diff.transpose();

							deriv_p += 0.25 * L(J) * (term1 + term2 + term3);


							// Case J : 1, 2

							denom_diff = V.row(i2) - V.row(j1);

							cross_term = TJ.cross(denom_diff);

							deriv_p += 0.25 * (L(J) / L(I)) * std::pow(cross_term.norm(), alpha) * std::pow(denom_diff.norm(), -1 * beta) * (V.row(i1) - V.row(i2));

						}
						
					}
				}

			}

		}

		// put values in the derivative vector for the derivative with respect to p (which is a vector) so it corresponds to 3 elements
		Deriv(3 * p) = deriv_p(0);
		Deriv(3 * p + 1) = deriv_p(1);
		Deriv(3 * p + 2) = deriv_p(2);
	}

	
}



---
File: /src/read_curve_file.cpp
---

#include "read_curve_file.h"  
#include <cmath>  
#include <Eigen/Dense>  
#include <iostream>  

#include <vector>
#include <fstream>
#include <string>
#include <sstream>

#include <iterator>

void read_curve_file(
    const std::string& input,
    Eigen::MatrixX3d& V,
    Eigen::MatrixXi& E)
{

    // First, we find the number of vertices and number of edges by looping through lines in the obj file
    std::ifstream file(input);

    int v_ind = 0;
    int e_ind = 0;

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            // if the line starts with 'v' (corresponding to a vert)
            if (line[0] == 'v') {
                v_ind++;
            }

            // if the line starts with 'l' (corresponding to an edge)
            if (line[0] == 'l') {
                e_ind++;
            }


        }
    }

    file.close();


    // set the sizes of V and E
    V.resize(v_ind, 3);
    E.resize(e_ind, 2);


    // now, we will go through the file again and actually set the elements of V and E to the values in the file
    std::ifstream file2(input);

    v_ind = 0;
    e_ind = 0;

    if (file2.is_open()) {
        std::string line;
        while (std::getline(file2, line)) {
            // if the line starts with 'v' (corresponding to a vert)
            if (line[0] == 'v') {
                line.erase(0, 1);
                std::istringstream iss(line);

                std::vector<std::string> tokens;

                // split up string by "space" delimeter
                std::copy(std::istream_iterator<std::string>(iss),
                    std::istream_iterator<std::string>(),
                    std::back_inserter(tokens));


                V(v_ind, 0) = std::stod(tokens[0]);
                V(v_ind, 1) = std::stod(tokens[1]);
                V(v_ind, 2) = std::stod(tokens[2]);

                v_ind++;
            }

            // if the line starts with 'l' (corresponding to an edge)
            if (line[0] == 'l') {
                line.erase(0, 1);
                std::istringstream iss(line);

                std::vector<std::string> tokens;

                // split up string by "space" delimeter
                std::copy(std::istream_iterator<std::string>(iss),
                    std::istream_iterator<std::string>(),
                    std::back_inserter(tokens));

                E(e_ind, 0) = std::stoi(tokens[0]) - 1;
                E(e_ind, 1) = std::stoi(tokens[1]) - 1;

                e_ind++;
            }


        }
    }

    file.close();


}


---
File: /src/repulsion_iteration.cpp
---

#include "repulsion_iteration.h"  

#include "loss_derivative.h"  
#include "build_A_bar.h"  
#include "read_curve_file.h"  

#include "build_weights.h"

#include "constraint_derivative.h"

#include <iostream>
#include <cmath>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>


void repulsion_iteration(
	const double& alpha,
	const double& beta,
	const double& a_const,
	const double& b_const,
	const double& threshold,
	const int& max_iters,
	const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& Ac,
	const std::vector<std::vector<int>>& E_adj,
	const Eigen::VectorXd& lengths,
	Eigen::MatrixX3d& V)
{
	
	std::cout << "Iteration in Progress" << std::endl;


	// number of points and number of edges
	int pt_num = V.rows();
	int edge_num = E.rows();

	// We want to constrain the length of edges connected to vertices
	// with 3 or more edges connected to the vertex.
	// Otherwise the tangent point energy will casue the graph will explode

	// number of verts with 3 or more edges
	int three_way = 0;
	// Count the number of constrain edges for graph embeddings
	for (size_t i = 0; i < pt_num; i++)
	{
		if (E_adj[i].size() > 2) {
			three_way += E_adj[i].size();
		}
	}

	// vector containing the indices of the edges whose lengths we want to fix
	Eigen::VectorXi constr_edges;
	constr_edges.resize(three_way);

	int n = 0;
	for (size_t i = 0; i < pt_num; i++)
	{
		// if there are more than 2 edges adjacent to this edge
		if (E_adj[i].size() > 2) {
			for (int I_ind = 0; I_ind < E_adj[i].size(); I_ind++) {

				int I = E_adj[i][I_ind];
				constr_edges(n) = I;
				n++;

			}			
		}
	}

	// total length of the original mesh before any of the iterations
	// it will be used later for the "total length constraint"
	double l0 = lengths.sum();

	// Get all the edge lengths
	Eigen::VectorXd L;
	L.resize(edge_num);
	for (size_t i = 0; i < edge_num; i++)
	{
		L(i) = (V.row(E(i, 0)) - V.row(E(i, 1))).norm();
	}

	Eigen::MatrixX3d T;
	T.resize(edge_num, 3);
	// for each edge I, we have a corresponding T_I vector (representing the tangent vector)
	for (size_t i = 0; i < edge_num; i++)
	{
		T.row(i) = (V.row(E(i, 1)) - V.row(E(i, 0))) / L(i);
	}

	// Build a weight matrix as described in the paper. It will be used to build
	// a matrix A_bar which is used for computing Sobolev inner products
	Eigen::MatrixXd W;
	Eigen::MatrixXd W0;
	build_weights(alpha, beta, E, Ac, V, T, L, W, W0);

	// Construct matrix A_bar
	Eigen::MatrixXd A_bar;
	A_bar.resize(pt_num * 3, pt_num * 3);
	build_A_bar(pt_num, E, Ac, T, L, W, W0, A_bar);

	
	// build derivative of the loss function
	Eigen::RowVectorXd Deriv;
	Deriv.setZero(3 * pt_num);
	loss_derivative(alpha, beta, E, Ac, E_adj, V, T, L, Deriv);

	// Make constraint derivative matrix
	Eigen::MatrixXd C;
	constraint_derivative(E, V, L, E_adj, constr_edges, C);

	int k = C.rows(); // number of constraints 

	Eigen::MatrixXd Z;
	Z.setZero(k, k);


	// SOLVE FOR UNKNOWN in equation
	// Left*Unknown = Right

	// the first 3*pt_num elements of "Unknown" give the descent direction

	// make the matrix on the left in the equation
	Eigen::MatrixXd Left;
	Left.resize(k + pt_num * 3, k + pt_num * 3);
	Left << A_bar, C.transpose(), C, Z;

	// make the vector on the right in the equation we are solving
	Eigen::VectorXd Right;
	Right.setZero(k + pt_num * 3);
	for (size_t i = 0; i < pt_num * 3; i++)
	{
		Right(i) = Deriv(i);
	}

	Eigen::VectorXd Unknown;
	
	Unknown = Left.fullPivLu().solve(Right);
	// we use full PivLu because partial pivoting sometimes results in "unknown" having 
	// undefined values

	// "g_tilde" from the paper is the first pt_num * 3 elements of Unknown
	
	
	// The step size is t. We start it at 1 and use backtracking line search to reduce it
	// Note that starting at 1 works because we will first normalize the descent direction and the derivative
	double t = 1;

	// Matrix to store the updated vertex position so we can compare it to the original
	Eigen::MatrixX3d V_new;
	V_new.resize(V.rows(), 3);

	// the descent direction as computed from the equation where we solved for "Unknown"
	Eigen::VectorXd descent_dir;
	descent_dir.resize(pt_num * 3);
	for (size_t i = 0; i < pt_num; i++)
	{
		descent_dir(3 * i + 0) = -1 * Unknown(3 * i + 0);
		descent_dir(3 * i + 1) = -1 * Unknown(3 * i + 1);
		descent_dir(3 * i + 2) = -1 * Unknown(3 * i + 2);
	}
	
	// normalize the gradient and descent direction before doing backtracking line search
	descent_dir = descent_dir/descent_dir.norm();
	Deriv = Deriv / Deriv.norm(); 

	// backtracking line search
	while (1) {

		// first take a step with step size t and store new vertex positions in V_new
		for (size_t i = 0; i < pt_num; i++)
		{
			V_new(i, 0) = V(i, 0) - t * Unknown(3 * i + 0);
			V_new(i, 1) = V(i, 1) - t * Unknown(3 * i + 1);
			V_new(i, 2) = V(i, 2) - t * Unknown(3 * i + 2);
		}

		// this vector will store the value of each constraint function at the new vertices	
		Eigen::VectorXd constraint_vals;
		constraint_vals.resize(k);

		// compute new edge lengths
		Eigen::VectorXd L_new;
		L_new.resize(L.rows());
		for (size_t i = 0; i < edge_num; i++)
		{
			L_new(i) = (V_new.row(E(i, 0)) - V_new.row(E(i, 1))).norm();
		}

		// compute new tangent vectors
		Eigen::MatrixX3d T_new;
		T_new.resize(edge_num, 3);
		for (size_t i = 0; i < edge_num; i++)
		{
			T_new.row(i) = (V_new.row(E(i, 1)) - V_new.row(E(i, 0))) / L_new(i);
		}


		// We project the step we have taken onto the constraint space in the following loop
		// This is done by solving an iterative equation
		// We terminate either when the constriants are below "threshold" or when we have done a max number of iterations
		// The max number of iterations is here so that we don't waste time making the constraint error small when 
		// just decreasing the step size will work better
		int iter = 0;
		while (iter < max_iters) {
			// the first constraint is that the total length doesn't change
			// note that l0 is the length of the mesh when it was initialized
			constraint_vals(0) = l0 - L_new.sum();

			// the remaining constraints are that the
			for (size_t i = 1; i < k; i++)
			{
				constraint_vals(i) = lengths(constr_edges(i - 1)) - L_new(constr_edges(i - 1));
			}

			if (constraint_vals.norm() <= threshold) {
				break;
			}

			// to do the projection, we solve a matrix equation
			// Left * Unknonw2 = Right2

			// "Left" is the same as the one from the previous equation

			// Now, construct the "Right2" vector
			Eigen::VectorXd Right2;
			Right2.setZero(k + 3 * pt_num);
			for (size_t i = 0; i < k; i++)
			{
				Right2(i + 3 * pt_num) = -1 * constraint_vals(i);
			}

			Left.resize(k + pt_num * 3, k + pt_num * 3);
			Left << A_bar, C.transpose(), C, Z;

			// Solve the matrix equation and update the vertex positions accordingly
			Eigen::VectorXd Unknown2;
			//Unknown2 = Left.lu().solve(Right2);
			Unknown2 = Left.fullPivLu().solve(Right2);
			for (size_t i = 0; i < pt_num; i++)
			{
				V_new(i, 0) += Unknown2(3 * i + 0);
				V_new(i, 1) += Unknown2(3 * i + 1);
				V_new(i, 2) += Unknown2(3 * i + 2);
			}

			iter++;

			// recompute lengths and tangents based on the new vertex positions before the next iter
			for (size_t i = 0; i < edge_num; i++)
			{
				L_new(i) = (V_new.row(E(i, 0)) - V_new.row(E(i, 1))).norm();
			}
			for (size_t i = 0; i < edge_num; i++)
			{
				T_new.row(i) = (V_new.row(E(i, 1)) - V_new.row(E(i, 0))) / L_new(i);
			}
		}


		// for backtracking line search, we want f_delta < "right side"
		// f_delta = f(x + tdelta(x)) where delta(x) is the direction of descent
		// right side = f(x) + a_const * t * the dot product between the gradient and the direction of descent
		// where f is the energy we want to minimize

		double f_delta = 0;
		double right_side = descent_dir.dot((Eigen::VectorXd)Deriv);

		right_side *= t * a_const;

		// we still need to add f(x) to "right_side"
		// and add f(x + t*delta(x)) to f_delta

		// we do this in the following nested for loop
		for (size_t I = 0; I < edge_num; I++)
		{
			for (size_t J_ind = 0; J_ind < Ac[I].size(); J_ind++)
			{
				int J = Ac[I][J_ind];

				for (size_t i = 0; i < 2; i++)
				{
					for (size_t j = 0; j < 2; j++)
					{
						Eigen::Vector3d d = V.row(E(I, i)) - V.row(E(J, j));
						double d_norm = d.norm();

						double cross_norm = (T.row(I)).cross(d).norm();

						right_side += L(I) * L(J) * 0.25 * std::pow(cross_norm, alpha) * std::pow(d_norm, -1*beta);



						Eigen::Vector3d d_new = V_new.row(E(I, i)) - V_new.row(E(J, j));
						double d_norm_new = d_new.norm();

						double cross_norm_new = (T_new.row(I)).cross(d_new).norm();

						f_delta += L_new(I) * L_new(J) * 0.25 * std::pow(cross_norm_new, alpha) * std::pow(d_norm_new, -1 * beta);
					}
				}
			}
		}

		if (f_delta <= right_side && constraint_vals.norm() < threshold) {
			V = V_new;
			break;
		}

		// decrease the step size
		t = t * b_const;
	}


	// in updating our vertex positions, we might have moved the entire mesh
	// we fix this by translating the mesh so the barycenter is 0

	// compute barycenter
	Eigen::Vector3d x0_new(0, 0, 0); 
	for (size_t i = 0; i < edge_num; i++)
	{
		x0_new += L(i) * 0.5 * (V.row(E(i, 0)) + V.row(E(i, 1)));
	}
	x0_new = x0_new / L.sum();

	// subtract off barycenter
	for (size_t i = 0; i < pt_num; i++)
	{
		V.row(i) = V.row(i) - (Eigen::RowVector3d)x0_new;
	}
	
	std::cout << "Iteration Complete" << std::endl;

}


---
File: /main.cpp
---

#include "init_curve.h"
#include "repulsion_iteration.h"
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <iostream>

#include <vector>

#include "read_curve_file.h"

// This project was made by Nathan Henry
// The main file was adapted from the registration assignment in CSC419
// PLEASE READ THE README.pdf

int main(int argc, char* argv[])
{

    // ----------------------------------------------------------------------------------------------------------------
    // THESE ARE THE PARAMETERS THAT CAN BE TUNED -- more details in the README pdf

    // copy the file path of the curve you want to use the program on
    // since I wrote a custom function to parse the file, IT ONLY ACCEPTS OBJ FILES!!!!!!!!
    std::string curve_file_path = "C:/Users/Nathan Henry/Desktop/Geometry_processing/nathan-henry-repulsive-curves/data/petersen.obj";

    // -------------

    // alpha and beta are parameters that change the nature of the tangent-point energy
    // (See Figure 4 in the paper for intuition on what they do)
    // values of 3 and 6 are recommended
    double alpha = 3;
    double beta = 6;

    // -------------

    double a_const = 0.01;
    double b_const = 0.9;
    // These are 2 parameters used in my implementation of line search
    // a_const must be between 0 and 0.5
    // b_const must be between 0 and 1
   
    // decreasing a_const gives preference to larger step sizes
    // though, these steps will have a tendency to overshoot

    // increasing b_const can increase the step sizes in a way where it won't overshoot as much if it does overshoot
    // however, increasing b_const results in each step taking a longer time to compute

    // empirically, I have found that for the examples I tried a_const = 0.01 and b_const = 0.9 work well

    
    // -------------

    int max_iters = 20;
    // THIS IS NOT THE MAXIMUM NUMBER OF DESCENT STEPS
    // for each descent step, we will have to project the result of the step back onto constrained configuration space
    // this projection is iterative, and the maximum number is something we can set here
    // usually it should only take 3 iterations, but this is set to 20 to deal with edge cases
    // on much bigger meshes, it should be set to be higher (ex. 50)
    // on very small meshes (less than 10 vertices), I have found that it's faster if you set max_iters to 10
    // though, this is not an absolute rule

    // -------------

    double threshold = 0.001;
    // This is the maximum absolute value of the constraint function
    // In this case, it means that the sum of squares of the total length deviation with the length
    // deviation of each constrained edge must be less than 0.001
    
    // ----------------------------------------------------------------------------------------------------------------


    Eigen::MatrixXd OVX, VX;
    Eigen::MatrixXi FX, E;
    Eigen::MatrixX3d OV, V;

    std::vector<std::vector<int>> Ac;
    std::vector<std::vector<int>> E_adj;
    Eigen::VectorXd lengths;


    // Load mesh which shows x,y,z axis
    igl::read_triangle_mesh(
        ("C:/Users/Nathan Henry/Desktop/Geometry_processing/nathan-henry-repulsive-curves/data/x_y_z_axis.obj"), OVX, FX);

    // I made a cpp file to read OBJ files which contain curves (i.e. no faces)
    read_curve_file(argc > 1 ? argv[1] : curve_file_path, OV, E);

    V.resize(OV.rows(), 3);

    double factor = OV.maxCoeff() / OVX.maxCoeff();
    OVX = OVX * factor / 1.25; // scale axis so it's a bit smaller than the curve



    bool show_samples = true;
    igl::opengl::glfw::Viewer viewer;
    const int xid = viewer.selected_data_index;
    viewer.append_mesh();

    std::cout << R"(
  [space]  toggle animation
  h        to a single step with fractional Sobolev descent
  r        reset the curve to its initial state
  p        hide points and edges of curve
)";

    // predefined color
    const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
    const auto& set_points = [&]()
    {
        if (show_samples)
        {
            // show verts and edges
            viewer.data_list[xid].set_points(V, (1. - (1. - orange.array()) * .8));
            viewer.data_list[xid].set_edges(V, E, Eigen::RowVector3d(1, 1, 1));
        }
        else
        {
            viewer.data_list[xid].clear_points();
            viewer.data_list[xid].clear_edges();
        }
    };
    const auto& reset = [&]()
    {
        V = OV;
        Ac.clear();
        E_adj.clear();

        // based on E and V, compute Ac, E_adj, and lengths (of edges in the curve's initial state)
        // Ac and E_adj as well as their motivation are described in the init_curve.cpp file
        init_curve(E, V, Ac, E_adj, lengths);

        set_points();
    };
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool
    {
        if (viewer.core().is_animating)
        {
            // do a sobolev descent step and update V
            repulsion_iteration(alpha, beta, a_const, b_const, threshold, max_iters, E, Ac, E_adj, lengths, V);

            set_points();
        }
        return false;
    };
    viewer.callback_key_pressed =
        [&](igl::opengl::glfw::Viewer&, unsigned char key, int)->bool
    {
        switch (key)
        {
        case ' ':
            viewer.core().is_animating ^= 1;
            break;
        case 'H':
        case 'h':
            // do a sobolev descent step and update V
            repulsion_iteration(alpha, beta, a_const, b_const, threshold, max_iters, E, Ac, E_adj, lengths, V);
            set_points();
            break;
        case 'M':
        case 'm':
        case 'P':
        case 'p':
            show_samples ^= 1;
            set_points();
            break;
        case 'R':
        case 'r':
            reset();
            break;
        case 'S':
        case 's':
        default:
            return false;
        }
        return true;
    };

    viewer.data_list[xid].set_mesh(OVX, FX);
    viewer.data_list[xid].set_colors(orange);

    reset();
    viewer.core().is_animating = false;
    viewer.data().point_size = 10;
    viewer.launch();

    return EXIT_SUCCESS;
}

Directory Structure:

└── ./
    └── src
        ├── lib
        │   ├── energyCalculations.js
        │   ├── graphDrawing.js
        │   ├── graphstate.js
        │   ├── graphUtils.js
        │   ├── innerProduct.js
        │   ├── interaction.js
        │   ├── optimization.js
        │   └── stores.js
        └── routes
            └── +page.svelte



---
File: /src/lib/energyCalculations.js
---

// src/lib/energyCalculations.js
import * as math from 'mathjs';
import { get } from 'svelte/store';
import { config } from '$lib/stores';

let logging = false;

export function calculateEdgeProperties(vertices, edges) {
	const edgeLengths = [];
	const edgeTangents = [];
	const edgeMidpoints = [];

	for (const edge of edges) {
		const v1 = vertices[edge[0]];
		const v2 = vertices[edge[1]];

		const dx = v2[0] - v1[0];
		const dy = v2[1] - v1[1];
		const length = Math.sqrt(dx * dx + dy * dy);
		edgeLengths.push(length);

		const unitTangent = length > 0 ? [dx / length, dy / length] : [0, 0];
		edgeTangents.push(unitTangent);

		const midpoint = [
			isNaN(v1[0]) || isNaN(v2[0]) ? 0 : (v1[0] + v2[0]) / 2,
			isNaN(v1[1]) || isNaN(v2[1]) ? 0 : (v1[1] + v2[1]) / 2
		];
		edgeMidpoints.push(midpoint);

		if (logging) {
			console.log(
				`Edge [${edge[0]}, ${edge[1]}]: length = ${length}, tangent = ${unitTangent}, midpoint = ${midpoint}`
			);
		}
	}

	if (logging) {
		console.log('Edge lengths:', edgeLengths);
		console.log('Unit tangents:', edgeTangents);
		console.log('Midpoints:', edgeMidpoints);
	}

	return { edgeLengths, edgeTangents, edgeMidpoints };
}

export function tangentPointKernel(p, q, T, alpha, beta) {
	// Ensure inputs are properly converted to matrices
	const p_ = math.matrix(p);
	const q_ = math.matrix(q);
	const T_ = math.matrix(T);
	const epsilon = get(config).epsilonKernel;

	// Calculate the difference vector
	const diff = math.subtract(p_, q_);
	const diffNorm = math.norm(diff) + epsilon; // Prevent division by zero
	const cross2D = T_.get([0]) * diff.get([1]) - T_.get([1]) * diff.get([0]); // 2D cross product (determinant)

	const numerator = Math.pow(Math.abs(cross2D), alpha);
	const denominator = Math.pow(diffNorm, beta);
	const result = numerator / denominator;

	// Check for NaN or Infinity
	if (!isFinite(result)) {
		console.warn('Invalid kernel result:', result, 'from inputs:', p, q, T, 'with cross2D:', cross2D, 'diffNorm:', diffNorm);
		return 0;
	}

	return result;
}

export function calculateDisjointEdgePairs(edges) {
	const numEdges = edges.length;
	const disjointPairs = [];

	for (let i = 0; i < numEdges; i++) {
		disjointPairs[i] = [];
		for (let j = 0; j < numEdges; j++) {
			if (i === j) continue;

			const edge1 = edges[i];
			const edge2 = edges[j];

			if (
				edge1[0] !== edge2[0] &&
				edge1[0] !== edge2[1] &&
				edge1[1] !== edge2[0] &&
				edge1[1] !== edge2[1]
			) {
				disjointPairs[i].push(j);
			}
		}
	}
	console.log('Calculated disjointPairs:', disjointPairs);
	return disjointPairs;
}

export function calculateDiscreteKernel(vertices, edges, edgeTangents, alpha, beta, disjointPairs) {
	const numEdges = edges.length;
	const kernelMatrix = math.zeros(numEdges, numEdges);

	if (!disjointPairs || !Array.isArray(disjointPairs) || disjointPairs.length === 0) {
		console.warn('No disjoint pairs found, returning zero kernel matrix');
		return kernelMatrix;
	}

	for (let i = 0; i < numEdges; i++) {
		if (!disjointPairs[i]) {
			console.warn(`No disjoint pairs for edge ${i}`);
			continue;
		}

		for (const j of disjointPairs[i]) {
			if (i < edges.length && j < edges.length) {
				let sum = 0;
				const combinations = [
					[vertices[edges[i][0]], vertices[edges[j][0]]],
					[vertices[edges[i][0]], vertices[edges[j][1]]],
					[vertices[edges[i][1]], vertices[edges[j][0]]],
					[vertices[edges[i][1]], vertices[edges[j][1]]]
				];

				for (const [p, q] of combinations) {
					sum += tangentPointKernel(p, q, edgeTangents[i], alpha, beta);
				}
				kernelMatrix.set([i, j], sum / 4);
				kernelMatrix.set([j, i], sum / 4); // Keep symmetry!
			} else {
				console.warn(
					'Invalid edge index:',
					i,
					j,
					'disjointPairs:',
					disjointPairs,
					'edges.length',
					edges.length
				);
			}
		}
	}
	return kernelMatrix;
}

export function calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs) {
	const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);
	const kernelMatrix = calculateDiscreteKernel(
		vertices,
		edges,
		edgeTangents,
		alpha,
		beta,
		disjointPairs
	);

	let totalEnergy = 0;
	const numEdges = edges.length;

	for (let i = 0; i < numEdges; i++) {
		for (const j of disjointPairs[i]) {
			if (i < edges.length && j < edges.length) {
				const kernelValue = kernelMatrix.get([i, j]);
					totalEnergy += kernelValue * edgeLengths[i] * edgeLengths[j];
			}
		}
	}
	return totalEnergy / 2; // Divide by 2 because of symmetry
}

function calculateAnalyticalDifferential(vertices, edges, alpha, beta, disjointPairs) {
    const numVertices = vertices.length;
    const differential = Array(numVertices).fill().map(() => [0, 0]);
    const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);

    for (let p = 0; p < numVertices; p++) {
        let deriv_p = [0, 0];
        const adjacentEdges = edges.filter(([v0, v1]) => v0 === p || v1 === p);

        for (const edgeI of adjacentEdges) {
            const I = edges.indexOf(edgeI);
            const i = edgeI[0] === p ? 0 : 1;
            const i1 = edgeI[i];
            const i2 = edgeI[(i + 1) % 2];
            const l_I = edgeLengths[I];
            const T_I = edgeTangents[I];

            for (const J of disjointPairs[I]) {
                const l_J = edgeLengths[J];
                const T_J = edgeTangents[J];

                for (let j = 0; j < 2; j++) {
                    const j1 = edges[J][j];
                    const p_i1 = vertices[i1];
                    const p_i2 = vertices[i2];
                    const p_j1 = vertices[j1];

                    // Terms from loss_derivative.cpp adapted for 2D
                    const cross_term = [
                        (p_i2[0] - p_j1[0]) * T_I[1] - (p_i2[1] - p_j1[1]) * T_I[0],
                        (p_i1[0] - p_j1[0]) * T_I[1] - (p_i1[1] - p_j1[1]) * T_I[0]
                    ];
                    const cross_norm = Math.sqrt(cross_term[0] * cross_term[0] + cross_term[1] * cross_term[1]);
                    const denom_diff_i1_j1 = [p_i1[0] - p_j1[0], p_i1[1] - p_j1[1]];
                    const denom_diff_i2_j1 = [p_i2[0] - p_j1[0], p_i2[1] - p_j1[1]];
                    const denom_norm_i1_j1 = Math.sqrt(denom_diff_i1_j1[0] * denom_diff_i1_j1[0] + denom_diff_i1_j1[1] * denom_diff_i1_j1[1]);
                    const denom_norm_i2_j1 = Math.sqrt(denom_diff_i2_j1[0] * denom_diff_i2_j1[0] + denom_diff_i2_j1[1] * denom_diff_i2_j1[1]);

                    // Analytical derivative terms (simplified for 2D)
                    const term1 = (1 - alpha) * Math.pow(l_I, -alpha - 1) * [p_i1[0] - p_i2[0], p_i1[1] - p_i2[1]] * Math.pow(cross_norm, alpha) * Math.pow(denom_norm_i1_j1, -beta);
                    // Add more terms as needed from loss_derivative.cpp
                    deriv_p[0] += 0.25 * l_J * term1[0];
                    deriv_p[1] += 0.25 * l_J * term1[1];
                }
            }
        }
        differential[p] = deriv_p;
    }
    return differential;
}

export function calculateDifferential(vertices, edges, alpha, beta, disjointPairs) {
    const method = get(config).differentialMethod;
    if (method === 'finiteDifference') {
        return calculateDifferentialFiniteDifference(vertices, edges, alpha, beta, disjointPairs);
    } else if (method === 'analytical') {
        return calculateAnalyticalDifferential(vertices, edges, alpha, beta, disjointPairs);
    } else {
        throw new Error('Unknown method for differential calculation');
    }
}

function calculateDifferentialFiniteDifference(vertices, edges, alpha, beta, disjointPairs) {
    // Existing finite difference implementation
    const h = get(config).finiteDiffH;
    const numVertices = vertices.length;
    const differential = [];

    const E_original = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

    for (let i = 0; i < numVertices; i++) {
        differential[i] = [0, 0];
        for (let dim = 0; dim < 2; dim++) {
            const vertices_perturbed = vertices.map((v) => [...v]);
            vertices_perturbed[i][dim] += h;
            const E_perturbed = calculateDiscreteEnergy(
                vertices_perturbed,
                edges,
                alpha,
                beta,
                disjointPairs
            );
            differential[i][dim] = (E_perturbed - E_original) / h;
        }
    }
    console.log('Computed differential:', differential);
    return differential;
}

function calculateL2Gradient(vertices, edges, alpha, beta, disjointPairs) {
	const h = get(config).finiteDiffH; // Use config finiteDiffH
	const numVertices = vertices.length;
	const gradient = [];

	const originalEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

	for (let i = 0; i < numVertices; i++) {
		gradient[i] = [0, 0];

		for (let j = 0; j < 2; j++) {
			const originalValue = vertices[i][j];
			vertices[i][j] += h;

			const newEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

			gradient[i][j] = (newEnergy - originalEnergy) / h;

			vertices[i][j] = originalValue;
		}
	}

	return gradient;
}



---
File: /src/lib/graphDrawing.js
---

// src/lib/graphDrawing.js
import * as math from 'mathjs';
import { drawArrow } from '$lib/graphUtils';
import { calculateDifferential } from '$lib/energyCalculations';
import { canvasTransform, subvertices } from '$lib/stores';
import { get } from 'svelte/store';

export function drawGraph(
	ctx,
	width,
	height,
	vertices,
	edges,
	edgeProps,
	kernelMatrix,
	alpha,
	beta,
	disjointPairs
) {
	const { offsetX, offsetY, zoom } = get(canvasTransform);
	const subs = get(subvertices);

	ctx.save();
	ctx.clearRect(0, 0, width, height);

	ctx.scale(zoom, zoom);
	ctx.translate(offsetX / zoom, offsetY / zoom);

	drawEdges(ctx, vertices, edges, kernelMatrix);
	drawVertices(ctx, vertices, edges, alpha, beta, disjointPairs);
	drawMidpoints(ctx, edges, edgeProps);
	drawSubvertices(ctx, subs);

	ctx.restore();
}

function drawEdges(ctx, vertices, edges, kernelMatrix) {
	const { zoom } = get(canvasTransform);

	if (kernelMatrix && math.isMatrix(kernelMatrix) && kernelMatrix.size()[0] === edges.length) {
		const maxKernelValue = math.max(kernelMatrix) || 1;

		edges.forEach((edge, i) => {
			const totalKernel = edges.reduce((sum, _, j) => {
				return i !== j ? sum + kernelMatrix.get([i, j]) : sum;
			}, 0);
			const avgKernel = totalKernel / (edges.length - 1 || 1);
			const normalizedValue = avgKernel / maxKernelValue;

			const blue = Math.round(255 * (1 - normalizedValue));
			const red = Math.round(255 * normalizedValue);

			ctx.strokeStyle = `rgb(${red}, 0, ${blue})`;
			ctx.lineWidth = (1 + normalizedValue * 4) / zoom;

			ctx.beginPath();
			ctx.moveTo(vertices[edge[0]][0], vertices[edge[0]][1]);
			ctx.lineTo(vertices[edge[1]][0], vertices[edge[1]][1]);
			ctx.stroke();
		});
	} else {
		ctx.strokeStyle = 'black';
		ctx.lineWidth = 1 / zoom;
		edges.forEach((edge) => {
			ctx.beginPath();
			ctx.moveTo(vertices[edge[0]][0], vertices[edge[0]][1]);
			ctx.lineTo(vertices[edge[1]][0], vertices[edge[1]][1]);
			ctx.stroke();
		});
	}
}

function drawVertices(ctx, vertices, edges, alpha, beta, disjointPairs) {
    const { offsetX, offsetY, zoom } = get(canvasTransform);
    const gradient = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);

    vertices.forEach((vertex, i) => {
        ctx.beginPath();
        ctx.arc(vertex[0], vertex[1], 5 / zoom, 0, 2 * Math.PI);
        ctx.fillStyle = 'blue';
        ctx.fill();

        const screenX = vertex[0] * zoom + offsetX;
        const screenY = vertex[1] * zoom + offsetY;

        ctx.save();
        ctx.setTransform(1, 0, 0, 1, 0, 0);

        ctx.fillStyle = 'black';
        ctx.font = '12px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'bottom';
        ctx.fillText(i.toString(), screenX, screenY - 10);

        const gradX = -gradient[i][0];
        const gradY = -gradient[i][1];
        const magnitude = Math.sqrt(gradX * gradX + gradY * gradY) || 1e-6;
        const screenGradX = gradX * zoom;
        const screenGradY = gradY * zoom;
        const screenMagnitude = Math.sqrt(screenGradX * screenGradX + screenGradY * screenGradY);
        const arrowLength = 20;
        const dirX = screenGradX / screenMagnitude;
        const dirY = screenGradY / screenMagnitude;
        ctx.strokeStyle = 'purple';
        ctx.lineWidth = 2;
        drawArrow(ctx, screenX, screenY, dirX, dirY, arrowLength);

        ctx.restore();
    });
}

function drawMidpoints(ctx, edges, edgeProps) {
    const { offsetX, offsetY, zoom } = get(canvasTransform);

    if (!edgeProps || !edgeProps.edgeMidpoints || edgeProps.edgeMidpoints.length !== edges.length) {
        return;
    }

    edges.forEach((edge, i) => {
        const midpoint = edgeProps.edgeMidpoints[i];
        const length = edgeProps.edgeLengths[i];
        const tangent = edgeProps.edgeTangents[i];

        if (!midpoint || !length || !tangent) return;

        ctx.fillStyle = 'red';
        ctx.beginPath();
        ctx.arc(midpoint[0], midpoint[1], 2 / zoom, 0, 2 * Math.PI);
        ctx.fill();

        const screenX = midpoint[0] * zoom + offsetX;
        const screenY = midpoint[1] * zoom + offsetY;

        ctx.save();
        ctx.setTransform(1, 0, 0, 1, 0, 0);

        ctx.fillStyle = 'blue';
        ctx.font = '10px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'bottom';
        ctx.fillText(`${i}L ${length.toFixed(6)}`, screenX, screenY - 15);
        ctx.fillText(`T ${tangent.map((t) => t.toFixed(1))}`, screenX, screenY - 5);
        ctx.fillText(`${edge[0]}, ${edge[1]}`, screenX, screenY + 15);

        ctx.strokeStyle = 'green';
        ctx.lineWidth = 1.5;
        const tangentScreenX = tangent[0] * zoom;
        const tangentScreenY = tangent[1] * zoom;
        const tangentMagnitude = Math.sqrt(tangentScreenX * tangentScreenX + tangentScreenY * tangentScreenY) || 1;
        const dirX = tangentScreenX / tangentMagnitude;
        const dirY = tangentScreenY / tangentMagnitude;
        drawArrow(ctx, screenX, screenY, dirX, dirY, 20);

        ctx.restore();
    });
}

function drawSubvertices(ctx, subvertices) {
	const { zoom } = get(canvasTransform);

	ctx.fillStyle = 'black';
	for (const subvertex of subvertices) {
		const [x, y] = subvertex.position; // Use precomputed position from updateKernelState
		ctx.beginPath();
		ctx.arc(x, y, 4 / zoom, 0, 2 * Math.PI);
		ctx.fill();
	}
}

export function drawKernelMatrix(kernelCanvas, kernelMatrix) {
	if (!kernelCanvas || !kernelMatrix || !math.isMatrix(kernelMatrix)) return;

	const ctx = kernelCanvas.getContext('2d');
	const size = kernelMatrix.size();
	const padding = 50;
	const maxWidth = Math.min(window.innerWidth / 2, 800);
	const maxContentWidth = maxWidth - padding * 2;

	const cellSize = Math.min(
		50,
		Math.floor(maxContentWidth / size[1]),
		Math.floor(maxContentWidth / size[0])
	);

	const canvasWidth = Math.min(maxWidth, padding * 2 + size[1] * cellSize);
	const canvasHeight = padding * 2 + size[0] * cellSize;

	kernelCanvas.width = canvasWidth;
	kernelCanvas.height = canvasHeight;
	kernelCanvas.style.width = `${canvasWidth}px`;
	kernelCanvas.style.height = `${canvasHeight}px`;

	ctx.clearRect(0, 0, canvasWidth, canvasHeight);

	const maxVal = math.max(kernelMatrix) || 1;

	for (let i = 0; i < size[0]; i++) {
		for (let j = 0; j < size[1]; j++) {
			const value = kernelMatrix.get([i, j]);
			const intensity = value / maxVal;
			const r = Math.round(255 - intensity * 255);
			const g = Math.round(255 - intensity * 255);
			const b = 255;
			const minOpacity = 0.1;
			const opacity = minOpacity + (1 - minOpacity) * intensity;
			ctx.fillStyle = `rgba(${r}, ${g}, ${b}, ${opacity})`;
			ctx.fillRect(padding + j * cellSize, padding + i * cellSize, cellSize, cellSize);
		}
	}

	ctx.strokeStyle = 'black';
	ctx.lineWidth = 1;
	for (let i = 0; i <= size[0]; i++) {
		ctx.beginPath();
		ctx.moveTo(padding, padding + i * cellSize);
		ctx.lineTo(padding + size[1] * cellSize, padding + i * cellSize);
		ctx.stroke();
	}
	for (let j = 0; j <= size[1]; j++) {
		ctx.beginPath();
		ctx.moveTo(padding + j * cellSize, padding);
		ctx.lineTo(padding + j * cellSize, padding + size[0] * cellSize);
		ctx.stroke();
	}

	ctx.fillStyle = 'black';
	ctx.font = '12px Arial';
	ctx.textAlign = 'center';
	for (let j = 0; j < size[1]; j++) {
		ctx.fillText(j.toString(), padding + j * cellSize + cellSize / 2, padding - 10);
	}

	ctx.textAlign = 'right';
	for (let i = 0; i < size[0]; i++) {
		ctx.fillText(i.toString(), padding - 5, padding + i * cellSize + cellSize / 2);
	}
}


---
File: /src/lib/graphstate.js
---

// src/lib/graphState.js
import {
	calculateEdgeProperties,
	calculateDisjointEdgePairs,
	calculateDiscreteKernel,
	calculateDiscreteEnergy
} from '$lib/energyCalculations';
import { generateSubvertices } from '$lib/graphUtils';
import { subvertices } from '$lib/stores';

export function initializeKernelState(vertices, edges, alpha, beta) {
	const disjointPairs = calculateDisjointEdgePairs(edges);
	const edgeProps = calculateEdgeProperties(vertices, edges);
	const kernelMatrix = calculateDiscreteKernel(
		vertices,
		edges,
		edgeProps.edgeTangents,
		alpha,
		beta,
		disjointPairs
	);
	const discreteEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
	const newSubvertices = generateSubvertices(vertices, edges);

	subvertices.set(newSubvertices);

	return {
		kernelMatrix,
		discreteEnergy,
		edgeProps,
		disjointPairs
	};
}

export function updateKernelState(vertices, edges, alpha, beta, disjointPairs) {
	const edgeProps = calculateEdgeProperties(vertices, edges);
	const kernelMatrix = calculateDiscreteKernel(
		vertices,
		edges,
		edgeProps.edgeTangents,
		alpha,
		beta,
		disjointPairs
	);
	const discreteEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
	const newSubvertices = generateSubvertices(vertices, edges);

	subvertices.set(newSubvertices);

	return {
		kernelMatrix,
		discreteEnergy,
		edgeProps
	};
}


---
File: /src/lib/graphUtils.js
---

// src/lib/graphUtils.js
import { get } from 'svelte/store';
import { config } from './stores';

export function getColor(vertexEnergy, vertexData) {
	const maxEnergy = Math.max(...vertexData.map((v) => v.energy), 1e-9);
	const normalizedEnergy = vertexEnergy / maxEnergy;
	const r = Math.floor(255 * normalizedEnergy);
	const g = 0;
	const b = Math.floor(255 * (1 - normalizedEnergy));
	return `rgb(${r}, ${g}, ${b})`;
}

export function getRadius(vertexEnergy, vertexData) {
	const maxEnergy = Math.max(...vertexData.map((v) => v.energy), 1e-9);
	const normalizedEnergy = vertexEnergy / maxEnergy;
	return 2 + 8 * normalizedEnergy;
}

export function drawArrow(ctx, fromX, fromY, dirX, dirY, length) {
	const headLength = 7;
	const headAngle = Math.PI / 6;

	const toX = fromX + dirX * length;
	const toY = fromY + dirY * length;

	ctx.beginPath();
	ctx.moveTo(fromX, fromY);
	ctx.lineTo(toX, toY);
	ctx.stroke();

	const angle = Math.atan2(dirY, dirX);
	ctx.beginPath();
	ctx.moveTo(toX, toY);
	ctx.lineTo(
		toX - headLength * Math.cos(angle - headAngle),
		toY - headLength * Math.sin(angle - headAngle)
	);
	ctx.moveTo(toX, toY);
	ctx.lineTo(
		toX - headLength * Math.cos(angle + headAngle),
		toY - headLength * Math.sin(angle + headAngle)
	);
	ctx.stroke();
}

export function generateSubvertices(vertices, edges) {
	const subvertices = [];
	const { subvertexGap } = get(config);

	for (const edge of edges) {
		const [v1, v2] = edge;
		const p1 = vertices[v1];
		const p2 = vertices[v2];
		const dx = p2[0] - p1[0];
		const dy = p2[1] - p1[1];
		const edgeLength = Math.sqrt(dx * dx + dy * dy);
		const subvertexCount = Math.max(1, Math.floor(edgeLength / subvertexGap));

		for (let i = 1; i <= subvertexCount; i++) {
			const t = i / (subvertexCount + 1);
			const x = p1[0] + t * dx;
			const y = p1[1] + t * dy;
			subvertices.push({ position: [x, y], edge });
		}
	}
	return subvertices;
}

export function generateRandomGraph(width, height) {
	const numVertices = Math.floor(Math.random() * 10) + 5;
	const vertices = [];
	const edges = [];

	for (let i = 0; i < numVertices; i++) {
		vertices.push([Math.random() * width, Math.random() * height]);
	}

	const edgeSet = new Set();
	for (let i = 0; i < numVertices; i++) {
		for (let j = i + 1; j < numVertices; j++) {
			if (Math.random() < 0.1) {
				const edge = [i, j];
				const edgeString = `${i}-${j}`;
				if (!edgeSet.has(edgeString)) {
					edges.push(edge);
					edgeSet.add(edgeString);
				}
			}
		}
	}

	const subvertices = generateSubvertices(vertices, edges);
	return { vertices, edges, subvertices };
}

export function generateBipartiteGraph(width, height) {
	const vertices = [];
	const edges = [];

	for (let i = 0; i < 3; i++) {
		vertices.push([width * 0.3, height * (0.2 + i * 0.3)]);
	}

	for (let i = 0; i < 3; i++) {
		vertices.push([width * 0.7, height * (0.2 + i * 0.3)]);
	}

	for (let i = 0; i < 3; i++) {
		for (let j = 0; j < 3; j++) {
			edges.push([i, 3 + j]);
		}
	}

	const subvertices = generateSubvertices(vertices, edges);
	return { vertices, edges, subvertices };
}

export function generate2x3BipartiteGraph(width, height) {
	const vertices = [];
	const edges = [];

	for (let i = 0; i < 2; i++) {
		vertices.push([width * 0.3, height * (0.3 + i * 0.4)]);
	}

	for (let i = 0; i < 3; i++) {
		vertices.push([width * 0.7, height * (0.2 + i * 0.3)]);
	}

	for (let i = 0; i < 2; i++) {
		for (let j = 0; j < 3; j++) {
			edges.push([i, 2 + j]);
		}
	}

	const subvertices = generateSubvertices(vertices, edges);
	return { vertices, edges, subvertices };
}


---
File: /src/lib/innerProduct.js
---

// src/lib/innerProduct.js
import * as math from 'mathjs';
import {
	calculateEdgeProperties,
	calculateDisjointEdgePairs,
	tangentPointKernel
} from '$lib/energyCalculations';
import { get } from 'svelte/store';
import { config } from '$lib/stores';

/**
 * Builds the weight matrices W and W0 used in energy calculations.
 * @param {number} alpha - Energy parameter.
 * @param {number} beta - Energy parameter.
 * @param {Array<Array<number>>} edges - Array of edges, where each edge is an array of two vertex indices.
 * @param {Array<Array<number>>} disjointPairs - Array of disjoint edge pairs.
 * @param {Array<Array<number>>} vertices - Array of vertex coordinates.
 * @param {Array<Array<number>>} edgeTangents - Array of edge tangent vectors.
 * @param {Array<number>} edgeLengths - Array of edge lengths.
 * @returns {{W: math.Matrix, W0: math.Matrix}} - The weight matrices W and W0.
 */
function build_weights(alpha, beta, edges, disjointPairs, vertices, edgeTangents, edgeLengths) {
	const edge_num = edges.length;
	const s = (beta - 1) / alpha; // Fractional order s from the paper
	const sigma = s - 1; // Correct sigma = s - 1 (Section 4.2.2)

	const W = math.zeros(edge_num, edge_num);
	const W0 = math.zeros(edge_num, edge_num);

	for (let I = 0; I < disjointPairs.length; I++) {
		for (const J of disjointPairs[I]) {
			let elt1 = 0;
			let elt2 = 0;

			for (let a = 0; a < 2; a++) {
				for (let b = 0; b < 2; b++) {
					const i = edges[I][a];
					const j = edges[J][b];
					const p = vertices[i];
					const q = vertices[j];
					const diff = math.subtract(p, q);
					const epsilon = get(config).epsilonStability; // Use config epsilon
					let diff_norm = math.norm(diff) + epsilon; // Prevent division by zero

					const term1 = 1 / Math.pow(diff_norm, 2 * sigma + 1);
					elt1 += term1;

					// Constants for B0 calculation
					const alph = 2;
					const bet = 4;

					const cross = math.det([diff, edgeTangents[I]]); // 2D cross product
					const cross_norm = Math.abs(cross); // Use absolute value for 2D

					const k_numerator = Math.pow(cross_norm, alph);
					const k_denominator = Math.pow(diff_norm, bet);
					const k = k_numerator / k_denominator;

					const term2 = k / Math.pow(diff_norm, 2 * sigma + 1);
					elt2 += term2;
				}
			}

			const w_ij_factor = 0.25 * edgeLengths[I] * edgeLengths[J];
			W.set([I, J], w_ij_factor * elt1);
			W0.set([I, J], w_ij_factor * elt2);
		}
	}

	return { W, W0 };
}

/**
 * Calculates the low-order term matrix B0.
 * @param {Array<Array<number>>} vertices - Array of vertex coordinates.
 * @param {Array<Array<number>>} edges - Array of edges.
 * @param {math.Matrix} W0 - The weight matrix W0.
 * @returns {math.Matrix} - The low-order term matrix B0.
 */
function calculateLowOrderTerm(vertices, edges, W0) {
	const numVertices = vertices.length;
	const B0 = math.zeros(numVertices, numVertices);
	const disjointEdges = calculateDisjointEdgePairs(edges);

	for (let I = 0; I < edges.length; I++) {
		for (const J of disjointEdges[I]) {
			const w_IJ_0 = W0.get([I, J]);

			for (let a = 0; a < 2; a++) {
				for (let b = 0; b < 2; b++) {
					const i_a = edges[I][a];
					const i_b = edges[I][b];
					const j_a = edges[J][a];
					const j_b = edges[J][b];

					B0.set([i_a, i_b], B0.get([i_a, i_b]) + 0.25 * w_IJ_0);
					B0.set([j_a, j_b], B0.get([j_a, j_b]) + 0.25 * w_IJ_0);
					B0.set([i_a, j_b], B0.get([i_a, j_b]) - 0.25 * w_IJ_0);
					B0.set([j_a, i_b], B0.get([j_a, i_b]) - 0.25 * w_IJ_0);
				}
			}
		}
	}
	console.log('Low order term B0:', B0.toArray());
	return B0;
}

/**
 * Calculates the high-order term matrix B.
 * @param {Array<Array<number>>} vertices - Array of vertex coordinates.
 * @param {Array<Array<number>>} edges - Array of edges.
 * @param {math.Matrix} W - The weight matrix W.
 * @param {Array<number>} edgeLengths - Array of edge lengths.
 * @param {Array<Array<number>>} edgeTangents - Array of edge tangent vectors.
 * @returns {math.Matrix} - The high-order term matrix B.
 */
function calculateHighOrderTerm(vertices, edges, W, edgeLengths, edgeTangents) {
	const numVertices = vertices.length;
	const B = math.zeros(numVertices, numVertices);
	const disjointEdges = calculateDisjointEdgePairs(edges);

	for (let I = 0; I < edges.length; I++) {
		for (const J of disjointEdges[I]) {
			const l_I = edgeLengths[I];
			const l_J = edgeLengths[J];
			const T_I = edgeTangents[I];
			const T_J = edgeTangents[J];
			const w_IJ = W.get([I, J]);
			const dot_TI_TJ = math.dot(T_I, T_J); // 2D dot product

			for (let a = 0; a < 2; a++) {
				for (let b = 0; b < 2; b++) {
					const sign = Math.pow(-1, a + b);
					const i_a = edges[I][a];
					const i_b = edges[I][b];
					const j_a = edges[J][a];
					const j_b = edges[J][b];

					const val_1 = (sign * w_IJ) / (l_I * l_I);
					const val_2 = (sign * w_IJ) / (l_J * l_J);
					const val_3 = (sign * w_IJ * dot_TI_TJ) / (l_I * l_J);

					B.set([i_a, i_b], B.get([i_a, i_b]) + val_1);
					B.set([j_a, j_b], B.get([j_a, j_b]) + val_2);
					B.set([i_a, j_b], B.get([i_a, j_b]) - val_3);
					B.set([j_a, i_b], B.get([j_a, i_b]) - val_3);
				}
			}
		}
	}

	console.log('High order term B:', B.toArray());
	return B;
}

/**
 * Calculates the discrete inner product matrix A and its components.
 * @param {Array<Array<number>>} vertices - Array of vertex coordinates.
 * @param {Array<Array<number>>} edges - Array of edges.
 * @param {number} alpha - Energy parameter.
 * @param {number} beta - Energy parameter.
 * @returns {{A_reg: math.Matrix, B0: math.Matrix, B: math.Matrix}} - The regularized inner product matrix A_reg and component matrices B0 and B.
 */
export function calculateDiscreteInnerProduct(vertices, edges, alpha, beta) {
	console.log('Calculating discrete inner product with alpha:', alpha, 'beta:', beta);
	const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);

	const { W, W0 } = build_weights(
		alpha,
		beta,
		edges,
		calculateDisjointEdgePairs(edges),
		vertices,
		edgeTangents,
		edgeLengths
	);
	const B = calculateHighOrderTerm(vertices, edges, W, edgeLengths, edgeTangents);
	const B0 = calculateLowOrderTerm(vertices, edges, W0);

	let A = math.add(B0, B);
	console.log('Initial A Matrix (B0 + B):', A.toArray());
	const A_reg = A;

	return { A_reg, B0, B };
}

export function build_A_bar_2D(alpha, beta, vertices, edges) {
	const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);
	const disjointEdges = calculateDisjointEdgePairs(edges);
	const numVertices = vertices.length;

	const { W, W0 } = build_weights(
		alpha,
		beta,
		edges,
		disjointEdges,
		vertices,
		edgeTangents,
		edgeLengths
	);
	const B = calculateHighOrderTerm(vertices, edges, W, edgeLengths, edgeTangents);
	const B0 = calculateLowOrderTerm(vertices, edges, W0);
	let A = math.add(B, B0);

	// Regularization to ensure invertibility (use config epsilon)
	const epsilon = get(config).epsilonStability;
	const reg = math.multiply(epsilon, math.identity(numVertices));
	A = math.add(A, reg);

	// Build 2D block-diagonal A_bar
	const A_bar = math.zeros(2 * numVertices, 2 * numVertices);
	A_bar.subset(math.index(math.range(0, numVertices), math.range(0, numVertices)), A);
	A_bar.subset(
		math.index(math.range(numVertices, 2 * numVertices), math.range(numVertices, 2 * numVertices)),
		A
	);

	return { A_bar, B, B0 };
}

export function computePreconditionedGradient(
	vertices,
	edges,
	edgeTangents,
	alpha,
	beta,
	differential
) {
	console.log('Computing preconditioned gradient with differential:', differential);
	const numVertices = vertices.length;
	const { A_bar } = build_A_bar_2D(alpha, beta, vertices, edges);
	const differentialFlat = differential.flat();

	let gradFlat;
	try {
		gradFlat = math.lusolve(A_bar, differentialFlat);
		console.log('Preconditioned gradient (flat):', gradFlat.toArray());
	} catch (e) {
		console.error('Linear solve failed:', e);
		throw new Error('Failed to compute preconditioned gradient due to singular matrix');
	}

	const gradArray = gradFlat.toArray();
	const grad = [];
	for (let i = 0; i < numVertices; i++) {
		grad[i] = [gradArray[i * 2], gradArray[i * 2 + 1]];
	}
	console.log('Gradient:', grad);
	return grad;
}



---
File: /src/lib/interaction.js
---

// src/lib/interaction.js
import { get, writable } from 'svelte/store';
import { canvasTransform } from '$lib/stores';

export function setupInteractions(canvas, vertices, updateFn) {
	let draggingVertex = null;
	let dragOffsetX = 0;
	let dragOffsetY = 0;
	let isDraggingCanvas = false; // For panning the canvas
	let lastMouseX = 0;
	let lastMouseY = 0;
	let isSpacePressed = false; // Track Space key state
	let isCtrlPressed = false; // Track Ctrl key state

	const MIN_ZOOM = 0.1; // Minimum zoom level
	const MAX_ZOOM = 5.0; // Maximum zoom level

	function getWorldCoords(screenX, screenY) {
		const { offsetX, offsetY, zoom } = get(canvasTransform);
		return {
			worldX: (screenX - offsetX) / zoom,
			worldY: (screenY - offsetY) / zoom
		};
	}

	function getScreenCoords(worldX, worldY) {
		const { offsetX, offsetY, zoom } = get(canvasTransform);
		return {
			screenX: worldX * zoom + offsetX,
			screenY: worldY * zoom + offsetY
		};
	}

	function handleMouseDown(event) {
		const rect = canvas.getBoundingClientRect();
		const mouseX = event.clientX - rect.left;
		const mouseY = event.clientY - rect.top;

		// Check if Space is pressed for canvas dragging (panning)
		if (isSpacePressed) {
			isDraggingCanvas = true;
			lastMouseX = mouseX;
			lastMouseY = mouseY;
			return; // Skip vertex dragging if panning
		}

		// Check for vertex dragging (existing behavior)
		const { zoom } = get(canvasTransform);
		for (let i = 0; i < vertices.length; i++) {
			const [vx, vy] = vertices[i];
			const { screenX, screenY } = getScreenCoords(vx, vy);
			const distance = Math.sqrt((mouseX - screenX) ** 2 + (mouseY - screenY) ** 2);
			if (distance <= 10 / zoom) {
				// Increased click radius for better hit detection
				draggingVertex = i;
				const worldCoords = getWorldCoords(mouseX, mouseY);
				dragOffsetX = vx - worldCoords.worldX;
				dragOffsetY = vy - worldCoords.worldY;
				break;
			}
		}
	}

	function handleMouseMove(event) {
		const rect = canvas.getBoundingClientRect();
		const mouseX = event.clientX - rect.left;
		const mouseY = event.clientY - rect.top;

		if (isDraggingCanvas && isSpacePressed) {
			// Panning: Update canvas offset with throttling for performance
			const dx = mouseX - lastMouseX;
			const dy = mouseY - lastMouseY;
			canvasTransform.update((transform) => ({
				...transform,
				offsetX: transform.offsetX + dx,
				offsetY: transform.offsetY + dy
			}));
			lastMouseX = mouseX;
			lastMouseY = mouseY;
			requestAnimationFrame(updateFn); // Use requestAnimationFrame for smoother updates
			return; // Skip vertex dragging if panning
		}

		if (draggingVertex !== null) {
			// Vertex dragging (optimized for performance)
			const worldCoords = getWorldCoords(mouseX, mouseY);
			let newX = worldCoords.worldX + dragOffsetX;
			let newY = worldCoords.worldY + dragOffsetY;

			// Keep vertices within canvas bounds (optional, adjust as needed)
			const { zoom } = get(canvasTransform);
			newX = Math.max(0, Math.min(canvas.width / zoom, newX));
			newY = Math.max(0, Math.min(canvas.height / zoom, newY));

			vertices[draggingVertex] = [newX, newY];
			requestAnimationFrame(updateFn); // Use requestAnimationFrame for smoother updates
		}
	}

	function handleMouseUp() {
		draggingVertex = null;
		isDraggingCanvas = false;
	}

	function handleMouseWheel(event) {
		if (isCtrlPressed) {
			event.preventDefault();
			const rect = canvas.getBoundingClientRect();
			const mouseX = event.clientX - rect.left;
			const mouseY = event.clientY - rect.top;

			// Calculate zoom factor
			const zoomFactor = event.deltaY < 0 ? 1.1 : 0.9; // Zoom in/out by 10%
			canvasTransform.update((transform) => {
				const newZoom = Math.min(MAX_ZOOM, Math.max(MIN_ZOOM, transform.zoom * zoomFactor));
				const worldCoords = getWorldCoords(mouseX, mouseY);
				const newOffsetX = mouseX - worldCoords.worldX * newZoom;
				const newOffsetY = mouseY - worldCoords.worldY * newZoom;
				return {
					offsetX: newOffsetX,
					offsetY: newOffsetY,
					zoom: newZoom
				};
			});

			requestAnimationFrame(updateFn); // Use requestAnimationFrame for smoother updates
		}
	}

	function handleKeyDown(event) {
		if (event.code === 'Space') {
			isSpacePressed = true;
		} else if (event.code === 'ControlLeft' || event.code === 'ControlRight') {
			isCtrlPressed = true;
		}
	}

	function handleKeyUp(event) {
		if (event.code === 'Space') {
			isSpacePressed = false;
			isDraggingCanvas = false; // Stop panning when Space is released
		} else if (event.code === 'ControlLeft' || event.code === 'ControlRight') {
			isCtrlPressed = false;
		}
	}

	// Add event listeners
	canvas.addEventListener('mousedown', handleMouseDown);
	canvas.addEventListener('mousemove', handleMouseMove);
	canvas.addEventListener('mouseup', handleMouseUp);
	canvas.addEventListener('mouseleave', handleMouseUp);
	canvas.addEventListener('wheel', handleMouseWheel, { passive: false });
	document.addEventListener('keydown', handleKeyDown);
	document.addEventListener('keyup', handleKeyUp);

	// Return a cleanup function
	return () => {
		canvas.removeEventListener('mousedown', handleMouseDown);
		canvas.removeEventListener('mousemove', handleMouseMove);
		canvas.removeEventListener('mouseup', handleMouseUp);
		canvas.removeEventListener('mouseleave', handleMouseUp);
		canvas.removeEventListener('wheel', handleMouseWheel);
		document.removeEventListener('keydown', handleKeyDown);
		document.removeEventListener('keyup', handleKeyUp);
	};
}



---
File: /src/lib/optimization.js
---

// src/lib/optimization.js
import {
    calculateDifferential,
    calculateEdgeProperties,
    calculateDiscreteEnergy
} from '$lib/energyCalculations';
import { computePreconditionedGradient } from '$lib/innerProduct';
import * as math from 'mathjs';
import { get } from 'svelte/store';
import { config } from '$lib/stores';
import { updateKernelState } from '$lib/graphState'; // Import to update subvertices

// Configuration toggles
const usePreconditioned = false;
const applyProjectConstraints = false;
const applyBarycenter = true;

function l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths) {
    const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);
    const gradient = differential.map(([dx, dy]) => [dx, dy]);
    
    const stepSize = get(config).l2StepSize;
    
    const newVertices = vertices.map((vertex, i) => [
        vertex[0] - stepSize * gradient[i][0],
        vertex[1] - stepSize * gradient[i][1]
    ]);

    let constrainedVertices = [...newVertices];
    if (applyProjectConstraints) {
        constrainedVertices = projectConstraints(constrainedVertices, edges, initialEdgeLengths);
    }
    
    if (applyBarycenter) {
        enforceBarycenter(constrainedVertices);
    }

    return constrainedVertices;
}

function preconditionedGradientDescentStep(
    vertices,
    edges,
    alpha,
    beta,
    disjointPairs,
    initialEdgeLengths
) {
    const { edgeTangents } = calculateEdgeProperties(vertices, edges);
    
    const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);
    
    let gradient;
    try {
        gradient = computePreconditionedGradient(
            vertices,
            edges,
            edgeTangents,
            alpha,
            beta,
            differential
        );
    } catch (e) {
        console.warn('Preconditioned gradient failed, falling back to L2 gradient:', e);
        return l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths);
    }

    const d = gradient.map(([gx, gy]) => [-gx, -gy]);
    const d_norm = Math.sqrt(d.flat().reduce((sum, val) => sum + val * val, 0)) || 1;
    const d_normalized = d.map(([dx, dy]) => [dx / d_norm, dy / d_norm]);

    const differentialFlat = differential.flat();
    const dFlat = d_normalized.flat();
    const slope = math.dot(differentialFlat, dFlat);
    
    const E_old = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
    
    const a_const = get(config).aConst;
    const b_const = get(config).bConst;
    const max_line_search = get(config).maxLineSearch;
    let t = get(config).tauInitial;
    const stepSize = get(config).precondStepSize;

    let vertices_new = [...vertices];
    for (let i = 0; i < max_line_search; i++) {
        vertices_new = vertices.map((vertex, idx) => [
            vertex[0] + t * d_normalized[idx][0] * stepSize,
            vertex[1] + t * d_normalized[idx][1] * stepSize
        ]);
        
        if (applyProjectConstraints) {
            vertices_new = projectConstraints(vertices_new, edges, initialEdgeLengths);
        }
        
        if (applyBarycenter) {
            enforceBarycenter(vertices_new);
        }
        
        const E_new = calculateDiscreteEnergy(vertices_new, edges, alpha, beta, disjointPairs);
        console.log(`Line search iteration ${i}: t=${t}, E_new=${E_new}, E_old=${E_old}, condition=${E_old + a_const * t * slope}`);
        if (E_new <= E_old + a_const * t * slope) {
            console.log('Line search converged at t=', t, 'with energy reduction:', E_old - E_new);
            return vertices_new;
        }
        
        t *= b_const;
    }

    console.warn('Line search did not converge, using smallest step size');
    return vertices_new;
}

function projectConstraints(
    vertices,
    edges,
    initialEdgeLengths,
    maxIterations = get(config).maxConstraintIterations,
    tolerance = get(config).constraintTolerance
) {
    const projectedVertices = vertices.map(v => [...v]);
    
    for (let iter = 0; iter < maxIterations; iter++) {
        let maxError = 0;
        
        for (let i = 0; i < edges.length; i++) {
            const [v1Idx, v2Idx] = edges[i];
            const v1 = projectedVertices[v1Idx];
            const v2 = projectedVertices[v2Idx];
            
            const dx = v2[0] - v1[0];
            const dy = v2[1] - v1[1];
            const currentLength = Math.sqrt(dx * dx + dy * dy) + get(config).epsilonStability;

            const targetLength = initialEdgeLengths[i];
            const error = (currentLength - targetLength) / targetLength;
            
            maxError = Math.max(maxError, Math.abs(error));
            
            if (Math.abs(error) > tolerance) {
                const correction = error / 2;
                const scaleFactor = 1 - correction;
                
                const midX = (v1[0] + v2[0]) / 2;
                const midY = (v1[1] + v2[1]) / 2;
                
                const halfDx = dx / 2;
                const halfDy = dy / 2;
                
                projectedVertices[v1Idx][0] = midX - halfDx * scaleFactor;
                projectedVertices[v1Idx][1] = midY - halfDy * scaleFactor;
                projectedVertices[v2Idx][0] = midX + halfDx * scaleFactor;
                projectedVertices[v2Idx][1] = midY + halfDy * scaleFactor;
            }
        }
        
        if (maxError < tolerance) {
            console.log(`Constraint projection converged after ${iter + 1} iterations`);
            break;
        }
    }
    
    return projectedVertices;
}

function enforceBarycenter(vertices, options = {}) {
    const {
        targetBarycenter = null,
        centerAtOrigin = false
    } = options;
    
    const currentBarycenter = [0, 0];
    const n = vertices.length;
    
    for (const vertex of vertices) {
        currentBarycenter[0] += vertex[0] / n;
        currentBarycenter[1] += vertex[1] / n;
    }
    
    let target;
    if (centerAtOrigin) {
        target = [0, 0];
    } else if (targetBarycenter) {
        target = targetBarycenter;
    } else {
        target = currentBarycenter;
    }
    
    const dx = target[0] - currentBarycenter[0];
    const dy = target[1] - currentBarycenter[1];
    
    for (const vertex of vertices) {
        vertex[0] += dx;
        vertex[1] += dy;
    }
    
    return vertices;
}

export function gradientDescentStep(
    vertices,
    edges,
    alpha,
    beta,
    disjointPairs,
    initialEdgeLengths
) {
    return usePreconditioned
        ? preconditionedGradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths)
        : l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths);
}

export function createOptimizer(
    vertices,
    edges,
    alpha,
    beta,
    disjointPairs,
    maxIterations,
    onUpdate,
    initialEdgeLengths
) {
    if (typeof onUpdate !== 'function') throw new Error('onUpdate must be a function');

    let currentIteration = 0;
    let intervalId = null;
    let lastEnergy = null;
    let stuckCounter = 0;
    
    const applyPerturbation = get(config).applyPerturbation;
    if (applyPerturbation) {
        lastEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
    }

    const optimizer = {
        step: () => {
            if (currentIteration < maxIterations) {
                const newVertices = gradientDescentStep(
                    vertices,
                    edges,
                    alpha,
                    beta,
                    disjointPairs,
                    initialEdgeLengths
                );
                vertices.forEach((v, i) => {
                    v[0] = newVertices[i][0];
                    v[1] = newVertices[i][1];
                });
                
                // Update kernel state including subvertices
                updateKernelState(vertices, edges, alpha, beta, disjointPairs);
                
                if (applyPerturbation) {
                    const newEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
                    const energyChange = newEnergy - lastEnergy;
                    
                    if (Math.abs(energyChange) < get(config).minEnergyChange) {
                        stuckCounter++;
                        
                        if (stuckCounter > get(config).maxStuckIterations) {
                            console.log('Optimizer stuck, applying random perturbation');
                            applyRandomPerturbation(vertices, get(config).perturbationScale);
                            stuckCounter = 0;
                        }
                    } else {
                        stuckCounter = 0;
                    }
                    
                    lastEnergy = newEnergy;
                }
                
                currentIteration++;
                onUpdate();
            } else {
                optimizer.stop();
            }
        },
        start: () => {
            currentIteration = 0;
            stuckCounter = 0;
            
            if (applyPerturbation) {
                lastEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
            }
            
            if (!intervalId) intervalId = setInterval(optimizer.step, 20);
        },
        stop: () => {
            if (intervalId) {
                clearInterval(intervalId);
                intervalId = null;
            }
        }
    };

    return optimizer;
}

function applyRandomPerturbation(vertices, scale) {
    for (const vertex of vertices) {
        vertex[0] += (Math.random() - 0.5) * scale;
        vertex[1] += (Math.random() - 0.5) * scale;
    }
}


---
File: /src/lib/stores.js
---

// src/lib/stores.js
import { writable, derived } from 'svelte/store';

export const config = writable({
	epsilonStability: 1e-7,
	epsilonKernel: 1e-6,
	finiteDiffH: 1e-4,
	constraintTolerance: 1e-7,
	tauInitial: 1.0,
	aConst: 0.1,
	bConst: 0.5,
	maxLineSearch: 20,
	differentialMethod: 'finiteDifference',
    precondStepSize: 20,
    l2StepSize: 100000,
    applyPerturbation: false,
    subvertexGap: 50 // New config: desired gap distance between subvertices (pixels)
});

export const vertices = writable([]); // Supervertices
export const edges = writable([]);
export const subvertices = writable([]); // New store for subvertices
export const kernelData = writable({
	kernelMatrix: null,
	discreteEnergy: 0,
	edgeProps: { edgeLengths: [], edgeTangents: [], edgeMidpoints: [] },
	disjointPairs: []
});
export const energyChange = writable(0);
export const previousEnergy = writable(0);

// Store for canvas transformations
export const canvasTransform = writable({
	offsetX: 0,
	offsetY: 0,
	zoom: 1.0
});

export const discreteEnergy = derived(kernelData, ($kernelData) => $kernelData.discreteEnergy);


---
File: /src/routes/+page.svelte
---

<script>
	import Controls from '$lib/Controls.svelte';
	import { onMount, onDestroy } from 'svelte';
	import { drawGraph, drawKernelMatrix } from '$lib/graphDrawing';
	import { createOptimizer } from '$lib/optimization';
	import {
		generateRandomGraph,
		generateBipartiteGraph,
		generate2x3BipartiteGraph
	} from '$lib/graphUtils';
	import { setupInteractions } from '$lib/interaction';
	import { initializeKernelState, updateKernelState } from '$lib/graphState';
	import {
		vertices,
		edges,
		subvertices,
		kernelData,
		energyChange,
		previousEnergy,
		discreteEnergy,
		config,
		canvasTransform
	} from '$lib/stores';
	import { get } from 'svelte/store';

	let graphCanvas;
	let kernelCanvas;
	let graphCtx;
	let optimizer;
	let cleanupInteractions = () => {};
	let isOptimizing = false;
	let graphType = 'bipartite';
	const width = 700;
	const height = 700;
	let alpha = 3;
	let beta = 6;
	const maxIterations = 1000;
	let initialEdgeLengths = [];

	onMount(() => {
		graphCtx = graphCanvas.getContext('2d');
		regenerateGraph();
	});

	onDestroy(() => {
		if (optimizer) optimizer.stop();
		cleanupInteractions();
	});

	function updateVisualization() {
		const updatedKernel = updateKernelState(
			$vertices,
			$edges,
			alpha,
			beta,
			$kernelData.disjointPairs
		);
		$kernelData = { ...updatedKernel, disjointPairs: $kernelData.disjointPairs };
		$energyChange = $discreteEnergy - $previousEnergy;
		$previousEnergy = $discreteEnergy;

		drawGraph(
			graphCtx,
			width,
			height,
			$vertices,
			$edges,
			updatedKernel.edgeProps,
			updatedKernel.kernelMatrix,
			alpha,
			beta,
			$kernelData.disjointPairs
		);

		drawKernelMatrix(kernelCanvas, updatedKernel.kernelMatrix);
	}

	function regenerateGraph() {
		$previousEnergy = 0;
		$energyChange = 0;
		let newVertices, newEdges, newSubvertices;
		if (graphType === 'random') {
			({
				vertices: newVertices,
				edges: newEdges,
				subvertices: newSubvertices
			} = generateRandomGraph(width, height));
		} else if (graphType === 'bipartite') {
			({
				vertices: newVertices,
				edges: newEdges,
				subvertices: newSubvertices
			} = generate2x3BipartiteGraph(width, height));
		}
		$vertices = newVertices;
		$edges = newEdges;
		$subvertices = newSubvertices;

		const initialKernel = initializeKernelState($vertices, $edges, alpha, beta);
		$kernelData = initialKernel;
		$previousEnergy = initialKernel.discreteEnergy;

		initialEdgeLengths = initialKernel.edgeProps.edgeLengths;

		optimizer = createOptimizer(
			$vertices,
			$edges,
			alpha,
			beta,
			initialKernel.disjointPairs,
			maxIterations,
			updateVisualization,
			initialEdgeLengths
		);

		cleanupInteractions();
		cleanupInteractions = setupInteractions(graphCanvas, $vertices, updateVisualization);

		canvasTransform.set({
			offsetX: 0,
			offsetY: 0,
			zoom: 1.0
		});
		updateVisualization();
	}

	function startOptimization() {
		if (optimizer && !isOptimizing) {
			console.log('Start optimization clicked');
			optimizer.start();
			isOptimizing = true;
		}
	}

	function stopOptimization() {
		if (optimizer && isOptimizing) {
			console.log('Stop optimization clicked');
			optimizer.stop();
			isOptimizing = false;
		}
	}

	function singleStep() {
		if (optimizer) {
			console.log('Single step clicked');
			optimizer.step();
		}
	}

	function updateAlphaBeta() {
		const updatedKernel = updateKernelState(
			$vertices,
			$edges,
			alpha,
			beta,
			$kernelData.disjointPairs
		);
		$kernelData = { ...updatedKernel, disjointPairs: $kernelData.disjointPairs };
		updateVisualization();
	}

	function updateConfig() {
		updateVisualization();
	}

	function getEnergyChangeColor() {
		return $energyChange < 0 ? 'green' : 'red';
	}
</script>

<div class="visualization-container">
	<div
		class="controls"
		style="display: flex; flex-direction: column; gap: 10px; top: 10px; left: 10px; z-index: 10;"
	>
		<button on:click={regenerateGraph}>Regenerate Graph</button>
		<button on:click={isOptimizing ? stopOptimization : startOptimization}>
			{isOptimizing ? 'Stop Optimization' : 'Start Optimization'}
		</button>
		<button on:click={singleStep}>Single Step</button>
		<div class="energy-value">
			<p>Discrete Energy: {$discreteEnergy.toFixed(4)}</p>
			<p style="color: {getEnergyChangeColor()}">Energy Change: {$energyChange.toFixed(4)}</p>
		</div>
		<Controls on:update={updateVisualization} />
	</div>
	<div class="graph-section">
		<div class="graph-container" style="position: relative; width: {width}px; height: {height}px;">
			<canvas bind:this={graphCanvas} {width} {height} style="position: absolute; top: 0; left: 0;"
			></canvas>
		</div>
	</div>
</div>

<style>
	.visualization-container {
		display: flex;
		flex-direction: row;
		gap: 20px;
	}
	.graph-section {
		flex: 1;
	}
</style>



Above is a cpp implementation of the paper by someone on gh, and then mine in js.

My implementation is for 2d. But generalize the code to have a choice whether to have a 2d or 3d canvas. Note that you will need to have to update kernel function to work for 2d and 3d (3d is cross product, whereas 2d is just a det), A_bar (for 3d) is 3|V | × 3|V |  matrix, instead of 2, and so on so forth. Also generalize graph creation, and other funcitons on your discretion. and of course Review the code from cpp and paper snippet carefully and apply it to 3d part of js. you must provide full code without simplifications.