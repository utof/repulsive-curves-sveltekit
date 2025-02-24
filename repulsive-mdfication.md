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



let me remind you this important part of the paper:
The discrete differential is then simply the partial derivatives of this energy with respect to the coordinates of all the curve vertices:

$$
d\hat{\mathcal{E}}^\alpha_\beta\big|_{\gamma} = \begin{bmatrix}
\partial \mathcal{E}^\alpha_\beta / \partial \gamma_1 & \cdots & \partial \mathcal{E}^\alpha_\beta / \partial \gamma_{|V|}
\end{bmatrix} \in \mathbb{R}^{3|V|}.
$$

These derivatives can be evaluated via any standard technique (*e.g.*, by hand, or using symbolic or automatic differentiation).

### 5.2 Discrete Inner Product

As in the smooth setting, we define our inner product matrix as a sum $A = B + B^0$ of high-order and low-order terms $B, B^0 \in \mathbb{R}^{|V| \times |V|}$ (as defined below).  For $\mathbb{R}^3$-valued functions, we also define a corresponding $3|V| \times 3|V|$ matrix

$$
\bar{A} =
\begin{bmatrix}
A & & \\
& A & \\
& & A
\end{bmatrix}. \qquad (19)
$$

Mirroring Equation 8, the discrete (fractional) Sobolev gradient $g \in \mathbb{R}^{3|V|}$ is then defined as the solution to the matrix equation

$$
\bar{A} g = d \hat{\mathcal{E}}_\beta^\alpha. \qquad (20)
$$

### 5.2.1 Discrete Derivative Operator.

For each edge \( I \in E \) we approximate the derivative \( \mathcal{D}u \) of a function \( u: M \to \mathbb{R} \) (Equation 11) via the finite difference formula \( \frac{1}{l_I}(u_{i_2} - u_{i_1})T_I \), where \( u_i \) denotes the value of \( u \) sampled at vertex \( i \). The corresponding derivative matrix \( D \in \mathbb{R}^{3|E| \times |V|} \) can be assembled from local \( 3 \times 2 \) matrices

\[
D_I = \frac{1}{l_I}
\begin{bmatrix}
-T_I & T_I
\end{bmatrix}.
\]

### 5.2.2 Discrete High-Order Term.

We approximate the high-order part of the inner product \( \langle \! \langle B_\sigma u, v \rangle \! \rangle \) as

\[
u^T B v = \sum_{I, J \in E, I \cap J = \emptyset} w_{IJ} \langle D_I u[I] - D_J u[J], D_I v[I] - D_J v[J] \rangle, \quad (21)
\]

where the weights \( w_{IJ} \) arise from applying trapezoidal quadrature to the denominator in Equation 25:

\[
w_{IJ} := \frac{1}{4} l_I l_J \sum_{i \in I} \sum_{j \in J} \frac{1}{|\gamma_i - \gamma_j|^{2\sigma + 1}}.
\]

The entries of the corresponding Gram matrix \( B \in \mathbb{R}^{|V| \times |V|} \) are obtained by differentiating Equation 21 with respect to the entries of \( u \) and \( v \). More explicitly, starting with the zero matrix one can build \( B \) by making the following increments for all pairs of disjoint edges \( I \cap J = \emptyset \), and all pairs of values \( a, b \in \{1, 2\} \):

\[
\begin{aligned}
B_{i_a i_b} += & (-1)^{a+b} w_{IJ} / l_I^2, & B_{i_a j_b} -= & (-1)^{a+b} w_{IJ} \langle T_I, T_J \rangle / (l_I l_J), \\
B_{j_a j_b} += & (-1)^{a+b} w_{IJ} / l_J^2, & B_{j_a i_b} -= & (-1)^{a+b} w_{IJ} \langle T_J, T_I \rangle / (l_J l_I). 
\end{aligned}
\]

### 5.2.3 Discrete Low-Order Term.

To discretize the low-order term \( B_\sigma^0 \) (Section 4.2.3), we use a different discrete weight

\[
w_{IJ}^0 := \frac{1}{4} l_I l_J \sum_{i \in I} \sum_{j \in J} \frac{k_4^2(\gamma_i, \gamma_j, T_I)}{|\gamma_i - \gamma_j|^{2\sigma + 1}},
\]

and define a matrix \( B^0 \in \mathbb{R}^{|V| \times |V|} \), given by the relationship

\[
u^T B^0 v = \sum_{I, J \in E, I \cap J = \emptyset} w_{IJ}^0 (u_I - u_J)(v_I - v_J).
\]

Following a similar derivation as above, this matrix can be constructed via the following increments:

\[
\begin{aligned}
B_{i_a i_b}^0 += & \frac{1}{4} w_{IJ}^0, & B_{i_a j_b}^0 -= & \frac{1}{4} w_{IJ}^0, \\
B_{j_a i_b}^0 -= & \frac{1}{4} w_{IJ}^0, & B_{j_a j_b}^0 += & \frac{1}{4} w_{IJ}^0.
\end{aligned}
\]

Above is a cpp implementation of the paper by someone on, and then mine in js.

Note that my implementation is for 2d only. Think thru and implement the needed step(s) in the js implementation and svelte5 based on the paper and the cpp code. You may refactor and optimize my code as well. Trust the given text and codes more than the paper to avoid vision and text incompatibilities.
Other known issues and tasks include:

Despite my attempts at implementation, even the simple 2x3 bipartite graph fails to minimize the energy, being stuck at the point of almost untangling. Figure out the issue. Even L2 gradient fails. For example a lot of times graphs (when random or 3x3 bipartite) get aligned along a long edge and then kinda get stuck there. I have disabled analytical difference because it's buggy but don't worry about it.
