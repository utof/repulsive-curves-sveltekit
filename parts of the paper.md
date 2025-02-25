let me remind you this important part of the paper:

Notation. In the discrete setting, we will model any collection of curves and loops (including several curves meeting at a common point) as a graph $G = (V,E)$ with vertex coordinates $\gamma : V \rightarrow \mathbb{R}^3$ (Figure 8); we use $|V|$ and $|E|$ to denote the number of vertices and edges, resp. For each edge $I \in E$ with endpoints $i_1, i_2$, we use

$\ell_I := |\gamma_{i_1} - \gamma_{i_2}|$, $T_I := (\gamma_{i_2} - \gamma_{i_1})/\ell_I$, and $x_I := (\gamma_{i_1} + \gamma_{i_2})/2$

to denote the edge length, unit tangent, and midpoint, resp. For any quantity $u: V \rightarrow \mathbb{R}$ on vertices we use $u_I := (u_{i_1} + u_{i_2})/2$ to denote the average value on edge $I = (i_1, i_2)$, and $u[I] := \begin{bmatrix} u_{i_1} & u_{i_2} \end{bmatrix}^T$ to denote the $2 \times 1$ column vector storing the values at its endpoints. Finally, we refer to any pair $(T,x) \in \mathbb{R}^6$ as a tangent-point.

5.1 Discrete Energy

Since the tangent-point energy is infinite for polygonal curves [Strzelecki and von der Mosel 2017, Figure 2.2], we assume that $\gamma$ is inscribed in some (unknown) smooth curve, and apply numerical quadrature to the smooth energy $\mathcal{E}_{\beta}^{\alpha}$. The resulting discrete energy then approximates the energy of any sufficiently smooth curve passing through the vertices $\gamma_i$. We start by integrating $k_{\beta}^{\alpha}$ over all pairs of edges:

$\sum_{I \in E} \sum_{J \in E} \int_{\overline{I}} \int_{\overline{J}} k_{\beta}^{\alpha}(\gamma(x), \gamma(y), T_I) dx_{\gamma} dy_{\gamma}.$ (16)

Here $\overline{I}$ denotes the interval along edge $I$. As given, this expression is ill-defined since two edges with a common endpoint contribute infinite energy. One idea is to instead use a term involving the curvature of the circle passing through the three distinct endpoints (in the spirit of Equation 1). However, such terms would contribute nothing to the energy in the limit of regular refinement (Figure 9) - hence, we simply omit neighboring edge pairs. Applying the (2D) trapezoidal rule to Equation 16 then yields a discrete energy

$\hat{\mathcal{E}}_{\beta}^{\alpha}(\gamma) = \sum_{I, J \in E, I \cap J = \emptyset} (\hat{k}_{\beta}^{\alpha})_{IJ} \ell_I \ell_J,$ (17)

where $\hat{k}$ is the discrete kernel

$(\hat{k}_{\beta}^{\alpha})_{IJ} := \frac{1}{4} \sum_{i \in I} \sum_{j \in J} k_{\beta}^{\alpha}(\gamma_i, \gamma_j, T_I)$. (18)

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