# Elements of Functional Analysis
This section collects basic knowledge of weak formulations and functional spaces to understand the math behind the finite element method {cite}`salsa, Gazzola_A3, quarteroni_2016`.

Let $\Omega\subset\mathbb{R}^{d}$ with $d=2,3$ and $\partial \Omega$ its boundary.

## $L^2$ and Sobolev spaces

````{prf:definition}
:label: L2-def

Let $L^2(\Omega)$ be an Hilbert space of square integrable functions, i.e. 

$$
u\in L^2(\Omega)\Longleftrightarrow \int_\Omega |u|^2\,d\Omega <\infty
$$

Since $L^2$ is an Hilbert space, an inner product and an induced norm can be endowed

$$
(u, \, v) = \int_\Omega u\cdot v\,d\Omega \qquad \qquad
\|u\|_{L^2}^2 = \int_\Omega |u|^2\,d\Omega
$$
````

The same space can be easily extended to vector quantities and complex functions (the conjugate is introduced).

Let us introduced the notion of multi-index derivative  $D^{{\alpha}}$: given ${\alpha}\in\mathbb{N}^d$ of order $p=\sum_{i=1}^d \alpha_i$, the multi-index derivative of a function is defined as
\begin{equation}
    D^{{\alpha}}u = \frac{\partial^pu}{\partial x_1^{\alpha_1} \dots \partial x_d^{\alpha_d}},
\end{equation}
in which $d$ is usually 2 or 3 {cite}`Gazzola_A3`.

````{prf:definition}
:label: Hp-def

Let $\mathcal{H}^p(\Omega)$ be the Sobolev space of order $p$, i.e. 

$$
\mathcal{H}^p(\Omega) = \left\{u\in L^2(\Omega)\,:\,\int_\Omega \left|D^{{\alpha}}u\right|^2\,d\Omega <\infty \right\}
$$
````

All these derivatives are meant to be in the weak/distribution sense {cite}`Gazzola_A3`.

## Useful formulas
````{prf:theorem}
:label: Gauss-theo
Let $\mathbf{u}\in[\mathcal{H}^1(\Omega)]^d$ a vector function, the Gauss (divergence) theorem states:

$$
\int_\Omega \nabla \cdot \mathbf{u}\,d\Omega = \int_{\partial\Omega}\mathbf{u}\cdot \mathbf{n}\,d\sigma
$$
````

From this theorem, the following formulas can be derived
\begin{equation}
\int_\Omega (\nabla \cdot \mathbf{u}) \,p\,d\Omega = 
- \int_\Omega \mathbf{u} \cdot \nabla p\,d\Omega + \int_{\partial\Omega}(\mathbf{u}\cdot \mathbf{n})\,p\,d\sigma
\end{equation}
given $p\in\mathcal{H}^1(\Omega)$.

````{prf:proof}
Recalling that (using index/Einstein notation)

$$
\frac{\partial}{\partial x_i}\left(u_i\,p\right) = 
\frac{\partial\mathbf{u}_i}{\partial x_i} \,p + u_i \frac{\partial p}{\partial x_i}
$$
we can rewrite the formula as

$$
\int_\Omega \left[(\nabla \cdot \mathbf{u}) \,p + \mathbf{u} \cdot \nabla p\right] \,d\Omega = 
\int_\Omega \nabla \cdot (p\mathbf{u})\,d\Omega = \int_{\partial\Omega}(\mathbf{u}\cdot \mathbf{n})\,p\,d\sigma
$$
which is {prf:ref}`Gauss theorem <Gauss-theo>`.
````