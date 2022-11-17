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

Let us introduced the notion of multi-index derivative  $D^{\bm{\alpha}}$: given $\bm{\alpha}\in\mathbb{N}^n$ of order $p=\sum_{i=1}^n \alpha_i$, the multi-index derivative of a function is defined as
\begin{equation}
    D^{\bm{\alpha}}u = \frac{\partial^pu}{\partial x_1^{\alpha_1} \dots \partial x_n^{\alpha_n}},
\end{equation}
in which $n$ is usually 2 or 3 {cite}`Gazzola_A3`.