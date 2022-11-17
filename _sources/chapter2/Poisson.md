# Poisson

The first problem that will be analysed is an elliptic equation, i.e. the Poisson problem
\begin{equation}
\left\{
    \begin{array}{ll}
        -\nabla \cdot \left(k\nabla T) = q''' & \mathbf{x}\in\Omega\\
        T = T_D & \mathbf{x}\in\Gamma_{D}\\
        -k\frac{\partial T}{\partial \mathbf{n}} = g_N & \mathbf{x}\in\Gamma_{N}
    \end{array}
\right.
\end{equation}
which can describe the temperature distribution inside a material.

## Derivation of the weak formulation
Let $\mathcal{V}\subset\mathcal{H}^1, \mathcal{V}_0\subset\mathcal{H}^1$ the trial and test space defined as

$$
\mathcal{V} = \left\{v\in\mathcal{H}^1:\;\left. v\right|_{\Gamma_D} = T_D\right\}\qquad 
\mathcal{V}_0 = \left\{v\in\mathcal{H}^1:\;\left. v\right|_{\Gamma_D} = 0\right\}
$$
Let $\varphi\in\mathcal{V}_0$ be the test function, then the strong problem is multiplied in scalar sense by it

$$
-\int_\Omega\nabla \cdot \left(k\nabla T\right)\,\varphi\,d\Omega = \int_\Omega q'''\,\varphi\,d\Omega
$$
we can apply the integration by parts to get the weak formulation
\begin{equation}
\int_\Omega k\nabla T\cdot\nabla\varphi\,d\Omega = \int_\Omega q'''\,\varphi\,d\Omega + \int_{\Gamma_N}g_N\,\varphi\,d\sigma
\end{equation}

## Algebraic system
The weak formulation can be then reduced to a linear system of equation by introducing the finite dimensional representation of the spaces $\left(T \simeq \sum_nt_n\phi_n\right)$

$$
A\mathbf{t} = \mathbf{f}
$$
It's worth noticing that the matrix $A$ is symmetric, thus a conjugate gradient solver can be used.
