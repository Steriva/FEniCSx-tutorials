# Steady Navier-Stokes

The Navier-Stokes equations describes the flows of incompressible viscous fluid 
\begin{equation}
\left\{
    \begin{array}{ll}
        \left(\mathbf{u}\cdot \nabla\right)\mathbf{u}-\nu\Delta \mathbf{u}+\nabla p = 0 & \mbox{in }\Omega\\
        \nabla\cdot \mathbf{u} = 0 & \mbox{in }\Omega\\
        \mathbf{u}=\mathbf{u}_D & \mbox{on }\Gamma_D\\
        \nu\frac{\partial \mathbf{u}}{\partial \mathbf{n}}-p\mathbf{n}=g\mathbf{n} & \mbox{on }\Gamma_N
    \end{array}
\right.
\end{equation}
This problem has a saddle point structure, making its numerical solution non-trivial, in addiction the convective term $\left(\mathbf{u}\cdot \nabla\right)\mathbf{u}$ is non-linear. Therefore, the numerical solver should include a Newton or fixed-point methods (or others).

Before entering into the details of the Galerkin problem let us derive the weak formulation.
Let $\mathcal{V}\subset[\mathcal{H}^1]^d,\; \mathcal{V}_0[\subset\mathcal{H}^1]^d$ the velocity trial and test spaces defined as

$$
\mathcal{V} = \left\{\mathbf{v}\in[\mathcal{H}^1(\Omega)]^d:\;\left. \mathbf{v}\right|_{\Gamma_D} = \mathbf{u}_D\right\}\qquad 
\mathcal{V}_0 = \left\{\mathbf{v}\in[\mathcal{H}^1(\Omega)]^d:\;\left. \mathbf{v}\right|_{\Gamma_D} = \mathbf{0}\right\}
$$
and let $\mathcal{Q}= L^2(\Omega)$ the pressure trial and test space. The momentum equation can be multiplied by the test function $\mathbf{v}\in\mathcal{V}_0$ and the integration by parts is applied

$$
\int_\Omega \left[\left(\mathbf{u}\cdot \nabla\right)\mathbf{u}\right]\cdot \mathbf{v}\,d\Omega+\int_\Omega \nu\nabla \mathbf{u}\cdot \nabla\mathbf{v}\,d\Omega -\int_\Omega p \nabla\cdot \mathbf{v}\,d\Omega - \int_{\partial\Omega} \left(\nu\frac{\partial \mathbf{u}}{\partial \mathbf{n}}-p\mathbf{n}\right)\cdot \mathbf{v}\,d\sigma=0
$$

whereas the continuity equation is multiplied by $q\in\mathcal{Q}$

$$
-\int_\Omega q \nabla\cdot \mathbf{u}\,d\Omega=0
$$

Imposing the boundary conditions, the weak formulation reads: *find $(\mathbf{u}, p)\in\mathcal{V}\times\mathcal{Q}$ s.t.*
\begin{equation}
\int_\Omega \left[\left(\mathbf{u}\cdot \nabla\right)\mathbf{u}\right]\cdot \mathbf{v}\,d\Omega+\int_\Omega \nu\nabla \mathbf{u}\cdot \nabla\mathbf{v}\,d\Omega -\int_\Omega p \nabla\cdot \mathbf{v}\,d\Omega - \int_\Omega q \nabla\cdot \mathbf{u}\,d\Omega = \int_{\Gamma_N} g\mathbf{n}\cdot \mathbf{v}\,d\sigma \qquad \forall (\mathbf{v}, q)\in\mathcal{V}_0\times\mathcal{Q}
\end{equation}

## Derivation of the linear system
When the finite dimensional spaces are introduced an important remark should be made. The Galerkin problem has a stable solution $(\mathbf{u}_h,p_h)$ if the finite dimensional spaces are *inf-sup* compatible. In fact, there exists a connection between the finite dimensional functional space of velocity and pressure referred to as the Taylor-Hood compatible spaces.
In order to have a stable solution {cite}`quarteroni_2016`, the most common couple is given by parabolic FE $P2$ for the velocity and linear finite element for pressure $P1$.

Let us consider the finite dimensional representation of the spaces (using Taylor-Hood elements), the correspondent non-linear system of equations results in
\begin{equation}
\left[
\begin{array}{cc}
A+C(\mathbf{U}) & B^T \\ B & 0
\end{array}
\right]\cdot
\left[
\begin{array}{c}
\mathbf{U} \\ \mathbf{P}
\end{array}
\right] = 
\left[
\begin{array}{c}
\mathbf{F} \\ \mathbf{0}
\end{array}
\right]
\end{equation}

Different techniques can be used to linearised this system and get a solution. The treatment of the non-linearity is an important concern in the solution of the Navier-Stokes equations.

## Treatment of the non-linear term
Let us define the Reynolds number

$$
Re = \frac{\mbox{Inertia forces}}{\mbox{Viscous forces}}=\frac{U_0\cdot D_h}{\nu}
$$
given $U_0$ a characteristic velocity, $D_h$ is the hydraulics diameter and it represents a characteristic length and $\nu$ is the kinematic viscosity. This dimensionless quantity measures the importance of inertia forces with respect of viscous forces, when the latter is stronger the flow is said to be **laminar** whereas as the numerator increases advection becomes dominant and the flow becomes **turbulent**.

A steady solution exists only for values of $Re<3000$, otherwise turbulence makes the flow unsteady.

```{prf:remark}
:label: my-remark

A steady solution can be obtained when introducing the Reynolds-Averaged Navier-Stokes (RANS) equations, in which turbulence is completely modelled {cite}`davidson2015turbulence`.

In this book, only laminar flow is considered.
```