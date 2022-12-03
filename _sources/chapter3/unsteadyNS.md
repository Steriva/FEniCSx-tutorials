# Unsteady Navier-Stokes

The Navier-Stokes equations describes the flows of incompressible viscous fluid 
\begin{equation}
\left\{
    \begin{array}{ll}
        \displaystyle\frac{\partial \mathbf{u}}{\partial t} +\left(\mathbf{u}\cdot \nabla\right)\mathbf{u}-\nu\Delta \mathbf{u}+\nabla p = 0 & \mbox{in }\Omega, \;t>0\\
        \nabla\cdot \mathbf{u} = 0 & \mbox{in }\Omega, \;t>0\\
        \mathbf{u}(\mathbf{x}, 0)=\mathbf{u}_0 & \mbox{in }\Omega\\
        \mathbf{u}=\mathbf{u}_D & \mbox{on }\Gamma_D, \;t>0\\
        \nu\frac{\partial \mathbf{u}}{\partial \mathbf{n}}-p\mathbf{n}=g\mathbf{n} & \mbox{on }\Gamma_N, \;t>0
    \end{array}
\right.
\end{equation}

Similar consideration made for the stationary case can be made even for the unsteady one. We are dealing with a non-linear parabolic problem, here we will present of the most famous algorithm to solve the unsteady Navier-Stokes equations, i.e. the incremental Chorin-Theman algorithm {cite}`Chorin1968`.

Let $\mathcal{V}\subset[\mathcal{H}^1]^d,\; \mathcal{V}_0[\subset\mathcal{H}^1]^d$ the velocity trial and test spaces defined as

$$
\mathcal{V} = \left\{\mathbf{v}\in[\mathcal{H}^1(\Omega)]^d:\;\left. \mathbf{v}\right|_{\Gamma_D} = \mathbf{u}_D\right\}\qquad 
\mathcal{V}_0 = \left\{\mathbf{v}\in[\mathcal{H}^1(\Omega)]^d:\;\left. \mathbf{v}\right|_{\Gamma_D} = \mathbf{0}\right\}
$$
and let $\mathcal{Q}= L^2(\Omega)$ the pressure trial and test space. The continuity and momentum equation can be multiplied by the test functions $\mathbf{v}\in\mathcal{V}_0,\;q\in\mathcal{Q}$, so the weak formulation reads: *find $(\mathbf{u}(t), p(t))\in\mathcal{V}\times\mathcal{Q}$ s.t.*
\begin{equation}
\int_\Omega \frac{\partial \mathbf{u}}{\partial t} \cdot \mathbf{v}\,d\Omega +\int_\Omega \left[\left(\mathbf{u}\cdot \nabla\right)\mathbf{u}\right]\cdot \mathbf{v}\,d\Omega+\int_\Omega \nu\nabla \mathbf{u}\cdot \nabla\mathbf{v}\,d\Omega -\int_\Omega p \nabla\cdot \mathbf{v}\,d\Omega - \int_\Omega q \nabla\cdot \mathbf{u}\,d\Omega = \int_{\Gamma_N} g\mathbf{n}\cdot \mathbf{v}\,d\sigma \qquad \forall (\mathbf{v}, q)\in\mathcal{V}_0\times\mathcal{Q}
\end{equation}
## Chorin-Theman algorithm
From the weak formulation a discrete system of non-linear ODEs can be derived and theoretically solved, the dimension of the system usually to large preventing a straightforward numerical solution. Therefore, the state-of-the-art for solving the Navier-Stokes system consists of fractional step (or projection) methods.

The pioneering work of {cite}`Chorin1968` was the inspiration for all successive developments. Here we will present a strategy quite similar to the original one, useful when laminar flows are simulated. The discrete time algorithm is made of 2 main steps, followed by the update one.

### Prediction
The continuity and the momentum balance from the NS equations in strong form are separated. At first, the momentum solved at time $t^{n+1}$

\begin{equation}
\left\{
    \begin{array}{ll}
        \displaystyle\frac{\tilde{\mathbf{u}}-\mathbf{u}^n}{\Delta t} +\left(\mathbf{u}^n\cdot \nabla\right)\tilde{\mathbf{u}}-\nu\Delta\tilde{\mathbf{u}}+\nabla p^n = 0 & \mbox{in }\Omega \\
        \mathbf{u}=\mathbf{u}_D & \mbox{on }\Gamma_D\\
        \nu\frac{\partial \mathbf{u}}{\partial \mathbf{n}}=g\mathbf{n} & \mbox{on }\Gamma_N
     \end{array}
\right.
\end{equation}
given ${\mathbf{u}}^n$ and $p^n$ to be the solution at time $t^n$.

The weak formulation of this problem can be derived in a similar way used for the steady NS equations.

### Pressure projection
The velocity field $\tilde{\mathbf{u}}$ obtained at the prediction step is not divergence-free. The core of fractional step method is the splitting of the differential operator: the remaining terms of the NS momentum are

$$
 \frac{\mathbf{u}^{n+1}-\tilde{\mathbf{u}}}{\Delta t} +\nabla (p^{n+1}-p^n) = 0
 $$
If these two equations are summed, the NS momentum equation discretised using an implicit Euler scheme (with a semi-implicit treatment of the non-linear term) is retrieved. 

Let us impose $\nabla \cdot \mathbf{u}^{n+1}=0$, thus the following Poisson problem is obtained

\begin{equation}
\left\{
    \begin{array}{ll}
        \nabla\cdot\mathbf{u}^n =\Delta t\,\Delta \delta p & \mbox{in }\Omega \\
        \frac{\partial p}{\partial \mathbf{n}}=0 & \mbox{on }\Gamma_D\\
        p=0 & \mbox{on }\Gamma_N
     \end{array}
\right.
\end{equation}

The weak formulation of this problem can be derived very easily.
### Correction
Once the pressure field is obtained and used to make the velocity field divergence-free, the update of the both fields is

$$
\mathbf{u}^{n+1} = \tilde{\mathbf{u}} -\Delta t\nabla \delta p = 0\qquad p^{n+1}=p^n+\delta p
$$

## Direct numerical simulation of turbulent flows
As the importance of the advection term $(\mathbf{u}\cdot\nabla)\mathbf{u}$ (i.e. inertia) increases over the viscous term $\nu\Delta\mathbf{u}$ (i.e. dissipation), the flow has a transition from laminar to turbulent.


$$
\frac{L_D}{L_0} \sim Re^{-3/4}\qquad 
\mathbf{U_D}{U_0}\sim Re^{-1/4}\qquad
\mathbf{t_D}{t_0}\sim Re^{-1/2}\qquad
Re = \frac{U_0L_0}{\nu}
$$