# Unsteady Navier-Stokes

The Navier-Stokes equations describes the flows of incompressible viscous fluid 
\begin{equation}
\left\{
    \begin{array}{ll}
        \frac{\partial \mathbf{u}}{\partial t} +\left(\mathbf{u}\cdot \nabla\right)\mathbf{u}-\nu\Delta \mathbf{u}+\nabla p = 0 & \mbox{in }\Omega, \;t>0\\
        \nabla\cdot \mathbf{u} = 0 & \mbox{in }\Omega, \;t>0\\
        \mathbf{u}(\cdot, t)=\mathbf{u}_0 & \mbox{in }\Omega,, \;t=0\\
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


## Direct numerical simulation of turbulent flows

$$
\frac{L_D}{L_0} \sim Re^{-3/4}\qquad 
\mathbf{U_D}{U_0}\sim Re^{-1/4}\qquad
\mathbf{t_D}{t_0}\sim Re^{-1/2}\qquad
Re = \frac{U_0L_0}{\nu}
$$