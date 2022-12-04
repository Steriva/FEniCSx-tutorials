# Steady Navier-Stokes + Energy equation (Boussinesq)

In the previous tutorials, the dependence on the temperature was not considered assuming that the temperature itself can be considered as a passive scalar advected by the flow. However, if the temperature variation is sufficiently high the thermophysical properties changes. This makes necessary the introduction of the energy balance to the Navier-Stokes equations, furthermore the incompressibility assumption may not be valid anymore.

If the temperature variations is sufficiently low, the flow can still be considered incompressible introducing the Boussinesq approximation which results in a linear dependence of the density on the temperature

$$
\rho \simeq \rho_0\,\left[1+\beta(T-T_{ref})\right]
$$

In the end, the steady Navier-Stokes equations and the energy equation under the Boussinesq approximation result in
\begin{equation}
\left\{
\begin{array}{ll}
    \nabla \cdot \mathbf{u} =0& in\;\Omega\\
    \displaystyle \left(\mathbf{u}\cdot \nabla\right)\mathbf{u}= \nu\Delta \mathbf{u}-\nabla p -\beta(T-T_{ref})\mathbf{g} & in\;\Omega\\ 
    \mathbf{u}\cdot \nabla T = \alpha \Delta T& in\;\Omega\\& \\
    \mathbf{u} = \mathbf{u}_D,\; T=T_D & on\;\Gamma_{D}\\
    \mathbf{u} = \mathbf{0},\;-k\frac{\partial T}{\partial \mathbf{n}}=q & on\;\Gamma_{w}\\
    \nu\frac{\partial \mathbf{u}}{\partial \mathbf{n}}-p\mathbf{n}=g\mathbf{n},\;
    \frac{\partial T}{\partial \mathbf{n}}=0  & \mbox{on }\Gamma_N
\end{array}
\right.
\end{equation}
This approximation is used to study steady buoyancy-driven flows.

## Derivition of the weak formulation

The problem we want to face is non-linear, whose weak formulation reads:
\begin{equation}
    \begin{split}
        \int_\Omega \left(\mathbf{u}\cdot \nabla\right)\mathbf{u}\cdot \mathbf{v}\,d\Omega + \nu \int_\Omega\nabla \mathbf{u}\cdot \nabla \mathbf{v}\,d\Omega -\int_\Omega p(\nabla\cdot\mathbf{v})\,d\Omega- \beta\int_\Omega(T-T_{ref})\mathbf{g}\cdot \mathbf{v}\,d\Omega &=0\\
        -\int_\Omega q(\nabla\cdot\mathbf{u})\,d\Omega &=0\\
        \int_\Omega \mathbf{u}\cdot \nabla T \, \phi\,d\Omega+\int_\Omega \alpha\nabla T \cdot \nabla\phi\,d\Omega&=0
    \end{split}
\end{equation}

Thus a proper non-linear solver should be implemented