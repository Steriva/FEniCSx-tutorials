# Steady Stokes

The Stokes equations describes the flows of incompressible viscous fluid when the importance of inertia is low with respect to viscous forces
\begin{equation}
\left\{
    \begin{array}{ll}
        \nu\Delta \mathbf{u}-\nabla p = 0 & \mbox{in }\Omega\\
        \nabla\cdot \mathbf{u} = 0 & \mbox{in }\Omega\\
        \mathbf{u}=\mathbf{u}_D & \mbox{on }\Gamma_D\\
        \nu\frac{\partial \mathbf{u}}{\partial \mathbf{n}}-p\mathbf{n}=g\mathbf{n} & \mbox{on }\Gamma_N
    \end{array}
\right.
\end{equation}