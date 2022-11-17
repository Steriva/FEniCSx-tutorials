# Madelung Transformation

The *Incompressible Schrodinger Flow* is a novel numerical technique, developed in {cite}`Chern2016` and {cite}`Chern2017`, able to describe inviscid fluids behaviour using the Schrodinger equation. This section gives some basic concepts on the analogy between hydrodynamics and quantum mechanics {cite}`Madelung1926`, onto which ISF is based, along with the description of the solution algorithm.

The Schrodinger equation is able to describe the time evolution of this physical state, subjected to a certain potential field $p$
\begin{equation}
    i\hbar \frac{\partial\psi}{\partial t} = -\frac{\hbar^2}{2}\Delta \psi + p \psi.
\end{equation}
This PDE is composed by linear operators for the wave function and it can be classified as parabolic. In the context of ISF, the wave function is used to describe inviscid flows, the probabilistic interpretation of the quantum mechanics is no more used, instead the ISF is set in the analogy between hydrodynamics and quantum mechanics.

This PDE is composed by linear operators for the wave function and it can be classified as parabolic. In the context of ISF, the wave function is used to describe inviscid flows, the probabilistic interpretation of the quantum mechanics is no more used, instead the ISF is set in the analogy between hydrodynamics and quantum mechanics, proposed by Madelung in 1926. The state of the system can be described by a 2 components wave function, i.e. $\Psi=[\psi_1,\psi_2]^T\in\mathbb{C}^2$, which is the solution of the following Schrodinger equation, i.e.
\begin{equation}
    i\hbar\frac{\partial\Psi}{\partial t} =- \frac{\hbar^2}{2}\Delta \Psi+p\Psi\longleftrightarrow
    i\hbar\frac{\partial}{\partial t}\left[\begin{array}{cc} \psi_1 \\ \psi_2\end{array}\right] = -\frac{\hbar^2}{2}\Delta \left[\begin{array}{cc} \psi_1 \\ \psi_2\end{array}\right]+p\left[\begin{array}{cc} \psi_1 \\ \psi_2\end{array}\right].
\end{equation}
```{note}
If the flow field is described by a single wave function wave function $\psi\in\mathbb{C}$, the associated velocity field is always characterised by a null vorticity, which is not generally true
```
It can be shown that the real and the imaginary part of this equation are linked to the continuity and momentum equation of the incompressible Euler equations
\begin{equation}
    \left\{
    \begin{array}{l}
        \nabla\cdot \vec{u} = 0 \\
        \displaystyle\frac{\partial\vec{u}}{\partial t}+\left(\vec{u}\cdot \nabla\right)\vec{u} = -\nabla p
    \end{array}
    \right.
\end{equation}
This translation is known as the Madelung transformation, in particular for $2-$component wave function the velocity can be computed as
\begin{equation}
    \vec{u} = \hbar\mbox{Re}\left\{-\left(\Psi^T\right)^*i\,\nabla\Psi\right\}.
\end{equation}
The viscosity, thus the dissipation, is not considered and the only existing external force is given by the potential $-\nabla p$. The spatial operators on the wave function $\Psi$, i.e. $-\frac{\hbar^2}{2}\Delta$ and $p$, have a physical interpretation in the Euler equations on $\vec{u}$, i.e. the advection $\left(\vec{u}\cdot \nabla\right)\vec{u}$ and the potential forces $-\nabla p$.

In the end, the incompressibility constraint $\nabla \cdot \vec{u}=0$ can be equivalently imposed to the wave function as
\begin{equation}
    \mbox{Re}\left\{\left(\Psi^T\right)^*i\,\Delta\Psi\right\} = 0.
    \label{eqn:SchrodingerIncCoinstraint}
\end{equation}