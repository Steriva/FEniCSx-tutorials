# ISF algorithm

Let $\Omega$ be the spatial domain, $\partial \Omega$ its boundary, usually composed by three different components $\Gamma_{in}\cup \Gamma_w\cup \Gamma_o$ and let $\mathcal{T} = [0, T]$ be the time interval considered. Given an initial condition $\Psi^0$, the {Schrodinger equation} with the incompressibility constraint, can be solved with a first order time splitting method

```{note}
A similar solution strategy is also employed in CFD, a very important example is the pioneering work of {cite}`Chorin1968`.
```

At each time step $t_n\rightarrow t_{n+1}$, the algorithm, in the space-continuous case solve the following problems:

- Schrodinger problem
- Normalisation
- Poisson equation (*pressure projection*)
- Phase Shift and velocity update

## Schrodinger problem
The prediction step consists of a free-particle Schrodinger equation
\begin{equation}
    \frac{\tilde{\Psi}-\Psi^n}{\Delta t} = \frac{i \hbar}{2}\Delta \tilde{\Psi}
\end{equation}
with the associated boundary conditions (a plane wave at the inlet and homogeneous Neumann boundaries for the others), i.e. 
\begin{equation}
     \left\{
        \begin{array}{ll}
            \displaystyle\tilde{\Psi} = e^{i\left(\vec{k}\cdot \vec{x}-\omega t\right)}\cdot \left[c_1, c_2\right]^T & \text{ on }\Gamma_{in}\\
            \nabla\tilde{\Psi}\cdot \vec{n}=0 & \text{ on }\Gamma_w\cup \Gamma_o
        \end{array}
        \right.
    \end{equation}
with $c_1,c_2 \in\mathbb{C}$ such that $|c_1|^2+|c_2|^2=1$.
## Normalisation
The wave function needs also to be normalised to 1, since this property is not necessarily conserved in the time-discrete case, i.e.
\begin{equation}
    \tilde{\Psi} \longleftarrow\frac{\tilde{\Psi}}{\;||\tilde{\Psi}||_2} =
    \frac{\tilde{\Psi}}{\sqrt{\tilde{\psi}_1^*\tilde{\psi}_1+\tilde{\psi}_2^*\tilde{\psi}_2}}  .
\end{equation}
## Pressure Projection
The correction phase involves a Poisson problem to enforce the incompressibility constraint. The unknown is $\varphi$, which is referred to as \textit{pressure}, even though its units of measure are $[m^2/s]$. This quantity has been introduced following the idea developed in {cite}`Chern2017`, which provides a complete discussion on this choice. Thus, the following Poisson problem is solved
\begin{equation}
        \left\{
    \begin{array}{ll}
       \Delta \varphi = \nabla \cdot \tilde{\vec{u}} & \mbox{ in }\Omega \\
       \varphi = 0 & \mbox{ on } \Gamma_{in}\\
        \nabla\varphi\cdot \vec{n}=0 & \mbox{ on } \Gamma_w\\ 
        \nabla\varphi\cdot \vec{n}=g & \mbox{ on } \Gamma_o
    \end{array}
    \right. 
    \label{eqn: Poisson-Strong}
\end{equation}
in which $\tilde{\vec{u}}$ is the velocity field computed from $\tilde{\Psi}$.

This phase is called *pressure projection*
```{note}
This idea is very similar to the standard CFD splitting method, in which a Poisson problem is solved to enforce the incompressibility.
```

## Phase Shift and Velocity Update
