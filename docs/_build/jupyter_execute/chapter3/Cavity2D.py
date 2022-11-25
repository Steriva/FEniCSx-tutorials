#!/usr/bin/env python
# coding: utf-8

# # Cavity (2D)
# This notebook implements a steady incompressible Navier-Stokes solver for the cavity problem.
# 
# The problem is strong form reads:
# \begin{equation}
# \left\{
# \begin{array}{ll}
#     \nabla \cdot \mathbf{u} =0& in\;\Omega\\
#     \displaystyle \left(\mathbf{u}\cdot \nabla\right)\mathbf{u}= \nu\Delta \mathbf{u}-\nabla p & in\;\Omega\\ & \\
#     \mathbf{u} = \mathbf{i} & on\;\Gamma_{lid}\\
#     \mathbf{u} = \mathbf{0} & on\;\partial\Omega\setminus\Gamma_{lid}
# \end{array}
# \right.
# \end{equation}

# In[6]:


import tqdm
import numpy as np

# Mesh generation
import dolfinx
from mpi4py import MPI
from dolfinx import mesh
from dolfinx.io import gmshio, XDMFFile
from dolfinx import fem
from dolfinx.fem import (Function, FunctionSpace, dirichletbc, locate_dofs_topological, 
                         form, assemble_scalar, locate_dofs_geometrical)
import ufl
from ufl import grad, div, nabla_grad, dx, inner, dot
from petsc4py import PETSc
from dolfinx.mesh import (CellType, GhostMode, create_rectangle, locate_entities_boundary)

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import cm

plt.rcParams.update({
  "text.usetex": True,
  "font.family": "serif"
})

rcParams['text.latex.preamble'] = r'\usepackage{amssymb} \usepackage{amsmath} \usepackage{amsthm} \usepackage{mathtools}'


# The problem we want to face is non-linear, whose weak formulation reads:
# \begin{equation}
# \int_\Omega \left(\mathbf{u}\cdot \nabla\right)\mathbf{u}\cdot \mathbf{v}\,d\Omega + \nu \int_\Omega\nabla \mathbf{u}\cdot \nabla \mathbf{v}\,d\Omega -\int_\Omega p(\nabla\cdot\mathbf{v})\,d\Omega -\int_\Omega q(\nabla\cdot\mathbf{u})\,d\Omega=0
# \end{equation}

# In[7]:


# Function to mark x = 0, x = 1 and y = 0
def noslip_boundary(x):
    return np.logical_or(np.logical_or(np.isclose(x[0], 0.0),
                                       np.isclose(x[0], 1.0)),
                         np.isclose(x[1], 0.0))

# Function to mark the lid (y = 1)
def lid(x):
    return np.isclose(x[1], 1.0)

# Lid velocity
def lid_velocity_expression(x):
    return np.stack((np.ones(x.shape[1]), np.zeros(x.shape[1])))

class cavity_NS():
    def __init__(self, N):

        # Domain
        self.domain = create_rectangle(MPI.COMM_WORLD,[np.array([0, 0]), np.array([1, 1])], [N, N], CellType.triangle, GhostMode.none)
        self.gdim = self.domain.geometry.dim
        self.fdim = self.gdim - 1

        # Functional Spaces
        self.P2 = ufl.VectorElement("Lagrange", self.domain.ufl_cell(), 2)
        self.P1 = ufl.FiniteElement("Lagrange", self.domain.ufl_cell(), 1)
        self.mixEl = ufl.MixedElement(self.P2, self.P1)
        self.W = FunctionSpace(self.domain, self.mixEl)

        # Test and trial functions: monolithic
        (self.v, self.q) = ufl.TestFunctions(self.W)
        self.dup = ufl.TrialFunction(self.W)
        self.up  = fem.Function(self.W)
        (self.u, self.p) = ufl.split(self.up)

    def parameters(self, Re):

        # Computing the viscosity
        self.nuValue = 1. / Re
        print('Re = {:.0f}'.format(Re)+' and nu = {:.2e}'.format(self.nuValue)+' [m2/s]')
        self.nu = fem.Constant(self.domain, PETSc.ScalarType(self.nuValue))

        # Setting Boundary Conditions
        self.lid_velocity = Function(self.W.sub(0).collapse()[0])
        self.lid_velocity.interpolate(lid_velocity_expression)
        self.ft_lid = locate_entities_boundary(self.domain, self.fdim, lid)
        self.dofs_lid = locate_dofs_topological((self.W.sub(0), self.W.sub(0).collapse()[0]), self.fdim, self.ft_lid)
        self.bc_lid = dirichletbc(self.lid_velocity, self.dofs_lid, self.W.sub(0))

        self.no_slip  = Function(self.W.sub(0).collapse()[0])
        self.ft_walls = locate_entities_boundary(self.domain, self.fdim, noslip_boundary)
        self.dofs_walls = locate_dofs_topological((self.W.sub(0), self.W.sub(0).collapse()[0]), self.fdim, self.ft_walls)
        self.bc_w = dirichletbc(self.no_slip, self.dofs_walls, self.W.sub(0))

        self.zero_p = Function(self.W.sub(1).collapse()[0])
        self.zero_p.x.set(0.0)
        self.dofs_p = locate_dofs_geometrical((self.W.sub(1), self.W.sub(1).collapse()[0]), 
                                              lambda x: np.isclose(x.T, [0, 0, 0]).all(axis=1))
        self.bc_p = dirichletbc(self.zero_p, self.dofs_p, self.W.sub(1))

        self.bcs = [self.bc_lid, self.bc_w, self.bc_p]

    def create_snes_solution(self) -> PETSc.Vec:  # type: ignore[no-any-unimported]
        """
        Create a petsc4py.PETSc.Vec to be passed to petsc4py.PETSc.SNES.solve.

        The returned vector will be initialized with the initial guess provided in `self._solution`.
        """
        x = self._solution.vector.copy()
        with x.localForm() as _x, self._solution.vector.localForm() as _solution:
            _x[:] = _solution
        return x

    def update_solution(self, x: PETSc.Vec) -> None:  # type: ignore[no-any-unimported]
        """Update `self._solution` with data in `x`."""
        x.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
        with x.localForm() as _x, self._solution.vector.localForm() as _solution:
            _solution[:] = _x
            
    def obj_fun(self, snes: PETSc.SNES, x: PETSc.Vec) -> np.float64:
            """Compute the norm of the residual."""
            self.F_fun(snes, x, self._obj_vec)
            return self.b.norm()  # type: ignore[no-any-return]

    def F_fun(self, snes: PETSc.SNES, x: PETSc.Vec, F_vec: PETSc.Vec) -> None:
            """Assemble the residual."""
            self.update_solution(x)
            with F_vec.localForm() as F_vec_local:
                F_vec_local.set(0.0)
            fem.petsc.assemble_vector(F_vec, self._F)
            dolfinx.fem.apply_lifting(F_vec, [self._J], [self.bcs], x0=[x], scale=-1.0)
            F_vec.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
            fem.set_bc(F_vec, self.bcs, x, -1.0)

    def J_fun( self, snes: PETSc.SNES, x: PETSc.Vec, J_mat: PETSc.Mat, P_mat: PETSc.Mat) -> None:
            """Assemble the jacobian."""
            J_mat.zeroEntries()
            fem.petsc.assemble_matrix(J_mat, self._J, self.bcs, diagonal=1.0)
            J_mat.assemble()

    def assemble(self, maxIter):

        # Variational forms
        self.F = (self.nu * inner(grad(self.u), grad(self.v)) * ufl.dx
                + inner(grad(self.u) * self.u, self.v) * ufl.dx
                - inner(self.p, ufl.div(self.v)) * ufl.dx
                + inner(div(self.u), self.q) * ufl.dx)
        self.J = ufl.derivative(self.F, self.up, self.dup)

        self._F = form(self.F)
        self._J = form(self.J)

        # Create matrix and vector
        self._solution = self.up
        self._obj_vec = fem.petsc.create_vector(self._F)
        self.b = fem.petsc.create_vector(self._F)
        self.A = fem.petsc.create_matrix(self._J)

        # Solver settings
        self.solver = PETSc.SNES().create(self.domain.comm)
        self.solver.setTolerances(max_it=maxIter)
        self.solver.getKSP().setType("preonly")
        self.solver.getKSP().getPC().setType("lu")
        self.solver.getKSP().getPC().setFactorSolverType("mumps")

        self.solver.setObjective(self.obj_fun)
        self.solver.setFunction(self.F_fun, self.b)
        self.solver.setJacobian(self.J_fun, J=self.A, P=None)
        # self.solver.setMonitor(lambda _, it, residual: print(it, residual))

    def solve(self):
        up_copy = self.create_snes_solution()
        self.solver.solve(None, up_copy)
        self.update_solution(up_copy)
        return self._solution


# The following function is used to extract the values of the velocity components.

# In[8]:


def extract2D(problem, N, u_sol):
    grid = np.linspace(0, 1, N)
    ux = np.zeros((N,N))
    uy = np.zeros((N,N))

    for ii in range(N):
        points = np.zeros((3, N))
        points[0, :] = grid[ii]
        points[1, :] = grid

        bb_tree = dolfinx.geometry.BoundingBoxTree(problem.domain, problem.domain.topology.dim)
        cells = []
        points_on_proc = []
        cell_candidates = dolfinx.geometry.compute_collisions(bb_tree, points.T)
        colliding_cells = dolfinx.geometry.compute_colliding_cells(problem.domain, cell_candidates, points.T)
        for i, point in enumerate(points.T):
            if len(colliding_cells.links(i))>0:
                points_on_proc.append(point)
                cells.append(colliding_cells.links(i)[0])
        xPlot = np.array(points_on_proc, dtype=np.float64)

        ux[ii, :] = u_sol.sub(0).eval(xPlot, cells).flatten()
        uy[ii, :] = u_sol.sub(1).eval(xPlot, cells).flatten()
    return ux.T, uy.T, grid


# Let's solve the problem

# In[9]:


N = 64
cavity = cavity_NS(N)

ReVec = np.array([20, 100, 400, 1000, 2000, 2500])

u_list = []

for ii in range(len(ReVec)):
    cavity.parameters(ReVec[ii])
    cavity.assemble(20)
    up_sol = cavity.solve()
    (u_sol, p_sol) = (up_sol.sub(0).collapse(), up_sol.sub(1).collapse())
    u_list.append(u_sol)


# Let's make some 2D streamline plots

# In[13]:


fig = plt.figure(figsize = (12,8))

for ii in range(len(ReVec)):
    ux, uy, grid = extract2D(cavity, N, u_list[ii])
    X, Y = np.meshgrid(grid, grid)
    color = np.sqrt(((uy+2)/2)*2 + ((ux+2)/2)*2)
    plt.subplot(2,3,ii+1)
    plt.streamplot(X, Y, ux, uy, color=np.sqrt(ux**2+uy**2), linewidth=1, cmap='turbo', density=2)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xlabel(r'$x$', fontsize = 15)
    plt.ylabel(r'$y$', fontsize = 15)
    plt.title(r'$Re = {:.0f}'.format(ReVec[ii])+'$', fontsize = 15)
    plt.tight_layout()

