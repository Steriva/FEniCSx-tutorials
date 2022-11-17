# Installation notes

The full instructions to install FEniCSx are reported on the [Github page of Dolfinx](https://github.com/FEniCS/dolfinx#installation). The easiest way is through *conda*.

Assuming to have correctly conda installed on your machine, the command to install the libraries is the following:
- For real mode:
```console
conda create -n envName -c conda-forge python=3.10 fenics-dolfinx petsc=*=real* mpich pyvista matplotlib tqdm meshio
````
- For complex mode:
```console
conda create -n envName -c conda-forge python=3.10 fenics-dolfinx petsc=*=complex* mpich pyvista matplotlib tqdm meshio
```

The former allows only real PDEs to be solved (e.g., heat equation), whereas the latter includes the generalised version of inner product for complex functions (e.g., Schrodinger equation).

Another possibility consists in using the Google Colab version provided by {cite}`FEMonColab`.