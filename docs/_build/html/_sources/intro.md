# Introduction

This book collects some basic theory and tutorials for the Finite Element library [FEniCSx](https://fenicsproject.org).

The book is organised in different chapters, according to the different physical problem analysed. Typically, the weak formulation is derived and the calculations are reported.

Moreover, a detailed analysis of the implementation is presented, along with some installation notes for [FEniCSx](https://fenicsproject.org) using *conda* on a computer (Linux or Mac) or the use in Google Colab {cite}`FEMonColab`.


In the book, the index notation is used which is briefly explained here
```{prf:assumption}
The main property to remember is that repeated indeces substitues summations

$$ 
\sum_i u_iv_i = u_i v_i
$$
Vectors as $\mathbf{a}$ are indicated with $a_i$, matrices $A$ with $A_{ij}$. 

Here are some differential operators:

$$ 
\nabla u \longleftrightarrow \frac{\partial u}{\partial x_i}\qquad 
\nabla \cdot \mathbf{u} \longleftrightarrow \frac{\partial u_i}{\partial x_i}\qquad 
\Delta \cdot u \longleftrightarrow \frac{\partial u}{\partial x_i^2}\qquad
\Delta \cdot \mathbf{u} \longleftrightarrow \frac{\partial u_j}{\partial x_i^2}
$$

```

```{tableofcontents}
```
