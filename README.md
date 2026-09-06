
<h1 style="color:teal;">Finite volume solver for Reactive Sedimentation</h1>

This Github repository contains the source files of a solver written in Matlab designed to approximate the reactive sedimentation model (in one spatial dimension) from [Bürger, Careaga & Diehl (2021)](https://academic.oup.com/imamat/article-abstract/86/3/514/6278612?redirectedFrom=fulltext), with first, second or third order of accuracy. 

You are welcome to use this software, to elaborate further simulation tests and to extend and implement the reaction terms to include more general activated sludge models such as the ASM1. However, we kindly ask you to to acknowledge the use of this software by citing the (current) paper:

$\color{blue}\texttt{(Current version)}$
J. Barajas-Calonge, J. Careaga, L.M. Villada. **Invariant-region-preserving high-order schemes for a model of reactive sedimentation**, 
*arXiv preprint arXiv:To set*, 2026, https://doi.org/To-set

--------------------------

The benchmark partial differential equation (PDE), in the one-dimensional case,  solved by this software is the following first-order system of $`n`$ equations:


The numerical scheme employed to solve the PDE combines a variety of ingredients:
- High-order time approximations for ....
- High-order polynomial reconstructions and maximum-principle limiters.


### Repositories and organization

The software is organized in three subfolders:
- **src**: Contains all the ".m" source codes, ...
- **post-processing**: This is a repository where compressed auxiliary files ...
- **examples**: 


## Code structure

Scripts in **scr** subfolder:
- **rs_const_functions.m**: contains all constituve functions of the model
- **rs_defaults.m**: fill in optional fields and compute derived constants
- **rs_diffusive_flux.m**: compute the diffusive flux approximation operator...
- **rs_ghost.m**: create ghost cells for computing boundary conditions
- **rs_initial.m**: compute cell averages of the initial conditions
- **rs_reconstruct.m**: compute second-order MUSCL and third-order CWENO reconstructions

## Running the program


```console
to do
```

## Output data



## Authorship

This Matlab-based solver has been developed by **Julio Careaga** (https://github.com/juliocareaga/) and **Juan Barajas-Calonge** (https://github.com/juanbarajascalonge/).




