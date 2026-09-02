
<h1 style="color:teal;">Finite volume solver for Reactive Sedimentation</h1>

This Github repository contains the source files of a solver written in Matlab tailored to approximate the reactive sedimentation up to order $`\mathcal{O}(\Delta z^3)`$ of accuracy in one spatial dimension. 

To acknowledge the use of this software, we kindly ask you to cite the (current) paper:

$\color{blue}\texttt{(Current version)}$
J.D. Barajas, J. Careaga, L.M. Villada. **Invariant-region-preserving high-order schemes for a model of reactive sedimentation**, 
*arXiv preprint arXiv:XXXX.XXXXX*, 2026, https://doi.org/XXXXXXXXXXX


The main structure of the code and a number of subroutines in this program were partially based on the paper
**SDIMA_MOL**


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


## Running the program


```console
to do
```

## Output data



## Authorship

This Matlab-based solver has been developed by **Julio Careaga** (https://github.com/juliocareaga/) and **Juan-David Barajas-Calogne** (https://github.com/juanbarajascalogne/).




