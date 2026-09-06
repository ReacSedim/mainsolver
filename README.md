
<h1 style="color:teal;">Finite volume solver for Reactive Sedimentation</h1>

This Github repository contains the source files of a finite volume solver written in Matlab designed to approximate the reactive sedimentation model (in one spatial dimension) from [Bürger, Careaga & Diehl (2021)](https://academic.oup.com/imamat/article-abstract/86/3/514/6278612?redirectedFrom=fulltext), with first, second or third order of accuracy in space and time. For the case of the second and third-order, MUSCL and central WENO reconstructions are implemented. 

You are welcome to use this software, to elaborate further simulation tests and to extend and implement the reaction terms to include more general activated sludge models such as the ASM1. However, we kindly ask you to to acknowledge the use of this software by citing the (current) paper:

$\color{blue}\texttt{(Current version)}$
J. Barajas-Calonge, J. Careaga, L.M. Villada. **Invariant-region-preserving high-order schemes for a model of reactive sedimentation**, 
*arXiv preprint arXiv:To set*, 2026, https://doi.org/To-set



--------------------------

### Model equations

The benchmark partial differential equation (PDE), in the one-dimensional case,  solved by this software is the following first-order system of $`n_{\boldsymbol{c}}+n_{\boldsymbol{s}}+1`$ equations:


$$
\begin{aligned}
\frac{\partial u}{\partial t}
&+ \frac{\partial}{\partial z}
\left[
\left(q(z,t)+\gamma(z)\left(v_{\rm hs}(u)-\partial_z\mathcal{B}(u)\right)\right)u
\right]
= b_u(\boldsymbol{c},\boldsymbol{s},z,t),
\\
\frac{\partial \boldsymbol{c}}{\partial t}
&+ \frac{\partial}{\partial z}
\left[
\left(q(z,t)+\gamma(z)\left(v_{\rm hs}(u)-\partial_z\mathcal{B}(u)\right)\right)
\boldsymbol{c}
\right]
= \boldsymbol{b}_{\boldsymbol{c}}(\boldsymbol{c},\boldsymbol{s},z,t),
\\
\frac{\partial \boldsymbol{s}}{\partial t}
&+ \frac{\partial}{\partial z}
\left[
\left(q(z,t)-\gamma(z)
\frac{\left(v_{\rm hs}(u)-\partial_z\mathcal{B}(u)\right)u}{\rho-u}\right)
\boldsymbol{s}
\right]
= \boldsymbol{b}_{\boldsymbol{s}}(\boldsymbol{c},\boldsymbol{s},z,t).
\end{aligned}
$$

where $`\boldsymbol{c} = (c_1,c_2,...,c_{n_{\boldsymbol{c}}})`$ is the vector of solid components, $`\boldsymbol{s} = (s_1,s_2,...,s_{n_{\boldsymbol{s}}})`$ is the vector substrates, and $`u`$ is the total concentration of solids. In the default example implemented in this program, the activated sludge reduced biokinetic model, the vectors are:

$$
\begin{aligned}
\boldsymbol{c} &= \left(X_{\rm OHO},X_{\rm U}\right)
\\
\boldsymbol{s} &= \left(S_{\rm N2},S_{\rm S},S_{\rm NO3}\right)
\end{aligned}
$$

where <br>
$`\color{purple}X_{\rm OHO}`$ is the concentration of heterotrophic organics<br>
$`\color{teal!80!blue}X_{\rm U}`$ is the concentration of undegradable matter<br>
$`S_{\rm NO3}`$ is the concentration of nitrate substrate<br>
$`S_{\rm S}`$ is the concentration of readily biodegradable substrate<br>
$`S_{\rm N2}`$ is the concentration of nitrogen substrate


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




