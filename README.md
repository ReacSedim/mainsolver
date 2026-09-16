
<h1 style="color:teal;">Finite volume solver for Reactive Sedimentation</h1>

This Github repository contains the source files of a finite volume solver written in Matlab designed to approximate the reactive sedimentation model (in one spatial dimension) from [Bürger, Careaga & Diehl (2021)](https://academic.oup.com/imamat/article-abstract/86/3/514/6278612?redirectedFrom=fulltext), with first, second or third order of accuracy in space and time, respectively. For the case of the second and third-order, MUSCL and central WENO reconstructions are implemented. This model of reactive sedimentation is an extension to the so-called **Bürger-Diehl model** for secondary settling tanks in simulations of wastewater treatment processes.

This repository is intended to provide an open source solver for the simulation of reactive and non-reactive sedimentation processes, and also to offer a technical higher-order numerical scheme for specialists in applied mathematics. You are welcome to use this software, to elaborate further simulation tests and to extend and implement the reaction terms to include more general activated sludge models such as the [ASM1](https://iwaponline.com/ebooks/book/96/Activated-Sludge-Models-ASM1-ASM2-ASM2d-and-ASM3). However, we kindly ask you to to acknowledge the use of this software by citing the (current) paper:

$\color{blue}\texttt{(Current version)}$
J. Barajas-Calonge, J. Careaga, L.M. Villada. **Invariant-region-preserving high-order schemes for a model of reactive sedimentation**, 
*ArXiv preprint arXiv:2609.06846*, 2026, [https://arxiv.org/abs/2609.06846](https://arxiv.org/abs/2609.06846)


--------------------------

### 1. Model equations

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
\color{#076D82}\frac{\partial \boldsymbol{c}}{\partial t}
&\color{#076D82}+ \frac{\partial}{\partial z}
\left[
\left(q(z,t)+\gamma(z)\left(v_{\rm hs}(u)-\partial_z\mathcal{B}(u)\right)\right)
\boldsymbol{c}
\right]
= \boldsymbol{b}_{\boldsymbol{c}}(\boldsymbol{c},\boldsymbol{s},z,t),
\\
\color{#822E07}\frac{\partial \boldsymbol{s}}{\partial t}
&\color{#822E07}+ \frac{\partial}{\partial z}
\left[
\left(q(z,t)-\gamma(z)
\frac{\left(v_{\rm hs}(u)-\partial_z\mathcal{B}(u)\right)u}{\rho-u}\right)
\boldsymbol{s}
\right]
= \boldsymbol{b}_{\boldsymbol{s}}(\boldsymbol{c},\boldsymbol{s},z,t).
\end{aligned}
$$

where $`\boldsymbol{c}`$ is the vector of $`{n_{\boldsymbol{c}}}`$ solid components, $`\boldsymbol{s}`$ is the vector $`{n_{\boldsymbol{s}}}`$ substrates, and $`u`$ is the total concentration of solids. In the default example implemented in this program, the [activated sludge reduced biokinetic model](https://www.sciencedirect.com/science/article/abs/pii/S0098135416301338), the vectors are:

$$
\begin{aligned}
\boldsymbol{c} &= \left(X_{\rm OHO},X_{\rm U}\right)\qquad\text{and}\qquad
\boldsymbol{s} = \left(S_{\rm NO3},S_{\rm S},S_{\rm N2}\right),
\end{aligned}
$$

where by components, we have: <br><br>
$`\qquad\color{#076D82}X_{\rm OHO}\text{ is the concentration of heterotrophic organics}`$,<br>
$`\qquad\color{#076D82}X_{\rm U}\text{ is the concentration of undegradable matter}`$,<br>
$`\qquad\color{#822E07}S_{\rm NO3}\text{ is the concentration of nitrate substrate}`$,<br>
$`\qquad\color{#822E07}S_{\rm S}\text{ is the concentration of readily biodegradable substrate}`$,<br>
$`\qquad\color{#822E07}S_{\rm N2}\text{ is the concentration of nitrogen substrate}`$.


The numerical scheme employed to solve the PDE combines a variety of ingredients:
- High-order time approximations for ....
- High-order polynomial reconstructions and maximum-principle limiters.

## Code structure

### Repositories and organization

The software is organized in three subfolders:
- **src**: Contains all the ".m" source codes
- **post-processing**: This is a repository where compressed auxiliary files ...
- **results**: This folder is intended to store all the results obtained from running the simulations; therefore, it is initially empty. 

### Scripts

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




