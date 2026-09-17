
<h1 style="color:teal;">Finite volume solver for Reactive Sedimentation</h1>

This Github repository contains the source files of a finite volume solver written in Matlab designed to approximate the reactive sedimentation model (in one spatial dimension) from [Bürger, Careaga & Diehl (2021)](https://academic.oup.com/imamat/article-abstract/86/3/514/6278612?redirectedFrom=fulltext), with first, second or third order of accuracy in space and time, respectively. For the case of the second and third-order, MUSCL and central WENO reconstructions are implemented. This model of reactive sedimentation is an extension to the so-called **Bürger-Diehl model** for secondary settling tanks in simulations of wastewater treatment processes.

This repository is intended to provide an open source solver for the simulation of reactive and non-reactive sedimentation processes, and also to offer a technical higher-order numerical scheme for specialists in applied mathematics and numerical analysis. You are welcome to use this software, to elaborate further simulation tests and to extend and implement the reaction terms to include more general activated sludge models such as the [ASM1](https://iwaponline.com/ebooks/book/96/Activated-Sludge-Models-ASM1-ASM2-ASM2d-and-ASM3). However, we kindly ask you to to acknowledge the use of this software by citing the (current) paper:

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

where by components, we have:
- $`\qquad\color{#076D82}X_{\rm OHO}\text{ is the concentration of heterotrophic organics}`$
- $`\qquad\color{#076D82}X_{\rm U}\text{ is the concentration of undegradable matter}`$
- $`\qquad\color{#822E07}S_{\rm NO3}\text{ is the concentration of nitrate substrate}`$
- $`\qquad\color{#822E07}S_{\rm S}\text{ is the concentration of readily biodegradable substrate}`$
- $`\qquad\color{#822E07}S_{\rm N2}\text{ is the concentration of nitrogen substrate}`$


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


## Input data

The input data is passed through a structure, which may contain the following fields:

| Field   |  Description  |
|---               |---                                                           |
| **par.HH**       | vessel height [m]                                            |
| **par.BB**       | vessel depth  [m]                                            |
| **par.rho**      | density of solids [kg/m^3]                                   |
| **par.rhoL**     | density of liquid [kg/m^3]                                   |
| **par.gg**       | gravity [m/s^2]                                              |
|---               |---                                                           |
| hindered-settling velocity                                                      |
|---               |---                                                           |
| **par.v0**       | hindered settling velocity at zero concentration [m/s]       |
| **par.utilde**   | parameter used in the Diehl flux function [kg/m^3]           |
| **par.eta**      | parameter used in the Diehl flux function [-]                |
| **par.umax**     | maximum solids concentration [kg/m^3]                        |
|---               |---                                                           |
| compression function                                                            |
|---               |---                                                           |
| **par.beta**     | parameter used in the compression function [-]               |
| **par.uc**       | critical concentration                                       |
| **par.diffScheme** | approximation type of the diffusion/compression terms      |
| denitrification kinetics                                                        |
| **par.Y**        | yield constant in the reduced biokinetic model               |
| **par.Ybar**     | additional yield constant in the reduced biokinetic model    |
| **par.fp**       | portion that decays to non-degradable organics               |
| **par.bdec**     | decay rate of heterotrophic organisms [1/s]                  |
| **par.mumax**    | maximum growth rate [1/s]                                    |
| **par.k1**       | saturation constant [kg/m^3]                                 |
| **par.k2**       | saturation constant [kg/m^3]                                 |
| **par.sigmac**   | stoichiometric matrix for the solid components [-]           |
| **par.sigmas**   | stoichiometric matrix for the substrates       [-]           |
|---               |---                                                           |
| bulk flows: batch                                                               |
|---               |---                                                           |
| **par.feedMode** | batch or continuous                                          |
| **par.cfrac**    | fraction of initial feed (taken constant; vector)            |
| **par.sf_val**   | feed concentration of substrates (taken constant; vector)    |
|---               |---                                                           |
| numerical boundary conditions                                                   |
|---               |---                                                           |
| **par.bc**       |  'transmissive'; gamma == 1 at every face and ghost cells    |
|---               |---                                                           |
| initial condition                                                               |
|---               |---                                                           |
| **par.icType**   | type of initial condition: 'jump' or 'bump'                  |
| **par.bumpP**    | power used in smooth initial condition                       |
| **par.z0**       | z-value used in 'bump' [m]                                   |
| **par.Rm**       | radious used in 'bump' [m]                                   |
| **par.ubar**     | u-value used in 'bump' [kg/m^3]                              |
| **par.uhat**     | u-value used in 'bump' [kg/m^3]                              |
| **par.sbar**     | s-value used in 'bump' [kg/m^3]                              |
| **par.shat**     | s-value used in 'bump' [kg/m^3]                              |
|---               |---                                                           |
| numerical parameters                                                            |
|---               |---                                                           |
| **par.Tfinal**   | final simulation time                                        |
| **par.CFL**      | CFL constant                                                 |
| **par.alphaM**   | Muscl parameter                                              |
| **par.wenoPower**| power used in the CWENO method, in the S_K factors           |
| **par.posEps**   | epsilon used in the CWENO method, in the S_K factors         |
| **par.useLim**   | Zhang-Shu limiter ('true' or 'false')                        |
| **par.nsnap**    | number of evenly saved snapshots of the solution             |
| **par.nmon**     | number of evenly saved min/max monitored values              |
| **par.dtRule**   | type of time step dt: 'cfl' or 'speed'                       |
| **par.sameDt**   | same rule (w_1 = 1/6) for all three reconstructions          |
|---               |---                                                           |


## Running the program

If you already have the parameters loaded in, for instance, the structure 'par', and consider the number of cells N, then you run:
```console
cd src
SOL = main_solver(N,par)
```
| Field   |  Description  |
|--- |---    |
| **SOL.N** | Number of cells |
| **SOL.T** | ending time |
| **SOL.dt** | time step |
| **SOL.nSteps** | total number of time iterations |
| **SOL.recon** | reconstruction type: 'fo', 'muscl', 'cweno3' |
| **SOL.rk** | time approximation: 'euler', 'ssprk2', 'ssprk3' |
| **SOL.useLim** | Zhang & Shu limiter: 0 ('false'), 1 ('true')  |
| **SOL.diffScheme** | diffusion approximation: 'none', d-method-2, d-method-4 |
|                    |                                  b-method-2, b-method-4 |
| **SOL.umax** | maximum solids concentration |
| **SOL.z** | mesh nodes in the z-coordinate |
| **SOL.tsnap** | time point snapshots |
| **SOL.usnap** | saved snapshots of vector solution u, size nsnap×N |
| **SOL.csnap** | saved snapshots of vector solution c, size nc×nsnap×N |
| **SOL.ssnap** | saved snapshots of vector solution s, size ns×nsnap×N |
| **SOL.tw** | monitored times for min/max values, size nmon |
| **SOL.minu** | monitored min values of u, size nmon |
| **SOL.maxu** | monitored max values of u, size nmon |
| **SOL.minc** | monitored min values of c, size nmon | 
| **SOL.mins** | monitored min values of s, size nmon |
| **SOL.sumc** | monitored sum of c, size nmon   |
| **SOL.U** | vector of solutions u, last time iteration, size 1×N |
| **SOL.C** | vector of solutions c, last time iteration, size nc×N |
| **SOL.S** | vector of solutions s, last time iteration, size ns×N | 
| **SOL.glob** | structure containing the max and min values of the variables |
|--- |--- |

## Authorship

This Matlab-based solver has been developed by **Julio Careaga** (https://github.com/juliocareaga/) and **Juan Barajas-Calonge** (https://github.com/juanbarajascalonge/).



