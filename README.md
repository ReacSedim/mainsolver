
<h1 style="color:teal;">Finite volume solver for Reactive Sedimentation</h1>

This Github repository contains the source files of a finite volume solver written in Matlab designed to approximate the reactive sedimentation model (in one spatial dimension) from [Bürger, Careaga & Diehl (2021)](https://academic.oup.com/imamat/article-abstract/86/3/514/6278612?redirectedFrom=fulltext), with first, second or third order of accuracy in space and time, respectively. For the case of the second and third-order, MUSCL and central WENO reconstructions are implemented. This model of reactive sedimentation is an extension to the so-called **Bürger-Diehl model** for secondary settling tanks in simulations of wastewater treatment processes.

This repository is intended to provide an open source solver for the simulation of reactive and non-reactive sedimentation processes, and also to offer a technical higher-order numerical scheme for specialists in applied mathematics and numerical analysis. You are welcome to use this software, to elaborate further simulation tests and to extend and implement the reaction terms to include more general activated sludge models such as the [ASM1](https://iwaponline.com/ebooks/book/96/Activated-Sludge-Models-ASM1-ASM2-ASM2d-and-ASM3). However, we kindly ask you to to acknowledge the use of this software by citing the (current) paper:

$\color{blue}\texttt{(Current version)}$
J. Barajas-Calonge, J. Careaga, L.M. Villada. **Invariant-region-preserving high-order schemes for a model of reactive sedimentation**, 
*ArXiv preprint arXiv:2609.06846*, 2026, [https://arxiv.org/abs/2609.06846](https://arxiv.org/abs/2609.06846)


--------------------------

## 1. Model equations

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

where $`u(z,t)`$ is the total concentration of solids, $`\boldsymbol{c}(z,t)`$ is the vector of $`{n_{\boldsymbol{c}}}`$ solid components, $`\boldsymbol{s}(z,t)`$ is the vector $`{n_{\boldsymbol{s}}}`$ substrates, and $`u`$ is the total concentration of solids. In the default example implemented in this program, the [activated sludge reduced biokinetic model](https://www.sciencedirect.com/science/article/abs/pii/S0098135416301338), the vectors are:

$$
\begin{aligned}
\color{#076D82}\boldsymbol{c} & {\color{#076D82}= \left(X_{\rm OHO},X_{\rm U}\right)}\qquad\text{and}\qquad{\color{#822E07}
\boldsymbol{s} = \left(S_{\rm NO3},S_{\rm S},S_{\rm N2}\right)},
\end{aligned}
$$

where by components, we have:
- $`\color{#076D82}X_{\rm OHO}\text{ is the concentration of heterotrophic organics}`$
- $`\color{#076D82}X_{\rm U}\text{ is the concentration of undegradable matter}`$
- $`\color{#822E07}S_{\rm NO3}\text{ is the concentration of nitrate substrate}`$
- $`\color{#822E07}S_{\rm S}\text{ is the concentration of readily biodegradable substrate}`$
- $`\color{#822E07}S_{\rm N2}\text{ is the concentration of nitrogen substrate}`$


The numerical scheme employed to solve the PDE combines a variety of ingredients:
- High-order time approximations for ....
- High-order polynomial reconstructions and maximum-principle limiters.

--------------------------
## 2. Code structure

### Repositories and organization

The software is organized in three subfolders:
- **src**: Contains all the ".m" source codes
- **post-processing**: This is a repository where compressed auxiliary files ...
- **results**: This folder is intended to store all the results obtained from running the simulations; therefore, it is initially empty. 

### Scripts

There are currently 7 source codes in folder **scr**, these are the following:
- **main_solver.m**: Corresponds to the main source code of the solver, contains the time and spatial loops, calls all the other functions and routines and produces the output structure.
- **rs_const_functions.m**: Contains all constitutive functions of the model plus some additional functions used in the numerical approximation.
- **rs_defaults.m**: fill in optional and default fields in the parameters data structure and compute some derived constants. The input in this function is the parameters input data structure and the output is also a parameters input data structure. It is intended to be called before running the main_solver.
- **rs_diffusive_flux.m**: computes the approximation of the diffusive flux.
- **rs_ghost.m**: create ghost cells for computing boundary conditions in the extended stencil of the polynomial reconstructions.
- **rs_initial.m**: computes cell averages of the initial conditions. New initial conditions for u, c and s need to be implemented in this routine.
- **rs_reconstruct.m**: computes the second-order MUSCL and third-order CWENO reconstructions.

### Running the program

If you already have the parameters loaded in, for instance, the structure 'par', and consider the number of cells N, then you run:
```matlab
cd src
% load the parameters in the structure par
par = rs_defaults(par);
sol = main_solver(N,par);
```

--------------------------
## 3. Authorship

This Matlab-based solver has been developed by **Julio Careaga** (https://github.com/juliocareaga/) and **Juan Barajas-Calonge** (https://github.com/juanbarajascalonge/).


<br>
<br>
<br>

--------------------------
### 🔵 Input data structure

The input data is passed through a structure, that we call here par, which may contain the following fields:

| Field   |  Description  |
|:---              |:---                                                      |
| par.HH       | vessel height [m]                                            |
| par.BB       | vessel depth  [m]                                            |
| par.rho      | density of solids [kg/m^3]                                   |
| par.rhoL     | density of liquid [kg/m^3]                                   |
| par.gg       | gravity [m/s^2]                                              |
| par.v0       | hindered settling velocity at zero concentration [m/s]       |
| par.utilde   | parameter used in the Diehl flux function [kg/m^3]           |
| par.eta      | parameter used in the Diehl flux function [-]                |
| par.umax     | maximum solids concentration [kg/m^3]                        |
| par.beta     | parameter used in the compression function [-]               |
| par.uc       | critical concentration                                       |
| par.diffScheme | approximation type of the diffusion/compression terms      |
| par.Y        | yield constant in the reduced biokinetic model               |
| par.Ybar     | additional yield constant in the reduced biokinetic model    |
| par.fp       | portion that decays to non-degradable organics               |
| par.bdec     | decay rate of heterotrophic organisms [1/s]                  |
| par.mumax    | maximum growth rate [1/s]                                    |
| par.k1       | saturation constant [kg/m^3]                                 |
| par.k2       | saturation constant [kg/m^3]                                 |
| par.sigmac   | stoichiometric matrix for the solid components [-]           |
| par.sigmas   | stoichiometric matrix for the substrates       [-]           |
| par.feedMode | batch or continuous                                          |
| par.cfrac    | fraction of initial feed (taken constant; vector)            |
| par.sf_val   | feed concentration of substrates (taken constant; vector)    |
| par.bc       |  'transmissive'; gamma == 1 at every face and ghost cells    |
| par.icType   | type of initial condition: 'jump' or 'bump'                  |
| par.bumpP    | power used in smooth initial condition                       |
| par.z0       | z-value used in 'bump' [m]                                   |
| par.Rm       | radious used in 'bump' [m]                                   |
| par.ubar     | u-value used in 'bump' [kg/m^3]                              |
| par.uhat     | u-value used in 'bump' [kg/m^3]                              |
| par.sbar     | s-value used in 'bump' [kg/m^3]                              |
| par.shat     | s-value used in 'bump' [kg/m^3]                              |
| par.Tfinal   | final simulation time                                        |
| par.CFL      | CFL constant                                                 |
| par.alphaM   | Muscl parameter                                              |
| par.wenoPower| power used in the CWENO method, in the S_K factors           |
| par.posEps   | epsilon used in the CWENO method, in the S_K factors         |
| par.useLim   | Zhang-Shu limiter ('true' or 'false')                        |
| par.nsnap    | number of evenly saved snapshots of the solution             |
| par.nmon     | number of evenly saved min/max monitored values              |
| par.dtRule   | type of time step dt: 'cfl' or 'speed'                       |
| par.sameDt   | same rule (w_1 = 1/6) for all three reconstructions          |


--------------------------
### 🟢 Output data structure
The output data is obtained as a structure, that we call here sol, which contains the following fields:

| Field  |  Description  |
|:--- |:---    |
| sol.N          | Number of cells                                                        |
| sol.T          | ending time                                                            |
| sol.dt         | time step                                                              |
| sol.nSteps     | total number of time iterations                                        |
| sol.recon      | reconstruction type: 'fo', 'muscl', 'cweno3'                           |
| sol.rk         | time approximation: 'euler', 'ssprk2', 'ssprk3'                        |
| sol.useLim     | Zhang & Shu limiter: 0 ('false'), 1 ('true')                           |
| sol.diffScheme | diffusion approximation: 'none', d-method-2, d-method-4                |
|                    |                                  b-method-2, b-method-4            |
| sol.umax       | maximum solids concentration                                           |
| sol.z          | mesh nodes in the z-coordinate                                         |
| sol.tsnap      | time point snapshots                                                   |
| sol.usnap      | saved snapshots of vector solution u, size nsnap×N                     |
| sol.csnap      | saved snapshots of vector solution c, size nc×nsnap×N                  |
| sol.ssnap      | saved snapshots of vector solution s, size ns×nsnap×N                  |
| sol.tw         | monitored times for min/max values, size nmon                          |
| sol.minu       | monitored min values of u, size nmon                                   |
| sol.maxu       | monitored max values of u, size nmon                                   |
| sol.minc       | monitored min values of c, size nmon                                   | 
| sol.mins       | monitored min values of s, size nmon                                   |
| sol.sumc       | monitored sum of c, size nmon                                          |
| sol.U          | vector of solutions u, last time iteration, size 1×N                   |
| sol.C          | vector of solutions c, last time iteration, size nc×N                  |
| sol.S          | vector of solutions s, last time iteration, size ns×N                  | 
| sol.glob       | structure containing the max and min values of the variables           |





