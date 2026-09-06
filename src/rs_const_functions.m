%%==================================================================================================
%%==================================================================================================
%%  CONSTITUTIVE FUNCTIONS:
%%  
%%  The functions in this file can also be programmed in separate files, mind that in that case
%%  you have to pass the right arguments, including the parameters par.
%%==================================================================================================
%%==================================================================================================


%global par.rs_vhs par.rs_psi par.rs_dfun par.rs_bfun par.growth_mu par.rs_reaction

bnd0 = @(u) max(0,min(par.umax,u));
bndc = @(u) min(max(u, par.uc), par.umax);

%%==================================================================================================
%% Hindered settling function (Diehl type)
%%==================================================================================================

auxvhs   = @(u) (u/par.utilde).^par.eta;
par.vinf = par.v0/(1 + auxvhs(par.umax));     % => v_hs(u_max) = 0

par.rs_vhs   = @(u,par) max(0, par.v0./(1 + auxvhs(bnd0(u))) - par.vinf);

%%==================================================================================================
%% Compression function D(u)
%%==================================================================================================

par.kap  = par.rho*par.beta/(par.gg*(par.rho - par.rhoL));   % d(u) = kap*v_hs(u)
par.abet = 1/par.eta;                                        % a = 1/eta
par.cbet = (par.utilde/par.eta)*pi/sin(pi/par.eta);          % (utilde/eta)*B(a,1-a)

%% par.rs_psi  Primitive of 1/(1+(s/utilde)^eta), used to evaluate D in closed form.
par.rs_psi = @(s,par) par.cbet*betainc((auxvhs(s))./(1 + auxvhs(s)), par.abet, 1 - par.abet);

par.Psic = par.rs_psi(par.uc,par);

par.rs_dfun = @(u,par) (u > par.uc).*par.kap.*(par.v0*(par.rs_psi(bndc(u),par) - par.Psic) - par.vinf*(bndc(u) - par.uc) );

%%==================================================================================================
%% Compression-related function B(u)
%%==================================================================================================

par.rs_bfun = @(u,par) (u > par.uc).*par.kap.*( (par.v0 - par.vinf)*log(bndc(u)/par.uc) ...
                    - (par.v0/par.eta)*log((1 + auxvhs(bndc(u)))/(1 + auxvhs(par.uc))));



%%==================================================================================================
%% Function mu (bacterial growth):
par.growth_mu = @(ss1, ss2) par.mumax .* (ss1./(par.k1 + ss1)) .* (ss2./(par.k2 + ss2));

%% Vector of nonlinear reaction functions
par.rs_reaction = @(C,S,par) [C(1,:).* par.growth_mu(S(1,:), S(2,:));  C(1,:).*par.bdec];

%%==================================================================================================

