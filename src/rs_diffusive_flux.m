function [Jd,Ef] = rs_diffusive_flux(U,par,uL,uR,Vf)
%RS_DIFFUSIVE_FLUX  Diffusive flux J_{j+1/2} at the N+1 faces, and the
%   compression velocity E_{j+1/2} used by the equation for c.
%
%   [JD,EF] = RS_DIFFUSIVE_FLUX(U,PAR,UL,UR,VF) returns the two row vectors of
%   length N+1.  UL, UR are the reconstructed traces of u on each side of every
%   face and VF is the convective velocity of (3.4).
%
%   Face m (m = 1,...,N+1) is z_{j+1/2} with j = m-1, and cell j sits at array
%   index j+2, so the four-point stencil is made of indices m,...,m+3.
%
%   PAR.diffScheme selects the approximation:
%
%     'none'   J = 0 
%
%     'c2'     Eq. (3.6).  Two-point central difference of D
%
%     'ho'     Eqs. (3.7)-(3.8).  Fourth-order staggered difference of D built on
%              the auxiliary cell-centre values ubar_j = u_j - (u_{j+1}-2u_j+u_{j-1})/24
%
%
%     'bhat'   Eq. (3.10).  
%

 
    N = numel(U) - 4;  m = 1:N+1;

    if strcmpi(par.diffScheme,'none')
        Jd = zeros(1,N+1);  Ef = zeros(1,N+1);  return
    end

    switch lower(par.diffScheme)
    case 'd-method-2'
        Dc = par.rs_dfun(U,par);                        
        Jd = par.gam.*( Dc(m+2) - Dc(m+1) )/par.dz;
        Ef = quotientE(Jd,uL,uR,par);                  
    case 'd-method-4'        
        ut = aux_reconstruction(U,N);                         
        Dp = par.rs_dfun(ut,par);                             
        Jd = par.gam.*( Dp(m) - 27*Dp(m+1) + 27*Dp(m+2) - Dp(m+3) )/(24*par.dz);
        Ef = quotientE(Jd,uL,uR,par);                          
    case 'b-method-2'
        Jd = 0.0;
        Bp = par.rs_bfun(U,par);
        Ef = par.gam.*(Bp(m+2) - Bp(m+1) )/par.dz;
    case 'b-method-4'
        Jd = 0.0;                             
        ut = aux_reconstruction(U,N); 
        Bp = par.rs_bfun(ut,par);
        Ef = par.gam.*( Bp(m) - 27*Bp(m+1) + 27*Bp(m+2) - Bp(m+3) )/(24*par.dz);         
    otherwise
        error('rs:diff','Unknown diffusive discretization ''%s''.',par.diffScheme);
    end
end

%--------------------------------------------------------------------------
function ut = aux_reconstruction(U,N)
% (3.7): fourth-order cell-centre values from the cell averages.
    ut = U;  jj = 2:N+3;
    ut(jj) = U(jj) - (U(jj+1) - 2*U(jj) + U(jj-1))/24;
end

function Ef = quotientE(Jd,uL,uR,par)
% (3.14): E = J/ubar with ubar the mean of the two traces.  J vanishes wherever
% D does, and the cut at u_c prevents a division by ~0.
    Ef = Jd./max(0.5*(uL + uR),par.uc);
end
