function [U,C,S] = rs_initial(N,par)
%RS_INITIAL  Cell averages of the initial data.
%
%   PAR.icType = 'bump'  smooth data of Section 7.1 (Examples 1 and 3a)
%   PAR.icType = 'jump'  discontinuous data of Section 7.2 (Examples 2 and 3b)
    dz = par.L/N;
    nC = numel(par.cfrac);
    U  = zeros(1,N+4);  C = zeros(nC,N+4);


    switch lower(par.icType)
    
    case 'bump'
        zc = par.zmin + ((1:N)-0.5)*dz;
        nS = numel(par.sbar);  
        S  = zeros(nS,N+4);
        ua = cellavg(@(z) par.ubar + par.uhat*bump(z,par), zc, dz);
        U(3:N+2) = ua;
        
        for i = 1:nC, C(i,3:N+2) = par.cfrac(i)*ua; end
        
        for k = 1:nS
            S(k,3:N+2) = cellavg(@(z) par.sbar(k) + par.shat(k)*bump(z,par), zc, dz);
        end

    case 'jump'
        zl = par.zmin + (0:N-1)*dz;  zr = zl + dz;  zj = par.zjump;
        nS = numel(par.sf_val);  
        S  = zeros(nS,N+4);
        ua = avg_lin(zl,zr,zj,3.8,1.6)/dz;
        U(3:N+2) = ua;
        
        for i = 1:nC, C(i,3:N+2) = par.cfrac(i)*ua; end
        
        S(1,3:N+2) = 0.006*(min(zr,zj)-min(zl,zj))/dz;     % 0.006 if z <  zj
        S(2,3:N+2) = avg_lin(zl,zr,zj,0.12,-0.12*zj)/dz;   % 0.12(z-zj) if z >= zj
        S(3,3:N+2) = 0.006*(max(zr,zj)-max(zl,zj))/dz;     % 0.006 if z >= zj

    otherwise
        error('rs:ic','Unknown initial condition ''%s''.',par.icType);
    end


    [U,C,S] = rs_ghost(U,C,S);
end

%--------------------------------------------------------------------------
function b = bump(z,par)
    r    = (z-par.z0)/par.Rm;  
    b    = zeros(size(z));  
    m    = abs(r) < 1;
    b(m) = (1-r(m).^2).^par.bumpP;
end

function m = cellavg(fn,zc,h)
% 5-point Gauss-Legendre rule 
    gn = [-0.9061798459386640,-0.5384693101056831,0, ...
           0.5384693101056831, 0.9061798459386640]/2;
    gw = [ 0.2369268850561891, 0.4786286704993665,0.5688888888888889, ...
           0.4786286704993665, 0.2369268850561891]/2;
    m = zeros(size(zc));
    for q = 1:numel(gn), m = m + gw(q)*fn(zc+gn(q)*h); end
end

function I = avg_lin(zl,zr,zj,a,b)     % integral of (a z + b) chi_{z >= zj}
    lo = max(zl,zj);  
    hi = max(zr,zj);
    I  = 0.5*a*(hi.^2 - lo.^2) + b*(hi - lo);
end
