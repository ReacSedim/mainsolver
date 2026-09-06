function SOL = main_solver(N,par)
%RS_SOLVE  Run the scheme of (3.13) + Section 5 on a mesh of N cells.
%
%   SOL = RS_SOLVE(N,PAR), with PAR produced by RS_DEFAULTS, returns
%
%     SOL.U, SOL.C, SOL.S   final cell averages on the N interior cells
%     SOL.z                 cell centers
%     SOL.dt, SOL.nSteps    time step and number of steps actually taken
%     SOL.glob              struct with the space-time extrema of (7.4):
%                             minu, maxu, minc, mins, sumc
%                           where sumc = max_{j,n} |u^n_j - 1^t c^n_j|
%     SOL.tw, SOL.minu, ... the same monitors sampled at PAR.nmon instants
%     SOL.t, SOL.u, ...     PAR.nsnap snapshots for the 3D plots
%
%   Time step.  PAR.dtRule selects
%     'cfl'   the CFL condition (6.4) with lambda/w_1 of Theorem 6.4, plus the
%             parabolic term ||d||_inf/h^2 of the diffusive case:
%
%               dt = CFL * w_1 / ( (||q_f||_inf + max{v_0,L_v})/h
%                                  + ||d||_inf/h^2 + K_reac ),
%
%             with K_reac = max{Lambda_{c,1}R, Lambda_{s,1}R, Lambda_{c,2}L_R}.
%             The number of steps is then rounded up so that N_T*dt = T exactly.
%     'speed' dt = CFL * w_1 * h / max wave speed, the latter estimated from the
%             initial data.  
%
    
    dz     = par.L/N;
    par.dz = dz;
    
    if strcmpi(par.recon,'fo') && ~isfield_true(par,'sameDt')
        par.w1 = 1;
    else
        par.w1 = par.wlgl(1);
    end

    % faces z_{j+1/2} = -H + j*h, j = 0,...,N
    zf = par.zmin + (0:N)*dz;
    par.gam = ones(1,N+1);
    
    if strcmpi(par.bc,'wall'), par.gam([1 end]) = 0; end     % zero-flux walls
    
    par.faceUp   = (zf <= 0);                                % faces with q = -q_e
    par.faceDown = (zf > 0);                                 % faces with q = -q_e
    
    par.jfeed  = 2 + round(par.H/dz) + 1;                    % array index of cell j_f

    [U,C,S] = rs_initial(N,par);


    rs_const_functions;
    %% Upwind operator:    
    par.rs_upw = @(a,bl,br) max(a,0).*bl + min(a,0).*br;    

    %%=============================================================================
    %% CFL constants (begins):
    if strcmpi(par.diffScheme,'none')
        par.dmax = 0;
    else
        par.dmax = par.kap*par.rs_vhs(par.uc,par);                 % ||d||_inf (v_hs decreases)
    end
    % ---- constants of the CFL condition (6.4) ------------------------------------
    if strcmpi(par.feedMode,'batch')
        par.qfmax = 0;
    else
        par.qfmax = max(par.qf_val);
    end
    % Lipschitz constant for vhs:
    uLv    = linspace(0,par.umax,20001);
    par.Lv = max(abs(diff(par.rs_vhs(uLv,par)))./diff(uLv));

    % Reaction terms related constants:
    Rbar_c = [0, par.bdec ; 0, par.bdec];            % n_c x n_r
    Rbar_s = [par.umax*par.mumax/par.k1, 0 ; ...     % n_s x n_r
              par.umax*par.mumax/par.k2, 0 ; ...
              0,                   0];
    Lc1 = max(sum(max(-par.sigmac,0).*Rbar_c, 2));   % Lambda_{c,1} R
    Ls1 = max(sum(max(-par.sigmas,0).*Rbar_s, 2));   % Lambda_{s,1} R
    Lc2 = sum(sum(max(par.sigmac,0)));               % Lambda_{c,2} = 1 + f_p
    par.Kreac = max([Lc1, Ls1, Lc2*par.mumax]);
    
    %% ~~~ CFL constants (end) ~~~
    %%=============================================================================
    
    % ---- time step -----------------------------------------------------------------
    switch lower(par.dtRule)
    case 'cfl'    
        dtc = par.CFL*par.w1/( (par.qfmax + max(par.v0,par.Lv))/dz ...
                               + par.cflDiffFactor*par.dmax/dz^2 + par.Kreac );
    case 'speed'
        dtc = par.CFL*par.w1*dz/max_speed(U,par);
    otherwise
        error('rs:dt','Unknown dtRule ''%s''.',par.dtRule);
    end
    nSteps = ceil(par.T/dtc);
    dt     = par.T/nSteps;

    ia    = 3:N+2;
    ksnap = unique(round(linspace(0,nSteps,min(par.nsnap,nSteps+1))));
    kmon  = unique([0:min(4,nSteps), round(linspace(0,nSteps,min(par.nmon,nSteps+1)))]);

    SOL = struct();
    SOL.N = N;  SOL.T = par.T;  SOL.dt = dt;  SOL.nSteps = nSteps;
    SOL.recon = par.recon;  SOL.rk = par.rk;  SOL.useLim = par.useLim;
    SOL.diffScheme = par.diffScheme;  SOL.umax = par.umax;
    SOL.z = par.zmin + ((1:N)-0.5)*dz;

    K  = numel(ksnap);  M = numel(kmon);
    nC = size(C,1);    nS = size(S,1);
    SOL.tsnap = zeros(1,K);
    SOL.usnap = zeros(K,N);
    SOL.csnap = zeros(nC,K,N);  
    SOL.ssnap = zeros(nS,K,N);
    SOL.tw    = zeros(1,M);
    SOL.minu  = zeros(1,M);  SOL.maxu = zeros(1,M);
    SOL.minc  = zeros(1,M);  SOL.mins = zeros(1,M);  SOL.sumc = zeros(1,M);

    m = monitor(U,C,S,ia);
    g = m;                                     
    is = 1;  im = 1;
    
    [SOL,is] = store(SOL,is,0,U,C,S,ia);
    [SOL,im] = store_mon(SOL,im,0,m);


    t = 0;
    for n = 1:nSteps
        [U,C,S] = rs_step(U,C,S,t,dt,par);
        t = t + dt;

        m = monitor(U,C,S,ia);
        g = [min(g(1),m(1)), max(g(2),m(2)), min(g(3),m(3)), min(g(4),m(4)), max(g(5),m(5))];

        if any(ksnap == n), [SOL,is] = store(SOL,is,t,U,C,S,ia); end
        if any(kmon  == n), [SOL,im] = store_mon(SOL,im,t,m);    end
    end

    SOL = trim(SOL,is,im);
    SOL.U = U(ia);  SOL.C = C(:,ia);  SOL.S = S(:,ia);
    SOL.glob = struct('minu',g(1),'maxu',g(2),'minc',g(3),'mins',g(4),'sumc',g(5));
end

%==========================================================================
% Section 5
%==========================================================================
function [U,C,S] = rs_step(U,C,S,t,dt,par)
    switch lower(par.rk)
    case 'euler'                                             % forward Euler
        [LU,LC,LS] = rs_spatial(U,C,S,t,par);
        U = U + dt*LU;  
        C = C + dt*LC;  
        S = S + dt*LS;
        [U,C,S] = rs_ghost(U,C,S);

    case 'ssprk2'                                            % RK2
        [LU,LC,LS] = rs_spatial(U,C,S,t,par);
        U1 = U + dt*LU;  
        C1 = C + dt*LC;  
        S1 = S + dt*LS;
        %----------------------------------------------
        [U1,C1,S1] = rs_ghost(U1, C1, S1);
        [LU,LC,LS] = rs_spatial(U1,C1,S1,t + dt,par);
        U = 0.5*U + 0.5*(U1 + dt*LU);
        C = 0.5*C + 0.5*(C1 + dt*LC);
        S = 0.5*S + 0.5*(S1 + dt*LS);
        [U,C,S] = rs_ghost(U,C,S);
    case 'ssprk3'                                            % RK3
        [LU,LC,LS] = rs_spatial(U,C,S,t,par);
        U1 = U + dt*LU;  
        C1 = C + dt*LC;  
        S1 = S + dt*LS;
        [U1,C1,S1] = rs_ghost(U1,C1,S1);
        %----------------------------------------------
        [LU,LC,LS] = rs_spatial(U1,C1,S1,t+dt,par);
        U2 = 0.75*U + 0.25*(U1+dt*LU);
        C2 = 0.75*C + 0.25*(C1+dt*LC);
        S2 = 0.75*S + 0.25*(S1+dt*LS);
        [U2,C2,S2] = rs_ghost(U2,C2,S2);
        %----------------------------------------------
        [LU,LC,LS] = rs_spatial(U2,C2,S2,t+0.5*dt,par);
        U = U/3 + (2/3)*(U2+dt*LU);
        C = C/3 + (2/3)*(C2+dt*LC);
        S = S/3 + (2/3)*(S2+dt*LS);
        [U,C,S] = rs_ghost(U,C,S);
        %----------------------------------------------

    otherwise
        error('rs:rk','Unknown time integrator ''%s''.',par.rk);
    end
end

%==========================================================================
%  MONITORS (7.4) AND STORAGE
%==========================================================================
function m = monitor(U,C,S,ia)
    u = U(ia);  c = C(:,ia);  s = S(:,ia);
    m = [min(u), max(u), min(c(:)), min(s(:)), max(abs(u-sum(c,1)))];
end

function [SOL,is] = store(SOL,is,t,U,C,S,ia)
    SOL.tsnap(is)   = t;
    SOL.usnap(is,:) = U(ia);
    SOL.csnap(:,is,:) = reshape(C(:,ia),size(C,1),1,numel(ia));
    SOL.ssnap(:,is,:) = reshape(S(:,ia),size(S,1),1,numel(ia));
    is = is + 1;
end

function [SOL,im] = store_mon(SOL,im,t,m)
    SOL.tw(im)   = t;
    SOL.minu(im) = m(1);  SOL.maxu(im) = m(2);
    SOL.minc(im) = m(3);  SOL.mins(im) = m(4);  SOL.sumc(im) = m(5);
    im = im + 1;
end

function SOL = trim(SOL,is,im)
    a = 1:is-1;  b = 1:im-1;
    SOL.tsnap = SOL.tsnap(a);      SOL.usnap = SOL.usnap(a,:);
    SOL.csnap = SOL.csnap(:,a,:);  SOL.ssnap = SOL.ssnap(:,a,:);
    SOL.tw   = SOL.tw(b);
    SOL.minu = SOL.minu(b);  SOL.maxu = SOL.maxu(b);
    SOL.minc = SOL.minc(b);  SOL.mins = SOL.mins(b);  SOL.sumc = SOL.sumc(b);
end

%==========================================================================
%  AUXILIARIES
%==========================================================================
function s = max_speed(U,par)
    
    ua = U(3:end-2);
    ug = linspace(max(0,min(ua)-0.1), min(par.umax,max(ua)+0.1), 400);
    qm = par.qfmax;
    f  = @(x) x.*par.rs_vhs(x,par);
    dx = 1e-6;
    sU = max(abs((f(ug + dx) - f(ug - dx))/(2*dx) + qm));
    sC = max(abs(par.rs_vhs(ug,par) + qm));
    sS = max(abs(par.rho*qm - f(ug))./(par.rho-ug));
    s  = max([sU,sC,sS,1e-14]);
end

function tf = isfield_true(p,f)
    tf = isfield(p,f) && ~isempty(p.(f)) && p.(f);
end


%==========================================================================
%  SPATIAL APPROXIMATION
%==========================================================================

function [LU,LC,LS] = rs_spatial(U,C,S,t,par)
%RS_SPATIAL  Semi-discrete operators 

    
    N  = numel(U)-4;  ia = 3:N+2;  dz = par.dz;
    nC = size(C,1);   nS = size(S,1);
    LU = zeros(1,N+4);  
    LC = zeros(nC,N+4);  
    LS = zeros(nS,N+4);

    [nu,nc,ns] = rs_reconstruct(U,C,S,par);

    % traces at the N+1 faces z_{j+1/2}, j = 0,...,N
    jl = 2:N+2;  
    jr = 3:N+3;
    uL = nu(3,jl);  
    uR = nu(1,jr);
    
    % --- feed terms, and eq (3.1) 
    auxt = find(t >= par.tbreak, 1, 'last');
    qf  = par.qf_val(auxt);
    qu  = par.qu_val(auxt); 
    qe  = qf - qu;
    uf  = par.uf_val(auxt);,
    cf  = par.uf_val(auxt)*par.cfrac;
    sf  = par.sf_val;
    qfc = qu*ones(1,N+1);  
    qfc(par.faceUp) = -qe; % (3.1)

    % --- convective velocity of u (3.4) ----------------------------------
    Vf = qfc + par.gam.*par.rs_vhs(uR,par);

    % --- diffusive flux J_{j+1/2} and compression velocity E_{j+1/2};
    [Jd,Ef] = rs_diffusive_flux(U,par,uL,uR,Vf);

    % --- convective flux of u and balance (3.13a) ------------------------ 
    Fu = par.rs_upw(Vf - Ef,uL,uR);
    LU(ia) = -(Fu(2:end) - Fu(1:end-1))/dz; % + (Jd(2:end) - Jd(1:end-1))/dz;

    % --- fluxes of c -----------------------------------------------------
    if par.cFluxSplit && par.dmax > 0
        % Split form
        [phiL,phiR] = fractions(nc,uL,uR,jl,jr,nC,C,U);
        for i = 1:nC
            Fc       = par.rs_upw(Vf, nc{i}(3,jl), nc{i}(1,jr)) + par.rs_upw(-Jd, phiL(i,:), phiR(i,:));
            LC(i,ia) = -(Fc(2:end) - Fc(1:end-1))/dz;
        end
    else
        % single upwind flux with velocity V - E.
        for i = 1:nC
            Fc = par.rs_upw(Vf - Ef, nc{i}(3,jl), nc{i}(1,jr));
            LC(i,ia) = -(Fc(2:end) - Fc(1:end-1))/dz;
        end
    end

    % --- fluxes of s (3.11): velocity rho q - F + J ----------------------
    Af   = par.rho*qfc - Fu;
    denL = max(par.rho - uL,1e-10);  
    denR = max(par.rho - uR,1e-10);
    for k = 1:nS
        Fs = par.rs_upw(Af, ns{k}(3,jl)./denL, ns{k}(1,jr)./denR);
        LS(k,ia) = -(Fs(2:end) - Fs(1:end-1))/dz;
    end

    % --- feed: discrete delta in the cell j_f ----------------------------
    if qf ~= 0
        jf = par.jfeed;
        LU(jf)   = LU(jf)   + qf*uf/dz;
        LC(:,jf) = LC(:,jf) + qf*cf/dz;
        LS(:,jf) = LS(:,jf) + qf*sf/dz;
    end

    % --- reaction: G = 3 Legendre-Gauss-Lobatto quadrature (3.9) ---------
    if strcmpi(par.recon,'fo')
        rq = par.rs_reaction(C(:,ia),S(:,ia),par);   
    elseif strcmpi(par.recon,'muscl')
        rq = par.rs_reaction(C(:,ia),S(:,ia),par);   
    else
        Cg = flatnodes(nc,ia);
        Sg = flatnodes(ns,ia);
        rr = par.rs_reaction(Cg,Sg,par);
        rq = par.wlgl(1)*rr(:,1:N) + par.wlgl(2)*rr(:,N+1:2*N) + par.wlgl(3)*rr(:,2*N+1:3*N);
    end
    
    RC = par.sigmac*rq;
    LC(:,ia) = LC(:,ia) + RC;
    LS(:,ia) = LS(:,ia) + par.sigmas*rq;
    LU(ia)   = LU(ia)   + sum(RC,1);
end

%--------------------------------------------------------------------------
function [phiL,phiR] = fractions(nc,uL,uR,jl,jr,nC,C,U)
% Percentages c^{(i)}/u at the two sides of every face.  
    m    = numel(jl);
    phiL = zeros(nC,m);  phiR = zeros(nC,m);
    okL  = uL > 0;  
    okR  = uR > 0;
    
    uc_L = U(jl);   
    uc_R = U(jr);            
    
    avL  = uc_L > 0;  
    avR  = uc_R > 0;
    
    for i = 1:nC
        cL = nc{i}(3,jl);  
        cR = nc{i}(1,jr);
        
        aL = C(i,jl);      
        aR = C(i,jr);
        
        phiL(i,:) = 1/nC;                    
        s = ~okL &  avL;   
        phiL(i,s)   = aL(s)./uc_L(s);
        phiL(i,okL) = cL(okL)./uL(okL);
        phiR(i,:) = 1/nC;
        s = ~okR &  avR;   
        phiR(i,s) = aR(s)./uc_R(s);
        phiR(i,okR) = cR(okR)./uR(okR);
    end
end

%--------------------------------------------------------------------------
function G = flatnodes(cellNodes,ia)
% Stack the three nodal values of every component into one row block, in the
% order [node 1 of all cells, node 2 of all cells, node 3 of all cells].
    m = numel(cellNodes);  n = numel(ia);
    G = zeros(m,3*n);
    for i = 1:m
        Y = cellNodes{i}(:,ia);
        G(i,:) = reshape(Y.',1,[]);
    end
end









