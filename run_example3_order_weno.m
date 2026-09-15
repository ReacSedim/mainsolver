function RES = run_example3_order(outdir)
%RUN_EXAMPLE3  Example 3 of Section 7.3: order of convergence, D > 0.

    if nargin < 1 || isempty(outdir), outdir = 'results'; end
    addpath(fullfile(fileparts(mfilename('fullpath')),'src'));
    addpath(fullfile(fileparts(mfilename('fullpath')),'post-processing'));
    if ~exist(outdir,'dir'), mkdir(outdir); end

    nameoutfile = "ex3_order_weno.mat"; %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DATA_EX3_ORDER  Example 3 of Section 7.3, first part: accuracy test WITH
    %   compression (beta > 0, so D is not identically zero).
    %
    %   The data are exactly those of Example 1 (batch regime, transmissive
    %   boundaries, C^5 compactly supported initial datum); the ONLY change is
    %   beta > 0 and the choice of u_c below.
    %
    %   CHOICE OF u_c -- the crux of this test.  The coefficient d has a JUMP at
    %   u = u_c (from 0 to kap*v_hs(u_c) > 0), so D is only Lipschitz there and
    %   D(u(.,t)) loses a derivative wherever u = u_c.  If the range of the solution
    %   contains u_c, the observable order saturates at about two no matter how fine
    %   the mesh.  Since here u_0 = ubar + uhat*psi(z) takes values in [5,7] and
    %   min u(.,t) stays close to ubar, u_c = 4 is taken: the diffusion term is
    %   ACTIVE on the whole domain and D is smooth along the solution.  RUN_EXAMPLE3_ORDER
    %   warns if min u < u_c in any run.
    %
    %   Defines the structure 'datos'; read by run_example3_order.m.    
    %% ---- geometry and physical constants ---------------------------------
    datos.HH     =    1.0;    % H     [m]
    datos.BB     =    3.0;    % B     [m]
    datos.rho    = 1050.0;    % rho   [kg/m^3]
    datos.rhoL   =  998.0;    % rho_L [kg/m^3]
    datos.gg     =    9.81;   % g     [m/s^2]    
    %% ---- hindered-settling velocity --------------------------------------
    datos.v0     = 1.76e-3;
    datos.utilde = 3.87;
    datos.eta    = 3.58;
    datos.umax   = 30.0;    
    %% ---- compression ON ---------------------------------------------------
    %  sigma_e(u) = beta*chi_{u>=u_c}*(u-u_c)  =>  sigma_e'(u) = beta,
    %  d(u) = rho v_hs(u) sigma_e'(u)/(g(rho-rho_L)) = kap v_hs(u) chi_{u>=u_c},
    %  D(u) = int_{u_c}^u d(s) ds .
    datos.beta   = 0.2;        % beta  [m^2/s^2]
    datos.uc     = 4.0;        % u_c   [kg/m^3]   see the note above    
    %  Diffusive flux used by BOTH schemes in Table 7.3; the reference solution uses
    %  datos.refDiff below.  'c2' -> (3.6),  'ho' -> (3.7)+(3.8).
    datos.diffScheme = 'c2';    
    %% ---- denitrification kinetics ----------------------------------------
    datos.Y        = 0.67;
    datos.Ybar     = 0.172216;
    datos.fp       = 0.2;
    datos.bdec     = 6.94e-6;
    datos.mumax    = 5.56e-5;
    datos.k1       = 5.0e-4;
    datos.k2       = 2.0e-2;    
    datos.sigmac   = [ 1.0, -1.0 ; ...
                       0.0,  datos.fp ];
    datos.sigmas   = [ -datos.Ybar , 0.0 ; ...
                       -1.0/datos.Y, 1.0-datos.fp ; ...
                        datos.Ybar , 0.0 ];    
    %% ---- batch regime, transmissive boundaries ---------------------------
    datos.feedMode  = 'batch';
    datos.cfrac     = [5/7; 2/7];
    datos.sf_val    = [0;0;0];
    datos.bc        = 'transmissive';    
    %% ---- initial datum, same as Example 1 --------------------------------
    datos.icType    = 'bump';
    datos.bumpP     = 6;
    datos.z0        = 1.0;
    datos.Rm        = 1.0;
    datos.ubar      = 5.0;
    datos.uhat      = 2.0;
    datos.sbar      = [0.50, 1.00, 0.10];
    datos.shat      = [0.05, 0.10, 0.01];
    datos.tbreak    = [0];   % s
    datos.qf_val    = [0];   % m/s
    datos.qu_val    = [0];   % m/s
    datos.uf_val    = [0];   % u_f [kg/m^3]    
    %% ---- numerical parameters --------------------------------------------
    datos.Tfinal    = 400.0;   % T = 0.1 h
    datos.CFL       = 0.9;
    datos.alphaM    = 1.4;     %% <-- 1.4
    datos.wenoPower = 2;
    datos.posEps    = 0.0;
    datos.nsnap     = 2;
    datos.nmon      = 2;    
    %  With diffusion the CFL condition (6.4) carries the parabolic term
    %  ||d||_inf/h^2, so the number of steps grows like N^2 and dominates the cost.
    datos.dtRule = 'cfl';
    datos.sameDt = true;    
    %% ---- convergence test -------------------------------------------------
    %  COST WARNING.  N_ref = 4*640 = 2560 with a parabolic time-step restriction is
    %  of the order of 3.7e5 steps.  Start from datos.Nlist = [40 80 160] and
    %  refFactor = 4 for a first check; the published table uses the values below.
    datos.Nlist     = [80 160 320 640];
    datos.refFactor = 4;    
    %% ---- reference solution (the same one for both schemes) --------------
    %datos.refDiff = 'ho';      % IRP-CWENO3 + SSPRK3 + the fourth-order flux (3.8)
    datos.refDiff = 'D-method-2';      % IRP-CWENO3 + SSPRK3 + the fourth-order flux (3.8)
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                  
    datos = rs_defaults(datos);
    Nlist = datos.Nlist(:).';
    Nref  = datos.refFactor * Nlist(end);

    if any(mod(Nref,Nlist) ~= 0)
        error('rs:mesh','N_ref = %d is not a multiple of every N in the list.',Nref);
    end

    schemes = { 'IRP-CWENO3-D-4', 'cweno3',  'ssprk3', true, 'D-method-4' ; ...
                'IRP-CWENO3-B-4', 'cweno3',  'ssprk3', true, 'B-method-4'};

    fprintf('\n=== Example 3: accuracy test WITH compression, t = %.3g h ===\n', datos.T/3600);
    fprintf('    beta=%.3g  u_c=%.3g\n', datos.beta,datos.uc);
    fprintf('    N = [%s],  N_ref = %d\n\n', num2str(Nlist), Nref);

    % ---- reference solution, cached ------------------------------------
    rname = fullfile(outdir,'ex3_Nref2560_D4.mat');    
    S   = load('ex3_Nref2560_D4.mat');  
    ref = S.ref;
    fprintf('  reference read from %s (%d steps)\n',rname,ref.nSteps);
    
    check_uc(ref,datos,'reference');

    % ---- mesh sweep ----------------------------------------------------
    nv  = 1 + size(datos.sigmac,1) + size(datos.sigmas,1);
    RES = struct('Nlist',Nlist,'Nref',Nref,'names',{schemes(:,1).'});
        
    for im = 1:size(schemes,1)        
        par   = setscheme(datos, schemes{im,2}, schemes{im,3}, schemes{im,4}, schemes{im,5});
        Error = zeros(numel(Nlist),nv);
        
        for k = 1:numel(Nlist)
            tic;  
            soln = main_solver(Nlist(k),par);  
            cpu  = toc;            
            Error(k,:) = rs_errors(soln, ref, Nref, datos.L);            
            check_uc(soln, datos, sprintf('N=%d',Nlist(k)));
            
            fprintf('  %-11s N=%5d  dt=%9.3e  steps=%7d  cpu=%6.1f s  e_u=%.4e\n',...
                    schemes{im,1}, Nlist(k), soln.dt, soln.nSteps, cpu, Error(k,1));
        end
        Order = zeros(size(Error));
        
        for iv = 1:nv
            Order(:,iv) = rs_order(Error(:,iv)); 
        end
        RES.E{im}     = Error;  
        RES.order{im} = Order;
    end
    
    
    rs_save(fullfile(outdir,nameoutfile),'RES',RES);
    fprintf('\n  -> %s\n\n', fullfile(outdir,nameoutfile));
end


%--------------------------------------------------------------------------
function p = setscheme(p,recon,rk,useLim,diffScheme)
    p.recon      = recon;  
    p.rk         = rk;  
    p.useLim     = useLim;  
    p.diffScheme = diffScheme;
    p.limitULow  = true;               % D ~= 0, see RS_RECONSTRUCT
    p            = rs_defaults(p);     % refresh ||d||_inf for this flux choice
end

function check_uc(sol,par,label)
    if sol.glob.minu < par.uc
        warning('rs:uc',['[%s] min u = %.4f < u_c = %.4f: the solution crosses ' ...
                'the jump of d and the observed order will be limited to ~2.'], ...
                label,sol.glob.minu,par.uc);
    end
end
