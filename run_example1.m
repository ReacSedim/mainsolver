function RES = run_example1(outdir)
%RUN_EXAMPLE1  Example 1 of Section 7.1: accuracy test without compression.
%
%   RES = RUN_EXAMPLE1            results into ./results
%   RES = RUN_EXAMPLE1(OUTDIR)    results into OUTDIR

    if nargin < 1 || isempty(outdir), outdir = 'results'; end
    
    addpath(fullfile(fileparts(mfilename('fullpath')),'src'));
    addpath(fullfile(fileparts(mfilename('fullpath')),'post-processing'));
    
    if ~exist(outdir,'dir'), mkdir(outdir); end

    %% PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hr = 3600.0;
    %% ---- geometry and physical constants (Section 7) ---------------------
    datos.HH     =    1.0;      % H     [m]
    datos.BB     =    3.0;      % B     [m]
    datos.rho    = 1050.0;      % rho   [kg/m^3]
    datos.rhoL   =  998.0;      % rho_L [kg/m^3]
    datos.gg     =    9.81;     % g     [m/s^2]
    %% ---- hindered-settling velocity --------------------------------------
    datos.v0     = 1.76e-3;     % v_0     [m/s]
    datos.utilde = 3.87;        % utilde  [kg/m^3]
    datos.eta    = 3.58;        % eta     [-]
    datos.umax   = 30.0;        % u_max   [kg/m^3]
    %% ---- compression: OFF in Example 1 -----------------------------------
    datos.beta   = 0.0;         % 
    datos.uc     = datos.umax;  % 
    datos.diffScheme = 'none';  % 
    %% ---- denitrification kinetics (Section 7) ----------------------------
    datos.Y      = 0.67;        % 
    datos.Ybar   = 0.172216;    % 
    datos.fp     = 0.2;         % 
    datos.bdec   = 6.94e-6;     % b       [1/s]
    datos.mumax  = 5.56e-5;     % mu_max  [1/s]
    datos.k1     = 5.0e-4;      % kappa_1 [kg/m^3]
    datos.k2     = 2.0e-2;      % kappa_2 [kg/m^3]
    datos.sigmac = [ 1.0,  -1.0 ; ...
                     0.0,   datos.fp ];  
    datos.sigmas = [ -datos.Ybar ,  0.0 ; ...
                     -1.0/datos.Y,  1.0 - datos.fp ; ...
                      datos.Ybar ,  0.0 ];
    %% ---- bulk flows: batch -----------------------------------------------
    datos.feedMode  = 'batch';    % batch or continuous
    datos.cfrac     = [5/7; 2/7]; %   
    datos.sf_val    = [0; 0; 0];  % 
    %% ---- numerical boundary conditions -----------------------------------
    datos.bc        = 'transmissive';  %  'transmissive': gamma == 1 at every face and ghost cells
    %% ---- initial condition, Section 7.1 --------------------------------------
    datos.icType    = 'bump';              %
    datos.bumpP     = 6;                   % C^5 datum
    datos.z0        = 1.0;                 % z_0  [m]
    datos.Rm        = 1.0;                 % R    [m]
    datos.ubar      = 5.0;                 % [kg/m^3]
    datos.uhat      = 2.0;                 % [kg/m^3]
    datos.sbar      = [0.50, 1.00, 0.10];  % [kg/m^3]
    datos.shat      = [0.05, 0.10, 0.01];  % [kg/m^3]
    %% ---- numerical parameters --------------------------------------------
    datos.Tfinal    = 400.0;               % T = 0.1 h  [s]
    datos.CFL       = 0.9;                 % CFL constant
    datos.alphaM    = 1.4;                 % Muscl parameter
    datos.wenoPower = 2;                   %
    datos.posEps    = 0.0;                 %
    datos.useLim    = true;                % Zhang-Shu Limiter
    datos.nsnap     = 2;                   %
    datos.nmon      = 2;                   %
    datos.dtRule    = 'speed';             %
    datos.sameDt    = true;                % same rule (w_1 = 1/6) for all three reconstructions
    %% ---- mesh sweep -------------------------------------------------------
    datos.Nlist     = [160 320 640 1280 2560];    % N_l = 10*2^l, l = 4,...,8
    datos.refFactor = 4;                          % N_ref = 4*2560 = 10240    
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    datos = rs_defaults(datos);  %% set default parameters
    Nlist = datos.Nlist(:).';
    Nref  = datos.refFactor*Nlist(end);

    schemes = { 'First-order', 'fo',     'euler',  false ; ...
                'MUSCL',       'muscl',  'ssprk2', false ; ...
                'IRP-CWENO3',  'cweno3', 'ssprk3', true  };
    fprintf('\n')
    fprintf('=== Example 1: accuracy test, D = 0, t = %.3g h ===\n',datos.T/hr);
    fprintf('    N = [%s],  N_ref = %d\n\n',num2str(Nlist),Nref);

    % ---- reference solution (IRP-CWENO3 + SSPRK3) ----------------------
    load('ex1_Nref10240_t400.mat') %% load reference solution, structure: ref
    tic;  
    fprintf('  reference N=%d: %d steps, %.1f s\n',Nref, ref.nSteps, toc);
    % ---- mesh sweep ----------------------------------------------------
    nv  = 1 + size(datos.sigmac,1) + size(datos.sigmas,1); % number of variables
    RES = struct('Nlist',Nlist,'Nref',Nref,'names',{schemes(:,1).'});
   
    for im = 1:size(schemes,1)
        
        %--------------------- Reconstruct    Time-scheme    Limiter ------ 
        par = setscheme(datos, schemes{im,2}, schemes{im,3}, schemes{im,4});
        Error  = zeros(numel(Nlist),nv);
        
        for k = 1:numel(Nlist)
            tic;  
            %% COMPUTE THE NUMERICAL SOLUTION:
            soln = main_solver(Nlist(k),par); 
            cpu  = toc;
            Error(k,:) = rs_errors(soln, ref, Nref,datos.L);
            fprintf('  %-11s N=%5d  steps=%7d  cpu=%6.1f s  e_u=%.4e\n', schemes{im,1}, Nlist(k), soln.nSteps, cpu,Error(k,1));
        end
        
        %% Order of convergence:
        Order = zeros(size(Error));
        for iv = 1:nv
            Order(:,iv) = rs_order(Error(:,iv)); 
        end
        RES.E{im}     = Error; 
        RES.order{im} = Order;
    end

    rs_print_errors(RES, sprintf('Table 7.1 -- Example 1: L^1 errors (7.2) and orders (7.3), t = %g s',datos.T));
    save(fullfile(outdir,'ex1_results.mat'),'RES');
    fprintf('  -> %s\n\n',fullfile(outdir,'ex1_results.mat'));
end
%--------------------------------------------------------------------------
function p = setscheme(p,recon,rk,useLim)
    p.recon     = recon;  
    p.rk        = rk;  
    p.useLim    = useLim;
    p.limitULow = false;   % D = 0.
end


