function RES = run_example2(Thours,outdir,figData)
%RUN_EXAMPLE2  Example 2 of Section 7.2: invariant-region preservation, D = 0.
%
%   RES = RUN_EXAMPLE2                 T = 3 h and T = 6 h, results into ./results
%   RES = RUN_EXAMPLE2(THOURS)         one or several final times, in hours
%   RES = RUN_EXAMPLE2(THOURS,OUTDIR)

    if nargin < 1 || isempty(Thours),  Thours  = [3 6]; end
    if nargin < 2 || isempty(outdir),  outdir  = 'results'; end
    if nargin < 3 || isempty(figData), figData = true;  end
    
    addpath(fullfile(fileparts(mfilename('fullpath')),'src'));
    addpath(fullfile(fileparts(mfilename('fullpath')),'post-processing'));    
    if ~exist(outdir,'dir'), mkdir(outdir); end

    %% PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hr = 3600;   % hour                   
    %% ---- geometry and physical constants ---------------------------------
    datos.HH     =    1.0;   % H     [m]
    datos.BB     =    3.0;   % B     [m]
    datos.AA     =  400.0;   % A     [m^2]  
    datos.rho    = 1050.0;   % rho   [kg/m^3]
    datos.rhoL   =  998.0;   % rho_L [kg/m^3]
    datos.gg     =    9.81;  % g     [m/s^2]    
    %% ---- hindered-settling velocity --------------------------------------
    datos.v0     = 1.76e-3;  % v_0   [m/s]
    datos.utilde = 3.87;     % utilde[kg/m^3]
    datos.eta    = 3.58;     % eta   [-]
    datos.umax   = 30.0;     % u_max [kg/m^3]    
    %% ---- compression: OFF in Example 2 -----------------------------------
    datos.beta   = 0.0;
    datos.uc     = datos.umax;
    datos.diffScheme = 'none';
    %% ---- bulk velocities and feed (piecewise constant in t) --------------
    datos.feedMode  = 'continuous';                     % batch or continuous
    datos.tbreak    = [0 2 4 7]*hr;                     % s
    datos.qf_val    = [1.1250 0.3250 1.6250 1.6250]/hr; % m/s    ([450, 130, 650]   / 400)/3600
    datos.qu_val    = [0.0750 0.2500 0.0875 0.1250]/hr; % m/s    ([30, 100,  35, 50]/ 400)/3600
    datos.uf_val    = [1.0    0.5    3.0    4.0];       % u_f [kg/m^3]
    datos.cfrac     = [5/7; 2/7];                       % c_f = u_f*(5/7,2/7)^t
    datos.sf_val    = [6.0e-3; 9.0e-4; 0.0];            % s_f [kg/m^3]    
    %% ---- denitrification kinetics (Section 7) ----------------------------
    datos.Y         = 0.67;      %
    datos.Ybar      = 0.172216;  %
    datos.fp        = 0.2;       %
    datos.bdec      = 6.94e-6;   % b       [1/s]
    datos.mumax     = 5.56e-5;   % mu_max  [1/s]
    datos.k1        = 5.0e-4;    % kappa_1 [kg/m^3]
    datos.k2        = 2.0e-2;    % kappa_2 [kg/m^3]    
    datos.sigmac    = [ 1.0, -1.0 ; ...
                        0.0,  datos.fp ];
    datos.sigmas    = [ -datos.Ybar , 0.0 ; ...
                        -1.0/datos.Y, 1.0 - datos.fp ; ...
                         datos.Ybar , 0.0 ];    
    %% ---- zero-flux walls at z = -H and z = B -----------------------------
    datos.bc        = 'wall';    
    %% ---- initial datum, discontinuous at z = 0.5 m -----------------------
    datos.icType    = 'jump';    %
    datos.zjump     = 0.5;       %
    %% ---- numerical parameters --------------------------------------------
    datos.Tfinal    = 9*hr;      % T = 9 h
    datos.CFL       = 0.9;       % safety factor on the CFL condition (6.4)
    datos.alphaM    = 1.4;       % alpha_M in (1,2), Section 4.1
    datos.wenoPower = 2;         % exponent in the weights (4.6)
    datos.posEps    = 1.0e-12;   % lower bound used by the limiters of Section 4.3
    datos.dtRule    = 'cfl';     % 
    datos.nsnap     = 200;       % instants stored for the 3D plots
    datos.nmon      = 200;       % instants stored for the monitors (7.4)    
    %% ---- mesh sweep -----------------------------------------
    datos.Nlist = [40 80 160 320 640];
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    datos = rs_defaults(datos);
    Nlist = datos.Nlist(:).';

    variants = { 'CWENO3',      false ; ...    % no limiters
                 'IRP-CWENO3',  true  };       % limiters of Section 4.3

    RES   = struct('Nlist',Nlist,'Thours',Thours,'names',{variants(:,1).'});
    RES.M = zeros(numel(variants(:,1)), numel(Thours), numel(Nlist),5);

    for it = 1:numel(Thours)
        for iv = 1:size(variants,1)
        
            par = datos;  
            par.recon  = 'cweno3';  
            par.rk     = 'ssprk3';
            par.T      = Thours(it) * 3600;
            par.useLim = variants{iv,2};  
            par.limitULow = false;
            
            for k = 1:numel(Nlist)
                tic;  
                %% COMPUTE THE NUMERICAL SOLUTION:  
                soln = main_solver(Nlist(k),par);  
                cpu  = toc;                
                G    = soln.glob;
                RES.M(iv,it,k,:) = [G.minu G.maxu G.minc G.mins G.sumc];
                
                %% Printing:
                fprintf(['  t = %gh  %-11s N = %4d  steps = %7d  dt = %.4f s  cpu = %6.1f s\n' ...
                         '        min u = %+10.3e  max u = %9.4f  min c = %+10.3e  ' ...
                         'min s = %+10.3e  max|u-1^t c| = %.3e\n'], ...
                         Thours(it),variants{iv,1},Nlist(k),soln.nSteps,soln.dt,cpu, ...
                         G.minu, G.maxu, G.minc, G.mins, G.sumc);
            end
        end
    end

    % ---- extra runs at T = 9 h behind Figures 7.4 and 7.5 --------------
    if figData
        fprintf('\n  data for Figures 7.4 and 7.5 (T = 9 h)\n');
        
        runs = { 'first_order', 'fo',     'euler',  false, 160 ; ...
                 'muscl',       'muscl',  'ssprk2', false, 160 ; ...
                 'irp_cweno3',  'cweno3', 'ssprk3', true , 160 ; ...
                 'irp_cweno3',  'cweno3', 'ssprk3', true , 640 };
                 
        for k = 1:size(runs,1)
            par        = datos;  
            par.recon  = runs{k,2};  
            par.rk     = runs{k,3};
            par.useLim = runs{k,4};
            par.T      = 9*3600;
            par.limitULow = false;  
            
            tic;
            %% COMPUTE THE NUMERICAL SOLUTION:  
            soln = main_solver(runs{k,5},par);
            cpu  = toc;
            
            fprintf('    %-12s N=%4d  steps=%7d  cpu=%6.1f s\n', runs{k,1}, runs{k,5}, soln.nSteps,cpu);
            rs_save(fullfile(outdir, sprintf('ex2_%s_N%d.mat', runs{k,1}, runs{k,5})), 'SOL', soln);
        end
    end

    for it = 1:numel(Thours)
        Mit = reshape(RES.M(:,it,:,:),[size(RES.M,1) numel(Nlist) 5]);
        rs_print_extrema(Nlist,Mit,RES.names, sprintf('Table 7.2 -- Example 2, D = 0, t = %g h',Thours(it)));
    end
    rs_save(fullfile(outdir,'ex2_table_7_2.mat'),'RES',RES);
    fprintf('  -> %s\n\n',fullfile(outdir, 'ex2_table_7_2.mat'));
end


