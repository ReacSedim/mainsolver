function par = rs_defaults(d)
%RS_DEFAULTS  Fill in optional fields and compute derived constants.
%
%   PAR = RS_DEFAULTS(D) takes the raw data structure D produced by one of the
%   files in examples/ and returns the parameter structure PAR used by
%   RS_SOLVE.  It adds
%     * the geometry shortcuts  L, zmin;
%     * the constants of  v_hs  and  D  (Section 7 of the paper);
%     * the constants entering the CFL condition (6.4);
%     * default values for every optional numerical switch.
%
%   Fields of D (see examples/data_*.m for complete, commented instances):
%
%   Geometry and physics
%     HH, BB          tank height above / depth below the feed level [m]
%     rho, rhoL, gg   solid density, liquid density, gravity
%     v0, utilde, eta, umax        parameters of v_hs, Eq. (7) of Section 7
%     beta, uc        compression parameters; beta = 0 switches D off
%
%   Kinetics (Section 7)
%     Y, Ybar, fp, bdec, mumax, k1, k2, sigmac, sigmas
%
%   Bulk flows.  feedMode = 'batch'      -> q = 0, no feed (Examples 1 and 3a)
%                feedMode = 'continuous' -> piecewise-constant-in-time data
%                                           tbreak, qf_val, qu_val, uf_val,
%                                           cfrac, sf_val (Examples 2 and 3b)
%
%   Initial data.  icType = 'bump' -> u_0 = ubar + uhat*psi(z), Section 7.1
%                  icType = 'jump' -> discontinuous data of Section 7.2
%
%   Numerics
%     recon       'fo'    | 'muscl'      | 'cweno3'      reconstruction, Section 4
%     rk          'euler' | 'ssprk2'     | 'ssprk3'      time stepping, Section 5
%     useLim      true/false             |               scaling limiters, Section 4.3
%     diffScheme  'none'  | 'D-method-2' | 'D-method-4' 
%                         | 'B-method-2' | 'B-method-4' 
%     bc          'wall'  | 'transmissive'            gamma at the two end faces
%     dtRule      'cfl'   | 'speed'                   see RS_SOLVE
%     CFL, alphaM, wenoPower, posEps, Tfinal
%
%   Part of the code accompanying
%     J. Barajas-Calonge, J. Careaga, L.-M. Villada, 
%     "Invariant-region-preserving high-order schemes for a model of reactive sedimentation".
%---------------------------------------------------------------------------------------------------
    par = d;
    % ---- geometry ------------------------------------------------------
    par.H    = d.HH;  
    par.B    = d.BB;
    par.L    = d.HH + d.BB;
    par.zmin = -d.HH;
    par.T    = d.Tfinal;

    % ---- G = 3 Legendre-Gauss-Lobatto rule on [-1/2,1/2], Section 3.1 ---
    par.wlgl = [1/6, 2/3, 1/6];
    par.w1   = par.wlgl(1);

    % ---- optional switches ---------------------------------------------
    par = setdef(par,'recon'     ,'cweno3');
    par = setdef(par,'rk'        ,'ssprk3');
    par = setdef(par,'useLim'    ,true    );
    par = setdef(par,'diffScheme','none'  );
    par = setdef(par,'bc'        ,'wall'  );
    par = setdef(par,'feedMode'  ,'continuous');
    par = setdef(par,'icType'    ,'jump'  );
    par = setdef(par,'dtRule'    ,'cfl'   );
    par = setdef(par,'CFL'       ,0.9     );
    par = setdef(par,'alphaM'    ,1.4     );
    par = setdef(par,'wenoPower' ,2       );
    par = setdef(par,'posEps'    ,0.0     );
    par = setdef(par,'nsnap'     ,200     );
    par = setdef(par,'nmon'      ,200     );
    par = setdef(par,'beta'      ,0.0     );
    par = setdef(par,'uc'        ,d.umax  );

    % The scaling limiter (4.14) bounds p^u from above by u_max and each
    % p^{c,(i)} from below by 0; by (4.11) this also gives p^u >= 0.  With
    % compression that implication is lost (the flux of u carries +J while the
    % flux of c carries -E c), so p^u >= 0 is imposed explicitly as well.
    par = setdef(par,'limitULow',~strcmpi(par.diffScheme,'none'));

    % Form of the numerical flux of c.  false = (3.10) as published; true = the
    % split form of RS_SPATIAL, which restores 1^t c = u to round-off when the
    % compression term is active.  Default false so that the runners reproduce
    % the published tables.
    par = setdef(par,'cFluxSplit',false);

    % Extra constraint on the scaling limiter needed by the invariant-region
    % result WITH compression: p^{c,i}_j(z)/p^u_j(z) <= betaPsi * c^{(i)}_j/u_j
    % wherever u_j > u_c.  Only meaningful together with the split flux.
    par = setdef(par,'betaPsi' ,2);
    par = setdef(par,'limitPsi',par.cFluxSplit && ~strcmpi(par.diffScheme,'none'));

    % Multiplier of ||d||_inf/h^2 in the CFL condition.  1 reproduces the time
    % step used for the published tables; 2*betaPsi is the value for which the
    % invariant-region theorem with compression is proved.
    par = setdef(par,'cflDiffFactor',1);

   
end

%--------------------------------------------------------------------------
function p = setdef(p,f,v)
    if ~isfield(p,f) || isempty(p.(f)), p.(f) = v; end
end


