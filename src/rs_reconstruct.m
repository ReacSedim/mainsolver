function [nu,nc,ns] = rs_reconstruct(U,C,S,par)
%RS_RECONSTRUCT  Polynomial reconstructions of Section 4, evaluated at the
%   G = 3 Legendre-Gauss-Lobatto nodes of every cell.
%
%   [NU,NC,NS] = RS_RECONSTRUCT(U,C,S,PAR) returns
%     NU     3 x (N+4)          rows = p^u_j(z_{j-1/2}), p^u_j(z_j), p^u_j(z_{j+1/2})
%     NC     {n_c} of 3 x (N+4) same for each p^{c,(i)}_j
%     NS     {n_s} of 3 x (N+4) same for each p^{s,(k)}_j
%
%   PAR.recon selects the reconstruction:
%     'fo'      First order scheme
%
%     'muscl'   Section 4.1.  
%               
%     'cweno3'  Section 4.2.  

    nc = cell(numel(C(:,1)),1);  ns = cell(numel(S(:,1)),1);
    nC = numel(nc);              nS = numel(ns);

    switch lower(par.recon)
    case 'fo'
        nu = repmat(U,3,1);
        for i=1:nC, nc{i} = repmat(C(i,:),3,1); end
        for k=1:nS, ns{k} = repmat(S(k,:),3,1); end
        return                                          % already IRP: nothing to limit

    case 'muscl'
        br = minmod_branch(U,par.alphaM);                              % (4.1)
        nu = muscl_nodes(U,br,par.alphaM);
        for i=1:nC, nc{i} = muscl_nodes(C(i,:),br,par.alphaM); end     % (4.2)
        for k=1:nS, ns{k} = muscl_nodes(S(k,:),minmod_branch(S(k,:),par.alphaM),par.alphaM); end  % (4.3)

    case 'cweno3'
        w  = cweno_weights(U,par);                                     % (4.6)-(4.7)
        nu = cweno_nodes(U,w);                                         % (4.5)
        for i=1:nC, nc{i} = cweno_nodes(C(i,:),w); end                 % (4.9)
        for k=1:nS, ns{k} = cweno_nodes(S(k,:),cweno_weights(S(k,:),par)); end % (4.10)

    otherwise
        error('rs:recon','Unknown reconstruction ''%s''.',par.recon);
    end

    if ~par.useLim, return; end
    ep = par.posEps;

    % ---- theta^solid 
    th = ones(1,numel(U));
    dU = max(nu,[],1) - U;                       % M^u_j - u_j  via (4.25)
    id = dU > 0;
    th(id) = min(th(id), (par.umax-U(id))./dU(id));

    for i = 1:nC
        dC = C(i,:) - min(nc{i},[],1);           % c^{(i)}_j - m^{c,i}_j
        id = dC > 0;
        th(id) = min(th(id), (C(i,id)-min(ep,C(i,id)))./dC(id));
    end

    th = max(0,min(1,th));

    nu = bsxfun(@plus,U,bsxfun(@times,th,bsxfun(@minus,nu,U)));            % (4.12)
    for i = 1:nC
        nc{i} = bsxfun(@plus,C(i,:), ...
                bsxfun(@times,th,bsxfun(@minus,nc{i},C(i,:))));            % (4.13)
    end

    % ---- theta^fluid 
    for k = 1:nS
        dS  = S(k,:) - min(ns{k},[],1);
        thk = ones(1,numel(U));
        id  = dS > 0;
        thk(id) = min(1, (S(k,id)-min(ep,S(k,id)))./dS(id));
        thk = max(0,min(1,thk));
        ns{k} = bsxfun(@plus,S(k,:), ...
                bsxfun(@times,thk,bsxfun(@minus,ns{k},S(k,:))));           % (4.17)
    end
end

%==========================================================================
%  MUSCL
%==========================================================================
function br = minmod_branch(Q,aM)
    M = numel(Q);  br = zeros(1,M);  j = 2:M-1;
    a = aM*(Q(j)    - Q(j-1));  
    b = 0.5*(Q(j+1) - Q(j-1));  
    c = aM*(Q(j+1)  - Q(j));
    sa = sign(a);  
    ok = (sa==sign(b)) & (sa==sign(c)) & (sa~=0);
    
    A  = abs(a);  
    B  = abs(b);  
    Cc = abs(c);
    
    bj = zeros(1,numel(j));
    
    bj(ok & A<=B  & A<=Cc) = 1;
    bj(ok & B< A  & B<=Cc) = 2;
    bj(ok & Cc< A & Cc< B) = 3;
    br(j) = bj;
end

function nodes = muscl_nodes(Q,br,aM)
    M  = numel(Q);  
    sl = zeros(1,M);  
    j = 2:M-1;
    d1 = aM*(Q(j)    - Q(j-1));  
    d2 = 0.5*(Q(j+1) - Q(j-1));  
    d3 = aM*(Q(j+1)  - Q(j));
    
    b  = br(j);
    sl(j) = (b == 1).*d1 + (b == 2).*d2 + (b == 3).*d3;
    nodes = [Q - 0.5*sl; Q; Q + 0.5*sl];
    
end

%==========================================================================
%  CWENO3
%==========================================================================
function w = cweno_weights(Q,par)
    M  = numel(Q);  j = 2:M-1;
    dL = Q(j)-Q(j-1);  dR = Q(j+1)-Q(j);
    ep = par.dz^2;
    SL = dL.^2;  SR = dR.^2;
    SC = (13/12)*(dR-dL).^2 + 0.25*(dR+dL).^2;
    aL = 0.25./(ep+SL).^par.wenoPower;
    aC = 0.50./(ep+SC).^par.wenoPower;
    aR = 0.25./(ep+SR).^par.wenoPower;
    as = aL+aC+aR;
    w = zeros(3,M);
    w(1,j) = aL./as;  w(2,j) = aC./as;  w(3,j) = aR./as;
end

function nodes = cweno_nodes(Q,w)
    M  = numel(Q);  j = 2:M-1;
    dL = Q(j) - Q(j-1);  
    dR = Q(j+1) - Q(j);
    A  = w(1,j).*dL + w(3,j).*dR + 0.5*w(2,j).*(dR + dL);
    B  = w(2,j).*(dR-dL);
    Q0 = Q(j);
    nodes = zeros(3,M);                 % columns 1 and M are never read
    nodes(1,j) = Q0 - 0.5*A + B/6;
    nodes(2,j) = Q0 - B/12;
    nodes(3,j) = Q0 + 0.5*A + B/6;
end
