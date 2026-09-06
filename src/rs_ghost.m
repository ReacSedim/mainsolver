function [U,C,S] = rs_ghost(U,C,S)
%RS_GHOST  Ghost cells for boundary conditions
    N = numel(U)-4;
    U(1:2)   = U(3);          U(N+3:N+4)   = U(N+2);
    C(:,1:2) = C(:,[3 3]);    C(:,N+3:N+4) = C(:,[N+2 N+2]);
    S(:,1:2) = S(:,[3 3]);    S(:,N+3:N+4) = S(:,[N+2 N+2]);
end
