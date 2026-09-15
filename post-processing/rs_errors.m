function [E,O] = rs_errors(sol,ref,Nref,L)
%RS_ERRORS  Approximate L^1 errors

    N = sol.N;  M = round(Nref/N);
    if M*N ~= Nref
        error('rs:mesh','N_ref = %d is not a multiple of N = %d.',Nref,N);
    end
    rU = mean(reshape(ref.U,M,N),1);
    rC = zeros(size(ref.C,1),N);  rS = zeros(size(ref.S,1),N);
    for i=1:size(ref.C,1), rC(i,:) = mean(reshape(ref.C(i,:),M,N),1); end
    for k=1:size(ref.S,1), rS(k,:) = mean(reshape(ref.S(k,:),M,N),1); end

    h = L/N;
    Q  = [sol.U ; sol.C ; sol.S];
    Qr = [rU    ; rC    ; rS   ];
    E  = h*sum(abs(Q-Qr),2).';
    O  = [];
end
