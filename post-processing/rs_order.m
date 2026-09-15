function o = rs_order(e)
%RS_ORDER  Observed convergence rates
    e = e(:);  o = NaN(size(e));  o(2:end) = log2(e(1:end-1)./e(2:end));
end
