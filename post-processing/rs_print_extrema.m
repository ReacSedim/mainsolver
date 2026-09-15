function rs_print_extrema(Nlist,M,names,header)
%RS_PRINT_EXTREMA  Print the space-time extrema (7.4)

    fprintf('\n  %s\n',header);
    for ib = 1:numel(names)
        fprintf('\n  --- %s ---\n',names{ib});
        fprintf('  %6s %14s %10s %14s %14s %13s\n', 'N','min u','max u','min c','min s','max|u-1^t c|');
        fprintf('  %s\n',repmat('-',1,76));
        for k = 1:numel(Nlist)
            m = squeeze(M(ib,k,:));
            fprintf('  %6d %14s %10.4f %14s %14s %13.2e\n', ...
                    Nlist(k),zstr(m(1)),m(2),zstr(m(3)),zstr(m(4)),m(5));
        end
    end
    fprintf('\n');
end

function s = zstr(x)
    if x == 0, s = '0.0'; else, s = sprintf('%+.4e',x); end
end
