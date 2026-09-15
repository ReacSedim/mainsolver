function rs_print_errors(R,title,cols)
%RS_PRINT_ERRORS  Print L^1 errors and observed orders

    vn = {'u','c^(1)','c^(2)','s^(1)','s^(2)','s^(3)'};
    if nargin < 3 || isempty(cols), cols = 1:size(R.E{1},2); end

    fprintf('\n  %s\n',title);
    for im = 1:numel(R.names)
        fprintf('\n  --- %s ---\n',R.names{im});
        fprintf('  %6s','N');
        
        for iv = cols, fprintf(' | %11s %5s',['e_' vn{iv}],'R'); end
        
        fprintf('\n  %s\n',repmat('-',1,6+numel(cols)*20));
        
        for k = 1:numel(R.Nlist)
            fprintf('  %6d',R.Nlist(k));
            for iv = cols
                o = R.order{im}(k,iv);
                
                if isnan(o), os = '--'; else, os = sprintf('%.2f',o); end
                
                fprintf(' | %11.4e %5s',R.E{im}(k,iv),os);
            end
            fprintf('\n');
        end
    end
    fprintf('\n');
end
