function rs_save(fname,vname,value)
%RS_SAVE  save() wrapper that produces a MAT-file readable by both MATLAB and
%   Octave (Octave needs the explicit '-v7' flag for that).
    eval([vname ' = value;']);                                       
    if exist('OCTAVE_VERSION','builtin') ~= 0
        save('-v7',fname,vname);
    else
        save(fname,vname);
    end
end
