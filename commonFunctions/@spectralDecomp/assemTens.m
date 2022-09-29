function out=assemTens(obj,stress)
    %% Assembly of second order tensor
    % member function of spectralDecomp
    % Robin Pfefferkorn 23.9.2019
    
    %changelog
    %----------------------------------------------------------
    %v1.0 (23.9.2019) 
    %first release
    %
    %v1.1 (21.10.20) 
    %added version without perturbation using exact eigen vectors and
    %adjusted tangent expressions in case of identical eigenvalues
    
    
    %% computations
    switch obj.algorithmNr
        %------------------------------------------------------------------
        %perturbation, exact eigenvectors
        %------------------------------------------------------------------
        case {0,1}
            eigenBases = obj.eigenBases;
            out = stress(1)*eigenBases{1} + stress(2)*eigenBases{2} + stress(3)*eigenBases{3};
        otherwise
            error('algorithm not implemented')
    end
    
end