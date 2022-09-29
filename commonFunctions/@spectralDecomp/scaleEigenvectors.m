function obj=scaleEigenvectors(obj,factors)
    %% Scaling of eigenvectors with given factors
    % member function of spectralDecomp
    % Robin Pfefferkorn 23.9.2019
    
    %changelog
    %----------------------------------------------------------
    %v1.0 (23.7.2020) 
    %first release
    %
    %v1.1.....
    
    
    %% checks
    if numel(factors)~=3 || ~isfloat(factors)
        error('input must have 3 float entries')
    end
    
%     if any(abs(factors-1)>0.01)
%         keyboard
%     end
        
    %% computations
    obj.eigenVectors = obj.eigenVectors*diag(factors);
    
    % bases
    factBases = [factors(1)^2, factors(2)^2, factors(3)^2, ...
        factors(1)*factors(2), factors(2)*factors(3), factors(1)*factors(3)];
    
    bases = obj.eigenBases;   
    for i=1:numel(bases)
        bases{i} = bases{i}*factBases(i);
    end
    obj.eigenBases = bases;
    obj.eigenBasesVoigt = obj.eigenBasesVoigt*diag(factBases);
end