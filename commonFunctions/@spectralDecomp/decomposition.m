function obj=decomposition(obj,b,mapVoigt)
    %% Spectral Decomposition Algorithms
    % member function of spectralDecomp
    % Robin Pfefferkorn 23.9.2019
    %
    % INPUT:
    % b...right or left cauchy green tensor (depending on spatial or
    % material formulation)
    % mapVoigt...vector containing voigt map information
    %
    %changelog
    %----------------------------------------------------------
    %v1.0 (23.9.2019) 
    %first release
    %
    %v1.1 (21.10.20) 
    %added version without perturbation using exact eigen vectors and
    %adjusted tangent expressions in case of identical eigenvalues
    
    %% check input
    if size(b,1)~=3 || size(b,2)~=3
        error('only 3x3 matrices supported')
    end
    
    temp=nargin;
    if temp==3
        voigt=true;
        if isrow(mapVoigt)
            mapVoigt = mapVoigt';
        elseif ~iscolumn(mapVoigt)
            error('mapVoigt must be a column vector')
        end
    elseif temp==2
        voigt=false;
    end
    
    %% compute eigenvalues and vectors
    switch obj.algorithmNr
        %------------------------------------------------------------------
        %Perturbation
        %------------------------------------------------------------------
        case 0
            %matlab decomp
            [nEig,lam2] = eig(b);
            lam = sqrt(diag(lam2));
    
            %perturbations
            tol = obj.perturbThreshold*max(abs(lam));
            pertubed = tol;
            if abs(lam(1)-lam(2)) < tol
                lam(1) = (1-pertubed)*lam(1);
                lam(2) = (1+pertubed)*lam(2);
                lam(3) = 1/((1-pertubed)*(1+pertubed))*lam(3);
            elseif abs(lam(1)-lam(3)) < tol
                lam(1) = (1-pertubed)*lam(1);
                lam(2) = 1/((1-pertubed)*(1+pertubed))*lam(2);
                lam(3) = (1+pertubed)*lam(3);
            elseif abs(lam(2)-lam(3)) < tol
                lam(1) = 1/((1-pertubed)*(1+pertubed))*lam(1);
                lam(2) = (1-pertubed)*lam(2);
                lam(3) = (1+pertubed)*lam(3);
            end
        %------------------------------------------------------------------
        %exact eigenvectors
        %------------------------------------------------------------------
        case 1  
            %matlab decomp
            [nEig,lam2] = eig(b);
            lam = sqrt(diag(lam2));
        otherwise
            error('algorithm not implemented')
    end
    
    %% compute bases and save results
    switch obj.algorithmNr
        %------------------------------------------------------------------
        %perturbation, exact eigenvectors (eigenvectors known)
        %------------------------------------------------------------------
        case {0,1}
            %save results
            obj.principalStretch = lam;
            obj.eigenVectors = nEig;

            %bases
            bases = cell(6,1);
            bases{1} = nEig(:,1)*nEig(:,1).';
            bases{2} = nEig(:,2)*nEig(:,2).';
            bases{3} = nEig(:,3)*nEig(:,3).';
            bases{4} = nEig(:,1)*nEig(:,2).'+nEig(:,2)*nEig(:,1).';
            bases{5} = nEig(:,2)*nEig(:,3).'+nEig(:,3)*nEig(:,2).';
            bases{6} = nEig(:,1)*nEig(:,3).'+nEig(:,3)*nEig(:,1).';
            obj.eigenBases = bases;
            
            %bases - voigt
            if voigt
                basesVoigt = zeros(numel(mapVoigt),6);
                basesVoigt(:,1) = bases{1}(mapVoigt);
                basesVoigt(:,2) = bases{2}(mapVoigt);
                basesVoigt(:,3) = bases{3}(mapVoigt);
                basesVoigt(:,4) = bases{4}(mapVoigt);
                basesVoigt(:,5) = bases{5}(mapVoigt);
                basesVoigt(:,6) = bases{6}(mapVoigt);
                obj.eigenBasesVoigt = basesVoigt;
            end
        otherwise
            error('algorithm not implemented')
    end
    
end