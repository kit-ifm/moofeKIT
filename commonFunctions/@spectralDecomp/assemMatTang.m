function out=assemMatTang(obj,stress,eMat,factStress)
    %% Assembly of material tangent tensor
    % member function of spectralDecomp
    % always in voigt notation
    % Robin Pfefferkorn 23.9.2019
    
    %changelog
    %----------------------------------------------------------
    %v1.0 (23.9.2019) 
    %first release
    %
    %v1.1 (21.10.20) 
    %added version without perturbation using exact eigen vectors and
    %adjusted tangent expressions in case of identical eigenvalues
    
    %% input
    if nargin<=3
        factStress = -2;
    end
    
    
    %% computations
    lam = obj.principalStretch;
    
    %bases  
    basesVoigt = obj.eigenBasesVoigt;
    M11_v = basesVoigt(:,1);
    M22_v = basesVoigt(:,2);
    M33_v = basesVoigt(:,3);
    M12_v = basesVoigt(:,4);
    M23_v = basesVoigt(:,5);
    M13_v = basesVoigt(:,6);
    
    %% principal stress contributions
    if obj.spatialFormulation
        principalStr = factStress*stress;
    else
        principalStr = factStress*(lam.^-2).*stress;
    end
    
    %% geometric contributions
    switch obj.algorithmNr
        %------------------------------------------------------------------
        %Perturbation
        %------------------------------------------------------------------
        case 0
            if obj.spatialFormulation
                %geometric parts
                geom = zeros(3,1);
                geom(1) = (stress(2)*lam(1)^2-stress(1)*lam(2)^2)/(lam(2)^2-lam(1)^2);
                geom(2) = (stress(3)*lam(2)^2-stress(2)*lam(3)^2)/(lam(3)^2-lam(2)^2);
                geom(3) = (stress(3)*lam(1)^2-stress(1)*lam(3)^2)/(lam(3)^2-lam(1)^2);
            else
                %geometric parts
                geom = zeros(3,1);
                geom(1) = (stress(2)-stress(1))/(lam(2)^2-lam(1)^2);
                geom(2) = (stress(3)-stress(2))/(lam(3)^2-lam(2)^2);
                geom(3) = (stress(3)-stress(1))/(lam(3)^2-lam(1)^2);
            end
        %------------------------------------------------------------------
        %exact eigenvectors
        %------------------------------------------------------------------
        case 1
            threshold = obj.perturbThreshold;
            if obj.spatialFormulation
                %geometric parts
                geom = zeros(3,1);
                geom(1) = switchGeomPartSpatial(lam,stress,eMat,1,2,threshold);
                geom(2) = switchGeomPartSpatial(lam,stress,eMat,2,3,threshold);
                geom(3) = switchGeomPartSpatial(lam,stress,eMat,1,3,threshold);
            else
                %geometric parts
                geom = zeros(3,1);
                geom(1) = switchGeomPartMaterial(lam,stress,eMat,1,2,threshold);
                geom(2) = switchGeomPartMaterial(lam,stress,eMat,2,3,threshold);
                geom(3) = switchGeomPartMaterial(lam,stress,eMat,1,3,threshold);
            end
        otherwise
            error('algorithm not implemented')
    end
    
    %% assembly
    switch obj.algorithmNr
        %------------------------------------------------------------------
        %perturbation, exact eigenvectors (eigenvectors known)
        %------------------------------------------------------------------
        case {0,1}
            out = ...
                eMat(1,1)*(M11_v*M11_v') + eMat(1,2)*(M11_v*M22_v') + eMat(1,3)*(M11_v*M33_v') + ...
                eMat(2,1)*(M22_v*M11_v') + eMat(2,2)*(M22_v*M22_v') + eMat(2,3)*(M22_v*M33_v') + ...
                eMat(3,1)*(M33_v*M11_v') + eMat(3,2)*(M33_v*M22_v') + eMat(3,3)*(M33_v*M33_v') + ...
                principalStr(1)*(M11_v*M11_v')+...
                principalStr(2)*(M22_v*M22_v')+...
                principalStr(3)*(M33_v*M33_v')+...
                geom(1)*(M12_v*M12_v')+...
                geom(2)*(M23_v*M23_v')+...
                geom(3)*(M13_v*M13_v');
        otherwise
            error('algorithm not implemented')
    end
    
end
function out = switchGeomPartSpatial(lam,stress,eMat,a,b,threshold)
    if abs(lam(a)-lam(b))/lam(a)>threshold
        %distinct eigenvalues
        out = (stress(b)*lam(a)^2-stress(a)*lam(b)^2)/(lam(b)^2-lam(a)^2);
    else
        %coincident eigenvalues
        out = 1/2*(eMat(b,b)-eMat(a,b))-stress(a);
    end
end
function out = switchGeomPartMaterial(lam,stress,eMat,a,b,threshold)
    if abs(lam(a)-lam(b))/lam(a)>threshold
        %distinct eigenvalues
        out = (stress(b)-stress(a))/(lam(b)^2-lam(a)^2);
    else
        %coincident eigenvalues
        out = 1/2*(eMat(b,b)-eMat(a,b));
    end
end