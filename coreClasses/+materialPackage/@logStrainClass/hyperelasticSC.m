function [W,S,Sv,CCv,errMat] = hyperelasticSC(materialObject,mapVoigtObject,C)
    errMat = 0; %no error
    %returns potential energy, stress and material derivatives for
    %hyperelastic materials in SC formulation
    
    %input:
    %-C: right Cauchy-Green
    mapVoigt = mapVoigtObject.mapVoigt;
    dimension = mapVoigtObject.dimension;
    
    % unity tensor
    I = eye(3);
    
    %material data
    mu = materialObject.mu;
    kappa = materialObject.kappa;
    
    %kinematics
    C = 0.5*(C+C.');
    
    %settings spectral decomposition
    specDecomp = spectralDecomp(materialObject.specDecompAlgo,false);
    specDecomp.perturbThreshold = materialObject.perturbThreshold;
    
    %spectral decomposition
    specDecomp.decomposition(C,mapVoigt);
    lam = specDecomp.principalStretch;
    J = lam(1)*lam(2)*lam(3);
    c = J^2;
    if any(lam>1e4) || any(lam<1e-4)
        errMat = 1;
    end
    
    %strain energy and derivatives
    W = 1/8*kappa*log(c)^2 + mu*sum((log(lam)-1/6*log(c)).^2);
    %derivatives
    DW_logL = 2*mu*(log(lam)-1/6*log(c));
    DW_logc = 1/4*kappa*log(c) - 1/3*mu*sum(log(lam)-1/6*log(c));
    
    %stress tensor (2.PK)
    s = (lam.^-2).*(DW_logL + 2*DW_logc);
    S = specDecomp.assemTens(s);
    
    % elasticity tensor
    eVol = kappa*ones(3,3);
    eIso = 2*mu*(I-1/3*ones(3,3));
    ETot = ((lam.^-2)*(lam.^-2).').*(eIso+eVol);
    CCv = specDecomp.assemMatTang(s,ETot);
       
    Sv = S(mapVoigt);
end