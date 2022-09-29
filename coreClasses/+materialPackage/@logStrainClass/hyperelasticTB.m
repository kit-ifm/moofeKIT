function [W,tau,tauv,ccMatv,errMat] = hyperelasticTB(materialObject,mapVoigtObject,F)
    errMat = 0; %no error
    %returns potential energy, stress and material derivatives for
    %hyperelastic materials in TB (tau and left cauchy green) formulation
    
    %input:
    %-F: deformation gradient always as 3x3 matrix
    %-matData: struct containing all material constants and further info
    %-mapVoigt: voigt map. should be [1,5,9,4,8,7] (corrsponds to
    % [11,22,33,12,23,13])
    mapVoigt = mapVoigtObject.mapVoigt;
    dimension = mapVoigtObject.dimension;
    
    %Robin Pfefferkorn, 8.4.2021
   
    % determinant and cofactor
    b = F*F.';
    
    %unitiy tensors
    I = eye(3);
    
    %material data
    mu = materialObject.mu;
    kappa = materialObject.kappa;
    
    %settings spectral decomposition
    specDecomp = spectralDecomp(materialObject.specDecompAlgo,true);
    specDecomp.perturbThreshold = materialObject.perturbThreshold;
    
    %spectral decomposition
    specDecomp.decomposition(0.5*(b+b.'),mapVoigt);
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
    
    %stress tensor (kirchhoff)
    tauA = DW_logL + 2*DW_logc;
    tau = specDecomp.assemTens(tauA);
    
    % simple form of symmetric tangent
    eVol = kappa*ones(3,3);
    eIso = 2*mu*(I-1/3*ones(3,3));
    eTot = eIso+eVol;
    ccMatv = specDecomp.assemMatTang(tauA,eTot);

    tauv = tau(mapVoigt);
end
