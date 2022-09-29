function [W,P,Pv,AAv,errMat] = hyperelasticPF(materialObject,mapVoigtObject,F)
    errMat = 0; %no error
    %returns potential energy, stress and material derivatives for
    %hyperelastic materials in PF formulation
    
    %input:
    %-F: deformation gradient    
    mapVoigt = mapVoigtObject.mapVoigt;
    dimension = mapVoigtObject.dimension;
    
    %material data
    mu = materialObject.mu;
    kappa = materialObject.kappa;
    
    %kinematics
    Finv = F\eye(3);
    b = F*F.';%left cauchy green
    
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
    Tau = specDecomp.assemTens(tauA);
    
    % simple form of symmetric tangent
    eVol = kappa*ones(3,3);
    eIso = 2*mu*(eye(3)-1/3*ones(3,3));
    eTot = eIso+eVol;
    ccMatv = specDecomp.assemMatTang(tauA,eTot);
    
    %stress tensor (1.PK)
    P = Tau*Finv.';
    
    % transformation tangent
    trafoMat = materialObject.specialMapDot(dimension,Finv.');
    AAv = trafoMat*ccMatv*trafoMat.' + materialObject.specialMapDot(dimension,Finv*Tau*Finv.');
    
    Pv = P(mapVoigt);
end