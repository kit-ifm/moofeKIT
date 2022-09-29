function [W,P,Pv,AAv,errMat] = hyperelasticPF(materialObject,mapVoigtObject,F)
    errMat = 0; %no error
    %returns potential energy, stress and material derivatives for
    %hyperelastic materials in PF formulation
    mapVoigt = mapVoigtObject.mapVoigt;
    dimension = mapVoigtObject.dimension;
    
    %material data
    mu = materialObject.mue;
    alpha = materialObject.alpha;
    beta = materialObject.beta;
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
    if any(lam>1e4) || any(lam<1e-4)
        errMat = 1;
    end
    lamBar = J^(-1/3)*lam;
    
    %strain energy
    Wvol = kappa*beta^(-2)*(beta*log(J)+J^(-beta)-1);
    Wiso = sum(sum( (mu./alpha).*(lamBar.^alpha-1) ));
    W = Wvol + Wiso;
    
    % stress tensor (kirchhoff)
    tauIsoA = sum(mu.*(lamBar.^alpha - 1/3*sum(lamBar.^alpha,1)) ,2);
    tauVol = (kappa/beta)*(1-J^(-beta));
    tauA = tauIsoA + tauVol;
    Tau = specDecomp.assemTens(tauA);
    %stress tensor (1.PK)
    P = Tau*Finv.';
    
    % simple form of symmetric tangent
    eVol = kappa*J^(-beta-1)*ones(3);
    eIso = zeros(3);
    for i=1:numel(alpha)
        eIso = eIso + mu(i)*alpha(i)*(diag(lamBar.^alpha(i))-1/3*(lamBar.^alpha(i)+(lamBar.').^alpha(i)) ...
            +1/9*sum(lamBar.^alpha(i)));
    end
    eTot = eIso+eVol;
    ccMatv = specDecomp.assemMatTang(tauA,eTot);
    
    % transformation to PF framework
    trafoMat = materialObject.specialMapDot(dimension,Finv.');
    AAv = trafoMat*ccMatv*trafoMat.' + materialObject.specialMapDot(dimension,Finv*Tau*Finv.');

    Pv = P(mapVoigt);
end