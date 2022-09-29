function [W,S,Sv,CCv,errMat] = hyperelasticSC(materialObject,mapVoigtObject,C)
    errMat = 0; %no error
    %returns potential energy, stress and material derivatives for
    %hyperelastic materials in SC formulation
    
    %input:
    %-C: right Cauchy-Green
    mapVoigt = mapVoigtObject.mapVoigt;
    dimension = mapVoigtObject.dimension;
    
    mu = materialObject.mue;
    alpha = materialObject.alpha;
    beta = materialObject.beta;
    kappa = materialObject.kappa;
    
    %settings spectral decomposition
    specDecomp = spectralDecomp(materialObject.specDecompAlgo,false);
    specDecomp.perturbThreshold = materialObject.perturbThreshold;
    
    %spectral decomposition
    specDecomp.decomposition(C,mapVoigt);
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
    
    % stress tensor
    tauIsoA = sum(mu.*(lamBar.^alpha - 1/3*sum(lamBar.^alpha,1)) ,2);
    tauVol = (kappa/beta)*(1-J^(-beta));
    tauA = tauIsoA + tauVol;
    s = (lam.^-2).*tauA;
    S = specDecomp.assemTens(s);
    
    % simple form of symmetric tangent
    eVol = kappa*J^(-beta-1)*ones(3);
    eIso = zeros(3);
    for i=1:numel(alpha)
        eIso = eIso + mu(i)*alpha(i)*(diag(lamBar.^alpha(i))-1/3*(lamBar.^alpha(i)+(lamBar.').^alpha(i)) ...
            +1/9*sum(lamBar.^alpha(i)));
    end
    ETot = ((lam.^-2)*(lam.^-2).').*(eIso+eVol);
    CCv = specDecomp.assemMatTang(s,ETot);
       
    Sv = S(mapVoigt);
end