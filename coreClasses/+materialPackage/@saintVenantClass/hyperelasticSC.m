function [W,S,Sv,CCv,errMat] = hyperelasticSC(materialObject,mapVoigtObject,C)
    errMat = 0; %no error
    %returns potential energy, stress and material derivatives for
    %hyperelastic materials in SC formulation
    
    %input:
    %-C: right Cauchy-Green
    mapVoigt = mapVoigtObject.mapVoigt;
    dimension = mapVoigtObject.dimension;
    
    % determinant and cofactor
    G = 1/2*wedge(C,C);
    
    % unity tensor
    I = eye(3);
    Iv = I(mapVoigt);
    
    %material parameters
    mu = materialObject.mu;
    lam = materialObject.lambda;
    
    %additional kinematics
    trC = trace(C);
    
    %strain energy and derivatives
    W = 1/4*(lam/2+mu)*(trC-3)^2+mu*(trC-3)-mu/2*(trace(G)-3);
    DC_W = (1/2*(lam/2+mu)*(trC-3) + mu)*I;
    DG_W = -mu/2*I;
    DCC_W = 1/2*(lam/2+mu)*(Iv*Iv.');
    
    %stress tensor (2.PK)
    S = 2*(DC_W + wedge(DG_W,C));
    
    %material derivative
    CCv = 4*(DCC_W + materialObject.specialMapCrossSymmetric(dimension,DG_W));
       
    Sv = S(mapVoigt);
end