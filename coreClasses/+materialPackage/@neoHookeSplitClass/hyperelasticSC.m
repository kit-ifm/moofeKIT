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
    c = det(C);
    
    % voigt notation
%     Cv = C(mapVoigt);
    Gv = G(mapVoigt);
    
    % unity tensor
    I = eye(3);
    Iv = I(mapVoigt);
    
    %material parameters
    mu = materialObject.mu;
    kappa = materialObject.kappa;
    
    %additional kinematics
    trC = trace(C);
    
    %strain energy and derivatives
    W = (mu/2*(c^(-1/3)*trC-3) + kappa/4*((c-1)-log(c)));
    DC_W = mu/2*c^(-1/3)*I;
    Dc_W = kappa/4*(1-c^(-1)) - mu/6*c^(-4/3)*trC;
    
    %stress tensors (2.PK)
    S = 2*(DC_W + Dc_W*G);
    
    % second derivatives of strain energy
    DCc_W_v = -mu/6*c^(-4/3)*Iv;
    DcC_W_v = -mu/6*c^(-4/3)*Iv;
    Dcc_W = kappa/4*c^(-2) + 2/9*mu*c^(-7/3)*trC;
    
    % elasticity tensor
    CCv = 4*(DCc_W_v*Gv.' + Gv*DcC_W_v.' + Dcc_W*(Gv*Gv.') + Dc_W*materialObject.specialMapCrossSymmetric(dimension,C));
       
    Sv = S(mapVoigt);
end