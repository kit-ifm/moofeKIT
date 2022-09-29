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
    Gv = G(mapVoigt);
    
    % unity tensor
    I = eye(3);
    
    %material parameters
    mu = materialObject.mu;
    lam = materialObject.lambda;
    
    %additonal kinematics
    trC = C(1,1)+C(2,2)+C(3,3);
    
    %strain energy and derivatives
    W = mu/2*(trC-3) + lam/8*(log(c))^2 - mu/2*log(c);
    DC_W = mu/2*I;
    Dc_W = (lam/4*log(c)-mu/2)*c^-1;
    Dcc_W = lam/4*c^-2 - (lam/4*log(c)-mu/2)*c^-2;
    
    %stress tensor (2.PK)
    S = 2*(DC_W + Dc_W*G);
    
    %material derivative
    CCv = 4*(Dcc_W*(Gv*Gv.') + Dc_W*materialObject.specialMapCrossSymmetric(dimension,C));
       
    Sv = S(mapVoigt);
end