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
    
    %material data
    aMR = materialObject.c1;
    bMR = materialObject.c2;
    cMR = materialObject.c;
    dMR = 2*(aMR+2*bMR);
    
    %additonal kinematics
    trC = C(1,1)+C(2,2)+C(3,3);
    trG = G(1,1)+G(2,2)+G(3,3);
    
    %strain energy and derivatives
    W = aMR*(trC-3) + bMR*(trG-3) - dMR/2*log(c) + cMR/2*(sqrt(c)-1)^2;
    DC_W = aMR*I;
    DG_W = bMR*I;
    Dc_W = cMR/2*(1-1/sqrt(c)) - dMR/2*c^-1;
    Dcc_W = cMR/4*c^(-3/2) + dMR/2*c^-2;
    
    %stress tensors (1.PK)
    S = 2*(DC_W + wedge(DG_W,C) + Dc_W*G);
    
    %material derivative
    CCv = 4*(materialObject.specialMapCrossSymmetric(dimension,DG_W) + Dcc_W*(Gv*Gv.') + Dc_W*materialObject.specialMapCrossSymmetric(dimension,C));

    Sv = S(mapVoigt);
end