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
    c = det(b);
    
    %unitiy tensors
    I = eye(3);
    Iv = I(mapVoigt);
    if dimension == 2
        IIsymv = diag([1,1,0.5]);
    else
        IIsymv = diag([1,1,1,0.5,0.5,0.5]);
    end
    
    %material data
    aMR = materialObject.c1;
    bMR = materialObject.c2;
    cMR = materialObject.c;
    dMR = 2*(aMR+2*bMR);
    
    %additonal kinematics
    g = 1/2*wedge(b,b);
    trb = b(1,1)+b(2,2)+b(3,3);
    trg = g(1,1)+g(2,2)+g(3,3);
    
    %strain energy and derivatives
    W = aMR*(trb-3) + bMR*(trg-3) -dMR*log(sqrt(c)) + cMR/2*(sqrt(c)-1)^2;
    Dc_W = 1/2 * (cMR*(1-1/sqrt(c)) - dMR/c);
    Dcc_W = 1/2 * (cMR/2*c^(-3/2) + dMR*c^-2);
    
    %stress tensors (kirchhoff)
    tau = 2*(aMR*b + bMR*wedge(g,I) + c*Dc_W*I);
    
    %material derivative
    ccMatv = 4*(bMR*materialObject.mapVoigtCross(dimension,g) + c*((Dcc_W*c+Dc_W)*(Iv*Iv')-Dc_W*IIsymv));

    tauv = tau(mapVoigt);
end
