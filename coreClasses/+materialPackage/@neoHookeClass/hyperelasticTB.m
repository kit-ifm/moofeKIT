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
    
    % determinant and cofactor
    b = F*F.';
    c = det(b);
   
    %unitiy tensor
    I = eye(3);
    Iv = I(mapVoigt);
    if dimension==2
        IIsymv = diag([1,1,0.5]);
    else
        IIsymv = diag([1,1,1,0.5,0.5,0.5]);
    end
    
    %material parameters
    mu = materialObject.mu;
    lam = materialObject.lambda;
    
    %strain energy and derivatives
    W = mu/2*(b(1,1)+b(2,2)+b(3,3)-3) + lam/2*(1/2*log(c))^2 - mu/2*log(c);
    Db_W = mu/2*I;
    Dc_W = lam/4*log(c)*c^(-1) - mu/2*c^(-1);
    Dcc_W = lam/4*(c^(-2)-log(c)*c^(-2)) + mu/2*c^(-2);
    
    %stress tensors (kirchhoff)
    tau = 2*(Db_W*b + c*Dc_W*I);
    
    %material derivative
    ccMatv = 4*(c*((Dcc_W*c+Dc_W)*(Iv*Iv')-Dc_W*IIsymv));

    tauv = tau(mapVoigt);
end
