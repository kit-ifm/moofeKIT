function [W,P,Pv,AAv,errMat] = hyperelasticPF(materialObject,mapVoigtObject,F)
    errMat = 0; %no error
    %returns potential energy, stress and material derivatives for
    %hyperelastic materials in PF formulation
    
    %input:
    %-F: deformation gradient    
    mapVoigt = mapVoigtObject.mapVoigt;
    dimension = mapVoigtObject.dimension;
    
    % determinant and cofactor
    H = 1/2*wedge(F,F);
    J = det(F);
    % voigt notation
    Hv = H(mapVoigt);    
    
    %material parameters
    mu = materialObject.mu;
    lam = materialObject.lambda;
    
    %strain energy and derivatives
    W = mu/2*(trace(F.'*F)-3) + lam/2*(log(J))^2 - mu*log(J);
    DF_W = mu*F;
    DJ_W = lam*log(J)*J^(-1) - mu/J;
    DFF_W = mu*eye(dimension^2);
    DJJ_W = lam*(J^(-2)-log(J)*J^(-2)) + mu*J^(-2);
    
    %stress tensors (1.PK)
    P = DF_W+DJ_W*H;
    
    %material derivative
    AAv = DFF_W + DJJ_W*(Hv*Hv.') + DJ_W*materialObject.specialMapCross(dimension,F);

    Pv = P(mapVoigt);
end