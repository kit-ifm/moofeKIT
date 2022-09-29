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
    kappa = materialObject.kappa;
    
    %strain energy and derivatives
    W = (mu/2*(J^(-2/3)*materialObject.dotdot(F,F)-3) + kappa/2*(1/2*(J^2-1)-log(J)));
    DF_W = mu*J^(-2/3)*F;
    DJ_W = kappa/2*(J-J^-1) - mu/3*J^(-5/3)*materialObject.dotdot(F,F);
    
    %stress tensors (1.PK)
    P = DF_W+DJ_W*H;
    
    % second derivatives of strain energy
    DFJ_W = -2/3*mu*J^(-5/3)*F;
    DJJ_W = kappa/2*(1+J^(-2)) + 5*mu/9*J^(-8/3)*materialObject.dotdot(F,F);
    DJF_W = -2/3*mu*J^(-5/3)*F;
    
    % voigt maps
    DFJ_W_v = DFJ_W(mapVoigt);
    DJF_W_v = DJF_W(mapVoigt);
    
    % elasticity tensor
    AAv = mu*J^(-2/3)*eye(dimension^2) + DJJ_W*(Hv*Hv.') + DJ_W*materialObject.specialMapCross(dimension,F) + DFJ_W_v*Hv.'+Hv*DJF_W_v.';
    
    Pv = P(mapVoigt);
end