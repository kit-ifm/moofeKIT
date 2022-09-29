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
    
    %unitiy tensors
    I = eye(3);
    Iv = I(mapVoigt);
    if dimension == 2
        IIsymv = diag([1,1,0.5]);
    else
        IIsymv = diag([1,1,1,0.5,0.5,0.5]);
    end
    
    %material parameters
    mu = materialObject.mu;
    kappa = materialObject.kappa;
    
    %additonal kinematics
    trb = b(1,1)+b(2,2)+b(3,3);
    
    %strain energy and derivatives
    W = mu/2*(c^(-1/3)*trb-3) + kappa/2*(0.5*(c-1)-log(sqrt(c)));
    Db_W = mu/2*c^(-1/3)*I;
    Dc_W = kappa/4*(1-c^(-1)) - mu/6*c^(-4/3)*trb;
    
    %stress tensors (kirchhoff)
    tau = 2*(Db_W*b + c*Dc_W*I);
    
    % 2nd derivatives of strain energy
    Dbc_W = -mu/6*c^(-4/3)*I*b;%includes push forward
    Dcc_W = kappa/4*c^(-2) + 2*mu/9*c^(-7/3)*trb;
    
    %voigt notation of spatial derivatives
    DbcW_v = Dbc_W(mapVoigt);
    
    % elasticity tensor
    ccMatv = 4*(c*(DbcW_v*Iv'+Iv*DbcW_v.') + c*((Dcc_W*c+Dc_W)*(Iv*Iv') - Dc_W*IIsymv));
    
    tauv = tau(mapVoigt);
end
