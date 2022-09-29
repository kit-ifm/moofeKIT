function [W,tau,tauv,ccMatv,errMat] = hyperelasticTB(materialObject,mapVoigtObject,F)
    errMat = 0; %no error
    %returns potential energy, stress and material derivatives for
    %hyperelastic materials in TB (tau and left cauchy green) formulation
    
    %input:
    %-F: deformation gradient always as 3x3 matrix
    mapVoigt = mapVoigtObject.mapVoigt;
    dimension = mapVoigtObject.dimension;
    
    %Robin Pfefferkorn, 8.4.2021
   
    %unitiy tensor
    I = eye(3);

    %material parameters
    mu = materialObject.mu;
    lam = materialObject.lambda;
    
    %additional kinematics
    E = 1/2*(F.'*F-I);
    trE = E(1,1)+E(2,2)+E(3,3);
    
    %strain energy and derivatives
    W = lam/2*trE^2+mu*trace(E*E);
    
    %stress tensor (kirchhoff)
    tau = F*(lam*trE*I+2*mu*E)*F.';
    
    %material tangent tensor
    if dimension==2
        CCMatv = [lam+2*mu lam 0; lam lam+2*mu 0; 0 0 mu];
    else
        CCMatv = [...
            lam+2*mu lam lam 0 0 0;
            lam lam+2*mu lam 0 0 0;
            lam lam lam+2*mu 0 0 0;
            0 0 0 mu 0 0;
            0 0 0 0 mu 0;
            0 0 0 0 0 mu];
    end
    
    trafoMatF = materialObject.trafoMatCC(dimension,F);
    ccMatv = trafoMatF*CCMatv*trafoMatF.';

    tauv = tau(mapVoigt);
end
