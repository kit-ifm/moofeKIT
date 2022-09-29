function [W,P,Pv,AAv,errMat] = hyperelasticPF(matierialObject,mapVoigtObject,F)
    errMat = 0; %no error
    %returns potential energy, stress and material derivatives for
    %hyperelastic materials in PF formulation

    mapVoigt = mapVoigtObject.mapVoigt;
    dimension = mapVoigtObject.dimension; 
    
    % cofactor
    H = 1/2*wedge(F,F);

    %material parameters
    mu = matierialObject.mu;
    lam = matierialObject.lambda;
    
    %additional kinematics
    trC = matierialObject.dotdot(F,F);
    
    %strain energy and derivatives
    W = 1/4*(lam/2+mu)*(trC-3)^2+mu*(trC-3)-mu/2*(matierialObject.dotdot(H,H)-3);
    DF_W = (lam/2+mu)*(trC-3)*F+2*mu*F;
    DH_W = -mu*H;
    DFF_W = ((lam/2+mu)*(trC-3)+2*mu)*eye(3^2) + 2*(lam/2+mu)*(F(mapVoigt)*F(mapVoigt).');
    
    %stress tensor (1.PK)
    P = DF_W + wedge(DH_W,F);
    
    %material derivative
    AAv = DFF_W - mu*matierialObject.specialMapCross(dimension,F)*matierialObject.specialMapCross(dimension,F) + matierialObject.specialMapCross(dimension,DH_W);
    
    %reduction in 2D (must always be computed in 3D due to
    %transformation)
    if dimension==2
        AAv = AAv([1,2,4,7],[1,2,4,7]);
    end
    
    Pv = P(mapVoigt);
end