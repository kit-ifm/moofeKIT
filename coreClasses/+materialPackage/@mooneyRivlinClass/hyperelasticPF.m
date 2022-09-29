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
    
    %material data
    aMR = materialObject.c1;
    bMR = materialObject.c2;
    cMR = materialObject.c;
    dMR = 2*(aMR+2*bMR);
    
    %strain energy and derivatives
    W = (aMR*(materialObject.dotdot(F,F)-3) + bMR*(materialObject.dotdot(H,H)-3) - dMR*log(J) + cMR/2*(J-1)^2);
    DF_W = 2*aMR*F;
    DH_W = 2*bMR*H;
    DJ_W = cMR*(J-1) - dMR/J;
    DJJ_W = cMR + dMR/J^2;
    
    %stress tensors (1.PK)
    P = DF_W+wedge(DH_W,F)+DJ_W*H;
    
    %material derivative (in 3D!)
    AAv = 2*aMR*eye(9) + 2*bMR*materialObject.specialMapCross(dimension,F)*materialObject.specialMapCross(dimension,F)...
        + DJJ_W*(H(mapVoigt)*H(mapVoigt).') + DJ_W*materialObject.specialMapCross(dimension,F) + materialObject.specialMapCross(dimension,DH_W);
    %reduction in 2D (must always be computed in 3D due to
    %transformation)
    if dimension==2
        AAv = AAv([1,2,4,7],[1,2,4,7]);
    end
    
    Pv = P(mapVoigt);
end