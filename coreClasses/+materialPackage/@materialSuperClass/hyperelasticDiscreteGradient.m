function [WN1,Salg,Salgv,CCalg,errMat] = hyperelasticDiscreteGradient(materialObject,mapVoigtObject,CN,CN1,activeDisGra)
    %computes the algorithmic stresses and corresponding tangent
    %Robin Pfefferkorn 13.4.2021
    %
    %input:
    %-CN, CN1: right Cauchy-Green at time tN and tN+1
    mapVoigt = mapVoigtObject.mapVoigt;
    dimension = mapVoigtObject.dimension;
    
    %evaluation at mid point of strains
    CMid = 0.5*(CN+CN1);
    [~,SMid,SMidv,CCMid,~] = hyperelasticSC(materialObject,mapVoigtObject,CMid);
    
    %evaluation at endpoint of strains (always needed for Epot)
    [WN1,~,SN1v,~,errMat] = hyperelasticSC(materialObject,mapVoigtObject,CN1);
    
    if activeDisGra
        %potential energy at time n
        WN = hyperelasticSC(materialObject,mapVoigtObject,CN);
        
        %matrix M
        M = CN1-CN;
        Mv = M(mapVoigt);
        Msymv = Mv;
        Msymv(dimension+1:end) = Msymv(dimension+1:end)*2;
        MddM = materialObject.dotdot(M,M);
        
        %projection factor
        factDG = (2*(WN1-WN)-materialObject.dotdot(SMid,M))/MddM;
        
        %identity tensor
        if dimension==2
            IIsymv = diag([1,1,0.5]);
        else
            IIsymv = diag([1,1,1,0.5,0.5,0.5]);
        end
        
        Salg = SMid + factDG*M;
        Salgv = Salg(mapVoigt);
        CCalg = CCMid + 2*( factDG*(2*IIsymv-4*(Mv*Mv.')/MddM) + Mv/MddM*(2*(SN1v-SMidv)-1/2*CCMid*Msymv).');
    else
        %mixed integrator or small difference between CN1 and CN
        Salg = SMid;
        Salgv = SMidv;
        CCalg = CCMid;
    end
end
