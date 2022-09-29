classdef logStrainClass < materialPackage.materialSuperClass
    properties (Dependent = true)
        lambda
        mu
    end
    properties
        E = 0;
        nu = 0;
        kappa = 0;
        specDecompAlgo
        perturbThreshold
    end
    methods % FIXME mapVoigt -> Klasse -> abstrakt
        [Wpot,P,P_v,FF,errMat] = hyperelasticPF(obj,F,mapVoigt,dimension);
        [Wpot,S,S_v,CC,errMat] = hyperelasticSC(obj,C,mapVoigt,dimension);
        [Wpot,Tau,Tau_v,CC,errMat] = hyperelasticTB(obj,F,mapVoigt,dimension);
    end
    %% SET- und GET-Methoden
    methods
        function out = get.mu(obj)
            [~,out] = getLame(obj);
        end
    end    
end