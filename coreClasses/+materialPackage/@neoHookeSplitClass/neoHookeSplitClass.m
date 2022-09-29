classdef neoHookeSplitClass < materialPackage.materialSuperClass
    properties (Dependent = true)
        lambda
        mu
    end
    properties
        E = 0;
        nu = 0;
        kappa = 0;
    end
    methods % FIXME mapVoigt -> Klasse -> abstrakt
        [Wpot,P,P_v,FF,errMat] = hyperelasticPF(obj,continuumObject,F);
        [Wpot,S,S_v,CC,errMat] = hyperelasticSC(obj,continuumObject,C);
        [Wpot,Tau,Tau_v,CC,errMat] = hyperelasticTB(obj,continuumObject,F);
    end
    %% SET- und GET-Methoden
    methods
        function out = get.mu(obj)
            [~,out] = getLame(obj);
        end
    end    
end