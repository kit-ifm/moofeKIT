classdef ogdenClass < materialPackage.materialSuperClass
    properties 
        mue
        alpha
        beta
        kappa
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
        function set.mue(obj,inp)
            assert(isscalar(inp) && isnumeric(inp),'Input must be scalar numeric value');
            obj.mue = inp;
        end
        function set.alpha(obj,inp)
            assert(isscalar(inp) && isnumeric(inp),'Input must be scalar numeric value');
            obj.alpha = inp;
        end
        function set.beta(obj,inp)
            assert(isscalar(inp) && isnumeric(inp),'Input must be scalar numeric value');
            obj.beta = inp;
        end
        function set.kappa(obj,inp)
            assert(isscalar(inp) && isnumeric(inp),'Input must be scalar numeric value');
            obj.kappa = inp;
        end
    end    
end