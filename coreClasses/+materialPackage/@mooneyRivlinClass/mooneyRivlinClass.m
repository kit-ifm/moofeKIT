classdef mooneyRivlinClass < materialPackage.materialSuperClass
    properties 
        c1
        c2
        c
    end
    methods % FIXME mapVoigt -> Klasse -> abstrakt
        [Wpot,P,P_v,FF,errMat] = hyperelasticPF(obj,continuumObject,F);
        [Wpot,S,S_v,CC,errMat] = hyperelasticSC(obj,continuumObject,C);
        [Wpot,Tau,Tau_v,CC,errMat] = hyperelasticTB(obj,continuumObject,F);
    end
    %% SET- und GET-Methoden
    methods
        function set.c1(obj,inp)
            assert(isscalar(inp) && isnumeric(inp),'Input must be scalar numeric value');
            obj.c1 = inp;
        end
        function set.c2(obj,inp)
            assert(isscalar(inp) && isnumeric(inp),'Input must be scalar numeric value');
            obj.c2 = inp;
        end
        function set.c(obj,inp)
            assert(isscalar(inp) && isnumeric(inp),'Input must be scalar numeric value');
            obj.c = inp;
        end
    end    
end