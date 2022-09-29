classdef saintVenantClass < materialPackage.materialSuperClass
    properties (Dependent = true)
        lambda  % first lame parameter
        mu      % second lame parameter
        G       % shear modulus
    end
    properties
        E = 0;   % Youngs modulus
        nu = 0;  % Poisson ratio
    end
    methods
        [Wpot,P,P_v,FF,errMat] = hyperelasticPF(obj,continuumObject,F);
        [Wpot,S,S_v,CC,errMat] = hyperelasticSC(obj,continuumObject,C);
        [Wpot,Tau,Tau_v,CC,errMat] = hyperelasticTB(obj,continuumObject,F);
        % TODO Umrechnung E, nu, lambda, mue
    end
    %% SET- und GET-Methoden
    methods
        function set.E(obj,inp)
            assert(isscalar(inp) && isnumeric(inp),'Input must be scalar numeric value');
            obj.E = inp;
        end
        function set.nu(obj,inp)
            assert(isscalar(inp) && isnumeric(inp),'Input must be scalar numeric value');
            obj.nu = inp;
        end
        function out = get.lambda(obj)
            out = getLame(obj);
        end
        function out = get.mu(obj)
            [~,out] = getLame(obj);
        end
        function out = get.G(obj)
            out = obj.mu;
        end
        function set.G(obj,inp)
            obj.mu = inp;
        end
        function materialDatabase(obj,type)
            [obj.E,obj.nu,obj.rho,~,~] = obj.database(type);
        end
    end
end