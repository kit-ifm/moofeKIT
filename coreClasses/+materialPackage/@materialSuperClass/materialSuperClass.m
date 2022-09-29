classdef (Abstract) materialSuperClass < handle 
    properties
        name
        rho     % mass density
    end
    methods (Abstract = true)
        [Wpot,P,P_v,FF,errMat] = hyperelasticPF(obj,mapVoigtObject,F);        
        [Wpot,S,S_v,CC,errMat] = hyperelasticSC(obj,mapVoigtObject,C);
        [Wpot,Tau,Tau_v,CC,errMat] = hyperelasticTB(obj,mapVoigtObject,F);
    end
    methods (Static = true)
        [E,nu,rho,G,err] = database(type);
    end
    methods
        [Wpot,Salg,Salg_v,CCalg,errMat] = hyperelasticDiscreteGradient(obj,CN,CN1,dimension,mapVoigt,activeDisGra);
    end
    methods
        function [E,nu] = getParameter(obj)
            if ~isempty(obj.lambda) && ~isempty(obj.mu)
                E  = ((3*obj.lambda+2*obj.mu)*obj.mu)/(obj.lambda+obj.mu);
                nu = obj.lambda/(2*(obj.lambda+obj.mu));
            else
                E  = 0;
                nu = 0;
            end
        end
        function [lambda,mu] = getLame(obj)
            if ~isempty(obj.E) && ~isempty(obj.nu)
                lambda = obj.nu*obj.E/((1+obj.nu)*(1-2*obj.nu));
                mu     = obj.E/(2*(1+obj.nu));
            else
                lambda = 0;
                mu     = 0;
            end
        end
    end
    methods (Static = true)
        function out = specialMapCrossSymmetric(dimension,A)
            out = zeros(6,6);
            out(2,1) = A(3,3);
            out(3,1) = A(2,2);
            out(3,2) = A(1,1);
            out(4,3) = -A(2,1);
            out(5,1) = -A(3,2);
            out(6,4) = 0.5*A(3,2);
            out(5,4) = 0.5*A(3,1);
            out(6,2) = -A(3,1);
            out(6,5) = 0.5*A(2,1);
            out = (out + out.');
            out(4,4) = -0.5*A(3,3);
            out(5,5) = -0.5*A(1,1);
            out(6,6) = -0.5*A(2,2);
            if dimension == 2
                out = out([1,2,4],[1,2,4]);
            end
        end
        function out = specialMapCross(dimension,F)
            out = [0, F(3,3), F(2,2),      0,-F(3,2),      0,      0,-F(2,3),      0;
                   0,      0, F(1,1),      0,      0,-F(3,1),      0,      0,-F(1,3);
                   0,      0,      0,-F(2,1),      0,      0,-F(1,2),      0,      0;
                   0,      0,      0,      0, F(3,1),      0,-F(3,3),      0, F(2,3);
                   0,      0,      0,      0,      0,      0,      0,-F(1,1), F(1,2);
                   0,      0,      0,      0,      0,      0, F(3,2), F(2,1),-F(2,2);
                   0,      0,      0,      0,      0,      0,      0, F(1,3),      0;
                   0,      0,      0,      0,      0,      0,      0,      0,      0;
                   0,      0,      0,      0,      0,      0,      0,      0,      0];
            out = out+out.';
            if dimension==2
                out = out([1,2,4,7],[1,2,4,7]);
            end
        end
        function out = trafoMatCC(dimension,F)
            %returns transformation matrix to transform material tangent to
            %spatial one (in voigt notation)
            if dimension==2
                out = [...
                    F(1,1)^2, F(1,2)^2, 2*F(1,1)*F(1,2);
                    F(2,1)^2, F(2,2)^2, 2*F(2,1)*F(2,2);
                    F(1,1)*F(2,1), F(1,2)*F(2,2), F(1,1)*F(2,2)+F(1,2)*F(2,1)];
            else
                out = [...
                    F(1,1)^2, F(1,2)^2, F(1,3)^2, 2*F(1,1)*F(1,2), 2*F(1,2)*F(1,3), 2*F(1,1)*F(1,3);
                    F(2,1)^2, F(2,2)^2, F(2,3)^2, 2*F(2,1)*F(2,2), 2*F(2,2)*F(2,3), 2*F(2,1)*F(2,3);
                    F(3,1)^2, F(3,2)^2, F(3,3)^2, 2*F(3,1)*F(3,2), 2*F(3,2)*F(3,3), 2*F(3,1)*F(3,3);
                    F(1,1)*F(2,1), F(1,2)*F(2,2), F(1,3)*F(2,3), F(1,1)*F(2,2)+F(1,2)*F(2,1), F(1,2)*F(2,3)+F(1,3)*F(2,2), F(1,1)*F(2,3)+F(1,3)*F(2,1);
                    F(2,1)*F(3,1), F(2,2)*F(3,2), F(2,3)*F(3,3), F(2,1)*F(3,2)+F(2,2)*F(3,1), F(2,2)*F(3,3)+F(2,3)*F(3,2), F(2,1)*F(3,3)+F(2,3)*F(3,1);
                    F(1,1)*F(3,1), F(1,2)*F(3,2), F(1,3)*F(3,3), F(1,1)*F(3,2)+F(1,2)*F(3,1), F(1,2)*F(3,3)+F(1,3)*F(3,2), F(1,1)*F(3,3)+F(1,3)*F(3,1)];
            end
        end
        function out = mapVoigtCross(dimension,A)
            out = ...
                [0,A(3,3),A(2,2),           0,     -A(2,3),           0;
                0,     0,A(1,1),           0,           0,     -A(1,3);
                0,     0,     0,     -A(1,2),           0,           0;
                0,     0,     0,-0.25*A(3,3), 0.50*A(1,3), 0.50*A(2,3);
                0,     0,     0,           0,-0.25*A(1,1), 0.50*A(1,2);
                0,     0,     0,           0,           0,-0.25*A(2,2)];
            out = (out + out.');
            if dimension == 2
                out = out([1,2,4],[1,2,4]);
            end
        end
        function out = specialMapDot(dimension,F)
            out = [F(1,1),      0,      0, F(2,1),      0, F(3,1),      0,      0,      0;
                0, F(2,2),      0,      0, F(3,2),      0, F(1,2),      0,      0;
                0,      0, F(3,3),      0,      0,      0,      0, F(2,3), F(1,3);
                F(1,2),      0,      0, F(2,2),      0, F(3,2),      0,      0,      0;
                0, F(2,3),      0,      0, F(3,3),      0, F(1,3),      0,      0;
                F(1,3),      0,      0, F(2,3),      0, F(3,3),      0,      0,      0;
                0, F(2,1),      0,      0, F(3,1),      0, F(1,1),      0,      0;
                0,      0, F(3,2),      0,      0,      0,      0, F(2,2), F(1,2);
                0,      0, F(3,1),      0,      0,      0,      0, F(2,1), F(1,1)];
            if dimension==2
                out = out([1,2,4,7],[1,2,4,7]);
            end
        end
        function out = trace3d(b)
            out = b(1,1) + b(2,2) + b(3,3);
        end
        function out = dotdot(A,B)
            %double contraction (Holzapfel (2000) Page 14 (1.93))
            out = sum(sum(A.*B));
        end
    end        
end