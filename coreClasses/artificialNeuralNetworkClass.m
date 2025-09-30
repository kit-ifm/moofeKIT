classdef artificialNeuralNetworkClass < handle
    %ARTIFICIALNEURALNETWORKCLASS ANN class for using ANN nets, e.g. for
    %incorporating tensorflow trained nets
    properties
        layer = struct('weight',[],'bias',[]);
        activation% computeActivationFunction('sigmoid'); % 'softPlus'
        alpha = 5e2;
        readData = struct('locationFolderANN','','biasNameCell','','weightNameCell','');
        normalizationFactor = 1;
        A = 0;
    end

    methods
        function obj = artificialNeuralNetworkClass()
            % constructor
        end
        function obj = read(obj,typeData)
            if strcmpi(typeData,'bias')
                typeDataCell = obj.readData.biasNameCell;
            elseif strcmpi(typeData,'weight')
                typeDataCell = obj.readData.weightNameCell;
            else
                error('data type not available')
            end
            for ii = 1:numel(typeDataCell)
                obj.layer(ii).(typeData) = readmatrix(strcat(obj.readData.locationFolderANN,typeDataCell{ii},'.txt'),'Delimiter',' ','ConsecutiveDelimitersRule','join');
                obj.layer(ii).(typeData) = readmatrix(strcat(obj.readData.locationFolderANN,typeDataCell{ii},'.txt'),'Delimiter',' ','ConsecutiveDelimitersRule','join');
            end
        end
        function obj = computeActivationFunction(obj,nameOfFunction)
            if strcmpi(nameOfFunction,'softPlus')
                obj.activation.function = @(x) log(1+exp(x));
                obj.activation.derivative = @(x) exp(x)./(1+exp(x));
                obj.activation.secondDerivative = @(x) exp(x)./((1+exp(x)).^2);
            elseif strcmpi(nameOfFunction,'sigmoid')
                obj.activation.function = @(x) 1./(1+exp(-x));
                obj.activation.derivative = @(x) x.*(1 - x);
            end
        end
        function W = computeEnergyANN(obj,C,G,c)
            w1 = obj.layer(1).weight;
            w2 = obj.layer(2).weight;
            b1 = obj.layer(1).bias;
            b2 = obj.layer(2).bias;
            I = trace(C);
            II = trace(G);
            J = sqrt(c);
            Jstar = -sqrt(c);
            x = [I;II;J;Jstar];
            h1 = w1*x + b1;
            h2 = obj.activation.function(h1);
%             Wdev = w2*h2 + b2;
            Wdev = w2*h2;
%             Wvol = obj.alpha*(J-1)^2;
            Wvol = (J + 1/J - 2 )^2;
            Wstress = -obj.A*(J-1);
            Wenergy = b2;
            W = obj.normalizationFactor*(Wdev + Wstress + Wenergy) + Wvol;
        end
        function W = computeEnergyANNCGJ(obj,C,G,J)
            w1 = obj.layer(1).weight;
            w2 = obj.layer(2).weight;
            b1 = obj.layer(1).bias;
            b2 = obj.layer(2).bias;
            I = trace(C);
            II = trace(G);
            Jstar = -J;
            x = [I;II;J;Jstar];
            h1 = w1*x + b1;
            h2 = obj.activation.function(h1);
%             Wdev = w2*h2 + b2;
            Wdev = w2*h2;
%             Wvol = obj.alpha*(J-1)^2;
            Wvol = (J + 1/J - 2 )^2;
            Wstress = -obj.A*(J-1);
            Wenergy = b2;
            W = obj.normalizationFactor*(Wdev + Wstress + Wenergy) + Wvol;
            %% TODO: for Anne's model
            % W = obj.normalizationFactor*(Wdev + Wstress + Wenergy + Wvol);
        end        
        function [dW_I, dW_II, dW_J, dW_Jstar] = computeDiffEnergyANN(obj,C,G,c)
            w1 = obj.layer(1).weight;
            w2 = obj.layer(2).weight;
            b1 = obj.layer(1).bias;
            I = trace(C);
            II = trace(G);
            J = sqrt(c);
            Jstar = -J;
            x = [I;II;J;Jstar];
            h1 = w1*x + b1;
            dh2dh1 = diag(obj.activation.derivative(h1));
            dWvol_J = 2*(J + 1 / J - 2)*(1-J^(-2)); %=2*(J-J^-3-2+2J^-2)
            dWvol_Jstar = 0;
            dW_x = obj.normalizationFactor*w2*(dh2dh1*w1);
            dW_x(3) = dW_x(3) + dWvol_J + -obj.A*obj.normalizationFactor;
            dW_x(4) = dW_x(4) + dWvol_Jstar;
            dW_I = dW_x(1);
            dW_II = dW_x(2);
            dW_J = dW_x(3);
            dW_Jstar = dW_x(4);
        end
        function [dW_I, dW_II, dW_J] = computeDiffEnergyANNJ(obj,C,G,J)
            w1 = obj.layer(1).weight;
            w2 = obj.layer(2).weight;
            b1 = obj.layer(1).bias;
            I = trace(C);
            II = trace(G);
            Jstar = -J;
            x = [I;II;J;Jstar];
            h1 = w1*x + b1;
            dh2dh1 = diag(obj.activation.derivative(h1));
            dWvol_J = 2*(J + 1 / J - 2)*(1-J^(-2)); %=2*(J-J^-3-2+2J^-2)
            dWvol_Jstar = 0;
            dW_x = obj.normalizationFactor*w2*(dh2dh1*w1);
            dW_x(3) = dW_x(3) + dWvol_J + -obj.A*obj.normalizationFactor;
            dW_x(4) = dW_x(4) + dWvol_Jstar;
            dW_I = dW_x(1);
            dW_II = dW_x(2);
            dW_J = dW_x(3);
            dW_Jstar = dW_x(4);
            dW_J = dW_J + dW_Jstar*(-1);
        end
        function [dW_I, dW_II, dW_J] = computeDiffEnergyANNIIIJ(obj,I,II,J)
            w1 = obj.layer(1).weight;
            w2 = obj.layer(2).weight;
            b1 = obj.layer(1).bias;
            Jstar = -J;
            x = [I;II;J;Jstar];
            h1 = w1*x + b1;
            dh2dh1 = diag(obj.activation.derivative(h1));
            dWvol_J = 2*(J + 1 / J - 2)*(1-J^(-2)); %=2*(J-J^-3-2+2J^-2)
            dWvol_Jstar = 0;
            dW_x = obj.normalizationFactor*w2*(dh2dh1*w1);
            dW_x(3) = dW_x(3) + dWvol_J + -obj.A*obj.normalizationFactor;
            dW_x(4) = dW_x(4) + dWvol_Jstar;
            dW_I = dW_x(1);
            dW_II = dW_x(2);
            dW_J = dW_x(3);
            dW_Jstar = dW_x(4);
            dW_J = dW_J + dW_Jstar*(-1);
        end
        function [dW_I, dW_II, dW_III] = computeDiffEnergyInvariantANN(obj,I,II,III)
            w1 = obj.layer(1).weight;
            w2 = obj.layer(2).weight;
            b1 = obj.layer(1).bias;
            J = sqrt(III);
            Jstar = -J;
            x = [I;II;J;Jstar];
            h1 = w1*x + b1;
            dh2dh1 = diag(obj.activation.derivative(h1));
            dWvol_J = 2*(J + 1 / J - 2)*(1-J^(-2)); %=2*(J-J^-3-2+2J^-2)
            dWvol_Jstar = 0;
            dW_x = obj.normalizationFactor*w2*(dh2dh1*w1);
            dW_x(3) = dW_x(3) + dWvol_J + -obj.A*obj.normalizationFactor;
            dW_x(4) = dW_x(4) + dWvol_Jstar;
            dW_I = dW_x(1);
            dW_II = dW_x(2);
            dW_J = dW_x(3);
            dW_Jstar = dW_x(4);
            dJ_c = 1/2*(III)^(-1/2);
            dJstar_c = -dJ_c;
            dW_III = dW_J*dJ_c + dW_Jstar*dJstar_c;            
        end        
        function [d2W_I_I] = compute2ndDiffEnergyANN(obj,C,G,c)
            w1 = obj.layer(1).weight;
            w2 = obj.layer(2).weight;
            b1 = obj.layer(1).bias;
            I = trace(C);
            II = trace(G);
            J = sqrt(c);
            Jstar = -J;
            x = [I;II;J;Jstar];
            h1 = w1*x + b1;
            d2W_I_I = obj.normalizationFactor*w1'*diag(w2'.*obj.activation.secondDerivative(h1))*w1;
%             d2Wvol_J_J = 2*obj.alpha;            
            d2Wvol_J_J = 2*(1+3*J^(-4)-4*J^-3);
            %     d2Wvol_J_Jstar = 0;
            %     d2Wvol_Jstar_J = 0;
            %     d2Wvol_Jstar_Jstar = 0;
            d2W_I_I(3,3) = d2W_I_I(3,3) + d2Wvol_J_J;
            %     d2W_I_I(3,4) = d2W_I_I(3,4) + d2Wvol_J_Jstar;
            %     d2W_I_I(4,3) = d2W_I_I(4,3) + d2Wvol_Jstar_J;
            %     d2W_I_I(4,4) = d2W_I_I(4,4) + d2Wvol_Jstar_Jstar;
        end        
        function [d2W_I_I] = compute2ndDiffEnergyANNCGJ(obj,C,G,J)
            w1 = obj.layer(1).weight;
            w2 = obj.layer(2).weight;
            b1 = obj.layer(1).bias;
            I = trace(C);
            II = trace(G);
            Jstar = -J;
            x = [I;II;J;Jstar];
            h1 = w1*x + b1;
            d2W_I_I = obj.normalizationFactor*w1'*diag(w2'.*obj.activation.secondDerivative(h1))*w1;
%             d2Wvol_J_J = 2*obj.alpha;            
            d2Wvol_J_J = 2*(1+3*J^(-4)-4*J^-3);
            %     d2Wvol_J_Jstar = 0;
            %     d2Wvol_Jstar_J = 0;
            %     d2Wvol_Jstar_Jstar = 0;
            d2W_I_I(3,3) = d2W_I_I(3,3) + d2Wvol_J_J;
            %     d2W_I_I(3,4) = d2W_I_I(3,4) + d2Wvol_J_Jstar;
            %     d2W_I_I(4,3) = d2W_I_I(4,3) + d2Wvol_Jstar_J;
            %     d2W_I_I(4,4) = d2W_I_I(4,4) + d2Wvol_Jstar_Jstar;
        end        
        function [dW_C, dW_G, dW_c, dW_I, dW_II, dW_J, dW_Jstar, dJ_c, dJstar_c] = computeDiffEnergyCGcANN(obj,C,G,c)
            [dW_I, dW_II, dW_J, dW_Jstar] = computeDiffEnergyANN(obj,C,G,c);
            dJ_c = 1/2*(c)^(-1/2);
            dJstar_c = -dJ_c;
            dW_c = dW_J*dJ_c + dW_Jstar*dJstar_c;
            dIC_C = eye(3);
            dIIC_G = eye(3);
            dW_C = dW_I*dIC_C;
            dW_G = dW_II*dIIC_G;
        end
        function [dW_I, dW_II, dW_c, dW_J, dW_Jstar, dJ_c, dJstar_c] = computeDiffEnergyIIIcANN(obj,C,G,c)
            [dW_I, dW_II, dW_J, dW_Jstar] = computeDiffEnergyANN(obj,C,G,c);
            dJ_c = 1/2*(c)^(-1/2);
            dJstar_c = -dJ_c;
            dW_c = dW_J*dJ_c + dW_Jstar*dJstar_c;
        end
        function [dWvec] = computeDiffEnergyANNVec(obj,C,G,c)
            [dW_I, dW_II, dW_J, dW_Jstar] = computeDiffEnergyANN(obj,C,G,c);
            dWvec = [dW_I, dW_II, dW_J, dW_Jstar];
        end
        function [dWvec] = computeDiffEnergyANNVecc(obj,C,G,c)
            [dW_I, dW_II, dW_c] = computeDiffEnergyIIIcANN(obj,C,G,c);
            dWvec = [dW_I, dW_II, dW_c];
        end
    end
end
