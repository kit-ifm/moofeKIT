function plot(obj,varargin)
% Plot Object

% set function-handle for determination of the color data
for i = 1:numel(obj)
    obj(i).PLOT_DATA.detColorDataFcn = @detColorData;
end

plotObject(obj,varargin) %common/post/plotObject

end

%% Function to determine color data
function colorData = detColorData(obj,options)
%obj: your Object
%options: structure with all plot-options (see common/plot/ploObject/checkInputArguments.m)
colorData = [];
if options.temp
    switch options.time
        case 'REF'
            colorData =  obj.QN1(:,obj.DIM+1);
        case 'N'
            colorData =  obj.QN(:,obj.DIM+1);
        case 'N1'
            if strcmp(obj.DESIGNATOR,'generic')
                if strcmp(obj.TAU,'entr')
                    colorData = obj.PLOTData;
                elseif strcmp(obj.TAU,'engy')
                    colorData = obj.PLOTData;
                elseif strcmp(obj.TAU,'temp')
                    colorData =  obj.QN1(:,obj.DIM+1);
                end
            else
                colorData =  obj.QN1(:,obj.DIM+1);
            end
    end
elseif options.inv1
    DESIGNATOR = obj.DESIGNATOR;
    if strcmp(DESIGNATOR,'generic')
        fdof = obj.NODESDOF(:,obj.DIM+obj.ADDFIELDS);
        maxElAct = max(fdof);
        internal=residual(obj,'J');
        rHam = full(sparse(vertcat(internal.structure(:).edofE),1,vertcat(internal.structure(:).Re),maxElAct,1));
        kHam = sparse(vertcat(internal.structure(:).pI), vertcat(internal.structure(:).pJ), vertcat(internal.structure(:).pK),maxElAct,maxElAct);
        colorData = zeros(size(fdof,1),1);
        colorData(:) = kHam(fdof,fdof)\rHam(fdof,1);
        [max(colorData) min(colorData)]
    else
        error('..')
    end
elseif options.stress~=0
    DESIGNATOR = obj.DESIGNATOR;
    DIM = obj.DIM;
    
    %TODO: implement all different in the same way as 'genericCascadeSC'
    if strcmp(DESIGNATOR,'genericCascadeSC') || strcmp(DESIGNATOR,'generic')
        fdof = obj.NODESDOF(:,obj.DIM+obj.ADDFIELDS);
        maxElAct = max(fdof);
        internal=residual(obj,'stress',options.stress);
        rHam = full(sparse(vertcat(internal.structure(:).edofE),1,vertcat(internal.structure(:).Re),maxElAct,1));
        kHam = sparse(vertcat(internal.structure(:).pI), vertcat(internal.structure(:).pJ), vertcat(internal.structure(:).pK),maxElAct,maxElAct);
        colorData = zeros(size(fdof,1),1);
        colorData(:) = kHam(fdof,fdof)\rHam(fdof,1);
    else
        ansatzfunction = obj.Ansatzfunction;
        EDOF = obj.EDOF;
        edof = obj.GLOBALEDOF;
        matName = obj.MAT.name;
        NAll = obj.SHAPEF.N;
        dNrAll = obj.SHAPEF.dNr;
        NGP = obj.NGP;
        I = eye(DIM);
        zw = struct('EDOF',[],'Spann',[],'pI',[],'pJ',[],'pK',[]);
        wp = obj.SHAPEF.wp;
        QREF = obj.QREF;
        QN1 = obj.QN1;
        KAPPA0 = [];
        BETA = [];
        GAMMA = [];
        ALPHA = [];
        ALPHA0 = [];
        THETA0 = [];
        mue = [];
        mu = [];
        lambda = [];
        beta = [];
        theta0 = [];
        c1 = [];
        c2 = [];
        c = [];
        alphaMat = [];
        wkalpha = [];
        numElements = numel(EDOF);
        switch lower(matName)
            case 'ogdenthermo'
                KAPPA0 = obj.MAT.kappa0;
                BETA = obj.MAT.beta;
                GAMMA = obj.MAT.gamma;
                ALPHA = obj.MAT.alpha;
                ALPHA0 = obj.MAT.alpha0;
                THETA0 = obj.MAT.theta0;
                mue = obj.MAT.mue;
            case 'neohookethermo'
                mu = obj.MAT.mu;
                lambda = obj.MAT.lambda;
                beta = obj.MAT.beta;
                theta0 = obj.MAT.theta0;
                % MATERIAL = []; %added by rp (parfor doesn't work without)
            case 'mooneyrivlinthermo'
                c1 = obj.MAT.c1;
                c2 = obj.MAT.c2;
                c = obj.MAT.c;
                beta = obj.MAT.beta;
                theta0 = obj.MAT.theta0;
            case 'linearthermo'
                lambda = obj.MAT.lambda;
                mu = obj.MAT.G;
                alphaMat = obj.MAT.alpha;
                wkalpha = obj.MAT.wkalpha;
                theta0 = obj.MAT.theta0;
            otherwise
                error('Material model is currently not implemented.')
        end
        % Reshape data & stress for hr element
        shape=[1 7 9 4 2 8 6 5 3]';
        Nlin = ones(NGP,0);
        SIGMA_FN1 = zeros(numElements,0);
        SIGMA_HN1 = zeros(numElements,0);
        SIGMA_JN1 = zeros(numElements,0);
        if strcmpi(DESIGNATOR,'mmhrSplit')
            % Linear shape functions
            Nlin = obj.SHAPEF_LM.N;
            % Conjungated stress
            SIGMA_FN1 = obj.SIGMA_FN1;
            SIGMA_HN1 = obj.SIGMA_HN1;
            SIGMA_JN1 = obj.SIGMA_JN1;
        end
        parfor j = 1:numElements
            dNr = [];
            N = [];
            [edofH1,edofH2] = expandEdof(EDOF{j});
            zw(j).pI = edofH1';
            zw(j).pJ = edofH2';
            zw(j).EDOF = EDOF{j}';
            edofLocal = EDOF{j};
            numDOFs = numel(edof{j});
            if strcmpi(ansatzfunction,'Lagrange')
                N = NAll;
                dNr = dNrAll;
            elseif strcmpi(ansatzfunction,'HNURBS')
                N = NAll{j};
                dNr = dNrAll{j};
            end
            edRef = QREF(edofLocal,1:DIM)';
            edN1 = QN1(edofLocal,1:DIM)';
            n1Temp = QN1(edofLocal,DIM+1)';
            J = QREF(edofLocal,1:DIM)'*dNr';
            J1 = QN1(edofLocal,1:DIM)'*dNr';
            nen = size(N,2);
            switch lower(DESIGNATOR)
                case 'wedgesplit'
                    switch lower(matName)
                        case 'mooneyrivlinthermo'
                            Spann = zeros(size(N,2),1);
                            Me = zeros(nen);
                            for k = 1:NGP
                                tempSpann = 0; %clear for each iteration
                                indx = DIM*k-(DIM-1):DIM*k;
                                detJ = det(J(:,indx)');
                                detJ1 = det(J1(:,indx)');
                                if detJ < 10*eps
                                    error('Jacobi determinant equal or less than zero.')
                                end
                                dNx = (J(:,indx)')\dNr(indx,:);
                                % Deformation gradient, co-factor & Jacobian determinant
                                FN1 = edN1*dNx';
                                HN1 = 0.5*wedge(FN1,FN1);
                                JN1 = det(FN1);
                                % Temperature
                                thetaN1 = N(k,:)*n1Temp';
                                % First derivatives
                                DPsi_FN1 = 2*c1*JN1^(-2/3)*FN1;
                                DPsi_HN1 = 3*c2*JN1^(-2)*(sum(sum(HN1.*HN1)))^(1/2)*HN1;
                                DPsi_JN1 = -2/3*c1*JN1^(-5/3)*sum(sum(FN1.*FN1)) - 2*c2*JN1^(-3)*(sum(sum(HN1.*HN1)))^(3/2) + c*(JN1-1) - DIM*beta*(thetaN1-theta0)*c;
                                % First Piola Kirchhoff stress tensor
                                PN1 = DPsi_FN1 + wedge(DPsi_HN1,FN1) + DPsi_JN1*HN1;
                                % Cauchy stress
                                sigma = 1/JN1*PN1*FN1';
                                tempSpann = selectStress(sigma,options.stress,DIM);
                                Spann   = Spann + N(k,:)'*tempSpann*detJ1*wp(k);
                                Me = Me + (N(k,:)'*N(k,:))*detJ1*wp(k);
                            end
                            zw(j).Spann = Spann;
                            zw(j).pK = Me(:);
                        otherwise
                            error('Material model is currently not implemented.')
                    end
                case 'mmhrsplit'
                    switch lower(matName)
                        case 'mooneyrivlinthermo'
                            Spann = zeros(nen,1);
                            Me = zeros(nen);
                            for k = 1:NGP
                                tempSpann = 0; %clear for each iteration
                                indx = DIM*k-(DIM-1):DIM*k;
                                detJ = det(J(:,indx)');
                                detJ1 = det(J1(:,indx)');
                                if detJ < 10*eps
                                    error('Jacobi determinant equal or less than zero.')
                                end
                                dNx = (J(:,indx)')\dNr(indx,:);
                                % Deformation gradient, co-factor & Jacobian determinant
                                FN1 = edN1*dNx';
                                HN1 = 0.5*wedge(FN1,FN1);
                                JN1 = det(FN1);
                                % Conjungated stress
                                SIGMA_FN1e = zeros(3);
                                SIGMA_HN1e = zeros(3);
                                for m = 1:size(Nlin,2)
                                    SIGMA_FN1e= SIGMA_FN1e + Nlin(k,m)*reshape(SIGMA_FN1(j,ones(1,9)*(m-1)*9+shape'),3,3);
                                    SIGMA_HN1e= SIGMA_HN1e + Nlin(k,m)*reshape(SIGMA_HN1(j,ones(1,9)*(m-1)*9+shape'),3,3);
                                end
                                SIGMA_JN1e  = SIGMA_JN1(j,:);
                                % 1. PK, cauchy stress & van Mises stress
                                PN1 = JN1^(-1/3)*SIGMA_FN1e + JN1^(-2/3)*wedge(SIGMA_HN1e,FN1) + (SIGMA_JN1e-1/3*JN1^(-4/3)*sum(sum(FN1.*SIGMA_FN1e))-2/3*JN1^(-5/3)*sum(sum(HN1.*SIGMA_HN1e)))*HN1;
                                sigma = 1/JN1*PN1*FN1';
                                tempSpann = selectStress(sigma,options.stress,DIM);
                                Spann   = Spann + N(k,:)'*tempSpann*detJ1*wp(k);
                                Me = Me + (N(k,:)'*N(k,:))*detJ1*wp(k);
                            end
                            zw(j).EDOF = edofLocal(:);
                            zw(j).Spann = Spann;
                            zw(j).pK = Me(:);
                        otherwise
                            error('Material model is currently not implemented.')
                    end
                case 'disp'
                    switch lower(matName)
                        case 'ogdenthermo'
                            Spann = zeros(size(N,2),1);
                            Me = zeros(nen);
                            for k = 1:NGP
                                tempSpann = 0; %clear for each iteration
                                indx = DIM*k-(DIM-1):DIM*k;
                                detJ = det(J(:,indx)');
                                detJ1 = det(J1(:,indx)');
                                if detJ < 10*eps
                                    error('Jacobi determinant equal or less than zero.')
                                end
                                dNx = (J(:,indx)')\dNr(indx,:);
                                FN1 = edN1*dNx';
                                detF = det(FN1);
                                % Cauchy-Green tensor
                                CN1 = FN1'*FN1;
                                % eigenvalue decomposition
                                [NeigN1,VeigN1] = eig(CN1);
                                lambda1N1 = sqrt(VeigN1(1,1));
                                lambda2N1 = sqrt(VeigN1(2,2));
                                lambda3N1 = sqrt(VeigN1(3,3));
                                tol = 10^(-14)*max([norm(lambda1N1) norm(lambda2N1) norm(lambda3N1)]);
                                pertubed = 10^(-9);
                                if norm(lambda1N1-lambda2N1) < tol
                                    lambda1N1 = (1+pertubed)*lambda1N1;
                                    lambda2N1 = (1-pertubed)*lambda2N1;
                                    lambda3N1 = 1/((1+pertubed)*(1-pertubed))*lambda3N1;
                                elseif norm(lambda1N1-lambda3N1) < tol
                                    lambda1N1 = (1+pertubed)*lambda1N1;
                                    lambda3N1 = (1-pertubed)*lambda3N1;
                                    lambda2N1 = 1/((1+pertubed)*(1-pertubed))*lambda2N1;
                                elseif norm(lambda2N1-lambda3N1) < tol
                                    lambda2N1 = (1+pertubed)*lambda2N1;
                                    lambda3N1 = (1-pertubed)*lambda3N1;
                                    lambda1N1 = 1/((1+pertubed)*(1-pertubed))*lambda1N1;
                                end
                                JN1 = lambda1N1*lambda2N1*lambda3N1;
                                lambda1N1Tilde = JN1^(-1/3)*lambda1N1;
                                lambda2N1Tilde = JN1^(-1/3)*lambda2N1;
                                lambda3N1Tilde = JN1^(-1/3)*lambda3N1;
                                % 2.PK
                                pN1 = KAPPA0*(N(k,:)*n1Temp'/THETA0)*BETA^(-1)*(1/JN1-JN1^(-BETA-1))-...
                                    3*ALPHA0*KAPPA0*JN1^(GAMMA-1)*(N(k,:)*n1Temp'-THETA0);
                                SVolN1 = pN1*JN1*I/CN1;
                                S1N1 = 1/lambda1N1^2*(sum((mue*(N(k,:)*n1Temp'/THETA0).*lambda1N1Tilde.^(ALPHA))) - ...
                                    1/3*(sum((mue*(N(k,:)*n1Temp'/THETA0).*lambda1N1Tilde.^(ALPHA))+...
                                    (mue*(N(k,:)*n1Temp'/THETA0).*lambda2N1Tilde.^(ALPHA))+...
                                    (mue*(N(k,:)*n1Temp'/THETA0).*lambda3N1Tilde.^(ALPHA)))));
                                S2N1 = 1/lambda2N1^2*(sum((mue*(N(k,:)*n1Temp'/THETA0).*lambda2N1Tilde.^(ALPHA))) - ...
                                    1/3*(sum((mue*(N(k,:)*n1Temp'/THETA0).*lambda1N1Tilde.^(ALPHA))+...
                                    (mue*(N(k,:)*n1Temp'/THETA0).*lambda2N1Tilde.^(ALPHA))+...
                                    (mue*(N(k,:)*n1Temp'/THETA0).*lambda3N1Tilde.^(ALPHA)))));
                                S3N1 = 1/lambda3N1^2*(sum((mue*(N(k,:)*n1Temp'/THETA0).*lambda3N1Tilde.^(ALPHA))) - ...
                                    1/3*(sum((mue*(N(k,:)*n1Temp'/THETA0).*lambda1N1Tilde.^(ALPHA))+...
                                    (mue*(N(k,:)*n1Temp'/THETA0).*lambda2N1Tilde.^(ALPHA))+...
                                    (mue*(N(k,:)*n1Temp'/THETA0).*lambda3N1Tilde.^(ALPHA)))));
                                SIsoN1 = S1N1*NeigN1(:,1)*NeigN1(:,1)'+S2N1*NeigN1(:,2)*NeigN1(:,2)'+S3N1*NeigN1(:,3)*NeigN1(:,3)';
                                S = SVolN1 + SIsoN1;
                                % Cauchy stress
                                sigma = 1/detF*FN1*S*FN1';
                                tempSpann = selectStress(sigma,options.stress,DIM);
                                Spann   = Spann + N(k,:)'*tempSpann*detJ1*wp(k);
                                Me = Me + (N(k,:)'*N(k,:))*detJ1*wp(k);
                            end
                            zw(j).Spann = Spann;
                            zw(j).pK = Me(:);
                        case 'neohookethermo'
                            Spann = zeros(size(N,2),1);
                            Me = zeros(nen);
                            for k = 1:NGP
                                tempSpann=0; %clear for each iteration
                                indx = DIM*k-(DIM-1):DIM*k;
                                detJ = det(J(:,indx)');
                                detJ1 = det(J1(:,indx)');
                                if detJ < 10*eps
                                    error('Jacobi determinant equal or less than zero.')
                                end
                                dNx = (J(:,indx)')\dNr(indx,:);
                                FN1 = edN1*dNx';
                                detF = det(FN1);
                                % Cauchy-Green tensor
                                CN1 = FN1'*FN1;
                                JN1 = (det(CN1))^0.5;
                                % Inverse strain tensor
                                CInvN1 = CN1\I;
                                % Delta Temp
                                thetaN1 = N(k,:)*n1Temp';
                                nuN1 = thetaN1 - theta0;
                                % Second PK
                                SN1 = mu*I + (lambda*(log(JN1)+JN1^2-JN1)-mu-3*beta*lambda*nuN1*(1/JN1+JN1-1/JN1*log(JN1)))*CInvN1;
                                % Cauchy stress
                                sigma = 1/detF*FN1*SN1*FN1';
                                tempSpann = selectStress(sigma,options.stress,DIM);
                                Spann   = Spann + N(k,:)'*tempSpann*detJ1*wp(k);
                                Me = Me + (N(k,:)'*N(k,:))*detJ1*wp(k);
                            end
                            zw(j).Spann = Spann;
                            zw(j).pK = Me(:);
                        case 'mooneyrivlinthermo'
                            Spann = zeros(size(N,2),1);
                            Me = zeros(nen);
                            for k = 1:NGP
                                tempSpann = 0; %clear for each iteration
                                indx = DIM*k-(DIM-1):DIM*k;
                                detJ = det(J(:,indx)');
                                detJ1 = det(J1(:,indx)');
                                if detJ < 10*eps
                                    error('Jacobi determinant equal or less than zero.')
                                end
                                dNx = (J(:,indx)')\dNr(indx,:);
                                FN1 = edN1*dNx';
                                detF = det(FN1);
                                % Right Cauchy-Green Tensor
                                CN1 = FN1'*FN1;
                                % Inverse strain tensor
                                CInvN1 = CN1\I;
                                % Invariant
                                I1N1 = trace(CN1);
                                I3N1 = (det(CN1))^0.5;
                                % Delta Tempy
                                thetaN1 = N(k,:)*n1Temp';
                                nuN1 = thetaN1 - theta0;
                                zw1 = -2*c1-4*c2;
                                % First derivatives of strain energy function
                                Psi_C = c1*I + c2*(I1N1*I-CN1) + 1/2*(zw1+DIM*beta*nuN1*zw1*I3N1^(-1))*CInvN1;
                                % Second PK
                                SN1 = 2*Psi_C;
                                % Cauchy stress
                                sigma = 1/detF*FN1*SN1*FN1';
                                tempSpann = selectStress(sigma,options.stress,DIM);
                                Spann   = Spann + N(k,:)'*tempSpann*detJ1*wp(k);
                                Me = Me + (N(k,:)'*N(k,:))*detJ1*wp(k);
                            end
                            zw(j).Spann = Spann;
                            zw(j).pK = Me(:);
                        case 'linearthermo'
                            Spann = zeros(size(N,2),1);
                            Me = zeros(nen);
                            DMat = [lambda+2*mu lambda lambda 0 0 0; lambda lambda+2*mu lambda 0 0 0; lambda lambda lambda+2*mu 0 0 0; 0 0 0 mu 0 0; 0 0 0 0 mu 0; 0 0 0 0 0 mu];
                            Dthermal = alphaMat*(3*lambda+2*mu) *  [1 1 1 0 0 0]';
                            uN1 = edN1-edRef;
                            for k = 1:NGP
                                indx = DIM*k-(DIM-1):DIM*k;
                                detJ = det(J(:,indx)');
                                detJ1=det(J1(:,indx)');
                                if detJ < 10*eps
                                    error('Jacobi determinant equal or less than zero.')
                                end
                                dNx = (J(:,indx)')\dNr(indx,:);
                                B = zeros(6,numDOFs);
                                B(1,1:DIM+1:end)=dNx(1,:);
                                B(2,2:DIM+1:end)=dNx(2,:);
                                B(3,3:DIM+1:end)=dNx(3,:);
                                B(4,1:DIM+1:end)=dNx(2,:);
                                B(4,2:DIM+1:end)=dNx(1,:);
                                B(5,2:DIM+1:end)=dNx(3,:);
                                B(5,3:DIM+1:end)=dNx(2,:);
                                B(6,1:DIM+1:end)=dNx(3,:);
                                B(6,3:DIM+1:end)=dNx(1,:);
                                pos = 1:numDOFs;
                                pos(4:4:end) = [];
                                % Delta Temp
                                thetaN1 = N(k,:)*n1Temp';
                                nuN1 = thetaN1 - theta0;
                                % sigma
                                sigma_v = DMat*B(1:6,pos)*uN1(:) - Dthermal*(1-wkalpha*nuN1)*nuN1;
                                sigma = [ sigma_v(1)  sigma_v(4)  sigma_v(6) ; sigma_v(4) sigma_v(2) sigma_v(5)  ; sigma_v(6) sigma_v(5) sigma_v(3) ];
                                tempSpann = selectStress(sigma,options.stress,DIM);
                                Spann   = Spann + N(k,:)'*tempSpann*detJ1*wp(k);
                                Me = Me + (N(k,:)'*N(k,:))*detJ1*wp(k);
                            end
                            zw(j).EDOF = edofLocal(:);
                            zw(j).Spann = Spann;
                            zw(j).pK = Me(:);
                        otherwise
                            error('Material model is currently not implemented.')
                    end
                otherwise
                    error('Element type is currently not implemented..')
            end
        end
        maxE = max(vertcat(zw(:).EDOF));
        zwRe = full(sparse(vertcat(zw(:).EDOF),1,vertcat(zw(:).Spann),maxE,1));
        fdof = 1:maxE;
        fdof(zwRe==0) = [];
        zwKe = sparse(vertcat(zw(:).pI), vertcat(zw(:).pJ), vertcat(zw(:).pK),maxE,maxE);
        colorData = zeros(maxE,1);
        colorData(fdof,1) = zwKe(fdof,fdof)\zwRe(fdof,1);
        
    end
end
end