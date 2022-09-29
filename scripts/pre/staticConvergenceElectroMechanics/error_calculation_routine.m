function out = error_calculation_routine(obj,analyticalVariablesCell)
%%Routine%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functionname:     error_calculation_routine.m
% Creation date:    22.12.2016
% Creator:          A.Janz
% Discription:      Routine calculates the L_2 error of the quantities            
% Based on:         --

% Modifications:    R.Pfefferkorn 10.5.2017, P. Kinon 09.03.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load objects
shapeFunctionObject = obj.shapeFunctionObject;
meshObject = obj.meshObject;
mixedFEObject = obj.mixedFEObject;

%% Shape functions
N = shapeFunctionObject.N;
dNrAll = shapeFunctionObject.dNr;
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
dimension = obj.dimension;

% Displacement
edof = meshObject.edof;
numberOfElements = size(edof,1);
qR = obj.qR;
numberOfInternalNodes = size(mixedFEObject.shapeFunction.N, 2);

% Electric Displacement
% FIXME different shapefunctions for mixed fields
N_D = mixedFEObject.shapeFunction.N;
N_C = mixedFEObject.shapeFunction.N;
N_G = mixedFEObject.shapeFunction.N;
N_c =  mixedFEObject.shapeFunction.N;
% numberOfNodes_D = size(N_D,2);
numberOfNodes_C = size(N_C,2);
% numberOfNodes_G = size(N_G,2);
% numberOfNodes_c = size(N_c,2);

% Output Variable
out(numberOfElements) = struct('QN1',[],'ana',[],'delta',[],'phi_analytical',[],'delta_phi',[],'DN1',[],'delta_D',[],'D_analytical',[],...
    'CN1',[],'delta_C',[],'C_analytical',[],'GN1',[],'delta_G',[],'G_analytical',[],'cN1',[],'delta_c',[],'c_analytical',[],'lambdaCN1',[],'delta_lambdaC',[],'lambdaC_analytical',[],...
    'lambdaGN1',[],'delta_lambdaG',[],'lambdaG_analytical',[],'lambdacN1',[],'delta_lambdac',[],'lambdac_analytical',[]);
for e = 1:numberOfElements
    edofLocal = edof(e,:);
    edofLocalReduced = edofLocal(1:numberOfNodes_C);
    dNr = dNrAll;
    
    save_QN1 = 0;
    save_anaDef_q = 0;
    save_delta_q = 0;
    
    save_delta_phi = 0;
    save_phi_analytical = 0;
    save_phiN1 = 0;
    
    save_DN1 = 0;
    save_delta_D = 0;
    save_D_analytical = 0;

    if strcmpi(obj.elementDisplacementType,'mixedSC')
        save_CN1 = 0;
        save_delta_C = 0;
        save_C_analytical = 0;
        save_GN1 = 0;
        save_delta_G = 0;
        save_G_analytical = 0;
        save_cN1 = 0;
        save_delta_c = 0;
        save_c_analytical = 0;

        save_lambdaCN1 = 0;
        save_delta_lambdaC = 0;
        save_lambdaC_analytical = 0;
        save_lambdaGN1 = 0;
        save_delta_lambdaG = 0;
        save_lambdaG_analytical = 0;
        save_lambdacN1 = 0;
        save_delta_lambdac = 0;
        save_lambdac_analytical = 0;
    end
    
    % analytical quantities
    ana_q = analyticalVariablesCell.q(edofLocal,1:dimension)';
    ana_Phi = analyticalVariablesCell.q(edofLocal,dimension+1)';
    ana_D = analyticalVariablesCell.D(edofLocalReduced,1:dimension)';
%     ana_D = analyticalVariablesCell.D(edofLocal,1:dimension)';
    if strcmpi(obj.elementDisplacementType,'mixedSC')
        ana_C = analyticalVariablesCell.C(edofLocalReduced,:)';
        ana_G = analyticalVariablesCell.G(edofLocalReduced,:)';
        ana_c = analyticalVariablesCell.c(edofLocalReduced,:)';
        ana_lambdaC = analyticalVariablesCell.lambdaC(edofLocalReduced,:)';
        ana_lambdaG = analyticalVariablesCell.lambdaG(edofLocalReduced,:)';
        ana_lambdac = analyticalVariablesCell.lambdac(edofLocalReduced,:)';
    end

    % numerical quantities 
    % extract dofs for element
    qN1 = obj.qN1(edofLocal,:);
    edAlphaN1 = obj.mixedFEObject.qN1(e,:);
    extractedDN1 = edAlphaN1(1:3*numberOfInternalNodes).';
    if strcmpi(obj.elementDisplacementType,'mixedSC')
        extractedCN1v = edAlphaN1(3*numberOfInternalNodes+1:9*numberOfInternalNodes).';
        extractedGN1v = edAlphaN1(9*numberOfInternalNodes+1:15*numberOfInternalNodes).';
        extractedcN1 = edAlphaN1(15*numberOfInternalNodes+1:16*numberOfInternalNodes).';
        extractedLambdaCN1v = edAlphaN1(16*numberOfInternalNodes+1:22*numberOfInternalNodes).';
        extractedLambdaGN1v = edAlphaN1(22*numberOfInternalNodes+1:28*numberOfInternalNodes).';
        extractedLambdacN1 = edAlphaN1(28*numberOfInternalNodes+1:29*numberOfInternalNodes).';
    end
    J = qR(edofLocal,1:dimension)'*dNr';
    
    % Run through all Gauss points
    for k = 1:numberOfGausspoints
        indx = dimension*k-(dimension-1):dimension*k;
        detJ = det(J(:,indx)');
        if detJ < 10*eps
            error('Jacobi determinant equal or less than zero.')
        end
        
        Phi = qN1(:,1:dimension)'*N(k,:)';
        Phi_ana = ana_q*N(k,:)';
        delta_q = Phi-Phi_ana;        
        % Summation over all Gauss Points
        save_QN1 = save_QN1 + Phi'*Phi*detJ*gaussWeight(k);
        save_anaDef_q = save_anaDef_q + Phi_ana'*Phi_ana*detJ*gaussWeight(k);
        save_delta_q = save_delta_q + delta_q'*delta_q*detJ*gaussWeight(k);
        
        phi = qN1(:,dimension+1)'*N(k,:)';
        phi_ana = ana_Phi*N(k,:)';
        delta_phi = phi-phi_ana;
        % Summation over all Gauss Points
        save_phiN1 = save_phiN1 + phi'*phi*detJ*gaussWeight(k);
        save_phi_analytical = save_phi_analytical + phi_ana'*phi_ana*detJ*gaussWeight(k);
        save_delta_phi = save_delta_phi + delta_phi'*delta_phi*detJ*gaussWeight(k);
        
        % Electric Displacements
        DN1 = reshape(extractedDN1, 3, []) * N_D(k, :)';
        D_ana = ana_D*N_D(k,:)';
%         D_ana = ana_D*N(k,:)';
        delta_D = DN1-D_ana;
        % Summation over all Gauss Points
        save_DN1 = save_DN1 + DN1'*DN1*detJ*gaussWeight(k);
        save_delta_D = save_delta_D + delta_D'*delta_D*detJ*gaussWeight(k);
        save_D_analytical = save_D_analytical + D_ana'*D_ana*detJ*gaussWeight(k);

        if strcmpi(obj.elementDisplacementType,'mixedSC')
            % Independent strains/conjungated stresses
            CN1v = reshape(extractedCN1v, 6, []) * N_C(k, :)';
            C_anav = ana_C*N_C(k,:)';
            CN1 = [CN1v(1) CN1v(4) CN1v(6); CN1v(4) CN1v(2) CN1v(5); CN1v(6) CN1v(5) CN1v(3)];
            C_ana = [C_anav(1) C_anav(4) C_anav(6); C_anav(4) C_anav(2) C_anav(5); C_anav(6) C_anav(5) C_anav(3)];
            delta_C = CN1-C_ana;

            GN1v = reshape(extractedGN1v, 6, []) * N_G(k, :)';
            GN1 = [GN1v(1) GN1v(4) GN1v(6); GN1v(4) GN1v(2) GN1v(5); GN1v(6) GN1v(5) GN1v(3)];
            G_anav = ana_G*N_G(k,:)';
            G_ana = [G_anav(1) G_anav(4) G_anav(6); G_anav(4) G_anav(2) G_anav(5); G_anav(6) G_anav(5) G_anav(3)];
            delta_G = GN1-G_ana;

            cN1 = extractedcN1.' * N_c(k, :)';
            c_ana = ana_c*N_c(k,:)';
            delta_c = cN1 - c_ana;

            % Independent strains/conjungated stresses
            lambda_CN1v = reshape(extractedLambdaCN1v, 6, [])*N_C(k, :)';
            lambda_CN1 = voigtToMatrix(lambda_CN1v, 'stress');
            lambdaC_anav = ana_lambdaC*N_C(k,:)';
            lambdaC_ana = [lambdaC_anav(1) lambdaC_anav(4) lambdaC_anav(6); lambdaC_anav(4) lambdaC_anav(2) lambdaC_anav(5); lambdaC_anav(6) lambdaC_anav(5) lambdaC_anav(3)];
            delta_lambdaC = lambda_CN1-lambdaC_ana;

            lambda_GN1v = reshape(extractedLambdaGN1v, 6, [])  * N_G(k, :)';
            lambda_GN1 = [lambda_GN1v(1) lambda_GN1v(4) lambda_GN1v(6); lambda_GN1v(4) lambda_GN1v(2) lambda_GN1v(5); lambda_GN1v(6) lambda_GN1v(5) lambda_GN1v(3)];
            lambdaG_anav = ana_lambdaG*N_G(k,:)';
            lambdaG_ana = [lambdaG_anav(1) lambdaG_anav(4) lambdaG_anav(6); lambdaG_anav(4) lambdaG_anav(2) lambdaG_anav(5); lambdaG_anav(6) lambdaG_anav(5) lambdaG_anav(3)];
            delta_lambdaG = lambda_GN1-lambdaG_ana;

            lambda_cN1  = extractedLambdacN1.' * N_c(k, :)';
            lambdac_ana = ana_lambdac*N_c(k,:)';
            delta_lambdac = lambda_cN1 - lambdac_ana;

            save_CN1 =  save_CN1 +  sum(sum(CN1.*CN1))*detJ*gaussWeight(k);
            save_delta_C = save_delta_C +sum(sum(delta_C.*delta_C))*detJ*gaussWeight(k);
            save_C_analytical = save_C_analytical + sum(sum(C_ana.*C_ana))*detJ*gaussWeight(k);
            save_GN1 =  save_GN1 +  sum(sum(GN1.*GN1))*detJ*gaussWeight(k);
            save_delta_G = save_delta_G +sum(sum(delta_G.*delta_G))*detJ*gaussWeight(k);
            save_G_analytical = save_G_analytical + sum(sum(G_ana.*G_ana))*detJ*gaussWeight(k);
            save_cN1 =  save_cN1 +  cN1*cN1*detJ*gaussWeight(k);
            save_delta_c = save_delta_c + delta_c*delta_c*detJ*gaussWeight(k);
            save_c_analytical = save_c_analytical + c_ana*c_ana*detJ*gaussWeight(k);

            save_lambdaCN1 =  save_lambdaCN1 +  sum(sum(lambda_CN1.*lambda_CN1))*detJ*gaussWeight(k);
            save_delta_lambdaC = save_delta_lambdaC +sum(sum(delta_lambdaC.*delta_lambdaC))*detJ*gaussWeight(k);
            save_lambdaC_analytical = save_lambdaC_analytical + sum(sum(lambdaC_ana.*lambdaC_ana))*detJ*gaussWeight(k);
            save_lambdaGN1 =  save_lambdaGN1 +  sum(sum(lambda_GN1.*lambda_GN1))*detJ*gaussWeight(k);
            save_delta_lambdaG = save_delta_lambdaG +sum(sum(delta_lambdaG.*delta_lambdaG))*detJ*gaussWeight(k);
            save_lambdaG_analytical = save_lambdaG_analytical + sum(sum(lambdaG_ana.*lambdaG_ana))*detJ*gaussWeight(k);
            save_lambdacN1 =  save_lambdacN1 + lambda_cN1*lambda_cN1*detJ*gaussWeight(k);
            save_delta_lambdac = save_delta_lambdac + delta_lambdac*delta_lambdac*detJ*gaussWeight(k);
            save_lambdac_analytical = save_lambdac_analytical + lambdac_ana*lambdac_ana*detJ*gaussWeight(k);
        end
    end
    
    % Output Variables
    out(e).QN1 = save_QN1;
    out(e).ana = save_anaDef_q;
    out(e).delta = save_delta_q;
    
    out(e).phiN1 = save_phiN1;
    out(e).phi_analytical = save_phi_analytical;
    out(e).delta_phi = save_delta_phi;

    out(e).DN1 = save_DN1;
    out(e).delta_D = save_delta_D;
    out(e).D_analytical = save_D_analytical;

    if strcmpi(obj.elementDisplacementType,'mixedSC')
        out(e).CN1 = save_CN1;
        out(e).delta_C = save_delta_C;
        out(e).C_analytical = save_C_analytical;
        out(e).GN1 = save_GN1;
        out(e).delta_G = save_delta_G;
        out(e).G_analytical = save_G_analytical;
        out(e).cN1 = save_cN1;
        out(e).delta_c = save_delta_c;
        out(e).c_analytical = save_c_analytical;

        out(e).lambdaCN1 = save_lambdaCN1;
        out(e).delta_lambdaC = save_delta_lambdaC;
        out(e).lambdaC_analytical = save_lambdaC_analytical;
        out(e).lambdaGN1 = save_lambdaGN1;
        out(e).delta_lambdaG = save_delta_lambdaG;
        out(e).lambdaG_analytical = save_lambdaG_analytical;
        out(e).lambdacN1 = save_lambdacN1;
        out(e).delta_lambdac = save_delta_lambdac;
        out(e).lambdac_analytical = save_lambdac_analytical;
    end
end  
end