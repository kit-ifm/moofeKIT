function out = displacementSCMooneyRivlinDiscreteGradient(obj,varargin)
%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% out = wedgeSC_mooneyRivlin_discreteGradient(obj,'PropertyName',PropertyValue)
%
% Description
%
% Mooney Rivlin strain-energy function, algorithmic stress formula, energs conserving.
%
% 20.09.2016 A.Janz

globalFullEdof = obj.globalFullEdof;
edof = obj.edof;
numberOfGausspoints = obj.numberOfGausspoints;
gaussWeight = obj.shapeFunctions.gaussWeight;
NAll = obj.shapeFunctions.N';
dNrAll = obj.shapeFunctions.dNr';
qR = obj.qR;
qN = obj.qN;
qN1 = obj.qN1;
numberOfElements = size(globalFullEdof,1);
numberOfDOFs = size(globalFullEdof,2);
dimension = obj.dimension;
a = obj.materialData.a;
b = obj.materialData.b;
c = obj.materialData.c;
d = obj.materialData.d;
I = eye(dimension);
out(numberOfElements) = struct('edofE',[],'Re',[],'ePot',[],'pI',[],'pJ',[],'pK',[]);
for e = 1:numberOfElements
    dNr = [];
    [edofH1,edofH2] = expandEdof(globalFullEdof(e,:));
    out(e).pI = edofH1';
    out(e).pJ = edofH2';
    out(e).edofE = double(globalFullEdof(e,:))';
    % Element routine
    Re = zeros(numberOfDOFs,1);
    Ke = zeros(numberOfDOFs);
    ePot = 0;
    edN1 = qN1(edof(e,:),1:dimension)';
    edN = qN(edof(e,:),1:dimension)';
    edN05 = 0.5*(edN+edN1);
    
    J = qR(edof(e,:),1:dimension)'*dNrAll;
    % Run through all Gauss points
    for k = 1:numberOfGausspoints
        indx = dimension*k-(dimension-1):dimension*k;
        detJ = det(J(:,indx)');
        if detJ < 10*eps
            error('Jacobi determinant equal or less than zero.')
        end
%         dNx = (J(:,indx)')\dNr(indx,:);
        dNx = (J(:,indx)')\(dNrAll(:,indx))';
        % Deformation gradient
        FN1 = edN1*dNx';
        FN = edN*dNx';
        FN05 = edN05*dNx';
        % B-matrix (midpoint configuration)
        BN05 = zeros(6,numberOfDOFs);
        BN05(1,1:3:end)=FN05(1,1)*dNx(1,:);
        BN05(1,2:3:end)=FN05(2,1)*dNx(1,:);
        BN05(1,3:3:end)=FN05(3,1)*dNx(1,:);
        BN05(2,1:3:end)=FN05(1,2)*dNx(2,:);
        BN05(2,2:3:end)=FN05(2,2)*dNx(2,:);
        BN05(2,3:3:end)=FN05(3,2)*dNx(2,:);
        BN05(3,1:3:end)=FN05(1,3)*dNx(3,:);
        BN05(3,2:3:end)=FN05(2,3)*dNx(3,:);
        BN05(3,3:3:end)=FN05(3,3)*dNx(3,:);
        BN05(4,1:3:end)=FN05(1,1)*dNx(2,:)+FN05(1,2)*dNx(1,:);
        BN05(4,2:3:end)=FN05(2,1)*dNx(2,:)+FN05(2,2)*dNx(1,:);
        BN05(4,3:3:end)=FN05(3,1)*dNx(2,:)+FN05(3,2)*dNx(1,:);
        BN05(5,1:3:end)=FN05(1,2)*dNx(3,:)+FN05(1,3)*dNx(2,:);
        BN05(5,2:3:end)=FN05(2,2)*dNx(3,:)+FN05(2,3)*dNx(2,:);
        BN05(5,3:3:end)=FN05(3,2)*dNx(3,:)+FN05(3,3)*dNx(2,:);
        BN05(6,1:3:end)=FN05(1,1)*dNx(3,:)+FN05(1,3)*dNx(1,:);
        BN05(6,2:3:end)=FN05(2,1)*dNx(3,:)+FN05(2,3)*dNx(1,:);
        BN05(6,3:3:end)=FN05(3,1)*dNx(3,:)+FN05(3,3)*dNx(1,:);
        BN05 = 2*BN05;
        % B-matrix (current configuration)
        BN1 = zeros(6,numberOfDOFs);
        BN1(1,1:3:end)=FN1(1,1)*dNx(1,:);
        BN1(1,2:3:end)=FN1(2,1)*dNx(1,:);
        BN1(1,3:3:end)=FN1(3,1)*dNx(1,:);
        BN1(2,1:3:end)=FN1(1,2)*dNx(2,:);
        BN1(2,2:3:end)=FN1(2,2)*dNx(2,:);
        BN1(2,3:3:end)=FN1(3,2)*dNx(2,:);
        BN1(3,1:3:end)=FN1(1,3)*dNx(3,:);
        BN1(3,2:3:end)=FN1(2,3)*dNx(3,:);
        BN1(3,3:3:end)=FN1(3,3)*dNx(3,:);
        BN1(4,1:3:end)=FN1(1,1)*dNx(2,:)+FN1(1,2)*dNx(1,:);
        BN1(4,2:3:end)=FN1(2,1)*dNx(2,:)+FN1(2,2)*dNx(1,:);
        BN1(4,3:3:end)=FN1(3,1)*dNx(2,:)+FN1(3,2)*dNx(1,:);
        BN1(5,1:3:end)=FN1(1,2)*dNx(3,:)+FN1(1,3)*dNx(2,:);
        BN1(5,2:3:end)=FN1(2,2)*dNx(3,:)+FN1(2,3)*dNx(2,:);
        BN1(5,3:3:end)=FN1(3,2)*dNx(3,:)+FN1(3,3)*dNx(2,:);
        BN1(6,1:3:end)=FN1(1,1)*dNx(3,:)+FN1(1,3)*dNx(1,:);
        BN1(6,2:3:end)=FN1(2,1)*dNx(3,:)+FN1(2,3)*dNx(1,:);
        BN1(6,3:3:end)=FN1(3,1)*dNx(3,:)+FN1(3,3)*dNx(1,:);
        BN1 = 2*BN1;
        % Right Cauchy-Green tensor
        CN1 = FN1'*FN1;
        CN = FN'*FN;
        CN05 = 0.5*(CN+CN1);            %averaged
        
        % Cofactor
        GN1 = 0.5*wedge(CN1,CN1);
        GN = 0.5*wedge(CN,CN);
        GN05 = 0.5*(GN+GN1);            %averaged
        GMid = 0.5*wedge(CN05,CN05);    %mid-configuration
        
        % Third invariant
        I3N1 = det(CN1);
        I3N = det(CN);
        I3N05 = 0.5*(I3N+I3N1);
        
        % Strain energy function
        ePot = ePot + (a*(trace(CN1)-3) + b*(trace(GN1)-3) - d*log(sqrt(I3N1)) + c/2*(sqrt(I3N1)-1)^2)*detJ*gaussWeight(k);
        
        % Derivative of the strain energy function
        Sigma_C = a*I;                       % S_C = 2*DW/DC
        Sigma_G = b*I;                       % S_G = 2*DW/DG
        CTest = norm(CN1-CN,'fro');
        if  CTest < 10^(-11)
            Sigma_I = -d/(2*I3N05)+c/2*(1-1/sqrt(I3N05));   % S_I =   2*DW/DI
            Sigma_I_I = d/(2*I3N05^2)+c/(2*(I3N05)^(3/2));    % S_I_I = 4*D2W/DI3^2
        else %Greenspan formula '84
            Sigma_I=  (-d*log(sqrt(I3N1)) + c/2*(sqrt(I3N1)-1)^2 - ( - d*log(sqrt(I3N)) + c/2*(sqrt(I3N)-1)^2)) /(I3N1-I3N) ;
            Sigma_I_I = ((d/I3N1 - (c*(I3N1^(1/2) - 1))/I3N1^(1/2))/(I3N - I3N1) + (c*(I3N^(1/2) - 1)^2 - c*(I3N1^(1/2) - 1)^2 - 2*d*log(I3N^(1/2)) + 2*d*log(I3N1^(1/2)))/(I3N - I3N1)^2);
        end
        
        % Second Piola Kirchhoff stress tensor
        SN05 = 2*(Sigma_C + wedge(Sigma_G,CN05) + Sigma_I*1/3*(wedge(CN05,CN05)+GN05));
        SN05_v = [SN05(1,1); SN05(2,2); SN05(3,3); SN05(1,2); SN05(2,3); SN05(1,3)];
        
        % Residual
        Re = Re + BN05'*1/2*SN05_v*detJ*gaussWeight(k);
        
        % Tangent
        
        % Derivative of wedge(Sigma_G,CN05)
        
        D=(Sigma_G);
        SecDiffOperator=[0          D(3,3)  D(2,2)      0        -D(3,2)  0        ;
            D(3,3)     0       D(1,1)      0         0       -D(3,1)  ;
            D(2,2)     D(1,1)  0           -D(2,1)	  0       0        ;
            0          0       -D(2,1)     -D(3,3)   D(3,1)  D(2,3)   ;
            -D(3,2)	    0       0           D(3,1)   -D(1,1)  D(1,2)   ;
            0          -D(3,1)	0           D(3,2)    D(2,1)  -D(2,2)] ;
        SecDiffOperator(4:6,4:6) = 0.5*SecDiffOperator(4:6,4:6);
        Kmat1 = SecDiffOperator;
        
        % Derivative of Sigma_I*1/3*(wedge(CN05,CN05)+GN05)  part I
        
        D=(2/3*(CN05+0.5*CN1));
        SecDiffOperator=[0          D(3,3)  D(2,2)      0        -D(3,2)  0        ;
            D(3,3)     0       D(1,1)      0         0       -D(3,1)  ;
            D(2,2)     D(1,1)  0           -D(2,1)	  0       0        ;
            0          0       -D(2,1)     -D(3,3)   D(3,1)  D(2,3)   ;
            -D(3,2)	    0       0           D(3,1)   -D(1,1)  D(1,2)   ;
            0          -D(3,1)	0           D(3,2)    D(2,1)  -D(2,2)] ;
        SecDiffOperator(4:6,4:6) = 0.5*SecDiffOperator(4:6,4:6);
        
        Kmat2 = Sigma_I*SecDiffOperator;
        
        % Derivative of Sigma_I*1/3*(wedge(CN05,CN05)+GN05)  part II
        GN05_v = [GN05(1,1); GN05(2,2); GN05(3,3); GN05(1,2); GN05(2,3); GN05(1,3)];
        GMid_v = [GMid(1,1); GMid(2,2); GMid(3,3); GMid(1,2); GMid(2,3); GMid(1,3)];
        GN1_v = [GN1(1,1); GN1(2,2); GN1(3,3); GN1(1,2); GN1(2,3); GN1(1,3)];
        Kmat3 = Sigma_I_I*((GN05_v+2*GMid_v)*GN1_v')*1/3;
        
        % Assembly of elasticity tensor
        ELA = (Kmat1 + Kmat2 +  Kmat3);
        
        A1 = dNx'*SN05*dNx*detJ*gaussWeight(k);
        MAT = zeros(numberOfDOFs);
        for g = 1:dimension
            MAT(g:dimension:numberOfDOFs,g:dimension:numberOfDOFs) = A1;
        end
        Ke = Ke + 0.5*BN05'*ELA*BN1*detJ*gaussWeight(k) + 0.5*MAT;
    end
    out(e).Re = Re;
    out(e).pK = Ke(:);
    out(e).ePot = ePot;
end
end