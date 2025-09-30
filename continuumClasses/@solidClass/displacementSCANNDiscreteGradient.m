function [rData, kData, elementEnergy, array] = displacementSCANN3DiscreteGradient(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTSCANNMOONEYRIVLIN3DISCRETEGRADIENT Element routine of class solidClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine  covering nonlinear
% mechanical processes employing a trained artificial neural network (ANN) 
% constitutive model which was trained with an isotropic, hyperelastic 
% Mooney-Rivlin ('MooneyRivlin') model.
% Implementation is due to work-conjugated 2nd PK-stress tensor and Cauchy-
% Green strain tensor ('SC').
% The routine is suitable for dynamic simulations where an energy-momentum 
% consistent integration scheme is used which is based on the conecpt of 
% discrete gradients in the sense of Gonzalez ('DiscreteGradient').
% In particular algorthmic derivatives for IC, IIC, J and J* are employed.
% -> does not work properly
%
% CALL
% displacementSCANNMooneyRivlinDiscreteGradient(obj,setupObject,computePostData)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
%
% REFERENCE
% https://doi.org/10.1016/j.cma.2018.01.013
%
% SEE ALSO
% displacementSCANNMooneyRivlinEndpoint
% displacementSCANNMooneyRivlinMidpoint
% displacementSCMooneyRivlinMidpoint,
% displacementSCMooneyRivlinDiscreteGradient
%
% CREATOR(S)
% Marlon Franke

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
numericalTangentObject = obj.numericalTangentObject;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof(e, :);

dimension = obj.dimension;

% aquire material data
a = materialObject.a;
b = materialObject.b;
c = materialObject.c;
d = materialObject.d;

% already trained ANN model
normalizeFactor = 1;
if 0
    b2 = -1.257220077514648438e+02;
    w1 =  [ 9.268754005432128906e+00 6.275703907012939453e-01 8.978516578674316406e+00 9.184016227722167969e+00 8.944575309753417969e+00 8.825480461120605469e+00 4.117710050195455551e-03 9.206143379211425781e+00;
        1.117671847343444824e+00 1.695705890655517578e+00 1.225380182266235352e+00 1.790273785591125488e+00 2.119191884994506836e+00 3.604783296585083008e+00 3.489471077919006348e-01 1.372677236795425415e-01;
        4.465130623430013657e-03 5.173204541206359863e-01 4.450271837413311005e-03 4.819013178348541260e-03 4.955890122801065445e-03 6.109093781560659409e-03 3.643061399459838867e+00 3.700215835124254227e-03;
        2.500126075744628906e+01 2.176816558837890625e+01 2.508148765563964844e+01 2.498847198486328125e+01 2.479306983947753906e+01 2.380420684814453125e+01 2.061729812622070312e+01 2.552124214172363281e+01];
    w2 = [  7.630357265472412109e+00;
        3.028271198272705078e+00;
        7.628017425537109375e+00;
        7.640446186065673828e+00;
        7.606583595275878906e+00;
        7.949126720428466797e+00;
        4.874103546142578125e+01;
        7.896135807037353516e+00];
    b1 = [  -3.441793203353881836e+00;
        1.242680740356445312e+01;
        -2.851085186004638672e+00;
        -4.938632488250732422e+00;
        -5.322002410888671875e+00;
        -1.029190444946289062e+01;
        9.544090270996093750e+00;
        -9.355480670928955078e-01];
else
    normalizeFactor = 1e3;
    w1 = [  -0.000000000000000000e+00 -0.000000000000000000e+00 -0.000000000000000000e+00 3.498095989227294922e+00;
            -0.000000000000000000e+00 -0.000000000000000000e+00 6.021448597311973572e-02 9.343179702758789062e+00;
            1.066247582435607910e+00 7.192744612693786621e-01 2.154079377651214600e-01 3.736973181366920471e-02;
            8.472843766212463379e-01 3.536386787891387939e-02 2.211225219070911407e-04 3.044792413711547852e+00;
            -0.000000000000000000e+00 -0.000000000000000000e+00 -0.000000000000000000e+00 3.267039537429809570e+00;
            -0.000000000000000000e+00 -0.000000000000000000e+00 -0.000000000000000000e+00 3.178473234176635742e+00;
            -0.000000000000000000e+00 -0.000000000000000000e+00 -0.000000000000000000e+00 3.943647623062133789e+00;
            -0.000000000000000000e+00 -0.000000000000000000e+00 -0.000000000000000000e+00 3.375712156295776367e+00];
    w1 = w1';
    b1 = [  -5.596243381500244141e+00 4.461702823638916016e+00 1.894991379231214523e-03 2.932721614837646484e+00 -5.644768238067626953e+00 -5.689310073852539062e+00 -5.542060375213623047e+00 -5.609594821929931641e+00]';
    w2 = [  -0.000000000000000000e+00 1.564468592405319214e-01 8.862167596817016602e-02 4.042990505695343018e-01 -0.000000000000000000e+00 -0.000000000000000000e+00 -0.000000000000000000e+00 -0.000000000000000000e+00]';
    b2 = -1.550519585609436035e+00;
end
alpha = 5e2;
% aquire the nodal values of the variables for the current element
edR = obj.qR(edof, 1:dimension).';
edN = obj.qN(edof, 1:dimension).';
edN1 = dofs.edN1;
edN05 = 0.5 * (edN + edN1);

% initialize residual & tangent
RX = rData{1};
KXX = kData{1, 1};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

% initialize flagDiscreteGradient
if ~flagNumericalTangent
    initializeFlagDiscreteGradient(numericalTangentObject, shapeFunctionObject.numberOfGausspoints, 4)
end
flagDiscreteGradient = numericalTangentObject.flagDiscreteGradient;

numberOfDOFs = dimension*size(edof, 2);

activation = computeActivationFunction();

netANN.w1 = w1;
netANN.w2 = w2;
netANN.b1 = b1;
netANN.b2 = b2;
netANN.activation = activation;
netANN.alpha = alpha;
netANN.normalizeFactor = normalizeFactor;

% Run through all Gauss points
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
    % Deformation gradient
    FN1 = edN1 * dN_X_I';
    FN = edN * dN_X_I';
    FN05 = edN05 * dN_X_I';
    % B-matrix (current configuration)
    BN1 = BMatrix(dN_X_I, FN1);
    BN05 = BMatrix(dN_X_I, FN05);
%     BN05 = 2 * BN05;
    % Right Cauchy-Green tensor
    CN1 = FN1' * FN1;
    CN = FN' * FN;
    CN05 = FN05' * FN05;
    CAlgo = 0.5 * (CN + CN1); %averaged
    % Cofactor
    GN1 = 0.5 * wedge(CN1, CN1); %cof(CN1)=det(CN1)*CN1^-T
    GN = 0.5 * wedge(CN, CN);
    GAlgo = 0.5 * (GN + GN1); %averaged
    GN05 = 0.5 * wedge(CN05, CN05); %mid-configuration

    % invariants for ANN model
    IN = trace(CN);
    IIN = trace(GN);
    IN1 = trace(CN1);
    IIN1 = trace(GN1);
    IN05 = trace(CN05);
    IIN05 = trace(GN05);
    cN = det(CN);
    JN = sqrt(cN);
    JstarN = -JN;
    cN1 = det(CN1);
    JN1 = sqrt(cN1);
    JstarN1 = -JN1;
    cN05 = det(CN05);
    JN05 = sqrt(cN05);
    JstarN05 = -JN05;

    % first derivative
    xN = [IN;IIN;JN;JstarN];
    xN1 = [IN1;IIN1;JN1;JstarN1];
    xN05 = [IN05;IIN05;JN05;JstarN05];

    h1N05 = w1'*xN05 + b1;

    WN1 = computeEnergyANN(IN1,IIN1,JN1,JstarN1,netANN);
    
    dh2dh1N05 = diag(activation.derivative(h1N05));
    dWvol_JN05 = 2*alpha*(JN05-1);
    dWvol_JstarN05 = 0;
    dW_xN05 = w2'*(dh2dh1N05*w1');
    dW_xN05(3) = dW_xN05(3) + dWvol_JN05;
    dW_xN05(4) = dW_xN05(4) + dWvol_JstarN05;
    dW_IN05 = dW_xN05(1);
    dW_IIN05 = dW_xN05(2);
    dW_JN05 = dW_xN05(3);
    dW_JstarN05 = dW_xN05(4);    

    deltaI = IN1 - IN;
    normDeltaI = abs(deltaI);
    deltaII = IIN1 - IIN;
    normDeltaII = abs(deltaII);
    deltaJ = JN1 - JN;
    normDeltaJ = abs(deltaJ);
    deltaJstar = JstarN1 - JstarN;
    normDeltaJstar = abs(deltaJstar);
    
    % discrete gradient
    if ~flagNumericalTangent
        discreteGradientThreshold = 1e-16; % FIXME: include into solidObject
        discreteGradientCondition_I = (normDeltaI > discreteGradientThreshold);
        discreteGradientCondition_II = (normDeltaII > discreteGradientThreshold);
        discreteGradientCondition_J = (normDeltaJ > discreteGradientThreshold);
        discreteGradientCondition_Jstar = (normDeltaJstar > discreteGradientThreshold);
%         discreteGradientCondition_D = setupObject.newton.step > 1;
        if (setupObject.newton.step > 1) && (~discreteGradientCondition_I || ~discreteGradientCondition_II || ~discreteGradientCondition_J || ~discreteGradientCondition_Jstar)
            disp(1)
        end
    else
        discreteGradientCondition_I = flagDiscreteGradient(k, 1);
        discreteGradientCondition_II = flagDiscreteGradient(k, 2);
        discreteGradientCondition_J = flagDiscreteGradient(k, 3);
        discreteGradientCondition_Jstar = flagDiscreteGradient(k, 4);
    end
    if discreteGradientCondition_I
        flagDiscreteGradient(k, 1) = 1;
        dW_IN1 = (computeEnergyANN(IN1,IIN,JN,JstarN,netANN)-computeEnergyANN(IN,IIN,JN,JstarN,netANN))/deltaI;
        dW_IN = (computeEnergyANN(IN1,IIN1,JN1,JstarN1,netANN)-computeEnergyANN(IN,IIN1,JN1,JstarN1,netANN))/deltaI;
        DW_I = 1/2*(dW_IN1+dW_IN);
    else
        flagDiscreteGradient(k, 1) = 0;
        DW_I = dW_IN05;
    end
    if discreteGradientCondition_II
        flagDiscreteGradient(k, 2) = 1;
        dW_IIN1 = (computeEnergyANN(IN1,IIN1,JN,JstarN,netANN)-computeEnergyANN(IN1,IIN,JN,JstarN,netANN))/deltaII;
        dW_IIN = (computeEnergyANN(IN,IIN1,JN1,JstarN1,netANN)-computeEnergyANN(IN,IIN,JN1,JstarN1,netANN))/deltaII;
        DW_II = 1/2*(dW_IIN1+dW_IIN);
    else
        flagDiscreteGradient(k, 2) = 0;
        DW_II = dW_IIN05;
    end    
    if discreteGradientCondition_J
        flagDiscreteGradient(k, 3) = 1;
        dW_JN1 = (computeEnergyANN(IN1,IIN1,JN1,JstarN,netANN)-computeEnergyANN(IN1,IIN1,JN,JstarN,netANN))/deltaJ;
        dW_JN = (computeEnergyANN(IN,IIN,JN1,JstarN1,netANN)-computeEnergyANN(IN,IIN,JN,JstarN1,netANN))/deltaJ;
        DW_J = 1/2*(dW_JN1+dW_JN);
    else
        flagDiscreteGradient(k, 3) = 0;
        DW_J = dW_JN05;
    end    
    if discreteGradientCondition_Jstar
        flagDiscreteGradient(k, 4) = 1;
        dW_JstarN1 = (computeEnergyANN(IN1,IIN1,JN1,JstarN,netANN)-computeEnergyANN(IN1,IIN1,JN1,JstarN,netANN))/deltaJstar;
        dW_JstarN = (computeEnergyANN(IN,IIN,JN,JstarN1,netANN)-computeEnergyANN(IN,IIN,JN,JstarN,netANN))/deltaJstar;
        DW_Jstar = 1/2*(dW_JstarN1+dW_JstarN);
    else
        flagDiscreteGradient(k, 4) = 0;
        DW_Jstar = dW_JstarN05;
    end    

    dIC_CN05 = eye(dimension);
    dIIC_GN05 = eye(dimension);
    dJ_cN05 = 1/2*(cN05)^(-1/2);
    dJstar_cN05 = -dJ_cN05;
      
    DW_C = DW_I*dIC_CN05;
    DW_G = DW_II*dIIC_GN05;
    DW_c = DW_J*dJ_cN05 + DW_Jstar*dJstar_cN05;

%     SN05 = 2*(DW_C + wedge(DW_G, CN05) + DW_c*GN05);
%     SN05 = 2*(DW_C + wedge(DW_G, CAlgo) + DW_c * 1 / 2 * (wedge(CAlgo, CAlgo) + GAlgo));
     SN05 = 2*(DW_C + wedge(DW_G, CAlgo) + DW_c * 1 / 3 * (wedge(CAlgo, CAlgo) + GAlgo));
    SN05_v = [SN05(1, 1); SN05(2, 2); SN05(3, 3); SN05(1, 2); SN05(2, 3); SN05(1, 3)];

    if ~computePostData
        % Residual
        RX = RX + BN05' * SN05_v * detJ * gaussWeight(k);
        % Strain energy
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + WN1 * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        PN1 = FN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FN1) * PN1 * FN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    if ~flagNumericalTangent
        numericalTangentObject.flagDiscreteGradient = flagDiscreteGradient;
    end
    rData{1} = RX;
    kData{1, 1} = KXX;
end
end

function D = secDiffOperator(A)
D = zeros(6, 6);
D(1, 2) = A(3, 3);
D(1, 3) = A(2, 2);
D(1, 5) = -0.5 * (A(2, 3) + A(3, 2));
D(2, 1) = A(3, 3);
D(2, 3) = A(1, 1);
D(2, 6) = -0.5 * (A(1, 3) + A(3, 1));
D(3, 1) = A(2, 2);
D(3, 2) = A(1, 1);
D(3, 4) = -0.5 * (A(1, 2) + A(2, 1));
D(4, 3) = -0.5 * (A(1, 2) + A(2, 1));
D(4, 4) = -0.5 * A(3, 3);
D(4, 5) = 0.25 * (A(1, 3) + A(3, 1));
D(4, 6) = 0.25 * (A(2, 3) + A(3, 2));
D(5, 1) = -0.5 * (A(2, 3) + A(3, 2));
D(5, 4) = 0.25 * (A(1, 3) + A(3, 1));
D(5, 5) = -0.5 * A(1, 1);
D(5, 6) = 0.25 * (A(1, 2) + A(2, 1));
D(6, 2) = -0.5 * (A(1, 3) + A(3, 1));
D(6, 4) = 0.25 * (A(2, 3) + A(3, 2));
D(6, 5) = 0.25 * (A(1, 2) + A(2, 1));
D(6, 6) = -0.5 * A(2, 2);
end
function [activation] = computeActivationFunction()
activation.function = @(x) log(1+exp(x));
activation.derivative = @(x) exp(x)./(1+exp(x));
activation.secondDerivative = @(x) exp(x)./((1+exp(x)).^2);
% activation.function = @(x) 1./(1+exp(-x));
% activation.derivative = @(x) x.*(1 - x);
end
function W = computeEnergyANN(I,II,J,Jstar,netANN)
w1 = netANN.w1;
w2 = netANN.w2;
b1 = netANN.b1;
b2 = netANN.b2;
activation = netANN.activation;
alpha = netANN.alpha;
normalizeFactor = netANN.normalizeFactor;
x = [I;II;J;Jstar];
h1 = w1'*x + b1;
h2 = activation.function(h1);
Wdev = w2'*h2 + b2;
Wvol = alpha*(J-1)^2;
W = normalizeFactor*Wdev + Wvol;
end

%     dW_I = zeros(4,1);
%     for indexI = 1:4
%         for aIndex = 1:8
%             dW_I(indexI,1) = w2(aIndex)*activation.derivative(h1(aIndex))*w1(indexI,aIndex);
%         end
%     end
%     d2W_I_I = zeros(4);
%     for indexI = 1:4
%         for indexJ = 1:4
%             for aIndex = 1:8
%                 d2W_I_I(indexI,indexJ) = d2W_I_I(indexI,indexJ) + w2(aIndex)*activation.secondDerivative(h1(aIndex))*w1(indexI,aIndex)*w1(indexJ,aIndex);
%             end
%         end
%     end
%     h = 1e-6;
%     d2WNum = zeros(numel(x));
%     for ii = 1:numel(x)
% %         x = [IC;IIC;J;Jstar];
%         xNum = x;
%         xNum(ii) = xNum(ii) + h;
%         h1Num = w1'*xNum + b1;
%         dh2dh1Num = diag(activation.derivative(h1Num));
%         dWvol_J = 2*alpha*(J-1);
%         dW_INum = w2'*(dh2dh1Num*w1');
%         dW_INum(3) = dW_INum(3) + dWvol_J;
%         d2WNum(:,ii) = (dW_INum-dW_I)/h;
%     end

