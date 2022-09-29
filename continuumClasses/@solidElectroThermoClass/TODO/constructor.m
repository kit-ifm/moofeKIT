function obj = constructor(obj,mesh,varargin)
%% Constructur of the SolidElectroThermo class
%
% 08.04.2019 M. Franke 

if isstruct(mesh) && isfield(mesh,'SolidElectroThermo') && numel(mesh.SolidElectroThermo) > 0
    for i = 1:numel(mesh.SolidElectroThermo)
        meshZw = mesh.SolidElectroThermo{i};
        obj(i).DESIGNATOR = meshZw.DESIGNATOR;
        obj(i).Ansatzfunction = meshZw.Ansatzfunction;
        obj(i).DIM = meshZw.DIM;
        obj(i).ADDFIELDS = meshZw.ADDFIELDS;
        obj(i).NGP = meshZw.NGP;
        obj(i).REFPHI = meshZw.REFPHI;
        obj(i).REFTEMP = meshZw.REFTEMP;
        obj(i).numericalTangent = meshZw.numericalTangent;
        obj(i).ORDER = meshZw.ORDER;
        obj(i).ORDER_MIXED_VAR = meshZw.ORDER_MIXED_VAR;
        obj(i).ORDER_MIXED_VAR_EL = meshZw.ORDER_MIXED_VAR_EL;
        obj(i).EDOF = cell(1,size(meshZw.ELEMENTS.EDOF,1));
        for ii = 1:size(meshZw.ELEMENTS.EDOF,1)
            obj(i).EDOF{ii} = meshZw.ELEMENTS.EDOF(ii,:);
        end
        %% Electrical and mechanical prestress
        if numel(meshZw.QREF) == 0
            if numel(obj.REFPHI) == 1
                obj(i).QREF = [meshZw.NODES ones(size(meshZw.NODES,1),1)*meshZw.REFPHI(:) ones(size(meshZw.NODES,1),1)*meshZw.REFTEMP(:)];
            else
                obj(i).QREF = [meshZw.NODES meshZw.REFPHI(:) meshZw.REFTEMP(:)];
            end
        else
            obj(i).QREF= [meshZw.QREF ones(size(meshZw.NODES,1),1)*meshZw.REFPHI(:) ones(size(meshZw.NODES,1),1)*meshZw.REFTEMP(:)];
        end
        if numel(meshZw.PHI) ==0
            obj(i).QN = [meshZw.NODES ones(size(meshZw.NODES,1),1)*meshZw.REFPHI(:) ones(size(meshZw.NODES,1),1)*meshZw.REFTEMP(:)];
            obj(i).QN1 = [meshZw.NODES ones(size(meshZw.NODES,1),1)*meshZw.REFPHI(:) ones(size(meshZw.NODES,1),1)*meshZw.REFTEMP(:)];
        else
            obj(i).QN = [meshZw.NODES meshZw.PHI(:) meshZw.TEMP(:)];
            obj(i).QN1 = [meshZw.NODES meshZw.PHI(:) meshZw.TEMP(:)];
        end
        obj(i).VN = [meshZw.V zeros(size(meshZw.NODES,1),1), zeros(size(meshZw.NODES,1),1)];
        obj(i).GRAVITY = mesh.Config.GRAVITY;
        obj(i).MAT = meshZw.Material;
        obj(i).STAGGERED = mesh.Config.STAGGERED;
        obj(i).IDENTIFIER = meshZw.IDENTIFIER;
        obj(i).WEIGHTS = meshZw.WEIGHTS;
        obj(i).HCONTROL = meshZw.HCONTROL;
        obj(i).STMatrix = meshZw.STMatrix;
        obj(i).ELEMENTS = meshZw.ELEMENTS;
        obj(i).KNOTS = meshZw.KNOTS;
        obj(i).ORDER = meshZw.ORDER;
        if isstruct(meshZw.levelVec)
            obj(i).levelVec = meshZw.levelVec.lv;
        else
            obj(i).levelVec = meshZw.levelVec;
        end
        if isstruct(meshZw.levelStruct)
            obj(i).levelStruct = meshZw.levelStruct.lv;
        else
            obj(i).levelStruct = meshZw.levelStruct;
        end
        
        % Element routine
        obj(i).ROUTINE = selectRoutine(obj(i),mesh.Config.INTEGRATOR);
        % Shape functions
        if ~strcmp(obj(i).Ansatzfunction,'Lagrange') && ~strcmp(obj(i).Ansatzfunction,'HNURBS')
            error('Type of shape function is currently not implemented')
        end
        structure = struct('DIM',[],'NGP',[],'Element',[],'HCONTROL',[],'Ansatz',obj(i).Ansatzfunction);
        structure.DIM = obj(i).DIM;
        structure.NGP = obj(i).NGP;
        structure.ORDER = meshZw.ORDER;
        structure.KNOTS = meshZw.KNOTS;
        structure.nders = 1;
        structure.NEN = meshZw.ELEMENTS.nen;
        structure.INC = meshZw.ELEMENTS.INC;
        structure.ID = meshZw.ELEMENTS.ID;
        structure.IDPLOT = meshZw.ELEMENTS.IDPLOT;
        structure.IEN = meshZw.ELEMENTS.IEN;
        structure.EDOF = obj(i).EDOF;
        structure.QREF = obj(i).QREF;
        structure.calcDelta = false;
        structure.calcNablaGradient = false;
        
        structure.WEIGHTS = obj(i).WEIGHTS;
        if strcmpi(obj(i).Ansatzfunction,'Lagrange')
            structure.Element = size(obj(i).EDOF{1},2);
            numEle = numel(obj(i).EDOF);
        elseif strcmp(obj(i).Ansatzfunction,'HNURBS')
            structure.Element = [];
            structure.HCONTROL = obj(i).HCONTROL;
            structure.levelVec = obj(i).levelVec;
            structure.levelStruct = obj(i).levelStruct;
            structure.NODESDOF = obj(i).NODESDOF;
            [EDOF,IDSTRUCT,IDTOPOLOGIE,HCONTROL] = reCalcArrays2(structure);
            obj(i).IDTOPOLOGIE = IDTOPOLOGIE;
            obj(i).HCONTROL = HCONTROL;
            structure.HCONTROL = obj(i).HCONTROL;
            obj(i).EDOF = EDOF;
            obj(i).IDSTRUCT = IDSTRUCT;
            structure.EDOF = obj(i).EDOF;
            structure.IDSTRUCT = obj(i).IDSTRUCT;
            structure.IDTOPOLOGIE = obj(i).IDTOPOLOGIE;
            structure.KNOTS = meshZw.KNOTS;
            structure.STMatrix = obj(i).STMatrix;
            numEle = numel(obj(i).EDOF);
        end
        obj(i).SHAPEF = shapeFunctions(structure);
        if strcmpi(obj(i).Ansatzfunction,'Lagrange')
            structure = struct('DIM',[],'NGP',[],'Element',[],'Ansatz',obj(i).Ansatzfunction);
            structure.DIM = obj(i).DIM;
            structure.NGP = 1;
            structure.Element = size(obj(i).EDOF{1},2);
            zw = shapeFunctions(structure);
            obj(i).SHAPEF.dNr0 = zw.dNr;
        end
        
        %% HR/ HW element
        % Interpolationorder and Shapefunctions
        if strcmpi(obj(i).DESIGNATOR,'mmhwPF') || strcmpi(obj(i).DESIGNATOR,'mmhwSC') || strcmpi(obj(i).DESIGNATOR,'mmhwCascadeSC') || strcmpi(obj(i).DESIGNATOR,'mmhwClassicSC')
            if obj.DIM == 3
                if size(obj.EDOF{1},2) == 4 || size(obj.EDOF{1},2) == 10   % Tetrahedron
                    % Deformation Gradient F
                    if obj.ORDER_MIXED_VAR(1)==2
                        Nodes_F = 10;
                    elseif obj.ORDER_MIXED_VAR(1)==1
                        Nodes_F = 4;
                    elseif obj.ORDER_MIXED_VAR(1)==0
                        Nodes_F = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_F;
                    obj(i).SHAPEF_F = shapeFunctions(structure);
                    % Co-Factor H
                    if obj.ORDER_MIXED_VAR(2)==2
                        Nodes_Cof = 10;
                    elseif obj.ORDER_MIXED_VAR(2)==1
                        Nodes_Cof = 4;
                    elseif obj.ORDER_MIXED_VAR(2)==0
                        Nodes_Cof = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_Cof;
                    obj(i).SHAPEF_Cof = shapeFunctions(structure);
                    % Determinant J
                    if obj.ORDER_MIXED_VAR(3)==2
                        Nodes_Det = 10;
                    elseif obj.ORDER_MIXED_VAR(3)==1
                        Nodes_Det = 4;
                    elseif obj.ORDER_MIXED_VAR(3)==0
                        Nodes_Det = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_Det;
                    obj(i).SHAPEF_Det = shapeFunctions(structure);
                    % Electrical displacement (material) D
                    if obj.ORDER_MIXED_VAR_EL(1)==2
                        Nodes_D = 10;
                    elseif obj.ORDER_MIXED_VAR_EL(1)==1
                        Nodes_D = 4;
                    elseif obj.ORDER_MIXED_VAR_EL(1)==0
                        Nodes_D = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_D;
                    obj(i).SHAPEF_D = shapeFunctions(structure);
                    % Electrical displacement (spatial) D
                    if obj.ORDER_MIXED_VAR_EL(2)==2
                        Nodes_d = 10;
                    elseif obj.ORDER_MIXED_VAR_EL(2)==1
                        Nodes_d = 4;
                    elseif obj.ORDER_MIXED_VAR_EL(2)==0
                        Nodes_d = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_d;
                    obj(i).SHAPEF_d = shapeFunctions(structure);
                    
                elseif size(obj.EDOF{1},2) == 8 || size(obj.EDOF{1},2) == 20  % Brick ( linear or serendepity)
                    % Deformation Gradient F
                    if obj.ORDER_MIXED_VAR(1)==2
                        Nodes_F = 20;
                    elseif obj.ORDER_MIXED_VAR(1)==1
                        Nodes_F = 8;
                    elseif obj.ORDER_MIXED_VAR(1)==0
                        Nodes_F = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_F;
                    obj(i).SHAPEF_F = shapeFunctions(structure);
                    %Co-Factor H
                    if obj.ORDER_MIXED_VAR(2)==2
                        Nodes_Cof = 20;
                    elseif obj.ORDER_MIXED_VAR(2)==1
                        Nodes_Cof = 8;
                    elseif obj.ORDER_MIXED_VAR(2)==0
                        Nodes_Cof = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_Cof;
                    obj(i).SHAPEF_Cof = shapeFunctions(structure);
                    % Determinant J
                    if obj.ORDER_MIXED_VAR(3)==2
                        Nodes_Det = 20;
                    elseif obj.ORDER_MIXED_VAR(3)==1
                        Nodes_Det = 8;
                    elseif obj.ORDER_MIXED_VAR(3)==0
                        Nodes_Det = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_Det;
                    obj(i).SHAPEF_Det = shapeFunctions(structure);
                    % Electrical displacement (material) D
                    if obj.ORDER_MIXED_VAR_EL(1)==2
                        Nodes_D = 20;
                    elseif obj.ORDER_MIXED_VAR_EL(1)==1
                        Nodes_D = 8;
                    elseif obj.ORDER_MIXED_VAR_EL(1)==0
                        Nodes_D = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_D;
                    obj(i).SHAPEF_D = shapeFunctions(structure);
                    % Electrical displacement (spatial) d
                    if obj.ORDER_MIXED_VAR_EL(2)==2
                        Nodes_d = 20;
                    elseif obj.ORDER_MIXED_VAR_EL(2)==1
                        Nodes_d = 8;
                    elseif obj.ORDER_MIXED_VAR_EL(2)==0
                        Nodes_d = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_d;
                    obj(i).SHAPEF_d = shapeFunctions(structure);
                elseif size(obj.EDOF{1},2) == 27  % Quadratic Brick
                    % Deformation Gradient F
                    if obj.ORDER_MIXED_VAR(1)==2
                        Nodes_F = 27;
                    elseif obj.ORDER_MIXED_VAR(1)==1
                        Nodes_F = 8;
                    elseif obj.ORDER_MIXED_VAR(1)==0
                        Nodes_F = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_F;
                    obj(i).SHAPEF_F = shapeFunctions(structure);
                    % Co-Factor H
                    if obj.ORDER_MIXED_VAR(2)==2
                        Nodes_Cof = 27;
                    elseif obj.ORDER_MIXED_VAR(2)==1
                        Nodes_Cof = 8;
                    elseif obj.ORDER_MIXED_VAR(2)==0
                        Nodes_Cof = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_Cof;
                    obj(i).SHAPEF_Cof = shapeFunctions(structure);
                    % Determinant J
                    if obj.ORDER_MIXED_VAR(3)==2
                        Nodes_Det = 27;
                    elseif obj.ORDER_MIXED_VAR(3)==1
                        Nodes_Det = 8;
                    elseif obj.ORDER_MIXED_VAR(3)==0
                        Nodes_Det = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_Det;
                    obj(i).SHAPEF_Det = shapeFunctions(structure);
                    % Electrical displacement (material) D
                    if obj.ORDER_MIXED_VAR_EL(1)==2
                        Nodes_D = 27;
                    elseif obj.ORDER_MIXED_VAR_EL(1)==1
                        Nodes_D = 8;
                    elseif obj.ORDER_MIXED_VAR_EL(1)==0
                        Nodes_D = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_D;
                    obj(i).SHAPEF_D = shapeFunctions(structure);
                    % Electrical displacement (spatial) d
                    if obj.ORDER_MIXED_VAR_EL(2)==2
                        Nodes_d = 27;
                    elseif obj.ORDER_MIXED_VAR_EL(2)==1
                        Nodes_d = 8;
                    elseif obj.ORDER_MIXED_VAR_EL(2)==0
                        Nodes_d = 1;
                    end
                    structure = struct('DIM',obj(i).DIM,'NGP',obj(i).NGP,'Element',[],'HCONTROL',[],'Ansatz','Lagrange','nders',0);
                    structure.Element = Nodes_d;
                    obj(i).SHAPEF_d = shapeFunctions(structure);
                else
                    error('Element not implemented yet')
                end
            else
                error('Dimension not implemented yet')
            end
        end
        % Set initial values for mixed fields
%         if strcmpi(obj(i).DESIGNATOR,'mmhwPF') && obj(i).DIM == 3
%             if strcmpi(obj(i).MAT.name,'mooneyrivlinelectro')
%                 obj(i).SIGMA_FN  = repmat([ones(numEle,3)*(2*obj(i).MAT.mu1+12*obj(i).MAT.mue), zeros(numEle,6)], [1,Nodes_F]);
%                 obj(i).SIGMA_FN1 = repmat([ones(numEle,3)*(2*obj(i).MAT.mu1+12*obj(i).MAT.mue), zeros(numEle,6)], [1,Nodes_F]);
%                 obj(i).SIGMA_HN  = repmat([ones(numEle,3)*2*(obj(i).MAT.mu2), zeros(numEle,6)], [1,Nodes_Cof]);
%                 obj(i).SIGMA_HN1 = repmat([ones(numEle,3)*2*(obj(i).MAT.mu2), zeros(numEle,6)], [1,Nodes_Cof]);
%                 obj(i).SIGMA_JN  = ones(numEle,Nodes_Det)*(-2*(obj(i).MAT.mu1+2*obj(i).MAT.mu2+6*obj(i).MAT.mue));
%                 obj(i).SIGMA_JN1 = ones(numEle,Nodes_Det)*(-2*(obj(i).MAT.mu1+2*obj(i).MAT.mu2+6*obj(i).MAT.mue));
%                 obj(i).SIGMA_dN  = zeros(numEle,3*Nodes_d);
%                 obj(i).SIGMA_dN1 = zeros(numEle,3*Nodes_d);
%                 obj(i).FN  = repmat([ones(numEle,3), zeros(numEle,6)], [1,Nodes_F]);
%                 obj(i).FN1 = repmat([ones(numEle,3), zeros(numEle,6)], [1,Nodes_F]);
%                 obj(i).HN  = repmat([ones(numEle,3), zeros(numEle,6)], [1,Nodes_Cof]);
%                 obj(i).HN1 = repmat([ones(numEle,3), zeros(numEle,6)], [1,Nodes_Cof]);
%                 obj(i).JNe  = ones(numEle,Nodes_Det);
%                 obj(i).JN1e = ones(numEle,Nodes_Det);
%                 obj(i).DN  =  zeros(numEle,3*Nodes_D);
%                 obj(i).DN1  = zeros(numEle,3*Nodes_D);
%                 obj(i).dN  =  zeros(numEle,3*Nodes_d);
%                 obj(i).dN1  = zeros(numEle,3*Nodes_d);
%             elseif strcmpi(obj(i).MAT.name,'neoHook')
%                 obj(i).SIGMA_FN  = repmat([ones(numEle,3)*2*(obj(i).MAT.c1), zeros(numEle,6)], [1,Nodes_F]);
%                 obj(i).SIGMA_FN1 = repmat([ones(numEle,3)*2*(obj(i).MAT.c1), zeros(numEle,6)], [1,Nodes_F]);
%                 obj(i).SIGMA_JN  = ones(numEle,Nodes_Det)*(-2*(obj(i).MAT.c1+2*obj(i).MAT.c2));
%                 obj(i).SIGMA_JN1 = ones(numEle,Nodes_Det)*(-2*(obj(i).MAT.c1+2*obj(i).MAT.c2));
%                 obj(i).FN  = repmat([ones(numEle,3), zeros(numEle,6)], [1,Nodes_F]);
%                 obj(i).FN1 = repmat([ones(numEle,3), zeros(numEle,6)], [1,Nodes_F]);
%                 obj(i).JNe  = ones(numEle,Nodes_Det);
%                 obj(i).JN1e = ones(numEle,Nodes_Det);
%                 obj(i).DN  = 4.8*ones(numEle,Nodes_D);
%                 obj(i).DN1  = 2.3*ones(numEle,Nodes_D);
%                 obj(i).dN  = 1.2*ones(numEle,Nodes_d);
%                 obj(i).dN1  = 3.7*ones(numEle,Nodes_d);
%             else
%                 error('material model not implementet yet')
%             end
%         elseif strcmpi(obj(i).DESIGNATOR,'mmhwClassicSC') && obj(i).DIM == 3
%             if strcmpi(obj(i).MAT.name,'mooneyrivlinelectro')
%                 obj(i).SIGMA_FN  = repmat([ones(numEle,3)*(obj(i).MAT.mu1/2), zeros(numEle,3)], [1,Nodes_F]);
%                 obj(i).SIGMA_FN1 = repmat([ones(numEle,3)*(obj(i).MAT.mu1/2), zeros(numEle,3)], [1,Nodes_F]);
%                 obj(i).SIGMA_HN  = repmat([ones(numEle,3)*(obj(i).MAT.mu2/2), zeros(numEle,3)], [1,Nodes_Cof]);
%                 obj(i).SIGMA_HN1 = repmat([ones(numEle,3)*(obj(i).MAT.mu2/2), zeros(numEle,3)], [1,Nodes_Cof]);
%                 obj(i).SIGMA_JN  = ones(numEle,Nodes_Det)*(-(obj(i).MAT.mu1+2*obj(i).MAT.mu2)/2);
%                 obj(i).SIGMA_JN1 = ones(numEle,Nodes_Det)*(-(obj(i).MAT.mu1+2*obj(i).MAT.mu2)/2);
%                 
%                 obj(i).FN  = repmat([ones(numEle,3), zeros(numEle,3)], [1,Nodes_F]);
%                 obj(i).FN1 = repmat([ones(numEle,3), zeros(numEle,3)], [1,Nodes_F]);
%                 obj(i).HN  = repmat([ones(numEle,3), zeros(numEle,3)], [1,Nodes_Cof]);
%                 obj(i).HN1 = repmat([ones(numEle,3), zeros(numEle,3)], [1,Nodes_Cof]);
%                 obj(i).JNe  = ones(numEle,Nodes_Det);
%                 obj(i).JN1e = ones(numEle,Nodes_Det);
%                 obj(i).DN  =  ones(numEle,3*Nodes_D)*0;
%                 obj(i).DN1  = ones(numEle,3*Nodes_D)*0;
%             else
%                 error('material model not implementet yet')
%             end
        if strcmpi(obj(i).DESIGNATOR,'mmhwCascadeSC') && obj(i).DIM == 3
            if strcmpi(obj(i).MAT.name,'mooneyrivlinelectrothermo')
                obj(i).SIGMA_FN  = repmat(zeros(numEle,6), [1,Nodes_F]);
                obj(i).SIGMA_FN1 = repmat(zeros(numEle,6), [1,Nodes_F]);
                obj(i).SIGMA_HN  = repmat([ones(numEle,3)*0.5*(obj(i).MAT.mu2)-1/6*(obj(i).MAT.mu1+2*obj(i).MAT.mu2), zeros(numEle,3)], [1,Nodes_Cof]);
                obj(i).SIGMA_HN1 = repmat([ones(numEle,3)*0.5*(obj(i).MAT.mu2)-1/6*(obj(i).MAT.mu1+2*obj(i).MAT.mu2), zeros(numEle,3)], [1,Nodes_Cof]);
                obj(i).SIGMA_JN  = ones(numEle,Nodes_Det)*0.5*(-(obj(i).MAT.mu1+2*obj(i).MAT.mu2));
                obj(i).SIGMA_JN1 = ones(numEle,Nodes_Det)*0.5*(-(obj(i).MAT.mu1+2*obj(i).MAT.mu2));
                
                obj(i).FN  = repmat([ones(numEle,3), zeros(numEle,3)], [1,Nodes_F]);
                obj(i).FN1 = repmat([ones(numEle,3), zeros(numEle,3)], [1,Nodes_F]);
                obj(i).HN  = repmat([ones(numEle,3), zeros(numEle,3)], [1,Nodes_Cof]);
                obj(i).HN1 = repmat([ones(numEle,3), zeros(numEle,3)], [1,Nodes_Cof]);
                obj(i).JNe  = ones(numEle,Nodes_Det);
                obj(i).JN1e = ones(numEle,Nodes_Det);
                obj(i).DN  =  ones(numEle,3*Nodes_D)*0;
                obj(i).DN1  = ones(numEle,3*Nodes_D)*0;
%             elseif strcmpi(obj(i).MAT.name,'mooneyrivlinelectrothermoTI')
%                 obj(i).SIGMA_FN  = repmat(zeros(numEle,6), [1,Nodes_F]);
%                 obj(i).SIGMA_FN1 = repmat(zeros(numEle,6), [1,Nodes_F]);
%                 obj(i).SIGMA_HN  = repmat([ones(numEle,3)*0.5*(obj(i).MAT.mu2)-1/6*(obj(i).MAT.mu1+2*obj(i).MAT.mu2), zeros(numEle,3)], [1,Nodes_Cof]);
%                 obj(i).SIGMA_HN1 = repmat([ones(numEle,3)*0.5*(obj(i).MAT.mu2)-1/6*(obj(i).MAT.mu1+2*obj(i).MAT.mu2), zeros(numEle,3)], [1,Nodes_Cof]);
%                 obj(i).SIGMA_JN  = ones(numEle,Nodes_Det)*0.5*(-(obj(i).MAT.mu1+2*obj(i).MAT.mu2));
%                 obj(i).SIGMA_JN1 = ones(numEle,Nodes_Det)*0.5*(-(obj(i).MAT.mu1+2*obj(i).MAT.mu2));
%                 
%                 obj(i).FN  = repmat([ones(numEle,3), zeros(numEle,3)], [1,Nodes_F]);
%                 obj(i).FN1 = repmat([ones(numEle,3), zeros(numEle,3)], [1,Nodes_F]);
%                 obj(i).HN  = repmat([ones(numEle,3), zeros(numEle,3)], [1,Nodes_Cof]);
%                 obj(i).HN1 = repmat([ones(numEle,3), zeros(numEle,3)], [1,Nodes_Cof]);
%                 obj(i).JNe  = ones(numEle,Nodes_Det);
%                 obj(i).JN1e = ones(numEle,Nodes_Det);
%                 obj(i).DN  =  ones(numEle,3*Nodes_D)*0;
%                 obj(i).DN1  = ones(numEle,3*Nodes_D)*0;
            else
                error('material model not implementet yet')
            end
        elseif strcmpi(obj(i).DESIGNATOR,'mmhwSC') && obj(i).DIM == 3
            if strcmpi(obj(i).MAT.name,'mooneyrivlinelectrothermo')
                if isempty(meshZw.DN)
                    obj(i).DN  =  ones(numEle,3*Nodes_D)*0;
                    obj(i).DN1  = ones(numEle,3*Nodes_D)*0;
                else
                    obj(i).DN  =  meshZw.DN;
                    obj(i).DN1  = meshZw.DN1;
                end
%             elseif strcmpi(obj(i).MAT.name,'mooneyrivlinelectrothermoTI')
%                 obj(i).DN  =  ones(numEle,3*Nodes_D)*0;
%                 obj(i).DN1  = ones(numEle,3*Nodes_D)*0;
            else
                error('material model not implementet yet')
            end
        end
    end
end

