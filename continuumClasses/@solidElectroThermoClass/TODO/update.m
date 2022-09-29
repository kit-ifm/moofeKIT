function obj = update(obj,varargin)
    %% Updates specified field
    %
    % Syntax
    %
    % update(obj,'PropertyName',PropertyValue)
    %
    % 09.01.2012 C.HESCH
    
    %% Check input
    newton = false;
    globalinput = false;
    field = '';
    input = [];
    inputglobal = [];
    for i = 1:size(varargin,2)
        if strcmp(varargin{i},'field')
            field = varargin{i+1};
        elseif strcmp(varargin{i},'input')
            input = varargin{i+1};
        elseif strcmp(varargin{i},'inputglobal')
            globalinput = true;
            inputglobal = varargin{i+1};
        elseif strcmp(varargin{i},'newton')
            newton = true;
        end
    end
    if (isempty(field) || (isempty(input) && isempty(inputglobal))) && ~ newton
        return
    end

%% Update
for i = 1:numel(obj)
    numElements = numel(obj(i).EDOF);
    edof = obj.GLOBALEDOF;
    if globalinput
        input = zeros(size(inputglobal));
        input(obj(i).NODESDOF) = - inputglobal(obj(i).NODESDOF) + obj(i).QN1;
    end
    if newton
        obj(i).QN1 = obj(i).QN1 - input(obj(i).NODESDOF);
        if  strcmpi(obj(i).DESIGNATOR,'mmhwPF')
            SIGMA_FN1   =  obj(i).SIGMA_FN1;
            SIGMA_HN1   =  obj(i).SIGMA_HN1;
            SIGMA_JN1   =  obj(i).SIGMA_JN1;
            SIGMA_dN1   =  obj(i).SIGMA_dN1;
            FN1         =  obj(i).FN1;
            HN1         =  obj(i).HN1;
            JN1e        =  obj(i).JN1e;
            DN1         =  obj(i).DN1;
            dN1         =  obj(i).dN1;
            StaConStruct = obj(i).StaConStruct;
            numElements = numel(obj(i).EDOF);
            edof = obj.GLOBALEDOF;
            YN1=zeros(numElements,size([StaConStruct.RY],1));
            SIGMA_YN1=zeros(numElements,size([StaConStruct.RSIGY],1));
            for j = 1:numElements
                edofLocal = edof{j};
                dofsLocalm=[1:4:size(edof{j},2) ; 2:4:size(edof{j},2) ; 3:4:size(edof{j},2)];
                dofsLocale=4:4:size(edof{j},2);
                edofLocalm = edofLocal(dofsLocalm(:));
                edofLocale = edofLocal(dofsLocale(:));
                YN1(j,:)   = [FN1(j,:), HN1(j,:), JN1e(j,:), dN1(j,:)]';
                SIGMA_YN1(j,:) = [SIGMA_FN1(j,:), SIGMA_HN1(j,:), SIGMA_JN1(j,:) SIGMA_dN1(j,:)]';
                
                SIGMA_YN1(j,:)  = SIGMA_YN1(j,:) + (StaConStruct(j).RbSIGY - StaConStruct(j).MSIGYX*input(edofLocalm) - StaConStruct(j).MSIGYPHI*input(edofLocale))';
                DN1(j,:)        = DN1(j,:) + (StaConStruct(j).RbD + StaConStruct(j).MbDX*input(edofLocalm) + StaConStruct(j).MbDPHI*input(edofLocale))';
                YN1(j,:)        = YN1(j,:) + (StaConStruct(j).RbY + StaConStruct(j).MbYD - StaConStruct(j).MbYX*input(edofLocalm) - StaConStruct(j).MbYPHI*input(edofLocale))';
                
                SIGMA_FN1(j,:) = SIGMA_YN1(j,1:size(SIGMA_FN1,2));
                SIGMA_HN1(j,:) = SIGMA_YN1(j,size(SIGMA_FN1,2)+(1:size(SIGMA_HN1,2)));
                SIGMA_JN1(j,:) = SIGMA_YN1(j,size(SIGMA_FN1,2)+size(SIGMA_HN1,2)+(1:size(SIGMA_JN1,2)));
                SIGMA_dN1(j,:) = SIGMA_YN1(j,size(SIGMA_FN1,2)+size(SIGMA_HN1,2)+size(SIGMA_JN1,2)+(1:size(SIGMA_dN1,2)));
                
                FN1(j,:) = YN1(j,1:size(SIGMA_FN1,2));
                HN1(j,:) = YN1(j,size(FN1,2)+(1:size(HN1,2)));
                JN1e(j,:) = YN1(j,size(FN1,2)+size(HN1,2)+(1:size(JN1e,2)));
                dN1(j,:) = YN1(j,size(FN1,2)+size(HN1,2)+size(JN1e,2)+(1:size(dN1,2)));
                
            end
            obj(i).SIGMA_FN1 = SIGMA_FN1;
            obj(i).SIGMA_HN1 = SIGMA_HN1;
            obj(i).SIGMA_JN1 = SIGMA_JN1;
            obj(i).SIGMA_dN1 = SIGMA_dN1;
            obj(i).FN1      = FN1;
            obj(i).HN1      = HN1;
            obj(i).JN1e     = JN1e;
            obj(i).dN1     = dN1;
            obj(i).DN1     = DN1;
        elseif  strcmpi(obj(i).DESIGNATOR,'mmhwSC')
            DN1         =  obj(i).DN1;
            StaConStruct = obj(i).StaConStruct;
            numElements = numel(obj(i).EDOF);
            edof = obj.GLOBALEDOF;
            for j = 1:numElements
                edofLocal = edof{j};
                dofsLocalm=[1:5:size(edof{j},2) ; 2:5:size(edof{j},2) ; 3:5:size(edof{j},2)];
                dofsLocale=4:5:size(edof{j},2);
                dofsLocalt=5:5:size(edof{j},2);
                edofLocalm = edofLocal(dofsLocalm(:));
                edofLocale = edofLocal(dofsLocale(:));
                edofLocalt = edofLocal(dofsLocalt(:));
%                 KDX*DX+KDP*DP+KDD*DD+KDT*DT = -RD
%                 DD = inv(KDD)*(-RD - KDX*DX - KDP*DP - KDT*DT)
                DN1(j,:)  = DN1(j,:) + ( inv(StaConStruct(j).KDD)*( -StaConStruct(j).RD + StaConStruct(j).KDX*input(edofLocalm) + StaConStruct(j).KDP*input(edofLocale) + StaConStruct(j).KDT*input(edofLocalt)) )';
            end
            obj(i).DN1     = DN1;
        elseif  strcmpi(obj(i).DESIGNATOR,'mmhwClassicSC') ||  strcmpi(obj(i).DESIGNATOR,'mmhwCascadeSC')
            SIGMA_FN1   =  obj(i).SIGMA_FN1;
            SIGMA_HN1   =  obj(i).SIGMA_HN1;
            SIGMA_JN1   =  obj(i).SIGMA_JN1;
            FN1         =  obj(i).FN1;
            HN1         =  obj(i).HN1;
            JN1e        =  obj(i).JN1e;
            DN1         =  obj(i).DN1;
            StaConStruct = obj(i).StaConStruct;
            numElements = numel(obj(i).EDOF);
            edof = obj.GLOBALEDOF;
            
            for j = 1:numElements
                edofLocal = edof{j};
                dofsLocalm=[1:4:size(edof{j},2) ; 2:4:size(edof{j},2) ; 3:4:size(edof{j},2)];
                dofsLocale=4:4:size(edof{j},2);
                edofLocalm = edofLocal(dofsLocalm(:));
                edofLocale = edofLocal(dofsLocale(:));
                                    
                DeltaEPSILONN1  = inv(StaConStruct(j).KSIGEPS)*(-StaConStruct(j).RSIG + StaConStruct(j).KSIGX*input(edofLocalm));
                DeltaDN1        = inv(StaConStruct(j).KDD)*(-StaConStruct(j).RD + StaConStruct(j).KDP*input(edofLocale)-StaConStruct(j).KDEPS*DeltaEPSILONN1);
                DeltaSIGMAN1    = inv(StaConStruct(j).KEPSSIG)*(-StaConStruct(j).REPS-StaConStruct(j).KEPSEPS*DeltaEPSILONN1-StaConStruct(j).KEPSD*DeltaDN1);
                              
                SIGMAN1(j,:)    = [SIGMA_FN1(j,:), SIGMA_HN1(j,:), SIGMA_JN1(j,:)]'  + DeltaSIGMAN1;
                EPSILONN1(j,:)  = [FN1(j,:), HN1(j,:), JN1e(j,:)]' + DeltaEPSILONN1;
                DN1(j,:)        = DN1(j,:)  + DeltaDN1';
                
                SIGMA_FN1(j,:) = SIGMAN1(j,1:size(SIGMA_FN1,2));
                SIGMA_HN1(j,:) = SIGMAN1(j,size(SIGMA_FN1,2)+(1:size(SIGMA_HN1,2)));
                SIGMA_JN1(j,:) = SIGMAN1(j,(size(FN1,2)+size(SIGMA_HN1,2))+1:end);
                
                FN1(j,:)  = EPSILONN1(j,1:size(FN1,2));
                HN1(j,:)  = EPSILONN1(j,size(FN1,2)+(1:size(HN1,2)));
                JN1e(j,:) = EPSILONN1(j,(size(FN1,2)+ size(HN1,2))+1:end);
            end
            obj(i).SIGMA_FN1 = SIGMA_FN1;
            obj(i).SIGMA_HN1 = SIGMA_HN1;
            obj(i).SIGMA_JN1 = SIGMA_JN1;
            obj(i).FN1       = FN1;
            obj(i).HN1       = HN1;
            obj(i).JN1e      = JN1e;
            obj(i).DN1       = DN1;
        end
    else
        if ~isproperty(obj(i),field)
            error('Specified field does not exist.');
        end
        obj(i).(field) = input(obj(i).NODESDOF);
    end
end
end
