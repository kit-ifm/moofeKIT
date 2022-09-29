function obj = update(obj,varargin)
%% Updates specific filed
%
% Syntax
% out = update(obj,'PropertyName',PropertyValue)
%
%
% Created: 			Mi, 10 Aug 2016
% Responsibilty: 	Mark Schiebl
% Editors: 
% gk534  gk604  gk612 
% Description: 
%----------------------------------------------------------------------
%
%
%

%% Check input
newton = false;
field = '';
input = [];
for i = 1:size(varargin,2)
    if strcmp(varargin{i},'field')
        field = varargin{i+1};
    elseif strcmp(varargin{i},'input')
        input = varargin{i+1};
    elseif strcmp(varargin{i},'newton')
        newton = true;
    end
end
if (isempty(field) || isempty(input)) && ~ newton
    return
end

%% Update
for i = 1:numel(obj)
    numElements = numel(obj(i).EDOF);
    edof = obj.GLOBALEDOF;
    if newton
        obj(i).QN1 = obj(i).QN1 - input(obj(i).NODESDOF);
        if strcmpi(obj(i).DESIGNATOR,'enheas') || strcmpi(obj(i).DESIGNATOR,'mmeas') || strcmpi(obj(i).DESIGNATOR,'mmeasimpr')
            ALPHAN1 = obj(i).ALPHAN1;
            ENHStruct = obj(i).ENHStruct;
            for j = 1:numElements
                edofLocal = edof{j};
                ALPHAN1(j,:) = ALPHAN1(j,:)-(ENHStruct(j).RSe-ENHStruct(j).KSe*input(edofLocal))';
            end
            obj(i).ALPHAN1 = ALPHAN1;
        end
        if  strcmpi(obj(i).DESIGNATOR,'mmhr') || strcmpi(obj(i).DESIGNATOR,'mmhrSplit')
            edof = obj(i).GLOBALEDOF;
            numElements = numel(obj(i).EDOF);
            SIGMA_FN1 =  obj(i).SIGMA_FN1;
            SIGMA_HN1 =  obj(i).SIGMA_HN1;
            SIGMA_JN1 =  obj(i).SIGMA_JN1;
            StaConStruct = obj(i).StaConStruct;
            for j = 1:numElements
                edofLocal = edof{j};
                SIGMAN1(j,:)   = [SIGMA_FN1(j,:), SIGMA_HN1(j,:), SIGMA_JN1(j,:)]';
                SIGMAN1(j,:)   = SIGMAN1(j,:)' + StaConStruct(j).KSe*input(edofLocal) - StaConStruct(j).RSe;
                SIGMA_FN1(j,:) = SIGMAN1(j,1:size(SIGMA_FN1,2));
                SIGMA_HN1(j,:) = SIGMAN1(j,size(SIGMA_FN1,2)+(1:size(SIGMA_HN1,2)));
                SIGMA_JN1(j,:) = SIGMAN1(j,end);
            end
            obj(i).SIGMA_FN1 = SIGMA_FN1;
            obj(i).SIGMA_HN1 = SIGMA_HN1;
            obj(i).SIGMA_JN1 = SIGMA_JN1;
        elseif (strcmpi(obj(i).DESIGNATOR,'mmhwCascadeSC') || strcmpi(obj(i).DESIGNATOR,'mmhwGenericCascadeSC'))
            SIGMA_FN1   =  obj(i).SIGMA_FN1;
            SIGMA_HN1   =  obj(i).SIGMA_HN1;
            SIGMA_JN1   =  obj(i).SIGMA_JN1;
            FN1         =  obj(i).FN1;
            HN1         =  obj(i).HN1;
            JN1e        =  obj(i).JN1e;
            StaConStruct = obj(i).StaConStruct;
            for j = 1:numElements
                edofLocal = edof{j};
                dofsLocalm=[1:4:size(edof{j},2) ; 2:4:size(edof{j},2) ; 3:4:size(edof{j},2)];
                dofsLocalt=4:4:size(edof{j},2);
                edofLocalm = edofLocal(dofsLocalm(:));
                edofLocalt = edofLocal(dofsLocalt(:));
                SIGMAN1(j,:)   = [SIGMA_FN1(j,:), SIGMA_HN1(j,:), SIGMA_JN1(j,:)]';
                EPSILONN1(j,:) = [FN1(j,:), HN1(j,:), JN1e(j,:)]';
                
                EPSILONN1(j,:)   = EPSILONN1(j,:) - (StaConStruct(j).RCon2 - StaConStruct(j).KCon2*input(edofLocalm))';
                SIGMAN1(j,:)     = SIGMAN1(j,:)   - (StaConStruct(j).RCon1 - StaConStruct(j).KCon1*StaConStruct(j).RCon2 +  StaConStruct(j).KCon1*StaConStruct(j).KCon2*input(edofLocalm) - StaConStruct(j).KEPSSIG\StaConStruct(j).KEPST*input(edofLocalt))';
                
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
            obj(i).FN1      = FN1;
            obj(i).HN1      = HN1;
            obj(i).JN1e     = JN1e;            
        end
    else
        if ~isproperty(obj(i),field)
            error('Specified field does not exist.');
        end
        obj(i).(field) = input(obj(i).NODESDOF);
    end
end
end
