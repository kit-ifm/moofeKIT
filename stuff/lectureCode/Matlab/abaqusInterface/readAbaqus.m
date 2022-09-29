function out = readAbaqus(obj,filename,varargin)
    % Import function for ABAQUS files (.inp).
    %
    % Syntax
    %
    % readAbaqus(obj,'filename')
    %
    % Description
    %
    % Reads the data from a specified Abaqus file and stores it in the nm.mesh
    % object.
    %
    % 05.04.2010 C.HESCH
    %
    % Change 19.05.2011
    % case '*Node' changed, use of the first column to identify the node-number.
    % C.Hesch
    
    %% Checking input
    if ~isa(obj, 'pre.abaqus') || ~ischar(filename)
        error('Please provide a string as first input and a mesh object as second input')
    end
    datei = fopen(filename,'r');
    if datei == -1
        error('No file with this name.')
    end
    dim = 3;
    for i = 1:size(varargin,2)
        if strcmp(varargin{i},'dim')
            if ~isnumeric(varargin{i+1})
                error('Please provide a numeric as input.')
            end
            dim = varargin{i+1};
        end
    end
    
    %% Definition of variables
    out = pre.abaqus;
    part = struct('Name',[],'Element',[],'Geometrie',[],'Edof',[]);
    instance = struct('Name',[],'Part',[],'Translation',[],'Rotation',[]);
    shell = struct('Thickness',[],'IntegrationPoints',[],'ElementSet',[],'Material',[]);
    section = struct('ElementSet',[],'Material',[]);
    nset = struct('Name',[],'Instance',[],'Nodes',[],'Part',[]);
    elset = struct('Name',[],'Instance',[],'Elements',[],'Part',[]);
    surfaces = struct('Name',[],'Elset',[],'Type',[],'Side',[]);
    laeufer = 0;
    nsetLaeufer = 0;
    elsetLaeufer = 0;
    partLaeufer = 0;
    shellLaeufer = 0;
    sectionLaeufer = 0;
    instanceLaeufer = 0;
    surfLaeufer = 0;
    partSwitch = false;
    
    %% Processing file
    tline = fgets(datei);
    while ~feof(datei)
        if size(tline,2) > 0 && strcmp(tline(1),'*')
            aline = textscan(tline, '%s', 'delimiter', ',');
        else
            aline{1} = {'ABC'};
        end
        
        switch aline{1}{1}
            
            case '*Part'
                
                partLaeufer = partLaeufer + 1;
                aline = textscan(tline, '%s %s', 'delimiter', '=');
                part(partLaeufer).Name = aline{2}{1};
                partSwitch = true;
                laeufer = laeufer + 1;
                tline = fgets(datei);
                
            case '*End Part'
                
                partSwitch = false;
                laeufer = laeufer + 1;
                tline = fgets(datei);
                
            case '*Node'
                
                laeufer = laeufer + 1;
                tline = fgets(datei);
                startNode = laeufer;
                while ~strncmp(tline,'*',1) && ~feof(datei)
                    laeufer = laeufer + 1;
                    tline = fgets(datei);
                end
                a = textscan(tline, '%s','delimiter',',');
                zw = dlmread(filename,',',[startNode 0 laeufer-1 dim]);
                part(partLaeufer).Geometrie(zw(:,1),:) = zw(:,2:end);
                clear zw
                
            case '*Element'
                
                aline = textscan(tline, '*%s %s','delimiter','=');
                part(partLaeufer).Element = aline{2}{1};
                laeufer = laeufer + 1;
                tline = fgets(datei);
                startNode = laeufer;
                a = textscan(tline, '%s','delimiter',',');
                rightCorner = size(a{1},1)-1;
                while ~strncmp(tline,'*',1) && ~feof(datei)
                    laeufer = laeufer + 1;
                    tline = fgets(datei);
                end
                zw = dlmread(filename,',',[startNode 0 laeufer-1 rightCorner]);
                if strcmp(aline{2},'C3D20R')
                    part(partLaeufer).Edof = [zw(1:2:end,2:end) zw(2:2:end,1:5)];
                else
                    part(partLaeufer).Edof = zw(:,2:end);
                end
                
            case '*Surface'
                
                if strcmp(aline{1}{2},'type=ELEMENT')
                    surfLaeufer = surfLaeufer + 1;
                    aline = textscan(tline, '%s','delimiter',',');
                    bline = textscan(aline{1}{3}, '%s','delimiter','=');
                    surfaces(surfLaeufer).Name = bline{1}{2};
                    bline = textscan(aline{1}{2}, '%s','delimiter','=');
                    if strcmp(bline{1}(1),'type')
                        surfaces(surfLaeufer).Type = bline{1}{2};
                    end
                    laeufer = laeufer + 1;
                    tline = fgets(datei);
                    aline = textscan(tline, '%s','delimiter',',');
                    surfaces(surfLaeufer).Elset = aline{1}{1};
                    surfaces(surfLaeufer).Side = aline{1}{2};
                    laeufer = laeufer + 1;
                    tline = fgets(datei);
                    while ~strncmp(tline,'*',1) && ~feof(datei)
                        surfLaeufer = surfLaeufer + 1;
                        surfaces(surfLaeufer).Name = surfaces(surfLaeufer-1).Name;
                        surfaces(surfLaeufer).Type = surfaces(surfLaeufer-1).Type;
                        aline = textscan(tline, '%s','delimiter',',');
                        surfaces(surfLaeufer).Elset = aline{1}{1};
                        surfaces(surfLaeufer).Side = aline{1}{2};
                        laeufer = laeufer + 1;
                        tline = fgets(datei);
                    end
                else
                    laeufer = laeufer + 1;
                    tline = fgets(datei);
                end
                
            case '*Shell Section'
                
                shellLaeufer = shellLaeufer + 1;
                aline = textscan(tline, '%s','delimiter',',');
                bline = textscan(aline{1}{2}, '%s','delimiter','=');
                shell(shellLaeufer).ElementSet = bline{1}{2};
                bline = textscan(aline{1}{3}, '%s','delimiter','=');
                shell(shellLaeufer).Material = bline{1}{2};
                laeufer = laeufer + 1;
                tline = fgets(datei);
                cline = textscan(tline, '%f','delimiter',',');
                shell(shellLaeufer).Thickness = cline{1}(1);
                shell(shellLaeufer).IntegrationPoints = cline{1}(2);
                laeufer = laeufer + 1;
                tline = fgets(datei);
                
            case '*Solid Section'
                
                sectionLaeufer = sectionLaeufer + 1;
                aline = textscan(tline, '%s','delimiter',',');
                bline = textscan(aline{1}{2}, '%s','delimiter','=');
                section(sectionLaeufer).ElementSet = bline{1}{2};
                bline = textscan(aline{1}{3}, '%s','delimiter','=');
                section(sectionLaeufer).Material = bline{1}{2};
                laeufer = laeufer + 1;
                tline = fgets(datei);
                
            case '*Instance'
                
                instanceLaeufer = instanceLaeufer + 1;
                aline = textscan(tline, '%s','delimiter',',');
                bline = textscan(aline{1}{2}, '%s','delimiter','=');
                instance(instanceLaeufer).Name = bline{1}{2};
                bline = textscan(aline{1}{3}, '%s','delimiter','=');
                instance(instanceLaeufer).Part = bline{1}{2};
                laeufer = laeufer + 1;
                tline = fgets(datei);
                if strncmp(tline,'*',1)
                    if dim == 2
                        instance(instanceLaeufer).Translation = [0 0];
                        instance(instanceLaeufer).Rotation = [1 0 0 0 0];
                    elseif dim == 3
                        instance(instanceLaeufer).Translation = [0 0 0];
                        instance(instanceLaeufer).Rotation = [1 0 0 0 0 0 0];
                    end
                else
                    if dim == 2
                        instance(instanceLaeufer).Translation = [0 0];
                        instance(instanceLaeufer).Rotation = [1 0 0 0 0];
                        if ~strncmp(tline,'*',1)
                            instance(instanceLaeufer).Translation = sscanf(tline, '%f, %f');
                        end
                        laeufer = laeufer + 1;
                        tline = fgets(datei);
                        if ~strncmp(tline,'*',1)
                            instance(instanceLaeufer).Rotation = sscanf(tline, '%f, %f, %f, %f, %f');
                        end
                    elseif dim == 3
                        instance(instanceLaeufer).Translation = [0 0 0];
                        instance(instanceLaeufer).Rotation = [1 0 0 0 0 0 0];
                        if ~strncmp(tline,'*',1)
                            instance(instanceLaeufer).Translation = sscanf(tline, '%f, %f, %f');
                        end
                        laeufer = laeufer + 1;
                        tline = fgets(datei);
                        if ~strncmp(tline,'*',1)
                            instance(instanceLaeufer).Rotation = sscanf(tline, '%f, %f, %f, %f, %f, %f, %f');
                        end
                    end
                end
                
            case '*Nset'
                
                nsetLaeufer = nsetLaeufer + 1;
                if partSwitch
                    nset(nsetLaeufer).Part = part(partLaeufer).Name;
                end
                aline = textscan(tline, '%s','delimiter',',');
                bline = textscan(aline{1}{2}, '%s','delimiter','=');
                nset(nsetLaeufer).Name = bline{1}{2};
                for i = 2:size(aline{1},1)
                    bline = textscan(aline{1}{i}, '%s','delimiter','=');
                    if strncmp(bline{1}(1),'instance',8)
                        nset(nsetLaeufer).Instance = bline{1}{2};
                    end
                end
                if strncmp(aline{1}(end),'generate',8)
                    laeufer = laeufer + 1;
                    tline = fgets(datei);
                    cline = textscan(tline, '%f','delimiter',',');
                    nset(nsetLaeufer).Nodes = (cline{1}(1):cline{1}(3):cline{1}(2))';
                else
                    laeufer = laeufer + 1;
                    tline = fgets(datei);
                    startNode = laeufer;
                    a = textscan(tline, '%s','delimiter',',');
                    rightCorner = size(a{1},1)-1;
                    while ~strncmp(tline,'*',1) && ~feof(datei)
                        laeufer = laeufer + 1;
                        tline = fgets(datei);
                    end
                    nodesMatrix = dlmread(filename,',',[startNode 0 laeufer-1 rightCorner]);
                    zw = nodesMatrix(:);
                    logikVektor = zw > 0;
                    nset(nsetLaeufer).Nodes = zw(logikVektor);
                end
                
            case '*Elset'
                
                elsetLaeufer = elsetLaeufer + 1;
                if partSwitch
                    elset(elsetLaeufer).Part = part(partLaeufer).Name;
                end
                aline = textscan(tline, '%s','delimiter',',');
                bline = textscan(aline{1}{2}, '%s','delimiter','=');
                elset(elsetLaeufer).Name = bline{1}{2};
                for i = 2:size(aline{1},1)
                    bline = textscan(aline{1}{i}, '%s','delimiter','=');
                    if strncmp(bline{1}(1),'instance',8)
                        elset(elsetLaeufer).Instance = bline{1}{2};
                    end
                end
                if strncmp(aline{1}(end),'generate',8)
                    laeufer = laeufer + 1;
                    tline = fgets(datei);
                    cline = textscan(tline, '%f','delimiter',',');
                    elset(elsetLaeufer).Elements = (cline{1}(1):cline{1}(3):cline{1}(2))';
                else
                    laeufer = laeufer + 1;
                    tline = fgets(datei);
                    startNode = laeufer;
                    a = textscan(tline, '%s','delimiter',',');
                    rightCorner = size(a{1},1)-1;
                    while ~strncmp(tline,'*',1) && ~feof(datei)
                        laeufer = laeufer + 1;
                        tline = fgets(datei);
                    end
                    nodesMatrix = dlmread(filename,',',[startNode 0 laeufer-1 rightCorner]);
                    zw = nodesMatrix(:);
                    logikVektor = zw > 0;
                    elset(elsetLaeufer).Elements = zw(logikVektor);
                end
                
            otherwise
                
                laeufer = laeufer + 1;
                tline = fgets(datei);
        end
    end
    fclose(datei);
    
    %% Merge data
    for i = 1:size(instance,2)
        for j = 1:size(part,2)
            if strcmp(instance(i).Part,part(j).Name)
                out(i).NAME = part(j).Name;
                out(i).ELEMENT = part(j).Element;
                out(i).QNODE = part(j).Geometrie;
                out(i).EDOF = part(j).Edof;
            end
        end
        % Node-sets
        for j = 1:size(nset,2)
            if size(elset(j).Instance,2) > 0
                if strcmp(nset(j).Instance,instance(i).Name)
                    out(i).SUBSETS = [out(i).SUBSETS struct('NAME',nset(j).Name,'type','NODESET','EDOF',nset(j).Nodes)];
                end
            else
                if strcmp(elset(j).Part,instance(i).Part)
                    out(i).SUBSETS = [out(i).SUBSETS struct('NAME',nset(j).Name,'type','NODESET','EDOF',nset(j).Nodes)];
                end
            end
        end
        % Element-sets
        for j = 1:size(elset,2)
            if size(elset(j).Instance,2) > 0
                if strcmp(elset(j).Instance,instance(i).Name)
                    out(i).SUBSETS = [out(i).SUBSETS struct('NAME',elset(j).Name,'type','ELSET','EDOF',out(i).EDOF(elset(j).Elements,:))];
                end
            else
                if strcmp(elset(j).Part,instance(i).Part)
                    out(i).SUBSETS = [out(i).SUBSETS struct('NAME',elset(j).Name,'type','ELSET','EDOF',out(i).EDOF(elset(j).Elements,:))];
                end
            end
        end
        % Surfaces
        SX = [1 2 3 4;
            5 8 7 6;
            1 5 6 2;
            2 6 7 3;
            3 7 8 4;
            4 8 5 1];
        for j = 1:size(surfaces,2)
            for k = 1:size(elset,2)
                if strcmp(surfaces(j).Elset,elset(k).Name)
                    if strcmp(elset(k).Instance,instance(i).Name)
                        foundSurf = false;
                        for m = 1:size(out(i).SUBSETS,2)
                            if strcmp(out(i).SUBSETS(m).NAME,surfaces(j).Name)
                                foundSurf = true;
                                switch out(i).ELEMENT
                                    case 'C3D8R'
                                        zw = out(i).EDOF(elset(k).Elements,:);
                                        % out(i).SUBSETS(m).EDOF = [out(i).SUBSETS(m).EDOF; zw(:,SX(str2double(surfaces(j).Side(2)),:))];
                                        out(i).SUBSETS(m).EDOF = [zw(:,SX(str2double(surfaces(j).Side(2)),:))];
                                    case 'C3D4'
                                        SX4nodes = [1 3 2 ;
                                            1 2 4 ;
                                            2 3 4 ;
                                            1 4 3 ];
                                        zw = out(i).EDOF(elset(k).Elements,:);
                                        out(i).SUBSETS(m).EDOF = [ zw(:,SX4nodes(str2double(surfaces(j).Side(2)),:))];
                                    case 'C3D20R'
                                        SX20nodes = [   1 2 3 4 9 10 11 12;
                                            5 8 7 6 16 15 14 13;
                                            1 5 6 2 17 13 18 9;
                                            2 6 7 3 18 14 19 10;
                                            3 7 8 4 19 15 20 11;
                                            4 8 5 1 20 16 17 12];
                                        zw = out(i).EDOF(elset(k).Elements,:);
                                        out(i).SUBSETS(m).EDOF = [out(i).SUBSETS(m).EDOF; zw(:,SX20nodes(str2double(surfaces(j).Side(2)),:))];
                                    case 'C3D10'
                                        SX10nodes = [   1 3 2 7 6 5;
                                            1 2 4 5 9 8;
                                            2 4 3 9 10 6;
                                            1 4 3 8 10 7];
                                        zw = out(i).EDOF(elset(k).Elements,:);
                                        out(i).SUBSETS(m).EDOF = [out(i).SUBSETS(m).EDOF; zw(:,SX10nodes(str2double(surfaces(j).Side(2)),:))];
                                    case 'CPE4'
                                        out(i).SUBSETS(m).EDOF = [out(i).SUBSETS(m).EDOF; out(i).EDOF(elset(k).Elements,:)];
                                    case 'CPS4R'
                                        SX2nodes = [   1 2;
                                            2 3;
                                            3 4;
                                            4 1];
                                        zw = out(i).EDOF(elset(k).Elements,:);
                                        out(i).SUBSETS(m).EDOF = [out(i).SUBSETS(m).EDOF; zw(:,SX2nodes(str2double(surfaces(j).Side(2)),:))];
                                    otherwise
                                        error('Type not implemented')
                                end
                            end
                        end
                        if ~foundSurf
                            switch out(i).ELEMENT
                                case 'C3D8R'
                                    zw = out(i).EDOF(elset(k).Elements,:);
                                    edofzw = zw(:,SX(str2double(surfaces(j).Side(2)),:));
                                case 'C3D4'
                                    SX4nodes = [    1 3 2 ;
                                        1 2 4 ;
                                        2 4 3 ;
                                        1 4 3 ];
                                    zw = out(i).EDOF(elset(k).Elements,:);
                                    edofzw = zw(:,SX4nodes(str2double(surfaces(j).Side(2)),:));
                                case 'C3D20R'
                                    SX20nodes = [   1 2 3 4 9 10 11 12;
                                        5 8 7 6 16 15 14 13;
                                        1 5 6 2 17 13 18 9;
                                        2 6 7 3 18 14 19 10;
                                        3 7 8 4 19 15 20 11;
                                        4 8 5 1 20 16 17 12];
                                    zw = out(i).EDOF(elset(k).Elements,:);
                                    edofzw = zw(:,SX20nodes(str2double(surfaces(j).Side(2)),:));
                                case 'C3D10'
                                    SX10nodes = [   1 3 2 7 6 5;
                                        1 2 4 5 9 8;
                                        2 4 3 9 10 6;
                                        1 4 3 8 10 7];
                                    zw = out(i).EDOF(elset(k).Elements,:);
                                    edofzw = zw(:,SX10nodes(str2double(surfaces(j).Side(2)),:));
                                case 'CPE4'
                                    edofzw = out(i).EDOF(elset(k).Elements,:);
                                case 'CPS4R'
                                    SX2nodes = [   1 2;
                                        2 3;
                                        3 4;
                                        4 1];
                                    zw = out(i).EDOF(elset(k).Elements,:);
                                    edofzw = zw(:,SX2nodes(str2double(surfaces(j).Side(2)),:));
                                otherwise
                                    error('Type not implemented')
                            end
                            out(i).SUBSETS = [out(i).SUBSETS struct('NAME',surfaces(j).Name,'type','SURFACE','EDOF',edofzw)];
                        end
                    end
                end
            end
        end
    end
    
%% Translate / rotate geometry (assembly operation in Abaqus)
for i = 1:size(out,2)
    out(i) = translate(out(i),instance(i).Translation);
    out(i) = rotate(out(i),instance(i).Rotation);
end
    
    %% Merge output
    list = [];
    for i = 1:size(obj,2)
        if size(obj(i).NAME,2) > 0
            list = [list i]; %#ok<AGROW>
        end
    end
    out = [obj(list) out];
end
