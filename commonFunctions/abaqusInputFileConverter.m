function out = abaqusInputFileConverter(filename,varargin)
% ABAQUSINPUTFILECONVERTER converts abaqus input files (*.inp) into moofeKIT format.
%
% CALL
% abaqusInputFileConverter(filename,varargin)
% filename: The first argument is the string of the file name of the abaqus
%           mesh, like e.g. 'XXX.inp'
% varargin: The second argument could be be for instance the dimension,
%           e.g. provide 'dimension',3
%
% CREATOR(S)
% Christian Hesch, Rannam Chaaban (esra) 05.04.2010/19.05.2011; Marlon Franke
% FIXME: needs to be rebuild

datei = fopen(filename,'r');
if datei == -1
    error('No file with this name.')
end
dimension = 3;
for i = 1:size(varargin,2)
    if strcmp(varargin{i},'dimension')
        if ~isnumeric(varargin{i+1})
            error('Please provide a numeric as input.')
        end
        dimension = varargin{i+1};
    end
end

%% variables
out = struct('meshName',[],'edof',[],'qR',[],'elementType',[],'subsets',[]);
part = struct('Name',[],'Element',[],'Geometrie',[],'Edof',[]);
instance = struct('Name',[],'Part',[],'Translation',[],'Rotation',[]);
shell = struct('Thickness',[],'IntegrationPoints',[],'ElementSet',[],'Material',[]);
section = struct('ElementSet',[],'Material',[]);
nset = struct('Name',[],'Instance',[],'Nodes',[],'Part',[]);
elset = struct('Name',[],'Instance',[],'Elements',[],'Part',[]);
surfaces = struct('Name',[],'Elset',[],'Type',[],'Side',[]);
runner = 0;
nsetRunner = 0;
elsetRunner = 0;
partRunner = 0;
shellRunner = 0;
sectionRunner = 0;
instanceRunner = 0;
surfRunner = 0;
partSwitch = false;
%% process input file
tline = fgets(datei);
while ~feof(datei)
    if size(tline,2) > 0 && strcmp(tline(1),'*')
        aline = textscan(tline, '%s', 'delimiter', ',');
    else
        aline{1} = {'ABC'};
    end
    switch aline{1}{1}
        case '*Part'
            partRunner = partRunner + 1;
            aline = textscan(tline, '%s %s', 'delimiter', '=');
            part(partRunner).Name = aline{2}{1};
            partSwitch = true;
            runner = runner + 1;
            tline = fgets(datei);
        case '*End Part'
            partSwitch = false;
            runner = runner + 1;
            tline = fgets(datei);
        case '*Node'
            runner = runner + 1;
            tline = fgets(datei);
            startNode = runner;
            while ~strncmp(tline,'*',1) && ~feof(datei)
                runner = runner + 1;
                tline = fgets(datei);
            end
            a = textscan(tline, '%s','delimiter',',');
            zw = dlmread(filename,',',[startNode 0 runner-1 dimension]);
            part(partRunner).Geometrie(zw(:,1),:) = zw(:,2:end);
            clear zw
        case '*Element'
            aline = textscan(tline, '*%s %s','delimiter','=');
            part(partRunner).Element = aline{2}{1};
            runner = runner + 1;
            tline = fgets(datei);
            startNode = runner;
            a = textscan(tline, '%s','delimiter',',');
            rightCorner = size(a{1},1)-1;
            while ~strncmp(tline,'*',1) && ~feof(datei)
                runner = runner + 1;
                tline = fgets(datei);
            end
            zw = dlmread(filename,',',[startNode 0 runner-1 rightCorner]);
            if strcmp(aline{2},'C3D20R')
                part(partRunner).Edof = [zw(1:2:end,2:end) zw(2:2:end,1:5)];
            else
                part(partRunner).Edof = zw(:,2:end);
            end
        case '*Surface'
            if strcmp(aline{1}{2},'type=ELEMENT')
                surfRunner = surfRunner + 1;
                aline = textscan(tline, '%s','delimiter',',');
                bline = textscan(aline{1}{3}, '%s','delimiter','=');
                surfaces(surfRunner).Name = bline{1}{2};
                bline = textscan(aline{1}{2}, '%s','delimiter','=');
                if strcmp(bline{1}(1),'type')
                    surfaces(surfRunner).Type = bline{1}{2};
                end
                runner = runner + 1;
                tline = fgets(datei);
                aline = textscan(tline, '%s','delimiter',',');
                surfaces(surfRunner).Elset = aline{1}{1};
                surfaces(surfRunner).Side = aline{1}{2};
                runner = runner + 1;
                tline = fgets(datei);
                while ~strncmp(tline,'*',1) && ~feof(datei)
                    surfRunner = surfRunner + 1;
                    surfaces(surfRunner).Name = surfaces(surfRunner-1).Name;
                    surfaces(surfRunner).Type = surfaces(surfRunner-1).Type;
                    aline = textscan(tline, '%s','delimiter',',');
                    surfaces(surfRunner).Elset = aline{1}{1};
                    surfaces(surfRunner).Side = aline{1}{2};
                    runner = runner + 1;
                    tline = fgets(datei);
                end
            else
                runner = runner + 1;
                tline = fgets(datei);
            end
        case '*Shell Section'
            shellRunner = shellRunner + 1;
            aline = textscan(tline, '%s','delimiter',',');
            bline = textscan(aline{1}{2}, '%s','delimiter','=');
            shell(shellRunner).ElementSet = bline{1}{2};
            bline = textscan(aline{1}{3}, '%s','delimiter','=');
            shell(shellRunner).Material = bline{1}{2};
            runner = runner + 1;
            tline = fgets(datei);
            cline = textscan(tline, '%f','delimiter',',');
            shell(shellRunner).Thickness = cline{1}(1);
            shell(shellRunner).IntegrationPoints = cline{1}(2);
            runner = runner + 1;
            tline = fgets(datei);
        case '*Solid Section'
            sectionRunner = sectionRunner + 1;
            aline = textscan(tline, '%s','delimiter',',');
            bline = textscan(aline{1}{2}, '%s','delimiter','=');
            section(sectionRunner).ElementSet = bline{1}{2};
            bline = textscan(aline{1}{3}, '%s','delimiter','=');
            section(sectionRunner).Material = bline{1}{2};
            runner = runner + 1;
            tline = fgets(datei);
        case '*Instance'
            instanceRunner = instanceRunner + 1;
            aline = textscan(tline, '%s','delimiter',',');
            bline = textscan(aline{1}{2}, '%s','delimiter','=');
            instance(instanceRunner).Name = bline{1}{2};
            bline = textscan(aline{1}{3}, '%s','delimiter','=');
            instance(instanceRunner).Part = bline{1}{2};
            runner = runner + 1;
            tline = fgets(datei);
            if strncmp(tline,'*',1)
                if dimension == 2
                    instance(instanceRunner).Translation = [0 0];
                    instance(instanceRunner).Rotation = [1 0 0 0 0];
                elseif dimension == 3
                    instance(instanceRunner).Translation = [0 0 0];
                    instance(instanceRunner).Rotation = [1 0 0 0 0 0 0];
                end
            else
                if dimension == 2
                    instance(instanceRunner).Translation = [0 0];
                    instance(instanceRunner).Rotation = [1 0 0 0 0];
                    if ~strncmp(tline,'*',1)
                        instance(instanceRunner).Translation = sscanf(tline, '%f, %f');
                    end
                    runner = runner + 1;
                    tline = fgets(datei);
                    if ~strncmp(tline,'*',1)
                        instance(instanceRunner).Rotation = sscanf(tline, '%f, %f, %f, %f, %f');
                    end
                elseif dimension == 3
                    instance(instanceRunner).Translation = [0 0 0];
                    instance(instanceRunner).Rotation = [1 0 0 0 0 0 0];
                    if ~strncmp(tline,'*',1)
                        instance(instanceRunner).Translation = sscanf(tline, '%f, %f, %f');
                    end
                    runner = runner + 1;
                    tline = fgets(datei);
                    if ~strncmp(tline,'*',1)
                        instance(instanceRunner).Rotation = sscanf(tline, '%f, %f, %f, %f, %f, %f, %f');
                    end
                end
            end
        case '*Nset'
            nsetRunner = nsetRunner + 1;
            if partSwitch
                nset(nsetRunner).Part = part(partRunner).Name;
            end
            aline = textscan(tline, '%s','delimiter',',');
            bline = textscan(aline{1}{2}, '%s','delimiter','=');
            nset(nsetRunner).Name = bline{1}{2};
            for i = 2:size(aline{1},1)
                bline = textscan(aline{1}{i}, '%s','delimiter','=');
                if strncmp(bline{1}(1),'instance',8)
                    nset(nsetRunner).Instance = bline{1}{2};
                end
            end
            if strncmp(aline{1}(end),'generate',8)
                runner = runner + 1;
                tline = fgets(datei);
                cline = textscan(tline, '%f','delimiter',',');
                nset(nsetRunner).Nodes = (cline{1}(1):cline{1}(3):cline{1}(2))';
            else
                runner = runner + 1;
                tline = fgets(datei);
                startNode = runner;
                a = textscan(tline, '%s','delimiter',',');
                rightCorner = size(a{1},1)-1;
                while ~strncmp(tline,'*',1) && ~feof(datei)
                    runner = runner + 1;
                    tline = fgets(datei);
                end
                nodesMatrix = dlmread(filename,',',[startNode 0 runner-1 rightCorner]);
                zw = nodesMatrix(:);
                logikVektor = zw > 0;
                nset(nsetRunner).Nodes = zw(logikVektor);
            end
        case '*Elset'
            elsetRunner = elsetRunner + 1;
            if partSwitch
                elset(elsetRunner).Part = part(partRunner).Name;
            end
            aline = textscan(tline, '%s','delimiter',',');
            bline = textscan(aline{1}{2}, '%s','delimiter','=');
            elset(elsetRunner).Name = bline{1}{2};
            for i = 2:size(aline{1},1)
                bline = textscan(aline{1}{i}, '%s','delimiter','=');
                if strncmp(bline{1}(1),'instance',8)
                    elset(elsetRunner).Instance = bline{1}{2};
                end
            end
            if strncmp(aline{1}(end),'generate',8)
                runner = runner + 1;
                tline = fgets(datei);
                cline = textscan(tline, '%f','delimiter',',');
                elset(elsetRunner).Elements = (cline{1}(1):cline{1}(3):cline{1}(2))';
            else
                runner = runner + 1;
                tline = fgets(datei);
                startNode = runner;
                a = textscan(tline, '%s','delimiter',',');
                rightCorner = size(a{1},1)-1;
                while ~strncmp(tline,'*',1) && ~feof(datei)
                    runner = runner + 1;
                    tline = fgets(datei);
                end
                nodesMatrix = dlmread(filename,',',[startNode 0 runner-1 rightCorner]);
                zw = nodesMatrix(:);
                logikVektor = zw > 0;
                elset(elsetRunner).Elements = zw(logikVektor);
            end
        otherwise
            runner = runner + 1;
            tline = fgets(datei);
    end
end
fclose(datei);
%% merge data
for i = 1:size(instance,2)
    for j = 1:size(part,2)
        if strcmp(instance(i).Part,part(j).Name)
            out(i).meshName = part(j).Name;
            out(i).elementType = part(j).Element;
            out(i).qR = part(j).Geometrie;
            out(i).edof = part(j).Edof;
        end
    end
    % node sets
    for j = 1:size(nset,2)
        if size(elset(j).Instance,2) > 0
            if strcmp(nset(j).Instance,instance(i).Name)
                out(i).subsets = [out(i).subsets struct('name',nset(j).Name,'type','nodeSet','edof',nset(j).Nodes)];
            end
        else
            if strcmp(elset(j).Part,instance(i).Part)
                out(i).subsets = [out(i).subsets struct('name',nset(j).Name,'type','nodeSet','edof',nset(j).Nodes)];
            end
        end
    end
    % element sets
    for j = 1:size(elset,2)
        if size(elset(j).Instance,2) > 0
            if strcmp(elset(j).Instance,instance(i).Name)
                out(i).subsets = [out(i).subsets struct('name',elset(j).Name,'type','elementSet','edof',out(i).edof(elset(j).Elements,:))];
            end
        else
            if strcmp(elset(j).Part,instance(i).Part)
                out(i).subsets = [out(i).subsets struct('name',elset(j).Name,'type','elementSet','edof',out(i).edof(elset(j).Elements,:))];
            end
        end
    end
    % surfaces
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
                    for m = 1:size(out(i).subsets,2)
                        if strcmp(out(i).subsets(m).name,surfaces(j).Name)
                            foundSurf = true;
                            switch out(i).elementType
                                case 'C3D8R'
                                    zw = out(i).edof(elset(k).Elements,:);
                                    out(i).subsets(m).edof = [out(i).subsets(m).edof; zw(:,SX(str2double(surfaces(j).Side(2)),:))];
                                case 'C3D4'
                                    SX4nodes = [1 3 2 ;
                                        1 2 4 ;
                                        2 3 4 ;
                                        1 4 3 ];
                                    zw = out(i).edof(elset(k).Elements,:);
                                    out(i).subsets(m).edof = [out(i).subsets(m).edof; zw(:,SX4nodes(str2double(surfaces(j).Side(2)),:))];
                                case 'C3D20R'
                                    SX20nodes = [   1 2 3 4 9 10 11 12;
                                        5 8 7 6 16 15 14 13;
                                        1 5 6 2 17 13 18 9;
                                        2 6 7 3 18 14 19 10;
                                        3 7 8 4 19 15 20 11;
                                        4 8 5 1 20 16 17 12];
                                    zw = out(i).edof(elset(k).Elements,:);
                                    out(i).subsets(m).edof = [out(i).subsets(m).edof; zw(:,SX20nodes(str2double(surfaces(j).Side(2)),:))];
                                case 'C3D10'
                                    SX10nodes = [   1 3 2 7 6 5;
                                        1 2 4 5 9 8;
                                        2 4 3 9 10 6;
                                        1 4 3 8 10 7];
                                    zw = out(i).edof(elset(k).Elements,:);
                                    out(i).subsets(m).edof = [out(i).subsets(m).edof; zw(:,SX10nodes(str2double(surfaces(j).Side(2)),:))];
                                case 'CPE4'
                                    out(i).subsets(m).edof = [out(i).subsets(m).edof; out(i).edof(elset(k).Elements,:)];
                                case 'CPS4R'
                                    SX2nodes = [   1 2;
                                        2 3;
                                        3 4;
                                        4 1];
                                    zw = out(i).edof(elset(k).Elements,:);
                                    out(i).subsets(m).edof = [out(i).subsets(m).edof; zw(:,SX2nodes(str2double(surfaces(j).Side(2)),:))];
                                otherwise
                                    error('Type not implemented')
                            end
                        end
                    end
                    if ~foundSurf
                        switch out(i).elementType
                            case 'C3D8R'
                                zw = out(i).edof(elset(k).Elements,:);
                                edofzw = zw(:,SX(str2double(surfaces(j).Side(2)),:));
                            case 'C3D4'
                                SX4nodes = [    1 3 2 ;
                                    1 2 4 ;
                                    2 4 3 ;
                                    1 4 3 ];
                                zw = out(i).edof(elset(k).Elements,:);
                                edofzw = zw(:,SX4nodes(str2double(surfaces(j).Side(2)),:));
                            case 'C3D20R'
                                SX20nodes = [   1 2 3 4 9 10 11 12;
                                    5 8 7 6 16 15 14 13;
                                    1 5 6 2 17 13 18 9;
                                    2 6 7 3 18 14 19 10;
                                    3 7 8 4 19 15 20 11;
                                    4 8 5 1 20 16 17 12];
                                zw = out(i).edof(elset(k).Elements,:);
                                edofzw = zw(:,SX20nodes(str2double(surfaces(j).Side(2)),:));
                            case 'C3D10'
                                SX10nodes = [   1 3 2 7 6 5;
                                    1 2 4 5 9 8;
                                    2 4 3 9 10 6;
                                    1 4 3 8 10 7];
                                zw = out(i).edof(elset(k).Elements,:);
                                edofzw = zw(:,SX10nodes(str2double(surfaces(j).Side(2)),:));
                            case 'CPE4'
                                edofzw = out(i).edof(elset(k).Elements,:);
                            case 'CPS4R'
                                SX2nodes = [   1 2;
                                    2 3;
                                    3 4;
                                    4 1];
                                zw = out(i).edof(elset(k).Elements,:);
                                edofzw = zw(:,SX2nodes(str2double(surfaces(j).Side(2)),:));
                            otherwise
                                error('Type not implemented')
                        end
                        out(i).subsets = [out(i).subsets struct('name',surfaces(j).Name,'type','surface','edof',edofzw)];
                    end
                end
            end
        end
    end
end
end