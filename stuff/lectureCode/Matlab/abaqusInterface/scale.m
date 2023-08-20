function obj = scale(obj,inp)
    % Scaling of an geometry object.
    %
    % Syntax
    %
    % scale(obj,scale)
    %
    % Description
    %
    % The function scales the geometry obj.NODES using the commited
    % scalar (uniform) oder vector (direction dependent).
    %
    % 03.01.2012 C.Hesch
    
    %% Check input
    inp = inp(:);
    DIM = size(obj(1).QNODE,2);
    
    if numel(inp) == 1
        inp = inp*ones(DIM,1);
    else
        for i = 1:numel(obj)
            if numel(inp) ~= DIM
                error(['Input must be of size ' num2str(DIM) ' or scalar'])
            end
        end
    end
    
    %% Scale
    for i = 1:numel(obj)
        if size(obj(i).QNODE,2) > 0
            for j = 1:DIM
                obj(i).QNODE(:,j) = obj(i).QNODE(:,j)*inp(j);
            end
        end
    end
end