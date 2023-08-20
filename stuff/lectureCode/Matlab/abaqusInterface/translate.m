function obj = translate(obj,translation)
    % Translation of an geometry object.
    %
    % Syntax
    %
    % translate(obj,translation)
    %
    % Description
    %
    % The function tranlates the geometry obj.NODES using the commited
    % translation vector.
    %
    % 03.01.2012 C.Hesch
    
    %% Check input
    translation = translation(:);
    if size(translation,1) ~= 2 && size(translation,1) ~= 3
        error('The translation vector MUST have 2 or 3 entries')
    end
    
    %% Translate
    for i = 1:size(obj,2)
        if size(obj(i).QNODE,2) > 0
            obj(i).QNODE = obj(i).QNODE + kron(translation,ones(1,size(obj(i).QNODE,1)))';
        end
    end
end