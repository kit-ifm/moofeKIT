function obj = rotate(obj,rotation)
    % Rotation of an geometry object.
    %
    % Syntax
    %
    % rotate(obj,rotation)
    %
    % Description
    %
    % The function rotates the geometry, i.e. obj.NODES using the commited
    % rotation vector. This vector consists of three elements for the position
    % vector of the first point of the rotation axis and three elements for the
    % position vector of the second point of the rotation axis. The last value
    % provides the rotation angle in degrees.
    %
    % 03.01.2012 C.Hesch
    
    %% Check input
    rotation = rotation(:);
    if numel(rotation) ~= 7
%         warning('MATLAB:notEnoughInputs','The rotation vector should have 7 entries, two dimensions are not implemented.')
    end
    
    %% Rotate
    if size(rotation,1) == 7
        a = rotation(1:3);
        b = rotation(4:6);
        e = (b-a)/(norm(b-a));
        eQuer = [0 -e(3) e(2);e(3) 0 -e(1);-e(2) e(1) 0];
        dreh = rotation(7)/360*2*pi;
        R = expm(eQuer*dreh);
        for i = 1:size(obj,2)
            if size(obj(i).QNODE,2) > 0
                obj(i).QNODE = (R*(obj(i).QNODE - kron(a,ones(1,size(obj(i).QNODE,1)))')')' + kron(a,ones(1,size(obj(i).QNODE,1)))';
            end
        end
    end
end