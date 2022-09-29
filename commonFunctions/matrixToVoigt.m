function Av = matrixToVoigt(A, type)
%% Converts matrix into voigt notation
    if strcmpi(type, 'stress')
        if size(A, 1) == 2 && size(A, 2) == 2
            Av = [A(1,1); A(2,2); A(1,2)];
        elseif size(A, 1) == 3 && size(A, 2) == 3
            Av = [A(1,1); A(2,2); A(3,3); A(1,2); A(2,3); A(1,3)];
        end
    elseif strcmpi(type, 'strain')
        if size(A, 1) == 2 && size(A, 2) == 2
            Av = [A(1,1); A(2,2); 2*A(1,2)];
        elseif size(A, 1) == 3 && size(A, 2) == 3
            Av = [A(1,1); A(2,2); A(3,3); 2*A(1,2); 2*A(2,3); 2*A(1,3)];
        end
    end
end