function A = voigtToMatrix(Av, type)
%% Converts voigt notation to matrix
    assert((size(Av,1)==1)||(size(Av,1)==3||size(Av,1)==6),'column size of Av should be either 3 or 6');
    assert((size(Av,2)==1),'row size of Av should be 1');
    if strcmpi(type, 'stress')
        if size(Av, 1) == 3
            A = [   Av(1), Av(3); 
                    Av(3), Av(2)];
        elseif size(Av, 1) == 6
            A = [   Av(1), Av(4), Av(6); 
                    Av(4), Av(2), Av(5); 
                    Av(6), Av(5), Av(3)];
        else
            A=Av;
        end
    elseif strcmpi(type, 'strain')
        if size(Av, 1) == 3
            A = [   Av(1), 0.5*Av(3); 
                    0.5*Av(3), Av(2)];
        elseif size(Av, 1) == 6
            A = [   Av(1), 0.5*Av(4),   0.5*Av(6); 
                0.5*Av(4), Av(2),       0.5*Av(5); 
                0.5*Av(6), 0.5*Av(5),   Av(3)];
        end
    end
end
