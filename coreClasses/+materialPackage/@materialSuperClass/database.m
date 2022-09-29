function [E,nu,rho,G,err] = database(type)
    % Database
    %
    % 03.01.20.12 C.Hesch
    
    err = false;
    switch lower(type)
        case 'academic'
            E   = 1;
            nu  = 0.35;
            rho = 1;
            G   = E/(2 + 2*nu);
            
        case 'steel'
            G   = 79.3e9;
            E   = 210e9;
            nu  = 0.28;
            rho = 7850;
            
        case 'aluminium'
            G   = 25.5e9;
            E   = 70e9;
            nu  = 0.33;
            rho = 2710;
            
        case 'copper'
            G   = 47.0e9;
            E   = 120e9;
            nu  = 0.33;
            rho = 8950;
            
        case 'titanium'
            G   = 41.4e9;
            E   = 105e9;
            nu  = 0.34;
            rho = 4500;
            
        case 'glass'
            G   = 26.2e9;
            E   = 70e9;
            nu  = 0.2;
            rho = 2600;
            
        case 'polyethylene'
            G   = 117e6;
            E   = 200e6;
            nu  = E/(2*G) - 1;
            rho = 925;
            
        case 'rubber'
            G   = 300e3;
            E   = 0.05e9;
            nu  = 0.5;
            rho = 960;
            
        case 'chroscielewski'
            E   = 2e7;
            nu  = 0.25;
            G   = E/(2 + 2*nu);
            rho = 1;
            
        otherwise
            out = [];
            err = true;
    end
end