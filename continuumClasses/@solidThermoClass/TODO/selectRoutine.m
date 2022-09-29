function out = selectRoutine(obj,integrator)
%% Select element routine
%
% Created: 			Mi, 10 Aug 2016
% Responsibilty: 	Mark Schiebl
% Editors: 
% g016464  gk534  gk538  gk604  gk612 
% Description: 
%----------------------------------------------------------------------
%
%
%

%% Check input
if numel(obj) ~= 1
    error('Please use only one object');
end
switch lower(obj.DESIGNATOR)
    case 'disp'
        switch obj.DIM
            case 3
                switch lower(integrator)
                    case 'endpoint'
                        switch lower(obj.MAT.name)
                            case 'linearthermo'
                                out = 'disp_linear_endPoint';
                            case 'neohookethermo'
                                out = 'disp_neoHooke_endPoint';
                            case 'mooneyrivlinthermo'
                                out = 'disp_mooneyRivlin_endPoint';
                            case 'ogdenthermo'
                                out = 'disp_ogden_endPoint';
                            otherwise
                                error('Endpoint 3D material not implemented.')
                        end
                    case 'midpoint'
                        switch lower(obj.MAT.name)
                            case 'neohookethermo'
                                out = 'disp_neoHooke_midPoint';
                            case 'mooneyrivlinthermo'
                                out = 'disp_mooneyRivlin_midPoint';
                            case 'ogdenthermo'
                                out = 'disp_ogden_midPoint';
                            otherwise
                                error('Midpoint 3D material not implemented.')
                        end
                    case 'mixed'
                        switch lower(obj.MAT.name)
                            case 'neohookethermo'
                                out = 'disp_neoHooke_mixed';
                            case 'ogdenthermo'
                                out = 'disp_ogden_mixed';
                            otherwise
                                error('Mixed 3D material not implemented.')
                        end
                    case 'discretegradient'
                        switch lower(obj.MAT.name)
                            case 'neohookethermo'
                                out = 'disp_neoHooke_discreteGradient';
                            case 'ogdenthermo'
                                out = 'disp_ogden_discreteGradient';
                            otherwise
                                error('Discrete gradient 3D material not implemented.')
                        end
                    otherwise
                        error('Time integration scheme is currently not implemented.')
                end
            otherwise
                error('Only 3D implemented.')
        end
    case 'enheas'
        switch obj.DIM
            case 3
                switch lower(integrator)
                    case 'midpoint'
                        switch lower(obj.MAT.name)
                            case 'ogdenthermo'
                                out = 'enheas_ogden_midPoint';
                            otherwise
                                error('Endpoint 3D material not implemented.')
                        end
                    case 'mixed'
                        switch lower(obj.MAT.name)
                            case 'ogdenthermo'
                                out = 'enheas_ogden_mixed';
                            otherwise
                                error('Endpoint 3D material not implemented.')
                        end
                        
                    case 'discretegradient'
                        switch lower(obj.MAT.name)
                            case 'ogdenthermo'
                                out = 'enheas_ogden_discreteGradient';
                            otherwise
                                error('Endpoint 3D material not implemented.')
                        end
                        
                    otherwise
                        error('Time integration scheme is currently not implemented.')
                end
            otherwise
                error('Only 3D implemented.')
        end
    case 'mmhr'
        switch obj.DIM
            case 3
                switch lower(integrator)
                    case 'endpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'mmhr_mooneyRivlin_endPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    otherwise
                        error('Time integration scheme is currently not implemented.')
                end
            otherwise
                error('Only 3D implemented.')
        end
    case 'dispcascadesc'
        switch obj.DIM
            case 3
                switch lower(integrator)
                    case 'endpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'dispCascadeSC_mooneyRivlin_endPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'midpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'dispCascadeSC_mooneyRivlin_midPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'discretegradient'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'dispCascadeSC_mooneyRivlin_discreteGradient';
                            otherwise
                                error('3D material not implemented.')
                        end
                    otherwise
                        error('Time integration scheme is currently not implemented.')
                end
            otherwise
                error('Only 3D implemented.')
        end
    case 'mmhwcascadesc'
        switch obj.DIM
            case 3
                switch lower(integrator)
                    case 'endpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'mmhwCascadeSC_mooneyRivlin_endPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'midpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'mmhwCascadeSC_mooneyRivlin_midPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'discretegradient'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'mmhwCascadeSC_mooneyRivlin_discreteGradient';
                            otherwise
                                error('3D material not implemented.')
                        end
                    otherwise
                        error('Time integration scheme is currently not implemented.')
                end
            otherwise
                error('Only 3D implemented.')
        end
    case 'mmhrsplit'
        switch obj.DIM
            case 3
                switch lower(integrator)
                    case 'endpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'mmhrSplit_mooneyRivlin_endPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'midpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'mmhrSplit_mooneyRivlin_midPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'mixed'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'mmhrSplit_mooneyRivlin_mixed';
                            otherwise
                                error('3D material not implemented.')
                        end
                    otherwise
                        error('Time integration scheme is currently not implemented.')
                end
            otherwise
                error('Only 3D implemented.')
        end
    case 'wedge'
        switch obj.DIM
            case 3
                switch lower(integrator)
                    case 'endpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'wedge_mooneyRivlin_endPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    otherwise
                        error('Time integration scheme is currently not implemented.')
                end
            otherwise
                error('Only 3D implemented.')
        end
    case 'wedgesc'
        switch obj.DIM
            case 3
                switch lower(integrator)
                    case 'endpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'wedgeSC_mooneyRivlin_endPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'midpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'wedgeSC_mooneyRivlin_midPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'discretegradient'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'wedgeSC_mooneyRivlin_discreteGradient';
                            otherwise
                                error('3D material not implemented.')
                        end
                    otherwise
                        error('Time integration scheme is currently not implemented.')
                end
            otherwise
                error('Only 3D implemented.')
        end
    case 'genericcascadesc'
        switch obj.DIM
            case 3
                switch lower(integrator)
                    case 'endpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'genericCascadeSC_mooneyRivlin_endPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'midpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'genericCascadeSC_mooneyRivlin_midPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'discretegradient'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'genericCascadeSC_mooneyRivlin_discreteGradient';
                                %                                 out = 'genericCascadeSC_mooneyRivlin_discreteGradientNumTang';
                                %                                  out = 'genericCascadeSC_mooneyRivlin_discreteGradientNumTang2';
                                %                                 out = 'genericCascadeSC_mooneyRivlin_discreteGradientTC';
                                %                                 out = 'genericCascadeSC_mooneyRivlin_discreteGradientTCNumTang';
                            otherwise
                                error('3D material not implemented.')
                        end
                    otherwise
                        error('Time integration scheme is currently not implemented.')
                end
            otherwise
                error('Only 3D implemented.')
        end
    case 'mmhwgenericcascadesc'
        switch obj.DIM
            case 3
                switch lower(integrator)
                    case 'endpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'mmhwGenericCascadeSC_mooneyRivlin_endPoint';
                                %                                 out = 'mmhwGenericCascadeSC_mooneyRivlin_endPointNumTang';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'midpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'mmhwGenericCascadeSC_mooneyRivlin_midPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'discretegradient'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'mmhwGenericCascadeSC_mooneyRivlin_discreteGradient';
                            otherwise
                                error('3D material not implemented.')
                        end
                    otherwise
                        error('Time integration scheme is currently not implemented.')
                end
            otherwise
                error('Only 3D implemented.')
        end
    case 'wedgesplit'
        switch obj.DIM
            case 3
                switch lower(integrator)
                    case 'endpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'wedgeSplit_mooneyRivlin_endPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'midpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'wedgeSplit_mooneyRivlin_midPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'mixed'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinthermo'
                                out = 'wedgeSplit_mooneyRivlin_mixed';
                            otherwise
                                error('3D material not implemented.')
                        end
                    otherwise
                        error('Time integration scheme is currently not implemented.')
                end
            otherwise
                error('Only 3D implemented.')
        end
    case 'generic'
        switch obj.DIM
            case 3
                switch lower(obj.TAU)
                    case 'entr'
                        switch lower(integrator)
                            case 'midpoint'
                                switch obj.L2Proj
                                    case true
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_entropy_L2_midPoint';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                    case false
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_entropy_midPoint';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                end
                            case 'tc'
                                switch obj.L2Proj
                                    case true
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_entropy_TC_L2';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                    case false
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_entropy_TC';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                end
                            case 'avgstrain'
                                switch obj.L2Proj
                                    case true
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                error('not implemented yet')
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                    case false
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_entropy_avgstrain';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                end
                        end
                    case 'temp'
                        switch lower(integrator)
                            case 'midpoint'
                                switch obj.L2Proj
                                    case true %real TC with L2 Proj
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_temp_midPoint_L2';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                    case false
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_temp_midPoint';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                end
                            case 'tc'
                                switch obj.L2Proj
                                    case true %real TC with L2 Proj
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_temp_TC_L2';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                    case false
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_temp_TC';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                end
                            case 'avgstrain'
                                switch obj.L2Proj
                                    case true
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                error('not implemented yet')
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                    case false
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_temp_avgstrain';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                end
                        end
                    case 'engy'
                        switch lower(integrator)
                            case 'midpoint'
                                switch obj.L2Proj
                                    case true
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_engy_midPoint_L2_v2';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                    case false
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_engy_midPoint';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                end
                            case 'tc'
                                switch obj.L2Proj
                                    case true
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_engy_TC_L2';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                    case false
                                        error('not implemented.')
                                end
                            case 'avgstrain'
                                switch obj.L2Proj
                                    case true
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                error('not implemented yet')
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                    case false
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_engy_avgstrain';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                end
                        end
                    case 'phi'
                        error('not implemented.')
                    otherwise
                        error('there is no such generic format')
                end
            otherwise
                error('Only 3D implemented.')
        end
    case 'generic_wedge'
        switch obj.DIM
            case 3
                switch lower(obj.TAU)
                    case 'temp'
                        switch lower(integrator)
                            case 'tc'
                                switch obj.L2Proj
                                    case true %real TC with L2 Proj
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_wedge_neoHooke_temp_TC_L2';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                    case false
                                        switch lower(obj.MAT.name)
                                            case 'neohookethermo'
                                                out='generic_neoHooke_temp_TC';
                                            otherwise
                                                error('3D material not implemented')
                                        end
                                    otherwise
                                        error('generic needs L2 Projection property')
                                end
                            otherwise
                                error('Time integration scheme is currently not implemented.')
                        end
                    otherwise
                        error('not implemented')
                end
                
            otherwise
                error('Only 3D implemented.')
        end
    otherwise
        error('Identifier not known.')
end
end
