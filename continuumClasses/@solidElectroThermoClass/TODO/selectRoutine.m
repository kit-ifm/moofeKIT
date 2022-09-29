function out = selectRoutine(obj,integrator)
%% Selects element routine
%
% 09.04.2019 M. Franke

%% Check input
if numel(obj) ~= 1
    error('Please use only one object');
end
switch lower(obj.DESIGNATOR)
%     case 'disp'
%         switch obj.DIM
%             case 3
%                 switch lower(integrator)
%                     case 'endpoint'
%                         switch lower(obj.MAT.name)
%                             case 'neohookeelectro'
%                                 out = 'disp_neoHooke_endPoint';
%                             otherwise
%                                 error('Endpoint 3D material not implemented.')
%                         end
%                     otherwise
%                         error('Time integration scheme is currently not implemented.')
%                 end
%             otherwise
%                 error('Only 3D implemented.')
%         end
    case 'mmhwsc'
        switch obj.DIM
            case 3
                switch lower(integrator)
                    case 'endpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinelectrothermo'
                                out = 'mmhwSC_mooneyRivlin_endPoint';
%                             case 'mooneyrivlinelectrothermoti'
%                                 out = 'mmhwSC_mooneyRivlinTI_endPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'midpoint'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinelectrothermo'
                                out = 'mmhwSC_mooneyRivlin_midPoint';
%                             case 'mooneyrivlinelectrothermoti'
%                                 out = 'mmhwSC_mooneyRivlinTI_midPoint';
                            otherwise
                                error('3D material not implemented.')
                        end
                    case 'discretegradient'
                        switch lower(obj.MAT.name)
                            case 'mooneyrivlinelectrothermo'
                                out = 'mmhwSC_mooneyRivlin_discreteGradient';
%                                 out = 'mmhwSC_mooneyRivlin_discreteGradientNumericalTangent';
%                             case 'mooneyrivlinelectrothermoti'
%                                 out = 'mmhwSC_mooneyRivlinTI_discreteGradient';
                            otherwise
                                error('3D material not implemented.')
                        end
                    otherwise
                        error('Time integration scheme is currently not implemented.')
                end
            otherwise
                error('Only 3D implemented.')
        end
%     case 'mmhwcascadesc'
%         switch obj.DIM
%             case 3
%                 switch lower(integrator)
%                     case 'endpoint'
%                         switch lower(obj.MAT.name)
%                             case 'mooneyrivlinelectro'
%                                 out = 'mmhwCascadeSC_mooneyRivlin_endPoint';
%                             case 'mooneyrivlinelectroti'
%                                 out = 'mmhwCascadeSC_mooneyRivlinTI_endPoint';
%                             otherwise
%                                 error('3D material not implemented.')
%                         end
%                     case 'midpoint'
%                         switch lower(obj.MAT.name)
%                             case 'mooneyrivlinelectro'
%                                 out = 'mmhwCascadeSC_mooneyRivlin_midPoint';
%                             case 'mooneyrivlinelectroti'
%                                 out = 'mmhwCascadeSC_mooneyRivlinTI_midPoint';
%                             otherwise
%                                 error('3D material not implemented.')
%                         end
%                     case 'discretegradient'
%                         switch lower(obj.MAT.name)
%                             case 'mooneyrivlinelectro'
%                                 out = 'mmhwCascadeSC_mooneyRivlin_discreteGradient';
%                             case 'mooneyrivlinelectroti'
%                                 out = 'mmhwCascadeSC_mooneyRivlinTI_discreteGradient';
%                             otherwise
%                                 error('3D material not implemented.')
%                         end
%                     otherwise
%                         error('Time integration scheme is currently not implemented.')
%                 end
%             otherwise
%                 error('Only 3D implemented.')
%         end
    otherwise
        error('Identifier not known.')
end
end
