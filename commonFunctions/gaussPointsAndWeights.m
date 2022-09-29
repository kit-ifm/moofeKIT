function [gaussPoints, gaussWeights] = gaussPointsAndWeights(dimension, numberOfGaussPoints, elementGeometryType)
%GAUSSPOINTSANDWEIGHTS Gauss points and weights in local coordinates.
%   This function returns a matrix with the local coordinates of all
%   Gauss points and a matrix with the corresponding Gauss weights.
%
%   CALL
%   [gaussPoints, gaussWeights] = gaussPointsAndWeights(dimension, numberOfGaussPoints, elementGeometryType)
%   dimension: number of dimensions
%   numberOfGaussPoints: number of Gauss points
%   elementGeometryType: element type
%   gaussPoints: matrix containing the local coordinates of all Gauss points.
%   gaussPoints = [xiGaussPoint1,   xiGaussPoint2,   ..., xiGaussPointNumberOfGaussPoints;
%                  etaGaussPoint1,  etaGaussPoint2,  ..., etaGaussPointNumberOfGaussPoints;
%                  zetaGaussPoint1, zetaGaussPoint2, ..., zetaGaussPointNumberOfGaussPoints];
%   gaussWeights: matrix containing the weights of all Gauss points.
%   gaussWeights = [weightGaussPoint1, weightGaussPoint2, ..., weightGaussPointNumberOfGaussPoints];
%
%   REFERENCE
%   -
%
%   CREATOR(S)
%   Felix Zaehringer, Marlon Franke

gaussPoints = zeros(dimension, numberOfGaussPoints);

if dimension == 1
    if strcmpi(elementGeometryType, 'oneDimensional')
        if mod(numberOfGaussPoints, 1) == 0
            [gaussPoints, gaussWeights] = standardGaussPointsAndWeights(numberOfGaussPoints, dimension);
        else
            error('Number of Gauss points is not implemented!');
        end
    else
        error('Element type is not implemented!');
    end
elseif dimension == 2
    if strcmpi(elementGeometryType, 'triangular')
        if numberOfGaussPoints == 1
            % sufficient for triangular elements with linear shape functions
            gaussPoints = [1 / 3; 1 / 3];
            gaussWeights = 1 / 2;
        elseif numberOfGaussPoints == 3
            % sufficient for triangular elements with quadratic shape functions
            gaussPoints(1, :) = [1 / 2, 0, 1 / 2];
            gaussPoints(2, :) = [1 / 2, 1 / 2, 0];

            % other possible evaluation points
            % gaussPoints(1, :) = [1/6, 2/3, 1/6];
            % gaussPoints(2, :) = [1/6, 1/6, 2/3];

            gaussWeights = [1 / 6, 1 / 6, 1 / 6];
        elseif numberOfGaussPoints == 4
            % sufficient for triangular elements with cubic shape functions
            gaussPoints(1, :) = [1 / 3, 1 / 5, 3 / 5, 1 / 5];
            gaussPoints(2, :) = [1 / 3, 1 / 5, 1 / 5, 3 / 5];
            gaussWeights = [-27 / 96, 25 / 96, 25 / 96, 25 / 96];
        elseif numberOfGaussPoints == 7
            gaussPoints(1, :) = [0.333333333333333, 0.736712498968435, 0.736712498968435, 0.237932366472434, 0.237932366472434, 0.025355134551932, 0.025355134551932];
            gaussPoints(2, :) = [0.333333333333333, 0.237932366472434, 0.025355134551932, 0.736712498968435, 0.025355134551935, 0.736712498968435, 0.237932366472434];
            gaussWeights = [0.375000000000000, 0.104166666666667, 0.104166666666667, 0.104166666666667, 0.104166666666667, 0.104166666666667, 0.104166666666667];
        elseif numberOfGaussPoints == 9
            % Sum of weights equal 1, order 13, degree of precision 7, taken from https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
            gaussPoints(1, :) = [0.124949503233232, 0.437525248383384, 0.437525248383384, 0.797112651860071, 0.797112651860071, 0.165409927389841, 0.165409927389841, 0.037477420750088, 0.037477420750088];
            gaussPoints(2, :) = [0.437525248383384, 0.124949503233232, 0.437525248383384, 0.165409927389841, 0.037477420750088, 0.797112651860071, 0.037477420750088, 0.797112651860071, 0.165409927389841];
            gaussWeights = [0.205950504760887, 0.205950504760887, 0.205950504760887, 0.063691414286223, 0.063691414286223, 0.063691414286223, 0.063691414286223, 0.063691414286223, 0.063691414286223];
        elseif numberOfGaussPoints == 13
            % Sum of weights equal 1, order 28, degree of precision 11, a rule from ACM TOMS algorithm #612
            gaussPoints(1, :) = [0.333333333333333, 0.479308067841923, 0.260345966079038, 0.260345966079038, 0.869739794195568, 0.0651301029022160, 0.0651301029022160, 0.638444188569809, 0.638444188569809, 0.312865496004875, 0.312865496004875, 0.0486903154253160, 0.0486903154253160];
            gaussPoints(2, :) = [0.333333333333333, 0.260345966079038, 0.479308067841923, 0.260345966079038, 0.0651301029022160, 0.869739794195568, 0.0651301029022160, 0.312865496004875, 0.0486903154253160, 0.638444188569809, 0.0486903154253160, 0.638444188569809, 0.312865496004875];
            gaussWeights = [-0.149570044467670, 0.175615257433204, 0.175615257433204, 0.175615257433204, 0.053347235608839, 0.053347235608839, 0.053347235608839, 0.077113760890257, 0.077113760890257, 0.077113760890257, 0.077113760890257, 0.077113760890257, 0.077113760890257];
        elseif numberOfGaussPoints == 28
            % Sum of weights equal 1, order 37, degree of precision 13, a rule from ACM TOMS algorithm #706
            gaussPoints(1, :) = [0.333333333333333, 0.948021718143423, 0.0259891409282883, 0.0259891409282883, 0.811424994704155, 0.0942875026479227, 0.0942875026479227, 0.0107264499655706, 0.494636775017215, 0.494636775017215, 0.585313234770972, 0.207343382614514, 0.207343382614514, 0.122184388599019, 0.438907805700491, 0.438907805700491, 0.677937654882590, 0.677937654882590, 0.0448416775891306, 0.0448416775891306, 0.277220667528279, 0.277220667528279, 0.858870281282636, 0.858870281282636, 0, 0, 0.141129718717364, 0.141129718717364];
            gaussPoints(2, :) = [0.333333333333333, 0.0259891409282883, 0.948021718143423, 0.0259891409282883, 0.0942875026479227, 0.811424994704155, 0.0942875026479227, 0.494636775017215, 0.0107264499655706, 0.494636775017215, 0.207343382614514, 0.585313234770972, 0.207343382614514, 0.438907805700491, 0.122184388599019, 0.438907805700491, 0.0448416775891306, 0.277220667528279, 0.677937654882590, 0.277220667528279, 0.677937654882590, 0.0448416775891306, 0, 0.141129718717364, 0.858870281282636, 0.141129718717364, 0.858870281282636, 0];
            gaussWeights = [0.0879773011622219, 0.00874431155373619, 0.00874431155373619, 0.00874431155373619, 0.0380815719939353, 0.0380815719939353, 0.0380815719939353, 0.0188554480561313, 0.0188554480561313, 0.0188554480561313, 0.0721596975447410, 0.0721596975447410, 0.0721596975447410, 0.0693291387055372, 0.0693291387055372, 0.0693291387055372, 0.0410563154292886, 0.0410563154292886, 0.0410563154292886, 0.0410563154292886, 0.0410563154292886, 0.0410563154292886, 0.00736238378330057, 0.00736238378330057, 0.00736238378330057, 0.00736238378330057, 0.00736238378330057, 0.00736238378330057];
        elseif numberOfGaussPoints == 37
            gaussPoints(1, :) = [0.333333333333333, 0.950275662924106, 0.0248621685379472, 0.0248621685379472, 0.171614914923835, 0.414192542538082, 0.414192542538082, 0.539412243677190, 0.230293878161405, 0.230293878161405, 0.772160036676533, 0.113919981661734, 0.113919981661734, 0.00908539994983535, 0.495457300025082, 0.495457300025082, 0.0622772903058870, 0.468861354847057, 0.468861354847057, 0.0220762896536244, 0.0220762896536244, 0.851306504174349, 0.851306504174349, 0.126617206172027, 0.126617206172027, 0.0186205228025210, 0.0186205228025210, 0.689441970728591, 0.689441970728591, 0.291937506468888, 0.291937506468888, 0.0965064812921592, 0.0965064812921592, 0.635867859433873, 0.635867859433873, 0.267625659273968, 0.267625659273968];
            gaussPoints(2, :) = [0.333333333333333, 0.0248621685379472, 0.950275662924106, 0.0248621685379472, 0.414192542538082, 0.171614914923835, 0.414192542538082, 0.230293878161405, 0.539412243677190, 0.230293878161405, 0.113919981661734, 0.772160036676533, 0.113919981661734, 0.495457300025082, 0.00908539994983535, 0.495457300025082, 0.468861354847057, 0.0622772903058870, 0.468861354847057, 0.851306504174349, 0.126617206172027, 0.0220762896536244, 0.126617206172027, 0.0220762896536244, 0.851306504174349, 0.689441970728591, 0.291937506468888, 0.0186205228025210, 0.291937506468888, 0.0186205228025210, 0.689441970728591, 0.635867859433873, 0.267625659273968, 0.0965064812921592, 0.267625659273968, 0.0965064812921592, 0.635867859433873];
            gaussWeights = [0.00800779955556480, 0.00800779955556480, 0.00800779955556480, 0.0468688989818216, 0.0468688989818216, 0.0468688989818216, 0.0465909401839765, 0.0465909401839765, 0.0465909401839765, 0.0310169433137964, 0.0310169433137964, 0.0310169433137964, 0.0107916127366313, 0.0107916127366313, 0.0107916127366313, 0.0321955342424316, 0.0321955342424316, 0.0321955342424316, 0.0154458342107016, 0.0154458342107016, 0.0154458342107016, 0.0154458342107016, 0.0154458342107016, 0.0154458342107016, 0.0178229899231787, 0.0178229899231787, 0.0178229899231787, 0.0178229899231787, 0.0178229899231787, 0.0178229899231787, 0.0370386836813846, 0.0370386836813846, 0.0370386836813846, 0.0370386836813846, 0.0370386836813846, 0.0370386836813846];
        else
            error('Number of Gauss points is not implemented!');
        end
    elseif strcmpi(elementGeometryType, 'quadrilateral')
        if mod(sqrt(numberOfGaussPoints), 1) == 0
            [gaussPoints, gaussWeights] = standardGaussPointsAndWeights(sqrt(numberOfGaussPoints), dimension);
        elseif numberOfGaussPoints == 5
            % Special 5 point integration according to simo, armero, taylor 1993
            gaussPoints(1, :) = [-0.774596669241483, 0.774596669241483, -0.774596669241483, 0.774596669241483, 0];
            gaussPoints(2, :) = [-0.774596669241483, -0.774596669241483, 0.774596669241483, 0.774596669241483, 0];
            gaussWeights = [0.555555555555556, 0.555555555555556, 0.555555555555556, 0.555555555555556, 16 / 9];
        else
            error('Number of Gauss points is not implemented!');
        end
    else
        error('Element type is not implemented!');
    end
elseif dimension == 3
    if strcmpi(elementGeometryType, 'tetrahedral')
        if numberOfGaussPoints == 1
            % sufficient for tetrahedral elements with linear shape functions
            gaussPoints = [1 / 4; 1 / 4; 1 / 4];
            gaussWeights = 1 / 6;
        elseif numberOfGaussPoints == 4
            a = 0.5854101966249685;
            b = 0.1381966011250105;
            gaussPoints = [a, b, b, b; b, b, b, a; b, b, a, b];
            gaussWeights = [1 / 4, 1 / 4, 1 / 4, 1 / 4] / 6;
        elseif numberOfGaussPoints == 5
            % sufficient for tetrahedral elements with quadratic shape functions
            gaussPoints(1, :) = [1 / 4, 1 / 6, 1 / 6, 1 / 6, 1 / 2];
            gaussPoints(2, :) = [1 / 4, 1 / 6, 1 / 6, 1 / 2, 1 / 6];
            gaussPoints(3, :) = [1 / 4, 1 / 6, 1 / 2, 1 / 6, 1 / 6];
            gaussWeights = [-4 / 5, 9 / 20, 9 / 20, 9 / 20, 9 / 20] / 6;
        elseif numberOfGaussPoints == 11
            gaussPoints(1, :) = [0.2500000000000000, 0.7857142857142857, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.1005964238332008, 0.3994035761667992, 0.3994035761667992, 0.3994035761667992, 0.1005964238332008, 0.1005964238332008];
            gaussPoints(2, :) = [0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.7857142857142857, 0.3994035761667992, 0.1005964238332008, 0.3994035761667992, 0.1005964238332008, 0.3994035761667992, 0.1005964238332008];
            gaussPoints(3, :) = [0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.7857142857142857, 0.0714285714285714, 0.3994035761667992, 0.3994035761667992, 0.1005964238332008, 0.1005964238332008, 0.1005964238332008, 0.3994035761667992];
            gaussWeights = [-0.0789333333333333, 0.0457333333333333, 0.0457333333333333, 0.0457333333333333, 0.0457333333333333, 0.1493333333333333, 0.1493333333333333, 0.1493333333333333, 0.1493333333333333, 0.1493333333333333, 0.1493333333333333] / 6;
        else
            error('Number of Gauss points is not implemented!');
        end
    elseif strcmpi(elementGeometryType, 'hexahedral')
        if mod(nthroot(numberOfGaussPoints, dimension), 1) == 0
            [gaussPoints, gaussWeights] = standardGaussPointsAndWeights(nthroot(numberOfGaussPoints, dimension), dimension);
        elseif numberOfGaussPoints == 4
            gaussPoints(1, :) = [0, 0, sqrt(2/3), -sqrt(2/3)];
            gaussPoints(2, :) = [sqrt(2/3), -sqrt(2/3), 0, 0];
            gaussPoints(3, :) = [-1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3)];
            gaussWeights = [2, 2, 2, 2];
        elseif numberOfGaussPoints == 6
            % sufficient for hexahedral elements with linear shape functions
            % Reference:
            % Irons 1971: Quadrature rules for brick based finite elements
            gaussPoints(1, :) = [1, -1, 0, 0, 0, 0];
            gaussPoints(2, :) = [0, 0, 1, -1, 0, 0];
            gaussPoints(3, :) = [0, 0, 0, 0, 1, -1];
            gaussWeights = [4 / 3, 4 / 3, 4 / 3, 4 / 3, 4 / 3, 4 / 3];
        elseif numberOfGaussPoints == 9
            % Special 9 point integration according to simo, armero, taylor 1993
            g2 = 0.774596669241483;
            w1 = 0.555555555555556;
            gaussPoints(1, :) = [-g2, g2, -g2, g2, -g2, g2, -g2, g2, 0];
            gaussPoints(2, :) = [-g2, -g2, g2, g2, -g2, -g2, g2, g2, 0];
            gaussPoints(3, :) = [-g2, -g2, -g2, -g2, g2, g2, g2, g2, 0];
            gaussWeights = [w1, w1, w1, w1, w1, w1, w1, w1, 32 / 9];
        elseif numberOfGaussPoints == 14
            % sufficient for hexahedral elements with quadratic shape functions
            % Reference:
            % Irons 1971: Quadrature rules for brick based finite elements
            gaussPoints(1, :) = [-sqrt(19/30), sqrt(19/30), 0, 0, 0, 0, -sqrt(19/33), sqrt(19/33), -sqrt(19/33), sqrt(19/33), -sqrt(19/33), sqrt(19/33), -sqrt(19/33), sqrt(19/33)];
            gaussPoints(2, :) = [0, 0, -sqrt(19/30), sqrt(19/30), 0, 0, -sqrt(19/33), -sqrt(19/33), sqrt(19/33), sqrt(19/33), -sqrt(19/33), -sqrt(19/33), sqrt(19/33), sqrt(19/33)];
            gaussPoints(3, :) = [0, 0, 0, 0, -sqrt(19/30), sqrt(19/30), -sqrt(19/33), -sqrt(19/33), -sqrt(19/33), -sqrt(19/33), sqrt(19/33), sqrt(19/33), sqrt(19/33), sqrt(19/33)];
            gaussWeights = [320 / 361, 320 / 361, 320 / 361, 320 / 361, 320 / 361, 320 / 361, 121 / 361, 121 / 361, 121 / 361, 121 / 361, 121 / 361, 121 / 361, 121 / 361, 121 / 361];
        else
            error('Number of Gauss points is not implemented!');
        end
    else
        error('Element type is not implemented!');
    end
else
    error('Number of dimensions is not implemented!');
end
end


function [gaussPoints, gaussWeights] = standardGaussPointsAndWeights(numberOfGaussPoints1D, dimension)
%STANDARDGAUSSPOINTSANDWEIGHTS Standard Gauss points and weights.
%   This function returns the standard Gauss points and Gauss weights for a
%   1D: line
%   2D: quadrilateral
%   3D: hexahedral
%   domain ([-1, 1]^dimension).
%
%   CREATOR(S)
%   Robin Pfefferkorn, Felix Zaehringer

n = 1:numberOfGaussPoints1D - 1;
y = n ./ (sqrt(4.*n.^2-1));
Tn = diag(y, -1) + diag(y, 1);
[eigVec, eigVal] = eig(Tn);

gaussPoints1D = ((1. + diag(eigVal)') - 1);
gaussWeights1D = ((eigVec(1, :).^2) * 2);

gaussPoints = gaussPoints1D;
gaussWeights = gaussWeights1D;

for i = 2:dimension
    gaussPoints = vertcat(repmat(gaussPoints, 1, numberOfGaussPoints1D), kron(gaussPoints1D, ones(1, numberOfGaussPoints1D^(i - 1))));
    gaussWeights = repmat(gaussWeights, 1, numberOfGaussPoints1D) .* kron(gaussWeights1D, ones(1, numberOfGaussPoints1D^(i - 1)));
end
end
