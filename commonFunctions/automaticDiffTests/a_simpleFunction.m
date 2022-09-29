% Generated by ADiMat 0.6.6-5530 (00419e1f)
% © 2001-2008 Andre Vehreschild <vehreschild@sc.rwth-aachen.de>
% © 2009-2018 Johannes Willkomm <johannes@johannes-willkomm.de>
% TU Darmstadt, 64289 Darmstadt, Germany
% Visit us on the web at http://www.adimat.de/
% Report bugs to johannes@johannes-willkomm.de
%
%                             DISCLAIMER
% 
% ADiMat was prepared as part of an employment at the Institute for Scientific Computing,
% RWTH Aachen University, Germany and at the Institute for Scientific Computing,
% TU Darmstadt, Germany and is provided AS IS. 
% NEITHER THE AUTHOR(S), THE GOVERNMENT OF THE FEDERAL REPUBLIC OF GERMANY
% NOR ANY AGENCY THEREOF, NOR THE RWTH AACHEN UNIVERSITY, NOT THE TU DARMSTADT,
% INCLUDING ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
% OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
% OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE
% WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
%
% Parameters:
%  - dependents=y
%  - independents=x, a, b, c
%  - inputEncoding=ISO-8859-1
%
% Functions in this file: a_simpleFunction, rec_simpleFunction,
%  ret_simpleFunction
%

function [a_x a_a a_b a_c nr_y] = a_simpleFunction(x, a, b, c, a_y)
   tmpca2 = x .^ 2;
   tmpca1 = tmpca2 * a;
   y = tmpca1 * b;
   nr_y = y;
   [a_tmpca1 a_tmpca2 a_x a_a a_b a_c] = a_zeros(tmpca1, tmpca2, x, a, b, c);
   if nargin < 5
      a_y = a_zeros1(y);
   end
   a_tmpca1 = adimat_adjsum(a_tmpca1, adimat_adjmultl(tmpca1, a_y, b));
   a_b = adimat_adjsum(a_b, adimat_adjmultr(b, tmpca1, a_y));
   a_tmpca2 = adimat_adjsum(a_tmpca2, adimat_adjmultl(tmpca2, a_tmpca1, a));
   a_a = adimat_adjsum(a_a, adimat_adjmultr(a, tmpca2, a_tmpca1));
   a_x = adimat_adjsum(a_x, adimat_adjred(x, 2 .* x.^1 .* a_tmpca2));
end

function y = rec_simpleFunction(x, a, b, c)
   tmpca2 = x .^ 2;
   tmpca1 = tmpca2 * a;
   y = tmpca1 * b;
   adimat_push(tmpca1, tmpca2, y, x, a, b, c);
end

function [a_x a_a a_b a_c] = ret_simpleFunction(a_y)
   [c b a x y tmpca2 tmpca1] = adimat_pop;
   [a_tmpca1 a_tmpca2 a_x a_a a_b a_c] = a_zeros(tmpca1, tmpca2, x, a, b, c);
   if nargin < 1
      a_y = a_zeros1(y);
   end
   a_tmpca1 = adimat_adjsum(a_tmpca1, adimat_adjmultl(tmpca1, a_y, b));
   a_b = adimat_adjsum(a_b, adimat_adjmultr(b, tmpca1, a_y));
   a_tmpca2 = adimat_adjsum(a_tmpca2, adimat_adjmultl(tmpca2, a_tmpca1, a));
   a_a = adimat_adjsum(a_a, adimat_adjmultr(a, tmpca2, a_tmpca1));
   a_x = adimat_adjsum(a_x, adimat_adjred(x, 2 .* x.^1 .* a_tmpca2));
end