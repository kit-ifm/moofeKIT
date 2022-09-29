% Generated by ADiMat 0.6.0-4867
% © 2001-2008 Andre Vehreschild <vehreschild@sc.rwth-aachen.de>
% © 2009-2013 Johannes Willkomm <johannes.willkomm@sc.tu-darmstadt.de>
% RWTH Aachen University, 52056 Aachen, Germany
% TU Darmstadt, 64289 Darmstadt, Germany
% Visit us on the web at http://www.adimat.de/
% Report bugs to adimat-users@lists.sc.informatik.tu-darmstadt.de
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
% Flags: BACKWARDMODE,  NOOPEROPTIM,
%   NOLOCALCSE,  NOGLOBALCSE,  NOPRESCALARFOLDING,
%   NOPOSTSCALARFOLDING,  NOCONSTFOLDMULT0,  FUNCMODE,
%   NOTMPCLEAR,  DUMP_XML,  PARSE_ONLY,
%   UNBOUND_ERROR
%
% Parameters:
%  - dependents=r
%  - independents=p
%  - inputEncoding=ISO-8859-1
%  - output-mode: plain
%  - output-file: ad_out/a_adimat_norm2_p2.m
%  - output-file-prefix: 
%  - output-directory: ad_out
% Generated by ADiMat 0.6.0-4867
% © 2001-2008 Andre Vehreschild <vehreschild@sc.rwth-aachen.de>
% © 2009-2013 Johannes Willkomm <johannes.willkomm@sc.tu-darmstadt.de>
% RWTH Aachen University, 52056 Aachen, Germany
% TU Darmstadt, 64289 Darmstadt, Germany
% Visit us on the web at http://www.adimat.de/
% Report bugs to adimat-users@lists.sc.informatik.tu-darmstadt.de
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
% Flags: BACKWARDMODE,  NOOPEROPTIM,
%   NOLOCALCSE,  NOGLOBALCSE,  NOPRESCALARFOLDING,
%   NOPOSTSCALARFOLDING,  NOCONSTFOLDMULT0,  FUNCMODE,
%   NOTMPCLEAR,  DUMP_XML,  PARSE_ONLY,
%   UNBOUND_ERROR
%
% Parameters:
%  - dependents=r
%  - independents=p
%  - inputEncoding=ISO-8859-1
%  - output-mode: plain
%  - output-file: ad_out/a_adimat_norm2_p2.m
%  - output-file-prefix: 
%  - output-directory: ad_out
%
% Functions in this file: a_adimat_norm2, rec_adimat_norm2,
%  ret_adimat_norm2
%

function [a_p nr_r] = a_adimat_norm2(x, p, a_r)
% r = norm(x, p);
   tmpda3 = 0;
   tmpda2 = 0;
   tmpda1 = 0;
   tmpca3 = 0;
   tmpca2 = 0;
   tmpca1 = 0;
   r = 0;
   answer = 0;
   a = 0;
   sa2 = 0;
   s = 0;
   tmpba1 = 0;
   if ischar(p)
      tmpba1 = 1;
      tmpba2 = 0;
      if strcmp(lower(p), 'fro')
         tmpba2 = 1;
         adimat_push1(tmpda3);
         tmpda3 = conj(x(:));
         adimat_push1(tmpda2);
         tmpda2 = x(:) .* tmpda3;
         adimat_push1(tmpda1);
         tmpda1 = sum(tmpda2);
         adimat_push1(r);
         r = sqrt(tmpda1);
      else
         error('Only "fro" is a valid string for p-norm computation currently.');
      end
      adimat_push1(tmpba2);
   else
      tmpba2 = 0;
      if isvector(x)
         tmpba2 = 1;
         tmpba3 = 0;
         if isinf(p)
            tmpba3 = 1;
            tmpba4 = 0;
            if p > 0
               tmpba4 = 1;
               adimat_push1(tmpda1);
               tmpda1 = abs(x);
               adimat_push1(r);
               r = max(tmpda1);
            else
               adimat_push1(tmpda1);
               tmpda1 = abs(x);
               adimat_push1(r);
               r = min(tmpda1);
            end
            adimat_push1(tmpba4);
         else
            tmpba4 = 0;
            if isreal(x) && mod(p, 2)==0
               tmpba4 = 1;
               adimat_push1(answer);
               answer = admGetPref('pnormEven_p_useAbs');
               tmpba5 = 0;
               if strcmp(answer, 'yes')
                  tmpba5 = 1;
                  adimat_push1(a);
                  a = abs(x);
               else
                  adimat_push1(a);
                  a = x;
               end
               adimat_push1(tmpba5);
            else
               adimat_push1(a);
               a = abs(x);
            end
            adimat_push(tmpba4, tmpca3);
            tmpca3 = 1 / p;
            adimat_push1(tmpca2);
            tmpca2 = a .^ p;
            adimat_push1(tmpca1);
            tmpca1 = sum(tmpca2);
            adimat_push1(r);
            r = tmpca1 .^ tmpca3;
         end
         adimat_push1(tmpba3);
      elseif ismatrix(x)
         tmpba2 = 2;
         tmpba3 = 0;
         if isinf(p)
            tmpba3 = 1;
            adimat_push1(a);
            a = abs(x);
            adimat_push1(sa2);
            sa2 = sum(a, 2);
            adimat_push1(r);
            r = max(sa2);
         elseif p == 2
            tmpba3 = 2;
            tmpba4 = 0;
            if issparse(x)
               tmpba4 = 1;
               adimat_push1(x);
               x = full(x);
            end
            adimat_push1(tmpba4);
            tmpba4 = 0;
            if isreal(x)
               tmpba4 = 1;
               adimat_push1(s);
               s = svd(x);
            else
               adimat_push1(s);
               [s] = adimat_svd(x);
            end
            adimat_push(tmpba4, r);
            r = max(s);
         elseif p == 1
            tmpba3 = 3;
            adimat_push1(a);
            a = abs(x);
            adimat_push1(sa2);
            sa2 = sum(a, 1);
            adimat_push1(r);
            r = max(sa2);
         else
            error('Derivatives of matrix-p-norm not implemented yet.');
         end
         adimat_push1(tmpba3);
      else
         error('Value is neither a matrix nor a vector!');
      end
      adimat_push1(tmpba2);
   end
   adimat_push1(tmpba1);
   nr_r = r;
   [a_tmpca3 a_tmpca2 a_tmpca1 a_p] = a_zeros(tmpca3, tmpca2, tmpca1, p);
   if nargin < 3
      a_r = a_zeros1(r);
   end
   tmpba1 = adimat_pop1;
   if tmpba1 == 1
      tmpba2 = adimat_pop1;
      if tmpba2 == 1
         r = adimat_pop1;
         a_r = a_zeros1(r);
         [tmpda1 tmpda2 tmpda3] = adimat_pop;
      else
         error('Only "fro" is a valid string for p-norm computation currently.');
      end
   else
      tmpba2 = adimat_pop1;
      if tmpba2 == 1
         tmpba3 = adimat_pop1;
         if tmpba3 == 1
            tmpba4 = adimat_pop1;
            if tmpba4 == 1
               r = adimat_pop1;
               a_r = a_zeros1(r);
               tmpda1 = adimat_pop1;
            else
               r = adimat_pop1;
               a_r = a_zeros1(r);
               tmpda1 = adimat_pop1;
            end
         else
            r = adimat_pop1;
            a_tmpca1 = adimat_adjsum(a_tmpca1, adimat_adjred(tmpca1, tmpca3 .* tmpca1.^(tmpca3 - 1) .* a_r));
            a_tmpca3 = adimat_adjsum(a_tmpca3, adimat_adjred(tmpca3, tmpca1.^tmpca3 .* log(tmpca1) .* a_r));
            a_r = a_zeros1(r);
            tmpca1 = adimat_pop1;
            a_tmpca2 = adimat_adjsum(a_tmpca2, a_sum(a_tmpca1, tmpca2));
            a_tmpca1 = a_zeros1(tmpca1);
            tmpca2 = adimat_pop1;
            a_p = adimat_adjsum(a_p, adimat_adjred(p, a.^p .* log(a) .* a_tmpca2));
            a_tmpca2 = a_zeros1(tmpca2);
            [tmpadjc2] = adimat_a_mrdivider(1, p, a_tmpca3);
            tmpca3 = adimat_pop1;
            a_p = adimat_adjsum(a_p, tmpadjc2);
            a_tmpca3 = a_zeros1(tmpca3);
            tmpba4 = adimat_pop1;
            if tmpba4 == 1
               tmpba5 = adimat_pop1;
               if tmpba5 == 1
                  a = adimat_pop1;
               else
                  a = adimat_pop1;
               end
               answer = adimat_pop1;
            else
               a = adimat_pop1;
            end
         end
      elseif tmpba2 == 2
         tmpba3 = adimat_pop1;
         if tmpba3 == 1
            r = adimat_pop1;
            a_r = a_zeros1(r);
            [sa2 a] = adimat_pop;
         elseif tmpba3 == 2
            r = adimat_pop1;
            a_r = a_zeros1(r);
            tmpba4 = adimat_pop1;
            if tmpba4 == 1
               s = adimat_pop1;
            else
               s = adimat_pop1;
            end
            tmpba4 = adimat_pop1;
            if tmpba4 == 1
               x = adimat_pop1;
            end
         elseif tmpba3 == 3
            r = adimat_pop1;
            a_r = a_zeros1(r);
            [sa2 a] = adimat_pop;
         else
            error('Derivatives of matrix-p-norm not implemented yet.');
         end
      else
         error('Value is neither a matrix nor a vector!');
      end
   end
end

function r = rec_adimat_norm2(x, p)
   tmpda3 = 0;
   tmpda2 = 0;
   tmpda1 = 0;
   tmpca3 = 0;
   tmpca2 = 0;
   tmpca1 = 0;
   r = 0;
   answer = 0;
   a = 0;
   sa2 = 0;
   s = 0;
   tmpba1 = 0;
   if ischar(p)
      tmpba1 = 1;
      tmpba2 = 0;
      if strcmp(lower(p), 'fro')
         tmpba2 = 1;
         adimat_push1(tmpda3);
         tmpda3 = conj(x(:));
         adimat_push1(tmpda2);
         tmpda2 = x(:) .* tmpda3;
         adimat_push1(tmpda1);
         tmpda1 = sum(tmpda2);
         adimat_push1(r);
         r = sqrt(tmpda1);
      else
         error('Only "fro" is a valid string for p-norm computation currently.');
      end
      adimat_push1(tmpba2);
   else
      tmpba2 = 0;
      if isvector(x)
         tmpba2 = 1;
         tmpba3 = 0;
         if isinf(p)
            tmpba3 = 1;
            tmpba4 = 0;
            if p > 0
               tmpba4 = 1;
               adimat_push1(tmpda1);
               tmpda1 = abs(x);
               adimat_push1(r);
               r = max(tmpda1);
            else
               adimat_push1(tmpda1);
               tmpda1 = abs(x);
               adimat_push1(r);
               r = min(tmpda1);
            end
            adimat_push1(tmpba4);
         else
            tmpba4 = 0;
            if isreal(x) && mod(p, 2)==0
               tmpba4 = 1;
               adimat_push1(answer);
               answer = admGetPref('pnormEven_p_useAbs');
               tmpba5 = 0;
               if strcmp(answer, 'yes')
                  tmpba5 = 1;
                  adimat_push1(a);
                  a = abs(x);
               else
                  adimat_push1(a);
                  a = x;
               end
               adimat_push1(tmpba5);
            else
               adimat_push1(a);
               a = abs(x);
            end
            adimat_push(tmpba4, tmpca3);
            tmpca3 = 1 / p;
            adimat_push1(tmpca2);
            tmpca2 = a .^ p;
            adimat_push1(tmpca1);
            tmpca1 = sum(tmpca2);
            adimat_push1(r);
            r = tmpca1 .^ tmpca3;
         end
         adimat_push1(tmpba3);
      elseif ismatrix(x)
         tmpba2 = 2;
         tmpba3 = 0;
         if isinf(p)
            tmpba3 = 1;
            adimat_push1(a);
            a = abs(x);
            adimat_push1(sa2);
            sa2 = sum(a, 2);
            adimat_push1(r);
            r = max(sa2);
         elseif p == 2
            tmpba3 = 2;
            tmpba4 = 0;
            if issparse(x)
               tmpba4 = 1;
               adimat_push1(x);
               x = full(x);
            end
            adimat_push1(tmpba4);
            tmpba4 = 0;
            if isreal(x)
               tmpba4 = 1;
               adimat_push1(s);
               s = svd(x);
            else
               adimat_push1(s);
               [s] = adimat_svd(x);
            end
            adimat_push(tmpba4, r);
            r = max(s);
         elseif p == 1
            tmpba3 = 3;
            adimat_push1(a);
            a = abs(x);
            adimat_push1(sa2);
            sa2 = sum(a, 1);
            adimat_push1(r);
            r = max(sa2);
         else
            error('Derivatives of matrix-p-norm not implemented yet.');
         end
         adimat_push1(tmpba3);
      else
         error('Value is neither a matrix nor a vector!');
      end
      adimat_push1(tmpba2);
   end
   adimat_push(tmpba1, answer, a, sa2, s, tmpda3, tmpda2, tmpda1, tmpca3, tmpca2, tmpca1, r, x, p);
end

function a_p = ret_adimat_norm2(a_r)
   [p x r tmpca1 tmpca2 tmpca3 tmpda1 tmpda2 tmpda3 s sa2 a answer] = adimat_pop;
   [a_tmpca3 a_tmpca2 a_tmpca1 a_p] = a_zeros(tmpca3, tmpca2, tmpca1, p);
   if nargin < 1
      a_r = a_zeros1(r);
   end
   tmpba1 = adimat_pop1;
   if tmpba1 == 1
      tmpba2 = adimat_pop1;
      if tmpba2 == 1
         r = adimat_pop1;
         a_r = a_zeros1(r);
         [tmpda1 tmpda2 tmpda3] = adimat_pop;
      else
         error('Only "fro" is a valid string for p-norm computation currently.');
      end
   else
      tmpba2 = adimat_pop1;
      if tmpba2 == 1
         tmpba3 = adimat_pop1;
         if tmpba3 == 1
            tmpba4 = adimat_pop1;
            if tmpba4 == 1
               r = adimat_pop1;
               a_r = a_zeros1(r);
               tmpda1 = adimat_pop1;
            else
               r = adimat_pop1;
               a_r = a_zeros1(r);
               tmpda1 = adimat_pop1;
            end
         else
            r = adimat_pop1;
            a_tmpca1 = adimat_adjsum(a_tmpca1, adimat_adjred(tmpca1, tmpca3 .* tmpca1.^(tmpca3 - 1) .* a_r));
            a_tmpca3 = adimat_adjsum(a_tmpca3, adimat_adjred(tmpca3, tmpca1.^tmpca3 .* log(tmpca1) .* a_r));
            a_r = a_zeros1(r);
            tmpca1 = adimat_pop1;
            a_tmpca2 = adimat_adjsum(a_tmpca2, a_sum(a_tmpca1, tmpca2));
            a_tmpca1 = a_zeros1(tmpca1);
            tmpca2 = adimat_pop1;
            a_p = adimat_adjsum(a_p, adimat_adjred(p, a.^p .* log(a) .* a_tmpca2));
            a_tmpca2 = a_zeros1(tmpca2);
            [tmpadjc2] = adimat_a_mrdivider(1, p, a_tmpca3);
            tmpca3 = adimat_pop1;
            a_p = adimat_adjsum(a_p, tmpadjc2);
            a_tmpca3 = a_zeros1(tmpca3);
            tmpba4 = adimat_pop1;
            if tmpba4 == 1
               tmpba5 = adimat_pop1;
               if tmpba5 == 1
                  a = adimat_pop1;
               else
                  a = adimat_pop1;
               end
               answer = adimat_pop1;
            else
               a = adimat_pop1;
            end
         end
      elseif tmpba2 == 2
         tmpba3 = adimat_pop1;
         if tmpba3 == 1
            r = adimat_pop1;
            a_r = a_zeros1(r);
            [sa2 a] = adimat_pop;
         elseif tmpba3 == 2
            r = adimat_pop1;
            a_r = a_zeros1(r);
            tmpba4 = adimat_pop1;
            if tmpba4 == 1
               s = adimat_pop1;
            else
               s = adimat_pop1;
            end
            tmpba4 = adimat_pop1;
            if tmpba4 == 1
               x = adimat_pop1;
            end
         elseif tmpba3 == 3
            r = adimat_pop1;
            a_r = a_zeros1(r);
            [sa2 a] = adimat_pop;
         else
            error('Derivatives of matrix-p-norm not implemented yet.');
         end
      else
         error('Value is neither a matrix nor a vector!');
      end
   end
end
% $Id: adimat_norm2.m 4281 2014-05-21 09:23:04Z willkomm $
