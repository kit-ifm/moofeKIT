% Generated by ADiMat 0.6.6-5530 (00419e1f)
% Copyright 2009-2013 Johannes Willkomm, Fachgebiet Scientific Computing,
% TU Darmstadt, 64289 Darmstadt, Germany
% Copyright 2001-2008 Andre Vehreschild, Institute for Scientific Computing,
% RWTH Aachen University, 52056 Aachen, Germany.
% Visit us on the web at http://www.adimat.de
% Report bugs to johannes@johannes-willkomm.de
%
%
%                             DISCLAIMER
%
% ADiMat was prepared as part of an employment at the Institute
% for Scientific Computing, RWTH Aachen University, Germany and is
% provided AS IS. NEITHER THE AUTHOR(S), THE GOVERNMENT OF THE FEDERAL
% REPUBLIC OF GERMANY NOR ANY AGENCY THEREOF, NOR THE RWTH AACHEN UNIVERSITY,
% INCLUDING ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY,
% EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY
% FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION OR
% PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
% PRIVATELY OWNED RIGHTS.
%
% Global flags were:
% FORWARDMODE -- Apply the forward mode to the files.
% NOOPEROPTIM -- Do not use optimized operators. I.e.:
%		 g_a*b*g_c -/-> mtimes3(g_a, b, g_c)
% NOLOCALCSE  -- Do not use local common subexpression elimination when
%		 canonicalizing the code.
% NOGLOBALCSE -- Prevents the application of global common subexpression
%		 elimination after canonicalizing the code.
% NOPRESCALARFOLDING -- Switch off folding of scalar constants before
%		 augmentation.
% NOPOSTSCALARFOLDING -- Switch off folding of scalar constants after
%		 augmentation.
% NOCONSTFOLDMULT0 -- Switch off folding of product with one factor
%		 being zero: b*0=0.
% FUNCMODE    -- Inputfile is a function (This flag can not be set explicitly).
% NOTMPCLEAR  -- Suppress generation of clear g_* instructions.
% UNBOUND_XML  -- Write list of unbound identifiers in XML format.
% DEPENDENCIES_XML  -- Write list of functions in XML format.
% UNBOUND_ERROR	-- Stop with error if unbound identifiers found (default).
% FUNCTION_LIST_XML	-- Write list of functions to XML file.
% VERBOSITYLEVEL=5
% AD_IVARS= edN1, k, dimension, J, dNrAll, DMat, gaussWeight, I
% AD_DVARS= Re

function [g_Re, Re]= g_gaussDisplacementSCSaintVenantEndpoint(g_edN1, edN1, g_k, k, g_dimension, dimension, g_J, J, g_dNrAll, dNrAll, g_DMat, DMat, g_gaussWeight, gaussWeight, g_I, I)
   indx= dimension* k- (dimension- 1): dimension* k; 
   g_tmp_J_00000= g_J(: , indx);
   tmp_J_00000= J(: , indx);
   g_tmp_dNrAll_00000= g_dNrAll(indx, : );
   tmp_dNrAll_00000= dNrAll(indx, : );
   g_dNX= adimat_g_mldivide((g_tmp_J_00000' ), (tmp_J_00000' ), g_tmp_dNrAll_00000, tmp_dNrAll_00000);
   dNX= (tmp_J_00000' )\ tmp_dNrAll_00000; 
   numberOfNodes= size(dNX, 2); 
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00024= call(@reshape, g_edN1, dimension, numberOfNodes);
   tmp_gaussDisplacementSCSaintVenantEndpoint_00024= reshape(edN1, dimension, numberOfNodes); 
   % Update detected: edN1= some_expression(edN1,...)
   g_edN1= g_tmp_gaussDisplacementSCSaintVenantEndpoint_00024;
   edN1= tmp_gaussDisplacementSCSaintVenantEndpoint_00024;
   g_tmp_J_00001= g_J(: , indx);
   tmp_J_00001= J(: , indx);
   detJ= det(tmp_J_00001' ); 
   g_detJ= detJ* call(@trace, inv(tmp_J_00001' )* g_tmp_J_00001' );
   if detJ< 10* eps
      error('Jacobi determinant equal or less than zero.')
   end
   g_FAkt= g_edN1* dNX' + edN1* g_dNX' ;
   FAkt= edN1* dNX' ; 
   % B-matrix
   numberOfDisplacementDofs= numberOfNodes* dimension; 
   numberOfSymmetricVoigtDofs= dimension+ (dimension^ 2- dimension)/ 2; 
   BAkt= zeros(numberOfSymmetricVoigtDofs, numberOfDisplacementDofs); 
   g_BAkt= g_zeros(size(BAkt));
   g_tmp_FAkt_00000= g_FAkt(1, 1);
   tmp_FAkt_00000= FAkt(1, 1);
   g_tmp_dNX_00000= g_dNX(1, : );
   tmp_dNX_00000= dNX(1, : );
   g_BAkt(1, 1: dimension: end)= g_tmp_FAkt_00000* tmp_dNX_00000+ tmp_FAkt_00000* g_tmp_dNX_00000;
   BAkt(1, 1: dimension: end)= tmp_FAkt_00000* tmp_dNX_00000; 
   g_tmp_FAkt_00001= g_FAkt(2, 1);
   tmp_FAkt_00001= FAkt(2, 1);
   g_tmp_dNX_00001= g_dNX(1, : );
   tmp_dNX_00001= dNX(1, : );
   g_BAkt(1, 2: dimension: end)= g_tmp_FAkt_00001* tmp_dNX_00001+ tmp_FAkt_00001* g_tmp_dNX_00001;
   BAkt(1, 2: dimension: end)= tmp_FAkt_00001* tmp_dNX_00001; 
   g_tmp_FAkt_00002= g_FAkt(3, 1);
   tmp_FAkt_00002= FAkt(3, 1);
   g_tmp_dNX_00002= g_dNX(1, : );
   tmp_dNX_00002= dNX(1, : );
   g_BAkt(1, 3: dimension: end)= g_tmp_FAkt_00002* tmp_dNX_00002+ tmp_FAkt_00002* g_tmp_dNX_00002;
   BAkt(1, 3: dimension: end)= tmp_FAkt_00002* tmp_dNX_00002; 
   g_tmp_FAkt_00003= g_FAkt(1, 2);
   tmp_FAkt_00003= FAkt(1, 2);
   g_tmp_dNX_00003= g_dNX(2, : );
   tmp_dNX_00003= dNX(2, : );
   g_BAkt(2, 1: dimension: end)= g_tmp_FAkt_00003* tmp_dNX_00003+ tmp_FAkt_00003* g_tmp_dNX_00003;
   BAkt(2, 1: dimension: end)= tmp_FAkt_00003* tmp_dNX_00003; 
   g_tmp_FAkt_00004= g_FAkt(2, 2);
   tmp_FAkt_00004= FAkt(2, 2);
   g_tmp_dNX_00004= g_dNX(2, : );
   tmp_dNX_00004= dNX(2, : );
   g_BAkt(2, 2: dimension: end)= g_tmp_FAkt_00004* tmp_dNX_00004+ tmp_FAkt_00004* g_tmp_dNX_00004;
   BAkt(2, 2: dimension: end)= tmp_FAkt_00004* tmp_dNX_00004; 
   g_tmp_FAkt_00005= g_FAkt(3, 2);
   tmp_FAkt_00005= FAkt(3, 2);
   g_tmp_dNX_00005= g_dNX(2, : );
   tmp_dNX_00005= dNX(2, : );
   g_BAkt(2, 3: dimension: end)= g_tmp_FAkt_00005* tmp_dNX_00005+ tmp_FAkt_00005* g_tmp_dNX_00005;
   BAkt(2, 3: dimension: end)= tmp_FAkt_00005* tmp_dNX_00005; 
   g_tmp_FAkt_00006= g_FAkt(1, 3);
   tmp_FAkt_00006= FAkt(1, 3);
   g_tmp_dNX_00006= g_dNX(3, : );
   tmp_dNX_00006= dNX(3, : );
   g_BAkt(3, 1: dimension: end)= g_tmp_FAkt_00006* tmp_dNX_00006+ tmp_FAkt_00006* g_tmp_dNX_00006;
   BAkt(3, 1: dimension: end)= tmp_FAkt_00006* tmp_dNX_00006; 
   g_tmp_FAkt_00007= g_FAkt(2, 3);
   tmp_FAkt_00007= FAkt(2, 3);
   g_tmp_dNX_00007= g_dNX(3, : );
   tmp_dNX_00007= dNX(3, : );
   g_BAkt(3, 2: dimension: end)= g_tmp_FAkt_00007* tmp_dNX_00007+ tmp_FAkt_00007* g_tmp_dNX_00007;
   BAkt(3, 2: dimension: end)= tmp_FAkt_00007* tmp_dNX_00007; 
   g_tmp_FAkt_00008= g_FAkt(3, 3);
   tmp_FAkt_00008= FAkt(3, 3);
   g_tmp_dNX_00008= g_dNX(3, : );
   tmp_dNX_00008= dNX(3, : );
   g_BAkt(3, 3: dimension: end)= g_tmp_FAkt_00008* tmp_dNX_00008+ tmp_FAkt_00008* g_tmp_dNX_00008;
   BAkt(3, 3: dimension: end)= tmp_FAkt_00008* tmp_dNX_00008; 
   g_tmp_FAkt_00009= g_FAkt(1, 1);
   tmp_FAkt_00009= FAkt(1, 1);
   g_tmp_dNX_00009= g_dNX(2, : );
   tmp_dNX_00009= dNX(2, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00000= g_tmp_FAkt_00009* tmp_dNX_00009+ tmp_FAkt_00009* g_tmp_dNX_00009;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00000= tmp_FAkt_00009* tmp_dNX_00009;
   g_tmp_FAkt_00010= g_FAkt(1, 2);
   tmp_FAkt_00010= FAkt(1, 2);
   g_tmp_dNX_00010= g_dNX(1, : );
   tmp_dNX_00010= dNX(1, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00001= g_tmp_FAkt_00010* tmp_dNX_00010+ tmp_FAkt_00010* g_tmp_dNX_00010;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00001= tmp_FAkt_00010* tmp_dNX_00010;
   g_BAkt(4, 1: dimension: end)= g_tmp_gaussDisplacementSCSaintVenantEndpoint_00000+ g_tmp_gaussDisplacementSCSaintVenantEndpoint_00001;
   BAkt(4, 1: dimension: end)= tmp_gaussDisplacementSCSaintVenantEndpoint_00000+ tmp_gaussDisplacementSCSaintVenantEndpoint_00001; 
   g_tmp_FAkt_00011= g_FAkt(2, 1);
   tmp_FAkt_00011= FAkt(2, 1);
   g_tmp_dNX_00011= g_dNX(2, : );
   tmp_dNX_00011= dNX(2, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00002= g_tmp_FAkt_00011* tmp_dNX_00011+ tmp_FAkt_00011* g_tmp_dNX_00011;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00002= tmp_FAkt_00011* tmp_dNX_00011;
   g_tmp_FAkt_00012= g_FAkt(2, 2);
   tmp_FAkt_00012= FAkt(2, 2);
   g_tmp_dNX_00012= g_dNX(1, : );
   tmp_dNX_00012= dNX(1, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00003= g_tmp_FAkt_00012* tmp_dNX_00012+ tmp_FAkt_00012* g_tmp_dNX_00012;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00003= tmp_FAkt_00012* tmp_dNX_00012;
   g_BAkt(4, 2: dimension: end)= g_tmp_gaussDisplacementSCSaintVenantEndpoint_00002+ g_tmp_gaussDisplacementSCSaintVenantEndpoint_00003;
   BAkt(4, 2: dimension: end)= tmp_gaussDisplacementSCSaintVenantEndpoint_00002+ tmp_gaussDisplacementSCSaintVenantEndpoint_00003; 
   g_tmp_FAkt_00013= g_FAkt(3, 1);
   tmp_FAkt_00013= FAkt(3, 1);
   g_tmp_dNX_00013= g_dNX(2, : );
   tmp_dNX_00013= dNX(2, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00004= g_tmp_FAkt_00013* tmp_dNX_00013+ tmp_FAkt_00013* g_tmp_dNX_00013;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00004= tmp_FAkt_00013* tmp_dNX_00013;
   g_tmp_FAkt_00014= g_FAkt(3, 2);
   tmp_FAkt_00014= FAkt(3, 2);
   g_tmp_dNX_00014= g_dNX(1, : );
   tmp_dNX_00014= dNX(1, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00005= g_tmp_FAkt_00014* tmp_dNX_00014+ tmp_FAkt_00014* g_tmp_dNX_00014;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00005= tmp_FAkt_00014* tmp_dNX_00014;
   g_BAkt(4, 3: dimension: end)= g_tmp_gaussDisplacementSCSaintVenantEndpoint_00004+ g_tmp_gaussDisplacementSCSaintVenantEndpoint_00005;
   BAkt(4, 3: dimension: end)= tmp_gaussDisplacementSCSaintVenantEndpoint_00004+ tmp_gaussDisplacementSCSaintVenantEndpoint_00005; 
   g_tmp_FAkt_00015= g_FAkt(1, 2);
   tmp_FAkt_00015= FAkt(1, 2);
   g_tmp_dNX_00015= g_dNX(3, : );
   tmp_dNX_00015= dNX(3, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00006= g_tmp_FAkt_00015* tmp_dNX_00015+ tmp_FAkt_00015* g_tmp_dNX_00015;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00006= tmp_FAkt_00015* tmp_dNX_00015;
   g_tmp_FAkt_00016= g_FAkt(1, 3);
   tmp_FAkt_00016= FAkt(1, 3);
   g_tmp_dNX_00016= g_dNX(2, : );
   tmp_dNX_00016= dNX(2, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00007= g_tmp_FAkt_00016* tmp_dNX_00016+ tmp_FAkt_00016* g_tmp_dNX_00016;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00007= tmp_FAkt_00016* tmp_dNX_00016;
   g_BAkt(5, 1: dimension: end)= g_tmp_gaussDisplacementSCSaintVenantEndpoint_00006+ g_tmp_gaussDisplacementSCSaintVenantEndpoint_00007;
   BAkt(5, 1: dimension: end)= tmp_gaussDisplacementSCSaintVenantEndpoint_00006+ tmp_gaussDisplacementSCSaintVenantEndpoint_00007; 
   g_tmp_FAkt_00017= g_FAkt(2, 2);
   tmp_FAkt_00017= FAkt(2, 2);
   g_tmp_dNX_00017= g_dNX(3, : );
   tmp_dNX_00017= dNX(3, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00008= g_tmp_FAkt_00017* tmp_dNX_00017+ tmp_FAkt_00017* g_tmp_dNX_00017;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00008= tmp_FAkt_00017* tmp_dNX_00017;
   g_tmp_FAkt_00018= g_FAkt(2, 3);
   tmp_FAkt_00018= FAkt(2, 3);
   g_tmp_dNX_00018= g_dNX(2, : );
   tmp_dNX_00018= dNX(2, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00009= g_tmp_FAkt_00018* tmp_dNX_00018+ tmp_FAkt_00018* g_tmp_dNX_00018;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00009= tmp_FAkt_00018* tmp_dNX_00018;
   g_BAkt(5, 2: dimension: end)= g_tmp_gaussDisplacementSCSaintVenantEndpoint_00008+ g_tmp_gaussDisplacementSCSaintVenantEndpoint_00009;
   BAkt(5, 2: dimension: end)= tmp_gaussDisplacementSCSaintVenantEndpoint_00008+ tmp_gaussDisplacementSCSaintVenantEndpoint_00009; 
   g_tmp_FAkt_00019= g_FAkt(3, 2);
   tmp_FAkt_00019= FAkt(3, 2);
   g_tmp_dNX_00019= g_dNX(3, : );
   tmp_dNX_00019= dNX(3, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00010= g_tmp_FAkt_00019* tmp_dNX_00019+ tmp_FAkt_00019* g_tmp_dNX_00019;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00010= tmp_FAkt_00019* tmp_dNX_00019;
   g_tmp_FAkt_00020= g_FAkt(3, 3);
   tmp_FAkt_00020= FAkt(3, 3);
   g_tmp_dNX_00020= g_dNX(2, : );
   tmp_dNX_00020= dNX(2, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00011= g_tmp_FAkt_00020* tmp_dNX_00020+ tmp_FAkt_00020* g_tmp_dNX_00020;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00011= tmp_FAkt_00020* tmp_dNX_00020;
   g_BAkt(5, 3: dimension: end)= g_tmp_gaussDisplacementSCSaintVenantEndpoint_00010+ g_tmp_gaussDisplacementSCSaintVenantEndpoint_00011;
   BAkt(5, 3: dimension: end)= tmp_gaussDisplacementSCSaintVenantEndpoint_00010+ tmp_gaussDisplacementSCSaintVenantEndpoint_00011; 
   g_tmp_FAkt_00021= g_FAkt(1, 1);
   tmp_FAkt_00021= FAkt(1, 1);
   g_tmp_dNX_00021= g_dNX(3, : );
   tmp_dNX_00021= dNX(3, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00012= g_tmp_FAkt_00021* tmp_dNX_00021+ tmp_FAkt_00021* g_tmp_dNX_00021;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00012= tmp_FAkt_00021* tmp_dNX_00021;
   g_tmp_FAkt_00022= g_FAkt(1, 3);
   tmp_FAkt_00022= FAkt(1, 3);
   g_tmp_dNX_00022= g_dNX(1, : );
   tmp_dNX_00022= dNX(1, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00013= g_tmp_FAkt_00022* tmp_dNX_00022+ tmp_FAkt_00022* g_tmp_dNX_00022;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00013= tmp_FAkt_00022* tmp_dNX_00022;
   g_BAkt(6, 1: dimension: end)= g_tmp_gaussDisplacementSCSaintVenantEndpoint_00012+ g_tmp_gaussDisplacementSCSaintVenantEndpoint_00013;
   BAkt(6, 1: dimension: end)= tmp_gaussDisplacementSCSaintVenantEndpoint_00012+ tmp_gaussDisplacementSCSaintVenantEndpoint_00013; 
   g_tmp_FAkt_00023= g_FAkt(2, 1);
   tmp_FAkt_00023= FAkt(2, 1);
   g_tmp_dNX_00023= g_dNX(3, : );
   tmp_dNX_00023= dNX(3, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00014= g_tmp_FAkt_00023* tmp_dNX_00023+ tmp_FAkt_00023* g_tmp_dNX_00023;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00014= tmp_FAkt_00023* tmp_dNX_00023;
   g_tmp_FAkt_00024= g_FAkt(2, 3);
   tmp_FAkt_00024= FAkt(2, 3);
   g_tmp_dNX_00024= g_dNX(1, : );
   tmp_dNX_00024= dNX(1, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00015= g_tmp_FAkt_00024* tmp_dNX_00024+ tmp_FAkt_00024* g_tmp_dNX_00024;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00015= tmp_FAkt_00024* tmp_dNX_00024;
   g_BAkt(6, 2: dimension: end)= g_tmp_gaussDisplacementSCSaintVenantEndpoint_00014+ g_tmp_gaussDisplacementSCSaintVenantEndpoint_00015;
   BAkt(6, 2: dimension: end)= tmp_gaussDisplacementSCSaintVenantEndpoint_00014+ tmp_gaussDisplacementSCSaintVenantEndpoint_00015; 
   g_tmp_FAkt_00025= g_FAkt(3, 1);
   tmp_FAkt_00025= FAkt(3, 1);
   g_tmp_dNX_00025= g_dNX(3, : );
   tmp_dNX_00025= dNX(3, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00016= g_tmp_FAkt_00025* tmp_dNX_00025+ tmp_FAkt_00025* g_tmp_dNX_00025;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00016= tmp_FAkt_00025* tmp_dNX_00025;
   g_tmp_FAkt_00026= g_FAkt(3, 3);
   tmp_FAkt_00026= FAkt(3, 3);
   g_tmp_dNX_00026= g_dNX(1, : );
   tmp_dNX_00026= dNX(1, : );
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00017= g_tmp_FAkt_00026* tmp_dNX_00026+ tmp_FAkt_00026* g_tmp_dNX_00026;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00017= tmp_FAkt_00026* tmp_dNX_00026;
   g_BAkt(6, 3: dimension: end)= g_tmp_gaussDisplacementSCSaintVenantEndpoint_00016+ g_tmp_gaussDisplacementSCSaintVenantEndpoint_00017;
   BAkt(6, 3: dimension: end)= tmp_gaussDisplacementSCSaintVenantEndpoint_00016+ tmp_gaussDisplacementSCSaintVenantEndpoint_00017; % Cauchy-Green tensor
   g_CAkt= g_FAkt' * FAkt+ FAkt' * g_FAkt;
   CAkt= FAkt' * FAkt; 
   % Green-Lagrange tensor
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00018= g_CAkt- g_I;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00018= CAkt- I;
   g_En1= 0.5* g_tmp_gaussDisplacementSCSaintVenantEndpoint_00018;
   En1= 0.5* tmp_gaussDisplacementSCSaintVenantEndpoint_00018; 
   g_tmp_En1_00000= g_En1(1, 1);
   tmp_En1_00000= En1(1, 1);
   g_tmp_En1_00001= g_En1(2, 2);
   tmp_En1_00001= En1(2, 2);
   g_tmp_En1_00002= g_En1(3, 3);
   tmp_En1_00002= En1(3, 3);
   g_tmp_En1_00003= g_En1(1, 2);
   tmp_En1_00003= En1(1, 2);
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00019= 2* g_tmp_En1_00003;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00019= 2* tmp_En1_00003;
   g_tmp_En1_00004= g_En1(3, 2);
   tmp_En1_00004= En1(3, 2);
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00020= 2* g_tmp_En1_00004;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00020= 2* tmp_En1_00004;
   g_tmp_En1_00005= g_En1(3, 1);
   tmp_En1_00005= En1(3, 1);
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00021= 2* g_tmp_En1_00005;
   tmp_gaussDisplacementSCSaintVenantEndpoint_00021= 2* tmp_En1_00005;
   g_tmp_gaussDisplacementSCSaintVenantEndpoint_00022= [g_tmp_En1_00000, g_tmp_En1_00001, g_tmp_En1_00002, g_tmp_gaussDisplacementSCSaintVenantEndpoint_00019, g_tmp_gaussDisplacementSCSaintVenantEndpoint_00020, g_tmp_gaussDisplacementSCSaintVenantEndpoint_00021];
   tmp_gaussDisplacementSCSaintVenantEndpoint_00022= [tmp_En1_00000, tmp_En1_00001, tmp_En1_00002, tmp_gaussDisplacementSCSaintVenantEndpoint_00019, tmp_gaussDisplacementSCSaintVenantEndpoint_00020, tmp_gaussDisplacementSCSaintVenantEndpoint_00021];
   g_En1_v= g_tmp_gaussDisplacementSCSaintVenantEndpoint_00022' ;
   En1_v= tmp_gaussDisplacementSCSaintVenantEndpoint_00022' ; 
   % Stresses
   tmp_gaussDisplacementSCSaintVenantEndpoint_00023= 1/ 2;
   g_DW1_v= tmp_gaussDisplacementSCSaintVenantEndpoint_00023* g_DMat* En1_v+ tmp_gaussDisplacementSCSaintVenantEndpoint_00023* DMat* g_En1_v;
   DW1_v= tmp_gaussDisplacementSCSaintVenantEndpoint_00023* DMat* En1_v; 
   % Residual
   g_tmp_gaussWeight_00000= g_gaussWeight(k);
   tmp_gaussWeight_00000= gaussWeight(k);
   g_Re= 2* g_BAkt' * DW1_v* detJ* tmp_gaussWeight_00000+ 2* BAkt' * g_DW1_v* detJ* tmp_gaussWeight_00000+ 2* BAkt' * DW1_v* g_detJ* tmp_gaussWeight_00000+ 2* BAkt' * DW1_v* detJ* g_tmp_gaussWeight_00000;
   Re= 2* BAkt' * DW1_v* detJ* tmp_gaussWeight_00000; 
   % Tangent
   % DW1 = [DW1_v(1) DW1_v(4) DW1_v(6); DW1_v(4) DW1_v(2) DW1_v(5); DW1_v(6) DW1_v(5) DW1_v(3)];
   % D2W1 = 1/4*DMat;
   % A1 = 2*dNX'*DW1*dNX*detJ*gaussWeight(k);
   % MAT = zeros(numberOfDOFs);
   % for g = 1:dimension
   %     MAT(g:dimension:numberOfDOFs,g:dimension:numberOfDOFs) = A1;
   % end
   % Ke = 4*BAkt'*D2W1*BAkt*detJ*gaussWeight(k) + MAT;
end