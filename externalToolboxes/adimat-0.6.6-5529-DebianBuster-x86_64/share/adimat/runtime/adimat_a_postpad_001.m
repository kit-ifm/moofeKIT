function a_c = adimat_a_postpad_001(a,b,c,a_z)
  dim = adimat_first_nonsingleton(a);
  [~, a_c] = adimat_a_postpad_1010(a,b,c,dim,a_z);
