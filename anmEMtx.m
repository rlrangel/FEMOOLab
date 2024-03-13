function E = anmEMtx(e)
% Generate material constituive matrix for a given element.
% The material constitutive matrix is assembled using the 
% element modulus of elasticity and poisson ratio.
% The matrix size and layout depends on the type of analysis.
% For plane stress and plane strain analysis, the size is (3,3).
% For axisymmetric analysis, the size is (4,4).
% Input arguments:
%  e:   index of element

 include_gblrefs;

 E_coeff = Elasticity(emat(e));          % get element elasticity coeff.
 P_ratio = Poisson(emat(e));             % get element poisson ration.

 if(analysis_type == AXISYMMETRIC)

  E = zeros(4,4);

  E(1,1) = 1.0 - P_ratio;
  E(1,2) = P_ratio;
  E(1,3) = P_ratio;
  E(1,4) = 0.0;
  E(2,1) = P_ratio;
  E(2,2) = 1.0 - P_ratio;
  E(2,3) = P_ratio;
  E(2,4) = 0.0;
  E(3,1) = P_ratio;
  E(3,2) = P_ratio;
  E(3,3) = 1.0 - P_ratio;
  E(3,4) = 0.0;
  E(4,1) = 0.0;
  E(4,2) = 0.0;
  E(4,3) = 0.0;
  E(4,4) = (1.0 - (2.0*P_ratio))/2.0;

  E = E * (E_coeff / ((1.0 + P_ratio) * (1.0 - (2.0*P_ratio))));

 elseif(analysis_type == PLANE_STRAIN)

  E = zeros(3,3);

  E(1,1) = 1.0 - P_ratio;
  E(1,2) = P_ratio;
  E(1,3) = 0.0;
  E(2,1) = P_ratio;
  E(2,2) = 1.0 - P_ratio;
  E(2,3) = 0.0;
  E(3,1) = 0.0;
  E(3,2) = 0.0;
  E(3,3) = (1.0 - (2.0*P_ratio))/2.0;

  E = E * (E_coeff / ((1.0 + P_ratio) * (1.0 - (2.0*P_ratio))));

 else                                    % plane stress analysis

  E = zeros(3,3);

  E(1,1) = 1.0;
  E(1,2) = P_ratio;
  E(1,3) = 0.0;
  E(2,1) = P_ratio;
  E(2,2) = 1.0;
  E(2,3) = 0.0;
  E(3,1) = 0.0;
  E(3,2) = 0.0;
  E(3,3) = (1.0 - P_ratio)/2.0;

  E = E * (E_coeff / (1.0 - (P_ratio*P_ratio)));
 end

end
