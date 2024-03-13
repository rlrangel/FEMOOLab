function [prc,thetap] = posPrincStress(str)
% Get principal stress components and orientation for a given stress tensor.
% Input arguments:
%  str:  stress tensor (sx, sy, txy) stored in a column vector.
% Output arguments:
%  prc:  principal stress components (s1, s2, taumax) stored in a column vector.
%  thetap: angle of normal of principal stress plane w.r.t x axis (angle is 
%          returned in radians from 0 to 180 degrees).

 prc = zeros(3,1);

 sx  = str(1);
 sy  = str(2);
 txy = str(3);

 center = (sx + sy)*0.5;
 deltas = (sx - sy)*0.5;
 radius = sqrt((deltas^2) + (txy*txy));

 prc(1) = center + radius;   % s1
 prc(2) = center - radius;   % s2
 prc(3) = radius;            % taumax

 if( abs(deltas) > 0.0 )
  thetap = 0.5 * atan2(txy,deltas);
 elseif( txy > 0.0 )
  thetap = pi / 4.0;
 elseif( txy < 0.0 )
  thetap = -pi / 4.0;
 else
  thetap = 0.0;
 end

 % Transform angle from 0 to 180 degrees
 if( thetap < 0.0 )
  thetap = pi + thetap;

end
