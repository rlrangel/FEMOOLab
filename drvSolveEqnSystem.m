function [D,F] = drvSolveEqnSystem(neq,neqfree,K,F,D)
% Partition and solve the system of equations.
% Input arguments:
%  neq:       total number of equations
%  neqfree:   number of equations of free d.o.f.'s
%  K:         global stiffness matrix
% Input/Output arguments:
%  F:         global force vector (right hand side)
%  D:         global displacement vector


 % Partition the coefficient matrix K, forcing vector F and unknown vector D:
 %  f --> free or natural B.C. (unknown) degree-of-freedown (numbered first)
 %  c --> constrainted (essential - known) B.C. degree-of-freedown
 Kff = K(1:neqfree,1:neqfree);                    % Extract Kff matrix 
 Kfc = K(1:neqfree,neqfree+1:neq);                % Extract Kfc matrix
 Kcf = K(neqfree+1:neq,1:neqfree);                % Extract Kcf matrix
 Kcc = K(neqfree+1:neq,neqfree+1:neq);            % Extract Kcc matrix
 Ff  = F(1:neqfree);                              % Extract Ff vector
 Dc  = D(neqfree+1:neq);                          % Extract Dc vector
 Fc  = F(neqfree+1:neq);                          % Extract Fc vector
  
 % Solve for Df
 Df = Kff \ (Ff - Kfc*Dc);

 % Reconstruct the global unknown vector D
 D = [Df
      Dc];

 % Recover forcing unknown values (reactions) at essential B.C.
 % It is assumed that the Fc vector currently stores combined nodal
 % loads applyed directly to fixed d.o.f's.
 % Superimpose computed reaction values to combined nodal loads,
 % with inversed direction, that were applyed directly to fixed d.o.f.'s.
 Fc = -Fc + Kcf*Df + Kcc*Dc;

 % Reconstruct the global forcing vector 
 F = [Ff
      Fc];
 
end
