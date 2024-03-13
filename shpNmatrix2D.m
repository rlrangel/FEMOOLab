function N = shpNmatrix2D(r,s)
% Evaluates 2D shape function matrix.
% Evaluates shape functions (in parametric coordinates) at point (r,s)
% for 2D elements.
% Input arguments:
%  r:    first parametric coordinate value
%  s:    second parametric coordinate value
% Output arguments:
%  N:    shape function matrix evaluated at given position

 include_gblrefs;

 N = zeros(1,nen);

 if(element_type == TRIA3)         % linear triangle element
 % 
 %                           s
 %                           ^
 %                           |
 %     
 %                         3 +
 %                           |\
 %                           | \
 %                           |  \
 %                           |   \ TRIA3
 %                           |    \
 %                         1 +-----+ 2  ----> r
 %
  N(1) = 1.0 - r - s;
  N(2) = r;
  N(3) = s;

 elseif(element_type == QUAD4)     % bilinear quadrilateral element
 % 
 %                         4 +---------------+ 3
 %                           |       s       |
 %                           |       ^       |
 %                           |       |       |
 %                           |       ----> r |
 %                           |               |
 %                           |     QUAD4     |
 %                           |               |
 %                         1 +---------------+ 2
 %
 %
  N(1) = 0.25*(1.0-r)*(1.0-s);
  N(2) = 0.25*(1.0+r)*(1.0-s);
  N(3) = 0.25*(1.0+r)*(1.0+s);
  N(4) = 0.25*(1.0-r)*(1.0+s);

 elseif(element_type == TRIA6)     % quadratic triangle element
 % 
 %                           s
 %                           ^
 %                           |
 %     
 %                         3 +
 %                           |\
 %                           | \
 %                         6 +  + 5
 %                           |   \ TRIA6
 %                           |    \
 %                         1 +--+--+ 2  ----> r
 %                              4
 %
  N(1) = 1.0 - 3.0*r - 3.0*s + 4.0*r*s + 2.0*r*r + 2.0*s*s;
  N(2) = 2.0*r*r - r;
  N(3) = 2.0*s*s - s;
  N(4) = 4.0*r - 4.0*r*r - 4.0*r*s;
  N(5) = 4.0*r*s;
  N(6) = 4.0*s - 4.0*r*s - 4.0*s*s;

 elseif(element_type == QUAD8)     % serendipity quadratic quadrilateral element
 % 
 %                                   7
 %                         4 +-------+-------+ 3
 %                           |       s       |
 %                           |       ^       |
 %                           |       |       |
 %                         8 +       ----> r + 6
 %                           |               |
 %                           |     QUAD8     |
 %                           |               |
 %                         1 +-------+-------+ 2
 %                                   5
 %
  N(5) = 0.50*(1.0-r*r)*(1.0-s);
  N(6) = 0.50*(1.0+r)*(1.0-s*s);
  N(7) = 0.50*(1.0-r*r)*(1.0+s);
  N(8) = 0.50*(1.0-r)*(1.0-s*s);
  N(1) = 0.25*(1.0-r)*(1.0-s) - 0.50*N(8) - 0.50*N(5);
  N(2) = 0.25*(1.0+r)*(1.0-s) - 0.50*N(5) - 0.50*N(6);
  N(3) = 0.25*(1.0+r)*(1.0+s) - 0.50*N(6) - 0.50*N(7);
  N(4) = 0.25*(1.0-r)*(1.0+s) - 0.50*N(7) - 0.50*N(8);
 end

end
