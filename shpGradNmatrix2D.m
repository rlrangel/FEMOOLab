function GradN = shpGradNmatrix2D(r,s)
% Evaluates derivatives of shape functions (in parametric coordinates).
% Input arguments:
%  r,s:  parametric coordinate values
% Output arguments:
%  GradN: Shape function gradient matrix evaluated at given position

 include_gblrefs;

 GradN = zeros(2,nen);

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
  GradN(1,1) = -1.0;
  GradN(2,1) = -1.0;
  GradN(1,2) =  1.0;
  GradN(2,2) =  0.0;
  GradN(1,3) =  0.0;
  GradN(2,3) =  1.0;

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
  GradN(1,1) = -0.25*(1.0 - s);
  GradN(2,1) = -0.25*(1.0 - r);
  GradN(1,2) =  0.25*(1.0 - s);
  GradN(2,2) = -0.25*(1.0 + r);
  GradN(1,3) =  0.25*(1.0 + s);
  GradN(2,3) =  0.25*(1.0 + r);
  GradN(1,4) = -0.25*(1.0 + s);
  GradN(2,4) =  0.25*(1.0 - r);

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
  GradN(1,1) =  4.0*r + 4.0*s - 3.0;
  GradN(2,1) =  4.0*r + 4.0*s - 3.0;
  GradN(1,2) =  4.0*r - 1.0;
  GradN(2,2) =  0.0;
  GradN(1,3) =  0.0;
  GradN(2,3) =  4.0*s - 1.0;
  GradN(1,4) =  4.0 - 8.0*r - 4.0*s;
  GradN(2,4) = -4.0*r;
  GradN(1,5) =  4.0*s;
  GradN(2,5) =  4.0*r;
  GradN(1,6) = -4.0*s;
  GradN(2,6) =  4.0 - 4.0*r - 8.0*s;

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
  GradN(1,1) = (2.0*r - 2.0*r*s - s*s + s) / 4.0;
  GradN(2,1) = (2.0*s - r*r - 2.0*r*s + r) / 4.0;
  GradN(1,2) = (2.0*r - 2.0*r*s + s*s - s) / 4.0;
  GradN(2,2) = (2.0*s - r*r + 2.0*r*s - r) / 4.0;
  GradN(1,3) = (2.0*r + 2.0*r*s + s*s + s) / 4.0;
  GradN(2,3) = (2.0*s + r*r + 2.0*r*s + r) / 4.0;
  GradN(1,4) = (2.0*r + 2.0*r*s - s*s - s) / 4.0;
  GradN(2,4) = (2.0*s + r*r - 2.0*r*s - r) / 4.0;
  GradN(1,5) = r*s - r;
  GradN(2,5) = (r*r - 1.0) / 2.0;
  GradN(1,6) = (1.0 - s*s) / 2.0;
  GradN(2,6) = -s - r*s;
  GradN(1,7) = -r - r*s;
  GradN(2,7) = (1.0 - r*r) / 2.0;
  GradN(1,8) = (s*s - 1.0) / 2.0;
  GradN(2,8) = r*s - s;
 end

end
