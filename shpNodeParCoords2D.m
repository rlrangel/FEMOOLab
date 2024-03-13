function pt = shpNodeParCoords2D
% Returns the parametric coordinates of the nodal points of an element.
% Output arguments:
%  pt:   array with pairs of element nodal paramatric coordinates

 include_gblrefs;

 pt = zeros(2,nen);

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
  pt(1,1) = 0.0;   pt(2,1) = 0.0;
  pt(1,2) = 1.0;   pt(2,2) = 0.0;
  pt(1,3) = 0.0;   pt(2,3) = 1.0;

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
  pt(1,1) = -1.0;   pt(2,1) = -1.0;
  pt(1,2) =  1.0;   pt(2,2) = -1.0;
  pt(1,3) =  1.0;   pt(2,3) =  1.0;
  pt(1,4) = -1.0;   pt(2,4) =  1.0;

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
  pt(1,1) = 0.0;   pt(2,1) = 0.0;
  pt(1,2) = 1.0;   pt(2,2) = 0.0;
  pt(1,3) = 0.0;   pt(2,3) = 1.0;
  pt(1,4) = 0.5;   pt(2,4) = 0.0;
  pt(1,5) = 0.5;   pt(2,5) = 0.5;
  pt(1,6) = 0.0;   pt(2,6) = 0.5;

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
  pt(1,1) = -1.0;   pt(2,1) = -1.0;
  pt(1,2) =  1.0;   pt(2,2) = -1.0;
  pt(1,3) =  1.0;   pt(2,3) =  1.0;
  pt(1,4) = -1.0;   pt(2,4) =  1.0;
  pt(1,5) =  0.0;   pt(2,5) = -1.0;
  pt(1,6) =  1.0;   pt(2,6) =  0.0;
  pt(1,7) =  0.0;   pt(2,7) =  1.0;
  pt(1,8) = -1.0;   pt(2,8) =  0.0;
 end

end
