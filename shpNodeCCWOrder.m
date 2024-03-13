function rv = shpNodeCCWOrder
% Returns the resequence vector of nodal indices of an element in CCW order.
% Output arguments:
%  rv:   resequence vector of nodal indices of an element in CCW order

 include_gblrefs;

 rv = zeros(nen,1);

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
  rv(1) = 1;
  rv(2) = 2;
  rv(3) = 3;

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
  rv(1) = 1;
  rv(2) = 2;
  rv(3) = 3;
  rv(4) = 4;

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
  rv(1) = 1;
  rv(2) = 4;
  rv(3) = 2;
  rv(4) = 5;
  rv(5) = 3;
  rv(6) = 6;

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
  rv(1) = 1;
  rv(2) = 5;
  rv(3) = 2;
  rv(4) = 6;
  rv(5) = 3;
  rv(6) = 7;
  rv(7) = 4;
  rv(8) = 8;
 end

end
