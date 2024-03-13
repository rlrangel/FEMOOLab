function [nne,IENe] = shpGetEdge(e,nod1,nod2)
% Get nodal incidence on a given 2D element edge (side).
% Input arguments:
%  e:      element number
%  nod1:   first node on element edge
%  nod2:   last node on element edge
% Output arguments:
%  nne:    number of nodes of target element edge
%  IENe:   local nodal connectivity array of target element edge

 include_gblrefs;

 nne = nedgen;                         % number of nodes of equiv.load vector
                                       % is equal to number of edge nodes

 IENe = zeros(nne,1);

 if(element_type == TRIA3)         % linear triangle element
 % 
 %                         3 +
 %                           |\
 %                           | \
 %                           |  \
 %                           |   \ TRIA3
 %                           |    \
 %                         1 +-----+ 2
 %
  if    ( nod1 == IEN(1,e) )
   if    ( nod2 == IEN(2,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
   elseif( nod2 == IEN(3,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
   end
  elseif( nod1 == IEN(2,e) )
   if    ( nod2 == IEN(3,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
   elseif( nod2 == IEN(1,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
   end
  elseif( nod1 == IEN(3,e) )
   if    ( nod2 == IEN(1,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
   elseif( nod2 == IEN(2,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
   end
  end

 elseif(element_type == QUAD4)     % bilinear quadrilateral element
 % 
 %                         4 +---------------+ 3
 %                           |               |
 %                           |               |
 %                           |               |
 %                           |     QUAD4     |
 %                           |               |
 %                           |               |
 %                         1 +---------------+ 2
 %
  if    ( nod1 == IEN(1,e) )
   if    ( nod2 == IEN(2,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
   elseif( nod2 == IEN(4,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
   end
  elseif( nod1 == IEN(2,e) )
   if    ( nod2 == IEN(3,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
   elseif( nod2 == IEN(1,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
   end
  elseif( nod1 == IEN(3,e) )
   if    ( nod2 == IEN(4,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
   elseif( nod2 == IEN(2,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
   end
  elseif( nod1 == IEN(4,e) )
   if    ( nod2 == IEN(1,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
   elseif( nod2 == IEN(3,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
   end
  end

 elseif(element_type == TRIA6)     % quadratic triangle element
 % 
 %                         3 +
 %                           |\
 %                           | \
 %                         6 +  + 5
 %                           |   \ TRIA6
 %                           |    \
 %                         1 +--+--+ 2
 %                              4
 %
  if    ( nod1 == IEN(1,e) )
   if    ( nod2 == IEN(2,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
    IENe(3) = IEN(4,e);
   elseif( nod2 == IEN(3,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
    IENe(3) = IEN(6,e);
   end
  elseif( nod1 == IEN(2,e) )
   if    ( nod2 == IEN(3,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
    IENe(3) = IEN(5,e);
   elseif( nod2 == IEN(1,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
    IENe(3) = IEN(4,e);
   end
  elseif( nod1 == IEN(3,e) )
   if    ( nod2 == IEN(1,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
    IENe(3) = IEN(6,e);
   elseif( nod2 == IEN(2,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
    IENe(3) = IEN(5,e);
   end
  end

 elseif(element_type == QUAD8)     % serendipity quadratic element
 % 
 %                                   7
 %                         4 +-------+-------+ 3
 %                           |               |
 %                           |               |
 %                         8 +     QUAD8     + 6
 %                           |               |
 %                           |               |
 %                         1 +-------+-------+ 2
 %                                   5
 %
  if    ( nod1 == IEN(1,e) )
   if    ( nod2 == IEN(2,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
    IENe(3) = IEN(5,e);
   elseif( nod2 == IEN(4,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
    IENe(3) = IEN(8,e);
   end
  elseif( nod1 == IEN(2,e) )
   if    ( nod2 == IEN(3,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
    IENe(3) = IEN(6,e);
   elseif( nod2 == IEN(1,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
    IENe(3) = IEN(5,e);
   end
  elseif( nod1 == IEN(3,e) )
   if    ( nod2 == IEN(4,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
    IENe(3) = IEN(7,e);
   elseif( nod2 == IEN(2,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
    IENe(3) = IEN(6,e);
   end
  elseif( nod1 == IEN(4,e) )
   if    ( nod2 == IEN(1,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
    IENe(3) = IEN(8,e);
   elseif( nod2 == IEN(3,e) )
    IENe(1) = nod1;
    IENe(2) = nod2;
    IENe(3) = IEN(7,e);
   end
  end

 end

end
