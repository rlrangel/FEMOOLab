function [ngp,w,gp] = gauss(type,order)
% Get gauss points coordinates in the parametric space of element and corresponding weights for a given quadrature type and order.
% Input arguments:
%  type:  gauss quadrature type (either line, triangle or quadrilateral)
%  order: quadrature order
% Output arguments:
%  ngp:   number of gauss points
%  w:     array of gauss point weights
%  gp:    array of gauss point parametric coordinates

 mlock;                   % prevent clearing this function from memory

 include_gblrefs;

 persistent qpts;
 persistent qwgt;

 if isempty(qpts)
  qpts = [.000000000000000 -.577350269189626 -.774596669241483 -.861136311594053
          .000000000000000  .577350269189626  .000000000000000 -.339981043584856
          .000000000000000  .000000000000000  .774596669241483  .339981043584856
          .000000000000000  .000000000000000  .000000000000000  .861136311594053];
 end

 if isempty(qwgt)
  qwgt = [2.00000000000000 1.000000000000000 0.555555555555555 0.347854845137454
          0.00000000000000 1.000000000000000 0.888888888888888 0.652145154862546
          0.00000000000000 0.000000000000000 0.555555555555555 0.652145154862546
          0.00000000000000 0.000000000000000 0.000000000000000 0.347854845137454];
 end

 % Initialize arrays of gauss point parametric coordenates and weights
 if(type == LINE_QUADRATURE)
  ngp = order;
 elseif(type == TRIA_QUADRATURE)
  if( order == 1 )
   ngp = 1;
  elseif( order == 2 || order == 3 )
   ngp = 3;
  elseif( order == 4 )
   ngp = 4;
  else
   ngp = 0;
  end
 elseif(type == QUAD_QUADRATURE)
  ngp = order*order;
 else
  ngp = 0;
 end

 if(ngp > 0)
  gp = zeros(2,ngp);
  w = zeros(1,ngp);
 end

 % Assemble arrays of gauss point coordinates and weights
 if(type == LINE_QUADRATURE)

  for i = 1:order
   gp(1,i) = qpts(i,order);
   gp(2,i) = 0.0;
   w(i)    = qwgt(i,order);
  end

 elseif(type == TRIA_QUADRATURE)

  if( order == 1 )
   gp(1,1) =  0.333333333333;
   gp(2,1) =  0.333333333333;
   w(1)    =  0.500000000000;
  elseif( order == 2 || order == 3 )
   gp(1,1) =  0.166666666666;
   gp(2,1) =  0.166666666666;
   w(1)    =  0.166666666666;
   gp(1,2) =  0.666666666666;
   gp(2,2) =  0.166666666666;
   w(2)    =  0.166666666666;
   gp(1,3) =  0.166666666666;
   gp(2,3) =  0.666666666666;
   w(3)    =  0.166666666666;
  elseif( order == 4 )
   gp(1,1) =  0.333333333333;
   gp(2,1) =  0.333333333333;
   w(1)    = -0.281250000000;
   gp(1,2) =  0.200000000000;
   gp(2,2) =  0.200000000000;
   w(2)    =  0.260416666666;
   gp(1,3) =  0.600000000000;
   gp(2,3) =  0.200000000000;
   w(3)    =  0.260416666666;
   gp(1,4) =  0.200000000000;
   gp(2,4) =  0.600000000000;
   w(4)    =  0.260416666666;
  end

 elseif(type == QUAD_QUADRATURE)

  npnts = 0;
  for i = 1:order
   for j = 1:order
    npnts = npnts + 1;
    gp(1,npnts) = qpts(i,order);
    gp(2,npnts) = qpts(j,order);
    w(npnts)    = qwgt(j,order) * qwgt(i,order);
   end
  end

 end

end
