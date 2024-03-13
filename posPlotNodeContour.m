function posPlotNodeContour
% Plots current active node contour response in current active figure.

include_gblrefs;

 XX = zeros(1,nen+1);
 YY = zeros(1,nen+1);

 % Get node resequecing vector in CCW order
 rv = shpNodeCCWOrder;

 % Display current active node contour response

 switch response_type
  case DX_CONTOUR
   contour = udispl;
  case DY_CONTOUR
   contour = vdispl;
  case SX_CONTOUR
   contour = sx_nodeextrap;
  case SY_CONTOUR
   contour = sy_nodeextrap;
  case TXY_CONTOUR
   contour = txy_nodeextrap;
  case S1_CONTOUR
   contour = s1_nodeextrap;
  case S2_CONTOUR
   contour = s2_nodeextrap;
  case TMAX_CONTOUR
   contour = tmax_nodeextrap;
 end

 for e = 1:nel
  IENe = IEN(rv,e);
  for j = 1:nen
   XX(j) = x(IENe(j));
   YY(j) = y(IENe(j));
   ZZ(j) = contour(IENe(j));
  end
  XX(nen+1) = x(IENe(1));
  YY(nen+1) = y(IENe(1));
  ZZ(nen+1) = contour(IENe(1));
  patch(XX,YY,ZZ);
  hold on;
 end

end
