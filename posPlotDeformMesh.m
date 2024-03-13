function posPlotDeformMesh(color)
% Plots deformed mesh in current active figure.

include_gblrefs;

 XX = zeros(1,nen+1);
 YY = zeros(1,nen+1);

 % Get node resequecing vector in CCW order
 rv = shpNodeCCWOrder;

 % Display deformed mesh
 for e = 1:nel
  IENe = IEN(rv,e);
  for j = 1:nen
   XX(j) = x(IENe(j)) + deform_fac*udispl(IENe(j));
   YY(j) = y(IENe(j)) + deform_fac*vdispl(IENe(j));
  end
  XX(nen+1) = x(IENe(1)) + deform_fac*udispl(IENe(1));
  YY(nen+1) = y(IENe(1)) + deform_fac*vdispl(IENe(1));
  plot(XX,YY,color);
  hold on;
 end

end
