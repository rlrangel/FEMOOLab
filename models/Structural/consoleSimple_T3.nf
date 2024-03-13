%HEADER
Neutral file created by Mtool

%HEADER.ANALYSIS
'PLANE_STRESS'

%HEADER.ANALYSIS.METHOD
'ISOPARAMETRIC'

%NODE
4

%NODE.COORD
4
1   10.0   0.0    0.0
2   0.0    0.0    0.0
3   10.0   10.0   0.0
4   0.0    10.0   0.0

%NODE.SUPPORT
2
2   1 1 0 0 0 0
4   1 0 0 0 0 0

%MATERIAL
1

%MATERIAL.ISOTROPIC
1
1   10000.0   0.25

%THICKNESS
1
1   1.0

%INTEGRATION.ORDER
1
1   1 1 1 1 1 1

%ELEMENT
2

%ELEMENT.T3
2
1   1 1 1   3 2 1	
2   1 1 1   4 2 3	

%LOAD.CASE.NODAL.FORCES
2
3   0.0 -10.0 0.0 0.0 0.0 0.0
1   0.0 -10.0 0.0 0.0 0.0 0.0

%END