%HEADER
Neutral file created by Mtool

%HEADER.ANALYSIS
'PLANE_STRESS'

%NODE
9

%NODE.COORD
9
1   10.0   0.0    0.0
2   10.0   5.0    0.0
3   5.0    0.0    0.0
4   10.0   10.0   0.0
5   5.0    5.0    0.0
6   0.0    0.0    0.0
7   5.0    10.0   0.0
8   0.0    5.0    0.0
9   0.0    10.0   0.0

%NODE.SUPPORT
5
6   1 1 0 0 0 0
3   0 1 0 0 0 0
1   0 1 0 0 0 0
8   1 0 0 0 0 0
9   1 0 0 0 0 0

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
8

%ELEMENT.T3
8
1   1 1 1   5 6 3
2   1 1 1   5 8 6
3   1 1 1   2 3 1
4   1 1 1   2 5 3
5   1 1 1   7 8 5
6   1 1 1   7 9 8
7   1 1 1   4 5 2
8   1 1 1   4 7 5

%LOAD.CASE.LINE.FORCE.UNIFORM
4
8   4 7   0   0.0  20.0 0.0
6   7 9   0   0.0  20.0 0.0
3   1 2   0   10.0 0.0  0.0
7   2 4   0   10.0 0.0  0.0

%END