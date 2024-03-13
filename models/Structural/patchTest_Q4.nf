%HEADER
Neutral file created by Mtool

%HEADER.ANALYSIS
'PLANE_STRESS'

%NODE
9

%NODE.COORD
9
1   0.0    10.0   0.000000
2   5.0    10.0   0.000000
3   10.0   10.0   0.000000
4   0.0    5.0    0.000000
5   0.0    0.0    0.000000
6   5.0    5.0    0.000000
7   10.0   5.0    0.000000
8   5.0    0.0    0.000000
9   10.0   0.0    0.000000

%NODE.SUPPORT
5
5   1 1 0 0 0 0
1   1 0 0 0 0 0
4   1 0 0 0 0 0
8   0 1 0 0 0 0
9   0 1 0 0 0 0

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
1   2 2 1 2 2 1

%ELEMENT
4

%ELEMENT.Q4
4
1   1 1 1   6 4 5 8
2   1 1 1   7 6 8 9
3   1 1 1   2 1 4 6
4   1 1 1   3 2 6 7

%LOAD.CASE.LINE.FORCE.UNIFORM
4
4   7 3   0   10.0 0.0  0.0
2   9 7   0   10.0 0.0  0.0
4   3 2   0   0.0  20.0 0.0
3   2 1   0   0.0  20.0 0.0

%END