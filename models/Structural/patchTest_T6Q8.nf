%HEADER
Neutral file created by Mtool

%HEADER.ANALYSIS
'PLANE_STRESS'

%NODE
23

%NODE.COORD
23
1    7.5    0.0    0.0
2    10.0   0.0    0.0
3    10.0   2.5    0.0
4    7.5    2.5    0.0
5    5.0    0.0    0.0
6    5.0    2.5    0.0
7    0.0    0.0    0.0
8    0.0    2.5    0.0
9    2.5    0.0    0.0
10   10.0   5.0    0.0
11   7.5    5.0    0.0
12   10.0   7.5    0.0
13   5.0    5.0    0.0
14   2.5    5.0    0.0
15   0.0    5.0    0.0
16   5.0    10.0   0.0
17   10.0   10.0   0.0
18   7.5    7.5    0.0
19   2.5    10.0   0.0
20   0.0    7.50   0.0
21   0.0    10.0   0.0
22   7.5    10.0   0.0
23   5.0    7.5    0.0

%NODE.SUPPORT
9
8    1 0 0 0 0 0
15   1 0 0 0 0 0
20   1 0 0 0 0 0
21   1 0 0 0 0 0
1    0 1 0 0 0 0
2    0 1 0 0 0 0
5    0 1 0 0 0 0
9    0 1 0 0 0 0
7    1 1 0 0 0 0

%MATERIAL
1

%MATERIAL.ISOTROPIC
1
1   10000.0   0.25

%THICKNESS
1
1   1.0

%INTEGRATION.ORDER
2
1   3 3 1 3 3 1
2   2 2 1 2 2 1

%ELEMENT
6

%ELEMENT.T6
4
1   1 1 1   5  1  2  3  10 4	
3   1 1 1   5  4  10 11 13 6	
4   1 1 1   13 11 10 12 17 18	
5   1 1 1   13 18 17 22 16 23	

%ELEMENT.Q8
2
2   1 1 2   13 14 15 8  7  9  5  6	
6   1 1 2   16 19 21 20 15 14 13 23	

%LOAD.CASE.LINE.FORCE.UNIFORM
4
1   2  10   0   10.0 0.0  0.0
4   10 17   0   10.0 0.0  0.0
5   17 16   0   0.0  20.0 0.0
6   16 21   0   0.0  20.0 0.0

%END