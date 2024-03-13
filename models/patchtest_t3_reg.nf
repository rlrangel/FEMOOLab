%HEADER
Neutral file created by FEMOOP program

%HEADER.VERSION
1-000  -  Aug 31 2000


%HEADER.ANALYSIS
'plane_stress'

%NODE
9

%NODE.COORD
9
1	10.000000	0.000000	0.000000
2	10.000000	5.000000	0.000000
3	5.000000	0.000000	0.000000
4	10.000000	10.000000	0.000000
5	5.000000	5.000000	0.000000
6	0.000000	0.000000	0.000000
7	5.000000	10.000000	0.000000
8	0.000000	5.000000	0.000000
9	0.000000	10.000000	0.000000

%NODE.SUPPORT
5
6	1	1	0	0	0	0
3	0	1	0	0	0	0
1	0	1	0	0	0	0
8	1	0	0	0	0	0
9	1	0	0	0	0	0

%MATERIAL
1

%MATERIAL.LABEL
1
1	'mat'

%MATERIAL.ISOTROPIC
1
1	10000.000000	0.250000

%THICKNESS
1
1	1.000000

%INTEGRATION.ORDER
5
1	2	2	1	2	2	1
2	1	1	1	1	1	1
3	2	2	1	2	2	1
4	2	2	1	2	2	1
5	2	2	1	2	2	1

%ELEMENT
8

%ELEMENT.T3
8
1	1	1	2	5	6	3
2	1	1	2	5	8	6
3	1	1	2	2	3	1
4	1	1	2	2	5	3
5	1	1	2	7	8	5
6	1	1	2	7	9	8
7	1	1	2	4	5	2
8	1	1	2	4	7	5

%LOAD
1
1	'Load_Case_1'

%LOAD.CASE
1

%LOAD.CASE.LINE.FORCE.UNIFORM
4
8	4	7	0	0.000000	20.000000	0.000000
6	7	9	0	0.000000	20.000000	0.000000
3	1	2	0	10.000000	0.000000	0.000000
7	2	4	0	10.000000	0.000000	0.000000

%END