%HEADER.ANALYSIS
'plane_conduction'

%HEADER.ANALYSIS.TYPE
'transient'

%HEADER.ANALYSIS.ALGORITHM
'crank_nicolson'

%HEADER.ANALYSIS.TIME.STEPS
0.01

%HEADER.ANALYSIS.MAXIMUM.STEPS
2000

%HEADER.ANALYSIS.PRINT.STEPS
20

%NODE
9

%NODE.COORD
9
1 0.000000	0.000000	0.000000
2 0.500000	0.000000	0.000000
3 1.000000	0.000000	0.000000
4 0.000000	0.500000	0.000000
5 0.500000	0.500000	0.000000
6 1.000000	0.500000	0.000000
7 0.000000	1.000000	0.000000
8 0.500000	1.000000	0.000000
9 1.000000	1.000000	0.000000

%MATERIAL
1

%MATERIAL.PROPERTY.DENSITY
1
1	1000.00

%MATERIAL.PROPERTY.THERMAL
1
1	1000.00	100.00

%THICKNESS
1
1	1.000000

%INTEGRATION.ORDER
1
1	2 2 1 2 2 1

%ELEMENT
4

%ELEMENT.Q4
4
1	1	1	1	1	2	5	4	
2	1	1	1	2	3	6	5	
3	1	1	1	4	5	8	7	
4	1	1	1	5	6	9	8	

%LOAD.CASE.NODAL.INITIAL.TEMPERATURE
9
1 300.0
2 300.0
3 300.0
4 300.0
5 300.0
6 300.0
7 300.0
8 300.0
9 300.0

%LOAD.CASE.NODAL.TEMPERATURE
8
1	400
2	400
3	400
4	400
6	400
7	400
8	400
9	400

%RESULT.CASE.STEP.NODAL.TEMPERATURE
1
5

%END
