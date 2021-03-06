% Main form-finding for a 3 strut prism known to be in equilibrium state
% Based on Yamamoto, Ohsaki, Fund and Gan papers

n = 6;
CONN = [1 1 2; ...
    2 3 4; ...
    3 5 6; ...
    4 1 3; ...
    5 1 5; ...
    6 1 6; ...
    7 2 3; ...
    8 2 4; ...
    9 2 5; ...
    10 3 6; ...
    11 4 6; ...
    12 4 5;];

q = [-1;-1;-1;1;1;1;1;1;1;1;1;1;];

coordinates = [159.42 138.82 0; ...
    -56.91 31.85 220.33; ...
    -45.79 206.23 0; ...
    155.17 72.02 219.99; ...
    -1.56 -5.19 0; ...
    14.35 235.56 220.39;];



% Second Truth Case - Gan et al
%{
n = 8;
CONN = [1 1 2; ...
    2 2 3; ...
    3 3 4;...
    4 1 4;...
    5 5 6; ...
    6 6 7; ...
    7 7 8; ...
    8 5 8; ...
    9 1 5; ...
    10 2 6; ...
    11 3 7; ...
    12 4 8; ...
    13 1 6; ...
    14 2 7; ...
    15 3 8; ...
    16 4 5;];
    

q = [1;1;1;1;1;1;1;1;1;1;1;1;-8;-3;-15;-7];

coordinates = [ -0.04699 -0.612765 -.3928; ...
0.100338 -.267531 .369638; ...
.462184 -.15565 .368403; ...
.502951 .101945 -.184040; ...
-.638767 -.233056 -.026372; ...
-.3347 .331275 .480091; ...
-.29177 .470376 -.063288; ...
-.015839 .364918 -.55163];

scaled_coordinates=[ 0 0 0; ...
    1437.34 -689.33 1688.4; ...
    2529.21 0 0;...
    2015.31 1808.41 0;...
    2261.04 -1106.38 1220.02; ...
    2062.49 902.53 2809.20; ...
    2517.07 1362.97 1185.58; ...
    66.22 1142.40 1839.52;];


q = [1.5331; 3.216; 4.9010; 1.6192; 2.835; 1.89; 1.77; 1.32; 3.42;...
    2.753; 1.65; 5.05; -3.16; -2.00; -3.712; -3.2504;];
    
%}


%{
% 2 dimensional tensegrity structure

CONN = [1 1 2;...
    2 2 3;...
    3 3 4;...
    4 1 3;...
    5 2 4;...
    6 1 4;...
    ];

coordinates = [0 1;...
    0 0;...
    1 0;...
    1 1];

q = [1;1;1;-1;1;-1];
n = 4;

%}
fitness = Nuwan_FDM(CONN, n, q, coordinates);