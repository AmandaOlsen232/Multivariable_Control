!Problem 1
clear all 
close all 
clc


A = [-3  5;
      7 -10];
B = [-1;
      3];

A*B

%%
!Problem 2 
clear all 
close all 
clc

A = [4 5 1; 
     3 7 10;
     1 0 1];
B = [1; 
     2;
     3];

A*B

%%
!Problem 3
clear all
close all
clc

dot_prod = dot([1,1], [1, -1])

x = [1;
     1]'
z = [1;
    -1]

x*z

%%
!Problem 4
close all
clear all
clc

x = [1;
     1]';
x_orth = [2;
          -2];

y = [-3;
     2;
     4]';
y_orth = [1;
          1;
          1/4];

z = [1;
    -1]';
z_orth =[2;
         2];

x*x_orth
y*y_orth
z*z_orth

%%
!Problem 5
clear all
close all
clc

A = [1 4 2;
     0 0 0;
     1 0 9];
B = [1 4;
     0 0;
     1 0];

rank(A)
rank(B)

%%
!Problem 6
close all 
clear all 
clc

% A1 = [2 3 5;
%       -4 2 3];
% 
% n_A1 = [-1/16;
%         -13/8;
%         1];
% test1 = A1*n_A1

% A2 = [1 0 1;
%       5 2 1;
%       1 2 2];
% test2 = null(A2)
% 
A3 = [2 1 1;
      1 1 0;
      1 0 1];
n_A3 = [-1; 1; 1];
test3 = A3*n_A3
% null(A3);
%%
!Problem 7
clear all
close all
clc

% A1 = [2 3 5;
%       -4 2 3];
% rank(A1)
% rref(A1)


% A2 = [1 0 1;
%       5 2 1;
%       1 2 2];
% rank(A2)

A3 = [2 1 1;
      1 1 0;
      1 0 1];

rank(A3)

rref(A3)

%%
!Problem 8
close all
clear all
clc

% Im_A2 = [1 0 1;
%          5 2 1;
%          1 2 2];
% null_A2 = [0;0;0];
% 
% Im_A2*null_A2;

Im_A3 = [2;
         1;
         1]
null_A3 = [-1;1;1]
dot(Im_A3, null_A3)

%%
!Problem 9
A2 = [1 0 1;
      5 2 1;
      1 2 2];
[V_2,D_2] = eig(A2);

A3 = [2 1 1;
      1 1 0;
      1 0 1];
[V_3,D_3] = eig(A3);

A4 = [4 5 1;
      3 7 10;
      1 0 1];
[V_4,D_4] = eig(A4);

A5 = [4 0 0;
      0 7 0;
      0 0 1];
[V_5,D_5] = eig(A5);

A6 = [2 1 1;
      0 5 10;
      0 0 8];
[V_6,D_6] = eig(A6)
%%
!Problem 10
close all
clear all 
clc

A2 = [1 0 1;
      5 2 1;
      1 2 2];

A3 = [2 1 1;
      1 1 0;
      1 0 1];

A4 = [4 5 1;
      3 7 10;
      1 0 1];

A5 = [4 0 0;
      0 7 0;
      0 0 1];

A6 = [2 1 1;
      0 5 10;
      0 0 8];

rank2 = rank(A2)
rank3 = rank(A3)
rank4 = rank(A4)
rank5 = rank(A5)
rank6 = rank(A6)

%%
!coding
clear all
close all
clc

A = [2 0 0; 2 2 2; 3 0 -1];
[V, D] = eig(A)