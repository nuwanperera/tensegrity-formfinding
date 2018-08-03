function [ fitness ] = tensegrity_FDM_fitness(n, CONN, q, coordinates )
%TENSEGRITY_FDM_FITNESS This function determines whether or not a
%tensegrity structure is self-stressed, thus at self-stressed equilibrium.
%   We use to force density method (Schek (1974), Fund (2008)) to determine
%   whether or not a tensegrity structure is at self-stressed equilibrium
%
% Input variables: 
% C - Connectivity matrix
% n - number of nodes
% q - force-density vector
% coordinates - x,y,z nodal coordinates of tensegrity structure
%
%
% Output variable:
% fitness - determines the fitness of the structure
% The output should minimize to zero for a self-stressed structure
[m,~] = size(CONN);
C = zeros(m,n);
% Populate Connectivity Matrix C
for j = 1:m
    r = CONN(j,2);
    C(j,r) = 1;
    s = CONN(j,3);
    C(j,s) = -1;
end

x = coordinates(:,1);
y = coordinates(:,2);
z = coordinates(:,3);
 
% Compute self-equilibrium matrix
A = [C'*diag(C*x); C'*diag(C*y); C'*diag(C*z);];

% Solve for Aq = 0 using SVD
[~,~,V] = svd(A);
new = V(:, end);
 
fitness = sum(abs(A*new));
 
end
