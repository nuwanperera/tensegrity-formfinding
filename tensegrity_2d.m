function [ output_args ] = tensegrity_2d( n, CONN, q, coordinates )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[m,~] = size(CONN);
C = zeros(m,n);
% Populate Connectivity Matrix
for j = 1:m
    r = CONN(j,2);
    C(j,r) = 1;
    s = CONN(j,3);
    C(j,s) = -1;
end

A = [C'*diag(C*x); C'*diag(C*y)];


end

