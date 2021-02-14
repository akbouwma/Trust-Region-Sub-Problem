%{
% Solve generalized eigenvalue problem.
% Authors: Alan Bouwman, Caleb Jacobs
% Date last modified: 02/14/2021
%}

% Clean workspace
clear;
close all;

% Load given matrices
A1 = load('Test_Matrices/A1.mat').Expression1;
A2 = load('Test_Matrices/A2.mat').Expression1;
B1 = sparse(load('Test_Matrices/B1.mat').Expression1);
B2 = sparse(load('Test_Matrices/B2.mat').Expression1);

% Solve generalized eigenvalue problems
[VA, DA] = eigs(A1, -A2, length(A1)); % Eigen system of (A1+lambda*A2)x = 0
[VB, DB] = eigs(B1, -B2, length(B1)); % Eigen system of (B1+lambda*B2)x = 0

% Create vectors of eigenvalues
lamA = diag(DA);
lamB = diag(DB);

% Plots of eigenvalues
figure(1)
scatter(real(lamA), imag(lamA), '*')
title('Eigenvalues for A1 \cdot x = \lambda A2 \cdot x')

figure(2)
scatter(real(lamB), imag(lamB), '.')
title('Eigenvalues for B1 \cdot x = \lambda B2 \cdot x')
