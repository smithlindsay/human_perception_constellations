function Ahat = human_expectations(A, beta)
% Inputs: NxN transition matrix A (row stochastic), and inverse temperature
% parameter (or acurracy of human learning) beta
%
% Output: NxN transition matrix Ahat (row stochastic) that represents human
% expectations of the transition matrix A.

N = size(A,1);

Ahat = (1-exp(-beta))*A/(eye(N) - exp(-beta)*A);