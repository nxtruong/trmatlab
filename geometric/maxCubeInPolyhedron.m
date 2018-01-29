function [L, R] = maxCubeInPolyhedron(A, B, solver, options)
%MAXCUBEINPOLYHEDRON Find the largest hypercube inside a polyhedron.
%   [L, R] = maxCubeInPolyhedron(A, B)
% Given a polyhedron P = {x | A*x <= B}, find the largest hypercube C =
% {x | L <= x <= L + R} that fits inside P. By largest, we are looking for
% C with the maximal R value.
%
% This function requires tropt.solveLP function, with a supported linear
% solver. A specific solver can be selected by an additional argument:
%   [L, R] = maxCubeInPolyhedron(A, B, SOLVER, OPTIONS)
% where the optional argument OPTIONS is the solver's option structure.
%
% If R is -inf, the problem is infeasible (P is empty for example); if R
% is +inf, the problem is unbounded (P is unbounded). If there is any other
% error, R is empty.
%
% See: tropt.solveLP
%
% (C) 2018 by Truong X. Nghiem.

% History:
% 2018-01-29 Truong started the function.

% Algorithm:
% It can be shown that C is inside P iff all vertices of C are inside P iff
% Ap*H + An*L <= B, where Ap = the non-negative elements of A = max(0, A),
% An = the non-positive elements of A = min(0, A), H = L + R. Note that A =
% Ap + An. So we solve the following LP:
%       maximize R
%       s.t. Ap*(L + 1*R) + An*L = A*L + Ap*1*R <= B, R >= 0

% Check arguments
narginchk(2, 4);
assert(isnumeric(A) && isnumeric(B), 'A and B must be numerical values.');
assert(isvector(B), 'B must be a vector.');
B = B(:);
[m, n] = size(A);
assert(length(B) == m, 'Dimensions of A and B mismatch.');

% Check solver and solver's options
if nargin > 2
    assert(ischar(solver), 'SOLVER must be a string.');
else
    solver = '';
end
if nargin > 3
    options = {options};
else
    options = {};
end

% Construct the LP and solve it
% The decision variable is x = [R; L], so
% min -[1, 0, ... , 0]*x
% s.t. [0; -inf; ...; -inf] <= x <= inf
%      [Ap*1, A] * x <= B

ApOne = sum(max(0, A), 2);
[xopt, ~, exitflag] = tropt.solveLP(solver, ...
    [-1; zeros(n,1)], [ApOne, A], B, [], [], ...
    [0; -inf(n,1)], [], options{:});

if exitflag == 0
    % Successful
    R = xopt(1);
    L = vec(xopt(2:end));
elseif exitflag == 1
    % Infeasible
    L = [];
    R = -inf;
elseif exitflag == 2
    % Unbounded
    L = [];
    R = inf;
else
    R = [];
    L = [];
end
end

