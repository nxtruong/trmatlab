function iss = issubset(A, B, randpntcheck)
%ISSUBSET Check if a hyper-rectangle is a subset of another.
%   iss = issubset(A, B)
%     returns true if A is a subset of B (A \subseteq B).
%
%   Both A and B must be of class hyperrect.
%   Both can be a single set or an array of sets (a union of sets).
%
%   An empty set A is a subset of every non-empty set B.
%   Any set A is not a subset of an empty set B.
%
%   This function does not check the bounding box, so either make sure that
%   bounding box of A is a subset of the bounding box of B, or use A <= B.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

assert(isa(A, 'hyperrect') && isa(B, 'hyperrect'),...
    'Both operands must be of class hyperrect.');

if B.dims == 0
    % if B is empty, no other set is its subset
    iss = false;
    return;
end

if A.dims == 0
    % if A is empty, it's a subset of every non-empty set.
    iss = true;
    return;
end

assert(A.dims == B.dims, 'Dimensions of both sets must agree.');

if nargin < 3, randpntcheck = false; end

iss = rawsubset(A.L, A.H, B.L, B.H, randpntcheck, true, true);

end

