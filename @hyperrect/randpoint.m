function p = randpoint(obj, n)
%RANDPOINT Returns a random point in the given hyperrectangle.
%   x = randpoint(H) returns one random point in H
%   X = randpoint(H, N) returns N random points in H as columns of X
%
%   H must be a single hyperrectangle.
%   If H is truly empty (i.e. dimension = 0) then x is empty.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

assert(length(obj) == 1, 'A single hyperrectangle polytope must be provided.');

if nargin < 2
    n = 1;
else
    assert(n >= 1, 'At least 1 random point must be generated (N >= 1).');
end

if obj.dims > 0
    x = rand(obj.dims, n);
    p = bsxfun(@times, obj.L, x) + bsxfun(@times, obj.H, 1-x);
else
    p = [];
end

end
