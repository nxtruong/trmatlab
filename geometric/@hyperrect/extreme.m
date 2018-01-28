function V = extreme(obj)
%EXTREME  Returns the extreme points (vertices).
%   This function returns the extreme points (vertices) of the given
%   single hyper-rectangle.
%       V = extreme(H)
%
%   H must be a single hyper-rectangle (can be unbounded).
%   The rows of V contain the extreme points.
%
%   NOTE: I think it would be better to store the extreme points in the
%   columns of V, but I followed the convention used in the polytope class
%   in the MPT toolbox.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

assert(length(obj) == 1, 'Only a single hyperrectangle can be provided.');

if obj.dims == 0
    V = [];
    return;
end

if obj.L(1) == obj.H(1)
    V = obj.L(1);
    n = 1;
else
    V = [obj.L(1); obj.H(1)];
    n = 2;
end

for k = 2:obj.dims
    if obj.L(k) == obj.H(k)
        % A plane in this dimension
        V(:,k) = obj.L(k);    % Add one column to the right c.t. dimension k
    else
        % An interval in this dimension
        V = repmat(V, 2, 1);
        V(1:n, k) = obj.L(k);
        V(n+1:end, k) = obj.H(k);
        n = 2*n;
    end
end

end