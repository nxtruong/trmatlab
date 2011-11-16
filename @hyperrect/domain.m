function obj = domain(obj, A, F, Q)
%DOMAIN Compute affine reverse image (domain of an affine transformation)
%
%       R = domain(H, A, F, Q)
%
%   This function computes a hyper-rectangle R contained in a
%   hyper-rectangle Q that is mapped to the given hyper-rectangle H using
%   an affine map (A, F):
%       R = {x0 in Q | x1 = A*x0 + F in H}
%   where A and F are both vectors of the same dimension as H and Q, and
%   A*x0 + F = [A(i)*x0(i) + F(i)]_i.
%
%   All entries of A must be non-zero.
%
%   H and Q can be a hyper-rectangle or a union of them.
%
%   Note that this function is meant to be fast, so it does not check its
%   arguments for validity (e.g. types, dimensions, etc) except that:
%   + A and F can be either row or column vectors.
%   + Q can be empty ([]) or ignored, in that case Q = R^n
%
%   Otherwise, you can always use R = intersect((H-F)*(1./A), Q).
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

hasQ = nargin >= 4 && ~isempty(Q);

if obj.dims == 0, return; end

N = size(obj.L, 2);
F = F(:);
A = A(:);

% Transform the sets

obj.L = bsxfun(@rdivide, bsxfun(@minus, obj.L, F), A);
obj.H = bsxfun(@rdivide, bsxfun(@minus, obj.H, F), A);

if hasQ
    % Take intersection
    obj = intersect(obj, Q);
end

end
