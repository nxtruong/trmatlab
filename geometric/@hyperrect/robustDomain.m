function obj = robustDomain(obj, A1, F1, A2, F2, Q)
%ROBUSTDOMAIN Compute robust affine reverse image (domain).
%
%       R = robustDomain(H, A1, F1, A2, F2, Q)
%
%   This function computes a hyper-rectangle R contained in a
%   hyper-rectangle Q that is robustly mapped to the given hyper-rectangle
%   H using affine maps (A1, F1) and (A2, F2):
%       R = {x0 in Q | (A1*x0 + F1) in H and (A2*x0 + F2) in H}
%   where A1, F1, A2, F2 are vectors of the same dimension as H and Q, and
%   Ak*x0 + Fk = [Ak(i)*x0(i) + Fk(i)]_i.
%
%   All entries of A1, A2 must be non-zero.
%
%   H and Q can be a hyper-rectangle or a union of them.
%
%   Note that this function is meant to be fast, so it does not check its
%   arguments for validity (e.g. types, dimensions, etc) except that:
%   + Ak and Fk can be either row or column vectors.
%   + Q can be empty ([]) or ignored, in that case Q = R^n
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

hasQ = nargin >= 4 && ~isempty(Q);

if obj.dims == 0, return; end

N = size(obj.L, 2);

% Transform the sets
A1 = A1(:);
F1 = F1(:);

L1 = bsxfun(@rdivide, bsxfun(@minus, obj.L, F1), A1);
H1 = bsxfun(@rdivide, bsxfun(@minus, obj.H, F1), A1);

A2 = A2(:);
F2 = F2(:);

L2 = bsxfun(@rdivide, bsxfun(@minus, obj.L, F2), A2);
H2 = bsxfun(@rdivide, bsxfun(@minus, obj.H, F2), A2);


% Compute the intersections of the corresponding subsets (hyperrectangles),
% that is intersect(domain(obj(k), A1, F1), domain(obj(k), A2, F2))
obj.L = max(L1, L2);
obj.H = min(H1, H2);

% Remove empty and thin subsets
empty = any(obj.L >= obj.H);
obj.L(:,empty) = [];
obj.H(:,empty) = [];

if isempty(obj.L), obj = hyperrect; return; end

if hasQ
    % Take intersection with Q
    obj = intersect(obj, Q);
end

end
