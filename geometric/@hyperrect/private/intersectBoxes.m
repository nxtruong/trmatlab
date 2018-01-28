function [L, H] = intersectBoxes(L1, H1, L2, H2)
%INTERSECTBOXES Intersection of two boxes, each can be a union of boxes.
%   [L, H] = intersectBoxes(L1, H2, L2, H2)
%   	returns the intersection (L, H) of two boxes (L1, H1) and (L2, H2)
%   	of the same dimension.  L1, H1, L2, H2 may contain multiple columns
%   	which represent a union of boxes.  The intersection L, H may be a
%   	union of boxes too, represented by columns.  All empty and thin
%   	sets will be removed.  If L is empty then the intersection is
%   	empty.
%
%   This function is a helper function and is supposed to be fast.  All
%   arguments must be provided and must be correct.  No checking is
%   performed.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

% The algorithm is to use bsxfun on 3-D arrays: (L1, H1) uses dimensions 1
% and 2, (L2, H2) uses dimensions 1 and 3.  The result will be unrolled to
% 2-D arrays of N1*N2 columns, where N1 and N2 are the numbers of columns
% of (L1, H1) and (L2, H2) respectively.

ndims = size(L1, 1);

L = reshape(bsxfun(@max, L1, reshape(L2, ndims, 1, [])), ndims, []);
H = reshape(bsxfun(@min, H1, reshape(H2, ndims, 1, [])), ndims, []);

emptyIdx = any(L >= H);
L(:, emptyIdx) = [];
H(:, emptyIdx) = [];

end

