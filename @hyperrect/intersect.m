function H1 = intersect(H1, H2)
%INTERSECT Intersection of two hyper-rectangles.
%
%   H = intersect(H1, H2)
%
%   Returns the intersection of two hyper-rectangles H1 and H2.  Either or
%   both of them can be hyper-rectangle arrays (considered as unions of
%   hyper-rectangles).
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

if H1.dims == 0 || H2.dims == 0
    H1 = hyperrect;
    return;
end

assert(H1.dims == H2.dims, 'Both hyper-rectangles must have the same dimension.');

% If the total number of elements is too large, we use the bounding boxes
% to remove unnecessary sets before actually computing the intersection
if H1.dims * size(H1.L, 2) * size(H2.L, 2) > 10000
    % Compute the bounding boxes of H1 and H2, then their intersection
    [Bl, Bh] = intersectBoxes(min(H1.L, [], 2), max(H1.H, [], 2),...
        min(H2.L, [], 2), max(H2.H, [], 2));

    % If these boxes are disjoint, H1 and H2 are disjoint
    if emptyBox
        H1 = hyperrect;
        return;
    end

    % In H1 and H2, remove those hyper-rectangles that do not intersect the box
    % B (the intersection must be a subset of this box).
    idx = any(bsxfun(@le, H1.H, Bl) | bsxfun(@ge, H1.L, Bh));
    H1.L(:, idx) = [];
    H1.H(:, idx) = [];
    
    idx = any(bsxfun(@le, H2.H, Bl) | bsxfun(@ge, H2.L, Bh));
    H2.L(:, idx) = [];
    H2.H(:, idx) = [];
end


% Compute the intersection
[H1.L, H1.H] = intersectBoxes(H1.L, H1.H, H2.L, H2.H);

if isempty(H1.L)
    H1 = hyperrect;
end

end