function obj = bounding_box(obj)
%BOUNDING_BOX Compute a bounding box for a given (union of) hyperrectangles
%   B = bounding_box(H)
%       Returns the bounding box for a given union of hyperrectangles.  If
%       H is a single hyper-rectangle, it is its bounding box.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

% Compute the max and min of the ends of the intervals in each
% dimension
obj = hyperrect(min(obj.L, [], 2), max(obj.H, [], 2));

end

