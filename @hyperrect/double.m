function [l,h] = double(obj)
%DOUBLE  Returns the intervals of a hyper-rectangle.
%
%   [l,h] = double(H)
%       l is the vector of the lower ends of the intervals, while h is the
%       vector of the higher ends; so H = [l(1),h(1)] x ... x [l(n),h(n)]
%
%   If H contains multiple hyperrectangles then l and h are matrices where
%   each column corresponds to one hyperrectangle.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

if obj.dims == 0
    l = [];
    h = [];
else
    l = obj.L;
    h = obj.H;
end

end
