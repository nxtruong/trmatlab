function p = analyticCenter(obj)
%ANALYTICCENTER  Computes the center point of the hyperrectangle.
%
%   x = analyticCenter(H)
%
%   Returns the center points of the hyperrectangles H in x.
%   If H is a union of multiple hyperrectangles, corresponding center
%   points are returned in the columns of x.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

assert(obj.dims > 0, 'The hyperrectangle must be not empty.');

p = (obj.L + obj.H)/2;

end
