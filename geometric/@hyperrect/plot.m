function plot(obj,varargin)
%PLOT  Plot 2-D or 3-D hyperrectangles.
%   See the function polytope/plot
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

assert(2 == obj.dims || 3 == obj.dims, 'Only 2-D or 3-D hyperrectangles can be plotted.');

plot(to_polytope(obj), varargin{:});
end
