function H = vertcat(varargin)
%VERTCAT Concatenates hyper-rectangles into an array (union)
%
%   H = vertcat(varargin)
%
%   Concatenates input hyper-rectangles and creates an array, which is the
%   union of those input hyper-rectangles.
%
%   Usage:
%       H = [H1; H2; H3];
%       AllH = [H; H4; H5];
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

H = horzcat(varargin{:});

end

