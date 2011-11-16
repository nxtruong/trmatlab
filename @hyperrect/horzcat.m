function obj = horzcat(varargin)
%HORZCAT Concatenates hyper-rectangles into an array (union)
%
%   H = horzcat(varargin)
%
%   Concatenates input hyper-rectangles and creates an array, which is the
%   union of those input hyper-rectangles.
%
%   Usage:
%       H = [H1 H2 H3];
%       AllH = [H H4 H5];
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

assert(all(cellfun('isclass', varargin, 'hyperrect')),...
    'Cannot concatenate hyperrect objects with non-hyperrect objects.');

% Obtain the dimensions
allDims = cellfun(@(x) x.dims, varargin);

% Remove 0-dimension (empty set)
emptyIdx = (allDims == 0);
varargin(emptyIdx) = [];

if isempty(varargin)
    % All sets are empty, so the result is empty
    obj = hyperrect;
    return;
end

assert(~any(diff(allDims(~emptyIdx))),...
    'All hyper-rectangles to be concatenated must have the same dimension.');

obj = varargin{1};
n = length(varargin);

for k = 2:n
    obj.L = [obj.L, varargin{k}.L];
    obj.H = [obj.H, varargin{k}.H];
end

end

