function obj = hyperrect(A, B)
%HYPERRECT Create a hyperrect object (axis-aligned hyper-rectangle)
% Constructor for a hyperrect class object.
%
% Examples
%   hyperrect([0.5 1.5; 1 2])
%       the rectangle [0.5 1.5] x [1 2] in 2D
%   hyperrect([0.5; 1], [1.5; 2])
%       same as above but providing the lower and upper-end values
%       separately in two vectors
%   hyperrect()
%       an empty hyper-rectangle
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

% Last update: 2012-08-05

error(nargchk(0, 2, nargin));

if nargin == 0
    % An empty hyperrect, can be used when objects are loaded from disk
    obj = init_fields;
    obj = class(obj, 'hyperrect');
    return;
end

obj = init_fields;

if nargin == 1
    if isa(A, 'hyperrect')  % used when objects are passed as arguments
        obj = A;
        return;
    end
    assert(size(A, 2) == 2, 'First argument must be 2-column matrix of thresholds.');
    assert(all(A(:,1) <= A(:,2)), 'Lower thresholds cannot be larger than upper thresholds.');
    
    obj.dims = size(A, 1);
    obj.L = A(:,1);
    obj.H = A(:,2);
else
    B = B(:);
    A = A(:);
    
    % Dimensions
    assert(length(B) == length(A), 'Mismatched dimensions.');
    
    if all(A <= B)
        % Create the hyperrect
        obj.L = A;
        obj.H = B;
        obj.dims = length(A);
    end
    % If lower thresholds are larger than upper thresholds, create
    % empty hyperrect
end

% class name tag
obj = class(obj, 'hyperrect');


end

function obj = init_fields()
% Initialize all fields to dummy values

% The columns of L and H hold the lower- and upper-thresholds of the atomic
% hyperrectangles.
obj.L = [];
obj.H = [];

obj.dims = 0;       % Number of dimensions
end