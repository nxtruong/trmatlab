function obj = mtimes(obj, A)
%MTIMES Multiplication between a hyper-rectangle and a scalar/vector/matrix
%   Three cases (H is a hyper-rectangle or a union of those):
%   1. H*c where c is a scalar: multiply each dimension of H with c.
%   2. H*v where v is a vector of the same dimension as H: multiply
%       dimension i of H with v(i).
%   3. H*A where A is a matrix of the same dimension as H: performs a
%       linear transformation of H; a polytope is returned.
%
%   Examples:
%       H1 = H*0.5;
%       H1 = H*[0.2; 0.3];
%       H1 = H*[0.3, 0; 0, 0.3];  % Same as above
%       P1 = H*[1 0; 1 1];  % A polytope is returned
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

if isa(obj, 'hyperrect') && isnumeric(A)
    vec = A;
elseif isa(A, 'hyperrect') && isnumeric(obj)
    vec = obj;
    obj = A;
else
    error('Unsupported multiplication of hyper-rectangles.');
end

% If obj is empty, return it immediately
if obj.dims == 0, return; end

if isscalar(vec)
    obj.L = obj.L * vec;
    obj.H = obj.H * vec;
    
elseif isvector(vec)
    % a vector of the same dimension
    assert(length(vec) == obj.dims,...
        'The vector must have the same dimension as the hyper-rectangle.');
    
    vec = vec(:);
    
    obj.L = bsxfun(@times, obj.L, vec);
    obj.H = bsxfun(@times, obj.H, vec);
    
else
    % A matrix
    obj = to_polytope(obj) * vec;
end

end

