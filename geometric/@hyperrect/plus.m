function obj = plus(obj1, obj2)
%PLUS Translate a hyperrectangle by a vector or a scalar.
%
%   M = H + v       (where v is a vector)
%       M is H translated by vector v, i.e. M = {x+v: x in H}
%
%   M = H + a       (where a is a scalar)
%       M is H translated by a in each dimension.
%
%   H can be an array of hyperrectangles.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

if isa(obj2,'hyperrect') && isnumeric(obj1)
    obj = obj2;
    vec = obj1;
elseif isa(obj1,'hyperrect') && isnumeric(obj2)
    obj = obj1;
    vec = obj2;
else
    error('Unsupported addition operation.');
end

% obj is a hyperrectangle, vec is the distance


if obj.dims == 0 || isempty(vec), return; end

% assert(isscalar(vec) || length(vec) == obj.dims,...
%     'The hyper-rectangle and the vector are not of the same dimension.');

if isscalar(vec)
    obj.L = obj.L + vec;
    obj.H = obj.H + vec;
else
    vec = vec(:);
    obj.L = bsxfun(@plus, obj.L, vec);
    obj.H = bsxfun(@plus, obj.H, vec);
end
    
end
