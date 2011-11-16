function obj = minus(obj1, obj2)
%MINUS Translate reversly a hyperrectangle by a vector or a scalar.
%
%   M = H - v       (where v is a vector)
%       M is H translated by vector v, i.e. M = {x-v: x in H}
%
%   M = H - a       (where a is a scalar)
%       M is H translated by -a in each dimension.
%
%   H can be an array of hyperrectangles.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

if isa(obj2,'hyperrect') && isa(obj1,'numeric')
    obj = obj2;
    vec = obj1;
elseif isa(obj1,'hyperrect') && isa(obj2,'numeric')
    obj = obj1;
    vec = obj2;
else
    error('The addition between two hyper-rectangles is not supported.');
end

obj = plus(obj, -vec);

end
