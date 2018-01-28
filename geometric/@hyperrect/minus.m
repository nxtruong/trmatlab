function obj = minus(obj1, obj2)
%MINUS Translate a hyperrectangle by a vector, or Minkowski difference.
%
%   M = H - v       (where v is a vector)
%       M is H translated by vector v, i.e. M = {x-v: x in H}
%
%   M = H - a       (where a is a scalar)
%       M is H translated by -a in each dimension.
%
%   M = H1 - H2     (where H1 and H2 are both single hyperrectangle)
%       Minkowski difference, the result is a hyperrectangle.
%
%   M = H - V       (where H: single hyperrect)
%       where V is a matrix of points in its columns, which represent the
%       vertices of a polytope. M is the Minkowski difference between the
%       two. M is a hyperrectangle.
%
%   M = H - P       (where H: single hyperrect, P: single polytope)
%       Minkowski difference; M is a hyperrectangle.
%
%   M = H - P
%       where H and P are as above, however either or both of them are
%       unions (arrays) of objects; then the objects are converted to
%       polytopes and MPT/polytope is used to compute the Minkowski
%       difference; M is a polytope (or polytope array).
%
% See also: hyperrect/plus, polytope/minus
%
% (C) 2011-2012 by Truong X. Nghiem (nghiem@seas.upenn.edu)

% History:
%   2012-05-10  Improved and added Minkowski difference.

% obj1 is always a hyperrect

if ~isfulldim(obj1)
    % Empty hyperrect
    obj = obj1;
    return;
end

if isa(obj1,'hyperrect') && ...
        isnumeric(obj2) && (isscalar(obj2) || isvector(obj2))
    obj = plus(obj1, -obj2);
elseif length(obj1) > 1 ||...
        ((isa(obj2,'hyperrect') || isa(obj2,'polytope')) && length(obj2) > 1)
    % If either obj1 or obj2 is a union then we need to use polytope
    % computation
    if isa(obj2, 'hyperrect')
        obj2 = to_polytope(obj2);
    end
    obj = minus(to_polytope(obj1), obj2);
else
    % So obj1 is single and obj2 is single
    if isa(obj2, 'polytope')
        BB = bounding_box(obj2, struct('Voutput', 1));
    elseif isa(obj2, 'hyperrect')
        BB = [obj2.L, obj2.H];
    elseif isnumeric(obj2)
        % obj2 is a matrix whose columns are the vertices of a polytope
        BB = [min(obj2, [], 2), max(obj2, [], 2)];
    else
        error('Unsupported class of the second object: %s', class(obj2));
    end
    
    % Check the dimension of obj1 and BB (obj2)
    assert(dimension(obj1) == size(BB, 1),...
        'Both objects must have the same dimension.');
    
    % New thresholds for the resulted hyperrect
    L = obj1.L - BB(:, 1);
    H = obj1.H - BB(:, 2);
    
    if any(L > H)
        % It is empty
        obj = hyperrect;
    else
        obj = hyperrect(L, H);
    end
end

end
