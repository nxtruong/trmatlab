function isin = isinside(obj, x0)
%ISINSIDE Check if a point is inside the given hyperrectangle(s)
%   isin = = isinside(H, x0)
%
%   H can be an array of hyperrectangles, which is the union of the
%   individual hyperrectangles.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

if obj.dims == 0
    isin = false;
    return;
end

assert(length(x0) == obj.dims, 'Incompatible dimensions.');

x0 = x0(:);

% Check if x0 is in one of the hyperrectangles
isin = any(all(bsxfun(@le, obj.L, x0) & bsxfun(@le, x0, obj.H)));

end

