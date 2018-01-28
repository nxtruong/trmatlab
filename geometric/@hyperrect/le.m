function iss = le(A, B)
%LE Check if a hyperrectangle set is a subset of another.
%
%   H1 <= H2
%       Returns true if H1 is a subset of H2; both of them can be an array
%       (union of hyper-rectangles).
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

if isa(A, 'polytope') && isa(B, 'hyperrect')
    iss = A <= to_polytope(B);
    return;
elseif isa(B, 'polytope') && isa(A, 'hyperrect')
    iss = to_polytope(A) <= B;
    return;
elseif isa(A, 'hyperrect') && isa(B, 'hyperrect')
    % Check dimensions
    if B.dims == 0
        % if B is empty, no other set is its subset
        iss = false;
        return;
    end
    
    if A.dims == 0
        % if A is empty, it's a subset of every non-empty set.
        iss = true;
        return;
    end
    
    assert(A.dims == B.dims, 'Dimensions of both sets must agree.');
    
    % Check the bounding box first
    bbA = bounding_box(A);
    bbB = bounding_box(B);
    iss = all(bbA.L >= bbB.L) && all(bbA.H <= bbB.H);
    if iss
        iss = issubset(A, B);
    end
else
    error('Unsupported classes of operands.');
end

end
