function iss = issubset(A, B, randpntcheck)
%ISSUBSET Check if a hyper-rectangle is a subset of another.
%   iss = issubset(A, B)
%     returns true if A is a subset of B (A \subseteq B).
%
%   Both A and B must be of class hyperrect.
%   Both can be a single set or an array of sets (a union of sets).
%
%   An empty set A is a subset of every non-empty set B.
%   Any set A is not a subset of an empty set B.
%
%   This function does not check the bounding box, so either make sure that
%   bounding box of A is a subset of the bounding box of B, or use A <= B.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

assert(isa(A, 'hyperrect') && isa(B, 'hyperrect'),...
    'Both operands must be of class hyperrect.');

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
ndims = A.dims;

if nargin < 3, randpntcheck = false; end

% Sort sets in A and B from large to small, to increase chance of breaking
% the loops sooner
[~, idxA] = sort(prod(A.H - A.L), 'descend');

if randpntcheck
    nrandpnts = min(size(B.L, 2)*ndims*10, round(10000/ndims));
    [~, idxB] = sort(prod(B.H - B.L), 'descend');
else
    [~, idxB] = sort(prod(B.H - B.L), 'descend');
end

B.H = B.H(:, idxB);
B.L = B.L(:, idxB);

iss = true;

for kA = idxA
    % For each A(k), check if it is a subset of B
    
    if randpntcheck
        % First, generate a number of random points in A(k) and check if they
        % all belong to B
        randpnts = rand(ndims, nrandpnts);
        for x = randpnts
            p = A.L(:, kA) .* x + A.H(:, kA) .* (1-x);
            if ~any(all(bsxfun(@le, B.L, p)) & all(bsxfun(@le, p, B.H)))
                iss = false;
                break;
            end
        end
        
        if ~iss
            % disp('Found by random points.');
            break;
        end
    end
    
    % Random point check passed, now we compute A(k) \ B(i) and check if
    % the final result is empty
    CL = A.L(:, kA);
    CH = A.H(:, kA);
    
    % Find sets in B that intersect (CL, CH)
    idxB = find(all(bsxfun(@lt, CL, B.H) & bsxfun(@gt, CH, B.L)));
    
    for kB = idxB
        [CL, CH] = regiondiff(CL, CH, B.L(:,kB), B.H(:, kB));
            
        % Remove thin sets from (CL, CH)
        if ~isempty(CL)
            thinidx = any(CL >= CH);
            CL(:, thinidx) = [];
            CH(:, thinidx) = [];
        end
            
        % Check if (CL, CH) is empty then we can stop
        if isempty(CL)
            break;
        end
    end
    
    % If (CL, CH) is empty then A(kA) is a subset of B, we continue;
    % otherwise, we can stop because A is not a subset of B
    if ~isempty(CL)
        iss = false;
        break;
    end
end

end

