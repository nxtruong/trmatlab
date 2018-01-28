function iss = rawsubset(AL, AH, BL, BH, randpntcheck, sortA, sortB)
%RAWSUBSET Check if a hyper-rectangle is a subset of another.
%   iss = rawsubset(AL, AH, BL, BH, randpntcheck, sortA, sortB)
%     returns true if A is a subset of B (A \subseteq B).
%
%   If sort* = true, * will be sorted in descending order of volumes; if
%   sort* = false, * should be already sorted.
%
%   This is an internal function and is not supposed to be used outside the
%   class hyperrect.
%
%   An empty set A is a subset of every non-empty set B.
%   Any set A is not a subset of an empty set B.
%
%   This function does not check the bounding box, so either make sure that
%   bounding box of A is a subset of the bounding box of B, or use A <= B.
%
% (C) 2012 by Truong X. Nghiem (nghiem@seas.upenn.edu)

ndims = size(AL, 1);

% Sort sets in A and B from large to small, to increase chance of breaking
% the loops sooner
if sortA
    [~, idxA] = sort(prod(AH - AL), 'descend');
else
    idxA = 1:size(AL, 2);
end

if randpntcheck
    nrandpnts = min(size(BL, 2)*ndims*10, round(10000/ndims));
end

if sortB
    [~, idxB] = sort(prod(BH - BL), 'descend');
    BH = BH(:, idxB);
    BL = BL(:, idxB);
end

iss = true;

for kA = idxA
    % For each A(k), check if it is a subset of B
    
    if randpntcheck
        % First, generate a number of random points in A(k) and check if they
        % all belong to B
        randpnts = rand(ndims, nrandpnts);
        for x = randpnts
            p = AL(:, kA) .* x + AH(:, kA) .* (1-x);
            if ~any(all(bsxfun(@le, BL, p)) & all(bsxfun(@le, p, BH)))
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
    CL = AL(:, kA);
    CH = AH(:, kA);
    
    % Find sets in B that intersect (CL, CH)
    idxB = find(all(bsxfun(@lt, CL, BH) & bsxfun(@gt, CH, BL)));
    
    for kB = idxB
        [CL, CH] = regiondiff(CL, CH, BL(:,kB), BH(:, kB));
            
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

