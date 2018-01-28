function [CL, CH] = regiondiff(AL, AH, BL, BH)
% C is the set difference between union of hyperrectangles A and a single
% hyperrectangle B

ndims = size(AL, 1);

% Find sets in A that do not intersect B
idx = ~(all(bsxfun(@lt, BL, AH) & bsxfun(@gt, BH, AL)));
CL = AL(:, idx);
CH = AH(:, idx);
AL(:, idx) = [];
AH(:, idx) = [];

nA = size(AL, 2);
if nA == 0, return; end

%% Vectorization
% The following code vectorizes the set diff operation, but it is slower
% than the simple loop code.

% For each dimension i, if AL(i, k) < BL(i) then 1 hyperrectangle is added;
% one more added if AH(i, k) > BH(i); the number can add up (to 2), or
% nothing at all (0).  So we can quickly compute the number of new sets.
% split(i, j) = 0/1/2 if dimension i of AL(:, j) will (not be split)/(be
% resized)/(be split to 2)
% split = bsxfun(@lt, AL, BL) + bsxfun(@gt, AH, BH);
% 
% 
% % Work in each dimension
% for k = 1:ndims
%     % For splitting to 2
%     idx = split(k, :) == 2;
%     if any(idx)
%         tmp = [bsxfun(@max, AL(1:k-1, idx), BL(1:k-1)); AL(k:end, idx)];
%         CL = [CL, tmp];
%         tmp(k, :) = BH(k);
%         CL = [CL, tmp];
%         tmp1 = [bsxfun(@min, AH(1:k-1, idx), BH(1:k-1)); AH(k:end, idx)];
%         tmp = tmp1;
%         tmp(k, :) = BL(k);
%         CH = [CH, tmp, tmp1];
%     end
%     
%     % For splitting to 1
%     idx = split(k, :) == 1;
%     if any(idx)
%         tmp = AL(k, idx);
%         tmp1 = BH(k) < AH(k, idx);
%         CL = [CL,...
%             [bsxfun(@max, AL(1:k-1, idx), BL(1:k-1));
%              tmp + (BH(k) - tmp).*tmp1;
%              AL(k+1:end, idx)]];
%         
%         CH = [CH,...
%             [bsxfun(@min, AH(1:k-1, idx), BH(1:k-1));
%              BL(k) + (AH(k, idx) - BL(k)).*tmp1;
%              AH(k+1:end, idx)]];
%     end
% end

% For each set in A, compute the difference between it and B
for k = 1:nA
    curL = AL(:, k);
    curH = AH(:, k);
    
    % Take the difference between (curL, curH) and
    % (BL,BH) by considering each dimension, create subsets
    for d = 1:ndims
        if curL(d) < BL(d)
            % split one rectangle
            CL = [CL, curL];
            CH = [CH, curH];
            CH(d, end) = BL(d);
            curL(d) = BL(d);
        end
        
        if curH(d) > BH(d)
            % split one rectangle
            CL = [CL, curL];
            CL(d, end) = BH(d);
            CH = [CH, curH];
            curH(d) = BH(d);
        end
    end
end

end
