function [ new_am, new_ap, new_bm, new_bp ] = lower_right_correction(...
    am, ap, bm, bp, m, p, n)
%LOWER_RIGHT_CORRECTION Compute the lower right correction after product.
%
% [ new_am, new_ap, new_bm, new_bp ] = ...
% lower_right_correction(am, ap, bm, bp, m, p, n)
% computes the updated symbols necessary to build the factors of the
% correction in the product of two matrices as Hankel with coefficients
% given from the new symbols.

nam = length(am);
nap = length(ap);
nbm = length(bm);
nbp = length(bp);

am = reshape(am, 1, nam);
ap = reshape(ap, 1, nap);
bm = reshape(bm, 1, nbm);
bp = reshape(bp, 1, nbp);

% Compute lower-right correction
if m >= p   % left factor
    if 1+m-p <= nam
        new_am = am(1+m-p:end);
    else
        new_am = 0;
    end
    new_ap = [zeros(1,max(0,1+m-p-nam)), am(min(1+m-p,nam):-1:2), ap];
else
    if 1+p-m <= nap
        new_ap = ap(1+p-m:end);
    else
        new_ap = 0;
    end
    new_am = [zeros(1,max(0,1+p-m-nap)),ap(min(1+p-m,nap):-1:2),am];
end
if p >= n   % right factor
    if 1+p-n <= nbm
        new_bm = bm(1+p-n:end);
    else
        new_bm = 0;
    end
    new_bp = [zeros(1,max(0,1+p-n-nbm)),bm(min(1+p-n,nbm):-1:2),bp];
else
    if 1+n-p <= nbp
        new_bp = bp(1+n-p:end);
    else
        new_bp = 0;
    end
    new_bm = [zeros(1,max(0,1+n-p-nbp)),bp(min(1+n-p,nbp):-1:2),bm];
end


end

