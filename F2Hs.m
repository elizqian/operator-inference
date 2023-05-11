function H = F2Hs(F)
    n = size(F,1);
    H = F * elimat(n) * nommat(n,n);
end