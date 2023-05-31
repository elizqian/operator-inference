function H = F2H(F)
% converts F matrix operator which operates on get_x_sq output (no
% Kronecker redundancy) to tensorized matrix operator H which operates on
% kron(x,x) output (with redundancy)
%
% INPUT
% F     w-by-(w*(w+1)/2) matrix operator that acts on get_x_sq.m output,
%
% OUTPUT
% H     w-by-w^2 tensorized matrix operator that acts on kron(x,x)
%
% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 17 June 2019

[n,~] = size(F);

[ii,jj,vv] = find(F);
iH = [];
jH = [];
vH = [];

bb = 0;
for i = 1:n     % loop thru chunks of F 
    cc = bb+n+1-i;
    sel = jj>bb & jj<=cc;   % find indices corresponding to said chunk
    itemp = ii(sel);
    jtemp = jj(sel);
    vtemp = vv(sel);
    for j = 1:length(jtemp)
        sj = jtemp(j)-bb;
        if sj==1 % then it's a diagonal term and doesn't need to be split
            iH = [iH; itemp(j)];
            jH = [jH; (i-1)*n+i+sj-1];
            vH = [vH; vtemp(j)];
        else % it's a cross term and needs to be split
            iH = [iH; itemp(j); itemp(j)];
            jH = [jH; (i-1)*n+i+sj-1; (i+sj-2)*n+i];
            vH = [vH; vtemp(j)/2; vtemp(j)/2];
        end
    end
    bb = cc;
end
H = sparse(iH,jH,vH,n,n^2);