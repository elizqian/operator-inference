function [A_interp, B_interp, F_interp] = spline_operators(ops, n, mus, mu_new)
    foo = mus > mu_new;
    idx = find(foo == 0, 1, "last");
    A_interp = ops{idx}.A + (mu_new - mus(idx)) * (ops{idx+1}.A - ops{idx}.A) / (mus(idx+1) - mus(idx));
    B_interp = ops{idx}.B + (mu_new - mus(idx)) * (ops{idx+1}.B - ops{idx}.B) / (mus(idx+1) - mus(idx));
    F_interp = ops{idx}.F + (mu_new - mus(idx)) * (ops{idx+1}.F - ops{idx}.F) / (mus(idx+1) - mus(idx));
end


% function [A_interp, B_interp, F_interp] = spline_operators(ops, n, mus, mu_new)
%     % Dimensions
%     M = length(ops);
%     l = size(ops{1}.B,2);
%     s = n * (n+1) / 2;
% 
%     % Unpack matrices 
%     A_all = zeros(n,n,M);
%     B_all = zeros(n,l,M);
%     F_all = zeros(n,s,M);
%     for i = 1:M
%         A_all(:,:,i) = ops{i}.A(1:n,1:n);
%         B_all(:,:,i) = ops{i}.B(1:n,1:l);
%         F_all(:,:,i) = ops{i}.F(1:n,1:s);
%     end
%     A_interp = spline(mus, A_all(:,:,:), mu_new);
%     B_interp = spline(mus, B_all(:,:,:), mu_new);
%     F_interp = spline(mus, F_all(:,:,:), mu_new);
% end
% 
% 
