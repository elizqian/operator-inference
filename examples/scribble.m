
clear; clc;
addpath('../')

n = 10;
H = zeros(n,n^2);
for i = 1:n
    tmp = randi([-5,5],n);
    tmp = reshape(tmp*tmp', 1, n^2);
    H(i,:) = tmp;
end


Dn = dupmat(n);
H;
F = H*Dn;

H1 = F2H(F);
H2 = F2Hs(F);

% nnz(H - H1)
% nnz(H - H2)
% nnz(H1 - H2)


CR_h = 0;
CR_f = 0;
CR_h_store = [0];
CR_f_store = [0];
for i = 1:n
    for j = 1:i
        for k = 1:j
            tmp1 = delta(j,k)*F(i,fidx(n,j,k)); fprintf("F(%d,%d,%d) = % 3.1f, ", [i,j,k,tmp1]);
            tmp1 = H(i,n*(k-1)+j); fprintf("H(%d,%d,%d) = % 3.1f, ", [i,j,k,tmp1]);

            tmp2 = delta(i,k)*F(j,fidx(n,i,k)); fprintf("F(%d,%d,%d) = % 3.1f, ", [j,i,k,tmp2]);
            tmp2 = H(j,n*(k-1)+i); fprintf("H(%d,%d,%d) = % 3.1f, ", [j,i,k,tmp2]);

            tmp3 = delta(i,j)*F(k,fidx(n,i,j)); fprintf("F(%d,%d,%d) = % 3.1f, ", [k,i,j,tmp3]);
            tmp3 = H(k,n*(j-1)+i); fprintf("H(%d,%d,%d) = % 3.1f\n", [k,i,j,tmp3]);
%             fprintf("F(%d,%d,%d) = % 3.1f, F(%d,%d,%d) = % 3.1f, F(%d,%d,%d) = % 3.1f ", [i,j,k,tmp1,tmp2])
            
            CR_h = CR_h + abs(H(i,n*(k-1)+j) + H(j,n*(k-1)+i) + H(k,n*(j-1)+i));
            CR_f = CR_f + abs(delta(j,k)*F(i,fidx(n,j,k)) + delta(i,k)*F(j,fidx(n,i,k)) + delta(i,j)*F(k,fidx(n,i,j)));

            CR_h_store = [CR_h_store, CR_h];
            CR_f_store = [CR_f_store, CR_f];
        end
    end
end
CR_h
CR_f

figure(1);
plot(CR_h_store, '-', DisplayName="H", LineWidth=4)
hold on;
plot(CR_f_store, DisplayName="F", LineStyle="--", LineWidth=2, Color='g')
hold off; grid on; legend(Location="southeast");


%% 
function del = delta(i,j)
    if i == j
        del = 1.0;
    else
        del = 0.5;
    end
end

function idx = fidx(n,j,k)
%     if j >= k
    idx = (n - k/2)*(k - 1) + j;
%     else
%         idx = (n - j/2)*(j - 1) + k;
%     end
end

