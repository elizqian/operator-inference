clear; close all; clc

N = 2;

a = 2;
b = 4;

H = [0 a/2 a/2 -b; -a b/2 b/2 0]


% check index-wise constraint
get_h = @(i,j,k) H(i,(k-1)*N+j);
derp = zeros(N,N,N);
for i = 1:N
    for j = 1:N
        for k = 1:N
            derp(i,j,k) = get_h(i,j,k) + get_h(j,i,k) + get_h(k,j,i);
        end
    end
end

% check inner product
viol = 0;
for i = 1:10
    randtest = rand(N,1);
    viol = viol + abs(randtest'*H*kron(randtest,randtest));
end
viol

V = rand(N,1);
Hhat = V'*H*kron(V,V)