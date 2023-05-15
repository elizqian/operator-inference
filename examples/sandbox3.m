clear; close all; clc

N = 3;

a = 2;
b = 3;
c = 4;
e = 6; 
f = 7;
g = 8;

H = [0 a/2 b/2 a/2 -c 0 b/2 0 -e; -a c/2 0 c/2 0 f/2 0 f/2 -g; -b 0 e/2 0 -f g/2 e/2 g/2 0]


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

%% project onto random basis
r = 2;
V = rand(N,r);
Hhat = V'*H*kron(V,V);

% check index-wise constraint
get_hhat = @(i,j,k) Hhat(i,(k-1)*r+j);
derphat = zeros(r,r,r);
for i = 1:r
    for j = 1:r
        for k = 1:r
            derphat(i,j,k) = get_hhat(i,j,k) + get_hhat(j,i,k) + get_hhat(k,j,i);
        end
    end
end

% check inner product
violhat = 0;
for i = 1:10
    randtest = rand(r,1);
    violhat = violhat + abs(randtest'*Hhat*kron(randtest,randtest));
end
violhat
