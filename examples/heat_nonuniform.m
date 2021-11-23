clear; close all

rng(13)

% set up non-uniform grid
dx1 = 0.05;
dx2 = 0.01;
dx3 = 0.1;
x = [0:dx1:0.3, 0.3+dx2:dx2:0.7, 0.8:dx3:1]';

% weights using trapezoid integration rule
w = zeros(size(x));
w(1:end-1) = w(1:end-1) + (x(2:end)-x(1:end-1))/2;
w(2:end) = w(2:end) + (x(2:end)-x(1:end-1))/2;
W = diag(w);

%% draw initial conditions from Gaussian random field
% set up Gaussian random field
s = 3/2;
lmax = 60;
l = 1:lmax;
a = (1./(l*pi)).^s;
X = sin(pi*l.*x).*a;    % this is a sin basis multiplied by decaying weights

num_draws = 1000;    
num_train = 500;
xi = randn(lmax,num_draws);
S0 = sqrt(2)*X*xi;              % to sample from GRF, multiply basis by random normals

%%
figure(1);
plot(x,S0(:,1:5),'-o','MarkerSize',3)
ax = gca;
ax.FontSize = 16; 
xlabel('$x$','interpreter','latex','fontsize',20)
title('Gaussian RF initial states','interpreter','latex','fontsize',20)



figure(2);
plot(x,sqrt(2)*sin(l(1:5)*pi.*x),'MarkerSize',3)
ax = gca;
ax.FontSize = 16; 
ylim([-2 2])
xlabel('$x$','interpreter','latex','fontsize',20)
title('Analytical principal components','interpreter','latex','fontsize',20)



%% evolve for timestep data
S = [];
Sdot = [];
for i = 1:num_draws
    [s, sdot] = timestep(S0(:,i),x);
    S = [S, s];
    Sdot = [Sdot, sdot];
end

r_vals = 1:10;
err = zeros(length(r_vals),4);
derr = zeros(length(r_vals),4);

%% unweighted POD
[U,Sig,V] = svd(S(:,1:num_train*10));
for rr = r_vals
    % get POD basis, project data
    r = r_vals(rr);
    Ur = U(:,1:r);
    Shat = Ur'*S(:,1:num_train*10);
    Shatdot = Ur'*Sdot(:,1:num_train*10);
    
    % LS op inf problem
    D = Shat';
    R = Shatdot';
    Ahat = (D\R)';
    
    ShatOI = [];
    for i = 1:num_draws
        s = timestepA(Ur'*S0(:,i),Ahat);
        ShatOI = [ShatOI, s];
    end
    SOI = Ur*ShatOI;
    
    err(rr,1) = mre(SOI(:,1:num_train*10),S(:,1:num_train*10),W);
    derr(rr,1) = mre(SOI(:,1:num_train*10),S(:,1:num_train*10));
    err(rr,2) = mre(SOI(:,num_train*10+1:end),S(:,num_train*10+1:end),W);
    derr(rr,2) = mre(SOI(:,num_train*10+1:end),S(:,num_train*10+1:end));
end

%% weighted POD
[Util,Sigtil,Vtil] = svd(W^0.5*S(:,1:num_train*10));
[boo,foo] = eig(S*S'*W);
UW = W^-0.5*Util;

for rr = r_vals
    % get POD basis, project data
    r = r_vals(rr);
    Ur = UW(:,1:r);
    Shat = Ur'*W*S(:,1:num_train*10);
    Shatdot = Ur'*W*Sdot(:,1:num_train*10);
    
    % LS op inf problem
    D = Shat';
    R = Shatdot';
    Ahat = (D\R)';
    
    ShatOI = [];
    for i = 1:num_draws
        s = timestepA(Ur'*W*S0(:,i),Ahat);
        ShatOI = [ShatOI, s];
    end
    SOI = Ur*ShatOI;
    
    err(rr,3) = mre(SOI(:,1:num_train*10),S(:,1:num_train*10),W);
    derr(rr,3) = mre(SOI(:,1:num_train*10),S(:,1:num_train*10));
    err(rr,4) = mre(SOI(:,num_train*10+1:end),S(:,num_train*10+1:end),W);
    derr(rr,4) = mre(SOI(:,num_train*10+1:end),S(:,num_train*10+1:end));
end
%%
figure(3);
plot(x,U(:,1:5).*(U(2,1:5)>0) - U(:,1:5).*(U(2,1:5)<0)) %this is to plot everything with a sign convention that matches the sin functions
ax = gca;
ax.FontSize = 16; 
xlabel('$x$','interpreter','latex','fontsize',20)
title('Euclidean POD modes','interpreter','latex','FontSize',20)


figure(4);
plot(x,UW(:,1:5).*(UW(2,1:5)>0) - UW(:,1:5).*(UW(2,1:5)<0))
ax = gca;
ax.FontSize = 16; 
xlabel('$x$','interpreter','latex','fontsize',20)
title('$L^2([0,1])$ POD modes','interpreter','latex','fontsize',20)

%%
figure(5); clf
semilogy(r_vals,err(:,1),'k:'); hold on
semilogy(r_vals,err(:,3))
ax = gca;
ax.FontSize = 16; 
grid on
xlabel('POD basis size $r$','interpreter','latex','fontsize',20)
ylabel('Mean relative error','interpreter','latex','fontsize',20)
ylim([3e-3 1])
legend({'Euclidean','$L^2([0,1])$'},'interpreter','latex','fontsize',18); legend boxoff
title('$L^2([0,1])$-norm training error','interpreter','latex','fontsize',20)


figure(51); clf
semilogy(r_vals,err(:,2),'k:'); hold on
semilogy(r_vals,err(:,4))
ax = gca;
ax.FontSize = 16; 
grid on
ylim([3e-3 1])
xlabel('POD basis size $r$','interpreter','latex','fontsize',20)
ylabel('Mean relative error','interpreter','latex','fontsize',20)
legend({'Euclidean','$L^2([0,1])$'},'interpreter','latex','fontsize',18); legend boxoff
title('$L^2([0,1])$-norm test error','interpreter','latex','fontsize',20)

%%
figure(6); clf
semilogy(r_vals,derr(:,1),'k:'); hold on
semilogy(r_vals,derr(:,3))
ax = gca;
ax.FontSize = 16; 
xlabel('POD basis size $r$','interpreter','latex','fontsize',16)
ylabel('Mean relative error','interpreter','latex','fontsize',16)
legend({'Euclidean','$L^2([0,1])$'},'interpreter','latex','fontsize',14); legend boxoff
title('Training error in Euclidean norm','interpreter','latex','fontsize',18)

figure(61); clf
semilogy(r_vals,derr(:,2),'k:'); hold on
semilogy(r_vals,derr(:,4))
ax = gca;
ax.FontSize = 16; 
xlabel('POD basis size $r$','interpreter','latex','fontsize',16)
ylabel('Mean relative error','interpreter','latex','fontsize',16)
legend({'Euclidean','$L^2([0,1])$'},'interpreter','latex','fontsize',14); legend boxoff
title('Test error in Euclidean norm','interpreter','latex','fontsize',18)

function out = mre(S,S_ref,W)
    if nargin == 2
        out = mean( colnorm(S-S_ref)./colnorm(S_ref) );
    else
        out = mean(colnorm(S-S_ref,W)./colnorm(S_ref,W));
    end
end

function nm = colnorm(S,W)
    if nargin == 1
        nm = sqrt(sum(S.^2));
    else
        nm = sqrt(sum((W^0.5*S).^2));
    end
end

function [s, sdot] = timestep(u,x)
dt = 1e-4;
s = [];
sdot = [];
for i = 1:100
    dx = x(2:end)-x(1:end-1);
    ux = (u(2:end) - u(1:end-1))./dx;
    udot = [0; (ux(2:end)-ux(1:end-1))./(x(3:end)-x(1:end-2)); 0];
    if mod(i,10)==0
        s = [s,u];
        sdot = [sdot, udot];
    end
    u = u + dt*udot;
end
end

function s = timestepA(u,A)
dt = 1e-4;
s = [];
for i = 1:100
    if mod(i,10) == 0
        s = [s,u];
    end
    u = u + dt*A*u;
end
end