function [R,ind] = ddt(X,dt,ddt_order)
% Implentation of several difference approximations to d/dt
% see % see https://en.wikipedia.org/wiki/Finite_difference_coefficient
%
%

[N,K] = size(X);
R = zeros(size(X));


switch ddt_order
    
    case {1,'FE'}   %         
%         R = ([X(2:end,:)-X(1:end-1,:); X(end,:)-X(end-1,:)])/dt;
        R = (X(:,2:end) - X(:,1:end-1))/dt;
        ind = 1:K-1;
        
    case {'BE'}
        R = (X(:,2:end) - X(:,1:end-1))/dt;
        ind = 2:K;
    case 12  % forward differencing
        R = ([X(2,:)-X(1,:); X(2:end,:)-X(1:end-1,:)])/dt;
        
    case 1  % forward differencing
        R = [( X(:,2:end)-X(:,1:end-1)), X(:,end)-X(:,end-1)]/dt;
        
    case 2 % central differencing
        R(1,:) = (-3*X(1,:) + 4*X(2,:) - X(3,:))/(2*dt);
        R(end,:) = (3*X(end,:) - 4*X(end-1,:) + X(end-2,:))/(2*dt);
        R(2:end-1,:) = (X(3:end,:) - X(1:end-2,:))/(2*dt);
        
    case 3
        error('3rd order implementation is not done')
    
    case 41 % uses mostly central differencing (same as Renee's code); % 4th order 5-point stencil from https://en.wikipedia.org/wiki/Five-point_stencil
        K = length(X);
        R(3:end-2,:) = (1/12*X(1:end-4,:) - 2/3*X(2:end-3,:) + 2/3*X(4:end-1,:) - 1/12*X(5:end,:))/dt;
        R(1:2,:) = (-25/12*X(1:2,:) + 4*X(2:3,:) - 3*X(3:4,:) + 4/3*X(4:5,:) - 1/4*X(5:6,:))/dt;
        R(K-1:K,:) = (25/12*X(K-1:K,:) - 4*X(K-2:K-1,:) +3*X(K-3:K-2,:) -4/3*X(K-4:K-3,:) + 1/4*X(K-5:K-4,:))/dt;
   
    case 42 % uses mostly backwards differencing
        R(5:end,:) = (25/12*X(5:end,:) - 4*X(4:end-1,:) + 3*X(3:end-2,:) - 4/3*X(2:end-3,:) + 1/4*X(1:end-4,:))/dt;
        R(1:4,:) = (-25/12*X(1:4,:) + 4*X(2:5,:) - 3*X(3:6,:) + 4/3*X(4:7,:) - 1/4*X(5:8,:))/dt;
    
    case 43 % uses mostly forwards differencing   
        R(end-3:end,:) = (25/12*X(end-3:end,:) - 4*X(end-4:end-1,:) + 3*X(end-5:end-2,:) - 4/3*X(end-6:end-3,:) + 1/4*X(end-7:end-4,:))/dt;
        R(1:end-4,:) = (-25/12*X(1:end-4,:) + 4*X(2:end-3,:) - 3*X(3:end-2,:) + 4/3*X(4:end-1,:) - 1/4*X(5:end,:))/dt;

end