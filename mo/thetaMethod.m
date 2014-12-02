function  error  = thetaMethod( dT,dX,tEnd,xEnd,leftBound,...
                                 rightBound,xBound,theta)
%This function comutes the error of a given theta method. For theta = 0,
%the theta method is the explicit Euler method. Foth theta = 0.5 it computes
%the error of the Crank Nicolson method, for theta = 1 of the implicit
%Euler method.
mu = dT/(dX^2);

t = 0:dT:tEnd;
x = 0:dX:xEnd;

U = meshgrid(x,t);
exSol = zeros(length(t),length(x));
error = zeros(length(t),length(x));

%Time boundary conditions.
U(:,1) = leftBound;
U(:,end) = rightBound;
%x boundary condition
U(1,:) = xBound;

%construct the step Matrix:
leftMat = toeplitz([(1 + 2*mu*theta) -mu*theta, zeros(1,length(x)-4)],...
             [(1 + 2*mu*theta) -mu*theta, zeros(1,length(x)-4)]);

rightMat = toeplitz([(1 - 2*mu*(1-theta)) mu*(1-theta), zeros(1,length(x)-4)],...
             [(1 - 2*mu*(1-theta)) mu*(1-theta), zeros(1,length(x)-4)]);
 
         
for n = 1:1:(length(t)-1)
    %compute the values for the next time step.
    U(n+1,2:(end-1)) = leftMat\(rightMat*U(n,2:end-1)')...
                       + mu*[U(n+1,1); zeros(length(x)-4,1); U(n+1,end)];
    %Exact solution:
    exSol(n,:) = exp(-pi^2*t(n))*sin(pi*x);
    %Error:
    error(n,:) = U(n,:) - exSol(n,:);
        
end
end

