%Solve the heat equation using crankNicolson.
clear all;



mu = 0.4;
dX = 1/20;
dT = mu * dX^2;

t = 0:dT:0.2;
x = 0:dX:1;

U = meshgrid(x,t);
exSol = zeros(length(t),length(x));
error = zeros(length(t),length(x));

%Time boundary conditions.
U(:,1) = 0;
U(:,end) = 0;
%x boundary condition
U(1,:) = sin(pi*x);
%U(1,:) = sin(5*pi*x/2);


%Crank-Nicolson =>
theta = 0.5;

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
    error(n,:) = abs(U(n,:) -exSol(n,:));
        
end

%plot solution
subplot(2,1,1)
mesh(x,t,U)
xlabel('x')
ylabel('time')
zlabel('function value')
%Plot the error
subplot(2,1,2)
mesh(x,t,error)
xlabel('x')
ylabel('time')
zlabel('abs(error)')
disp('the biggest error is:')
max(error(:))
