%This script attempts to find the error of the one dimensional, time dependent heat
%problem by using the explicit euler method.
clear all;
%Create a grid with t e [0.05, 0.5]. x e [0.1]
%Stability with mu < 0.5. mu = (delta t)/(delta x)^2
mu = 0.4;
dX = 1/20;
dT = mu * dX^2;

t = 0:dT:0.2;
x = 0:dX:1;

%Boundary conditions:
U = meshgrid(x,t);
exSol = zeros(length(t),length(x));
error = zeros(length(t),length(x));
%Left time boundary
U(:,1) = 0;
%Right time boundary.
U(:,end) = 0;
%u(x,0) = sin(5*pi*x)/2;
%U(1,:) = sin(pi*x);
U(1,:) = sin(pi*(x));


%Heat equation:
%u_t = u_xx
%Explicit Euler :
% U_(j+1)^n = U_j^n + mu *( U_(j+1)^n  - 2*U_j^n + U_(j-1)^n )


for n = 1:1:(length(t)-1)
    for j = 2:1:(length(x)-1)
        %Find the numerical solution:
        U(n+1,j) = U(n,j) + mu* ( U(n,j+1) - 2*U(n,j) + U(n,j-1));
        
        %Exact solution:
        exSol(n,j) = exp(-pi^2*t(n))*sin(pi*x(j));
        %Error:
        error(n,j) = abs(U(n,j) -exSol(n,j));
    end
end



%plot solution
subplot(1,2,1)
mesh(x,t,U)
xlabel('x')
ylabel('time')
zlabel('function value')
%Plot the error
subplot(1,2,2)
mesh(x,t,error)
xlabel('x')
ylabel('time')
zlabel('abs(error)')
disp('the biggest error is:')
max(error(:))
