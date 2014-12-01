%The script solves our  conveciton problem using the implicit euler method.
%The differential is approximated with a forward difference.
clear all;

mu = 0.01;
dX = 1/20;
dT = mu * dX^2;

t = 0:dT:0.01;
x = 0:dX:1;

U = meshgrid(x,t);
exSol = zeros(length(t),length(x));
error = zeros(length(t),length(x));
%Boundary condition:
%Time boundary conditions.
U(:,1) = 0;
U(:,end) = 0;

%x boundary:
%1.0 for x in [0.450, 0.550]
pos1 = floor(0.450/dX);
pos2 = floor(0.550/dX);
U(1,:) = 0;
U(1,pos1) = 1;
U(1,pos2) = 1;

%U(1,:) = 0;
%U(1,:) = sin(pi*x);


%Heat equation with convection:
%u_t = 50*ux + u_xx
%set dircetion:
dir = -1; %-1 or 1!!!
beta = 50*mu*(dX/2)*dir;

%construct the step Matrix:
stepMat = toeplitz([(1 + 2*mu + beta) -(mu + beta), zeros(1,length(x)-4)],...
             [(1 + 2*mu + beta) (-mu ), zeros(1,length(x)-4)]);

for n = 1:1:(length(t)-1)
    
    %compute the values for the next time step.
    U(n+1,2:(end-1)) = stepMat\U(n,2:(end-1))'...
    + [(mu )*U(n+1,1); zeros(length(x)-4,1); (mu + beta)*U(n+1,end)];
    
    %Exact solution:
    %exSol(n,:) = exp(-pi^2*t(n))*sin(pi*x);
  
end

%plot solution
mesh(x,t,U)
xlabel('x')
ylabel('time')
zlabel('function value')
title('forward difference')