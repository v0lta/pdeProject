%This script attempts to solve the convection problem.
%Create a grid with t e [0.05, 0.5]. x e [0.1]
%Stability with mu < 0.5. mu = (delta t)/(delta x)^2
mu = 0.01;
dX = 1/20;
beta = mu * dX/2;
dT = mu * dX^2;
%First requested plot
t = 0:dT:0.01;
%second requested plot.
%t = 0:deltaT:0.05;

x = 0:dX:1;

%Boundary conditions:
U = zeros(length(t),length(x));
U(:,1) = 0;
U(:,end) = 0;

%1.0 for x in [0.450, 0.550]
pos1 = floor(0.450/dX);
pos2 = floor(0.550/dX);
U(1,pos1) = 1;
U(1,pos2) = 1;

%U(1,:) = sin(pi*x);


%Heat equation with convection:
%u_t = 50*ux +  u_xx
%Explicit Euler :
% U_(j+1)^n = U_j^n + mu *( U_(j+1)^n  - 2*U_j^n + U_(j-1)^n )

for n = 1:1:(length(t)-1)
    for j = 2:1:(length(x)-1)
        U(n+1,j) = U(n,j) + mu * ( U(n,j+1) - 2*U(n,j) + U(n,j-1) ) ...
                    + 50*beta*( U(n,j+1) - U(n,j-1) );
    end
end

mesh(x,t,U)
%shading interp
%lighting phong
%surf(x,t,U)
xlabel('x')
ylabel('time')
zlabel('function value')
