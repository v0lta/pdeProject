%This script attempts solves the euqation for the second initial condition
%set and calculates the error in comparison with the exact solution.

%Create a grid with t e [0.05, 0.5]. x e [0.1]
%Stability with mu < 0.5. mu = (delta t)/(delta x)^2
mu = 0.3;
dX = 1/20;
dT = mu * dX^2;
%First requested plot
t = 0:dT:0.5;
%second requested plot.
%t = 0:deltaT:0.05;

x = 0:dX:1;

%Boundary conditions:
U = zeros(length(t),length(x));
%u(0,t) = 0;
U(:,1) = 0;
%u(1,t) = 1;
U(:,end) = 1;
%u(x,0) = sin(5*pi*x)/2;
%U(1,:) = sin(pi*x);
U(1,:) = sin(5*pi*x/2);


%Heat equation:
%u_t = u_xx
%Explicit Euler :
% U_(j+1)^n = U_j^n + mu *( U_(j+1)^n  - 2*U_j^n + U_(j-1)^n )

for n = 1:1:(length(t)-1)
    for j = 2:1:(length(x)-1)
        U(n+1,j) = U(n,j) + mu* ( U(n,j+1) - 2*U(n,j) + U(n,j-1));
    end
end

mesh(x,t,U)
%shading interp
%lighting phong
%surf(x,t,U)
xlabel('x')
ylabel('time')
zlabel('function value')
