

%Solve the third question:
h = [1/10 1/20 1/40 1/80];
k = [1/20, 1/40, 1/80, 1/160, 1/320];
tEnd = 0.5;
xEnd = 1;
leftBound = 0;
rightBound = 0;
theta = 0.5;



for a = 1:1:length(h)
    for b = 1:1:length(k)
        dT = h(a);
        dX = k(b);
        mu = dT/(dX^2);
        x = 0:dX:xEnd;
        xBound = sin(pi*x);
        errorTmp = thetaMethod( dT,dX,tEnd,xEnd,leftBound,...
                                 rightBound,xBound,theta);
        errorMat(a,b) = norm(errorTmp(end-1,:));
        
    end
end
figure(1);
bar3(errorMat);
xlabel('grid')% = [1/20, 1/40, 1/80, 1/160, 1/320]');
ylabel('time')% = [1/10 1/20 1/40 1/80]');
zlabel('error')
title('Crank-Nicolson')


theta = 1;
for a = 1:1:length(h)
    for b = 1:1:length(k)
        dT = h(a);
        dX = k(b);
        mu = dT/(dX^2);
        x = 0:dX:xEnd;
        xBound = sin(pi*x);
        errorTmp = thetaMethod( dT,dX,tEnd,xEnd,leftBound,...
                                 rightBound,xBound,theta);
        errorMat(a,b) = norm(errorTmp(end-1,:));
        
    end
end
figure(2);
bar3(errorMat);
xlabel('grid') %= [1/20, 1/40, 1/80, 1/160, 1/320]');
ylabel('time') %= [1/10 1/20 1/40 1/80]');
zlabel('error')
title('Implicit Euler');

%plot mu:
for a = 1:1:length(h)
    for b = 1:1:length(k)
        dT = h(a);
        dX = k(b);
        muMat(a,b) = dT/(dX^2);
    end
end
figure(3);
bar3(muMat);
xlabel('grid')% = [1/20, 1/40, 1/80, 1/160, 1/320]');
ylabel('time')% = [1/10 1/20 1/40 1/80]');
zlabel('magnitude')
title('\mu- Magnitude Plot');

