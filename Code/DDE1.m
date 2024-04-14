clear all;
global a b e p k alpha lambda tau m miu d r n

tau =  1;
lags = [tau tau tau];
tspan = [0 2200];

% 2D Figure
figure;
sol = dde23(@ddefun, lags, @history, tspan);
plot(sol.x, sol.y,LineWidth=2)
xlabel('Time t');
ylabel('Solution y');
legend('prey', 'susceptible predator', 'infected predator', 'Location', 'NorthEast');

% 3D Figure
figure;
a0=10
b0=5
c0=5
z0=[a0 b0 c0];
multisol = dde23(@ddefun, lags, z0, tspan);
y1 = multisol.y(1,:);
y2 = multisol.y(2,:);
y3 = multisol.y(3,:);
plot3(y1, y2, y3);
xlabel('prey');
ylabel('susceptible predator');
zlabel('infected predator');
title('Species Dynamics');
grid on;


function dydt = ddefun(t, y, Z)
  k = 50;
  alpha = 0.2;
  tau =  1;
  m = 0.1;
  miu = 0.1;
  d = 0.8;
  r = 0.5;
  n=1;
  b=5;
  p=0.7;
  e=0.4;

  ylag1 = Z(:,1);
  ylag2 = Z(:,2);
  ylag3 = Z(:,3);
  dydt = [r*y(1)*(1-y(1)/k-n/(y(1)+b))-m*y(1)*y(2)-p*m*y(1)*y(3); 
          e*m*y(1)*y(2)-alpha*y(2)*y(3)-miu*y(2);
          alpha*ylag2(2)*ylag3(3)+e*m*p*y(1)*y(3)-d*y(3)];
end


function s = history(t)
  s = ones(3, 1); 
  s(1) = 10;
  s(2) = 5;
  s(3) = 5;
end