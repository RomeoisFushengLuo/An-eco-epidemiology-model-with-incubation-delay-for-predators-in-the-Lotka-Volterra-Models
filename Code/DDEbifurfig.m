clear all;
global a b e p k alpha lambda tau m miu d r n

tau = 0.001;
lags = [tau tau tau];
tspan = [0 1000];


parameters = linspace(2.4, 4, 100);  
results = zeros(length(parameters), 1000);  

for i = 1:length(parameters)
    r = parameters(i);
    sol = dde23(@ddefun, lags, @history, tspan);
    results(i, :) = sol.y(1, :);  
end

figure;
plot(parameters, results(:, end), '.');
xlabel('r');
ylabel('Population');
title('Bifurcation diagram of the Logistic map');

function y = ddefun(t, Z)
  k = 50;
  alpha = 0.1;
  tau = 0.05;
  m = 0.01;
  miu = 0.1;
  d = 0.8;
  r = 0.5;
  n = 1;
  b = 5;
  p = 0.1;
  e = 0.3;

  ylag1 = Z(:,1);
  ylag2 = Z(:,2);
  ylag3 = Z(:,3);
  
  y = [r*ylag1*(1-ylag1/k-n/(ylag1+b))-m*ylag1*ylag2-p*m*ylag1*ylag3; 
       e*m*ylag1*ylag2-alpha*ylag2*ylag3-miu*ylag2;
       alpha*ylag2*ylag3+e*m*p*ylag1*ylag3-d*ylag3];
end

function s = history(t)
  s = ones(3, 1); 
  s(1) = 40;
  s(2) = 7;
  s(3) = 0.1;
end
