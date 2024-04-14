clear all
global b e p k alpha m miu d r n

k = 18.79;
alpha = 0.1;
m = .1;
miu = 0.1;
d = 0.2;
r = 0.95;
n=9;
b=10;
p=0.5;
e=.333333;

% Figure 1 
figure(1);
for a0=10
    b0=5
    c0=5
ts=[0 2000];
z0=[a0 b0 c0];

[t, z] = ode45('epidemic', ts , z0); 


s=plot3(z(:,1),z(:,2),z(:,3))
xlabel('prey'), ylabel('susceptible predator'),zlabel('infected predator')
hold on
end

% Figure 2
ts=[0 2000];
z0=[20 10 1];
[t, z] = ode45('epidemic', ts , z0); 

figure (2)
X=plot(t, z(:, 1),LineWidth=2)
legend(X,'prey');
hold on
Y=plot(t, z(:, 2),LineWidth=2)
legend(Y,'susceptible predator');
hold on
Z=plot(t, z(:, 3),LineWidth=2)
legend([X,Y,Z],'prey','susceptible predator','infected predator')
