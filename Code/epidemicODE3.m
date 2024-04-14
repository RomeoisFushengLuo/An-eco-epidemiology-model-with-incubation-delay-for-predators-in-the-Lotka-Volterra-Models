clear all
global b e p k alpha m miu d r n

k = 9;
alpha = 0.0002;
m = .01;
miu = 0.0031;
d = 0.99;
r = 0.1;
n=.1;
b=7;
p=0.01;
e=.8;

% Figure 1 
figure(1);
for a0=1.16:1.29
    b0=a0-.3
    c0=2-a0
ts=[0 50000];
z0=[a0 b0 c0];

[t, z] = ode45('epidemic', ts , z0); 


s=plot3(z(:,1),z(:,2),z(:,3))
xlabel('prey'), ylabel('susceptible predator'),zlabel('infected predator')

hold on
end

% Figure 2
ts=[0 50000];
z0=[1.2499 0.1 0.01];
[t, z] = ode45('epidemic', ts , z0); 

figure (2)
X=plot(t, z(:, 1))
legend(X,'prey');
hold on
Y=plot(t, z(:, 2))
legend(Y,'susceptible predator');
hold on
Z=plot(t, z(:, 3))
legend([X,Y,Z],'prey','susceptible predator','infected predator')
xlabel('t'),ylabel('population')