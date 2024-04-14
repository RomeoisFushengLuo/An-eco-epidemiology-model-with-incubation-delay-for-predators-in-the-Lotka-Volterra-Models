 d = 0.15;
x= linspace(10,30,1000);
x1=linspace(0.5,1,1000);
x2=linspace(1,1.5,1000);
x4 = ones(1,1000);
p2 = ones(1,1000);
p3 = ones(1,1000);
del = ones(1,1000);
p4 = zeros(1,1000);
for n = 1:1000
    %del(n)=sqrt((b+0.85*x(n))^2-4*x(n)*9);
    x4(n)= (p*miu+(alpha*r/m)*(1-miu/(e*m*x(n))-9/((miu/(e*m))+b)))/d;
    p2(n) = (0.9*x(n)-b)/2+(sqrt((b+0.9*x(n))^2-4*x(n)*9))/2;
    p3(n)=(0.9*x(n)-b)/2-(sqrt((b+0.9*x(n))^2-4*x(n)*9))/2;
end
z2=(e*m*p2-miu)/alpha;
z3=(e*m*p3-miu)/alpha;

figure;
plot(x,z2,'LineWidth',2)
hold on
plot(x,z3,'LineWidth',2)


figure;
plot(x4,z2,'Color',[0.5 0.7 1],'LineWidth',2)
hold on
plot(x4,z3,'LineWidth',2) 
hold on
plot(x1,p4,'LineWidth',3)
hold on
plot(x2,p4,'--','Color',[0.91 0.588 0.478],'LineWidth',3)

set(gca,'FontSize',13)
xlim([0.5 1.5])
ylim([0 2.6])
xlabel('R_0','fontsize',15);
ylabel('z','fontsize',15);