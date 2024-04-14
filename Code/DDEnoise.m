clear all;

global tau

tau = 0.001;
tspan = [0 10];
dt = 0.02;

sigma = 0.00001;

a0 = 10;
b0 = 4;
c0 = 2;

num_steps = ceil((tspan(2) - tspan(1)) / dt);
y = zeros(3, num_steps);
y(:, 1) = [a0; b0; c0];
timesteps = zeros(1, num_steps);
timesteps(1) = tspan(1);

for i = 1:num_steps-1
    t = tspan(1) + (i-1) * dt;
    timesteps(i+1) = t + dt;
    ylag1 = y(:, max(1, i - ceil(tau / dt)));
    ylag2 = y(:, max(1, i - ceil(tau / dt)));
    ylag3 = y(:, max(1, i - ceil(tau / dt)));
    
    dW = sigma * randn(3,1);
    
    dydt = ddefun(t, y(:, i), [ylag1, ylag2, ylag3]);
    
    y(:, i+1) = y(:, i) + dt * dydt + dW;
end

figure;
plot3(y(1,:), y(2,:), y(3,:));
xlabel('prey');
ylabel('susceptible predator');
zlabel('infected predator');
title('Species Dynamics');
grid on;

figure;
plot(timesteps, y(1,:), 'b-', 'LineWidth', 2); hold on;
plot(timesteps, y(2,:), 'r--', 'LineWidth', 2); hold on;
plot(timesteps, y(3,:), 'g-.', 'LineWidth', 2);
xlabel('Time');
ylabel('Population');
title('Species Dynamics over Time');
legend('Prey', 'Susceptible predator', 'Infected predator');
grid on;

function dydt = ddefun(t, y, Z)
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

    dydt = [r * y(1) * (1 - y(1) / k - n / (y(1) + b)) - m * y(1) * y(2) - p * m * y(1) * y(3);
            e * m * y(1) * y(2) - alpha * y(2) * y(3) - miu * y(2);
            alpha * ylag2(2) * ylag3(3) + e * m * p * y(1) * y(3) - d * y(3)];
end
