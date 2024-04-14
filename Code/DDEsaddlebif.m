clear all;
global a e p alpha lambda tau m miu d r n

tau = 10;
lags = [tau tau tau];
tspan = [0 200];


e_range = linspace(0.3, 0.6, 500);   
miu_range = linspace(0, 0.5, 500);   

bifurcation_data = zeros(length(e_range), length(miu_range));


for i = 1:length(e_range)
    for j = 1:length(miu_range)
        e = e_range(i);
        miu = miu_range(j);
        

        sol = dde23(@(t, y, Z) ddefun(t, y, Z, e, miu), lags, @history, tspan);
        
     
        y_final = sol.y(:, end);
        
       
        J = zeros(length(y_final));
        delta = 1e-6;
        for ii = 1:length(y_final)
            y_perturbed = y_final;
            y_perturbed(ii) = y_perturbed(ii) + delta;
            sol_perturbed = dde23(@(t, y, Z) ddefun(t, y, Z, e, miu), lags, @(t) y_perturbed, [0 tau]);
            J(:, ii) = (sol_perturbed.y(:, end) - y_final) / delta;
        end
        
        
        eig_values = eig(J);
        
       
        if any(real(eig_values) > 0)
            bifurcation_data(i, j) = 1;  
        end
    end
end


figure;
imagesc(miu_range, e_range, bifurcation_data);
colormap("parula");
xlabel('$$\mu$$', 'Interpreter', 'latex');
ylabel('e');
%title('Bifurcation Diagram (Interior Point Analysis)');
colorbar;


function dydt = ddefun(t, y, Z, e, miu)
  alpha = 0.2;
  m =0.1;
  b = 5;
  d=0.8;
  k=50;
  n=1;
  p=0.7;

  ylag1 = Z(:,1);
  ylag2 = Z(:,2);
  ylag3 = Z(:,3);
  dydt = [e*y(1)*(1-y(1)/k-n/(y(1)+b))-m*y(1)*y(2)-p*m*y(1)*y(3); 
          e*m*y(1)*y(2)-alpha*y(2)*y(3)-miu*y(2);
          alpha*ylag2(2)*ylag3(3)+e*m*p*y(1)*y(3)-d*y(3)];
end

function s = history(t)
  s = ones(3, 1); 
  s(1) = 10;
  s(2) = 4;
  s(3) = 2;
end
