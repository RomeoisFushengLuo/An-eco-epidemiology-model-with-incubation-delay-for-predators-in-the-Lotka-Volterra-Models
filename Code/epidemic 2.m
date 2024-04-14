function zdot=epidemic(t,z);
global b e p k alpha m miu d r n
zdot= [r*z(1)*(1-z(1)/k-n/(z(1)+b))-m*z(1)*z(2)-p*m*z(1)*z(3);e*m*z(1)*z(2)-alpha*z(2)*z(3)-miu*z(2);alpha*z(2)*z(3)+e*m*p*z(1)*z(3)-d*z(3)];