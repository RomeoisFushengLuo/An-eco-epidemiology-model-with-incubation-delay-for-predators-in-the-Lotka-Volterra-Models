syms a b e p k alpha lambda tau m miu d r n x y z
J = jacobian([r*x*(1-x/k-n/(x+b))-m*x*y-p*m*x*z;e*m*x*y-alpha*y*z-miu*y;a*y*z+e*m*p*x*z-d*z],[x,y,z])
[x,y,z] = solve(r*x*(1-x/k-n/(x+b))-m*x*y-p*m*x*z==0,e*m*x*y-alpha*y*z-miu*y==0,a*y*z+e*m*p*x*z-d*z==0)
% J1=subs(J, {sym('x'), sym('y'), sym('z')},{x(1), y(1),z(1)})
% EigenvaleofJ1=eig(J1)
% J2=subs(J, {sym('x'), sym('y'), sym('z')},{x(2), y(2),z(2)})
