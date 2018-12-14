function [dfdx,dfdp] = tf_der(t,x,p)
K = p(1);
a0 = p(2);
a1 = p(3);
a2 = p(4);
dfdx = [0 1; -a0/a2 -a1/a2];
dfdp = [zeros(1,4); 1/a2 -1/a2*x(1) -1/a2*x(2) 1/a2*(a0*x(1)+a1*x(2))];