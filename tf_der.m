function [dfdx,dfdp] = tf_der(t,x,p)
dfdx = -p(2:3);
dfdp = [1 -x(1) -x(2)];