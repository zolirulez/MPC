function Dx = tf_function(t,x,p)
Dx = -p(2:3).*x+p(1);