function Dx = tf_function(t,x,p)
K = p(1);
a0 = p(2);
a1 = p(3);
a2 = 1;
Dx = [0 1; -a0/a2 -a1/a2]*x + [0; K/a2];