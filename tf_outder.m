function [dgdx,dgdp] = tf_outder(x,p)
K = p(1);
a0 = p(2);
a1 = p(3);
a2 = p(4);
dgdx = [1 0];
dgdp = zeros(1,4);