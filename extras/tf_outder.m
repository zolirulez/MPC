function [dgdx,dgdp] = tf_outder(x,p)

dgdx = [1 0];
dgdp = zeros(1,3);