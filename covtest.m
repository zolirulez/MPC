QQQ = [1 0.5; 0.5 2];
qqq = chol(QQQ);
sss = qqq'*randn(2,10000);
cov(sss')