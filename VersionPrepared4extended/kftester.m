clearvars
system.A = [0.9 -1; 0 0.8];
system.B = [0.1; 2];
system.C = eye(2);
system.Cz = [1 0];
system.G = [0.1; 2];
noise.R = 0.1*eye(2);
noise.Q = 0.01*eye(1);
noise.S = zeros(1,2);
initial.x = zeros(2,1);
initial.P = eye(2);
kfType = 'stationary';
kf = KalmanFilter;
kf.initialize(system,noise,initial,kfType)
j = 500;
u = ones(1,j);
Lr = chol(noise.R);
y = initial.x + Lr*randn(2,1);
Q = cell(5,1);
for it = 1:j
    Q{it} = noise.Q;
end
[xf, x1, xj, z] = kf.outputPredictor(u,y,Q,j)