clearvars
system.A = [0.9 -1; 0 0.8];
system.B = [0.1 0.1; 2 2];
system.C = eye(2);
system.Cz = [1 0];
system.G = [0.1; 2];
noise.R = 0.1*eye(2);
noise.Q = 0.01*eye(1);
noise.S = 0.1*ones(1,2);
initial.x = zeros(2,1);
initial.P = eye(2);
kfType = 'timeinvariant';
kf = KalmanFilter;
horizon = 1000;
kf.initialize(system,noise,initial,kfType,horizon)
j = horizon;
u = ones(2,j);
Lr = chol(noise.R);
y = initial.x + Lr*randn(2,1);
Q = cell(5,1);
for it = 1:j
    Q{it} = noise.Q;
end
% WARNING: THE MARKOV PREDICTOR STILL DOES NOT WORK WELL
% [xf, x1, xj, z] = kf.outputPredictor(u,y,Q,j)
[xf, x1, z] = kf.markovPredictor(u,y)
% disp(['Final value should be ' num2str(kf.Cz*((eye(2)-kf.A)\kf.B))])