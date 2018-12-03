clearvars
fts = FourTankSystem;
gamma = [0.65; 0.55];
initials.m = zeros(4,1);
ODEoptions = [];
fts.initialize(gamma,initials,ODEoptions);
ssInput  = [250; 325; 0; 0];     
% fts.timestep([0; 1],ssInput);
initialGuess = [5000; 5000; 5000; 5000];  
fts.steadyState(ssInput,initialGuess)

% Simulation
Ts = 1;
itmax = 500;
t = 0:Ts:itmax*Ts;
for it = 1:length(t)-1
    fts.timestep([t(it); t(it+1)],ssInput);
end

% Plotting
figure(1)
plot(fts.record.t,fts.record.x)

% Linearization and discretization
fts.linearize;
fts.discretize(Ts);