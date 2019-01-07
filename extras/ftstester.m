clearvars
close all
fts = FourTankSystem;
gamma = [0.65; 0.4]; 
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

% Linearization and discretization
fts.linearize;
fts.discretize(Ts);

% Poles and zeros
pzmap(fts.tfd)
xlim([0.9 1.1])

%% Test for identification
% simulation around equilibrium
fts.m = fts.ms;
fts.record.t = [];
fts.record.x = [];
for it = 1:length(t)-1
    fts.timestep([t(it); t(it+1)],ssInput.*[1.05 1 0 0]);
end

y = fts.record.x'-fts.ms;
y = y(1,:)';

% Plotting
subplot(212)
plot(fts.record.t,y)

% ---------- DELETE THIS AFTER EVERYTHING IS ALIGHT
trafo = tf(4.09e-5,[1 0.059 0.0006484])
[y, t] = step(trafo);
% -------------------------------------

% Initialization of object for id.
idobj = IdentificationObject;
n.np = 3;
n.nx = 2;
idobj.initialize(n);
% Initialization of PI
pi = ParameterIdentifier;
Bounds.LowerBound = [];
Bounds.UpperBound = [];
ODEoptions = [];
OPToptions = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Jacobian','on',...
    'SpecifyObjectiveGradient',true,'FunctionTolerance',1e-10,...
    'StepTolerance',1e-10);
pi.initialize(idobj,n,Bounds,ODEoptions,OPToptions);

% Identification
ResidualArguments.t = t; %fts.record.t;
ResidualArguments.y = y;
ResidualArguments.x0 = [0; 0];
InitialParameters = [0.0001; 0.1; 0.001];
p = pi.identify(ResidualArguments,InitialParameters)
% fts.tfc(1,1).num{1}
% fts.tfc(1,1).den{1}