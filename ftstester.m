clearvars
close all
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

% The DISTURBANCES ARE NOT CONSIDERED; TODO

%% Test for identification
% simulation around equilibrium
fts.m = fts.ms;
fts.record.t = [];
fts.record.x = [];
for it = 1:length(t)-1
    fts.timestep([t(it); t(it+1)],ssInput.*[1.05 1 0 0]);
end
% Plotting
figure(2)
subplot(211)
plot(fts.record.t,fts.record.x)

y = fts.record.x'-fts.ms;
y = y(1,:)';

% Plotting
subplot(212)
plot(fts.record.t,y)

% ---------- DELETE THIS AFTER EVERYTHING IS ALIGHT
trafo = tf(2,[1 1 1]);
[y, t] = step(trafo);
% -------------------------------------

% Initialization
pi = ParameterIdentifier;
Model = @tf_function;
PartialDerivatives = @tf_der;
Output = @tf_out;
OutputPartialDerivatives = @tf_outder;
Bounds.LowerBound = [];
Bounds.UpperBound = [];
ODEoptions = [];
OPToptions = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Jacobian','on',...
    'SpecifyObjectiveGradient',true,'FunctionTolerance',1e-10,...
    'StepTolerance',1e-10);
pi.initialize(Model,PartialDerivatives,Output,OutputPartialDerivatives,...
    Bounds,ODEoptions,OPToptions);

% Identification
ResidualArguments.t = t; %fts.record.t;
ResidualArguments.y = y;
ResidualArguments.np = 4;
ResidualArguments.nx = 2;
ResidualArguments.x0 = [0; 0];
InitialParameters = [2; 1; 1; 2];
p = pi.identify(ResidualArguments,InitialParameters)
fts.tfc(1,1).num{1}
fts.tfc(1,1).den{1}