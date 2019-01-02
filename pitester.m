
clearvars

% Creating a transfer function
trafo = tf(1,[10 1]);
[y, t] = step(trafo);

% Initialization
pi = ParameterIdentifier;
Model = @tf_function;
PartialDerivatives = @tf_der;
Output = @tf_out;
OutputPartialDerivatives = @tf_outder;
Bounds.LowerBound = [];
Bounds.UpperBound = [];
ODEoptions = [];
optioptions = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Jacobian','on',...
    'SpecifyObjectiveGradient',true,'FunctionTolerance',1e-10,...
    'StepTolerance',1e-3);
pi.initialize(Model,PartialDerivatives,Output,OutputPartialDerivatives,...
    Bounds,ODEoptions,optioptions);

% Identification
ResidualArguments.t = t;
ResidualArguments.y = y;
ResidualArguments.np = 2;
ResidualArguments.nx = 1;
ResidualArguments.x0 = 0;
InitialParameters = [0.5; 0.5];
p = pi.identify(ResidualArguments,InitialParameters)
