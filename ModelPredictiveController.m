classdef ModelPredictiveController < handle
    % Model Predictive Controller. Handles constrains on input, input
    %   rate, output (for the latter soft constrains as well). The current
    %   version is not prepared for a time varying system.
    properties
        % Markov matrix
        Gamma 
        % Auxuliary matrices
        I0
        Lambda
        % Weights
        Wz
        Wu
        WDu
        % Final quadratic properties
        M       % Slope coefficient
        H       % Hessian
        g       % Gradient
        % Constraints
        A       % Constraint matrix
        bl      % Lower function bound
        bu      % Upper function bound
        l       % Lower independent bound
        u       % Upper independent bound
        % Memory
        U_1     % Previous value of u
        % Length
        nx      % Number of states
        nz      % Number of outputs
        nu      % Number of inputs
        j       % Length of horizon
        no      % Number of unconstrained objectives
        nscu    % Number of soft upper output constraints
        nscl    % Number of soft lower output constraints
    end
    methods
        function controlPrepare(mpc,c,ExtraFeatures,u_1)
            % Gradients
            for it = 1:mpc.no
                mpc.g(1+(it-1)*mpc.j:it*mpc.j) = mpc.M{it}*c{it};
            end
            % Constraints
            mpc.bl = [];
            mpc.bu = [];
            if mpc.useInputConstraints
                mpc.l(1:mpc.j*mpc.nu,1) = ExtraFeatures.Bounds.Umin;
                mpc.u(1:mpc.j*mpc.nu,1) = ExtraFeatures.Bounds.Umax;
            end
            if mpc.useInputRateConstraints
                mpc.bl = ExtraFeatures.Bounds.DUmin + mpc.I0*u_1;
                mpc.bu = ExtraFeatures.Bounds.DUmax + mpc.I0*u_1;
            end
            if mpc.useSoftOutputConstraints
                for it = 1:mpc.n.nscl
                    mpc.bl = [mpc.bl; ExtraFeatures.Bounds.Rmin{it}-ExtraFeatures.b];
                    mpc.bu = [mpc.bu; 1e12*ones(mpc.j*mpc.nz,1)];
                end
                for it = 1:mpc.n.nscu
                    mpc.bl = [mpc.bu; -1e12*ones(mpc.j*mpc.nz,1)];
                    mpc.bu = [mpc.bl; ExtraFeatures.Bounds.Rmax{it}-ExtraFeatures.b];
                end
            end
        end
        function U = controlCompute(mpc,c,ExtraFeatures,optioptions)
            mpc.controlPrepare(c,ExtraFeatures)
            U = quadprog(mpc.H,mpc.g',[mpc.A; -mpc.A],[mpc.bu -mpc.bl],...
                [],[],mpc.l,mpc.u,mpc.U_1,optioptions);
            mpc.U_1 = U;
        end
        function initialize(mpc,W,Ucoeff,horizon,n,ExtraFeatures,U_1)
            % Function help: initialization of MPC.
            %   Optimization problem is phrased as ||W(Ucoeff*U-c)||_2^2
            %   Soft optimization subproblem minimizes nonnegative
            %       excessions of output constraints
            %   ExtraFeatures include the Markov matrix for soft output
            %       constraints (if used), and booleans stating what kind
            %       of constraints are to be considered.
            %   Parameter U_1 is the previous value of u.
            
            % Previous value
            mpc.U_1 = U_1;
            % Descriptive numbers
            mpc.j = horizon;
            mpc.nu = n.nu;
            mpc.nx = n.nx;
            mpc.nz = n.nz;
            mpc.no = n.no;
            mpc.nscu = n.nscu;
            mpc.nscl = n.nscl;
            % Extra features
            mpc.useInputConstraints = ExtraFeatures.useInputConstraints;
            mpc.useInputRateConstraints = ExtraFeatures.useInputRateConstraints;
            mpc.useSoftOutputConstraints = ExtraFeatures.useSoftOutputConstraints;
            % Coefficient of independent variable in objective function
            mpc.Ucoeff = Ucoeff;
            % Initialization of matrices
            IN = eye(mpc.j);
            mpc.W = cell(mpc.no+(mpc.nscl+n.nscu));
            mpc.H = zeros(mpc.j*mpc.nu);
            mpc.M = cell(mpc.no+(mpc.nscl+n.nscu));
            mpc.g = zeros(mpc.j*(1+(mpc.nscl+n.nscu)),1);
            % Unconstrained objectives
            for it = 1:mpc.no
                % Weights
                mpc.W{it} = kron(IN,W{it});
                % Hessian
                mpc.H = mpc.H +...
                    (mpc.W{it}*mpc.Ucoeff{it})'*(mpc.W{it}*mpc.Ucoeff{it});
                % Gradient coefficients
                mpc.M{it} = -(mpc.W{it}*mpc.Ucoeff{it})'*mpc.W{it};
            end
            % Soft output constrained optimization
            if mpc.useSoftOutputConstraints
                mpc.e = ones(mpc.j*mpc.nu,1);
                for it = mpc.no+1:mpc.no+(mpc.nscl+n.nscu)
                    % Weights
                    mpc.W{it} = kron(IN,W{it});
                    % Hessian
                    mpc.H = blkdiag(mpc.H,mpc.W{it}'*mpc.W{it});
                    % Gradient coefficients
                    mpc.M{it} = mpc.W{it};
                    % Gradients
                    mpc.g(mpc.no*mpc.j+1+(it-1)*mpc.j:mpc.no*mpc.j+it*mpc.j)...
                        = mpc.M{it}*mpc.e;
                end
            end
            % Constraints
            mpc.A = [];
            mpc.l = zeros(mpc.j*(mpc.useInputConstraints*mpc.nu+...
                mpc.nz*(mpc.nscl+n.nscu)));
            mpc.u = mpc.l;
            if mpc.useInputRateConstraints
                mpc.I0 = [eye(mpc.nu); zeros(mpc.nu*(mpc.j-1),mpc.nu)];
                mpc.Lambda = kron(IN,eye(mpc.nu))-...
                    [zeros(mpc.nu,mpc.nu*mpc.j);...
                    [kron(eye(mpc.j-1),eye(mpc.nu))...
                    zeros(mpc.nu*(mpc.j-1),mpc.nu)]];
                mpc.A = mpc.Lambda;
            end
            if mpc.useSoftOutputConstraints
                mpc.Markov = ExtraFeatures.Markov;
                mpc.A = [mpc.A zeros(length(mpc.A),mpc.j*mpc.nz*(mpc.nscl+n.nscu));...
                    kron(ones(mpc.j*(mpc.nscl+n.nscu),1),mpc.Markov)...
                    blkdiag(eye(mpc.j*mpc.nz*mpc.nscl),-eye(mpc.j*mpc.nz*mpc.nscu))];
                mpc.l(mpc.j*mpc.nu+1:...
                    mpc.j*(mpc.nu+mpc.nz*(mpc.nscl+n.nscu)),1) = 0;
                mpc.u(mpc.j*mpc.nu+1:...
                    mpc.j*(mpc.nu+mpc.nz*(mpc.nscl+n.nscu)),1) = 1e12; % Inf
            end
        end
    end
end