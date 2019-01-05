classdef ModelPredictiveController < matlab.mixin.Copyable
    % Function designed by Zoltan Mark Pinter, DTU student behind s172040,
    %   for Model Predictive Control Course 02619, and his MSc thesis, 2019
    %
    % Model Predictive Controller.
    properties
        % Usage mode
        useInputConstraints
        useInputRateConstraints
        useSoftOutputConstraints
        % Optimization options
        optioptions
        % Auxuliary matrices
        Markov  % Markov matrix
        IN      % Kronecker matrix
        I0      % Picking 0th matrix
        Lambda  % Differentiation matrix
        e       % One-vector
        % Quadratic properties
        Ucoeff  % Coefficients of the input in optimization problem
        W       % Cell of weights
        M       % Cell of slope coefficients
        H       % Hessian
        g       % Gradient
        % Constraints
        A       % Constraint matrix
        bl      % Lower function bound
        bu      % Upper function bound
        l       % Lower independent bound
        u       % Upper independent bound
        % Memory
        u_1     % Previous value of u
        % Lengths
        j       % Length of horizon
        nx      % Number of states
        nz      % Number of outputs
        nu      % Number of inputs
        no      % Number of unconstrained objectives
        nscu    % Number of soft upper output constraints
        nscl    % Number of soft lower output constraints
    end
    methods
        function controlPrepare(mpc,c,ExtraFeatures)
            % Function help: time dependent part of initialization
            
            % Gradients
            mpc.g(1:mpc.j*mpc.nu) = 0;
            for it = 1:mpc.no
                mpc.g(1:mpc.j*mpc.nu) = mpc.g(1:mpc.j*mpc.nu) + ...
                    mpc.M{it}*c{it};
            end
            % Constraints
            mpc.bl = [];
            mpc.bu = [];
            if mpc.useInputConstraints
                mpc.l(1:mpc.j*mpc.nu,1) = ExtraFeatures.Bounds.Umin;
                mpc.u(1:mpc.j*mpc.nu,1) = ExtraFeatures.Bounds.Umax;
            end
            if mpc.useInputRateConstraints
                mpc.bl = ExtraFeatures.Bounds.DUmin + mpc.I0*mpc.u_1;
                mpc.bu = ExtraFeatures.Bounds.DUmax + mpc.I0*mpc.u_1;
            end
            if mpc.useSoftOutputConstraints
                for it = 1:mpc.nscl
                    mpc.bl = [mpc.bl; ExtraFeatures.Bounds.Rmin{it}-ExtraFeatures.b];
                    mpc.bu = [mpc.bu; 1e10*ones(mpc.j*mpc.nz,1)];
                end
                for it = 1:mpc.nscu
                    bl = mpc.bl;
                    bu = mpc.bu;
                    mpc.bl = [bl; -1e10*ones(mpc.j*mpc.nz,1)];
                    mpc.bu = [bu; ExtraFeatures.Bounds.Rmax{it}-ExtraFeatures.b];
                end
            end
        end
        function [u, U] = controlCompute(mpc,c,ExtraFeatures)
            % Calculates control input from the objective coefficients and
            %   the extra features
            
            mpc.controlPrepare(c,ExtraFeatures)
            if mpc.useInputConstraints || mpc.useInputRateConstraints ||...
                    mpc.useSoftOutputConstraints
                U = quadprog(mpc.H,mpc.g',[mpc.A; -mpc.A],[mpc.bu; -mpc.bl],...
                    [],[],mpc.l,mpc.u,mpc.u_1,mpc.optioptions);
            else
                U = -mpc.H\mpc.g;
            end
            U = U(1:mpc.nu*mpc.j);
            u = U(1:mpc.nu,1);
            mpc.u_1 = u;
        end
        function initialize(mpc,W,Ucoeff,n,ExtraFeatures,u_1,optioptions)
            % Function help: initialization of MPC.
            %   Optimization problem is phrased as ||W(Ucoeff*U-c)||_2^2
            %   Soft optimization subproblem minimizes nonnegative
            %       excessions of output constraints
            %   ExtraFeatures include the Markov matrix for soft output
            %       constraints (if used), and booleans stating what kind
            %       of constraints are to be considered.
            %   Parameter u_1 is the previous value of u.
            
            % Descriptive numbers
            mpc.no = n.no;
            mpc.nscu = n.nscu;
            mpc.nscl = n.nscl;
            % Previous value
            mpc.u_1 = u_1;
            % Extra features
            mpc.useInputConstraints = ExtraFeatures.useInputConstraints;
            mpc.useInputRateConstraints = ExtraFeatures.useInputRateConstraints;
            mpc.useSoftOutputConstraints = ExtraFeatures.useSoftOutputConstraints;
            % Coefficient of independent variable in objective function
            mpc.Ucoeff = Ucoeff;
            % Initialization of matrices
            mpc.W = cell(mpc.no+(mpc.nscl+mpc.nscu),1);
            mpc.H = zeros(mpc.j*mpc.nu);
            mpc.M = cell(mpc.no+(mpc.nscl+mpc.nscu),1);
            mpc.g = zeros(mpc.j*(mpc.nu+mpc.nz*(mpc.nscl+mpc.nscu)),1);
            % Unconstrained objectives
            for it = 1:mpc.no
                % Weights
                mpc.W{it} = kron(mpc.IN,W{it});
                % Hessian
                mpc.H = mpc.H +...
                    (mpc.W{it}*mpc.Ucoeff{it})'*(mpc.W{it}*mpc.Ucoeff{it});
                % Gradient coefficients
                mpc.M{it} = -(mpc.W{it}*mpc.Ucoeff{it})'*mpc.W{it};
            end
            % Soft output constrained optimization
            if mpc.useSoftOutputConstraints
                mpc.e = ones(mpc.j*mpc.nu,1);
                for it = mpc.no+1:mpc.no+(mpc.nscl+mpc.nscu)
                    % Weights
                    mpc.W{it} = kron(mpc.IN,W{it});
                    % Hessian
                    mpc.H = blkdiag(mpc.H,mpc.W{it}'*mpc.W{it});
                    % Gradient coefficients
                    mpc.M{it} = mpc.W{it};
                    % Gradients
                    mpc.g(mpc.nu*mpc.j*(it-mpc.no)+1:mpc.nu*mpc.j*(it-mpc.no+1))...
                        = mpc.M{it}*mpc.e;
                end
            end
            % Come change
            % Constraints
            mpc.A = [];
            mpc.l = zeros(mpc.j*(mpc.useInputConstraints*mpc.nu+...
                mpc.nz*(mpc.nscl+mpc.nscu)),1);
            mpc.u = mpc.l;
            if mpc.useInputRateConstraints
                mpc.A = mpc.Lambda;
            end
            if mpc.useSoftOutputConstraints
                mpc.Markov = ExtraFeatures.Markov;
                mpc.A = [mpc.A zeros(length(mpc.A),mpc.j*mpc.nz*(mpc.nscl+mpc.nscu));...
                    kron(ones((mpc.nscl+mpc.nscu),1),mpc.Markov)...
                    blkdiag(eye(mpc.j*mpc.nz*mpc.nscl),-eye(mpc.j*mpc.nz*mpc.nscu))];
                mpc.l(mpc.j*mpc.nu+1:...
                    mpc.j*(mpc.nu+mpc.nz*(mpc.nscl+mpc.nscu)),1) = 0;
                mpc.u(mpc.j*mpc.nu+1:...
                    mpc.j*(mpc.nu+mpc.nz*(mpc.nscl+mpc.nscu)),1) = 1e10;     % Inf
            end
            % Correcting nnumerical errors of the Hessian
            mpc.H = (mpc.H+mpc.H')/2;
            % Optimization options
            mpc.optioptions = optioptions;
        end
        function preinitialize(mpc,horizon,n)
            % Function help: recording system size and horizon and produce auxuilary matrices
            
            % Descriptive numbers
            mpc.j = horizon;
            mpc.nu = n.nu;
            mpc.nx = n.nx;
            mpc.nz = n.nz;
            % Auxuliary matrices
            mpc.IN = eye(mpc.j);
            mpc.I0 = [eye(mpc.nu); zeros(mpc.nu*(mpc.j-1),mpc.nu)];
            mpc.Lambda = kron(mpc.IN,eye(mpc.nu))-...
                [zeros(mpc.nu,mpc.nu*mpc.j);...
                [kron(eye(mpc.j-1),eye(mpc.nu))...
                zeros(mpc.nu*(mpc.j-1),mpc.nu)]];
        end
    end
end