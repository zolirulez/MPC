classdef KalmanFilter < handle
    % This Kalman filter is an object for solving general computations
    %   required by a Kalman Filter. Note that the object does not work
    %   like a single function, therefore for different applications
    %   different instances are needed (e.g. different predictor and
    %   filter.
    properties
        % Types: stationary, timevariant, timeinvariant
        kalmanFilterType
        % Signals
        x1   % Prediction for one time step
        xf   % Estimate
        xj   % Predictions for j time steps
        wf   % Estimate of process noise
        y    % Measurement
        z    % Output
        % Covariance and filtering matrices
        R    % Covariance of measurement noise
        Q    % Covariance of process noise
        S    % Covariance between process and measurement noise
        Pf   % State covariance at estimation
        P1   % State covariance at prediction for one time step
        Pj   % State covariance at prediction for j time steps
        P1s  % Stationary state covariance at prediction for one timestep
        Kx   % Kalman gain for state
        Kw   % Kalman gain for process noise
        Rz   % Covariance of output
        Ry   % Covariance of measurement
        % System matrices
        A    % State matrix
        G    % Noise input matrix
        B    % Input matrix
        C    % Measurement matrix
        Cz   % Output matrix
        obs  % Extended observability matrix
        obsn % Noise propagation matrix
        markov % Markov parameter matrix
        % Length
        nx   % Number of states
        nz   % Number of outputs
        nu   % Number of inputs
        j    % Length of horizon
    end
    methods
        function measurementUpdate(kf,y)
            % Innovation covariance
            Re = kf.C*kf.P1*kf.C' + kf.R;
            % Kalman filter
            if ~strcmp(kf.kalmanFilterType,'stationary')
                kf.Kx = kf.P1*kf.C'/Re;
                kf.Kw = kf.S/Re;
            end
            % Signals
            e = y - kf.C*kf.x1;
            kf.xf = kf.x1 + kf.Kx*e;
            kf.wf = kf.Kw*e;
            % Covariances
            kf.Pf = kf.P1 - kf.Kx*Re*kf.Kx';
            kf.Q = kf.Q - kf.Kw*Re*kf.Kw';
        end
        function timeUpdate(kf,u)
            % Signals
            kf.x1 = kf.A*kf.xf + kf.B*u(:,1) + kf.G*kf.wf;
            % Covariances
            kf.P1 = kf.A*kf.Pf*kf.A' + kf.G*kf.Q*kf.G'...
                - kf.A*kf.Kx*kf.S'*kf.G' - kf.G*kf.S*kf.Kx'*kf.A';
        end
        function statePrediction(kf,u,Q,j)
            % Signals
            kf.xj(:,1) = kf.A*kf.xf + kf.B*u(:,1);
            % Covariances
            kf.Pj{1} = kf.A*kf.Pf*kf.A' + kf.G*kf.Q*kf.G';
            for it = 2:j
                % Signals
                kf.xj(:,it) = kf.A*kf.xj(:,it-1) + kf.B*u(:,it);
                % Covariances
                kf.Pj{it} = kf.A*kf.Pj{it-1}*kf.A' + kf.G*Q{it}*kf.G';
            end
        end
        function outputPrediction(kf,j)
            % Signals
            kf.z(:,1) = kf.Cz*kf.xj(:,1);
            % Covariances
            kf.Rz{1} = kf.Cz*kf.Pj{1}*kf.Cz';
            for it = 2:j
                % Signals
                kf.z(:,it) = kf.Cz*kf.xj(:,it);
                % Covariances
                kf.Rz{it} = kf.Cz*kf.Pj{it}*kf.Cz';
            end
        end
        function measurementPrediction(kf)
            % Signals
            kf.y = kf.C*kf.xf;
            % Covariances
            kf.Ry = kf.C*kf.Pf*kf.C' + kf.R;
        end
        function markovPrediction(kf,u)
            u = reshape(u,kf.j*kf.nu,1);
            % Signal
            kf.z = kf.obs*kf.xf + kf.obsn*kf.wf + kf.markov*u;
        end
        function [xf, x1, xj, z] = outputPredictor(kf,u,y,Q,j)
            kf.measurementUpdate(y);
            if strcmp(kf.kalmanFilterType,'timevariant')
                kf.timeVariation(kf,t);
            end
            kf.timeUpdate(u);
            kf.statePrediction(u,Q,j);
            kf.outputPrediction(j);
            xf = kf.xf;
            x1 = kf.x1;
            xj = kf.xj;
            z = kf.z;
        end
        function [xf, x1, z] = markovPredictor(kf,u,y)
            kf.measurementUpdate(y);
            kf.timeUpdate(u);
            kf.markovPrediction(u);
            xf = kf.xf;
            x1 = kf.x1;
            z = kf.z;
        end
        function timeVariation(kf,t)
            % Function for the case of changing system matrices. The
            %     methodology is not clarified yet.
            
            kf.A = 0;
            kf.B = 0;
            kf.C = 0;
            kf.Cz = 0;
            kf.R = 0;
            kf.Q = 0;
            kf.S = 0;
        end
        function initialize(kf,system,noise,initial,kalmanFilterType,horizon)
            kf.kalmanFilterType = kalmanFilterType;
            kf.A = system.A;
            kf.B = system.B;
            kf.C = system.C;
            kf.Cz = system.Cz;
            kf.G = system.G;
            kf.R = noise.R;
            kf.Q = noise.Q;
            kf.S = noise.S;
            kf.x1 = initial.x;
            kf.P1 = initial.P;
            % Length of vectors
            kf.nu = size(kf.B,2);
            kf.nx = length(kf.A);
            kf.nz = size(kf.Cz,1);
            kf.j = horizon;
            if strcmp(kf.kalmanFilterType,'stationary')
                kf.stationaryInitialization();
            end
            if nargin>4 % horizon is preset
                kf.markovInitialization();
            end
        end
        function stationaryInitialization(kf)
            kf.P1s = dare(kf.A',kf.C',kf.G*kf.Q*kf.G',kf.R,kf.G*kf.S);
            Res = kf.C*kf.P1s*kf.C' + kf.R;
            kf.Kx = kf.P1s*kf.C'/Res; % TODO
            kf.Kw = kf.S/Res; % TODO
        end
        function markovInitialization(kf)
            obs = zeros(kf.nz*(kf.j+1),kf.nx);
            kf.markov = zeros(kf.nz*kf.j,kf.nu*kf.j);
            obs(1:kf.nz,:) = kf.Cz;
            for it = 1:kf.j
                % Extended observability matrix
                obs(it*kf.nz+1:(it+1)*kf.nz,:) =...
                    obs((it-1)*kf.nz+1:it*kf.nz,:)*kf.A;
            end
            H = obs(1:kf.j*kf.nz,:)*kf.B;
            for it = 1:kf.j
                % Markov parameter matrix
                kf.markov(end-it*kf.nz+1:end,...
                    end-it*kf.nu+1:end-(it-1)*kf.nu) = H(1:it*kf.nz,:); 
            end
            kf.obsn = obs(1:kf.j*kf.nz,:)*kf.G;
            kf.obs = obs(kf.nz+1:end,:);
        end
    end
end