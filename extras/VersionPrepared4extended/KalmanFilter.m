classdef KalmanFilter < handle
    % This Kalman filter is an object for solving general computations
    %   required by a Kalman Filter. Note that the object does not work
    %   like a single function, therefore for different applications
    %   different instances are needed (e.g. different predictor and
    %   filter.
    properties
        % Types: stationary, timevarying, timeinvariant
        kalmanFilterType
        % Domains: discrete, continuous-discrete
        kalmanFilterDomain
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
        Ac   % Continuous state matrix
        Gc   % Continuous noise input matrix
        Bc   % Continuous input matrix
        % Lengths
        nx   % Number of states
        % Independent variable
        t    % Time
        Ts   % Sampling time
    end
    methods
        function DxP = process(kf,t,xP,u)
            x = xP(1:kf.nx,:);
            P = xP(kf.nx+1:end,:);
            % State evolution
            Dx = kf.Ac*x + kf.Bc*u;
            % Covariance evolution
            DP = kf.Ac*P + P*kf.Ac' + kf.Gc*kf.Gc';
            % Derivative output
            DxP = [Dx; DP];
        end
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
        function timeUpdate(kf,u,t)
            if strcmp(kf.kalmanFilterDomain,'discrete')
                % State
                kf.x1 = kf.A*kf.xf + kf.B*u(:,1) + kf.G*kf.wf;
                % Covariance
                kf.P1 = kf.A*kf.Pf*kf.A' + kf.G*kf.Q*kf.G'...
                    - kf.A*kf.Kx*kf.S'*kf.G' - kf.G*kf.S*kf.Kx'*kf.A';
            else
                xP = [hx.xf; hx.Pf];
                [t, xP] = ode15s(kf.process,[t t+kf.Ts],xP,kf.ODEoptions,u(:,1));
                % State
                kf.x1 = xP(end,1:kf.nx);
                % Covariance
                kf.P1 = xP(end,kf.nx+1:end);
            end
        end
        function statePrediction(kf,u,Q,j,t)
            if strcmp(kf.kalmanFilterDomain,'discrete')
                % State
                kf.xj(:,1) = kf.A*kf.xf + kf.B*u(:,1);
                % Covariance
                kf.Pj{1} = kf.A*kf.Pf*kf.A' + kf.G*kf.Q*kf.G';
                for it = 2:j
                    % State
                    kf.xj(:,it) = kf.A*kf.xj(:,it-1) + kf.B*u(:,it);
                    % Covariance
                    kf.Pj{it} = kf.A*kf.Pj{it-1}*kf.A' + kf.G*Q{it}*kf.G';
                end
            else
                xP = [hx.xf; hx.Pf];
                [t, xP] = ode15s(kf.process,[t t+kf.Ts],xP,kf.ODEoptions,u(:,1));
                xPj(:,1) = xP(end,:)';
                for it = 2:j
                    [t, xP] = ode15s(kf.process,[t+(it-1)*kf.Ts t+it*kf.Ts],...
                        xP,kf.ODEoptions,u(:,it));
                    xPj(:,it) = xP(end,:)';
                end
                % Todo: give values to the states of the kalman filter.
                %   Make it for extended, because this way it has no sense.
                %   Also, resizing is NOT considered yet for the covariance
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
        function [xf, x1, xj, z] = outputPredictor(kf,u,y,Q,j)
            kf.measurementUpdate(y);
            kf.timeUpdate(u);
            if strcmp(kf.kalmanFilterType,'timevariant')
                kf.timeVariation(kf,t);
            end
            kf.statePrediction(u,Q,j);
            kf.outputPrediction(j);
            xf = kf.xf;
            x1 = kf.x1;
            xj = kf.xj;
            z = kf.z;
        end
        function timeVariation(kf,t) % not ready!
            % Function for the case of changing system matrices. The
            %     methodology is not clarified yet.
            
            % Furthermore, continuous version is not implemented yet
            kf.A = 0;
            kf.B = 0;
            kf.C = 0;
            kf.Cz = 0;
            kf.R = 0;
            kf.Q = 0;
            kf.S = 0;
        end
        function initialize(kf,system,noise,initial,kalmanFilterSpecs)
            kf.kalmanFilterType = kalmanFilterSpecs.kalmanFilterType;
            kf.kalmanFilterDomain = kalmanFilterSpecs.kalmanFilterDomain;
            if strcmp(kf.kalmanFilterDomain,'discrete')
                kf.A = system.A;
                kf.B = system.B;
                kf.G = system.G;
            else
                kf.Ac = system.Ac;
                kf.Bc = system.Bc;
                kf.Gc = system.Gc;
                kf.ODEoptions = kalmanFilterSpecs.ODEoptions;
            end
            kf.C = system.C;
            kf.Cz = system.Cz;
            kf.R = noise.R;
            kf.Q = noise.Q;
            kf.S = noise.S;
            kf.x1 = initial.x;
            kf.P1 = initial.P;
            kf.nx = length(kf.x1);
            if strcmp(kf.kalmanFilterType,'stationary')
                kf.stationaryInitialization();
            end
        end
        function stationaryInitialization(kf)
            kf.P1s = dare(kf.A',kf.C',kf.G*kf.Q*kf.G',kf.R,kf.G*kf.S);
            Res = kf.C*kf.P1s*kf.C' + kf.R;
            kf.Kx = kf.P1s*kf.C'/Res;
            kf.Kw = kf.S/Res;
        end
    end
end