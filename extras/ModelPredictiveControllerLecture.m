classdef ModelPredictiveControllerLecture < handle
    % Model Predictive Controller. Handles constrains on input, input
    %   rate, output (for the latter soft constrains as well). The current
    %   version is not prepared for a time varying system.
    properties
        mpcType
        % Markov matrix
        Gamma 
        % Auxuliary matrices
        I0
        Lambda
        % Weights
        Wz
        Wu
        WDu
        % Intermediate gradients
        Mz
        Mu
        MDu
        % Intermediate gradients
        gs
        gt
        % Final quadratic properties
        H
        g
        rho
        % Constraints
        A
        bl
        bu
        l
        u
        % DC control signal
        Uc
        % Memory
        U_1
        % Length
        nx   % Number of states
        nz   % Number of outputs
        nu   % Number of inputs
        j    % Length of horizon
    end
    methods
        function controlPrepare(mpc,b,R,refBounds,inputBounds,u_1)
            Umin = inputBounds.Umin;
            Umax = inputBounds.Umax;
            DUmin = inputBounds.DUmin;
            DUmax = inputBounds.DUmax;
            Rmin = refBounds.Rmin;
            Rmax = refBounds.Rmax;
            % Inputs from predictor
            c = R-b;
            % Intermediate gradients
            g1 = mpc.Mz*c + mpc.Mu*mpc.Uc + mpc.MDu*u_1;
            mpc.g = [g1; mpc.gs; mpc.gt];
            % Intermediate constants
%             rhoz = 0.5*(mpc.Wz*c)'*(mpc.Wz*c);
%             rhou = 0.5*(mpc.Wu*mpc.Uc)'*(mpc.Wu*mpc.Uc);
%             rhoDu = 0.5*(mpc.WDu*mpc.I0*u_1)'*(mpc.WDu*mpc.I0*u_1);
%             mpc.rho = rhoz + rhou + rhoDu;
            % Constraints
            mpc.l = [Umin; zeros(mpc.nz*mpc.j,1); zeros(mpc.nz*mpc.j,1)];
            mpc.u = [Umax; Inf*ones(mpc.nz*mpc.j,1); 1e10*ones(mpc.nz*mpc.j,1)];
            mpc.bl = [DUmin+mpc.I0*u_1; Rmin-b; -1e10*ones(mpc.nz*mpc.j,1)];
            mpc.bu = [DUmax+mpc.I0*u_1; 1e10*ones(mpc.nu*mpc.j,1); Rmax-b];
        end
        function U = controlCompute(mpc,b,R,refBounds,inputBounds,u_1,optioptions)
            mpc.controlPrepare(b,R,refBounds,inputBounds,u_1)
            U = quadprog(mpc.H,mpc.g',[mpc.A; -mpc.A],[mpc.bu -mpc.bl],...
                [],[],mpc.l,mpc.u,mpc.U_1,optioptions);
            mpc.U_1 = U;
        end
        function initialize(mpc,W,Gamma,Uc,horizon,n)
            mpc.j = horizon;
            mpc.nu = n.nu;
            mpc.nx = n.nx;
            mpc.nz = n.nz;
            % DC control signal
            mpc.Uc = Uc;
            mpc.U_1 = Uc;
            % Markov matrix
            mpc.Gamma = Gamma;
            % Auxuliary matrices
            mpc.I0 = [eye(mpc.nu); zeros(mpc.nu*(mpc.j-1),mpc.nu)];
            IN = eye(mpc.j);
            mpc.Lambda = kron(IN,eye(mpc.nu))-...
                [zeros(mpc.nu,mpc.nu*mpc.j);...
                [kron(eye(mpc.j-1),eye(mpc.nu))...
                zeros(mpc.nu*(horizon-1),mpc.nu)]];
            % Weights
            mpc.Wz = kron(IN,W.Wz);
            mpc.Wu = kron(IN,W.Wu);
            mpc.WDu = kron(IN,W.WDu);
            Ws1 = kron(IN,W.Ws1);
            Ws2 = kron(IN,W.Ws2);
            Wt1 = kron(IN,W.Wt1);
            Wt2 = kron(IN,W.Wt2);
            % Intermediate Hessians
            Hz = (mpc.Wz*mpc.Gamma)'*(mpc.Wz*mpc.Gamma);
            Hu = mpc.Wu'*mpc.Wu;
            HDu = (mpc.WDu*mpc.Lambda)'*(mpc.WDu*mpc.Lambda);
            % Intermediate gradients
            mpc.Mz = -(mpc.Wz*mpc.Gamma)'*mpc.Wz;
            mpc.Mu = -Hu;
            mpc.MDu = -(mpc.WDu*mpc.Lambda)'*mpc.WDu*mpc.I0;
            % Intermediate Hessians
            H1 = Hz + Hu + HDu;
            Hs = Ws2'*Ws2;
            Ht = Wt2'*Wt2;
            % Intermediate gradients
            mpc.gs = Ws1*ones(mpc.j*mpc.nu,1); 
            mpc.gt = Wt1*ones(mpc.j*mpc.nu,1);
            % Final Hessian
            mpc.H = blkdiag(H1,Hs,Ht);
            % Constraints
            mpc.A = [mpc.Lambda zeros(mpc.nu*mpc.j) zeros(mpc.nu*mpc.j);...
                mpc.Gamma eye(mpc.nu*mpc.j) zeros(mpc.nu*mpc.j);...
                mpc.Gamma zeros(mpc.nu*mpc.j) -eye(mpc.nu*mpc.j)];
        end
    end
end