classdef FourTankSystem < matlab.mixin.Copyable
    properties
        t
        ODEoptions
        record
        m
        gamma
        a
        A
        rho
        h
        g
        qi
        qo
        ms
        hs
        m0
        Ts
        ssc
        ssd
        tfc
        tfd
        Lrvv
        Lrww
    end
    methods
        function timestep(fts,t,inputs)
            % Function help: 
            
            inputs = inputs + fts.Lrww*randn(4,1);
            x = fts.m;
            [t, x] = ode15s(@fts.process,[t(1) t(2)],x,fts.ODEoptions,inputs);
            fts.record.t = [fts.record.t; t(1:end-1,:)];
            fts.record.x = [fts.record.x; x(1:end-1,:)];
            fts.m = x(end,:)';
            fts.record.meas = [fts.record.meas;...
                (fts.m/fts.A/fts.rho+fts.Lrvv*randn(4,1))'];
        end
        function Dx = process(fts,t,x,F)
            % Time
            fts.t = t;
            fts.m = x;
                        
            % Inflows
            fts.qi = [fts.gamma(1)*F(1);...     % Valve 1 to tank 1 [cm3/s]
                fts.gamma(2)*F(2);...           % Valve 2 to tank 2 [cm3/s]
                (1-fts.gamma(2))*F(2)+F(3);...  % Valve 2 to tank 3 [cm3/s]
                (1-fts.gamma(1))*F(1)+F(4)];    % Valve 1 to tank 4 [cm3/s]
            
            % Outflows
            fts.h = fts.m/(fts.rho*fts.A);         % Liquid level in each tank [cm]
            fts.qo = fts.a*sqrt(2*fts.g*fts.h);    % Outflow from each tank [cm3/s]
            
            % Differential equations, mass balances
            Dx = [fts.rho*(fts.qi(1)+fts.qo(3)-fts.qo(1));...   % Tank 1
                fts.rho*(fts.qi(2)+fts.qo(4)-fts.qo(2));...     % Tank 2
                fts.rho*(fts.qi(3)-fts.qo(3));...               % Tank 3
                fts.rho*(fts.qi(4)-fts.qo(4))];                 % Tank 4
        end
        function Dx = steadyStateProcess(fts,x,F)
            Dx = fts.process(0,x,F);
        end
        function initialize(fts,gamma,initials,ODEoptions)
            % Function help:
            
            fts.a = 1.2272;      % [cm2] Area of outlet pipe 
            fts.A = 380.1327;    % [cm2] Cross sectional area of tank
            fts.g = 981;         % [cm/s2] The acceleration of gravity
            fts.rho = 1.00;      % [g/cm3] Density of water
            fts.gamma = gamma;   % [-] Valve constants
            % Initial values
            fts.m0 = initials.m;
            fts.m = fts.m0;
            % Noise initialization to zero
            fts.Lrww = zeros(4);
            fts.Lrvv = zeros(4);
            % Solver options
            fts.ODEoptions = ODEoptions;
            % Records
            Record.t = [];
            Record.x = [];
            Record.meas = [];
            fts.record = Record;
        end
        function addnoise(fts,noise)
            % Function help:
            
            % DO WE NEED RHO FOR THE LINEAR SYSTEM?
            % WHAT IS THE CONNECTION BETWEEN G AND LRWW?
            % TODO
            % Noise
            fts.Lrvv = chol(noise.Rvv,'lower');
            if any(any(noise.Rww))
                fts.Lrww(3:4,3:4) = chol(noise.Rww(3:4,3:4),'lower');
            end
        end
        function reinitialize(fts)
            % Initialization of the equilibrium
            fts.m = fts.ms;
            % Initialization of the records
            fts.record.t = [];
            fts.record.x = [];
            fts.record.meas = [];
        end
        function steadyState(fts,ssInput,initialGuess)
            fts.ms = fsolve(@fts.steadyStateProcess,initialGuess,[],ssInput);
            fts.hs = fts.ms/(fts.A*fts.rho);
            % Initial values refill
            fts.m = fts.m0;
        end
        function linearize(fts)
            %T = (fts.A/fts.a)*sqrt(2*fts.hs/fts.g);
            T = sqrt(fts.ms/(fts.rho*fts.g/fts.A/2))/fts.a;
            System.A = [-1/T(1) 0 1/T(3) 0;0 -1/T(2) 0 1/T(4);0 0 -1/T(3) 0;0 0 0 -1/T(4)];
            System.B = [fts.rho*fts.gamma(1) 0;0 fts.rho*fts.gamma(2);...
                0 fts.rho*(1-fts.gamma(2)); fts.rho*(1-fts.gamma(1)) 0];
            System.G = [0 0; 0 0; fts.rho 0; 0 fts.rho];
            System.C = 1/(fts.rho*fts.A)*eye(4);
            System.Cz = System.C(1:2,:);
            % Noise
            if any(any(fts.Lrvv)) || any(any(fts.Lrww))
                System.Rww = fts.Lrww*fts.Lrww';
                System.Rvv = fts.Lrvv*fts.Lrvv';
                System.Rwv = zeros(4);
            end
            fts.ssc = System;
            % Transfer function
            fts.tfc = tf(ss(System.A,System.B,System.Cz,zeros(2,2)));
        end
        function discretize(fts,Ts)
            System = fts.ssc;
            System.Ts = Ts;
            if any(any(fts.Lrvv)) || any(any(fts.Lrww))
                System.Rww = fts.c2d_noise(System.A,System.G*fts.Lrww(3:4,3:4),System.Ts);
                System.G = eye(4);
            end
            [Ad,Bd] = c2d(System.A,System.B,System.Ts);
            System.A = Ad;
            System.B = Bd;
            fts.ssd = System;
            % Transfer Function
            fts.tfd = tf(ss(System.A,System.B,System.Cz,zeros(2,2),System.Ts));
        end
        function Qd = c2d_noise(fts,A,G,Ts)
            [nx,~] = size(G);
            M = [-A' G*G'; zeros(nx,nx) A];
            Phi = expm(M*Ts);
            Ad = Phi(nx+1:nx+nx,nx+1:nx+nx);
            Qd = Ad'*Phi(1:nx,nx+1:nx+nx);
        end
    end
end