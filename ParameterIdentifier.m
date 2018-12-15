classdef ParameterIdentifier < handle
    properties
        residualArguments
        f
        g
        Df_Dxp
        Dg_Dxp
        LowerBound
        UpperBound
        ODEoptions
        OPToptions
    end
    methods
        function p = identify(pi,ResidualArguments,InitialParameters)
            % Residual arguments has the fields: t,y,x0,np
            pi.residualArguments = ResidualArguments;
            pi.residualArguments.m = size(pi.residualArguments.y,1);
            p = lsqnonlin(@pi.Residual,InitialParameters,pi.LowerBound,pi.UpperBound,pi.OPToptions);
        end
        function [residual,jacobian] = Residual(pi,Parameters)
            t = pi.residualArguments.t;
            y = pi.residualArguments.y;
            nx = pi.residualArguments.nx;
            np = pi.residualArguments.np;
            x0 = pi.residualArguments.x0;
            m = pi.residualArguments.m;
            z0 = [x0; zeros(nx*np,1)]; 
            % Implementation of integration
            [t,z] = ode45(@(t,z) pi.ModelAndSensitivity(t,z,Parameters,nx,np),t,z0,pi.ODEoptions);
            x = z(:,1:nx);
            S = reshape(z(:,nx+1:end),m*nx,np); % THIS RESHAPE TODO
            % Residual
            residual = pi.g(x,Parameters)-y;
            % Jacobian 
            jacobian = zeros(m,np);
            [Dg_Dx,Dg_Dp] = pi.Dg_Dxp(x,Parameters);
            for it = 1:m
                jacobian(it,:) = reshape(Dg_Dx*S(1+(it-1)*nx:it*nx,:)+Dg_Dp,1,np);
            end
        end
        function Dz = ModelAndSensitivity(pi,t,z,p,nx,np)
            x = z(1:nx,1);                          % Unpack states
            s = z(nx+1:end,1);                      % Unpack sensitivities
            S = reshape(s,nx,np);                   % Sensitivities as a matrix
            Dx = pi.f(t,x,p);                       % Evaluate the model equations
            [Df_Dx,Df_Dp] = pi.Df_Dxp(t,x,p);       % Evaluate the derivatives
            DSp = Df_Dx*S + Df_Dp;
            Dz = [Dx; DSp(:)];                      % Return derivatives as a vector
        end
        function initialize(pi,Model,PartialDerivatives,Output,OutputPartialDerivatives,Bounds,ODEoptions,OPToptions)
            % input model and derivative model
            pi.f = Model;
            pi.Df_Dxp = PartialDerivatives;
            pi.g = Output;
            pi.Dg_Dxp = OutputPartialDerivatives;
            pi.LowerBound = Bounds.LowerBound;
            pi.UpperBound = Bounds.UpperBound;
            pi.ODEoptions = ODEoptions;
            pi.OPToptions = OPToptions;
        end
    end
end