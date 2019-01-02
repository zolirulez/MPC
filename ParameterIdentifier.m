classdef ParameterIdentifier < handle
    properties
        f
        g
        Df_Dxp
        Dg_Dxp
        LowerBound
        UpperBound
        ODEoptions
        optioptions
        args
        nx
        np
    end
    methods
        function p = identify(pi,ResidualArguments,InitialParameters)
            % Residual arguments has the fields: t,y,x0
            pi.args = ResidualArguments;
            pi.args.m = size(pi.args.y,1);
            p = lsqnonlin(@pi.Residual,InitialParameters,pi.LowerBound,pi.UpperBound,pi.optioptions);
        end
        function [residual,jacobian] = Residual(pi,Parameters)
            z0 = [pi.args.x0; zeros(pi.nx*pi.np,1)];
            % Implementation of integration
            [t,z] = ode15s(@(t,z) pi.ModelAndSensitivity(t,z,Parameters),pi.args.t,z0,pi.ODEoptions);
            x = z(:,1:pi.nx);
            try
                S = reshape(z(:,pi.nx+1:end),pi.args.m*pi.nx,pi.np); % THIS RESHAPE TODO
            catch
                disp('')
            end
            % Residual
            residual = pi.g(x,Parameters)-pi.args.y;
            % Jacobian
            jacobian = zeros(pi.args.m,pi.np);
            [Dg_Dx,Dg_Dp] = pi.Dg_Dxp(x,Parameters);
            for it = 1:pi.args.m
                try
                jacobian(it,:) =...
                    reshape(...
                    Dg_Dx*S(it:pi.args.m:end,:)+Dg_Dp,...
                    1,pi.np);
                catch
                    disp('itt')
                end
            end
        end
        function Dz = ModelAndSensitivity(pi,t,z,p)
            x = z(1:pi.nx,1);                          % Unpack states
            s = z(pi.nx+1:end,1);                      % Unpack sensitivities
            S = reshape(s,pi.nx,pi.np);                   % Sensitivities as a matrix
            Dx = pi.f(t,x,p);                       % Evaluate the model equations
            [Df_Dx,Df_Dp] = pi.Df_Dxp(t,x,p);       % Evaluate the derivatives
            DS = Df_Dx*S + Df_Dp;
            Dz = [Dx; DS(:)];                      % Return derivatives as a vector
        end
        function initialize(pi,IdentificationObject,n,Bounds,ODEoptions,optioptions)
            % Function help:
            
            % Dynamics, measurement, and their partial derivatives
            pi.f = @IdentificationObject.f;
            pi.g = @IdentificationObject.g;
            pi.Df_Dxp = @IdentificationObject.Df_Dxp;
            pi.Dg_Dxp = @IdentificationObject.Dg_Dxp;
            pi.nx = n.nx;
            pi.np = n.np;
            pi.LowerBound = Bounds.LowerBound;
            pi.UpperBound = Bounds.UpperBound;
            pi.ODEoptions = ODEoptions;
            pi.optioptions = optioptions;
        end
    end
end