classdef IdentificationObject < handle
    % Object for identifying continuous transfer functions with no zeros
    properties
        np
        nx
    end
    methods
        function Dx = f(idobj,t,x,p)
            Dx = [zeros(idobj.nx-1,1) eye(idobj.nx-1); -flip(p(2:end)')]*x...
                + [zeros(idobj.nx-1,1); p(1)];
        end
        function y = g(idobj,x,p)
            y = x(:,1);
        end
        function [dfdx,dfdp] = Df_Dxp(idobj,t,x,p)
            dfdx = [zeros(idobj.nx-1,1) eye(idobj.nx-1); -flip(p(2:end)')];
            dfdp = [zeros(idobj.nx-1,idobj.np); 1 -flip(x')];
        end
        function [dgdx,dgdp] = Dg_Dxp(idobj,x,p)
            dgdx = [1 zeros(1,idobj.nx-1)];
            dgdp = zeros(1,idobj.np);
        end
        function initialize(idobj,n)
            % Function help:

            idobj.np = n.np;
            idobj.nx = n.nx;
        end
    end
end