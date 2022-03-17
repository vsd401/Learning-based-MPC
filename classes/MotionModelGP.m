%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -
%------------------------------------------------------------------

classdef (Abstract) MotionModelGP < handle
%--------------------------------------------------------------------------
%   Abstract class for implementing Motion Model of Plant
%   Intended to be used with GP and NMPC classes
%
%   Please inherit this class and implement all the (Abstract) methods and
%   variables
%
%   xk+1 = fd(xk,uk) + Bd * ( d(zk) + w ),    
%
%       where: zk = [Bz_x*xk ; Bz_u*uk],
%              d ~ N(mean_d(zk),var_d(zk))
%              w ~ N(0,var_w)
%   
%--------------------------------------------------------------------------

    properties (Abstract, Constant)
        Bd    % <n,xx>  xk+1 = fd(xk,uk) + Bd*d(zk)
        Bz_x  % <yy,n>  z = [Bz_x*x;Bz_u*u]  
        Bz_u  % <ww,n>
        n     % <1>     number of outputs x(t)
        m     % <1>     number of inputs u(t)
        nd    % <1>     output dimension of d(z)
        nz    % <1>     dimension of z(t)
    end
    
    properties (SetAccess=private)
        discretization = 'ode2';
        
        d      % [E[d(z)] , Var[d(z)]] = d(z): disturbace model
        var_w  % measurement noise covariance matrix. w ~ N(0,var_w)
    end
    
    properties (SetAccess=public)
        % is_d_active = false
    end
    
    methods (Abstract)
        xdot = f (obj, x, u) %系统函数 
        %------------------------------------------------------------------
        %   Continuous time dynamics.
        %   out:
        %       xdot: <n,1> time derivative of x given x and u
        %------------------------------------------------------------------
        
        gradx = gradx_f(obj, x, u) %系统方程x_dot=f(x) 对x的偏导，也就是A矩阵
        %------------------------------------------------------------------
        %   Continuous time dynamics.
        %   out:
        %       gradx: <n,n> gradient of xdot w.r.t. x
        %------------------------------------------------------------------
        
        gradu = gradu_f(obj, x, u) %系统方程x_dot=f(x) 对u的偏导，也就是B矩阵
        %------------------------------------------------------------------
        %   Continuous time dynamics.
        %   out:
        %       gradu: <m,n> gradient of xdot w.r.t. u
        %------------------------------------------------------------------
    end
    
    methods
        function obj = MotionModelGP (d, var_w)
        %------------------------------------------------------------------
        %   object constructor
        %   args:
        %       d: evaluates nonlinear motion model mean and covariance 
        %          function [mean_d, var_d] = d(z),   with z = Bz*x
        %       var_w: <1> measurement noise covariance
        %   obs: if d or var_w are empty, then this function will set them
        %        to zero with the correct dimensions
        %------------------------------------------------------------------
            obj.d = d;
            if isempty(obj.d)  %如果是空的，那么它就是维度为nd的0矩阵
                mu_d  = @(z) zeros(obj.nd,1);
                var_d = @(z) zeros(obj.nd);
                obj.d = @(z) deal( mu_d(z), var_d(z) );
            end
            
            obj.var_w = var_w;
            if isempty(obj.var_w)
                obj.var_w = zeros(obj.nd);
            end
            
            %--------------------------------------------------------------
            % assert model
            %--------------------------------------------------------------
            assert(size(obj.Bz_x,1) + size(obj.Bz_u,1) == obj.nz, ...
                sprintf('obj.Bz_x and obj.Bz_u matrices should have %d columns in total, but have %d',obj.nz,size(obj.Bz_x,1) + size(obj.Bz_u,1)))
            assert(size(obj.Bd,2) == obj.nd, ...
                sprintf('obj.Bd matrix should have %d columns, but has %d',obj.nd,size(obj.Bd,2)))
            
            assert( all(size(obj.var_w)==[obj.nd,obj.nd]), ...
                sprintf('Variable var_w should have dimension %d, but has %d',obj.nd,size(obj.var_w,1)))
            assert(size(obj.Bd,1) == obj.n, ...
                sprintf('obj.Bd matrix should have %d rows, but has %d',obj.n,size(obj.Bd,1)))
            assert(size(obj.Bz_x,2) == obj.n || isempty(obj.Bz_x), ...
                sprintf('obj.Bz_x matrix should have %d columns, but has %d',obj.n,size(obj.Bz_x,1)))
            assert(size(obj.Bz_u,2) == obj.m || isempty(obj.Bz_u), ...
                sprintf('obj.Bz_u matrix should have %d columns, but has %d',obj.m,size(obj.Bz_u,1)))
            
            % validate given disturbance model
            ztest = [obj.Bz_x*zeros(obj.n,1) ; obj.Bz_u*zeros(obj.m,1)];
            [muy,vary] = obj.d(ztest);
            assert( size(muy,1)==obj.nd, ...
                sprintf('Disturbance model d evaluates to a mean value with wrong dimension. Got %d, expected %d',size(muy,1),obj.nd))
            assert( all(size(vary)==[obj.nd,obj.nd]), ...
                sprintf('Disturbance model d evaluates to a variance value with wrong dimension. Got %d, expected %d',size(vary,1),obj.nd))
        end
        
        
        function zk = z(obj, xk, uk)
        %------------------------------------------------------------------
        % select variables (xk,uk) -> z
        %------------------------------------------------------------------
            if ~isempty(obj.Bz_x)
                z_xk = obj.Bz_x * xk;  else, z_xk=[];
            end
            if ~isempty(obj.Bz_u)
                z_uk = obj.Bz_u * uk; else, z_uk = [];
            end
            zk = [ z_xk ; z_uk ];
        end
        
        
        function [xkp1, gradx_xkp1] = fd (obj, xk, uk, dt) % fd为下一步状态预测；gradx_fd为 I+A*det_t 
        %------------------------------------------------------------------
        %   Discrete time dynamics (ODE1 Euler approximation)
        %   args:
        %       xkp1: <n,1> state prediction (without disturbance model)
        %       grad_xkp1: <n,n> gradient of xkp1 w.r.t. xk
        %------------------------------------------------------------------            
            if strcmp(obj.discretization,'ode1')
                %-----------------
                % Fixed step ODE1
                %-----------------
                % calculate continous time dynamics
                xkp1 = xk + dt * obj.f(xk,uk);
                
            elseif strcmp(obj.discretization,'ode2')
                %-----------------
                % Fixed step ODE2 (developed by myself)
                %-----------------
                [~,xkp1] = ODE.ode2(@(t,x) obj.f(x,uk), xk, dt, dt);
                
            elseif strcmp(obj.discretization,'ode4')
                %-----------------
                % Fixed step ODE4 (developed by myself)
                %-----------------
                [~,xkp1] = ODE.ode4(@(t,x) obj.f(x,uk), xk, dt, dt);
                
            elseif strcmp(obj.discretization,'ode23')
                %-----------------
                % Variable step ODE23
                %-----------------
                opts_1 = odeset('Stats','off','RelTol',1e-1,'AbsTol',1e-1);
                [~,xkp1_23] = ode23( @(t,x) obj.f(x,uk), [0 dt], xk, opts_1);
                xkp1 = xkp1_23(end,:)';
            
            elseif strcmp(obj.discretization,'ode23')
                %-----------------
                % Variable step ODE23
                %-----------------
                opts_1 = odeset('Stats','off','RelTol',1e-1,'AbsTol',1e-1);
                [~,xkp1_23] = ode23( @(t,x) obj.f(x,uk), [0 dt], xk, opts_1);
                xkp1 = xkp1_23(end,:)';
                
            else
                error('Chosen discretization: %s is not yet implemented',obj.discretization);
            end
            
            % for now, gradient is being discretized using a simple ode1
            gradx_xdot = obj.gradx_f(xk,uk);
            gradx_xkp1 = eye(obj.n) + dt * gradx_xdot;
        end
        
        
        function [mu_xkp1,var_xkp1] = xkp1 (obj, mu_xk, var_xk, uk, dt)
        %------------------------------------------------------------------
        %   State prediction (motion model) using Extended Kalman Filter 
        %   approach
        %
        %       xk+1 = fd(xk,uk) + Bd * ( d(zk) + w ),  zk=Bz*xk
        %
        %------------------------------------------------------------------
            % calculate discrete time dynamics     预估预估预估预估预估预估预估预估
            [fd, gradx_fd] = obj.fd(mu_xk,uk,dt);  % fd为下一步状态预测；gradx_fd为 I+A*det_t 
            
            % calculate grad_{x,d,w} xk+1
            grad_xkp1 = [gradx_fd; obj.Bd'; obj.Bd']; % 对fest求偏导
            
            % select variables (xk,uk) -> z  z=[zx,zu] zx=Bzx*x zu=Bzu*u
            z = obj.z(mu_xk,uk);
           
            % evaluate disturbance
            [mu_d, var_d] = obj.d(z);  %当数据长度符合要求的时候，就用GP模型来获取高斯分布的均值和协方差，数据长度太短的时候，认为均值和协方差都为0
                
            % A) Mean Equivalent Approximation:
            var_x_d_w = blkdiag(var_xk, var_d, obj.var_w);
            
            % B) Taylor Approximation
            %--------------------------------------------------------------
            %  TODO
            %--------------------------------------------------------------
            
            % predict mean and variance (Extended Kalman Filter) 矫正矫正矫正矫正矫正矫正
            mu_xkp1  = fd  + obj.Bd * ( mu_d );
            var_xkp1 = grad_xkp1' * var_x_d_w * grad_xkp1; % zeros(obj.n);
        end
        
    end
end
