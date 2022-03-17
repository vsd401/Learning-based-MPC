%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -
%------------------------------------------------------------------

classdef NMPC < handle
    %------------------------------------------------------------------
    % Nonlinear MPC class
    % 
    %   solve:
    %
    %   MIN     { SUM_{k=0:N-1} fo(tk,xk,uk,r(t)) } + fend(tN,xN,r(tN))
    %
    %   s.t.    xk+1 = E[f(xk,uk)]
    %           h(xk,uk) == 0
    %           g(xk,uk) <= 0
    %
    %   where the motion model evaluates   [E[xk+1],Var[xk+1]] = f(xk,uk)
    %
    %   for [u0;...;uN-1; e1;...;eN]
    %
    %   where xk: state variables
    %         zk: selected state variables zk=Bd'*xk
    %         r(tk): trajectory   Ŀ��
    %         tk: current time
    %------------------------------------------------------------------
    
    properties
        % Optimizer settings
        maxiter = 200       % max solver iterations
        tol     = 1e-6      % optimizer tolerance
        N       = 30        % prediction horizon
        
        % Define optimization problem
        f       % @fun nonlinear dynamics:  [E[xk+1],Var[xk+1]] = f(xk,uk)
        h       % nonlinear equality constraint function 
        g       % nonlinear inequality constraint function
        B       % ���Բ���ʽ��ʽԼ���Ҳ�
        Beq   % ���Ե�ʽԼ���Ҳ�
        
        fo      % @fun nonlinear cost function�ɱ�����
        fend    % @fend nonlinear cost function for the final state�ն˳ɱ�����
        
        % Optimization dimension
        n   % dim(x) state space dimension ״̬��ά��
        m   % dim(u) input dimension ����ά��
        ne  % dim(e) additional optimization variables dimension
        dt  % time step size
        nh  % number of additional eq. constraints for every time step
        ng  % number of additional ineq. constraints for every time step
        nB  %���Բ���ʽԼ������
        nBeq%���Ե�ʽԼ������
        
        % save last optimal results computed, in order to use as initial guess
        uguess  % <m,N>  initial guess for inputs
        eguess  % <ne,N> initial guess for extra variables
    end
    
    properties(Access=private)
        lb    % lower bound constraints  lb <= vars
        ub    % upper bound constraints        vars <= ub
    end
    
    methods
        
        function obj = NMPC (f, B, Beq, h, g, u_lb, u_ub, n, m, ne, fo, fend, N, dt,bias)
        %------------------------------------------------------------------
        % MPC constructor
        %
        % args:
        %   f: motion model that evaluates  [E[xk+1],Var[xk+1]] = f(xk)
        % 
        % varargin:
        %   provideDynamicsGradient: <bool> when set to true, then the
        %   parameter f must return []
        %------------------------------------------------------------------
           % constraints
           obj.f  = f;
           obj.B = B;
           obj.Beq = Beq;
           obj.h = h;
           obj.g = g;
           % variable dimensions
           obj.n = n;
           obj.m = m;
           obj.ne = ne;
           % get size of additional constraints
           obj.nh = length(h(zeros(n,1),zeros(m,1),zeros(ne,1)));
           obj.ng = length(g(zeros(n,1),bias));
           obj.nB =  length(B);
           obj.nBeq =  length(Beq);
           % cost functions
           obj.fo = fo;
           obj.fend = fend;
           % optimizer parameters
           obj.N = N;  %Ԥ�ⲽ��
           obj.dt = dt; %���㲽��
           
           % set vector of initial guess for optimization
           obj.uguess = zeros(m,N);
           obj.eguess = zeros(ne,N);
           
           % define lower and upper bound constraints
            if ~isempty(u_lb)
                obj.lb = [repmat(u_lb,obj.N,1);    % repeat lower bound for all u0,...,uN-1
                          -Inf(obj.ne*obj.N,1)];             
            else
                obj.lb = [];
            end
            if ~isempty(u_ub)
                obj.ub = [repmat(u_ub,obj.N,1);     % repeat upper bound for all u0,...,uN-1
                          Inf(obj.ne*obj.N,1)];
            else
                obj.ub = [];
            end
        end
        
        
        function numvars = optSize(obj)
        %------------------------------------------------------------------
        % How many variables we need to optimize?
        %
        %   vars_opt = [x0; u0;...;uN-1; e1;...;eN]
        %------------------------------------------------------------------
            numvars = obj.N*obj.m + obj.N*obj.ne;
        end
        
        
        function [u_opt, e_opt,fval,existflag] = optimize(obj, x0, t0, r, UseParallel,bias)
        %------------------------------------------------------------------
        % Calculate first uptimal control input
        %------------------------------------------------------------------
            
            %-------- Set initial guess for optimization variables  -------
            varsguess = [obj.uguess(:); obj.eguess(:)]; % ��ʼֵ
            
            
            %------------------ Optimize  ---------------------------------
            assert(all(size(x0)==[obj.n,1]), 'x0 has wrong dimension!!')
            assert(numel(varsguess) == obj.optSize(), ...
                'There is something wrong with the code. Number of optimization variables does not match!' );
                
            % define cost and constr. functions, as a function only of the
            % optimazation variables. This is a prerequisite for the function fmincon
            costfun = @(vars) obj.costfun(vars, t0, r, x0); % Ŀ�꺯��  vars��Ҫ�Ż��ı���
            nonlcon = @(vars) obj.nonlcon(vars, t0, x0,bias);  %������Լ��
            [Aineq,Bineq,Aeq,Beq]  = obj.lcon( ); %����Լ��
            
            % define optimizer settings
            options = optimoptions('fmincon',...
                                   'Display','iter',...
                                   'Algorithm', 'interior-point',... % 'interior-point',... % 'sqp','interior-point'
                                   'SpecifyConstraintGradient',false,...
                                   'UseParallel',UseParallel,... %'ConstraintTolerance',obj.tol,... Ϊ��ʱ��������ò��з�ʽ�����ݶȡ�ͨ������ΪĬ��ֵfalse�����Խ���
                                   'MaxIterations',obj.maxiter);  %Ѱ������������
            
            % solve optimization problem
            [vars_opt,fval,existflag,~] = fmincon(costfun,varsguess,Aineq,Bineq,Aeq,Beq,obj.lb,obj.ub,nonlcon,options);
            %                      fmincom��Ŀ�꺯������ʼֵ�����Բ���ʽԼ��A��b�����Ե�ʽԼ��Aeq��beq������x���ϡ��½磬 �����Բ���ʽ�������Ż�������
            %         [x,fval]=fmincon(fun,   x0,           A,b,                            Aeq,beq,                         lb,ub ,                   nonlcon,              options)
            
            %------------------ Output results  ---------------------------
            % split variables since vars_opt = [x_opt; u_opt; e_opt]
            [u_opt, e_opt] = splitvariables(obj, vars_opt);
            
            % store current optimization results to use as initial guess for future optimizations
            obj.uguess = u_opt(:,[2:end,end]);
            obj.eguess = e_opt(:,[2:end,end]);
        end
        
    
        function [uvec, evec] = splitvariables(obj, vars)
        %------------------------------------------------------------------
        % args:
        %   vars: <optSize,1> optimization variables
        % out:
        %   uvec: <m,N>
        %   evec: <ne,N>
        %------------------------------------------------------------------
            % split variables
            uvec = vars( (1:obj.N*obj.m) );
            evec = vars( (1:obj.N*obj.ne) + length(uvec) );
            % reshape the column vector <m*N,1> to <m,N>
            uvec = reshape(uvec, obj.m, obj.N);
            % reshape the column vector <ne*N,1> to <ne,N>
            evec = reshape(evec, obj.ne, obj.N);
        end
        
        
        function [mu_xk,var_xk] = predictStateSequence(obj, mu_x0, var_x0, uk)
        %------------------------------------------------------------------
        % Propagate mean and covariance of state sequence, given control
        % input sequence.
        % out:
        %   mu_xk:  <n,N+1>
        %   var_xk: <n,n,N+1>
        %------------------------------------------------------------------
            mu_xk  = zeros(obj.n,obj.N+1);
            var_xk = zeros(obj.n,obj.n,obj.N+1);
            
            mu_xk(:,1)    = mu_x0;
            var_xk(:,:,1) = var_x0;
            
            for iN=1:obj.N      % [x1,...,xN]
                [mu_xk(:,iN+1),var_xk(:,:,iN+1)] = obj.f(mu_xk(:,iN),var_xk(:,:,iN),uk(:,iN));

                % % % % if sum(isnan(mu_xk),'all') || sum(isinf(mu_xk),'all')
                % % % %     error('%s','System dynamics evaluated to NaN or Inf')
                % % % % end

            end
        end

        
        function cost = costfun(obj, vars, t0, r, mu_x0)
        %------------------------------------------------------------------
        % Evaluate cost function for the whole horizon, given variables
        %------------------------------------------------------------------
            % split variables
            [uvec, evec] = obj.splitvariables(vars);
            var_x0 = zeros(obj.n);
            
            % calculate state sequence for given control input sequence and x0
            [mu_xvec,var_xvec] = obj.predictStateSequence(mu_x0, var_x0, uvec);
            
            uncertainty =   zeros(obj.n,obj.N+1);
            for  iN=1:obj.N+1      % i=0:N-1
                % add cost: fo=@(t,mu_x,var_x,u,e,r)
                uncertainty(:,iN) = 1.96 * sqrt(diag(var_xvec(:,:,iN)));  %��ÿ������var_xvec(:,:,i)�ĶԽ�����ȡ������Ϊ��������Ȼ�󿪸���
                
            end
                  figure(555)
                    subplot(4,1,1)
                    NN=linspace(0,obj.N+1,obj.N+1);                    
                    xconf =[NN NN(end:-1:1)];
                    yconf = [mu_xvec(1,:) + uncertainty(1,:)        mu_xvec(1,end:-1:1) - uncertainty(1,end:-1:1) ];
                    fill(xconf,yconf,'r','FaceColor',[1 0.8 0.8],'EdgeColor','none');%FaceColorΪ�����ɫ��EdgeColorΪ�߿���ɫ            
                    hold on
                    plot(NN,mu_xvec(1,:),'*')
                    title('x1')
                    hold off
                    
                    subplot(4,1,2)
                    yconf = [mu_xvec(2,:) + uncertainty(2,:)        mu_xvec(2,end:-1:1) - uncertainty(2,end:-1:1) ];
                    fill(xconf,yconf,'r','FaceColor',[1 0.8 0.8],'EdgeColor','none');%FaceColorΪ�����ɫ��EdgeColorΪ�߿���ɫ            
                    hold on
                    plot(NN,mu_xvec(2,:),'*')
                    title('x2')
                    hold off
                    
                    subplot(4,1,3)
                    yconf = [mu_xvec(3,:) + uncertainty(3,:)        mu_xvec(3,end:-1:1) - uncertainty(3,end:-1:1) ];
                    fill(xconf,yconf,'r','FaceColor',[1 0.8 0.8],'EdgeColor','none');%FaceColorΪ�����ɫ��EdgeColorΪ�߿���ɫ            
                    hold on
                   plot(NN,mu_xvec(3,:),'*')
                    title('x3')
                    hold off
                    
                    subplot(4,1,4)
                    yconf = [mu_xvec(4,:) + uncertainty(4,:)        mu_xvec(4,end:-1:1) - uncertainty(4,end:-1:1) ];
                    fill(xconf,yconf,'r','FaceColor',[1 0.8 0.8],'EdgeColor','none');%FaceColorΪ�����ɫ��EdgeColorΪ�߿���ɫ   
                    hold on
                    plot(NN,mu_xvec(4,:),'*')
                    xlabel('Prediction Step')
                    title('x4')
                    hold off
                    
                    
                    figure(666)
                    subplot(3,1,1)
                    NN=linspace(0,obj.N+1,obj.N+1);                    
                    xconf =[NN NN(end:-1:1)];
                    yconf = [mu_xvec(4,:) + uncertainty(4,:)        mu_xvec(4,end:-1:1) - uncertainty(4,end:-1:1) ];
                    fill(xconf,yconf,'r','FaceColor',[1 0.8 0.8],'EdgeColor','none');%FaceColorΪ�����ɫ��EdgeColorΪ�߿���ɫ            
                    hold on
                    plot(NN,mu_xvec(4,:),'*')
                    plot(NN,mu_xvec(4,:))
                    ylabel('vx/ m/s')
                    hold off
                    
                    subplot(3,1,2)
                    yconf = [mu_xvec(5,:) + uncertainty(5,:)        mu_xvec(5,end:-1:1) - uncertainty(5,end:-1:1) ];
                    fill(xconf,yconf,'r','FaceColor',[1 0.8 0.8],'EdgeColor','none');%FaceColorΪ�����ɫ��EdgeColorΪ�߿���ɫ            
                    hold on
                    plot(NN,mu_xvec(5,:),'*')
                    plot(NN,mu_xvec(5,:))
                    ylabel('vy/ m/s')
                    hold off
                    
                    subplot(3,1,3)
                    yconf = [mu_xvec(6,:) + uncertainty(6,:)        mu_xvec(6,end:-1:1) - uncertainty(6,end:-1:1) ];
                    fill(xconf,yconf,'r','FaceColor',[1 0.8 0.8],'EdgeColor','none');%FaceColorΪ�����ɫ��EdgeColorΪ�߿���ɫ            
                    hold on
                   plot(NN,mu_xvec(6,:),'*')
                    plot(NN,mu_xvec(6,:))
                    ylabel('yaw rate/ rad/s')
%                     axis([0,16,-2,2])
                    xlabel('Prediction horizon')
                    hold off
                    
                     figure(777)
                    subplot(1,1,1)
                    yconf = [mu_xvec(7,:) + uncertainty(7,:)        mu_xvec(7,end:-1:1) - uncertainty(7,end:-1:1) ];
                    fill(xconf,yconf,'r','FaceColor',[1 0.8 0.8],'EdgeColor','none');%FaceColorΪ�����ɫ��EdgeColorΪ�߿���ɫ            
                    hold on
                    plot(NN,mu_xvec(7,:),'*')
                    title('vx/ m/s')
                    hold off
            cost = 0;
            t = t0;
            for iN=1:obj.N      % i=0:N-1
                % add cost: fo=@(t,mu_x,var_x,u,e,r)
                cost = cost + obj.fo(t, mu_xvec(:,iN), var_xvec(:,:,iN), uvec(:,iN), evec(:,iN), r);
                
                % % % % if sum(isnan(cost),'all') || sum(isinf(cost),'all')
                % % % %     error('Cost function evaluated to NaN or Inf')
                % % % % end
                
                % update current time
                t = t + iN * obj.dt;
            end
            % final cost: fend=@(t,mu_x,var_x,e,r)
            cost = cost + obj.fend(t, mu_xvec(:,end), var_xvec(:,:,end), evec(:,iN), r);
            
            % normalize cost by horizon size
            cost = cost / (obj.N+1);
        end
        

        function [cineq,ceq] = nonlcon(obj, vars, t0, mu_x0,bias) %������Լ��
        % function [cineq,ceq,gradvars_cineq,gradvars_ceq] = nonlcon(obj, vars, t0, x0)
        %------------------------------------------------------------------
        % Evaluate nonlinear equality and inequality constraints
        % out:
        %   cineq = g(x,u) <= 0 : inequality constraint function
        %   ceq   = h(x,u) == 0 : equality constraint function
        %   gradx_cineq(x,u): gradient of g(x,u) w.r.t. x
        %   gradx_ceq(x,u):   gradient of h(x,u) w.r.t. x
        %------------------------------------------------------------------             
            % init outputs
            ceq = []; 
            cineq = [];
        
            % if there are no constraints, then there is nothing to do here
            if obj.nh==0 && obj.ng==0
                return
            end
            
            % init vectors to speedup calculations
            ceq_h   = zeros(obj.nh, obj.N);
            cineq_g = zeros(obj.ng, obj.N);
            
            % vars_size = obj.optSize();
            % gradvars_cineq = zeros(vars_size,obj.n);
        
            % split variables
            [uvec, evec] = obj.splitvariables(vars);
            var_x0 = zeros(obj.n);
            
            % calculate state sequence for given control input sequence and x0
            [mu_xvec,var_xvec] = obj.predictStateSequence(mu_x0, var_x0, uvec);
            

            t = t0;
            for iN=1:obj.N
                % append provided equality constraints(h==0)
                ceq_h(:,iN) = obj.h(mu_xvec(:,iN),uvec(:,iN), evec(:,iN));
                % provided inequality constraints (g<=0)
                cineq_g(:,iN) = obj.g(mu_xvec(:,iN),bias);
                t = t + iN * obj.dt;
            end

            ceq   = ceq_h(:);
            cineq = cineq_g(:);
        end
   
        
       function [Aineq,Bineq,Aeq,Beq] = lcon(obj) %����Լ��,Լ������������
        % function [Aineq,Aeq,Bineq,Beq] = lcon(obj, vars, t0, x0)
        %------------------------------------------------------------------
        % Evaluate linear equality and inequality constraints
        % out:
        %   Aineq * x<= 0 : inequality constraint function
        %   Aeq *x  = Beq : equality constraint function
        %------------------------------------------------------------------             
            % init outputs
            Aeq = []; 
            Aineq = [];
            Beq = []; 
            Bineq = [];
        
            % if there are no constraints, then there is nothing to do here
            if obj.nB==0 && obj.Beq==0
                return
            end
            Ain =  [-diag(ones(1,obj.nB*obj.N), 0)+ diag(ones(1,obj.nB*(obj.N-1)), obj.nB); diag(ones(1,obj.nB*obj.N), 0)- diag(ones(1,obj.nB*(obj.N-1)), obj.nB)];
            % init vectors to speedup calculations
           Aineq   =[Ain(1:obj.nB*(obj.N-1),:); Ain(obj.nB*obj.N+1:end-obj.nB,:)];
           Aeq = zeros(obj.nBeq*obj.N, obj.nBeq*obj.N);
           Bineq = [repmat(obj.B, (obj.N-1)*2,1)];
           Beq = zeros(obj.nBeq*obj.N, 1);
            
            

        end
        
    end
end


