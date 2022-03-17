%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -

%   Generate Animation for the main_singletrack.m script
%------------------------------------------------------------------


classdef SingleTrackAnimation < handle
    
    properties
        % object that contains track coordinates
        racetrack @ RaceTrack
        
        % variables that contain history of vehicle states and inputs to be ploted
        mu_x_pred_opt
        var_x_pred_opt
        u_pred_opt
        x_ref
        
        %variables that obstacles state info
       XY_o
       PSIt_o
       bias
       v_obs

        % Track animation handles
        h_fig
        h_ltrack
        h_rtrack
        h_track_center
        h_mu_x_pred_opt
        h_var_x_pred_opt
        h_x_ref
        h_x_trace
        h_car
        h_obs_cars
        h_obs_ellipse
        h_x_trace_obs
        
        
        % Scope handles
        h_scopex
        h_scopeu
        
        % covariance ellipse properties
        ell_npoints = 30     % number of points to make an ellipse
        ell_level = 2        % plot ell_level*sigma ellipse curves
        
        k   % current time step
        N   % horizon length
    end
    
    methods
        function obj = SingleTrackAnimation(racetrack, mu_x_pred_opt, var_x_pred_opt, u_pred_opt, x_ref, XY_o, PSIt_o, bias,v_obs)
            obj.racetrack       = racetrack;
            obj.mu_x_pred_opt   = mu_x_pred_opt;
            obj.var_x_pred_opt  = var_x_pred_opt;
            obj.u_pred_opt      = u_pred_opt;
            obj.x_ref                  = x_ref;
            
                    obj.XY_o          = XY_o;
                    obj.PSIt_o         =  PSIt_o;
                    obj.bias          = bias;
                    obj.v_obs       =v_obs;
            
            
            % get horizon length from inputs  
            obj.N = size(obj.mu_x_pred_opt,2);  
        end
        
        function initTrackAnimation(obj)
        % -----------------------------------------------------------------
        %   Init Animation of the vehicle driving the track. Please call
        %   updateTrackAnimation(k) to move forward with the animation.
        % -----------------------------------------------------------------
            obj.h_fig = figure('Color','w','Position',[68 128 1272 333]);
%             title('Adaptive Gaussian-Process MPC')
            hold on;
%             grid on;
%             axis equal;
%             axis([0+obj.mu_x_pred_opt(7,1,k),30+obj.mu_x_pred_opt(7,1,k),-5,5]);

            
            % -------------------------------------------------------------
            %   plot track asphalt
            % -------------------------------------------------------------
            n = length(obj.racetrack.track_l);
            v = [obj.racetrack.track_l(:,1:n)'; obj.racetrack.track_r(:,1:n)'];   % <n,2>
            f = [1:n-1 ; 2:n; n+2:2*n; n+1:2*n-1]';
            patch('Faces',f,'Vertices',v,'FaceColor',[0.5 0.5 0.5],'LineStyle', 'none')  %face的每一行定义了一个多边形，有多少行就定义多少个多边形；其中face的每一行是索引号，去Vetrices中索引顶点坐标
            %%例如
%          v=[  0.5000    3.0000          f=[1     2     7     6
%                 1.0000    3.0000                2     3     8     7
%                 1.5000    3.0000                3     4     9     8
%                 2.0000    3.0000                4     5    10     9]                                      f
%                 2.5000    3.0000                
%                 0.5000   -3.0000
%                 1.0000   -3.0000
%                 1.5000   -3.0000
%                 2.0000   -3.0000
%                 2.5000   -3.0000]
%  f中的第1行 表示 v中的第2行、第7行、第6行形成一个多边形的四个顶点坐标，
            
            
            %%
            % -------------------------------------------------------------
            %   plot track borders
            % -------------------------------------------------------------
            obj.h_ltrack = plot(obj.racetrack.track_l(1,:),obj.racetrack.track_l(2,:),'k','LineWidth',5);
            obj.h_rtrack = plot(obj.racetrack.track_r(1,:),obj.racetrack.track_r(2,:),'k','LineWidth',5);
            
            %% plot track center
            % -------------------------------------------------------------
%             obj.h_track_center = plot(obj.racetrack.track_c(1,:),obj.racetrack.track_c(2,:),'--','LineWidth',0.5);
             %% plot track marker
            % -------------------------------------------------------------
            obj.h_track_center = plot(obj.racetrack.track_l(1,:),obj.racetrack.track_l(2,:)./3,'w.--','LineWidth',3);
            obj.h_track_center = plot(obj.racetrack.track_r(1,:),obj.racetrack.track_r(2,:)./3,'w.--','LineWidth',3);
            % -------------------------------------------------------------
            %   setup state predictions
            % -------------------------------------------------------------
            k = 1;
            
            % trace vehicle path
            obj.h_x_trace = patch(obj.mu_x_pred_opt(1,1,k),obj.mu_x_pred_opt(2,1,k),obj.mu_x_pred_opt(3,1,k),...
                                  'EdgeColor','interp',...
                                  'Marker','*',...
                                  'MarkerSize',2,...
                                  'LineWidth',3,...
                                  'EdgeAlpha',0.5); 
            
            % plot car
            obj.h_car = patch('Faces',1:4,'Vertices',NaN(4,2),...
                        'EdgeColor','black',...
                        'FaceColor','none',...
                        'LineWidth',1);
             % plot obstacle cars
             for i = 1: size(obj.bias,2)*2
                 obj.h_x_trace_obs{i} = patch(obj.XY_o(1,1,k),obj.XY_o(2,1,k),obj.v_obs(1,1),...
                                      'EdgeColor','interp',...
                                      'Marker','none',...
                                      'LineWidth',1.5,...
                                      'EdgeAlpha',0.5); 
                   obj.h_obs_cars{i} = patch('Faces',1:4,'Vertices',NaN(4,2),...
                            'EdgeColor','black',...
                            'FaceColor','black',...
                            'LineWidth',1); 
                        
                    obj.h_obs_ellipse{i} = plot(NaN,NaN,...
                        '-g','LineWidth',0.75, 'Marker','none');
             end
            % reference trajectory for prdiction
            obj.h_x_ref = plot(NaN,NaN,...
                        '-','LineWidth',1.0, 'Marker','o',...
                        'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k', 'Color','k',...
                        'MarkerSize',2,...
                        'DisplayName','Reference trajectory',...
                        'XDataSource', 'obj.x_ref(1,:,obj.k)',...
                        'YDataSource', 'obj.x_ref(2,:,obj.k)');

            % optimal prediction means
            obj.h_mu_x_pred_opt = patch(obj.mu_x_pred_opt(1,:,k),obj.mu_x_pred_opt(2,:,k),obj.mu_x_pred_opt(3,:,k),...
                        'EdgeColor','interp','Marker','none',...
                        'LineWidth',1,...
                        'DisplayName','replanned optimial trajectory' );
                    
            % plot prediction covariance ellipses
            ell = NaN(obj.ell_level, obj.ell_npoints);
            for i=1:obj.N
                obj.h_var_x_pred_opt{i} = patch('Faces',1:obj.ell_npoints,'Vertices',ell','FaceColor',[1 0 0],'FaceAlpha',0.1,'LineStyle', 'none');
            end      
                             
                              
            % display legend, colorbar, etc.
            leg = legend([obj.h_mu_x_pred_opt,obj.h_x_ref],'Location','northeast');
            c = colorbar;
            c.Label.String = 'Vehicle predicted velocity [m/s]';
            caxis([1 30])
            
            % lock up axis limits
            xlim('manual')
            ylim('manual')
            drawnow
        end
        
        function initScope(obj)
        % -----------------------------------------------------------------
        %   Init Scope that shows vehicle prediction signals
        % -----------------------------------------------------------------
            obj.h_scopex = figure('Position',[-1006 86 957 808]);
            names = {'I-x','I-y','psi','V-vx','V-vy','psidot','track_dist'};
            angles = [0 0 1 0 0 1 0];
            for i=1:numel(names)
                subplot(4,2,i);
                p = plot(NaN);
                if angles(i)
                    p.YDataSource = sprintf('rad2deg(obj.mu_x_pred_opt(%d,:,obj.k))',i);
                else
                    p.YDataSource = sprintf('obj.mu_x_pred_opt(%d,:,obj.k)',i);
                end
                xlabel('Prediction horizon')
                grid on;
                title(names{i});
            end
            
            obj.h_scopeu = figure('Position',[-1879 93 867 795]);
            names = {'delta','T','Track vel'};
            angles = [1 0 0];
            for i=1:numel(names)
                subplot(2,2,i);
                p = plot(NaN);
                if angles(i)
                    p.YDataSource = sprintf('rad2deg(obj.u_pred_opt(%d,:,obj.k))',i);
                else
                    p.YDataSource = sprintf('obj.u_pred_opt(%d,:,obj.k)',i);
                end
                xlabel('Prediction horizon')
                grid on;
                title(names{i});
            end
        end
        
        function status = updateTrackAnimation(obj,k)
        % -----------------------------------------------------------------
        %   Update track animation with the current time step k. Beware
        %   that the vectors obj.mu_x_pred_opt and obj.x_ref must have the
        %   correct values at position (:,:,k)
        % -----------------------------------------------------------------
            status = 0;
            hold on
            if k < 1 || k > size(obj.mu_x_pred_opt,3) || any(isnan(obj.mu_x_pred_opt(:,:,k)),'all')
                return;
            end
           axis([-3+obj.mu_x_pred_opt(7,1,k),30+obj.mu_x_pred_opt(7,1,k),-5,5]);     
            vel = vecnorm(obj.mu_x_pred_opt(4:5,:,k));
            % update predicted trajectory
            obj.h_mu_x_pred_opt.XData = [obj.mu_x_pred_opt(1,:,k) 0];
            obj.h_mu_x_pred_opt.YData = [obj.mu_x_pred_opt(2,:,k) NaN];
            obj.h_mu_x_pred_opt.CData = [vel min(vel)];
            
            % update state covariances
            for i=1:obj.N
                mean = obj.mu_x_pred_opt(1:2,i,k);
                var  = obj.var_x_pred_opt(1:2,1:2,i,k);
                if isnan(var)
                    break
                else
                    if all(svd(var)>1e-6)
                        ell = sigmaEllipse2D(mean, var, obj.ell_level, obj.ell_npoints);
                    else
                        ell = repmat(mean,1,obj.ell_npoints);
                    end
                    obj.h_var_x_pred_opt{i}.Vertices = ell';
                end
            end
            
            % update projected reference
            obj.h_x_ref.XData = obj.x_ref(1,:,k);
            obj.h_x_ref.YData = obj.x_ref(2,:,k);
            
            % update trace
            veltrace = vecnorm(squeeze(obj.mu_x_pred_opt(4:5,1,1:k)));
            obj.h_x_trace.XData = [squeeze(obj.mu_x_pred_opt(1,1,1:k))' NaN];
            obj.h_x_trace.YData = [squeeze(obj.mu_x_pred_opt(2,1,1:k))' NaN];
            obj.h_x_trace.CData = [veltrace NaN];
            
            % update car
            carpos = obj.mu_x_pred_opt(1:2,1,k); %[0;0]
            psi    = obj.mu_x_pred_opt(3,1,k); %deg2rad(30);
            car_w = 1;
            car_l = 2;
            V_carpoints = [[car_l/2;car_w/2],[car_l/2;-car_w/2],[-car_l/2;-car_w/2],[-car_l/2;car_w/2]];
            I_carpoints = [cos(psi) -sin(psi);
                           sin(psi)  cos(psi)] * V_carpoints + carpos;
            obj.h_car.Vertices = I_carpoints';
            
             % update obstacles cars
             for i=1: size(obj.bias,2)
                carpos_obs = obj.XY_o(:,i,k); %[0;0]
                psi_obs    = obj.PSIt_o(:,i,k); %deg2rad(30);
                car_w = 1;
                car_l = 2;
                V_carpoints = [[car_l/2;car_w/2],[car_l/2;-car_w/2],[-car_l/2;-car_w/2],[-car_l/2;car_w/2]];
                I_carpoints_obs = [cos(psi_obs ) -sin(psi_obs );  sin(psi_obs )  cos(psi_obs )] * V_carpoints + carpos_obs; %+[cos(psi_obs ) -sin(psi_obs );  sin(psi_obs )  cos(psi_obs )]*[0; obj.bias(2,i)] ;%+ repmat( [cos(psi_obs ) -sin(psi_obs );  sin(psi_obs )  cos(psi_obs )]*[0; obj.bias(2,i)], [1,4]) ;
                obj.h_obs_cars{i}.Vertices= I_carpoints_obs';
                
                obj.h_x_trace_obs{i}.XData = [squeeze(obj.XY_o(1,i,1:k))' NaN];
                obj.h_x_trace_obs{i}.YData = [squeeze(obj.XY_o(2,i,1:k))' NaN];
                obj.h_x_trace_obs{i}.CData = [repmat(obj.v_obs(1,i),[1,k]) NaN];
                
                a1=car_l;b1=car_w;  %车长；车宽
                a=sqrt(a1^2+b1^2); b= a*b1/a1;
                t=linspace(0,2*pi,obj.N);
                xx=a*sin(t);
                yy=b*cos(t);
                phi=psi_obs;
                A=[cos(phi) -sin(phi);sin(phi) cos(phi)]*[xx;yy]+carpos_obs;
                obj.h_obs_ellipse{i}.XData = A(1,:);
                obj.h_obs_ellipse{i}.YData = A(2,:);
             end
                      xy_oa = obj.XY_o(:,:,k);  %虚拟约束for凸化非凸集
                      Dis_obs=squareform(pdist(xy_oa'))+10*tril(ones(nchoosek(size(obj.bias,2),2),nchoosek(size(obj.bias,2),2)),0);  %计算各个障碍车之间的距离   
                      Dis_min = min(min(Dis_obs));
                      [indx,indy]=find(Dis_min ==Dis_obs);
                      for j = 1:size(indx,1)
                          carpos_obs = (xy_oa(:,indx(j,1))+xy_oa(:,indy(j,1)))/2;
                          phi =atan((xy_oa(2,indx(j,1))-xy_oa(2,indy(j,1)))/(xy_oa(1,indx(j,1))-xy_oa(1,indy(j,1))));
    %                       xx =2*xx;
                          A=[cos(phi) -sin(phi);sin(phi) cos(phi)]*[xx;yy]+carpos_obs;
                          obj.h_obs_ellipse{i+j}.XData = A(1,:);
                          obj.h_obs_ellipse{i+j}.YData = A(2,:);
                      end
            % no error when updating graphics
            status = 1;
        end
        
        
        function updateScope(obj,k)
        % -----------------------------------------------------------------
        %   Update scope with signals from the current time step k. Beware
        %   that the vectors obj.mu_x_pred_opt and obj.u_pred_opt must have 
        %   the correct values at position (:,:,k)
        % -----------------------------------------------------------------
            obj.k = k;
            refreshdata(obj.h_scopex,'caller');
            refreshdata(obj.h_scopeu,'caller');
        end
        
        function recordvideo(obj, videoName, format, FrameRate)
        % -----------------------------------------------------------------
        %   Record video of the track animation
        % -----------------------------------------------------------------
            % video rec
            videoframes = struct('cdata',[],'colormap',[]);
            obj.initTrackAnimation();
            xlim manual
            ylim manual
            for k=1:size(obj.mu_x_pred_opt,3)
                status = obj.updateTrackAnimation(k);
                if status == 0
                    break;
                end
                videoframes(k) = getframe(obj.h_fig);
            end
            % -----------------------------------------------------------------
            %   Save video
            % -----------------------------------------------------------------
            writerObj = VideoWriter(videoName,format);
            writerObj.FrameRate = FrameRate;
            open(writerObj);
            % Write out all the frames.
            numberOfFrames = length(videoframes);
            for k=1:numberOfFrames 
               writeVideo(writerObj, videoframes(k));
            end
            close(writerObj);
            disp('Video saved successfully')
        end
    end
end