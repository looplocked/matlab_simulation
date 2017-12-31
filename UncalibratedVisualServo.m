classdef UncalibratedVisualServo < handle
    properties
        P
        uv_star
        Tcam
        camera
        URrobot
        jointpos
        Tf
        T0
        P0
        pf
        
        history
        
        niter
        fps
        
        verbose
        arglist
        axis
        
        lambda
        eterm
        uv_p
        
        depth
        depthest
        vel_p
        theta
        smoothing
        
        jacobian
    end
    
    methods
        function vs = UncalibratedVisualServo(cam, robot, varargin)
            vs.camera = cam;
            vs.URrobot = robot;
            vs.history = [];

            z = 3;
            opt.niter = [];
            opt.fps = 5;
            opt.posef = [];
            opt.pose0 = cam.T;
            opt.P = [];
            opt.targetsize = 0.5;       % dimensions of target
            opt.pstar = [];
            opt.axis = [];
            
            opt.eterm = 0.5;
            opt.lambda = 0.08;         % control gain
            opt.depth = [];
            opt.depthest = false;
            
            [opt,vs.arglist] = tb_optparse(opt, varargin);
            vs.niter = opt.niter;
            vs.fps = opt.fps;
            vs.verbose = opt.verbose;
            if ~isempty(opt.pose0)
                vs.T0 = SE3(opt.pose0);
                vs.P0 = robot.ikine(SE3(opt.pose0));
                vs.jointpos = robot.ikine(SE3(opt.pose0));
            end
            if ~isempty(opt.posef)
                vs.Tf = SE3(opt.posef);
            end

            % define feature points in XY plane, make vertices of a square
            if isempty(opt.P)
                targetpos = SE3(-0.7, -0.7, -1);
                opt.P = mkgrid(2, opt.targetsize, 'pose', targetpos);
            end
            vs.P = opt.P;
            vs.pf = opt.pstar;
            vs.axis = opt.axis;
            
            vs.lambda = opt.lambda;
            vs.eterm = opt.eterm;
            vs.theta = 0;
            vs.smoothing = 0.80;
            vs.depth = opt.depth;
            vs.depthest = opt.depthest;
        end
        
        function init(vs)
            if ~isempty(vs.pf)
                vs.uv_star = vs.pf
            else
                if ~isempty(vs.Tf)
                    vs.Tf = transl(0, 0, 1);
                    warning('setting Tf to default');
                end
                vs.uv_star = vs.camera.project(vs.P, 'Tcam', inv(vs.Tf))
            end
            
            Theta = [];
            F = [];
            vs.camera.T = vs.T0;
            preuv = vs.camera.project(vs.P);
            preuv = preuv(:);
            prejoint = vs.P0';
            for i = 1:6
                deltatheta = [0 0 0 0 0 0]';
                deltatheta(i) = 0.017;
                Theta = [Theta deltatheta];
                joint = prejoint + deltatheta;
                coor = vs.URrobot.fkine(joint);
                prejoint = joint;
                imageuv = vs.camera.project(vs.P, 'pose', coor);
                imageuv = imageuv(:);
                F = [F imageuv-preuv];
                preuv = imageuv;
            end
            
            vs.jacobian = F / Theta;
            vs.jointpos = prejoint';
            
            vs.camera.T = vs.URrobot.fkine(vs.jointpos);
            vs.Tcam = vs.camera.T;
            
            if 0
                vs.camera.clf()
                vs.camera.plot(vs.uv_star, '*')
                vs.camera.hold(true);
                vs.camera.plot(vs.P, 'Tcam', vs.T0, 'o');
                pause(2)
                vs.camera.hold(false);
                vs.camera.clf();
            end
            
            vs.camera.plot(vs.P);
            pause(1)
            
            figure
            if ~isempty(vs.axis)
                axis(vs.axis);
            end
            plot_sphere(vs.P, 0.05, 'b')
            lighting gouraud
            light
            vs.camera.plot_camera(vs.P, 'label', 'scale', 0.3);
            %hold on;
            
            %vs.jointpos = vs.URrobot.ikine(vs.T0);
            vs.URrobot.plot(vs.jointpos, 'tilesize', 1, 'jointdiam', 2, 'basewidth', 5);
            
            vs.vel_p = [];
            vs.uv_p = [];
            vs.history = [];
            hist.uv = preuv;
        end
        
        function broyden_update(vs)
            deltajoint = vs.jointpos - vs.history(:, end).jointpos;
            uv = vs.camera.plot(vs.P);
            uv = uv(:);
            deltaimage = uv - vs.history(:, end).uv;
            vs.jacobian = vs.jacobian + vs.lambda * (deltaimage - vs.jacobian * deltajoint')...
                * deltajoint / (deltajoint * deltajoint');
        end
        
        function rls = rls_update(vs)
        end
        
        function status = step(vs)
            status = 0;
            Zest = [];
            
            % compute the view
            uv = vs.camera.plot(vs.P);

            % optionally estimate depth
            if vs.depthest
                % run the depth estimator
                [Zest,Ztrue] = vs.depth_estimator(uv);
                if vs.verbose
                    Zest
                    Ztrue
                end
                vs.depth = Zest;
                hist.Ztrue = Ztrue(:);
                hist.Zest = Zest(:);
            end

            % compute image plane error as a column
            e = uv - vs.uv_star;   % feature error
            e = e(:);
        
            
            % compute the Jacobian
            if isempty(vs.depth)
                % exact depth from simulation (not possible in practice)
                pt = inv(vs.Tcam) * vs.P;
                J = vs.jacobian;
            elseif ~isempty(Zest)
                J = vs.camera.visjac_p(uv, Zest);
            else
                J = vs.camera.visjac_p(uv, vs.depth);
            end

            % compute the velocity of camera in camera frame
            try
                v = -vs.lambda * pinv(J) * e;
            catch
                status = -1;
                return
            end

            if vs.verbose
                fprintf('v: %.3f %.3f %.3f %.3f %.3f %.3f\n', v);
            end

            Td = SE3(trnorm(delta2tr(v)));
            %vs.Tcam = vs.Tcam .* Td; 
            vs.jointpos = vs.jointpos + v';
            vs.Tcam = vs.URrobot.fkine(vs.jointpos);
            vs.camera.T = vs.Tcam;
            hist.uv = uv(:);
            vel = tr2delta(Td);
            hist.vel = vel;
            hist.e = e;
            hist.en = norm(e);
            hist.jcond = cond(J);
            hist.Tcam = vs.Tcam;
            
            vs.URrobot.plot(vs.jointpos, 'tilesize', 1, 'jointdiam', 2, 'basewidth', 5);
            
            hist.jointpos = vs.jointpos;
            hist.jacobian = vs.jacobian;

            vs.history = [vs.history hist]; 

            vs.vel_p = vel;
            vs.uv_p = uv;
            
            vs.broyden_update();

            if norm(e) < vs.eterm,
                status = 1;
                return
            end
        end
        
        function run(vs, nsteps)
            vs.init();
            
            if nargin < 2
                nsteps = vs.niter;
            end
            ksteps = 0;
            while true
                ksteps = ksteps + 1;
                status = vs.step();
                
                drawnow
                pause(1/vs.fps)
                
                if status > 0
                    fprintf('completed on error tolerance\n');
                    break;
                elseif status < 0
                    fprintf('failed on error\n');
                    break;
                end
                
                if ~isempty(nsteps) && (ksteps > nsteps)
                    break;
                end
            end
            
            if status == 0
                fprintf('completed on iteration count\n');
            end
        end
        
        function [Zest, Ztrue] = depth_estimator(vs, uv)
            if isempty(vs.uv_p)
                Zest = [];
                Ztrue = [];
                return;
            end
            
            J = vs.camera.visjac_p(uv, 1);
            Jv = J(:, 1:3);
            Jw = J(:, 4:6);
            uv_d = uv(:) - vs.uv_p(:);
            B = uv_d - Jw*vs.vel_p(4:6);
            A = Jv * vs.vel_p(1:3);
            
            AA = zeros(numcols(uv), numcols(uv)/2);
            for i = 1:numcols(uv)
                AA(i*2-1:i*2, i) = A(i*2-1:i*2);
            end
            eta2 = AA\B;
            eta2 = A(1:2) \ B(1:2);
            vs.theta = (1-vs.smoothing) * 1./eta' + vs.smoothing * vs.theta;
            Zest = vs.theta;
            
            P_CT = inv(vs.Tcam) * vs.P;
            Ztrue = P_CT(3, :);
            
            if vs.verbose
                fprintf('depth %.4g, est depth %.4g, rls depth %.4g\n', ...
                    Ztrue, 1/eta, Zest);
            end
        end
        
        function plot_p(vs)
            if isempty(vs.history)
                return
            end
            clf
            hold on
            % image plane trajectory
            uv = [vs.history.uv]'; 
            % result is a vector with row per time step, each row is u1, v1, u2, v2 ...
            for i=1:numcols(uv)/2
                p = uv(:,i*2-1:i*2);    % get data for i'th point
                plot(p(:,1), p(:,2), 'b')
            end
            
            % mark the initial target shape
            plot_poly( reshape(uv(1,:), 2, []), 'o--');
            uv(end,:)
            
            % mark the final target shape
            if ~isempty(vs.uv_star)
                plot_poly(vs.uv_star, 'rh:', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
            else
                plot_poly( reshape(uv(end,:), 2, []), 'rh--', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            end
            axis([0 vs.camera.npix(1) 0 vs.camera.npix(2)]);
            daspect([1 1 1])
            set(gca, 'Ydir' , 'reverse');
            grid
            xlabel('u (pixels)');
            ylabel('v (pixels)');
            hold off
        end
        
        function plot_vel(vs)
            if isempty(vs.history)
                return
            end
            clf
            vel = [vs.history.vel]';
            plot(vel(:,1:3), '-')
            hold on
            plot(vel(:,4:6), '--')
            hold off
            ylabel('Cartesian velocity')
            grid
            xlabel('Time step')
            xaxis(length(vs.history));
            legend('v_x', 'v_y', 'v_z', '\omega_x', '\omega_y', '\omega_z')
        end
        
        function plot_camera(vs)
            if isempty(vs.history)
                return
            end
            clf
            T = [vs.history.Tcam];
            subplot(211)
            plot(T.tv');
            xaxis(length(vs.history));
            ylabel('Camera position')
            legend('X', 'Y', 'Z');
            grid
            
            subplot(212)
            plot(T.torpy)
            ylabel('Camera orientation (rad)')
            grid
            xlabel('Time step')
            xaxis(length(vs.history));
            legend('R', 'P', 'Y');
        end
        
        function plot_jcond(vs)
            if isempty(vs.history)
                return
            end
            clf

            Jcond = [vs.history.jcond];
            % Image Jacobian condition number vs time
            plot(Jcond);
            grid
            ylabel('Jacobian condition number')
            xlabel('Time step')
            xaxis(length(vs.history));
        end
        
        function plot_z(vs)
            if isempty(vs.history)
                return
            end
            clf
            Zest = [vs.history.Zest];
            Ztrue = [vs.history.Ztrue];
            plot(Ztrue', '-')
            hold on
            set(gca, 'ColorOrderIndex', 1);
            plot(Zest', '--')
            grid
            ylabel('Depth (m)')
            xlabel('Time step')
            xaxis(length(vs.history));
        end
        
        function plot_error(vs)
            if isempty(vs.history)
                return
            end
            clf
            e = [vs.history.e]';
            plot(e(:,1:2:end), 'r');
            hold on
            plot(e(:,2:2:end), 'b');
            hold off
            ylabel('Feature error (pixel)')
            grid
            xlabel('Time')
            xaxis(length(vs.history));
            legend('u', 'v');
        end
        
        function plot_all(vs, name, dev)
            if nargin < 3
                dev = '-depsc';
            end
            if nargin < 2
                name = [];
            end

            figure
            vs.plot_p();
            if ~isempty(name)
                print(gcf, dev, sprintf(name, '-p'));
            end

            figure
            vs.plot_vel();
            if ~isempty(name)
                print(gcf, dev, strcat(name, '-vel'));
            end

            figure
            vs.plot_camera();
            if ~isempty(name)
                print(gcf, dev, strcat(name, '-camera'));
            end
            figure
            vs.plot_error();
            if ~isempty(name)
                print(gcf, dev, strcat(name, '-error'));
            end
            
            % optional plots depending on what history was recorded
            if isfield(vs.history, 'Zest')
                figure
                vs.plot_z();
                if ~isempty(name)
                    print(gcf, dev, strcat(name, '-z'));
                end
            end
            
            if isfield(vs.history, 'jcond')
                figure
                vs.plot_jcond();
                if ~isempty(name)
                    print(gcf, dev, strcat(name, '-jcond'));
                end
            end
        end
        
        function display(vs)
            loose = strcmp( get(0, 'FormatSpacing'), 'loose');
            if loose
                disp(' ');
            end
            disp([inputname(1), ' = '])
            disp( char(vs) );
        end
        
        function s = char(vs)
            s = sprintf('Visual servo object: camera=%s\n  %d iterations, %d history', ...
                vs.camera.name, vs.niter, length(vs.history));
            s = strvcat(s, [['  P= '; '     '; '     '] num2str(vs.P)]);
            if 0
            s = strvcat(s, sprintf('  T0:'));
            s = strvcat(s, [repmat('      ', 4,1) num2str(vs.T0)]);
            s = strvcat(s, sprintf('  C*_T_G:'));
            s = strvcat(s, [repmat('      ', 4,1) num2str(vs.Tf)]);
            else
                if ~isempty(vs.T0)
                    s = strvcat(s, ['  C_T0:   ' trprint(vs.T0, 'fmt', ' %g', 'angvec')]);
                end
                if ~isempty(vs.Tf)
                    s = strvcat(s, ['  C*_T_G: ' trprint(vs.Tf, 'fmt', ' %g', 'angvec')]);
                end
            end
        end
    end
end
        