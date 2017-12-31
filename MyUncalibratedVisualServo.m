classdef MyUncalibratedVisualServo < handle
    properties
        P
        uv_star
        Tcam
        camera
        robot
        T0
        P0
        
        niter
        fps
        eterm
        lambda
        history
        
        image
        jacobian
        jointpos
        
        preimage
        prejacobian
        prejoint
        
        arglist
    end
    
    methods
        function vs = MyUncalibratedVisualServo(cam, urrobot, varargin)
            vs.camera = cam;
            vs.robot = urrobot;
            vs.history = [];
            
            opt.niter = [];
            opt.fps = 5;
            opt.pose0 = cam.T;
            opt.P = [];
            opt.target = [];
            opt.pstar = [];
            
            opt.eterm = 0.5;
            opt.lambda = 0.08;
            
            [opt, vs.arglist] = tb_optparse(opt, varargin);
            vs.niter = opt.niter;
            vs.fps = opt.fps;
            if ~isempty(opt.pose0)
                vs.T0 = SE3(opt.pose0);
                vs.P0 = urrobot.ikine(opt.pose0);
                vs.jointpos = urrobot.ikine(opt.pose0);
            end
            
            vs.P = opt.target;
            vs.lambda = opt.lambda;
            vs.eterm = opt.eterm;
            vs.uv_star = opt.pstar;
        end
        
        function [J0, newinitjoint, newinitimage] = initjacobian(vs)
%             theta = [];
%             F = [];
%             vs.camera.T = vs.T0;
%             uv_p = vs.camera.project(vs.P);
%             uv_p = uv_p(:);
%             joint_p = vs.P0';
%             for i = 1:6
%                 deltatheta = [0 0 0 0 0 0]';
%                 deltatheta(i) = 0.017;
%                 theta = [theta deltatheta];
%                 joint = joint_p + deltatheta;
%                 coor = vs.robot.fkine(joint);
%                 joint_p = joint;
%                 imageuv = vs.camera.project(vs.P, 'pose', coor);
%                 imageuv = imageuv(:);
%                 F = [F imageuv-uv_p];
%                 uv_p = imageuv;
%             end
%             J0 = F / theta;
%             newinitjoint = joint_p';
%             newinitimage = uv_p;
              vs.camera.T = vs.T0;
              uv_p = vs.camera.project(vs.P);
              pt = inv(vs.T0) * vs.P;
              J0 = vs.camera.visjac_p(uv_p, pt(3, :));
              uv_p = uv_p(:);
              newinitjoint = vs.P0;
              newinitimage = uv_p;
        end
        
        function init(vs)
            
            [vs.jacobian, vs.jointpos, vs.image] = vs.initjacobian();
            vs.prejacobian = vs.jacobian;
            vs.preimage = vs.image;
            vs.prejoint = vs.jointpos;
            
            vs.camera.T = vs.robot.fkine(vs.jointpos);
            vs.Tcam = vs.camera.T;
            
            vs.camera.plot(vs.P);
            pause(1)
            
            plot_sphere(vs.P, 0.05, 'b')
            lighting gouraud
            light
            vs.camera.plot_camera(vs.P, 'label', 'scale', 0.3);
            
            vs.robot.plot(vs.jointpos, 'tilesize', 1, 'jointdiam', 2, 'basewidth', 5);
            
            hist.uv = vs.image;
            e = vs.image - vs.uv_star(:);
            hist.e = e(:);
            hist.jointpos = vs.jointpos;
            hist.jacobian = vs.jacobian;
            vs.history = [vs.history hist];
        end
        
        function newjacobian = broyden_update(vs)
            deltaimage = vs.image - vs.preimage;
            deltajoint = vs.jointpos - vs.prejoint;
            newjacobian = vs.prejacobian + (deltaimage - vs.prejacobian * deltajoint')...
                * deltajoint / (deltajoint * deltajoint');
        end
        
        function status = step(vs)
            status = 0;
            
            e = vs.image - vs.uv_star(:);
            e = e(:);
            try
                v = -pinv(vs.jacobian) * e;
            catch
                status = -1;
                return
            end
            
            vs.jointpos = vs.prejoint + v';
            vs.robot.plot(vs.jointpos, 'tilesize', 1, 'jointdiam', 2, 'basewidth', 5);
            vs.Tcam = vs.robot.fkine(vs.jointpos);
            vs.camera.T = vs.Tcam;
            vs.camera.plot(vs.P);
            uv = vs.camera.project(vs.P);
            vs.image = uv(:);
            
            vs.jacobian = vs.broyden_update();
            
            vs.prejoint = vs.jointpos;
            vs.preimage = vs.image;
            vs.prejacobian = vs.jacobian;
            
            hist.e = e;
            hist.jointpos = vs.prejoint;
            hist.uv = vs.preimage;
            hist.jacobian = vs.prejacobian;
            vs.history = [vs.history hist];
            
            if norm(e) < vs.eterm
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
        
        function plot_p(vs)
        end
    end
end