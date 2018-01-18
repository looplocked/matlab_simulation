clear
clc
cam = CentralCamera('default');
T_C0 = SE3(1,1,-3)*SE3.Rz(0.6);
Cd_T_G = SE3(0, 0, 1);
pd = bsxfun(@plus, 200*[-1 -1 1 1; 1 -1 -1 1], cam.pp')
T_C0 = SE3(1,1,-3)*SE3.Rx(0.6);
%ibvs = UncalibratedVisualServo(cam, 'pose0', T_C0, 'pstar', pd)

T0 = SE3(-0.5, -0.5, 0.5) * SE3.Rx(-pi) * SE3.Rz(pi/4);
T1 = SE3(-0.7, -0.7, 0) * SE3.Rx(-pi) * SE3.Rz(-pi/4);
T2 = SE3(-0.8, -0.6, 0) * SE3.Rx(-pi) * SE3.Rz(pi/2);
P0 = [0.5514 -1.1683 1.3407 1.3984 1.5708 -0.2340]';
P1 = [0.6190 -0.3503 0.7929 1.1282 1.5708 -1.7372]';

T3 = SE3(-0.9, 0, 0.1) * SE3.Rx(-pi*5/4) * SE3.Ry(pi/16) * SE3.Rz(pi/4);

L1=Link('d',0.1273,'a',0,'alpha',1.570796327);
L2=Link('d',0,'a',-0.612,'alpha',0);
L3=Link('d',0,'a',-0.5723,'alpha',0);
L4=Link('d',0.163941,'a',0,'alpha', 1.570796327);
L5=Link('d',0.1157,'a',0,'alpha',-1.570796327);
L6=Link('d',0.0922,'a',0,'alpha',0);
robot = SerialLink([L1,L2,L3,L4,L5,L6], 'name', 'URRobot');

psize = 0.5;
ppos = SE3(-0.7, -0.7, -1);
ptarget = mkgrid(2, 0.5, 'pose', ppos);
uibvs = MyUncalibratedVisualServo(cam, robot, 'pose0', T3, 'pstar', pd, 'target', ptarget, 'lambda', 0.02, 'eterm', 0.5);
uibvs.run();
uibvs.plot_p();

% uibvs = UncalibratedVisualServo(cam, robot, 'pose0', T0, 'pstar', pd, 'lambda', 0.002, 'eterm', 0.5)
% uibvs.run()
% figure()
% uibvs.plot_p();
% uibvs = UncalibratedVisualServo(cam, robot, 'pose0', T1, 'pstar', pd);
% uibvs.run()
% uibvs.plot_p();
% uibvs = UncalibratedVisualServo(cam, robot, 'pose0', T2, 'pstar', pd, 'niter', 10);
% uibvs.run()
% uibvs.plot_p();

