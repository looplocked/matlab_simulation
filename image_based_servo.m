clear
clc
cam = CentralCamera('default');
T_C0 = SE3(1,1,-3)*SE3.Rz(0.6);
Cd_T_G = SE3(0, 0, 1);
pd = bsxfun(@plus, 200*[-1 -1 1 1; 1 -1 -1 1], cam.pp')
T_C0 = SE3(1,1,-3)*SE3.Rx(0.6);
%ibvs = UncalibratedVisualServo(cam, 'pose0', T_C0, 'pstar', pd)

T0 = SE3(-0.5, -0.5, 0.25) * SE3.Rx(-pi) * SE3.Rz(pi/4);
%T0 = SE3(-0.7, -0.7, 0) * SE3.Rx(-pi);

L1=Link('d',0.1273,'a',0,'alpha',1.570796327);
L2=Link('d',0,'a',-0.612,'alpha',0);
L3=Link('d',0,'a',-0.5723,'alpha',0);
L4=Link('d',0.163941,'a',0,'alpha', 1.570796327);
L5=Link('d',0.1157,'a',0,'alpha',-1.570796327);
L6=Link('d',0.0922,'a',0,'alpha',0);
robot = SerialLink([L1,L2,L3,L4,L5,L6], 'name', 'URRobot');

uibvs = UncalibratedVisualServo(cam, robot, 'pose0', T0, 'pstar', pd, 'lambda', 0.02, 'eterm', 0.5)
uibvs.run()
uibvs.plot_p();
uibvs = UncalibratedVisualServo(cam, robot, 'pose0', SE3(0, 0, -1)*SE3.Rz(1), 'pstar', pd);
uibvs.run()
uibvs.plot_camera
uibvs = UncalibratedVisualServo(cam, robot, 'pose0', SE3(0, 0, -1)*SE3.Rz(pi), 'pstar', pd, 'niter', 10);
uibvs.run()
uibvs.plot_camera