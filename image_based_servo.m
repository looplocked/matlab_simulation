cam = CentralCamera('default');
T_C0 = SE3(1,1,-3)*SE3.Rz(0.6);
Cd_T_G = SE3(0, 0, 1);
pd = bsxfun(@plus, 200*[-1 -1 1 1; -1 1 1 -1], cam.pp');
T_C0 = SE3(1,1,-3)*SE3.Rx(0.6);
%ibvs = UncalibratedVisualServo(cam, 'pose0', T_C0, 'pstar', pd)
 
T0 = SE3(-2.1, 0, -3) * SE3.Rz(5*pi/4);

uibvs = UncalibratedVisualServo(cam, 'pose0', T0, 'pstar', pd, 'lambda', 0.002, 'eterm', 0.5)
uibvs.run()
uibvs.plot_p();
uibvs = UncalibratedVisualServo(cam, 'pose0', SE3(0, 0, -1)*SE3.Rz(1), 'pstar', pd);
uibvs.run()
uibvs.plot_camera
uibvs = UncalibratedVisualServo(cam, 'pose0', SE3(0, 0, -1)*SE3.Rz(pi), 'pstar', pd, 'niter', 10);
uibvs.run()
uibvs.plot_camera