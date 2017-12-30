clear
clc
clf
L1=Link('d',0.1273,'a',0,'alpha',1.570796327);
L2=Link('d',0,'a',-0.612,'alpha',0);
L3=Link('d',0,'a',-0.5723,'alpha',0);
L4=Link('d',0.163941,'a',0,'alpha', 1.570796327);
L5=Link('d',0.1157,'a',0,'alpha',-1.570796327);
L6=Link('d',0.0922,'a',0,'alpha',0);
URRobot = SerialLink([L1,L2,L3,L4,L5,L6],'name', 'URRobot');

T0 = SE3(-0.5, -0.5, 0.5) * SE3.Rx(-pi) * SE3.Rz(pi/4);
T1 = SE3(-0.7, -0.7, 0) * SE3.Rx(-pi) * SE3.Rz(-pi/4);
T2 = SE3(-0.8, -0.6, 0) * SE3.Rx(-pi) * SE3.Rz(pi/2);
t = 0:0.56:2;
q0 = URRobot.ikine(T0);
q1 = URRobot.ikine(T1);
q = jtraj(q0, q1, t);

cam = CentralCamera('default');
jointpos1 = [pi/2,0.2,0.3,0.4,0.5,0.6];
jointpos2 = [0.1,0.2,0.3,0.4,0.5,0.6];
Q = [-pi/2 -pi/2 -pi*3/4 -pi/4 pi/2 0];
jointtraj = [(0:0.1:1)', (0:0.1:1)', (0:0.1:1)', (0:0.1:1)', (0:0.1:1)', (0:0.1:1)'];

pos = URRobot.ikine(T2)
URRobot.plot(q, 'tilesize', 1, 'jointdiam', 2, 'basewidth', 5);
cam.plot_camera();
cam.T = T0;
cam.T = T1;
%URRobot.plot(jointpos2);