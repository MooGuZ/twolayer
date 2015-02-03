% This script trying to read different kind of transformation data and
% train Transformation Pattern Model

fdata  = '/home/hxz244/TPModel/#data/Transformation/';
fmodel = '/home/hxz244/TPModel/#models/';

% Motion with 23 pattern kernel and 17 transformation kernel
load([fmodel,'20150129/B.mat']);
M = tpmodel([fdata,'motion/'],23,17,30,10,5,M);  save([fmodel,'20150202/motion.mat'],'M');
M = tpmodel([fdata,'motion/'],23,17,30,20,10,M); save([fmodel,'20150202/motion.mat'],'M');
M = tpmodel([fdata,'motion/'],23,17,30,50,15,M); save([fmodel,'20150202/motion.mat'],'M');
clear('M');

% Rotation with 23 pattern kernel and 17 transformation kernel
load([fmodel,'20150129/C.mat']);
M = tpmodel([fdata,'rotation/'],23,17,30,10,5,M);  save([fmodel,'20150202/rotation.mat'],'M');
M = tpmodel([fdata,'rotation/'],23,17,30,20,10,M); save([fmodel,'20150202/rotation.mat'],'M');
M = tpmodel([fdata,'rotation/'],23,17,30,50,15,M); save([fmodel,'20150202/rotation.mat'],'M');
clear('M');

% Rotation with 23 pattern kernel and 17 transformation kernel
load([fmodel,'20150129/D.mat']);
M = tpmodel([fdata,'scaling/'],23,17,30,10,5,M);  save([fmodel,'20150202/scaling.mat'],'M');
M = tpmodel([fdata,'scaling/'],23,17,30,20,10,M); save([fmodel,'20150202/scaling.mat'],'M');
M = tpmodel([fdata,'scaling/'],23,17,30,50,15,M); save([fmodel,'20150202/scaling.mat'],'M');
clear('M');