clc;
clear;
close all
load C.mat
load A_d.mat
load B_d.mat
load D.mat


Ts = 1e-3;
sys_d = ss(Ad,Bd,C1,0);
%% Question C Kalman Filter
sys = ss(Ad,[Bd Bd],C1,0,Ts);% State space model discrete  system.
QN = [(0.3/3)^2 0; 0 (0.1/3)^2]; %covariance matrix of Process Noise
RN = [(0.02/3)^2 0; 0 (0.01/3)^2];  %covariance matrix of Measurment Noise
NN = zeros(2,2); 
[KEST,L,P1] = kalman(sys,QN,RN,NN,'current'); % Kalman filter KEST - Estimates states L - Kalman gain 
eigen = eig(Ad - L*C1); % observer stability

%% Question D Referance Gain
clc
Qx = diag([1 1 1 150 1/700]);
Qu = [100 0;0 50];
N = 0;
K = dlqr(Ad,Bd,Qx,Qu,N); % LQR Gain calculation
% Kr1=-inv(C1*inv(eye(5)-Ad+(Bd*K))*Bd);
% Kr1=diag(Kr1)/10
Kr1=-inv(C1*inv(eye(5)-Ad+(Bd*K))*Bd);

Kr1=diag(Kr1)*1e-14;
%% Question D Integral Action
CI = [0 0 0 0 1];
A = [Ad,zeros(5,1);-CI, 1];
B = [Bd;zeros(1,2)];
C = [CI,zeros(1,1)];
D = D_1;
ctr =  ctrb(A,B);
Rank =  rank(ctr);
cond_number = cond(ctr);

sys_extend = ss(A,B,C,0,Ts);
% Calculate LQI controller
Qx_e = diag([1 1 1 1 100 1/700]);
Qu_e = [1000 0;0 500];

K2 = dlqr(A, B, Qx_e, Qu_e);
% Set variables for Simulink
K_p = K2(:,1:5);
K_I = K2(:,6)*150;
