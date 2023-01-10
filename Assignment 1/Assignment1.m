clear 
close 
clc
%% Given variables

R=1;
Ke=0.1;
Kt=0.1;
J1=1e-5;
J2=4e-5;
Bf=2e-3;
D1=20;
D2=2;

%% Question D

A= [0       0       0       1       0;
    0       0       0       0       1;
    0     D2/Bf   -D2/Bf    0       0;
    -D1/J1 D1/J1    0      (-Kt*Ke)/(R*J1) 0;
    D1/J2 -(D1+D2)/J2 D2/J2 0       0];

eigenValues = eig(A);

B= [0       0;
    0       0;
    0       1/Bf;
    Kt/(R*J1) 0;
    0       0];

C1=[0 1 0 0 0;
    0 0 0 0 1];

C2=[0 0 0 -Ke/R 0;
    0 D2/Bf -D2/Bf 0 0];

D1=zeros(2,2);

D2=[1/R 0;
    0 1/Bf];

ss1= ss(A,B,C1,D1);
ss2= ss(A,B,C2,D2);

z1=tzero(ss1);
z2=tzero(ss2);

pzmap(ss1,ss2)
%% Question E

B2=[0;
    0;
    0;
    Kt/(R*J1);
    0];

D3 = [1/R;
      0];

ss3= ss(A,B2,C2,D3);
tf3=tf(ss3);

z3=tzero(ss3);
pzmap(ss3)

%% Simulation
ASim=A;
BSim=B;
C1Sim=C1;
C2Sim=C2;
D1Sim=D1;
D2Sim=D2;
