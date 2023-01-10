clear 
clc
close

syms R D1 R_0 D1_0 real positive %non-zero positive and real symbols for R and D1.
assume(R ~= 0)
assume(D1 ~= 0)
assume(R_0 ~= 0)
assume(D1_0 ~= 0)

%Given parameters 
Ke=0.1;
Kt=0.1;
J1=1e-5;
J2=4e-5;
Bf=2e-3;
D2=2;


%Matrices A, B, C and D for the state-space modelling

A= [0       0       0       1       0;
    0       0       0       0       1;
    0     D2/Bf   -D2/Bf    0       0;
    -D1/J1 D1/J1    0      (-Kt*Ke)/(R*J1) 0;
    D1/J2 -(D1+D2)/J2 D2/J2 0       0];

B= [0       0;
    0       0;
    0       1/Bf;
    Kt/(R*J1) 0;
    0       0]; 

C1=[0 1 0 0 0;
    0 0 0 0 1];

C2=[0 0 0 -Ke/R 0;
    0 D2/Bf -D2/Bf 0 0];

D_1=zeros(2,2);

D_2=[1/R 0;
    0 1/Bf];

%Controllbility matrix and its reduced-row echelon form.
contr=[B A*B A^2*B A^3*B A^4*B];
rref_c=rref(contr);
Rank_contr=rank(contr);

%Observability matrices and their reduced-row echelon forms. 

obsvr1=[C1; C1*A; C1*A^2; C1*A^3; C1*A^4];
rref_o1=rref(obsvr1);
Rank_obs1=rank(obsvr1);

obsvr2=[C2; C2*A; C2*A^2; C2*A^3; C2*A^4];
rref_o2=rref(obsvr2);
Rank_obs2=rank(obsvr2);

[R_0, D1_0] = solve(det(obsvr1.'*obsvr1),[R D1]);

%% Question B

%Eigen values for the A matrix
eigen=eig(A);

SA=zeros(5,7);
stab=zeros(5,1);
Dec1=zeros(5,1);
Dec2=zeros(5,1);

%PBH test for stability and detectability
for i=1:length(eigen)
    SA=[A-eigen(i)*eye(5,5), B];
    stab(i)=rank(SA);

    DA1=[A-eigen(i)*eye(5,5); C1];
    Dec1(i)=rank(DA1);

    DA2=[A-eigen(i)*eye(5,5); C2];
    Dec2(i)=rank(DA2);
end
%% Question C

R=1;
D1=20;
Ke=0.1;
Kt=0.1;
J1=1e-5;
J2=4e-5;
Bf=2e-3;
D2=2;

A_new = [0       0       0       1       0;
    0       0       0       0       1;
    0     D2/Bf   -D2/Bf    0       0;
    -D1/J1 D1/J1    0      (-Kt*Ke)/(R*J1) 0;
    D1/J2 -(D1+D2)/J2 D2/J2 0       0];

B_new = [0       0;
    0       0;
    0       1/Bf;
    Kt/(R*J1) 0;
    0       0]; 

C2_new = [0 0 0 -Ke/R 0;
    0 D2/Bf -D2/Bf 0 0];

D1_new=zeros(2,2);

D2_new=[1/R 0;
    0 1/Bf];
contr_m = ctrb(A_new,B_new); % controllability using Matlab function
rref_c_m = rref(contr_m);
Rank_contr_m = rank(rref_c_m);  % Rank controllability
cond_contr = cond(contr_m); % condition number for controllability

obsvr1_m = obsv(A_new,C1); % observability using Matlab function
rref_o1_m = rref(obsvr1_m);
Rank_obs1_m = rank(rref_o1_m); % Rank observability
cond_obs1 = cond(obsvr1_m); % condition number for observability

obsvr2_m = obsv(A_new,C2_new); % observability using Matlab function
rref_o2_m = rref(obsvr2_m);
Rank_obs2_m = rank(rref_o2_m);  % Rank observability
cond_obs2 = cond(obsvr2_m); % condition number for observability


%% Question D

Ts =1/1000; % sampling interval[ms]
Ad = expm(A_new*Ts); % Converting continues to discrete time

%% Question E 

fun = @(t) expm(A_new*t)*B_new;
Bd = integral(fun,0,0.001,'ArrayValued', 1); % Converting continues to discrete time

%% Question F

% Case-1
ss1= ss(Ad,Bd,C1,D1_new); % State Space for discrete time
P1 = pole(ss1); % Pole for discrete time
Z1 = tzero(ss1);% Zeros for discrete time

figure(1)
zplane(Z1,P1)

contr_md = ctrb(Ad,Bd); % controllability using Matlab function
Rank_contr_md=rank(contr_md);  % Rank controllability

obsvr1_md = obsv(ss1); % observability using Matlab function
Rank_obsvr1_md=rank(obsvr1_md); % Rank observability

%Case-2
ss2= ss(Ad,Bd,C2_new,D2_new); % State Space for discrete time
P2 = pole(ss2); % Pole for discrete time
Z2 = tzero(ss2); % Zeros for discrete time

figure(2)
zplane(Z2,P2)

obsvr2_md = obsv(ss2); % observability using Matlab function
Rank_obsvr2_md=rank(obsvr2_md); % Rank observability

