clear all; close all; clc;
lc = lines(7);
if ~ exist(fullfile(pwd,'images'),'dir'), mkdir images; end

lt = @(s) clipboard('copy',latex(s));

%% Modelling of the DC Motor

syms v_a T_e R L K_e K_t D_1 D_2 J_1 J_2 B

ia = sym('i_a_',[1,2]);
phi1 = sym('phi_1_',[1,2]);
w1 = sym('omega_1_',[1,2]);
phi2 = sym('phi_2_',[1,2]);
w2 = sym('omega_2_',[1,2]);
phi3 = sym('phi_3_',[1,2]);

x = [ia; phi1; w1; phi2; w2; phi3];
u = [v_a; T_e];


eqs = [       v_a == ia(1)*R + L*ia(2) + K_e*w1(1) ; 
        J_1*w1(2) == K_t*ia(1) - D_1*(phi1(1)-phi2(1)) ;
        J_2*w2(2) == D_1*(phi1(1)-phi2(1)) - D_2*(phi2(1)-phi3(1)) ;
               0 == D_2*(phi2(1)-phi3(1)) - B*phi3(2) + T_e ;
         phi1(2) == w1(1) ;
         phi2(2) == w2(1) ];
     
eqs = lhs(eqs) - rhs(eqs);
     
%% Assignment M01

[Am,Bm] = get_state_space(eqs, x(:,2), x(:,1), u);

% lstring = sprintf('\t%s \n\t=\n \t%s \n\t*\n \t%s \n\t+\n \t%s \n\t*\n \t%s', ...
%                   latex(x(:,2)), latex(Am), latex(x(:,1)), latex(Bm), latex(u));
% clipboard('copy',lstring)

x = x(2:end,:);

eqs_s = subs(eqs, ia(2), 0);
ia_v = solve(eqs_s(1), ia(1));
eqs_s = subs(eqs_s, ia(1), ia_v);
eqs_s = eqs_s(2:end);

[Am,Bm] = get_state_space(eqs_s, x(:,2), x(:,1), u);

% Case 1
x(:,1);
C1 = jacobian(x(:,1), [phi2(1);w2(1)]).';
D1 = 0;

% Case 2
[Ca,Da] = get_state_space( subs(eqs(1),ia(2),0), ia(1), x(:,1), u);
[Cb,Db] = get_state_space( eqs(4)              , phi3(2), x(:,1), u);
C2 = [Ca;Cb];
D2 = [Da;Db];



%% Assignment M02


% Question a)

vars_sym   = [K_e, K_t, J_1,   J_2,   B,    D_1, D_2]; 
vars_value = [.1, .1, 1e-5, 4e-5, 2e-3, 20,  2];
repl = @(x) subs(x, vars_sym, vars_value);

Av = repl(Am);
Bv = repl(Bm);
C1v = repl(C1);
C2v = repl(C2);

W_r = ctrb_s(Av,Bv);
W_o1 = obsv_s(Av,C1v);
W_o2 = obsv_s(Av,C2v);

rref(W_r)
rref(W_o1)
rref(W_o2)


% Question b)

svd([Av;C2])        % PBH test with lambda=0, results in rank<=4 (rank deficient)


% Question c)

vars_sym   = [R, K_e, K_t, J_1,   J_2,   B,    D_1, D_2]; 
vars_value = [1, .1, .1, 1e-5, 4e-5, 2e-3, 20,  2];
repl = @(x) double(subs(x, vars_sym, vars_value));

Av = repl(Am);
Bv = repl(Bm);
C1v = repl(C1);
C2v = repl(C2);

W_r = ctrb(Av,Bv);
W_o1 = obsv(Av,C1v);
W_o2 = obsv(Av,C2v);

svd(W_r)
svd(W_o1)
svd(W_o2)

kappa = @(W) max(svd(W))/min(svd(W));

kappa(W_r)
kappa(W_o1)
kappa(W_o2)

eig(Av)


svd([Av;C2v])        % PBH test with lambda=0, results in rank<=4 (rank deficient)

% PBH test
rank([Av, Bv])
rank([Av;C1v])
rank([Av;C2v])



% Question e)

Ts = 1e-3;

syms s t
matrix = inv(eye(size(Av))*s-Av);
exp_At = ilaplace(matrix);

Bd = double(int(exp_At,t,0,Ts)*Bv);
Ad=expm(Av*Ts);

sst = ss(Av,Bv,C1v,0);
c2d(sst,Ts,'zoh');


% Question f)

svd(ctrb(Ad,Bd))
svd(obsv(Ad,C1v))
svd(obsv(Ad,C2v))

abs(eig(Ad))



%% Assignment M03


% Question 1)

N = Bd;


% Question 2)

% Calculate noise variaces
sigma2_va   = calcVariance(0.3, 0.997);
sigma2_Te   = calcVariance(0.1, 0.997);
sigma2_phi2 = calcVariance(0.02, 0.997);
sigma2_w2   = calcVariance(0.01, 0.997);

R = diag([sigma2_va sigma2_Te sigma2_phi2 sigma2_w2])


% Question 3)

% Calculate Kalman gain
sysmodel = ss(Ad, [Bd N], C1v, 0, Ts)
Qm = diag([sigma2_va,sigma2_Te]);
Rm = diag([sigma2_phi2 sigma2_w2]);
[kest,L_kalman,P] = kalman(sysmodel, Qm , Rm)

% Check observer stability
eig(Ad-L_kalman*C1v)



%% Question 4)

x0 = [0 0 0 0 0]';

C1v_I = [1 1];
C1v_fb = C1v_I*C1v;

% Set extended system model
Ad_e = [Ad  0*Ad*C1v_fb';
       -C1v_fb eye(size(C1v_fb,1))];
Bd_e = [Bd;
        0*(C1v_fb*Bd)];
C1v_e= [C1v_fb 0*(C1v_fb*C1v_fb')];

svd_pos = svd(ctrb(Ad_e,Bd_e));
cond_number = max(svd_pos)/min(svd_pos)


sysmodel_e = ss(Ad_e,Bd_e,C1v_e, 0, Ts);

% Test if system is stabilizable
[Abar,Bbar,Cbar,T,k] = ctrbf(Ad_e,Bd_e,C1v_e);
Nctrb = sum(k);
eig(Abar(1:end-Nctrb,1:end-Nctrb));

% Calculate LQI controller
Qx = diag([1 1 10 0.001 1 0.001]);
Qu = diag([2 2])*1e6;
[Klqr,S,e] = dlqr(Ad_e, Bd_e, Qx, Qu);

Klqr

% Set variables for Simulink
K_fb  = Klqr(:,1:length(Ad));
K_I = Klqr(:,length(Ad)+1:end)*20;
% K_I(:,1) = 0;



%% Help functions

function sigma2_val = calcVariance(ub, realization_ub)
    syms sigma2 x
    eq = realization_ub == int(1/sqrt(2*pi*sigma2)*exp(-(x^2)/(2*sigma2)),x,-ub,ub);
    sigma2_val = double(solve(eq,sigma2));
    mu = 0;
    pd = makedist('Normal',mu,abs( sqrt(sigma2_val)));
    v_aw = linspace(-ub*2,ub*2,100);
    y = pdf(pd,v_aw);
%     figure('Color','white');
%     plot(v_aw,y);
end

function mctrb = ctrb_s(A,B)
    mctrb = [];
    for i=1:size(A,2)
        mctrb = [mctrb A^(i-1)*B];
    end
end

function mobsv = obsv_s(A,C)
    mobsv = [];
    for i=1:size(A,2)
        mobsv = [mobsv; C*A^(i-1)];
    end
end

function [A,b] = get_state_space(eqs, xdot, x, u)
    A = -jacobian(eqs, xdot) \ jacobian(eqs, x);
    b = -jacobian(eqs, xdot) \ jacobian(eqs, u);
end