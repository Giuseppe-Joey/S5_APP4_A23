% Resolution de la problematique
% APP4 S5 GE

clc
clear
close all

% Matrice ABCD modele lineaire
A = [-0.018223 -0.088571 -9.78 0;
     -0.003038 -1.2563 0 1;
      0 0 0 1;
      0.0617 -28.075 0 -4.5937];
 
 B = [0 1.1962;
      0 -0.00120;
      0 0;
      7.84 -4.05];
  
  C = [1 0 0 0;
       0 57.296 0 0;
       0 0 57.296 0;
       0 0 0 57.296;
       0 -57.296 57.296 0];
   
   D = [0 0;
        0 0;
        0 0;
        0 0;
        0 0];

% fonctions de transfert du systeme
% 2 entrees -> 5 sorties DONC 10 FT
[num1, den] = ss2tf(A, B, C, D, 1);
[num2, den] = ss2tf(A, B, C, D, 2);
FTBO_0(1) = tf(num1(1,:), den);
FTBO_0(2) = tf(num1(2,:), den);
FTBO_0(3) = tf(num1(3,:), den);
FTBO_0(4) = tf(num1(4,:), den);

% FT pour angle profondeur et angle vol
% -> phase non-minimale <-
FTBO_0(5) = tf(num1(5,:), den);

% FT pour vitesse v selon propulsion aprop
FTBO_0(6) = tf(num2(1,:), den)

FTBO_0(7) = tf(num2(2,:), den);
FTBO_0(8) = tf(num2(3,:), den);
FTBO_0(9) = tf(num2(4,:), den);
FTBO_0(10) = tf(num2(5,:), den);

% trouve les valeurs propres du systeme
e = eig(A);

% voir feuille equation exam APP3
wd_0 = imag(e);
alpha_0 = real(e);

theta_0 = atan(wd_0./alpha_0);
zeta_0 = cos(theta_0);

wn_0 = abs(alpha_0./zeta_0);
ts_0 = 4./(zeta_0.*wn_0);
t_max_0 = abs(pi./wd_0);
M_max_0 = 100*exp(-pi./tan(abs(theta_0)));
ymax = exp(-(pi .* zeta_0)./(sqrt(1-zeta_0.^2))) + 1;

t = [0:0.1:500]';
u = ones(size(t));
y6 = lsim(FTBO_0(6), u, t);
y5 = lsim(FTBO_0(5), u, t);

figure(1)
hold on
step(FTBO_0(9))
stepinfo(FTBO_0(9))
hold off

figure(2)
hold on
step(FTBO_0(10))
stepinfo(FTBO_0(10))
hold off

% trouver les poles et les zeros
[p, z] = pzmap(FTBO_0(6));
tot = sum(p) - sum(z);

figure(3)
rlocus(FTBO_0(6))

c7 = [1.196 13.994 241.81 1138.95 3587.17 3334.02 -67.57];
roots(c7);

% montrer l'effet du gain Kv
t2 = [0:0.01:50]';
u2 = ones(size(t2));
boucle_Kv1 = feedback(FTBO_0(6), 0.5)
boucle1 = lsim(boucle_Kv1, u2, t2);

boucle_Kv2 = feedback(FTBO_0(6), 1.03)
boucle2 = lsim(boucle_Kv2, u2, t2);

figure(4)
hold on
rlocus(FTBO_0(6), 0.5)
rlocus(FTBO_0(6), 1.33)
legend('Gain Kv = 0.5','Gain Kv = 1.03')
hold off

Kv = 1.03;

% % boucle interne %
% B1 = B(:,1);
% C1 = C(1,:);
% B2 = B(:,2);
% Aa = A-B2*Kv*C1;
% Ba = [B1 B2];
% Ca = C;
% Da = D;
% [numa,dena] = ss2tf(Aa,Ba,Ca,Da,2);
% tf6a = tf(numa(1,:),dena)

% boucle interne %
B1 = B(:,1);
C1 = C(1,:);
B2 = B(:,2);
kv = 1.03;
Aa = A-B2*kv*C1
Ba = B1
Ca = C(5,:)
Da = zeros(1,1)
[numa,dena] = ss2tf(Aa,Ba,Ca,Da);
tf6a = tf(numa(1,:),dena)

% verification des marges
figure(5)
hold on
bode(FTBO_0(6)*Kv)
margin(num2(1,:), den)
hold off
[Gm,Pm,Wcg,Wcp] = margin(num2(1,:), den)

% reduction de la FT
[r, p, k] = residue(num2(1,:), den);

ratio = abs(r./real(p));

[numr, denr] = residue(r(3:4), p(3:4), k);
g0 = dcgain(num2(1,:), den);
gr = dcgain(numr, denr);
numr = numr * g0/gr;

FTBO_1R = tf(numr, denr)

% Lieu des racines en Gamma_d et Gamma
figure(6)
rlocus(FTBO_0(5) * tf6a)