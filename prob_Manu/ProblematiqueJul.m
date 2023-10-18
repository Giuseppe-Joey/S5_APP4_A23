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
[num1, den1] = ss2tf(A, B, C, D, 1);
[num2, den2] = ss2tf(A, B, C, D, 2);
FTBO_0(1) = tf(num1(1,:), den1);
FTBO_0(2) = tf(num1(2,:), den1);
FTBO_0(3) = tf(num1(3,:), den1);
FTBO_0(4) = tf(num1(4,:), den1);

% FT pour angle profondeur et angle vol
% -> phase non-minimale <-
FTBO_0(5) = tf(num1(5,:), den1);

% FT pour vitesse v selon propulsion aprop
FTBO_0(6) = tf(num2(1,:), den2)

FTBO_0(7) = tf(num2(2,:), den2);
FTBO_0(8) = tf(num2(3,:), den2);
FTBO_0(9) = tf(num2(4,:), den2);
FTBO_0(10) = tf(num2(5,:), den2);



% trouve les valeurs propres du systeme
e = eig(A)

% voir feuille equation exam APP3
wd_0 = imag(e)
alpha_0 = real(e)

theta_0 = atan(wd_0./alpha_0)
zeta_0 = cos(theta_0)

wn_0 = abs(alpha_0./zeta_0)
ts_0 = 4./(zeta_0.*wn_0)
t_max_0 = abs(pi./wd_0)
M_max_0 = 100*exp(-pi./tan(abs(theta_0)))
ymax = exp(-(pi .* zeta_0)./(sqrt(1-zeta_0.^2))) + 1

t = [0:0.1:500]';
u = ones(size(t));
y6 = lsim(FTBO_0(6), u, t);
y5 = lsim(FTBO_0(5), u, t);

boucle_Kv = feedback(FTBO_0(6), 1)

figure(1)
hold on
plot(t, y6)

% figure(2)
% hold on
% plot(t, y5)

% trouver les poles et les zeros
[p, z] = pzmap(FTBO_0(6))
tot = sum(p) - sum(z)

figure(3)
rlocus(FTBO_0(6))


s = roots([-1.196 -14.03 -241.81 -1138.95 -3587.17 -3333.95 68.03])  


% boucle interne %
B1 = B(:,1);
C1 = C(1,:);
B2 = B(:,2);
kv = 1.03;
Aa = A-B2*kv*C1;
Ba = B1;
Ca = C(5,:);
Da = zeros(1,1);
[numa,dena] = ss2tf(Aa,Ba,Ca,Da);
tf6a = tf(numa,dena)


% verification des marges
figure(5)
hold on
bode(tf6a)
margin(numa, dena)
hold off
[Gm,Pm,Wcg,Wcp] = margin(numa, dena)


% reduction de la FT
pzmap(tf6a)
[r, p, k] = residue(num2(1,:), den2);

ratio = abs(r./real(p))

[numr, denr] = residue(r(3:4), p(3:4), k);
g0 = dcgain(num2(1,:), den2);
gr = dcgain(numr, denr);
numr = numr * g0/gr;

FTBO_1R = tf(numr, denr)


%%Calcul du Kv
y6 = lsim(FTBO_1R, u, t);
figure(4)
hold on
plot(t, y6)

P_inter = roots([1.061^2 -1.465 (0.01406^2-0.04719)])

%rlocus(FTBO) -> rlocus(10*FTBO)
    % les stations changent pour facteur 10
% bode ftbo négale pas bode(10*ftbo)
    % phase reste pareil ligne de gain va changer (10* +haut)
figure(10)
bode(FTBO_1R)
hold on
bode(1.03*FTBO_1R)
% [Gm,Pm,Wcg,Wcp] = margin(FTBO_1R)
% [Gm,Pm,Wcg,Wcp] = margin(1.03*FTBO_1R)


%Comparaison lieux des racines
figure(3)
rlocus(FTBO_0(6))
hold on
rlocus(FTBO_1R)
legend('FTBO', 'FTBO reduite')

%Lieu des racine delata_d et delta
figure(6)
rlocus(FTBO_0(5)*tf6a)


%Design boucle interne

figure(7)
hold on
Gs = tf6a;
bode(Gs)
[num1, den1] = tfdata(Gs, 'v');

Kpos = num1(end)/den1(end)

margin(num1, den1)
hold off
[Gm,Pm,Wcg,Wcp] = margin(num1, den1)

Klim = 20*log10(Gm);
kp = 10^((Klim-6)/20);

%Modèle A2B2C2D2

C5 = C(5,:);
A2 = Aa - (Ba*kp*C5);
B2 = Ba*kp;
C2 = C5;
D2 = 0;

[num2, den2] = ss2tf(A2,B2,C2,D2);

tf6b = tf(num2,den2)


%Compensateur

PD = tf([kp kp],1);
PI = tf([kp kp],[1 0]);
PID = tf([kp kp kp],[1 0]);

FTBO_P = tf6b
FTBO_PD = series(PD, tf6a);
FTBO_PI = series(PI, tf6a);
FTBO_PID = series(PID, tf6a);

FTBF_P = feedback(FTBO_P,1)
FTBF_PD = feedback(FTBO_PD,1)
FTBF_PI = feedback(FTBO_PI,1)
FTBF_PID = feedback(FTBO_PID,1)

figure(8)
subplot(2,2,1)
step(FTBF_P)

subplot(2,2,2)
step(FTBF_PD)

subplot(2,2,3)
step(FTBF_PI)

subplot(2,2,4)
step(FTBF_PID)