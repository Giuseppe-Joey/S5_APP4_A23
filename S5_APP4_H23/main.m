%%  S5 - APP4 - PROBLEMATIQUE - MAIN.M
%   Auteur:     Giuseppe Lomonaco
%   CIP:        LOMG2301
%   Auteur:     Lucas Corrales
%   CIP:        CORL0701

%   Date:       2-MARS-2023
%   Modifications (Date - initiales - détails):


% JAMAIS UN BODE AVEC UNE FTBF!!!!!!!

%% DEBUT DE La problematique
clc
close all
clear all

% Initialisation
constantes % call le fichier des constantes




%% QUESTION a)

disp(['------QUESTION a)------']);

% Calcul des valeurs propres du systeme
val_propres_A = eig(A)

% calcul des carac temporelles
wn = abs(val_propres_A);
zeta = -real(val_propres_A)./wn;
wa = wn.*sqrt(1-zeta.^2);
phi = acos(zeta);
Mp = 100*exp(-pi./tan(phi));
ts = 4./(zeta.*wn);
tp = pi./wa;

% affichage des carac temporelles
disp(["Affichage des carac temporelles:"]);
disp(['wn = ', num2str(wn(end)), ' rad/s']);
disp(['zeta = ', num2str(zeta(end)), ' unites']);
disp(['wa = ', num2str(wa(end)), ' rad/s']);
disp(['phi = ', num2str(phi(end)), ' radian']);
disp(['Mp = ', num2str(Mp(end)), ' %']);
disp(['ts = ', num2str(ts(end)), ' s']);
disp(['tp = ', num2str(tp(end)), ' s']);






%% QUESTION b) - v/a_prop
disp(['------QUESTION a)------']);
[num,den] = ss2tf(A,B,C,D,2);   % le 2 signifie quon veut le 2e element de U soit aprop
v_aprop = tf(num(1, :), den)         % v/aprop
alpha_aprop = tf(num(2, :), den);         % alpha/aprop
teta_aprop = tf(num(3, :), den);         % teta / aprop
q_aprop = tf(num(4, :), den);         % q/aprop
gamma_aprop = tf(num(5, :), den);         % gamma/aprop


numerateur = num(1,:);
zeros = roots(num(1,:))
poles = roots(den)


p1 = poles(1)
p2 = poles(2)
p3 = poles(3)
p4 = poles(4)


z1 = zeros(1);
z2 = zeros(2);
z3 = zeros(3);


figure('Name', 'Rlocus de V/a_prop')
rlocus(v_aprop)

figure('Name', 'Step de V/a_prop')
step(v_aprop)


[num, den] = zp2tf(zeros, poles, 1);
TF_a_partir_zeros_poles_v_aprop = tf(num,den)
 

% section lucas
% identification de la fonction de transfert à partir des pôles et zéros
z = zeros;
p = poles;
k = 1;
sys_zpk = zpk(z, p, k)







%% QUESTION b) - phase non minimale
disp(['------QUESTION b)------']);
[num,den] = ss2tf(A,B,C,D,1);   % le 1 signifie quon veut le 1er element de U soit delta_c
v_delta_c = tf(num(1, :), den)         % v/delta_c

figure('Name', 'V sur Delta_C')
rlocus(v_delta_c)

figure('Name', 'V sur Delta_C')
bode(v_delta_c), grid




%% QUESTION c)
disp(['------QUESTION c)------']);








%% FEEDBACK d)
disp(['------QUESTION d)------']);

% EXEMPLE FEEDBACK
%num = [1 2];
%den = [1 3 2];
%G = tf(num, den)
%Kv = 0.5;
%H = tf(Kv, 1);
%T = feedback(G*H, 1);
%step(T)

%% question e)
disp(['------QUESTION e)------']);

G = v_aprop
Kv = 1.0262;

    figure('Name', 'Avec feedback')
    T = feedback(G, Kv);
    step(T)
    figure('Name', 'Sans feedback')
    step(G)
%     rlocus(G)
    legend
T
[num_T,den_T] = tfdata(T)
[A1,B1,C1,D1] = tf2ss(num_T{1},den_T{1})




%% question f)
disp(['------QUESTION f)------']);

figure('Name', 'Bode')
for i = 1:length(Kv)  
    bode(Kv(i)*G), grid
    margin(Kv(i)*G);
    hold on
end

% fprintf('Phase margin: %0.2f degrees\n', Pm);
% fprintf('Gain margin: %0.2f dB\n', Gm);
% fprintf('Omega G: %0.2f rad/s\n', Wog);
% fprintf('Omega P: %0.2f rad/s\n', wop);







%% question g)
disp(['------QUESTION g)------']);

%trouver les pôles dominant (les pôle qui affect le plus la fonction)
[R,P,K]=residue(G.numerator{1},G.denominator{1});

%trouver le poid des pôles
Cdom=abs(R)./abs(real(P))
 
%le poid pole avec la plus grande valeur=plus proche de l'axe imaginaire
%on veut toujours deux pole... systeme d'ordre 2

%reduction du systeme
[num,den]=residue(R(3:4),P(3:4),K);
 
TFR=tf(num,den);

g0 = dcgain(G);
g1 = dcgain(TFR);

G_simpl = (g0/g1)*TFR

[num_G_simpl,den_G_simpl] = tfdata(G_simpl)

Wn_G_simpl = sqrt(den_G_simpl{1}(3))
Zeta = 1

Kv1 = ((Wn_G_simpl^2)-0.047)/1.465
Kv2 = ((den_G_simpl{1}(2))-0.028)/2.122
 
figure('Name', 'comparaison avec nouvelle tf')
rlocus(G_simpl)
% step(G)
% hold on
% step(G_simpl)






%% question h)
disp(['------QUESTION h)------']);
Kv = 1.306;
FFFT = tf([2.122*Kv 1.465*Kv],[1 (2.122*Kv+0.028) (1.465*Kv+0.047)])

rlocus(T)
hold on
rlocus(FFFT)






%% question i)
disp(['------QUESTION i)------']);

FT_originale = v_aprop
FT_reduite = G_simpl


figure
rlocus(FT_reduite), grid
hold on
rlocus(FT_originale), grid
legend







%% question j)
disp(['------QUESTION j)------']);

[num,den] = ss2tf(A1,B1,C1,D1,1);   % le 1 signifie quon veut le 1e element de U (soit delta_c)
gamma_delta = tf(num(5), den)         % gamma/delta

figure('Name', 'Rlocus de gamma/delta');
rlocus(gamma_delta);

% figure('Name', 'Step de gamma/delta');
% step(gamma_delta);

% figure('Name', 'Bode de gamma/delta');
% bode(gamma_delta);









%% question k)
disp(['------QUESTION k)------']);

[num,den] = ss2tf(A1,B1,C1,D1,1);   % le 1 signifie quon veut le 1e element de U (soit delta_c)
gamma_delta = tf(num(5), den)         % gamma/delta

G1 = gamma_delta
%Kp = 1;

%13.9 = 20*log(Kp)

Kp = 10^(9/20)


    figure('Name', 'Kp')
    %T2 = feedback(G1*Kp, 1);
    %step(G1)
    bode(G1*Kp)
    margin(G1*Kp)
    legend


    [Gm,Pm,Wcg,Wcp] = margin(G1*Kp)
    Gm_dB = 20*log10(Gm)
    
% Calculer l'erreur en régime permanent
FT = G1*Kp
[num_FT,den_FT] = tfdata(FT);
% puisque s tend vers 0
Kpos = num_FT{1}(end)/den_FT{1}(end)
erp = 1/(1+Kpos)




disp(['']);
disp(['']);







%% question l)
disp(['------QUESTION l)------']);

% Obtain the magnitude response data
[mag] = bode(FT);
mag = squeeze(mag);

Kpos_bode = mag(1)
erp_bode = 1/(1+Kpos_bode)
disp(['']);
disp(['']);






%% question m)
disp(['------QUESTION m)------']);

    figure('Name', 'G1*Kp')
    T1 = feedback(G1*Kp, 1);
    step(T1)
    hold on
    step(T)
    legend


    FT_T1 = tf(T1)

[num_T1,den_T1] = tfdata(T1)
[A2,B2,C2,D2] = tf2ss(num_T1{1},den_T1{1})



disp(['']);
disp(['']);








%% question n)
disp(['------QUESTION n)------']);





disp(['']);
disp(['']);


