%%  S5 - APP4 - PROBLÉMATIQUE
%   Auteur:     Giuseppe Lomonaco
%   CIP:        LOMG2301
%   Auteur:     Isabelle Handfield
%   CIP:        HANI1401

%   Date de creation:                       10-Octobre-2023
%   Date de derniere modification:          16-Octobre-2023

%   DESCRIPTION: PROBLÉMATIQUE


clc
close all
clear all

% Chargement du fichier Annexe A
fprintf("*** Chargement du fichier 'Annexe_A' ***\n");
Annexe_A


%% a) ANALYSE DES CARACTÉRISTIQUES DYNAMIQUES DE L'AVION
%     À PÂRTIR DES VALEURS PROPRES DU SYSTEME ET DE LA RÉPONSE TEMPORELLE
%     DE LA VITESSE v DU SYSTEME SOUMIS A UN ÉCHELON SUR a_prop
%     DONNER, POUR CHAQUE MODE:         -le temps du premier pic;
%                                       -le dépassement maximum;
%                                       -le temps de stabilisation;
%                                       -le facteur d'amortissement;
%                                       -la période des oscillations amorties
%                                           et naturelles et vérifier ces
%                                           résultats à partir de la réponse
%                                           temporelle.

disp(" ")
disp("*** Question a) ***")

% calcul des valeurs propres du systeme 
val_propes = eig(A)
disp("On remarque qu'il y a deux paires de poles conj. complexes alors 2 modes dynamiques")

% on trouve la fonction de transfert pour le 2e element de la matrice U
% (soit aprop)
[num, den] = ss2tf(A, B, C, D, 2);      
v_sur_aprop = tf(num(1,:), den)

zeros = roots(num(1,:));
poles = roots(den);

% figure('Name', 'Question a)')
% pzmap(zeros, poles)

% figure('Name', 'Question a)')
% rlocus(v_sur_aprop)

% reduction de la FT avec residue
[R, P, K] = residue(num(1,:), den);
poid = abs(R) ./ abs(real(P));
disp('On remarque que les poles dominants sont les 3e et 4e, mode dominant')
disp(' ')
disp(' ')


% MODE 1
disp("*** MODE 1 (mode dominant) ***")

[numR, denR] = residue(R(3:4), P(3:4), K);
TF_temporaire = tf(numR, denR);

% ajustement du gain DC
gain0 = dcgain(v_sur_aprop);
gainR = dcgain(TF_temporaire);

numR = numR * (gain0/gainR);
TF_reduce_mode1 = tf(numR, denR)

% on trouve les zeros et les poles du mode dominant
zeros = roots(numR);
poles = roots(denR);

% figure('Name', 'Question a)')
% pzmap(zeros, poles)

% affichage de la reponse temporelle MODE 1
figure('Name', 'Question a)')
step(TF_reduce_mode1)
title("Réponse temporelle MODE 1")
grid on

% selon la formule du     wa = wn* sqrt(1-zeta^2)
wa = imag(poles);
wa = wa(1);

wn = abs(poles);
wn = wn(1);

zeta = -real(poles) / wn;
zeta = zeta(1);

phi = acos(zeta);
phi_degres = acosd(zeta);

ts = 4 / (zeta * wn);
tp = pi / wa;
Mp = 100 * exp(-pi/tan(phi));


% affichage des caracartéristiques temporelles
disp(["--------------------------------------------------------------------------"]);
disp(["Caracartéristiques temporelles *MODE 1* (poles qui ont le plus GRAND poid)"]);
disp(["--------------------------------------------------------------------------"])
disp(['wn   = ', num2str(wn), ' rad/s']);
disp(['zeta = ', num2str(zeta), ' unites']);
disp(['wa   = ', num2str(wa), ' rad/s']);
disp(['phi  = ', num2str(phi), ' radian']);
disp(['phi  = ', num2str(phi_degres), ' degrés']);
disp(['Mp   = ', num2str(Mp), ' %']);
disp(['ts   = ', num2str(ts), ' s']);
disp(['tp   = ', num2str(tp), ' s']);
disp(["--------------------------------------------------------------------------"])

fprintf("\n\n")
% *************************************************************************************







% MODE 2
disp("*** MODE 2 ***")
% on trouve les poles 
[num, den] = ss2tf(A, B, C, D, 2);
[R, P, K] = residue(num(1,:), den);
poid = abs(R) ./ abs(real(P))
disp('On garde maintenant les 1er et 2e poles pour faire le mode 2')
[numR, denR] = residue(R(1:2), P(1:2), K);

% correction du gain DC
gain0 = dcgain(num(1,:), den);
gainR = dcgain(numR, denR);
numR = numR * (gain0/gainR);
TF_reduce_mode2 = tf(numR, denR)

% affichage de la reponse temporelle MODE 2
figure('Name', 'Question a)')
step(TF_reduce_mode2)
title("Réponse temporelle MODE 2")
grid on

% on trouve les caracteristiques temporelles
poles = roots(denR);

wa = imag(poles);
wa = wa(1);

wn = abs(poles);
wn = wn(1);

zeta = -real(poles) / wn;
zeta = zeta(1);

phi = acos(zeta);
phi_degres = acosd(zeta);

Mp = 100 * exp(-pi / tan(phi));
ts = 4 / (zeta*wn);
tp = pi / wa;






% affichage des caracartéristiques temporelles
disp(["--------------------------------------------------------------------------"]);
disp(["Caracartéristiques temporelles *MODE 2* (poles qui ont le plus PETIT poid)"]);
disp(["--------------------------------------------------------------------------"])
disp(['wn   = ', num2str(wn), ' rad/s']);
disp(['zeta = ', num2str(zeta), ' unites']);
disp(['wa   = ', num2str(wa), ' rad/s']);
disp(['phi  = ', num2str(phi), ' radian']);
disp(['phi  = ', num2str(phi_degres), ' degrés']);
disp(['Mp   = ', num2str(Mp), ' %']);
disp(['ts   = ', num2str(ts), ' s']);
disp(['tp   = ', num2str(tp), ' s']);
disp(["--------------------------------------------------------------------------"]);


% on affiche quelques lignes vides pour separer de la prochaine section
fprintf("\n\n\n")














%% b) IDENTIFICATION DE LA FONCTION DE TRANSFERT À PHASE NON-MINIMALE
%       À PARTIR DES PÔLES ET DES ZÉROS

disp('*******************************************************************')
disp("*** Question b) ***")

disp("on desire gamma sur delta_c")

% on veut la premiere entree soit delta c
[num, den] = ss2tf(A, B, C, D, 1)
gamma_sur_deltaC = tf(num(5,:), den)        % on utilise la 5 ligne qui correspond a gamma

% figure('Name', 'Question b)')
% rlocus(gamma_sur_deltaC)
% title("Gamma sur DeltaC")

figure('Name', 'Question b)')
pzmap(gamma_sur_deltaC)
title("Gamma sur DeltaC")
axis([-0.05  0.05  -1    1   ])
disp("on verifie quil y a bel et bien un zero a droite de laxe imaginauire")




% % on utilise le 2e element de lentree soit a_prop
% [num, den] = ss2tf(A, B, C, D, 2);   
% v_sur_aprop = tf(num(1,:), den);
% 
% zeros = roots(num(1,:));
% poles = roots(den);
% 
% n = length(poles);
% m = length(zeros);
% fprintf("n = %d et m = %d \n", n, m)
% fprintf("Il y a '%d' asymptote(s) a l'infini", n-m)
% 
% [num, den] = zp2tf(zeros, poles, 1);
% v_sur_aprop_avec_zp2tf = tf(num,den);
% 
% figure('Name', 'Question b)')
% rlocus(v_sur_aprop_avec_zp2tf)
% title("FT a partir de zp2tf()")
% 
% % % identification de la fonction de transfert à partir des pôles et zéros
% % systeme_zpk = zpk(zeros, poles, 1)
% % 
% % figure('Name', 'Question b)')
% % rlocus(systeme_zpk)
% % title("FT a partir de zpk()")







fprintf("\n\n\n")







%% c) PRÉSENTATION DU LIEU DES RACINES FAIT À LA MAIN ENTRE a_prop ET v 
%       (AVEC ÉTAPES) ET VALIDATION MATLAB

disp('*******************************************************************')
disp("*** Question c) ***")
disp('**** dessiner a la main le lieu des racines et verification avec MATLAB ****')

figure('Name', 'Question c)')
rlocus(v_sur_aprop_avec_zp2tf)
title("FT a partir de zp2tf()")


fprintf("\n\n\n")









%% d) A PARTIR DU LIEU DE RACINES, EXPLICATION DE L'EFFET DE LA RÉTROACTION
%      Kv SUR LA STABILITÉ, LE TEMPS DE STABilisation ET LE DEPASSement MAX
 
% point de separation...demander a rlocus de retourner la valeur au gain
% -1.39
% [r, k] = rlocus(sys)
% r = rlocus(sys, k)
disp('*******************************************************************')
disp("*** Question d) ***")
disp('on fait un zoom sur le lieu des racines pour trouver Kv')

[num, den] = ss2tf(A, B, C, D, 2);   
v_sur_aprop = tf(num(1,:), den);

figure('Name', 'Question d)')
rlocus(v_sur_aprop)
axis([-1.5    0.5   -1  1])

Kv = 1.03;
disp(['On trouve sur le graph le gain max Kv = ', num2str(Kv), ' au point dintersection'])
fprintf("\n\n\n")
                                                                                                                                                                                            









%% e) CONCEPTION DE LA BOUCLE INTERNE À PARTIR DU LIEU DE RACINES: 
%       1) CHOIX DE Kv
%       2) CALCUL DU NOUVEAU MODELE n1(s)/d1(s) ET DU MODELE A1, B1, C1, D1
%       INCLUANT LA BOUCLE

disp('*******************************************************************')
disp("*** Question e) ***")

Kv = 1.03;

B1 = B(:,1);
B2 = B(:,2);

C1 = C(1,:);
C2 = C(2,:);
C3 = C(3,:);
C4 = C(4,:);
C5 = C(5,:);   

D5 = D(5,1);

%Nouvelles matrices incluant leffet de la boucle interne (voir prob 5 procedural 1)
A1 = A - B2*Kv*C1;  % on garde A - B2*Kv*C1 car A+B... ne donne pas de resultats coherents
% Aa = A + B2*Kv*C1;
B1 = B1;
C1 = C5;
D1 = D5;


disp(' ')
disp('*** La FT suivante a ete validee par les enseignants et est bonne')
[num, den] = ss2tf(A1, B1, C1, D1);
TF = tf(num, den)

figure('Name', 'Question e)')
rlocus(TF)

fprintf("\n\n\n")








%% f) VÉRIFICATION DES MARGES AVEC DIAGRAMMES DE BODE ET Kv;
%       COMMENTAIRE SUR LE SENS DES MARGES, UTILITÉ

disp('*******************************************************************')
disp("*** Question f) ***")

figure('Name', 'Question f)')
margin(TF)
grid on


figure('Name', 'Question f)')
step(TF)
grid on

fprintf("\n\n\n")











%% g) réduction de la FT entre a_prop et v, présentation de la nouvelle FT
%     réduite, fermeture de la boucle à la main sur cette FT, explication
%     de leffet de Kv sur les param standards (K, zeta, wn, tau) et la reponse

disp('*******************************************************************')
disp("*** Question g) ***")

[R, P, K] = residue(num, den);
poid = abs(R) ./ abs(real(P))
[numR, denR] = residue(R(3:4), P(3:4), K);

tf_tempo = tf(numR, denR);

% ajustement du gain dc
gain0 = dcgain(num, den);
gainR = dcgain(numR, denR);
numR = numR * (gain0/gainR);
TF_reduce = tf(numR, denR)


figure('Name', 'Question g)')
step(TF_reduce)
grid on


zeros = roots(numR);
poles = roots(denR);

wn = abs(poles);
wn = wn(1);

wa = imag(poles);
wa = wa(1);

zeta = -real(poles) / wn;
zeta = zeta(1);

phi = acos(zeta);

Mp = 100 * exp(-pi / (tan(phi)));
ts = 4 / (zeta*wn);
tp = pi / (wn*sqrt(1-(zeta^2)));



figure('Name', 'Question g)')

fprintf("\n\n\n")









%% h) calcul du Kv avec le système réduit et comparaison avec la valeur
%     trouvée ci-dessus en (e)

disp('*******************************************************************')
disp("*** Question h) ***")

Kv = 1.03;

FTBF = feedback(v_sur_aprop, Kv);
TF = tf([2.122*Kv    1.465*Kv],[1   (2.122*Kv+0.028)    (1.465*Kv+0.047)])

figure('Name', 'Question h)')
rlocus(FTBF)
hold on
rlocus(TF)
legend('FTBF', 'TF')

fprintf("\n\n\n")






%% i) comparaison des lieux des racines original et réduit de la FT entre
%     a_prop et v sur le même graphique

disp('*******************************************************************')
disp("*** Question i) ***")

[num, den] = ss2tf(A, B, C, D, 2);
v_sur_aprop = tf(num(1,:), den)


[R, P, K] = residue(num(1,:), den);
poid = abs(R) ./ abs(real(P))
[numR, denR] = residue(R(3:4), P(3:4), K);
tf_tempo = tf(numR, denR);

% ajustement du gain DC
gain0 = dcgain(v_sur_aprop);
gainR = dcgain(tf_tempo);
numR = numR * (gain0/gainR);

tf_reduite = tf(numR, denR)


figure('Name', 'Question i)')
rlocus(v_sur_aprop)
hold on
rlocus(tf_reduite)
title('Comparaison entre FT originale et FT reduite')
legend('FT v sur a_prop originale', 'FT reduite')

fprintf("\n\n\n")







%% j) lieu des racines entre gamma_d et gamma incluant la boucle interne
%     (i.e. n1/d1), effet de Kp sur la réponse

disp('*******************************************************************')
disp("*** Question j) ***")

[num,den] = ss2tf(A1,B1,C1,D1, 1);   % le 1 signifie quon veut le 1e element de U (soit delta_c)
gamma_sur_deltaC = tf(num(5), den)         % on veut 

figure('Name', 'Question j)');
rlocus(gamma_sur_deltaC);
title("Rlocus de gamma sur deltaC")

fprintf("\n\n\n")







%% k) design de la boucle externe: avec Bode, 
% (1) calculer Kp limite pour instabilité;
% (2) calculer Kp pour marges de 6 dB et 30 degrés;
% (3) mesurer l'erreur en régime permanent par simulation sur MATLAB

disp('*******************************************************************')
disp("*** Question k) ***")

[num,den] = ss2tf(A1,B1,C1,D1,1);   % le 1 signifie quon veut le 1e element de U (soit delta_c)
gamma_sur_delta = tf(num(5), den)         % gamma/delta

% Kp donne a 6 dB
Kp = 10^(6/20)


figure('Name', 'Question k)')
margin(gamma_sur_delta * Kp)
grid on


[GM,PM,wp,wg] = margin(gamma_sur_delta * Kp)
GM_dB = 20*log10(GM)

% Calcule de lerreur en régime permanent
FT = gamma_sur_delta * Kp
[num,den] = tfdata(FT);

% puisque s tend vers 0
K_pos = num{1}(end)/den{1}(end)
e_RP = 1/(1+K_pos)

fprintf("\n\n\n")








%% l) calcul de l'erreur en régime permanent à partir du gain statique tel
% que lu sur le diagramme de Bode, comparaison avec (3) de l’item ci-dessus


disp('*******************************************************************')
disp("*** Question l) ***")

% on obtient la magnitude
[mag] = bode(gamma_sur_delta * Kp);
mag = squeeze(mag);

Kpos_bode = mag(1)
erp_bode = 1/(1+Kpos_bode)

fprintf("\n\n\n")









%% m) calcul du nouveau modèle n2(s) / d2(s) et le modèle A2, B2, C2, D2

disp('*******************************************************************')
disp("*** Question m) ***")

FTBF = feedback(gamma_sur_delta * Kp, 1)

figure('Name', 'Question m)')
step(gamma_sur_delta)
hold on
step(FTBF)
legend('FTBO Gamma sur delta', 'FTBF Gamma sur delta')


[num,den] = tfdata(FTBF)
[A2,B2,C2,D2] = tf2ss(num{1},den{1})

fprintf("\n\n\n")










%% n) comparaison et discussion des réponses temporelles avec
%     compensateurs P, PD, PI et PID

disp('*******************************************************************')
disp("*** Question n) ***")

P = 0;

PD = 0;

PI = 0;

PID = 0;

fprintf("\n\n\n")












