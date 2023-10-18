%%  S5 - APP4 - PROBLÉMATIQUE
%   Auteur:     Giuseppe Lomonaco
%   CIP:        LOMG2301
%   Auteur:     Isabelle Handfield
%   CIP:        HANI1401

%   Date de creation:                       10-Octobre-2023
%   Date de derniere modification:          18-Octobre-2023

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
poles = val_propes;

% reduction de la FT avec residue
[R, P, K] = residue(num(1,:), den);
poid = abs(R) ./ abs(real(P));
disp('On remarque que les poles dominants sont les 3e et 4e, mode dominant')
disp(' ')
disp(' ')


% MODE 1
disp("*** MODE 1 - PHUGOIDE ***")

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

% affichage de la reponse temporelle MODE 1
figure('Name', 'Question a)')
step(TF_reduce_mode1)
title("Réponse temporelle MODE 1 - PHUGOIDE")
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

disp("on verifie quil y a bel et bien un zero a droite de laxe imaginauire avec pzmap")

figure('Name', 'Question b)')
pzmap(gamma_sur_deltaC)
title("Gamma sur DeltaC")
axis([-0.05  0.05  -1    1   ])

figure
rlocus(gamma_sur_deltaC)


fprintf("\n\n\n")







%% c) PRÉSENTATION DU LIEU DES RACINES FAIT À LA MAIN ENTRE a_prop ET v 
%       (AVEC ÉTAPES) ET VALIDATION MATLAB

disp('*******************************************************************')
disp("*** Question c) ***")
disp('**** dessiner a la main le lieu des racines et verification avec MATLAB ****')

figure('Name', 'Question c)')
rlocus(v_sur_aprop)
title("FT v sur aprop")
axis([-4    1     -8       8])


fprintf("\n\n\n")









%% d) A PARTIR DU LIEU DE RACINES, EXPLICATION DE L'EFFET DE LA RÉTROACTION
%      Kv SUR LA STABILITÉ, LE TEMPS DE STABilisation ET LE DEPASSement MAX
 
disp('*******************************************************************')
disp("*** Question d) ***")
disp('on fait un zoom sur le lieu des racines pour trouver Kv')

figure('Name', 'Question d)')
rlocus(v_sur_aprop)
axis([-1.5    0.5   -1  1])

Kv = 1.03;
disp(['On trouve sur le graph le gain Kv = ', num2str(Kv), ' au point dintersection'])
fprintf("\n\n\n")
                                                                                                                                                                                            









%% e) CONCEPTION DE LA BOUCLE INTERNE À PARTIR DU LIEU DE RACINES: 
%       1) CHOIX DE Kv
%       2) CALCUL DU NOUVEAU MODELE n1(s)/d1(s) ET DU MODELE A1, B1, C1, D1
%       INCLUANT LA BOUCLE

disp('*******************************************************************')
disp("*** Question e) ***")

% Nouvelles matrices incluant leffet de la boucle interne (voir prob 5 procedural 1)
A1 = A - B(:,2)*Kv*C(1,:)  % on garde A - B2*Kv*C1 car on demande retroaction negative
B1 = B(:,1)
C1 = C(5,:)
D1 = D(5,1)

[num, den] = ss2tf(A1, B1, C1, D1)
TF = tf(num, den)

figure('Name', 'Question e)')
rlocus(TF)

fprintf("\n\n\n")








%% f) VÉRIFICATION DES MARGES AVEC DIAGRAMMES DE BODE ET Kv;
%       COMMENTAIRE SUR LE SENS DES MARGES, UTILITÉ

disp('*******************************************************************')
disp("*** Question f) ***")

FTBF = feedback(v_sur_aprop, Kv)

figure('Name', 'Question f)')
margin(v_sur_aprop)
hold on
margin(FTBF)
grid on
legend('v sur aprop', 'v sur aprop avec Kv')


fprintf("\n\n\n")











%% g) réduction de la FT entre a_prop et v, présentation de la nouvelle FT
%     réduite, fermeture de la boucle à la main sur cette FT, explication
%     de leffet de Kv sur les param standards (K, zeta, wn, tau) et la reponse

disp('*******************************************************************')
disp("*** Question g) ***")

[num, den] = ss2tf(A, B, C, D, 2);
[R, P, K] = residue(num(1,:), den);
poid = abs(R) ./ abs(real(P))
[numR, denR] = residue(R(3:4), P(3:4), K);

% ajustement du gain dc
gain0 = dcgain(num(1,:), den);
gainR = dcgain(numR, denR);
numR = numR * (gain0/gainR);

TF_reduce = tf(numR, denR)
FTBF = feedback(TF_reduce, Kv)
[numR, denR] = tfdata(FTBF, 'v');

% calcul des caracteristiques temporelles
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



% affichage des caracartéristiques temporelles
disp(["--------------------------------------------------------------------------"]);
disp(["Caracartéristiques temporelles FTBF Reduite"]);
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



figure('Name', "Question g)")
step(TF_reduce)
title("Reponse impultionnelle v sur aprop reduite")
grid on

figure('Name', "Question g)")
step(FTBF)
title("Reponse impultionnelle v sur aprop reduite avec gain Kv")
grid on



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

figure('Name', 'Question i)')
rlocus(v_sur_aprop)
hold on
rlocus(TF_reduce)
title('Comparaison entre FT originale et FT reduite')
legend('FT v sur a_prop originale', 'FT reduite')
axis([-1.5    0.1     -1    1])

fprintf("\n\n\n")










%% j) lieu des racines entre gamma_d et gamma incluant la boucle interne
%     (i.e. n1/d1), effet de Kp sur la réponse

disp('*******************************************************************')
disp("*** Question j) ***")

[num,den] = ss2tf(A1,B1,C1,D1, 1);   % le 1 signifie quon veut le 1e element de U (soit delta_c)
gamma_sur_deltaC = tf(num(5), den)         % on veut 

figure('Name', 'Question j)');
rlocus(gamma_sur_deltaC);
title("Rlocus de gamma sur deltaC ZOOM")
axis([-4    0    -6     6])

figure('Name', 'Question j)');
rlocus(gamma_sur_deltaC);
title("Rlocus de gamma sur deltaC")

figure('Name', 'Question j)');
step(gamma_sur_deltaC);
title("Gamma sur deltaC FTBO")
grid on

FTBF_gamma_sur_deltac = feedback(gamma_sur_deltaC, Kv);

figure('Name', 'Question j)');
step(FTBF_gamma_sur_deltac);
title("Gamma sur deltaC FTBF")
grid on

fprintf("\n\n\n")










%% k) design de la boucle externe: avec Bode, 
% (1) calculer Kp limite pour instabilité;
% (2) calculer Kp pour marges de 6 dB et 30 degrés;
% (3) mesurer l'erreur en régime permanent par simulation sur MATLAB

disp('*******************************************************************')
disp("*** Question k) ***")

[num,den] = ss2tf(A1,B1,C1,D1, 1);   % le 1 signifie quon veut le 1e element de U (soit delta_c)
gamma_sur_deltaC = tf(num(5), den) 

num = num(5);
[GM,PM,Wp,Wg] = margin(num, den);

figure('Name', 'Question k)')
margin(num, den)
grid on

Kp_limite = 20*log10(GM)
Kp = 10^((Kp_limite-6)/20)
erp = 1/(1+Kp)








%% l) calcul de l'erreur en régime permanent à partir du gain statique tel
% que lu sur le diagramme de Bode, comparaison avec (3) de l’item ci-dessus


disp('*******************************************************************')
disp("*** Question l) ***")

figure('Name', 'Question l)')
margin(gamma_sur_deltaC)
grid on

Kpos = num(end) / den(end)
erp = 1/(1+Kpos)

fprintf("\n\n\n")











%% m) calcul du nouveau modèle n2(s) / d2(s) et le modèle A2, B2, C2, D2

disp('*******************************************************************')
disp("*** Question m) ***")

A2 = A1 - (B1*Kp*C(5,:))
B2 = B1*Kp
C2 = C(5,:)
D2 = 0

[num2, den2] = ss2tf(A2,B2,C2,D2);
TF2 = tf(num2, den2)

fprintf("\n\n\n")






%% n) comparaison et discussion des réponses temporelles avec
%     compensateurs P, PD, PI et PID

disp('*******************************************************************')
disp("*** Question n) ***")

% compensateur PD
num_PD = Kp * [1    1];
den_PD = [1];
G_PD = tf(num_PD, den_PD)

% compensateur PI
num_PI = Kp * [1    1];
den_PI = [1     0];
G_PI = tf(num_PI, den_PI)

% compensateur PID
num_PID = Kp * [1    1      1];
den_PID = [1        0];
G_PID = tf(num_PID, den_PID)

FTBF_P = feedback(TF2, 1);
FTBF_PD = feedback(G_PD*TF2, 1);
FTBF_PI = feedback(G_PI*TF2, 1);
FTBF_PID = feedback(G_PID*TF2, 1);


figure('Name', 'Question n)')
subplot(4,1,1)
step(FTBF_P)
title('FTBF avec compensateur P')
grid on

subplot(4,1,2)
step(FTBF_PD)
title('FTBF avec compensateur PD')
grid on

subplot(4,1,3)
step(FTBF_PI)
title('FTBF avec compensateur PI')
grid on

subplot(4,1,4)
step(FTBF_PID)
title('FTBF avec compensateur PID')
grid on
