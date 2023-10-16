%%  S5 - APP4 - PROBLÉMATIQUE
%   Auteur:     Giuseppe Lomonaco
%   CIP:        LOMG2301
%   Auteur:     Isabelle Handfield
%   CIP:        HANI1401

%   Date de creation:                       10-Octobre-2023
%   Date de derniere modification:          14-Octobre-2023

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
%     DONNER:       -le temps du premier pic;
%                   -le dépassement maximum;
%                   -le temps de stabilisation;
%                   -le facteur d'amortissement;
%                   -la période des oscillations amorties et naturelles et 
%                    vérifier ces résultats à partir de la réponse
%                    temporelle.

disp(" ")
disp("*** Question a) ***")

% calcul des valeurs propres du systeme 
val_propes = eig(A)

% on trouve la fonction de transfert pour le 2e element de la matrice U
% (soit aprop)
[num, den] = ss2tf(A, B, C, D, 2);      
v_sur_aprop = tf(num(1,:), den)

zeros = roots(num(1,:))
poles = roots(den)

figure('Name', 'Question a)')
pzmap(zeros, poles)



figure('Name', 'Question a)')
rlocus(v_sur_aprop)



% reduction de la FT avec residue
[R, P, K] = residue(num(1,:), den);
poid = abs(R) ./ abs(real(P))
[numR, denR] = residue(R(3:4), P(3:4), K);
TF_temporaire = tf(numR, denR)

% ajustement du gain DC
gain0 = dcgain(v_sur_aprop);
gainR = dcgain(TF_temporaire);

numR = numR * (gain0/gainR);
TF_reduce = tf(numR, denR)



zeros = roots(numR);
poles = roots(denR);

figure('Name', 'Question a)')
pzmap(zeros, poles)


figure('Name', 'Question a)')
step(TF_reduce)
title("Réponse temporelle à partir des carac temporelles")
grid on



% selon la formule du     wa = wn* sqrt(1-zeta^2)
wa = imag(poles)
wa = wa(1)

wn = abs(poles)
wn = wn(1)

zeta = -real(poles) / wn
zeta = zeta(1)

phi = acos(zeta);
phi_degres = acosd(zeta);

ts = 4 / (zeta * wn);
tp = pi / wa;
Mp = 100 * exp(-pi/tan(phi));


% affichage des caracartéristiques temporelles
disp(["---------------------------------------------"]);
disp(["Affichage des caracartéristiques temporelles:"]);
disp(["---------------------------------------------"])
disp(['wn   = ', num2str(wn), ' rad/s']);
disp(['zeta = ', num2str(zeta), ' unites']);
disp(['wa   = ', num2str(wa), ' rad/s']);
disp(['phi  = ', num2str(phi), ' radian']);
disp(['phi  = ', num2str(phi_degres), ' degrés']);
disp(['Mp   = ', num2str(Mp), ' %']);
disp(['ts   = ', num2str(ts), ' s']);
disp(['tp   = ', num2str(tp), ' s']);
disp(["---------------------------------------------"])

fprintf("\n\n\n")











%% b) IDENTIFICATION DE LA FONCTION DE TRANSFERT À PHASE NON-MINIMALE
%       À PARTIR DES PÔLES ET DES ZÉROS

disp("*** Question b) ***")

% on utilise le 2e element de lentree soit a_prop
[num, den] = ss2tf(A, B, C, D, 2)      
v_sur_aprop = tf(num(1,:), den)

zeros = roots(num(1,:))
poles = roots(den)

longueur_zeros = length(zeros);
longueur_poles = length(poles);

[num, den] = zp2tf(zeros, poles, 1);
TF_a_partir_zeros_poles_v_aprop = tf(num,den)


figure('Name', 'Question b)')
rlocus(TF_a_partir_zeros_poles_v_aprop)
title("FT a partir de zp2tf()")

% identification de la fonction de transfert à partir des pôles et zéros
systeme_zpk = zpk(zeros, poles, 1)

figure('Name', 'Question b)')
rlocus(systeme_zpk)
title("FT a partir de zpk()")

fprintf("\n\n\n")







%% c) PRÉSENTATION DU LIEU DES RACINES FAIT À LA MAIN ENTRE a_prop ET v 
%       (AVEC ÉTAPES) ET VALIDATION MATLAB

disp("*** Question c) ***")
disp('**** dessiner a la main le lieu des racines ****')
fprintf("\n\n\n")









%% d) EXPLICATION DE L'EFFET DE LA RÉTROACTION Kv SUR LA STABILITÉ,
%       LE TEMPS DE STAB ET LE DEPASS MAX A PARTIR DU LIEU DE RACINES

disp("*** Question d) ***")

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

disp("*** Question e) ***")

Kv = 1.03;

B1 = B(:,1)
B2 = B(:,2)

C1 = C(1,:)     
C2 = C(2,:)
C3 = C(3,:)
C4 = C(4,:)
C5 = C(5,:)     

%Nouvelles matrices incluant leffet de la boucle interne (voir prob 5procedural 1)
Aa = A - B2*Kv*C1;
% Aa = A + B2*Kv*C1;
Ba = [B1    B2];
Ca = C1;
Da = [0     0];

[num, den] = ss2tf(Aa, Ba, Ca, Da, 2)
TF = tf(num, den)

figure('Name', 'Question e)')
rlocus(TF)

fprintf("\n\n\n")








%% f) VÉRIFICATION DES MARGES AVEC DIAGRAMMES DE BODE ET Kv;
%       COMMENTAIRE SUR LE SENS DES MARGES, UTILITÉ



disp("*** Question f) ***")

figure('Name', 'Question f)')
margin(TF)
grid on


figure('Name', 'Question f)')
step(TF)
grid on

fprintf("\n\n\n")











%% g) RÉDUCTION DE LA FT  ENTRE a_prop ET v presentation de la nouvelle FT reduite
% fermeture de la boucle a la main explications de leffet de Kv sur les param
% standards (K, zeta, wn, tau) et la reponse

disp("*** Question g) ***")

[R, P, K] = residue(num, den);
poid = abs(R) ./ abs(real(P))
[numR, denR] = residue(R(3:4), P(3:4), K)

tf_tempo = tf(numR, denR)

% ajustement du gain dc
gain0 = dcgain(num, den);
gainR = dcgain(numR, denR);
numR = numR * (gain0/gainR);
TF_reduce = tf(numR, denR)


figure('Name', 'Question g)')
step(TF_reduce)
grid on


zeros = roots(numR)
poles = roots(denR)

wn = abs(poles)
wn = wn(1)
% wa = imag(poles)


zeta = -real(poles) / wn
zeta = zeta(1)

phi = acos(zeta)

Mp = 100 * exp(-pi / (tan(phi)))
ts = 4 / (zeta*wn)
tp = pi / (wn*sqrt(1-(zeta^2)))



figure('Name', 'Question g)')






%% h) 










%% i)









%% j)









%% k)










%% l)











%% m)












%% n)














