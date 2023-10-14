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

% calcul des valeurs propres du systeme 
val_propes = eig(A)

% on trouve la fonction de transfert pour le 2e element de la matrice U
% (soit aprop)
[num, den] = ss2tf(A, B, C, D, 2);      
v_sur_aprop = tf(num(1,:), den)

zeros = roots(num(1,:))
poles = roots(den)

figure
pzmap(zeros, poles)
grid on


figure
rlocus(v_sur_aprop)
grid on


% on peut trouver les ft sur chacun des elements de la sortie 
% alpha_sur_aprop = tf(num(2,:), den)
% teta_sur_aprop = tf(num(3,:), den)
% q_sur_aprop = tf(num(4,:), den)
% gamma_sur_aprop = tf(num(5,:), den)

% reduction de la FT avec residue
[R, P, K] = residue(num(1,:), den);
poid = abs(R) ./ abs(real(P))
[numR, denR] = residue(R(3:4), P(3:4), K);
TF_reduce_temporaire = tf(numR, denR)

% ajustement du gain DC
gain0 = dcgain(v_sur_aprop);
gainR = dcgain(TF_reduce_temporaire);

new_num = numR * (gain0/gainR);
TF_reduce = tf(new_num, denR)



zeros = roots(new_num);
poles = roots(denR);

figure
pzmap(zeros, poles)
grid on

figure('Name', 'Question a)')
step(TF_reduce)
title("Réponse temporelle à partir des carac temporelles")
grid on



% selon la formule du     wa = wn* sqrt(1-zeta^2)
wa = abs(imag(poles));
wa = wa(1);

% alpha est egal zeta * wn
alpha = abs(real(poles));
alpha = alpha(1);

wn = abs(poles);
wn = wn(1);

zeta = abs(alpha / wn);
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


% % fonction de transfert a partir des caract/ristiques temporelles
% num = [wn^2];
% den = [1    2*zeta*wn   wn^2];
% TF = tf(num, den)
% 
% figure('Name', 'Question a)')
% step(TF)
% title("Réponse temporelle à partir des carac temporelles")
% grid on







%% b) IDENTIFICATION DE LA FONCTION DE TRANSFERT À PHASE NON-MINIMALE
%       À PARTIR DES PÔLES ET DES ZÉROS












%% c) PRÉSENTATION DU LIEU DES RACINES FAIT À LA MAIN ENTRE a_prop ET v 
%       (AVEC ÉTAPES) ET VALIDATION MATLAB













%% d) EXPLICATION DE L'EFFET DE LA RÉTROACTION Kv SUR LA STABILITÉ,
%       LE TEMPS DE STAB ET LE DEPASS MAX A PARTIR DU LIEU DE RACINES

                                                                                                                                                                                            












%% e) CONCEPTION DE LA BOUCLE INTERNE À PARTIR DU LIEU DE RACINES: 
%       1) CHOIX DE Kv
%       2) CALCUL DU NOUVEAU MODELE n1(s)/d1(s) ET DU MODELE A1, B1, C1, D1
%       INCLUANT LA BOUCLE










%% f) VÉRIFICATION DES MARGES AVEC DIAGRAMMES DE BODE ET Kv;
%       COMMENTAIRE SUR LE SENS DES MARGES, UTILITÉ









%% g) 
