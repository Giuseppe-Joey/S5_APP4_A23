%%  S5 - APP4 - GUIDE MATLAB - CHAPITRE_17.M
%   Auteur:     Giuseppe Lomonaco
%   CIP:        LOMG2301
%   Auteur:     Isabelle Handfield
%   CIP:        HANI1401

%   Date de creation:                       10-Octobre-2023
%   Date de derniere modification:          10-Octobre-2023

%   DESCRIPTION: LIEU DES RACINES

clc
close all
clear all


%% Lieu des racines
num = [3 12 9]*100;
den = [1 23 206 902 1948 1680];
FT = tf(num,den);

[A, B, C, D] = tf2ss(num,den);
VE = ss(A, B, C, D);

figure
rlocus(num,den)
axis([-10, 2, -12, 12])

figure
rlocus(FT)









%% Ajout de pôles à des gains spécifiques
num = [3 12 9]*100;
den = [1 23 206 902 1948 1680];

figure
rlocus(num,den)
axis([-15, 2, -12, 12])
hold on
p = rlocus(num,den,1); % Calcul des pôles de la FTBF à K = 1
plot(real(p), imag(p), 'p') % Utilise le pentagramme
p = rlocus(num,den,3); % Calcul des pôles de la FTBF à K = 3
plot(real(p), imag(p), 's') % Utilise le carré

figure
rlocus(A, B, C, D)

figure
rlocus(VE)





%% Effet du gain
num = [3 12 9]*100;
den = [1 23 206 902 1948 1680];
FT1 = tf(num,den);
FT3 = 3*FT1;
figure
rlocus(FT1)
axis([-15, 2, -12, 12])
hold on
p = rlocus(FT1,3); % Tracé des pôles de FT1 à K = 3
plot(real(p), imag(p), 'p') % Utilise le pentagramme
rlocus(FT3)
p = rlocus(FT3,1); % Tracé des pôles de FT3 à K = 1
plot(real(p), imag(p), 's') % Utilise le carré








