%%  S5 - APP4 - GUIDE MATLAB - CHAPITRE_18_reponse_en_frequence.M
%   Auteur:     Giuseppe Lomonaco
%   CIP:        LOMG2301
%   Auteur:     Isabelle Handfield
%   CIP:        HANI1401

%   Date de creation:                       10-Octobre-2023
%   Date de derniere modification:          10-Octobre-2023

%   DESCRIPTION: REPONSE EN FREQUENCE

clc
close all
clear all





%% Bode sans arguments à gauche
num = [3 12 9]*200;
den = [1 23 206 902 1948 1680];
FT = tf(num,den);

[A, B, C, D] = tf2ss(num,den);
VE = ss(A, B, C, D);

figure
bode(num,den)

figure
bode(FT)

figure
bode(A, B, C, D)

figure
bode(VE)









%% Bode avec arguments à gauche
A = [0 1 0; 0 0 1; -2 -5 -1];
B = [0 2; 1 0; 1 -1]; % deux entrées (r=2)
C = [3 -1 -2; 0 1 1; 2 1 0]; % trois sorties (p=3)
D = zeros(3,2);

VE = ss(A, B, C, D);
[mag,pha,wVE] = bode(VE);

[p,r,nw] = size(mag);

magVE = reshape(mag,[p*r,nw])';
phaVE = reshape(pha,[p*r,nw])';

[magABCD, phaABCD, wABCD] = bode(A,B,C,D);

subplot(2,2,1)
loglog(wVE, magVE)
title('Bode avec forme objet (VE)')

subplot(2,2,3)
semilogx(wVE, phaVE)
xlabel('Fréquence (rad/s)')

subplot(2,2,2)
loglog(wABCD, magABCD)
title('Bode avec forme matricielle (ABCD)')

subplot(2,2,4)
semilogx(wABCD, phaABCD)
xlabel('Fréquence (rad/s)')








%% Bode avec logspace
num = [3 12 9]*200;
den = [1 23 206 902 1948 1680];
FT = tf(num,den);

[A, B, C, D] = tf2ss(num,den);
VE = ss(A, B, C, D);

w = logspace(-1,2,200);

bode(num,den,w)
bode(FT,w)
bode(VE,w)






%% Margin sans arguments à gauche
num = [3 12 9]*200;
den = [1 23 206 902 1948 1680];
FT = tf(num,den);
[A, B, C, D] = tf2ss(num,den);
VE = ss(A, B, C, D);

figure
margin(num,den)
grid on

figure
margin(FT)

figure
margin(A, B, C, D)

figure
margin(VE)








%% Margin avec arguments à gauche
num = [3 12 9]*200;
den = [1 23 206 902 1948 1680];
FT = tf(num,den);

[A, B, C, D] = tf2ss(num,den);
VE = ss(A, B, C, D);

margin(num,den)
[GM,PM,wp,wg] = margin(num,den);

disp(' ')
disp(['La marge de gain est ', num2str(GM),' à la fréquence ', num2str(wp), ' rad/s'])
disp(['La marge de phase est ', num2str(PM),' deg à la fréquence ', num2str(wg), ' rad/s'])

[GM,PM,wp,wg] = margin(FT);
disp(['La marge de gain est ', num2str(GM),' à la fréquence ', num2str(wp), ' rad/s'])
disp(['La marge de phase est ', num2str(PM),' deg à la fréquence ', num2str(wg), ' rad/s'])

[GM,PM,wp,wg] = margin(A,B,C,D);
disp(['La marge de gain est ', num2str(GM),' à la fréquence ', num2str(wp), ' rad/s'])
disp(['La marge de phase est ', num2str(PM),' deg à la fréquence ', num2str(wg), ' rad/s'])

[GM,PM,wp,wg] = margin(VE);
disp(['La marge de gain est ', num2str(GM),' à la fréquence ', num2str(wp), ' rad/s'])
disp(['La marge de phase est ', num2str(PM),' deg à la fréquence ', num2str(wg), ' rad/s'])










%% Effet du gain

margin(num,den)
[GM,PM,wp,wg] = margin(num,den);
disp(' ')
disp(['La marge de gain est ', num2str(GM),' à la fréquence ', num2str(wp), ' rad/s'])
disp(['La marge de phase est ', num2str(PM),' deg à la fréquence ', num2str(wg), ' rad/s'])

margin(10*num,den)
[GM,PM,wp,wg] = margin(10*num,den);
disp(' ')
disp(['La marge de gain est ', num2str(GM),' à la fréquence ', num2str(wp), ' rad/s'])
disp(['La marge de phase est ', num2str(PM),' deg à la fréquence ', num2str(wg), ' rad/s'])








