


clc
close all
clear all

% format short g

%% question a) - reduire le systeme a ses poles dominants

disp("*** Question a) ***")
num = [1        3       23];
den = [1        5       22      7       9];
disp("**Fonction de transfert originale")
G = tf(num, den)

figure
pzmap(G)



[R, P, K] = residue(num, den);
poid = abs(R) ./ abs(real(P))
% disp("Residus, poles et poid")
% [R, P, poid]

% on garge que les poles dominants qui m,aximise le test de dominance
[numR, denR] = residue(R(3:4), P(3:4), K);
TF_tempo = tf(numR, denR);


% correction du gain DC
gain0 = dcgain(num, den);
gainR = dcgain(numR, denR);
numR = numR * (gain0/gainR);


disp(' ')
disp("FT reduite")
TF_reduce = tf(numR, denR)

disp(' ')
disp(['Numerateur (reduit): ', num2str(numR)])
disp(['Denominateur (reduit): ', num2str(denR)])
fprintf("\n\n\n")








%% question b) - comparer le lieu des racines

disp("*** Question b) ***")

figure
rlocus(G)
hold on
rlocus(TF_reduce)
title("lieu de racines du sys original et du sys reduit")


figure
rlocus(G)
hold on
% p = rlocus(real(G), imag(G))
% plot(p, 'p')
hold on
rlocus(TF_reduce)
axis([-1    1   -5  5])
title("lieu de racines ZOOM")


fprintf("\n\n\n")








%% question c) - 

disp("*** Question c) ***")
disp("systeme a phase non minimale car il y a un zero du cote droit de laxe imaginaire")
fprintf("\n\n\n")









%% question d) - diagramme de bode

disp("*** Question d) ***")

% fonction originale
figure
bode(G)
title("FT originale")
grid on

figure
margin(G)
grid on

disp("margin value")
[GM_0, PM_0, wp_0, wg_0] = margin(G)


% fonction reduite
figure
bode(TF_reduce)
title("FT reduite")
grid on

figure
margin(TF_reduce)
grid on

[GM_R, PM_R, wp_R, wg_R] = margin(TF_reduce)

fprintf("\n\n\n")








%% question e) - reponse temporelle

disp("*** Question e) ***")
