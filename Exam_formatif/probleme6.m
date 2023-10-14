


clc
close all
clear all


%% question a)

disp("*** Question a) ***")
num = [1        3       23];
den = [1        5       22      7       9];
G = tf(num, den)

poles = roots(den)

disp("--> Deux paires de poles conjugu/s complexes donc 2 modes dynamiques!! <--")
disp("--------------------------------------------------------------------------- ")
disp(" ")




%% question b) - trouver Mp, tp et ts
disp("*** Question b) ***")

% on trouve zeta et wn grace aux poles
poles = roots(den)
wn = abs(poles)
zeta = -real(poles) ./ wn
wa = imag(poles)
zeta_wn = abs(real(poles))

% on trouve les autres valeurs a partir de zeta et wn
phi = acosd(zeta)
Mp = 100 * exp(-pi ./ tand(phi))
ts = 4 ./ (zeta .* wn)
tp = pi ./ wa








%% question c) - Lieu des racines
disp("*** Question c) ***")

% lieu de racines
figure
rlocus(G, 5000)
plot(real(poles), imag(poles), 'p')
hold on
rlocus(G)

ts1 = 4 ./ abs(real(poles))




% regle numro 5 - intersextion avec axe reel
zeros = roots(num)
n_moins_m = length(poles) - length(zeros)
sigma = (sum(poles) - sum(zeros)) / n_moins_m      % intersection avec laxe real
ts2 = 4 / abs(sigma)






