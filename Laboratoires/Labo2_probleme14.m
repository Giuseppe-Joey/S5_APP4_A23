

clc
close all
clear all


%% Probleme 14

% a)
num = [17.8885];
den = [1    2       0];

% soit la FTBO suivante
G = tf(num, den)









% b)
numGr = 0.9997 * [1     0.10];
denGr = 0.9997 * [1     0.02];
Ga = tf(numGa, denGa)








% c)







% d)









% e)
syms s
Kp = 0.9997;
Ki = 0.10;
numPI = Kp + Ki/s







