

clc
close all
clear all


%% probleme 13

% a)
num = [17.8885];
den = [1    2       0];

% soit la FTBO suivante
G = tf(num, den)





% b)
numGa = 1.5697 * [1     2.5483];
denGa = 1.5697 * [1     6.2787];
Ga = tf(numGa, denGa)







% c) 










% d)










% e)
Kp = 1;
Td = 0:0.1:100;
numGPD = [Kp         Kp * Td];
denGPD = [1];
GPD = tf(numGPD, denGPD)



