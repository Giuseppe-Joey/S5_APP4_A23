

clc
close all
clear all


%% Probleme 14

% a)
num = [17.8885];
den = [1    2       0];

% soit la FTBO suivante
G = tf(num, den)

figure
bode(G)
grid on






% b)
numGr = 0.9997 * [1     0.10];
denGr =  [1     0.02];
Gr = tf(numGr, denGr)


figure
margin(G)
hold on
margin(G*Gr)
grid on
legend("Original", "Compense RePh")





% c)







% d)









% % e)
% syms s
% Kp = 0.9997;
% Ki = 0.10;
% numPI = Kp + Ki/s
% 
% 





