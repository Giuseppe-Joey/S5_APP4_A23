%%  S5 - APP4 - LABORATOIRE
%   Auteur:     Giuseppe Lomonaco
%   CIP:        LOMG2301
%   Auteur:     Isabelle Handfield
%   CIP:        HANI1401

%   Date de creation:                       10-Octobre-2023
%   Date de derniere modification:          10-Octobre-2023

%   DESCRIPTION: 




clc
close all
clear all


%% Probleme 8 - Compensateur P, PD et avance de phase

% fonction de transfert
num = [4];
den = [1    2   0];
TF = tf(num, den);

% compensateur proportionnel de gain unitaire
Gc = 1;

% compensateur avance de phase 
num_AvPh = 4.68 * [1    2.9];
den_AvPh = 4.68 * [1    5.4];
Gc_AvPh = tf(num_AvPh, den_AvPh);

% compensateur PD
num_PD = 4.68 * [1    2.9];
den_PD = 4.68 * [5.4];
Gc_PD = tf(num_PD, den_PD);



% question a) 
FTBO = Gc * TF;
FTBF = feedback(FTBO, 1)
rlocus(FTBF)

[numFTBF, denFTBF] = tfdata(FTBF, 'v');
wn = sqrt(denFTBF(3))
zeta = denFTBF(2) / (wn*2)
ts = 4 / (zeta * wn)




% question b)
FTBO = Gc_AvPh * TF
% FTBF = feedback(FTBO, 1)
[numFTBO, denFTBO] = tfdata(FTBO, 'v');

[R, P, K] = residue(numFTBO, denFTBO);
poid = abs(R) ./ abs(real(P))
[numR, denR] = residue(R(2:3), P(2:3), K);
TF_Reduce = tf(numR, denR)
% FTBF = feedback(TF_Reduce, 1);
[numFTBF, denFTBF] = tfdata(TF_Reduce, 'v');

wn = sqrt(denFTBF(3))
zeta = denFTBF(2) / (wn*2)
ts = 4 / (zeta * wn)
