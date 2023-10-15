%Mode 2
clear all
close all
clc
Annexe_A

%Valeurs propres de la matrice A
valeurs_propes = eig(A);

%Fonction de transfert
[num, den] = ss2tf(A,B,C,D,2)


FT = tf(num(1,:), den)

%zéros et pôles
zeros = roots(num(1,:));
poles = roots(den);

%reduction de la fonction de transfert
[R, P, K] = residue(num(1,:), den);
ratio = abs(R./real(P));
[numN,denN] = residue(R(1:2), P(1:2), K);
NFT = tf(numN,denN);

gainO = dcgain(FT);
gainN = dcgain(NFT);

numR = numN *(gainO/gainN);
RFT = tf(numR,denN)

figure(1)
step(RFT)
grid on

polesN = roots(denN)
%Oscillations naturelles
wn = abs(polesN)
%Oscillations amorties
wa = imag(polesN)
%Facteur d'amortissement
zeta = abs(real(polesN)./wn)
%Phi
phi = acosd(zeta)
%Temps de stabilisation
ts = 4./(zeta.*wn)
%Temps de pic
tp = pi./wa
%Depassement maximal
Mp = 100*exp(-pi./tan(phi))
