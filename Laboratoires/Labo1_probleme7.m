%%  S5 - APP4 - LABORATOIRE
%   Auteur:     Giuseppe Lomonaco
%   CIP:        LOMG2301
%   Auteur:     Isabelle Handfield
%   CIP:        HANI1401

%   Date de creation:                       10-Octobre-2023
%   Date de derniere modification:          10-Octobre-2023

%   DESCRIPTION: REPONSE EN FREQUENCE



%% probleme 7 - Effet du retard
clc 
close all
clear all

Kp = 1;
wn = 1;
zeta = 0.1;

t = [0:0.1:50]';

num = [wn^2];
den = [1    2*zeta*wn   wn^2];
TF = tf(num, den);

Tr = [0.05  0.10    0.20    0.25];

u = ones(size(t));
ordre_pade = 6;



for i=1:length(Tr)
    num = [-Tr(i)   0];
    den = [1];

    tf_sys = tf(num, den);
    exp_tf = exp(tf_sys);
    [num1, den1] = pade(Tr(i), ordre_pade);
    
    FTBO = TF * tf(num1, den1);


    figure(1)
    rlocus(FTBO);
    FTBF = feedback(FTBO, 1);
    hold on


    figure(2)
    y = lsim(FTBF, u, t);
    plot(t, y)
    hold on

end








