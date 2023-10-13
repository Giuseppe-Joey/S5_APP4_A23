%%  S5 - APP4 - LABORATOIRE
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



%% probleme 6 - effet de lajout de zweros et de poles
Kp = 1;
wn = 1;
zeta = 0.1;

t = 0:0.1:50;

num = [wn^2];
den = [1    2*zeta*wn   wn^2];
TF = tf(num, den);


% a)
Td = [0.0   0.2     0.5     1.0     2.5     5.0 ];

for i = 1:length(Td)
    num_PD = [Td(i)  1];
    den_PD = [1];
    TF_PD = tf(num_PD, den_PD);

    FTBO = TF * TF_PD ;
    FTBF = feedback(FTBO, 1);

    figure(1)
    rlocus(FTBF) %rlocus ferme la boucle pour nous pour tout K de 0 a linfini
    hold on
    grid on
    
    % lsim ne ferme pas la boucle pour nous alors on ne verrais pas leffet
    % du compensateur 
    %     u = ones(size(t))
    %     yy = lsim()

    
    y = step(FTBF, t);

    figure(2)
    plot(t, y)
    hold on
    legend("Td=0.0", "Td = 0.2", "Td = 0.5", "Td = 1.0", "Td = 2.5", "Td = 5.0")
    grid on
end







% b)
tau = [0.0    0.01    0.05      0.1     0.2];
u = ones(size(tau));

for i = 1:length(tau)
    num = [1];
    den = [tau(i)   1];
    FTBO = tf(num, den);
    FTBF = feedback(FTBO, 1)


%     figure(3)
%     rlocus(FTBF)
%     p = rlocus(num, den, 1);
%     plot(real(p), imag(p), 'p')
%     hold on
% 
% 
%     figure(4)
%     y = lsim(num, den, u, t);
%     plot(y)
%     hold on
end





