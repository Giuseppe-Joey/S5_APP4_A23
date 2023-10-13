



clc
close all
clear all





%% Probleme 15

% a)
masse = 1;                  % en kg  
Force_moteur = 20.4;        % en Newton
b = 10;                     % coefficient de friction dynamique en N/(m/s)



% on doit faire une entree PWM de + ou - 1....
% Parameters
Fs = 1000;          % Sampling frequency (Hz)
T = 1;              % Total time (seconds)
DutyCycle = 0.4;    % Duty cycle (between 0 and 1)
Frequency = 100;    % PWM frequency (Hz)

% Create a time vector
t = 0:1/Fs:T;

% Generate PWM signal
pwm_signal = PWM(t, Frequency, DutyCycle);

% Plot the PWM signal
plot(t, pwm_signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('PWM Signal');
grid on;




% a)








% b)








% c)










