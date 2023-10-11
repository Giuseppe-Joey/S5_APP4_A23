%%  S5 - APP4 - PROBLEMATIQUE - ANNEXE_A.M
%   Auteur:     Giuseppe Lomonaco
%   CIP:        LOMG2301
%   Auteur:     Isabelle Handfield
%   CIP:        HANI1401

%   Date de creation:                       10-Octobre-2023
%   Date de derniere modification:          10-Octobre-2023

%   DESCRIPTION: fichier fourni en Annexe pour la problématique

clc
close all
clear all

% Modèle linéaire de la dynamique longitudinale d'un avion
A = [-0.018223  -0.088571   -9.78   0;
     -0.003038  -1.2563     0       1;
     0          0           0       1;
     0.0617     -28.075     0       -4.5937];
 
 B = [0     1.1962;
      0     -0.00120;
      0     0;
      7.84 -4.05];
  
 C = [1     0       0       0;
      0     57.296  0       0;
      0     0       57.296  0;
      0     0       0       57.296;
      0 -57.296     57.296  0];
      
 D = [0     0;
      0     0;
      0     0;
      0     0;
      0     0];
      
      
 delta_c = 1;
 a_prop = 1;
 v = 1;
 alpha = 1;
 teta = 1;
 q = 1;
 gamma = 1;
 
 
 u = [delta_c a_prop]';        % en degres et en fraction de la poussee maximale
 x = [v alpha teta q]';
 y = [v alpha teta q gamma]';
 
 


