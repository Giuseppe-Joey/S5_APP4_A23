%%  S5 - APP4 - PROBLEMATIQUE - ANNEXE_A.M
%   Auteur:     Giuseppe Lomonaco
%   CIP:        LOMG2301
%   Auteur:     Isabelle Handfield
%   CIP:        HANI1401

%   Date de creation:                       10-Octobre-2023
%   Date de derniere modification:          17-Octobre-2023

%   DESCRIPTION: fichier fourni en Annexe pour la problématique


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
      
      
%  delta_c;       % angle du gouvernail de profondeur
%  a_prop;        % fraction de la poussée maximale des moteurs (0 < aprop < 1)
%  v;             % vitesse de l<avion
%  alpha;         % angle d<attaque (entre l<axe longitudinal de lavion et la vitesse)
%  teta;          % angle de tanguage( entre l<axe longitudinal de lavion et le plan horizontal)
%  q;             % vitesse angulaire en tanguage de lavion (q = d teta / dt)
%  gamma;         % angle de vol de lavion ( gamma = teta - alpha): si gamma > 0, lavion monte, si gamma < 0 lavion descend


%  % Variables d'entrée (actionneurs)
%  u = [delta_c   a_prop]';        % en degres et en fraction de la poussee maximale

%  % variables d'états
%  x = [v     alpha   teta    q]';
 
%  % variables de sorties (capteurs)
%  y = [v     alpha   teta    q   gamma]';
 
 


