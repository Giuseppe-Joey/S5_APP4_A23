%%  S5 - APP4 - PROBLÉMATIQUE
%   Auteur:     Giuseppe Lomonaco
%   CIP:        LOMG2301
%   Auteur:     Isabelle Handfield
%   CIP:        HANI1401

%   Date de creation:                       10-Octobre-2023
%   Date de derniere modification:          10-Octobre-2023

%   DESCRIPTION: PROBLÉMATIQUE


clc
close all
clear all

% Chargement du fichier Annexe A
fprintf("*** Chargement du fichier 'Annexe_A' ***\n");
Annexe_A

%% a) ANALYSE DES CARACTÉRISTIQUES DYNAMIQUES DE L'AVION
%     À PÂRTIR DES VALEURS PROPRES DU SYSTEME ET DE LA RÉPONSE TEMPORELLE
%     DE LA VITESSE v DU SYSTEME SOUMIS A UN ÉCHELON SUR a_prop
%     DONNER:       -le temps du premier pic;
%                   -le dépassement maximum;
%                   -le temps de stabilisation;
%                   -le facteur d'amortissement;
%                   -la période des oscillations amorties et naturelles et 
%                    vérifier ces résultats à partir de la réponse
%                    temporelle.


val_propes = eig(A)

% - TEMPS DU PREMIER PIC 


% creation dune matrice identite de la meme taille que val propres
I = eye(size(val_propes))
den = det((val_propes.*I) - A)






%% b) IDENTIFICATION DE LA FONCTION DE TRANSFERT À PHASE NON-MINIMALE
%       À PARTIR DES PÔLES ET DES ZÉROS

% iu = 1;
% TS = 100;
% [num, den] = ss2tf(A, B, C, D, iu)
% tf = tf(num, den)









%% c) PRÉSENTATION DU LIEU DES RACINES FAIT À LA MAIN ENTRE a_prop ET v 
%       (AVEC ÉTAPES) ET VALIDATION MATLAB













%% d) EXPLICATION DE L'EFFET DE LA RÉTROACTION Kv SUR LA STABILITÉ,
%       LE TEMPS DE STAB ET LE DEPASS MAX A PARTIR DU LIEU DE RACINES

                                                                                                                                                                                            












%% e) CONCEPTION DE LA BOUCLE INTERNE À PARTIR DU LIEU DE RACINES: 
%       1) CHOIX DE Kv
%       2) CALCUL DU NOUVEAU MODELE n1(s)/d1(s) ET DU MODELE A1, B1, C1, D1
%       INCLUANT LA BOUCLE










%% f) VÉRIFICATION DES MARGES AVEC DIAGRAMMES DE BODE ET Kv;
%       COMMENTAIRE SUR LE SENS DES MARGES, UTILITÉ









%% g) 
