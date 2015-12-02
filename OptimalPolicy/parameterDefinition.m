% Purpose of this file:
% define parameters of the model

beta      = 0.99;                % discount factor
eta       = 0.5;                 % relative bargaining power (workers)
alpha     = 0.5;                 % matching parameter
xi        = 0.762;               % matching parameter (constant)
lambda    = 0.05;                % probability of match seperation
k         = 0.598;               % cost of posting vacancy
a         = 0.5;                 % workers unemployment benefit
z         = 1.0;                 % average productivity
phist     = 0.868;               % enforcement parameter (shock)
epsilon   = 6;                   % price elasticity paremeter
omega     = 0.75;                % 0.75; %price adjustment probability
Xst       = epsilon/(epsilon-1); % price ratio/markup
sigma     = 2;                   % risk aversion, comes from household utility function
intrst    = 1/beta;              % interest rate in steady state

% this file defines parameter values
   
paramValues = struct('eta',eta,'beta',beta,'alpha',alpha,'xi',xi,'lambda',lambda,'k',k,'a',a,'z',z,'phist',phist,'epsilon',epsilon,'omega',omega,'Xst',Xst,'sigma',sigma,'intrst',intrst);
    
clear eta beta alpha xi lambda k a z epsilon Xst intrst
save paramValues paramValues
