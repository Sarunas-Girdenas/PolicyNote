
%Steady state solution for Monacelli, Quadrini, Trigarri (2011
function [ys,check]=optimalPolicy_steadystate(ys,exe)
global M_ lgy_
global beta eta alpha xi lambda k a z phist omega sigma epsilon Sst Xst bst intrst qst ust Nst Cst Kst Zst G1s G2s G3s G4s G5s G6s G7s G8s G9s wst
if isfield(M_,'param_nbr') == 1
NumberOfParameters = M_.param_nbr;
for i = 1:NumberOfParameters
  paramname = deblank(M_.param_names(i,:));
  eval([ paramname ' = M_.params(' int2str(i) ');']);
end
check = 0;
end

load paramValues

beta      = paramValues.beta;    % discount factor
eta       = paramValues.eta;     % relative bargaining power (workers)
alpha     = paramValues.alpha;   % matching parameter
xi        = paramValues.xi;      % matching parameter (constant)
lambda    = paramValues.lambda;  % probability of match seperation
k         = paramValues.k;       % cost of posting vacancy
a         = paramValues.a;       % Workers unemployment benefit
z         = paramValues.z;       % average productivity
phist     = paramValues.phist;   % enforcement parameter (shock)
epsilon   = paramValues.epsilon; % price elasticity paremeter
omega     = paramValues.omega;   % price adjustment probability
Xst       = paramValues.Xst;     % price ratio/markup
sigma     = paramValues.sigma;   % risk aversion, comes from household utility function
intrst    = paramValues.intrst;  % steady state interest

%Here we solve for the steady state

[ std_Variables, std_Multipliers ] = steadyStateFunction(paramValues);

Sst = std_Variables.Sst;
bst = std_Variables.bst;
qst = std_Variables.qst;
Cst = std_Variables.Cst;
Nst = std_Variables.Nst;
Kst = std_Variables.Kst;
Zst = std_Variables.Zst;
ust = std_Variables.ust;
wst = std_Variables.wst;

G1s = std_Multipliers.G1s;
G2s = std_Multipliers.G2s;
G3s = std_Multipliers.G3s;
G4s = std_Multipliers.G4s;
G5s = std_Multipliers.G5s;
G6s = std_Multipliers.G6s;
G7s = std_Multipliers.G7s;
G8s = std_Multipliers.G8s;
G9s = std_Multipliers.G9s;

Kplus    = 0;
bplus    = 0;
Zplus    = 0;
inflplus = 0;
Cplus    = 0;
Splus    = 0;
Nminus   = 0;

G1       = 0;
G2       = 0;
G3       = 0;
G4       = 0;
G5       = 0;
G6       = 0;
G7       = 0;
G8       = 0;
G9       = 0;
S        = 0;
X        = 0;
b        = 0;
intr     = 0;
infl     = 0;
q        = 0;
N        = 0;
C        = 0;
phi      = 0;
K        = 0;
Z        = 0;
w        = 0;
u        = 0;

for iter = 1:length(M_.params)
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

if isfield(M_,'param_nbr') == 1

if isfield(M_,'orig_endo_nbr') == 1
NumberOfEndogenousVariables = M_.orig_endo_nbr;
else
NumberOfEndogenousVariables = M_.endo_nbr;
end
ys = zeros(NumberOfEndogenousVariables,1);
for i = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(i,:));
  eval(['ys(' int2str(i) ') = ' varname ';']);
end
else
ys=zeros(length(lgy_),1);
for i = 1:length(lgy_)
    ys(i) = eval(lgy_(i,:));
end
check = 0;
end
end