% this file computes steady state figures

% firstly, load parameter values

load paramValues

% construct the interval for phi in steady state

intervalLength = 50;

phistInt = linspace(0.25,0.75,intervalLength);

storeUnemployment = zeros(1,intervalLength);
storeWage         = zeros(1,intervalLength);


% Here we solve for the steady state variables

for i = 1:intervalLength
    
    % different values of phi
    
    paramValues.phist = phistInt(i);
    
    % solve steady state
    
    syms Sst bst qst Cst Nst

    eq1 = -Sst + paramValues.z/paramValues.Xst - paramValues.a - bst + paramValues.intrst*bst + (1-paramValues.lambda)*paramValues.beta*Sst - paramValues.eta*paramValues.beta*(qst^(-paramValues.alpha/(1-paramValues.alpha)))*(paramValues.xi^(paramValues.alpha/(1-paramValues.alpha)))*Sst;
    eq2 = -bst*paramValues.intrst + paramValues.phist*(1-paramValues.eta)*Sst;
    eq3 = -paramValues.k/qst + paramValues.intrst*bst + paramValues.beta*(1-paramValues.eta)*Sst;
    eq4 = -paramValues.z*Nst + Cst - paramValues.a*(1-Nst) + paramValues.k*(qst^(-1/(1-paramValues.alpha)))*(paramValues.xi^((2*paramValues.alpha-1)/(paramValues.alpha*(1-paramValues.alpha))))*(1-(1-paramValues.lambda)*Nst);
    eq5 = -Nst + (1-paramValues.lambda)*Nst + (qst^(-paramValues.alpha/(1-paramValues.alpha)))*(paramValues.xi^(paramValues.alpha/(1-paramValues.alpha)))*(1-(1-paramValues.lambda)*Nst);

    sol = solve(eq1,eq2,eq3,eq4,eq5);

    Sst  = double(sol.Sst(2));
    bst  = double(sol.bst(2));
    qst  = double(sol.qst(2));
    Cst  = double(sol.Cst(2));
    Nst  = double(sol.Nst(2));
    
    Kst = Nst*paramValues.z/(1-paramValues.omega*paramValues.beta);
    Zst = Nst/(paramValues.Xst*(1-paramValues.omega*paramValues.beta));
    ust = 1-(1-paramValues.lambda)*Nst;

    wst = (paramValues.z/paramValues.Xst)*paramValues.eta + paramValues.a*(1-paramValues.eta)+(paramValues.eta*paramValues.k*(paramValues.phist+paramValues.beta*(qst^((-paramValues.alpha)/(1-paramValues.alpha))))*(paramValues.xi^(paramValues.alpha/(1-paramValues.alpha))))/(qst*(paramValues.phist+paramValues.beta));

    % store results
    
    storeUnemployment(1,i) = ust;
    storeWage(1,i)         = wst;
    
end

% this is figure 2 from the paper

figure(1)
[AX,H1,H2] = plotyy(phistInt,storeUnemployment,phistInt,storeWage);
set(get(AX(1),'Ylabel'),'String','Unemployment') 
set(get(AX(2),'Ylabel'),'String','Wage')
set(AX,{'ycolor'},{'k';'k'})
xlabel('Value of Borrowing Parameter, \phi')
hold(AX(1))
plot([0.433 0.433], [0.06 0.075])


% compute figure 1 from the paper

clearvars -except paramValues

intervalLength = 50;

etaInt = linspace(0.25,0.75,intervalLength);

storeUnemployment = zeros(1,intervalLength);
storeborrowing    = zeros(1,intervalLength);


% Here we solve for the steady state variables

for i = 1:intervalLength
    
    % different values of phi
    
    paramValues.eta = etaInt(i);
    
    % solve steady state
    
    syms Sst bst qst Cst Nst

    eq1 = -Sst + paramValues.z/paramValues.Xst - paramValues.a - bst + paramValues.intrst*bst + (1-paramValues.lambda)*paramValues.beta*Sst - paramValues.eta*paramValues.beta*(qst^(-paramValues.alpha/(1-paramValues.alpha)))*(paramValues.xi^(paramValues.alpha/(1-paramValues.alpha)))*Sst;
    eq2 = -bst*paramValues.intrst + paramValues.phist*(1-paramValues.eta)*Sst;
    eq3 = -paramValues.k/qst + paramValues.intrst*bst + paramValues.beta*(1-paramValues.eta)*Sst;
    eq4 = -paramValues.z*Nst + Cst - paramValues.a*(1-Nst) + paramValues.k*(qst^(-1/(1-paramValues.alpha)))*(paramValues.xi^((2*paramValues.alpha-1)/(paramValues.alpha*(1-paramValues.alpha))))*(1-(1-paramValues.lambda)*Nst);
    eq5 = -Nst + (1-paramValues.lambda)*Nst + (qst^(-paramValues.alpha/(1-paramValues.alpha)))*(paramValues.xi^(paramValues.alpha/(1-paramValues.alpha)))*(1-(1-paramValues.lambda)*Nst);

    sol = solve(eq1,eq2,eq3,eq4,eq5);

    Sst  = double(sol.Sst(2));
    bst  = double(sol.bst(2));
    qst  = double(sol.qst(2));
    Cst  = double(sol.Cst(2));
    Nst  = double(sol.Nst(2));
    
    Kst = Nst*paramValues.z/(1-paramValues.omega*paramValues.beta);
    Zst = Nst/(paramValues.Xst*(1-paramValues.omega*paramValues.beta));
    ust = 1-(1-paramValues.lambda)*Nst;

    wst = (paramValues.z/paramValues.Xst)*paramValues.eta + paramValues.a*(1-paramValues.eta)+(paramValues.eta*paramValues.k*(paramValues.phist+paramValues.beta*(qst^((-paramValues.alpha)/(1-paramValues.alpha))))*(paramValues.xi^(paramValues.alpha/(1-paramValues.alpha))))/(qst*(paramValues.phist+paramValues.beta));

    % store results
    
    storeUnemployment(1,i) = ust;
    storeborrowing(1,i)    = bst;
    
end

figure(2)
[AX,H1,H2] = plotyy(etaInt,storeUnemployment,etaInt,storeborrowing);
set(get(AX(1),'Ylabel'),'String','Unemployment') 
set(get(AX(2),'Ylabel'),'String','Borrowing')
set(AX,{'ycolor'},{'k';'k'})
xlabel('Value of Bargaining Power, \eta')

% compute second figure in the two figures block

% first, choose internvals for eta, and two params in Taylor rule

% Taylor rule coefs
inflBegin = 1.5; 
inflEnd   = 3; 
uBegin    = 0.01;
uEnd      = 2;

% eta

infl_Int = linspace(1.5,3,11); 
un_Int   = linspace(0.01,-2,11);
eta_Int  = linspace(0.25,0.75,11);

pairs = zeros([ length(infl_Int)^3 3 ]);

i = 0;
    
for val1 = infl_Int
        
    for val2 = un_Int
        
        for val3 = eta_Int
            
            i = i + 1;
            
            pairs(i,1) = val1;
            pairs(i,2) = val2;
            pairs(i,3) = val3;
            
        end
            
    end
        
end

u_IntFig    = pairs(:,2);
infl_IntFig = pairs(:,1);
eta_IntFig  = pairs(:,3);

save u_IntFig u_IntFig
save eta_IntFig eta_IntFig
save infl_IntFig infl_IntFig





    
