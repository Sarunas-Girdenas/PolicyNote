 % Purpose: solve for steady state
 % INPUT:  struct of parameter values
 % OUTPUT: struct of steady state variables
 %         struct of steady state Lagrange Multipliers
 
 function [ std_Variables, std_Multipliers ] = steadyStateFunction(paramStruct)
 
    % Here we solve for the steady state variables

    syms Sst bst qst Cst Nst

    eq1 = -Sst + paramStruct.z/paramStruct.Xst - paramStruct.a - bst + paramStruct.beta*bst + (1-paramStruct.lambda)*paramStruct.beta*Sst - paramStruct.eta*paramStruct.beta*(qst^(-paramStruct.alpha/(1-paramStruct.alpha)))*(paramStruct.xi^(paramStruct.alpha/(1-paramStruct.alpha)))*Sst;
    eq2 = -bst*paramStruct.beta + paramStruct.phist*(1-paramStruct.eta)*Sst;
    eq3 = -paramStruct.k/qst + paramStruct.beta*bst + paramStruct.beta*(1-paramStruct.eta)*Sst;
    eq4 = -paramStruct.z*Nst + Cst - paramStruct.a*(1-Nst) + paramStruct.k*(qst^(-1/(1-paramStruct.alpha)))*(paramStruct.xi^((2*paramStruct.alpha-1)/(paramStruct.alpha*(1-paramStruct.alpha))))*(1-(1-paramStruct.lambda)*Nst);
    eq5 = -Nst + (1-paramStruct.lambda)*Nst + (qst^(-paramStruct.alpha/(1-paramStruct.alpha)))*(paramStruct.xi^(paramStruct.alpha/(1-paramStruct.alpha)))*(1-(1-paramStruct.lambda)*Nst);

    sol = solve(eq1,eq2,eq3,eq4,eq5);

    Sst  = double(sol.Sst(2));
    bst  = double(sol.bst(2));
    qst  = double(sol.qst(2));
    Cst  = double(sol.Cst(2));
    Nst  = double(sol.Nst(2));
    
    Kst = Nst*paramStruct.z/(1-paramStruct.omega*paramStruct.beta);
    Zst = Nst/(paramStruct.Xst*(1-paramStruct.omega*paramStruct.beta));
    ust = 1-(1-paramStruct.lambda)*Nst;
    wst = (paramStruct.z/paramStruct.Xst-bst)*paramStruct.eta + paramStruct.a*(1-paramStruct.eta)+paramStruct.eta*paramStruct.k*(paramStruct.phist+paramStruct.beta*qst^(-paramStruct.alpha/(1-paramStruct.alpha))*(paramStruct.xi^(paramStruct.alpha/(1-paramStruct.alpha))))/(qst*(paramStruct.phist+paramStruct.beta));
    
    % make a struct of steady state variables
    
    std_Variables = struct('Sst',Sst,'bst',bst,'qst',qst,'Cst',Cst,'Nst',Nst,'Kst',Kst,'Zst',Zst,'ust',ust,'wst',wst);
        
    % Solve for steady state Lagrange multipliers
    
    syms G7s G8s G9s G1s G2s G3s G4s G5s
    
    e1  = G7s*((paramStruct.epsilon-1)/paramStruct.epsilon)*paramStruct.z-G7s*(paramStruct.z/paramStruct.Xst)+G8s*paramStruct.z-paramStruct.a*G8s+paramStruct.beta*G8s*(1-paramStruct.lambda)*paramStruct.k*(qst^(-1/(1-paramStruct.alpha)))*(paramStruct.xi^((2*paramStruct.alpha-1)/(paramStruct.alpha*(1-paramStruct.alpha))))+G9s-paramStruct.beta*G9s*(1-paramStruct.lambda)+paramStruct.beta*G9s*(qst^(-paramStruct.alpha/(1-paramStruct.alpha)))*((paramStruct.xi^(paramStruct.alpha/(1-paramStruct.alpha))))*(1-paramStruct.lambda);
    e2  = (-paramStruct.alpha/(1-paramStruct.alpha))*G1s*paramStruct.eta*paramStruct.beta*(qst^(-1/(1-paramStruct.alpha)))*(paramStruct.xi^(paramStruct.alpha/(1-paramStruct.alpha)))*Sst-G3s*(paramStruct.k/(qst^2))+(1/(1-paramStruct.alpha))*G8s*paramStruct.k*(qst^((paramStruct.alpha-2)/(1-paramStruct.alpha)))*(paramStruct.xi^(1/(1-paramStruct.alpha)))*(1-(1-paramStruct.lambda)*Nst)+(paramStruct.alpha/(1-paramStruct.alpha))*G9s*(qst^(-1/(1-paramStruct.alpha)))*(paramStruct.xi^(paramStruct.alpha/(1-paramStruct.alpha)))*(1-(1-paramStruct.lambda)*Nst);
    e3  = -G1s*bst+G2s*bst-G3s*bst+G4s+G5s*(paramStruct.omega/(1-paramStruct.omega))*Kst+G7s*((paramStruct.epsilon-1)/paramStruct.epsilon)*paramStruct.omega*(paramStruct.epsilon-1)*Kst-G7s*paramStruct.omega*paramStruct.epsilon*Zst;
    e4  = (Cst^(1-paramStruct.sigma))-G4s*paramStruct.sigma*paramStruct.beta+G4s*paramStruct.sigma - G8s*Cst;
    e5  = G1s-G1s*(1-paramStruct.lambda)+G1s*paramStruct.eta*(qst^(-paramStruct.alpha/(1-paramStruct.alpha)))*(paramStruct.xi^(paramStruct.alpha/(1-paramStruct.alpha)))-(G2s*paramStruct.phist*(1-paramStruct.eta)/paramStruct.beta)-G3s*(1-paramStruct.eta);
    e6  = G1s*paramStruct.z + G7s*paramStruct.z*Nst;
    e7  = G2s - G3s;
    e8  = G1s*bst-G4s-G7s*((paramStruct.epsilon-1)/paramStruct.epsilon)*paramStruct.omega*Kst+G7s*paramStruct.omega*Zst;
  
    sln = solve(e1,e2,e3,e4,e5,e6,e7,e8);
    
    G7s = double(sln.G7s);
    G8s = double(sln.G8s);
    G9s = double(sln.G9s);
    G1s = double(sln.G1s);
    G2s = double(sln.G2s);
    G3s = double(sln.G3s);
    G4s = double(sln.G4s);
    G5s = double(sln.G5s);
    
    G6s = -G7s*(paramStruct.epsilon-1)/paramStruct.epsilon;
    
    % make a struct of steady state Lagrange multipliers
    
    std_Multipliers = struct('G1s',G1s,'G2s',G2s,'G3s',G3s,'G4s',G4s,'G5s',G5s,'G6s',G6s,'G7s',G7s,'G8s',G8s,'G9s',G9s);
 
 end