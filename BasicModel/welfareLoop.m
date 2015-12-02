% Purpose: load simulation results from dynare and calculate the welfare

% INPUT: inputFile   -> File with Dynare Output
%        paramValues -> Struct of Parameter Values

function [ welfareOut ] = welfareLoop( inputFile,paramValues )


	% Set parameter values from the input struct

	beta      = paramValues.beta;    % discount factor
	%eta       = paramValues.eta;     % relative bargaining power (workers)
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
	intrst    = paramValues.intrst;  % steady state interest rate

	% Here we solve for the steady state

	[ std_Variables, std_Multipliers ] = steadyStateFunction( paramValues );

	% Variables

	Sst = std_Variables.Sst;
	bst = std_Variables.bst;
	qst = std_Variables.qst;
	Cst = std_Variables.Cst;
	Nst = std_Variables.Nst;
	Kst = std_Variables.Kst;
	Zst = std_Variables.Zst;

	% Lagrange multipliers

	G1s = std_Multipliers.G1s;
	G2s = std_Multipliers.G2s;
	G3s = std_Multipliers.G3s;
	G4s = std_Multipliers.G4s;
	G5s = std_Multipliers.G5s;
	G6s = std_Multipliers.G6s;
	G7s = std_Multipliers.G7s;
	G8s = std_Multipliers.G8s;
	G9s = std_Multipliers.G9s;

	% Now calculate second moments
	% The following piece of code extracts matrix A (variables) in the same
	% order as variables are declared
	% Let's get the variables list and call it C
	C = inputFile.oo_.dr.order_var;
	% Now add it to our ghx matrix
	A_1 = [C inputFile.oo_.dr.ghx];
	% Now sort the matrix in ascending order
	A_2 = sortrows(A_1);
	% Now delete the first column (which is sorting column)
	A_2(:,1) = [];
	% Call the matrix which is now is exactly the same is in dynare output
	A    = transpose(A_2);
	n_pr = rows(A);
	% Here we take predetermined variables
	k1     = find(inputFile.oo_.dr.kstate(:,2) <= inputFile.M_.maximum_lag+1);
	klag   = inputFile.oo_.dr.kstate(k1, [1,2]);
	state  = [ inputFile.oo_.dr.order_var(klag(:,1)), klag(:,2)-inputFile.M_.maximum_lag-2 ];
	state1 = state(:,1);
	A_sm   = A(:,state1);

	% Let's get the variables list and call it C
	C   = inputFile.oo_.dr.order_var;
	% Now add it to our ghu matrix
	B_1 = [C inputFile.oo_.dr.ghu];
	% Now sort the matrix in ascending order
	B_2 = sortrows(B_1);
	% Now delete the first column (which is sorting column)
	B_2(:,1) = [];
	% Call the matrix which is now the same as in dynare output
	B        = transpose(B_2);
	% Taking predetermined variables
	B_sm     = B(:,state1);
	% computing covariance matrix for predetermined variables
	C1       = (beta^(1/2))*A_sm';
	% omega1=[0.01 0;0 0.01]; % Shock correlation matrix, we call it omega1 because omega already exists
	omega1   = 1; % since we have only one shock
	D_1      = (transpose(B_sm))*omega1*(B_sm)*(beta/(1-beta));
	S_1_0    = zeros(n_pr,n_pr);
	S_1      = D_1;
	H_1      = C1;

	while max(max( abs( S_1-S_1_0 ) )) > 10e-20
    
    	S_1_0 = S_1;
    	S_1   = S_1_0 + H_1 * S_1_0 * H_1';
    	H_1   = H_1^2;
    
	end

	%S_big = beta*A'*S_1*A+(beta/(1-beta))*B'*omega1*B;
    
    S_big = A'*S_1*A+B'*omega1*B;

    % S(1) X(2) b(3) intr(4) infl(5) q(6) phi(7) C(8) K(9) Z(10) N(11) Kplus(12) bplus(13) Zplus(14) inflplus(15) Cplus(16) Splus(17) Nminus(18) u(19)

	welfareOut = (-sigma*(Cst^(1-sigma))*S_big(8,8))+... 
				G1s*( -2*(z/Xst)*S_big(2,2)+2*(bst/intrst)*S_big(13,4)+2*(bst/intrst)*S_big(4,15)-2*(bst/intrst)*S_big(13,15)-2*(bst/intrst)*S_big(4,4)-2*(alpha/(1-alpha))*eta*beta*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))*Sst*S_big(17,6)+(alpha/(1-alpha))*(1/(1-alpha))*eta*beta*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))*Sst*S_big(6,6) )+...        
				G2s*( -2*(bst/intrst)*S_big(13,4)-2*(bst/intrst)*S_big(4,15)+2*(bst/intrst)*S_big(13,15)-2*phist*(1-eta)*Sst*S_big(17,7)+2*(bst/intrst)*S_big(4,4) )+...
				G3s*( 2*(k/qst)*S_big(6,6)+2*(bst/intrst)*S_big(13,4)+2*(bst/intrst)*S_big(4,15)-2*(bst/intrst)*S_big(13,15)-2*bst*intrst*S_big(4,4) )+...
				G4s*( (2/intrst)*S_big(4,4)-2*beta*S_big(15,15)+beta*sigma*(-sigma-1)*S_big(16,16)-sigma*beta*(sigma-1)*S_big(8,8)+2*beta*(sigma^2)*S_big(16,8) )+...
				G5s*( (omega/(1-omega))*Kst*(epsilon*(omega/(1-omega))+epsilon-2)*S_big(5,5)+2*(omega/(1-omega))*Kst*S_big(9,5) )+...
				G6s*( -omega*(epsilon-1)*(epsilon-2)*(Kst/intrst)*S_big(15,15)-2*(omega/intrst)*Kst*S_big(4,4)+2*omega*(epsilon-1)*(1/intrst)*Kst*S_big(4,15)+2*(omega/intrst)*Kst*S_big(4,12)-2*omega*((epsilon-1)/intrst)*Kst*S_big(12,15) )+...          
				G7s*( -omega*epsilon*((epsilon-1)/intrst)*Zst*S_big(15,15)-2*(omega/intrst)*Zst*S_big(4,4)+2*z*(Nst/Xst)*S_big(11,2)-2*z*(Nst/Xst)*S_big(2,2)+2*omega*(epsilon/intrst)*Zst*S_big(15,4)-2*omega*(epsilon/intrst)*Zst*S_big(14,15)+2*(omega/intrst)*Zst*S_big(14,4) )+...
				G8s*( (1/(1-alpha))*((alpha-2)/(1-alpha))*k*(qst^(-1/(1-alpha)))*(xi^((2*alpha-1)/(alpha*(1-alpha))))*S_big(6,6)*(1-(1-lambda)*Nst)-2*(1/(1-alpha))*k*(qst^(-1/(1-alpha)))*(xi^((2*alpha-1)/(alpha*(1-alpha))))*(1-lambda)*Nst*S_big(18,6) )+...
				G9s*( ((-alpha)/(1-alpha))*(1/(1-alpha))*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))*(1-(1-lambda)*Nst)*S_big(6,6)-2*(alpha/(1-alpha))*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))*(1-lambda)*Nst*S_big(18,6) );



end





