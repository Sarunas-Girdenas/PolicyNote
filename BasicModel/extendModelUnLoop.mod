// this is simpliefied dynare file

var S X b intr infl q phi C K Z N Kplus bplus Zplus inflplus Cplus Splus Nminus u w;

predetermined_variables b;

varexo shock;

parameters beta eta alpha xi lambda k a z phist omega sigma epsilon Sst Xst bst intrst qst Nst Cst Kst Zst ust taylInfl taylUn wst;


load infl_IntFig;
load u_IntFig;
load eta_IntFig;

@#for i in 1:1331
taylInfl           = infl_IntFig(@{i});
taylUn             = u_IntFig(@{i});
eta                = eta_IntFig(@{i});
stoch_simul( order = 1 );
save results_unemployment@{i} oo_ M_;
@#endfor


model(linear);

// variables for welfare

K(+1)    = Kplus;
b(+1)    = bplus;
Z(+1)    = Zplus;
infl(+1) = inflplus;
C(+1)    = Cplus;
S(+1)    = Splus;
N(-1)    = Nminus;

// eq1 policy, Taylor rule

intr = taylInfl * infl + taylUn * u;

// unemployment

u*ust = -(1-lambda)*Nst*N(-1);

// wage

wst*w=-((z/Xst)*X+bst*b)*eta-((alpha/(1-alpha))*eta*k*beta*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))*(phist+beta)+eta*k*(phist+beta)*(phist+beta*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))))*qst*q + eta*k*qst*((phist+beta)-(phist+beta*(qst^(-alpha/(1-alpha))))*(xi^(alpha/(1-alpha))))*phist*phi;

// eq2 surplus

Sst*S = -(z/Xst)*X-bst*b+(bst/intrst)*b(+1)-(bst/intrst)*intr+(bst/intrst)*infl(+1)+(1-lambda)*beta*Sst*S(+1)-eta*beta*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))*S(+1)*Sst+q*(alpha/(1-alpha))*eta*beta*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))*Sst;

// eq3 interest and borrowing

(bst/intrst)*b(+1) = (bst/intrst)*intr-(bst/intrst)*infl(+1)+phi*phist*(1-eta)*Sst+phist*(1-eta)*Sst*S(+1);

// eq4 borrowing and labour market

-(k/qst)*q = (bst/intrst)*b(+1)-(bst/intrst)*intr+(bst/intrst)*infl(+1)+beta*(1-eta)*Sst*S(+1);

// eq5 Euler equation

-intr = -sigma*C(+1)+sigma*C-infl(+1);

// 3 equations block for Phillips Curve

z*Nst*N = Kst*K + omega*(1/intrst)*intr(+1)*Kst - omega*(1/intrst)*(epsilon-1)*infl(+1)*Kst-omega*(1/intrst)*Kst*K(+1);

-Zst*Z = z*(Nst/Xst)*X - z*(Nst/Xst)*N + omega*(1/intrst)*intr(+1)*Zst-omega*(epsilon/intrst)*infl(+1)*Zst-omega*(1/intrst)*Z(+1)*Zst;

K = Z - (omega/(1-omega))*infl;

// end of 3 equations block

// eq6 Phillips Curve

// beta*infl(+1) = (1/omega)*(1-omega)*(1-omega*beta)*X + infl;

// end of 3 equations block

// eq7 N and consumption

z*Nst*N = Cst*C+a*Nst*N-(1/(1-alpha))*q*k*(qst^(-1/(1-alpha)))*(xi^((2*alpha-1)/(alpha*(1-alpha))))*(1-(1-lambda)*Nst)-k*(qst^(-1/(1-alpha)))*(xi^((2*alpha-1)/(alpha*(1-alpha))))*(1-lambda)*Nst*N(-1);

// eq8 Output and stuff

Nst*N = (1-lambda)*Nst*N(-1) - (alpha/(1-alpha))*q*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))*(1-(1-lambda)*Nst)-(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))*(1-lambda)*Nst*N(-1);

// eq9 phi equation

phi = 0.9*phi(-1) + shock;

end;

shocks;
var shock;
stderr 1;

end;

resid(1);

steady;

check(qz_zero_threshold=1e-15);

stoch_simul(irf=40,order = 1);

conditional_variance_decomposition=1;













































