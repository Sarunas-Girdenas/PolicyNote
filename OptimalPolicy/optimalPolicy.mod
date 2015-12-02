// this is simpliefied dynare file

var S X b intr infl q phi C K Z N Kplus bplus Zplus inflplus Cplus Splus Nminus G1 G2 G3 G4 G5 G6 G7 G8 G9 w u;

predetermined_variables b;

varexo shock;

parameters beta eta alpha xi lambda k a z phist omega sigma epsilon Sst Xst bst intrst qst Nst Cst Kst Zst G1s G2s G3s G4s G5s G6s G7s G8s G9s wst ust;

model(linear);

// Some variables for welfare

Kplus    = K(+1);
bplus    = b(+1);
Zplus    = Z(+1);
inflplus = infl(+1);
Cplus    = C(+1);
Splus    = S(+1);
Nminus   = N(-1);


//////////////////         BEGIN OPTIMAL POLICY

// 1. FOC w.r.t. N

G6*G6s*z = -G7*G7s*(z/Xst)+(G7s/Xst)*z*X+G8*G8s*z-a*G8*G8s+beta*G8s*G8(+1)*(1-lambda)*k*(qst^(-1/(1-alpha)))*(xi^((2*alpha-1)/(alpha*(1-alpha))))-(1/(1-alpha))*beta*G8s*(1-lambda)*q(+1)*k*(qst^(-1/(1-alpha)))*(xi^((2*alpha-1)/(alpha*(1-alpha))))+G9*G9s-beta*G9(+1)*G9s*(1-lambda)-beta*G9s*G9(+1)*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))*(1-lambda)-(alpha/(1-alpha))*beta*G9s*q(+1)*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))*(1-lambda);                                                                                                                      

// 2. FOC w.r.t. q

(alpha/(1-alpha))*G1*G1s*eta*beta*(qst^(-1/(1-alpha)))*(xi^(alpha/(1-alpha)))*Sst = (alpha/(1-alpha))*(1/(1-alpha))*G1s*eta*beta*q*(qst^(-1/(1-alpha)))*(xi^(alpha/(1-alpha)))*Sst-(alpha/(1-alpha))*G1s*eta*beta*(qst^(-1/(1-alpha)))*(xi^(alpha/(1-alpha)))*Sst*S(+1)-G3s*G3*(k/(qst^2))+2*G3s*(k/(qst^2))*q+(1/(1-alpha))*G8s*G8*k*(qst^((alpha-2)/(1-alpha)))*(xi^((2*alpha-1)/(alpha*(1-alpha))))*(1-(1-lambda)*Nst)+(1/(1-alpha))*((alpha-2)/(1-alpha))*G8s*k*q*(qst^((alpha-2)/(1-alpha)))*(xi^((2*alpha-1)/(alpha*(1-alpha))))*(1-(1-lambda)*Nst)-(1/(1-alpha))*G8s*k*(qst^((alpha-2)/(1-alpha)))*(xi^((2*alpha-1)/(alpha*(1-alpha))))*(1-lambda)*Nst*N(-1)+(alpha/(1-alpha))*G9s*G9*(qst^(-1/(1-alpha)))*(xi^(alpha/(1-alpha)))*(1-(1-lambda)*Nst)-(alpha/(1-alpha))*(1/(1-alpha))*G9s*q*(qst^(-1/(1-alpha)))*(xi^(alpha/(1-alpha)))*(1-(1-lambda)*Nst)-(alpha/(1-alpha))*G9s*(qst^(-1/(1-alpha)))*(xi^(alpha/(1-alpha)))*(1-lambda)*Nst*N(-1);

// 3. FOC w.r.t. pi

G1(-1)*G1s*bst = G1s*bst*intr(-1)-G1s*bst*b-G1s*bst*infl+G2(-1)*G2s*bst-G2s*bst*intr(-1)+G2s*bst*b+G2s*bst*infl-G3(-1)*G3s*bst+G3s*bst*intr(-1)-G3s*bst*b-G3s*bst*infl+G4s*G4(-1)-sigma*G4s*C+sigma*G4s*C(-1)-G4s*infl+G5*G5s*(omega/(1-omega))*Kst+G5s*(omega/(1-omega))*Kst*K+G5s*(omega/(1-omega))*Kst*((epsilon-1)+epsilon*omega/(1-omega))*infl-G6(-1)*G6s*omega*(epsilon-1)*Kst-G6s*omega*((epsilon-1)^2)*infl*Kst-G6s*omega*(epsilon-1)*K*Kst+G6s*omega*(epsilon-1)*Kst*intr(-1)-G7s*G7(-1)*omega*epsilon*Zst-G7s*omega*Zst*Z*epsilon-G7s*omega*(epsilon^2)*infl*Zst+G7s*omega*epsilon*Zst*intr(-1);

// 4. FOC w.r.t. C

-(1-sigma)*(Cst^(1-sigma))*C = -beta*sigma*G4s*G4-beta*(sigma^2)*G4s*C+beta*(sigma^2)*G4s*C(+1)+beta*sigma*G4s*infl(+1)+sigma*G4s*G4(-1)-(sigma^2)*G4s*C+(sigma^2)*G4s*C(-1)-G4s*sigma*infl-G8s*G8*Cst-G8s*Cst*C;

// 5. FOC w.r.t. S

-G1*G1s = -G1(-1)*G1s*(1-lambda)-(alpha/(1-alpha))*G1s*eta*q*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))+G1s*G1(-1)*eta*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))-G2s*G2(-1)*phist*((1-eta)/beta)-G2s*phist*phi(-1)*((1-eta)/beta)-G3(-1)*G3s*(1-eta);

// 6. FOC w.r.t. X

-G1s*G1*z = G7s*z*Nst*G7+G7s*z*Nst*N;

// 7. FOC w.r.t. b

-G1s*G1 = intr(-1)*(G1s-G2s+G3s)+infl*(-G1s+G2s-G3s)-G1s*G1(-1)+G2s*G2(-1)-G3s*G3(-1);

// 8. FOC w.r.t. i

-G1s*G1*(intrst/bst) = G1s*(intrst/bst)*b(+1)-G1s*(intrst/bst)*intr+G1s*(intrst/bst)*infl(+1)-G2s*G2*(intrst/bst)-G2s*(intrst/bst)*b(+1)-G2s*(intrst/bst)*infl(+1)+G2s*(intrst/bst)*intr+G3*G3s*(intrst/bst)+G3s*(intrst/bst)*b(+1)-G3s*(intrst/bst)*intr+G3s*(bst/intrst)*infl(+1)-G4s*G4*beta+G4s*beta*intr+G6s*G6*omega*(1/intrst)*Kst-G6s*omega*beta*Kst*intr+G6s*omega*beta*(epsilon-1)*infl(+1)*Kst+G6s*omega*beta*Kst*K(+1)+G7s*omega*epsilon*beta*infl(+1)*Zst-G7s*omega*beta*intr*Zst+G7s*omega*beta*Zst*Z(+1)+G7*G7s*omega*beta*Zst;

// 9. FOC w.r.t. K

-G5s*G5 = G5s*(omega/(1-omega))*infl+G6s*G6-G6(-1)*G6s*omega+G6s*omega*intr(-1)-G6s*omega*(epsilon-1)*infl;

// 10. FOC w.r.t. Z

G5s*G5*(epsilon/(epsilon-1)) = G7*G7s-G7s*G7(-1)*omega-G7s*omega*epsilon*infl+G7s*omega*intr(-1);

///////////////////          END OPTIMAL POLICY


// unemployment

u*ust = -(1-lambda)*Nst*N(-1);

// wage

wst*w=-((z/Xst)*X+bst*b)*eta-((alpha/(1-alpha))*eta*k*beta*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))*(phist+beta)+eta*k*(phist+beta)*(phist+beta*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))))*qst*q + eta*k*qst*((phist+beta)-(phist+beta*(qst^(-alpha/(1-alpha))))*(xi^(alpha/(1-alpha))))*phist*phi;

// eq1 policy, Taylor rule

// intr = taylInfl * infl + taylUn * u;

// eq2 surplus

Sst*S = -(z/Xst)*X-bst*b+(bst/intrst)*b(+1)-(bst/intrst)*intr+(bst/intrst)*infl(+1)+(1-lambda)*beta*Sst*S(+1)-eta*beta*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))*S(+1)*Sst+q*(alpha/(1-alpha))*eta*beta*(qst^(-alpha/(1-alpha)))*(xi^(alpha/(1-alpha)))*Sst;

// eq3 interest and borrowing

(bst/intrst)*b(+1) = (bst/intrst)*intr-(bst/intrst)*infl(+1)+phi*phist*(1-eta)*Sst+phist*(1-eta)*Sst*S(+1);

// eq4 borrowing and labour market

-(k/qst)*q = (bst/intrst)*b(+1)-(bst/intrst)*intr+(bst/intrst)*infl(+1)+beta*(1-eta)*Sst*S(+1);

// eq5 Euler equation

-intr = -sigma*C(+1)+sigma*C-infl(+1);

// eq6 Phillips Curve

//beta*infl(+1) = (1/omega)*(1-omega)*(1-omega*beta)*X+infl;

// 3 equations block for Phillips Curve

z*Nst*N = Kst*K + omega*(1/intrst)*intr*Kst - omega*(1/intrst)*(epsilon-1)*infl(+1)*Kst-omega*(1/intrst)*Kst*K(+1);

-Zst*Z = z*(Nst/Xst)*X - z*(Nst/Xst)*N + omega*(1/intrst)*intr*Zst-omega*(epsilon/intrst)*infl(+1)*Zst-omega*(1/intrst)*Z(+1)*Zst;

//-(omega/(1-omega))*infl = -K + Z;

K = Z - (omega/(1-omega))*infl;

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

// check;

check;

stoch_simul(irf=40,order = 1);

conditional_variance_decomposition=1;

save inputFile oo_ M_;













































