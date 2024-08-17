function dydt = RetroM_C(t,y,D,gamma,alpha,beta,theta,mu_1max,mu_2max,K_S,Y_XS,Y_PS,S0)
dydt=zeros(4,1);

X1 = y(1);
X2 = y(2);
P = y(3);
S = y(4);

mu1 = mu_1max*S/(K_S+S)*(1+gamma*P/S0/Y_PS);
mu2 = mu_2max*S/(K_S+S);

dydt(1) = mu1*X1*(1-theta)-D*X1;
dydt(2) = mu2*X2 + theta*mu1*X1-D*X2;
dydt(3) = (alpha*mu1+beta)*X1-D*P;
dydt(4) = D*(S0-S)-(mu1*X1+mu2*X2)/Y_XS-(alpha*mu1+beta)*X1/Y_PS;

end