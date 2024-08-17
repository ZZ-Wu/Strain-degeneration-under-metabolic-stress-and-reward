function y = eqn_for_numeric_Jam(x) % only for one-step retro-mutation
% parameters that can be modified from te input
alpha=x(5); beta=x(6); theta=x(7); gamma=x(8); D=x(9);
mu_1max = x(10); mu_2max = x(11); K_S=x(12); Y_XS=x(13); Y_PS=x(14); S0=x(15);

X1 = x(1); X2 = x(2); P = x(3); S = x(4);
y = zeros(1,4);

% mu1 = mu_max*S/(K_S+S)*(1+gamma*P/S0/Y_PS);
% mu2 = mu_max*S/(K_S+S);

y(1) = (mu_1max*S/(K_S+S)*(1+gamma*P/S0/Y_PS))*X1*(1-theta)-D*X1;
y(2) = (mu_2max*S/(K_S+S))*X2 + theta*(mu_1max*S/(K_S+S)*(1+gamma*P/S0/Y_PS))*X1-D*X2;
y(3) = (alpha*(mu_1max*S/(K_S+S)*(1+gamma*P/S0/Y_PS))+beta)*X1-D*P;
y(4) = D*(S0-S)-((mu_1max*S/(K_S+S)*(1+gamma*P/S0/Y_PS))*X1+(mu_2max*S/(K_S+S))*X2)/Y_XS-(alpha*(mu_1max*S/(K_S+S)*(1+gamma*P/S0/Y_PS))+beta)*X1/Y_PS;
end