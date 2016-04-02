function [ S_T,v_T ] = CIRProcess( S0,v_T0,mu,sigma,rho,kappa, theta,T,dt )
% Returns the stock price and volatility simulated using CIR process
S_T = zeros(1,T+1);
v_T = zeros(1,T+1);
v_T(1) = v_T0;
S_T(1) = S0;
for i=2:T+1
    dz1 = randn;
    n2 = randn;
    dz2 = (rho)*dz1 + sqrt(1-rho^2)*n2;
    S_T(i) = S_T(i-1)+ S_T(i-1)* ( mu*dt + sqrt(v_T(i-1)*dt)*dz1 );
    v_T(i) = v_T(i-1) + kappa*(theta- v_T(i-1))*dt + sigma*dz2*sqrt(v_T(i-1)*dt);
end

