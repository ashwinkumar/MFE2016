function [ CallPrice ] = HestonCall( S_T,K,v_T,r,rho,sigma,kappa,lambda,theta,t)
% Return the Call Price using Heston Model. 

    I1 = integral(@(u) getReal_Pj(u,S_T,K,v_T,r,rho,sigma,kappa,lambda,theta, t,1),0,500);
    I2 = integral(@(u) getReal_Pj(u,S_T,K,v_T,r,rho,sigma,kappa,lambda,theta, t,2),0,500);
    CallPrice = S_T*(0.5+ I1/pi) - K*exp(-r*t)*(0.5+I2/pi);
end

