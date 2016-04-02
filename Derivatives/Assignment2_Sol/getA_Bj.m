function [ A_j, B_j ] = getA_Bj(u,r,rho,sigma,kappa,lambda,theta, t, j )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
b_j = kappa + lambda -(j==1)*rho*sigma;
u_j = (j==1)*1/2 - (j==2)*1/2;

z1 = complex(-b_j,rho*sigma*u);
z2 = complex(-u.^2, 2*u_j*u);

d_j = sqrt(z1.^2- (sigma^2)*z2);


g_j = (-z1+ d_j)./(-z1 - d_j);

z3 = complex(0,r*u*t);
z5 = (kappa*theta/(sigma^2))*((d_j-z1)*t -2*log( (1- exp(d_j*t).*g_j)./(1-g_j) ));


A_j= z3 + z5;
B_j = 1/(sigma^2)*(d_j-z1).*(1- exp(d_j*t))./(1-g_j.*exp(d_j*t));
end

