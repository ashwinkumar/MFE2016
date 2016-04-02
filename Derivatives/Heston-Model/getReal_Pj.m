function [ real_P_j] = getReal_Pj(u,S_T,K,v_T,r,rho,sigma,kappa,lambda,theta, t, j )
  x_T = log(S_T);
 [A_j,B_j]=getA_Bj(u,r,rho,sigma,kappa,lambda,theta,t,j);
 
 z1 = complex(0,u*x_T);
 phi_j = exp(A_j + B_j.*v_T + z1);
 
 z2 = complex(0,-u*log(K));
 z3 = complex(0,u);
 
 real_P_j =  real(exp(z2).*phi_j./z3);
 %real_P_j(isnan(real_P_j)) = 0 ;

 end

