function [ x ,rnd_gen] = Euler_Discretization( x_0,coeff_dt, coeff_dw ,dt,T,type,y,jump_perc)
% Returns the simulated values of x using the Euler Discretization method.
% The process uses the coefficients of dt, dw ,x0, dependent variables y and jump points 
% Type can be 
% a)R - Euler Discretization process to simulate rate movements
% b)S1/S2/S - Euler Discretization process to simulate Stock Price
% movements
    if (~exist('y', 'var'))
        y = zeros(1,T+1);
    end
    if (~exist('jump_perc', 'var'))
         jump_perc = zeros(1,T+1);
    end
    x = zeros;
    x(1) = x_0;
    rnd_gen = zeros(1,T);
    for i=2:T+1
        if jump_perc(i)== 0
            [sigma,rnd_gen(:,i)]= get_sigma_x(coeff_dw,x(i-1),type,y(:,i));
            x(i) = x(i-1)+ dt*get_mu_x(coeff_dt,x(i-1),type,y(:,i)) + sqrt(dt)*sigma;
        else
            x(i) = x(i-1)*(1+jump_perc(1,i));
        end
    end
end

