function [ c ,rnd_gen] = get_sigma_x( coeff_dw,x,type,y)
% Returns the coefficient of dw and random number used to generate it
% in Euler discretization for different simulations. 
% Type can be 
% a)R - Euler Discretization process to simulate rate movements
% b)S1/S2/S - Euler Discretization process to simulate Stock Price
% movements
rnd_gen = zeros;
switch type
    % Call
    case 'R' 
         delta = coeff_dw(1);
         c= delta*sqrt(x)*randn;
         rnd_gen = randn;
    %European Put
    case 'S1' 
        sigma11 = coeff_dw(1);
        sigma12 = coeff_dw(2);
        %rnd_gen(1) = randn;
        rnd_gen = randn;
        c = sigma11*sqrt(x)*y(2) + sigma12*x*rnd_gen;
    case 'S2'
        sigma21 = coeff_dw(1);
        S1 = y(2);
        c = sigma21*(S1-x)*y(3);
    case 'S'
        c = coeff_dw(1)*randn;
    otherwise
        c = 0;
end

end

