function [ c ] = get_mu_x( coeff_dt,x,type,y)
% Returns the coefficient of dt in Euler discretization for different
% simulations. Type can be 
% a)R - Euler Discretization process to simulate rate movements
% b)S1/S2/S - Euler Discretization process to simulate Stock Price
% movements
switch type
    % Call
    case 'R' 
         alpha = coeff_dt(1);
         beta = coeff_dt(2);
         c= alpha*(beta - x);
    %European Put
    case 'S1'  
        c = y(1)*x;
    case 'S2'
        c = y(1)*x;
    case 'S'
        c = coeff_dt(1);
    otherwise
        c = 0;
end

end


