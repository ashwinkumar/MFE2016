function [ Option_Price, Exercise ] = get_Option_Price( B, Delta ,S0,K, Exercise,t,optionType)
% This function is used to calculate the option price for constant dividend
% yield.
    p = (S0-K) > (B +  Delta.*S0);
    p1 = abs(S0-K) > (B +  Delta.*S0);
    p2 = (K - S0) > (B +  Delta.*S0);
    switch optionType
        %European Call
        case 'EC' 
            Option_Price = B + Delta.*S0;
        %European Put
        case 'EP' 
            Option_Price = B + Delta.*S0;
        %European Binary
        case 'EB'
            Option_Price = B + Delta.*S0;
        %European Straddle
        case 'ES'
            Option_Price = B + Delta.*S0;
        %American Call
        case 'AC'  
            Option_Price = max(B + Delta.*S0, S0-K);
            Exercise(:,t+1) = p;
        %American Put
        case 'AP'  
            Option_Price = max(B + Delta.*S0, K-S0); 
            Exercise(:,t+1) = p2;
        %American Straddle
        case 'AS'  
            Option_Price = max(B + Delta.*S0, abs(S0-K));
            Exercise(:,t+1) = p1;
        otherwise
            warning('Unexpected pay off function.')
    end
    
    
end

