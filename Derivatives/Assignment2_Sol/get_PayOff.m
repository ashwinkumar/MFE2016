function [ PayOff ] = get_PayOff(StockPrice,K,optionType)
% Returns the payoff for different types of options

switch optionType
    %European Call
    case 'EC' 
        PayOff = max(StockPrice- K,0);
    %European Put
    case 'EP'  
        PayOff = max(K -StockPrice,0);
end
end
