function [ PayOff,Exercise ] = get_PayOffFn( StockPrice,K,T,optionType)
% Returns the PayOff function for different options at step T (Final) and
% if the option is exercised
% StockPrice    - Stock Price at given T
% T - current Step
% K     - Strike Price
% optionType - (American/European Call/Put/Straddle) (A/E C/P/S)
switch optionType
    %European Call
    case 'EC' 
        PayOff = max(StockPrice(:,T+1) - K,0);
        %Exercise = (T)*(StockPrice(:,T+1) >= K) + (large)*(K > StockPrice(:,T+1)) ;
        Exercise = StockPrice(:,T+1) >= K;
    %European Put
    case 'EP'  
        PayOff = max(K -StockPrice(:,T+1),0);
        %Exercise = (large)*(StockPrice(:,T+1) >= K) + (T)*(K > StockPrice(:,T+1)) ;
        Exercise = K >=StockPrice(:,T+1);

    %European Straddle
    case 'ES'
        PayOff = abs(K -StockPrice(:,T+1));
        %Exercise = large*ones(T+1,1);
        Exercise = ones(T+1,1);
    %European Binary
    case 'EB'  
        PayOff = StockPrice(:,T+1)>= K;
        Exercise = PayOff ;
    %American Call
    case 'AC' 
        PayOff = max(StockPrice(:,T+1) - K,0);
        %Exercise = (T)*(StockPrice(:,T+1) >= K) + (large)*(K > StockPrice(:,T+1)) ;
        Exercise = StockPrice(:,T+1) >= K;
    %American Put
    case 'AP'  
        PayOff = max(K -StockPrice(:,T+1),0);
        %Exercise = (large)*(StockPrice(:,T+1) > K) + (T)*(StockPrice(:,T+1) <=K) ;
        Exercise = K >=StockPrice(:,T+1);
    %American Straddle
    case 'AS'  
        PayOff = abs(K -StockPrice(:,T+1));
        Exercise = ones(T+1,1);
    otherwise
        warning('Unexpected pay off function.')
end
end

