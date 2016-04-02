function [ St, Knocked_Out,Knocked_In ] = GetSimulatedGBMStockPrice(S0,Sb,r,sigma,N,dt,n_sims,optionType,jump_perc)
% Returns the "n_sims" number of Stock Prices simulated using geometric 
%brownian motion and the boolean array of paths in which the stock prices 
%got knocked_out/knocked_in. The function takes the type of crossing 
%(crossing from below or above) and the jump points as parameters

switch optionType
    % Up and Out Call
    case 'U-O-C' 
        c = -1;
    % Up and In Call
    case 'U-I-C'
        c= 1;
    % Down and Out Put
    case 'D-O-P'  
        c = 1;
    % Down and In Put
    case 'D-I-P'
    otherwise
        c = 0;
end

if (~exist('jump_perc', 'var'))
    jump_perc = zeros(1,N+1);
end

a1 = r-sigma^2/2;
St = zeros(n_sims,N+1);
Knocked_Out = ones(n_sims,1);
Knocked_In = zeros(n_sims,1);

randn('seed',0);
for i=1:n_sims
    St(i,1)=S0;
    notknocked = 1;
    for j=1:N
        if jump_perc(1,j)== 0
            St(i,j+1) = St(i,j)*(exp(a1*dt + sigma*sqrt(dt)*randn));
        else
            St(i,j+1) = St(i,j)*(1+jump_perc(1,j));
        end
        if( ((St(i,j+1)-Sb)*c <0) && (notknocked))
            Knocked_Out(i)=0;
            Knocked_In(i)=1;
            notknocked = 0;
        end      
    end
end
end

