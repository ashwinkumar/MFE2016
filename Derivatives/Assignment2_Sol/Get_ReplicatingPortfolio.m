function [ V_T,B_T,Delta_T ] = Get_ReplicatingPortfolio(V0,S_T,r,dt,d1,N,cost,type)
% Return the replicating portfolio for the simulated stock price 
% V_T - Value of the replicating portfolio
% B_T - Vaule in bond
% Delta_T - Number of shares
V_T = zeros(1,N+1);
B_T = zeros(1,N+1);
switch(type)
    case 'C'
        Delta_T = normcdf(d1);
    case 'P'
        Delta_T = normcdf(d1)-1;
    otherwise
        warning('Unexpected option type.')
end
V_T(1) = V0;
B_T(1) = V_T(1)- Delta_T(1)*S_T(1);

for i=2:N+1
    delta_change = abs(Delta_T(i)-Delta_T(i-1));
    V_T(i) = Delta_T(i-1)*S_T(i) + B_T(i-1)*exp(r(1)*dt) - delta_change*S_T(i)*cost/10^4;
    B_T(i)=  V_T(i) - Delta_T(i)*S_T(i);
end

end

