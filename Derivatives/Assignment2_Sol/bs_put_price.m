function [ put_price ] = put_price(Stock_price,Strike_price, rate,T,sigma,yield)
%Returns the BS Put Option price using the given paremeters
%   Detailed explanation goes here
put_price = Strike_price.*exp(-rate.*T).*normcdf(-getd2(Stock_price, Strike_price, rate,T, sigma)) -Stock_price.*exp(-yield.*T).*normcdf(-getd1(Stock_price, Strike_price, rate,T, sigma));
end