function p = tTest(r,n)

% t= abs(r)./sqrt((1-r.^2)/n);
% p = 2*tcdf(-t,n);

p = betainc(1-r.^2,n/2-1,1/2);