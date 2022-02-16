function tau_new = tau_update(Pi,simga_sq,C,Theta,Y,W,Z)
q = size(C,2);
numerator = zeros(length(Y),q);
for ii = 1:q
    numerator(:,ii) = Pi(ii).*exp(-(Y-W*C(:,ii)-Z*Theta(:,ii)).^2./simga_sq);
end
SUM = sum(numerator,2);
denominator = kron(ones(1,q),SUM);
tau_new = numerator./denominator;