function theta_new = theta_update(Z,tau,Y,W,c)
T = diag(tau);

theta_new = (Z'*T*Z)^(-1)*Z*T*(Y-W*c);
end
