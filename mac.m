function mac=mac(Phi1,Phi2)
% This function calculates mac between phi1 and phi2
mac= (abs(Phi1'*Phi2))^2/((Phi1'*Phi1)*(Phi2'*Phi2));
end