function out = znormalise (value, mean, omega)

%Performs znormalization based on this equation:
% value-mean/omega

out = (value - mean)/omega;



end