function [z,J] = calcJac(funct,x)
%CALCJAC [nextState, Jacobian] = calcJac(functionHandle, state)
%   This function calculates the state evolution and jacobian of a function

%J = 
% [dz1/dx1 dz1/dx2 ... dz1/dxn
%  dz2/dx1 ...         dz2/dxn
%  ...                  ...
%  dzn/dx1 ...         dzn/dxn]

% So the Jacobian has shape [elements(observations), elements(states)]

z = funct(x);
cols = numel(x);
rows = numel(z);

J = zeros(rows,cols);
h = cols*eps;

for k = 1:cols
    xi = x;
    xi(k) = xi(k) + h*i; %x(i) with a complex part
    J(:,k) = imag(funct(xi))/h; %See differentiation by complex parts
    
end
end

