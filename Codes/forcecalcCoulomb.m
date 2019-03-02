function [force,f] = forcecalc(r, q, epsilon)
%UNTITLED5 Calculates forces for MD simulations, assumes perfect crstal boundaries.
%   Detailed explanation goes here
[rij, dir] = distancematrix(r, r);
size(rij);
size(q*transpose(q));
f = (q*transpose(q))./(rij.*rij+eps);
for i = 1:length(f)
    f(i,i) = 0;
end
[na,~, dim] = size(dir);
force(:,:,1) = f.*dir(:,:,1);
force(:,:,2) = f.*dir(:,:,2);
end

