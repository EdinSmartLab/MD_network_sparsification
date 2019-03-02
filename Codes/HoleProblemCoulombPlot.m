function [ output_args ] = HolePlacements()
%
 set(0,'DefaultAxesFontSize',24)
 set(groot,'defaultLineLineWidth',2)
epsilon = 1;
sigma = 1;
Fmax = 48*(-1/((23/7)^(1/6))^13 + 1/2/((23/7)^(1/6))^7);
r = Lattice2D(1, [40, 40], 'sq2');
r = GenerateEllipse(r, [19.75,19.75], 10,5);
q2 = CreateChargeList2(r);
q = CreateChargeList(r);
q = abs(q);
r = 1.56*r; %Moves close to equilibrium position
% r = r + normrnd(0, .1, size(r));
% natoms = length(r(:,1))
% for i = 1:natoms
%     r = max(max(r))*rand(natoms,2);
% end
% figure
% scatter(r(:,1), r(:,2))
% axis equal
%r = [0.5*r(:,1), 2*r(:,2)];
%scatter(r(:,1), r(:,2))
[Ar,Af, Au] = getACoulomb(r, q, epsilon);
max(max(Af));
min(min(Af));
% input('Continue?')
%Af = Af.*signForceCoul(r);

%[force, ~] = forcecalcCoulomb(r, q, epsilon);
% size(force);
% fx = force(:,:,1);
% fy = force(:,:,2);
sum(sum(Au));
k1 = find(q>0);
k2 = find(q<0);
r;
q;
k1;
k2;
rp1 = r(k1,:);
rp2 = r(k2,:);
figure(96)
scatter(rp1(:,1),rp1(:,2), 100, 'b+')
axis equal
hold on
scatter(rp2(:,1),rp2(:,2), 100, 'r*')
legend('Positive charges', 'Negative Charges')
dir =  '../figures4/'
name = 'AtomPositions';
str = 'CrystalSeparate';
set(figure(96), 'Position', [204   272   763   585])
% saveas(96, strcat(dir, name, str))
% saveas(96, strcat(dir, name, str), 'epsc')
% j1 = find(q2>0);
% j2 = find(q2<0);
% rp1 = r(j1,:);
% rp2 = r(j2,:);
% figure
% scatter(rp1(:,1),rp1(:,2), 100, 'b+')
% axis equal
% hold on
% scatter(rp2(:,1),rp2(:,2), 100, 'r*')


%% Sparisfy Spectral Potential Energy
% sparsifyTest(-Au, Ar, r)
% 
% sparsifyTest(fx, Ar, r)
% 
% sparsifyTest(fy, Ar, r)
% 
% sparsifyTest(sqrt(fx.*fx + fy.*fy), Ar, r)

% min(min(Af))
% max(max(Af))
% min(min(sqrt(fx.*fx + fy.*fy)))
% max(max(sqrt(fx.*fx + fy.*fy)))
% % 
% norm(sqrt(fx.*fx + fy.*fy)-Af)
fmaxSeparateTest(Af, Ar, r, q)


end





function [n] = nonzeroCount(A)
%%Counts the number of nonzero entries in symmetric matrix A excluding the diagonal.
n = 0;
[nx, ny] = size(A);
for i = 1:nx
    for j = i+1:ny
        if A(i,j) ~= 0
            n = n+1;
        end
    end
end

end

function [q] = CreateChargeList(r)
[nx, ny] = size(r);
q = zeros(nx,1);
for i=1:nx
    if mod(floor(2*r(i,1)),2) == 0
        q(i) = 1;
    else
        q(i) = -1;
    end
end
end

function [q] = CreateChargeList2(r)
[nx, ny] = size(r);
q = zeros(nx,1);
for i=1:nx
    if mod(i,2) == 0
        q(i) = 1;
    else
        q(i) = -1;
    end
end
end



function [sgn] =signForceCoul(q)
sgn = q*transpose(q);

end


function [] = fmaxSeparateTest(Af, Ar, r, q)
directory =  '../figures4/'
str = 'Crystal';
[~,dir] = distancematrix(r,r);
Afp = Af;
Afn = Af;
[n, ~] = size(Af)
for i = 1:n
    for j = i+1:n
        if sign(q(i)) == sign(q(j))
            Afp(i,j) = 0;
            Afp(j,i) = 0;
        else
            Afn(i,j) = 0;
            Afn(j,i) = 0;
        end
    end
end

Fmax = max(max(abs(Af)));
sgn = signForceCoul(q);


[rows, ~]  = find(abs(r(:,1)-29.75*1.56)<4*1.56);
[cols, ~]  = find(abs(r(:,2)-20*1.56)<4*1.56);
indices = intersect(rows,cols)

figure(100)
C2 = DistanceCut(Af, Ar, 2.5);
weightedGraphPlot(C2, r(:,:), sum(C2));
axis equal
saveas(100, 'graphCoulombCut')
saveas(100, 'graphCoulombCut', 'epsc')
% xlim([29 32])
% ylim([26 30])
figure(101)
weightedGraphPlot(C2(indices,indices), r(indices,:), sum(C2));
xlim([(19.75+8)*1.56,(19.75+12)*1.56])
ylim([18*1.56,22*1.56])
saveas(101, 'graphCoulombCutZoom', 'epsc')
saveas(101, 'graphCoulombCutZoom')
fs = sgn.*sparsify_spectral(Af, 1.0);
figure(200)
weightedGraphPlot(Af, r, sum(Af));
axis equal
saveas(200, 'graphCoulomb')
saveas(200, 'graphCoulomb', 'epsc')
figure(201)
weightedGraphPlot(Af(indices,indices), r(indices,:), sum(Af));
xlim([(19.75+8)*1.56,(19.75+12)*1.56])
ylim([18*1.56,22*1.56])
saveas(201, 'graphCoulombZoom', 'epsc')
saveas(201, 'graphCoulombZoom')
Fsp = zeros(size(Afp));
Fsn = sparsify_spectral(Afn, 1);
Fs = Fsp - Fsn;
figure(300)
weightedGraphPlot(Fs, r, sum(Fs));
axis equal
saveas(300, 'graphCoulombSparse')
saveas(300, 'graphCoulombSparse', 'epsc')
figure(301)
weightedGraphPlot(Fs(indices,indices), r(indices,:), sum(Fs));
xlim([(19.75+8)*1.56,(19.75+12)*1.56])
ylim([18*1.56,22*1.56])
axis equal
saveas(301, 'graphCoulombSparseZoom', 'epsc')
saveas(301, 'graphCoulombSparseZoom')
% save('HoleProblemCoulombRaPositive')
end