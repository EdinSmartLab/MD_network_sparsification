function [ output_args ] = HoleProblemCombined()
%
set(0,'DefaultAxesFontSize',24)
set(groot,'defaultLineLineWidth',2)

for i = 1:14
    epsilon = 1;
    sigma = 1;
    Fmax = 48*(-1/((23/7)^(1/6))^13 + 1/2/((23/7)^(1/6))^7);
    r = Lattice2D(1, [5*i, 5*i], 'sq2');
    if i < 20
        rhole = 5*i/4;
    else
        rhole = 5;
    end
    r = GenerateHole(r, [5*i/2-0.25,5*i/2-0.25], rhole);
    q = CreateChargeList(r);
    %q = abs(q);
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
    Af = Af.*signForceCoul(r);
    Af = abs(Af);
    str = 'Crystal';
    [~,dir] = distancematrix(r,r);
    Afp = Af;
    Afn = Af;
    [n, ~] = size(Af)
    for ii = 1:n
        for j = ii+1:n
            if sign(q(ii)) == sign(q(j))
                Afp(ii,j) = 0;
                Afp(j,ii) = 0;
            else
                Afn(ii,j) = 0;
                Afn(j,ii) = 0;
            end
        end
    end
    Fsp = sparsify_spectral(Afp, 1);
    Fsn = sparsify_spectral(Afn, 1);
    Fs = Fsp - Fsn;
    Fpc = DistanceCut(Afp, Ar, 15);
    Fnc = DistanceCut(Afn, Ar, 15);
    Fpc25 = sparsify_spectral(DistanceCut(Afp, Ar, 2.5),1);
    Fpc50 = sparsify_spectral(DistanceCut(Afp, Ar, 5),1);
    Fpc100 = sparsify_spectral(DistanceCut(Afp, Ar, 10),1);
    Fpc150 = sparsify_spectral(DistanceCut(Afp, Ar, 15),1);
    Fpc250 = sparsify_spectral(DistanceCut(Afp, Ar, 25),1);
    Fnc25 = sparsify_spectral(DistanceCut(Afn, Ar, 2.5),1);
    Fnc50 = sparsify_spectral(DistanceCut(Afn, Ar, 5),1);
    Fnc100 = sparsify_spectral(DistanceCut(Afn, Ar, 10),1);
    Fnc150 = sparsify_spectral(DistanceCut(Afn, Ar, 15),1);
    Fnc250 = sparsify_spectral(DistanceCut(Afn, Ar, 25),1);
    Fc = Fpc - Fnc;
    Fc25 = Fpc25 - Fnc25;
    Fc50 = Fpc50 - Fnc50;
    Fc100 = Fpc100 - Fnc100;
    Fc150 = Fpc150 - Fnc150;
    Fc250 = Fpc250 - Fnc250;
    edge(i) =  nonzeroCount(Fs);
    edgec(i) = nonzeroCount(Fc);

    edgec100(i) = nonzeroCount(Fc100);
    edgec150(i) = nonzeroCount(Fc150);
    edgec250(i) = nonzeroCount(Fc250);  
    NN(i) = length(r(:,1));
end

figure
plot(NN,edge)
hold on
plot(NN, edgec)
plot(NN, edgec100)
plot(NN, edgec150)
plot(NN, edgec250)
plot(NN, NN)
plot(NN, NN.*NN)
xlabel('Number of atoms')
ylabel('Number of edges')
legend('Sparsification', 'Thresholding', 'Combined r_c = 10', 'Combined r_c =15', 'Combined r_c =25', 'N', 'N^2')

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


function [sgn] =signForceCoul(q)
sgn = q*transpose(q);

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