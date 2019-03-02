function [A] = DistanceCut(A, Ar, rcut )
%Sparsifies A, based on a cutoff distance specified by Ar.
[nx,ny] = size(A);

for i = 1:nx
    for j = 1:ny
        if Ar(i,j)>rcut
            A(i,j) = 0;
        end
    end
end



end

