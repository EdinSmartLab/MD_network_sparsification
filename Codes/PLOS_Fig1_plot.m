clear;clc;close all
load Af.mat;
load Ar.mat;
load HoleProblemCoulombData.mat;
edgee = edge-edgee;
edgece = edgec-edgece;
[rows, ~]  = find(abs(r(:,1)-29*1.56)<3*1.56);
[cols, ~]  = find(abs(r(:,2)-20*1.56)<3*1.56);
indices = intersect(rows,cols);


C2 = DistanceCut(Af, Ar, 2.5);
Fs = sparsify_spectral(Af, 1.0);

figure(1);
weightedGraphPlot(Af, r, sum(Af));
axis equal
%print('-depsc', 'graphCoulomb.eps');close;

figure(2);
weightedGraphPlot(Af(indices,indices), r(indices,:), sum(Af));
axis equal
%print('-depsc', 'graphCoulomb_zoom.eps');close;

figure(3);
weightedGraphPlot(C2(:,:), r(:,:), sum(C2));
axis equal
%print('-depsc', 'graphCoulombCut.eps');close;

figure(4);
weightedGraphPlot(C2(indices,indices), r(indices,:), sum(C2));
axis equal
%print('-depsc', 'graphCoulombCut_zoom.eps');close;

figure(5);
weightedGraphPlot(Fs(:,:), r(:,:), sum(Fs));
axis equal
%print('-depsc', 'graphCoulombSparse.eps');close;

figure(6);
weightedGraphPlot(Fs(indices,indices), r(indices,:), sum(Fs));
axis equal
%print('-depsc', 'graphCoulombSparse_zoom.eps');close;


% C2 = DistanceCut(Af, Ar, 2.5);
% weightedGraphPlot(C2(indices,indices), r(indices,:), sum(C2));
% axis equal
% saveas(100, 'graphCoulombCut', 'epsc')
% % xlim([29 32])
% % ylim([26 30])
% xlim([27*1.56,31*1.56])
% ylim([18*1.56,22*1.56])
% saveas(100, 'graphCoulombCutZoom', 'epsc')
% fs = sgn.*sparsify_spectral(Af, 1.0);
% figure(200)
% weightedGraphPlot(Af(indices,indices), r(indices,:), sum(Af));
% axis equal
% saveas(200, 'graphCoulomb')
% saveas(200, 'graphCoulomb', 'epsc')
% xlim([27*1.56,31*1.56])
% ylim([18*1.56,22*1.56])
% % saveas(200, 'graphCoulombZoom', 'epsc')
% figure(300)
% weightedGraphPlot(Fs(indices,indices), r(indices,:), sum(Fs));
% axis equal
% saveas(300, 'graphCoulombSparse')
% saveas(300, 'graphCoulombSparse', 'epsc')
% xlim([27*1.56,31*1.56])
% ylim([18*1.56,22*1.56])
% saveas(300, 'graphCoulombSparseZoom', 'epsc')