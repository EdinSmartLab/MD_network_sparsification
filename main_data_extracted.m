clear;clc;close all;

%% figure 2

open('figure2new_a.fig');
a = get(gca,'Children');
xdata = get(a, 'XData');
ydata = get(a, 'YData');
zdata = get(a, 'ZData');
