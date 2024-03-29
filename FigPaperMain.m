clc
clear 
close all
h00=openfig('FigBFinal.fig','reuse');
ylabel('$$\log_{10}(\textrm{BAS})$$');

h0=openfig('Fig1aV2.fig','reuse');
h1 = openfig('Fig1a.fig','reuse'); % open figure

ax1 = gca; % get handle to axes of figure
h2 = openfig('Fig2_text.fig','reuse');
ax2 = gca;
% test1.fig and test2.fig are the names of the figure files which you would % like to copy into multiple subplots
h3 = figure; %create new figure
s1 = subplot(1,2,1); %create and get handle to the subplot axes
s2 = subplot(1,2,2);
fig1 = get(ax1,'children'); %get handle to all the children in the figure
fig2 = get(ax2,'children');
copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
copyobj(fig2,s2);