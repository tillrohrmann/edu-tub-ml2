fprintf('Loading Data\n')
load tkcca_toy_data

fprintf('Computing tkCCA\n')
[c,U,V] = tkcca_simple(x,y,lag,kappas);

fprintf('Plotting Results\n')
% plot results
figure(42),clf
subplot(2,3,1:3)
x_ = filter2(U,x,'valid');
y_=V'*y;
t=1:400;
plot(t(11:end-10),x_,'k-o',t,y_,'r-v')
xlim([270 310])
xlabel('Time [samples]','fontsize',6)
legend({'Canonical Component 1', 'Canonical Component 2'},'Location','SouthWest')
set(gca,'fontsize',6)

subplot(2,3,4), 
imagesc(U),colorbar,
ylabel('Channel','fontsize',6)
xlabel('\tau','fontsize',6)
set(gca,'xtick',1:5:(2*lag)+1,'xticklabel',-lag:5:lag,'fontsize',6)
title('Variate of X')
subplot(2,3,5), 
imagesc(V),colorbar
set(gca,'xtick',[],'fontsize',6)
ylabel('Channel','fontsize',6)
title('Variate of Y','fontsize',6)


subplot(2,3,6), 
plot(-lag:lag,c,'r')
xlabel('\tau','fontsize',6)
title('canonical correlogram','fontsize',6)
set(gca,'fontsize',6)

set(gcf,'paperunits','centimeters','papersize',[14 7],'paperposition',[0 0 14 7])
print('tkcca_example.pdf','-dpdf')

