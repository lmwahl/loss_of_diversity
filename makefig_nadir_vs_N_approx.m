load nadir_vs_N_data
newcolors = {'#F00','#F80','#0B0','#00F','#50F','#A0F'};
colororder(newcolors(1:4));
M=8;
L = 2;

subplot(1,2,1)
  loglog(Ns,Neffs,'-','linewidth',L);
hold on;
loglog(Ns,time2fix_logistic,'o','Markersize',M);
loglog(Ns,time2fix,'^','markersize',M);
xlabel('N','interpreter','latex');
set(gca,'fontsize',20);
set(gca,'XLim',[100 1e7])
ylabel('$\mathcal{N}_e^*$','interpreter','latex');


subplot(1,2,2);
loglog(Ns,fractioNs,'linewidth',L);
hold on;
semilogx(Ns,fractioNsaa,'o','Markersize',M);
xlabel('N','interpreter','latex');
ylabel('$\mathcal{N}_e^*/N$','interpreter','latex');
legend('s_b=.05','s_b=.1','s_b=.2','s_b=.5','location','NE');
set(gca,'fontsize',20);
set(gca,'XLim',[100 1e7])
set(gcf,'color','w');

