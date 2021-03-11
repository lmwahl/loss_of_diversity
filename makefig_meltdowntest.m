Ns = [ 20:10:100 150 200];
mprob = [ 0.99 0.98 .91 .57 .52 .21 .11 .03  .01 0 0];

% the data "mprob" above are prob(meltdown), obtained by running
% Ndiversity_zone_test with the following parameters:
nreps = 100;
Kmax =  2*Ninit;
mubsave = 1e-4;
sb = .05;
mud = 0.01;
sd = .005;
mun = 0.0;
nsteps = 10000;
burnin = 750;
Nwaveapprox = 187.05;
Nwave = 214.9;

sds = sqrt(mprob.*(1-mprob)./nreps);

semilogx(Ns,mprob,'o','markersize',8,'color',rgb('DarkGreen'),'linewidth',2);
hold on;
plot([Nwave Nwave],[0 1],'--','color',rgb('DarkGreen'),'linewidth',2);
errorbar(Ns, mprob, sds, 'Linestyle','none','color',rgb('ForestGreen'),'linewidth',2);
xlabel('population size, N');
ylabel('probability of meltdown');
text(0.78*Nwave,.95,'N_{c}','color',rgb('DarkGreen'),'fontsize',16);

hold on;
Nwaveapprox = 591.5;
Nwave = 679.6;
mubsave = 1e-5;
nsteps = 1e5;
nreps = 25;

Ns = [ 150 200 300 400 500 600 800];
mprob = [ 1 1 .96 .52 .08 .08  0];

sds = sqrt(mprob.*(1-mprob)./nreps);
semilogx(Ns,mprob,'bo','markersize',8,'linewidth',2);
plot([Nwave Nwave],[0 1],'b--','linewidth',2);
errorbar(Ns, mprob, sds, 'Linestyle','none','linewidth',2);
text(0.78*Nwave,.95,'N_{c}','color','b','fontsize',16);

Nwaveapprox = 93.5;
Nwave = 107.45;
mubsave = 4e-4;
nsteps = 1e5;
nreps = 100;

Ns = [ 10:10:70];
mprob = [ 1 .93 .42 .11 .01 0 0];

sds = sqrt(mprob.*(1-mprob)./nreps);
semilogx(Ns,mprob,'ro','markersize',8,'linewidth',2);
plot([Nwave Nwave],[0 1],'r--','linewidth',2);
errorbar(Ns, mprob, sds, 'Linestyle','none','linewidth',2);
hold off;
text(0.78*Nwave,.95,'N_{c}','color','r','fontsize',16);
axis([10 1000 0 1.]);


set(gca,'fontsize',18)
set(gca,'linewidth',2);
set(gcf,'color','w');
