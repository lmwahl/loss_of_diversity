load data/Ndiv_5k_sb2_mun005.mat   % data file
F = 16;
L = 2;
figure(1)
subplot(4,1,1);
plot(N-M,'color',rgb('palegreen'),'Linewidth',3);
hold on;
plot(M,'color','c','Linewidth',3);
plot(N,'k','Linewidth',1);
ylabel('N');
set(gca,'fontsize',F);
set(gca,'YLim',[0 8000])

subplot(4,1,2)
plot(Wbar./Wbar(1),'b','Linewidth',L);
ylabel('Fitness');
set(gca,'YLim',[.9 1.3]);
set(gca,'fontsize',F);

subplot(4,1,3);
grey = 0.75;
load data/Ndiv_5k_sb2_mun005_b.mat;
plot(pdiversity,'color',grey*[1 1 1]);
hold on;
load data/Ndiv_5k_sb2_mun005_c.mat;
plot(pdiversity,'color',grey*[1 1 1]);
load data/Ndiv_5k_sb2_mun005_d.mat;
plot(pdiversity,'color',grey*[1 1 1]);
load data/Ndiv_5k_sb2_mun005_e.mat;
plot(pdiversity,'color',grey*[1 1 1]);
plot(h,'m','linewidth',2);
fixtime;   % annoyingly, some of this code plots stuff we need
ylabel('Heterozygosity');
set(gca,'fontsize',F);
set(gca,'YLim',[ 0 1.1]);

load data/Ndiv_5k_sb2_mun005.mat

subplot(4,1,4);
plot(NeffH,'color',.8*[1 1 1]);
hold on;
h = ylabel('N_e')
xlabel('time (generations)');
set(gca,'fontsize',F);

load data/Ndiv_5k_sb2_mun005_b.mat
NeffH = 1/2/mun*pdiversity./(1-pdiversity);
plot(NeffH,'color',grey*[1 1 1]);
load data/Ndiv_5k_sb2_mun005_c.mat
NeffH = 1/2/mun*pdiversity./(1-pdiversity);
plot(NeffH,'color',grey*[1 1 1]);
load data/Ndiv_5k_sb2_mun005_d.mat
NeffH = 1/2/mun*pdiversity./(1-pdiversity);
plot(NeffH,'color',grey*[1 1 1]);
load data/Ndiv_5k_sb2_mun005_e.mat
NeffH = 1/2/mun*pdiversity./(1-pdiversity);
plot(NeffH,'color',grey*[1 1 1]);

plot(N,'k');
plot(Neffpred,'r','linewidth',2);
set(gca,'Yscale','linear');

set(gcf,'color','w');

uisetfont(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
filebasesweepeg = 'data/Ndiv_5k_sb2_mun005_';
pdivs = zeros(size(pdiversity));
neffhs = zeros(size(NeffH));
nsweeps = 50;
for iunique=1:nsweeps
  filename = [ filebasesweepeg int2str(iunique) '.mat'];
  load(filename);
  pdivs(iunique,:) = pdiversity;
  neffhs(iunique,:) = 1/2/mun*pdiversity./(1-pdiversity);
end
meanpdiv = mean(pdivs);
meanneffh = mean(neffhs);
stdpdiv = std(pdivs);
stdneffh = std(neffhs);

figure(3);
handle = fill([1:1:2000 2000:-1:1],[meanpdiv-stdpdiv fliplr(meanpdiv+stdpdiv)],.8*[1 1 1],'EdgeColor',.8*[1 1 1]);
hold on;
plot(h,'m','linewidth',3);
fixtime;
plot(burnin+[ttf:nsteps-burnin],hpred(ttf:nsteps-burnin),'m-','linewidth',3)
plot(burnin+ttf,hpred(ttf),'kp','markerfacecolor','k','markersize',14);
plot(meanpdiv,'k--','linewidth',2);
ylabel('Heterozygosity');
set(gca,'fontsize',F);
set(gca,'YLim',[ 0 1.1]);
xlabel('time (generations)');
set(gcf,'color','w');

