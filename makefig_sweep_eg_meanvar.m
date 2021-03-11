%%  This code got really annoying b/c I accidentally saved every variable in
% the workspace when I previously saved the data files.
% This means that various counters can be overwritten when we read in the
% datafile, so we had to use different variable names here.

% the variable 'h' is assumed and it is stored in some of the early
% files and never gets re-written.
% I'll re-write it when mu changes (using the recursion h = 1-(1-h+h/N)(1-mu)^2 ).
% Also note that N = 5000 is hard-coded.

F = 16;
L = 2;
subplot(1,3,1);
filebasesweepeg = 'data/Ndiv_5k_sb2_mun005_';
nsweepsa = 50;
for iunique=1:nsweepsa
  filename = [ filebasesweepeg int2str(iunique) '.mat'];
  load(filename);
  pdivsa(iunique,:) = pdiversity;
  neffhsa(iunique,:) = 1/2/mun*pdiversity./(1-pdiversity);
end
meanpdiva = mean(pdivsa);
meanneffha = mean(neffhsa);
stdpdiva = std(pdivsa);
stdneffha = std(neffhsa);
handle = fill([1:1:2000 2000:-1:1],[meanpdiva-stdpdiva fliplr(meanpdiva+stdpdiva)],.8*[1 1 1],'EdgeColor',.8*[1 1 1]);
hold on;
sb = .2; mun =0.005;
fixtime;
plot(h,'m','linewidth',3);
plot(burnin+[ttf:nsteps-burnin],hpred(ttf:nsteps-burnin),'m-','linewidth',3)
plot(burnin+ttf,hpred(ttf),'kp','markerfacecolor','k','markersize',14);
plot(meanpdiva,'k--','linewidth',2);
ylabel('Heterozygosity');
set(gca,'fontsize',F);
set(gca,'YLim',[ 0 1.05]);
xlabel('time (generations)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot(1,3,3)
filebasesweepegb = 'data/Ndiv_5k_sb2_mun0001_';
nsweepsb = 50;
for iuniqueb=1:nsweepsb
  filename = [ filebasesweepegb int2str(iuniqueb) '.mat'];
  load(filename);
  pdivsb(iuniqueb,:) = pdiversity;
  neffhsb(iuniqueb,:) = 1/2/mun*pdiversity./(1-pdiversity);
end
ikeep = find(pdivsb(:,10)~=0);
meanpdivb = mean(pdivsb(ikeep,:));
meanneffhb = mean(neffhsb(ikeep,:));
stdpdivb = std(pdivsb(ikeep,:));
stdneffhb = std(neffhsb(ikeep,:));

handle = fill([1:1:nsteps nsteps:-1:1],[meanpdivb-stdpdivb min(fliplr(meanpdivb+stdpdivb),1)],.8*[1 1 1],'EdgeColor',.8*[1 1 1]);
hold on;
sb = .2; mun = 0.0001;
h(1) = 0;
for is = 2:burnin
   h(is) = 1- (1-h(is-1)+h(is-1)/5000)*(1-mun).^2;
end	   
fixtime;
plot(h,'m','linewidth',3);
plot(burnin+[ttf:nsteps-burnin],hpred(ttf:nsteps-burnin),'m-','linewidth',3)
plot(burnin+ttf,hpred(ttf),'kp','markerfacecolor','k','markersize',14);
plot(meanpdivb,'k--','linewidth',2);
ylabel('Heterozygosity');
set(gca,'fontsize',F);
set(gca,'YLim',[ 0 1.05]);
set(gca,'XLim',[0 nsteps]);
xlabel('time (generations)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,2)
filebasesweepegc = 'data/Ndiv_5k_sb05_mun005_';
nsweepsc = 31;
clear pdivsc neffhsc
for iuniquec=1:nsweepsc
  filename = [ filebasesweepegc int2str(iuniquec) '.mat'];
  load(filename);
  pdivsc(iuniquec,:) = pdiversity;
  neffhsc(iuniquec,:) = 1/2/mun*pdiversity./(1-pdiversity);
end
ikeep = find(pdivsc(1:nsweepsc,10)~=0);
meanpdivc = mean(pdivsc(ikeep,:));
meanneffhc = mean(neffhsc(ikeep,:));
stdpdivc = std(pdivsc(ikeep,:));
stdneffhc = std(neffhsc(ikeep,:));

handle = fill([1:1:2000 2000:-1:1],[meanpdivc-stdpdivc min(fliplr(meanpdivc+stdpdivc),1)],.8*[1 1 1],'EdgeColor',.8*[1 1 1]);
hold on;
sb = 0.05; mun = 0.005;
h(1) = 0;
for is = 2:burnin
   h(is) = 1- (1-h(is-1)+h(is-1)/5000)*(1-mun).^2;
end	   
fixtime;
plot(h,'m','linewidth',3);
plot(burnin+[ttf:nsteps-burnin],hpred(ttf:nsteps-burnin),'m-','linewidth',3)
plot(burnin+ttf,hpred(ttf),'kp','markerfacecolor','k','markersize',14);
plot(meanpdivc,'k--','linewidth',2);
ylabel('Heterozygosity');
set(gca,'fontsize',F);
set(gca,'YLim',[ 0 1.05]);
set(gca,'XLim',[0 nsteps]);
xlabel('time (generations)');

set(gcf,'color','w');

