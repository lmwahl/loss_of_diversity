% currently takes about 160 seconds to run this
tic
muns = [ 1e-5 1e-4 1e-3 ];  % neutral mutation rates to test
for imun = 1:length(muns)
mun = muns(imun);

ss = [ .001 .002 .0025 .003 .004 .005 .007 .01 .02 .03 .04 .05 .07 .1  .2  .5];
Ns = [1e5 1e6 1e7 1e8];
%nadir = zeros(length(Ns),length(ss));

clear fractioNs fractioNsa

for iN = 1:length(Ns)
Ninit = Ns(iN);
for is = 1:length(ss)
s = ss(is);	   
p0 = 1/Ninit;
% first we compute the fixation time using the diffusion approximation
if p0>1e-6
  delx = p0/10;
 else
   delx = p0;
end
xs = p0:delx:1-delx;
e2ns = exp(-2*Ninit*s);
e2nsx = exp(-2*Ninit*s*xs);

num=(1-e2nsx).*(1 - exp(-2*Ninit*s*(1-xs)));
denom = s*xs.*(1-xs).*(1-e2ns);
term1 = trapz(xs,num./denom);

num2 = exp(-2*Ninit*s*p0-e2ns);
denom2 = (1-exp(-2*Ninit*s*p0))*s;

delx = p0/100;
x2s = delx:delx:p0;
num3 = (1-exp(-2*Ninit*s*x2s)).^2;
denom3 = x2s.*(1-x2s).*(1-e2ns).*exp(-2*Ninit*s*x2s);

term2 = num2/denom2*trapz(x2s,num3./denom3);

time2fix(iN,is) = term1 + term2;
time2fix_logistic(iN,is) = 2*log(Ninit-1)/s;
time2fix_exp(iN,is) = log(Ninit)/s;
  
clear hpred ntmp;
hpred(1) = 0;
ntmp(1) = 1;
ttf = round(time2fix(iN,is));
for i=2:ttf
  hpred(i) = 1 - (1-hpred(i-1) + hpred(i-1)/ntmp(i-1))*(1-mun)^2;
  ntmp(i) = ntmp(i-1) + s*ntmp(i-1)*(1-ntmp(i-1)/Ninit);
end
%napprox(iN,is) = 1/(mean(1./ntmp));
%napprox(iN,is) = mean(ntmp);
napprox(iN,is) = Ninit/2;
hpreda(1) = 0;
for i=2:ttf
  hpreda(i) = 1 - (1-hpreda(i-1) + hpreda(i-1)/napprox(iN,is))*(1-mun)^2;
end

nadir(iN,is) = hpred(ttf);
nadira(iN,is) = hpreda(ttf);
nadiraa(iN,is) = 2*mun*2*log(Ninit)/s;
p = (1-mun).^2;
K = 1/napprox(iN,is);
nadira_det(iN,is) = (1-p)*(1-(p*(1-K)).^(ttf-1))./(1 - p*(1-K));

%for i=ttf+1:nsteps-burnin
%   hpred(i) =  1 - (1-hpred(i-1) + hpred(i-1)/Ninit/(1+2*log(1+sb)))*(1-mun)^2;
%end

fractioNsaa(iN,is) = time2fix_logistic(iN,is)/Ninit;
end
end

Neffs = 1/2/(mun).*nadir./(1-nadir);
Neffsa = 1/2/(mun).*nadira./(1-nadira);
for is =1:length(ss)
  fractioNs(:,is) = Neffs(:,is)./Ns';
  fractioNsa(:,is) = Neffsa(:,is)./Ns';
end


newcolors = {'#F00','#F80','#0B0','#00F','#50F','#A0F'};
colororder(newcolors(1:length(Ns)));

M = 8; L = 2;

subplot(1,length(muns),imun)
loglog(ss,(flipud(Neffs))','-','linewidth',L);
xlabel('strength of sweep ($s_b$)','interpreter','latex');
ylabel('$\mathcal{N}_e^*$','interpreter','latex');
legend('N=10^8','N=10^7','N=10^6','N=10^5','location','NE');
set(gca,'XLim',[1e-3 .5])
set(gca,'XTick',[.001 .01 .1 .5]);
set(gca,'fontsize',16);


end % loop on mun
set(gcf,'color','w');

toc

