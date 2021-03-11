seed=16;   % seed 16 gives a nice clear example with 2 distinct sweeps
rng(seed);
Ninit = 2000;
Kmax =  2*Ninit;
mubsave = 1e-5;
sb = .05;
mud = 0.003;
sd = .02;
mun = 0.001;
burnin = 1500;
nsteps = burnin + 2000;
bentimes = []; benss = [];   %keep track of when beneficial mutations occur
% compute N_c, N critical, called Nwave here
if sb==sd
  Nwave = exp(mud/sd)/sqrt(6.4*mubsave*sb)
else
  i = floor(sb/sd);
  ivec = [0:1:i];
  pfixes = 2*(sb - sd*ivec);
  lambda = mud/sd;
  fi = exp(-lambda)*lambda.^ivec./factorial(ivec);
  Nsq = sd/3.2/exp(-lambda)/(mubsave*sum(pfixes.*fi.*(sb - sd*ivec)));
  Nwave = sqrt(Nsq)
end

clear Ntypes Ws N Wbar Wvar totaltypes shannon pdiversity k0 
% note that Ntypes is number of each type, not total number of types
% ICs below if we want to start monomorphic
Ntypes = Ninit;
Ws = exp(Ninit/Kmax);  % then growth will be 1
Wsave = Ws;
nmus = zeros(length(Ntypes),length(Ntypes));

mub=0;
mutotal = mub+mun+mud;
fben = mub/mutotal;
fbd = (mub+mud)/mutotal;

for istep=1:nsteps
  if istep==burnin
    mub=mubsave;
    mutotal = mub+mun+mud;
    fben = mub/mutotal;
    fbd = (mub+mud)/mutotal;
% seed a sweep (we won't do that for this figure)
%    Ntypes(end+1) = 1;
%    Ws(end+1) = Ws(end)*(1+sb); 
%    nmus(end+1,:) = nmus(end,:)+1;
%    nmus(:,end+1) = [ nmus(end,:) 0]';
  end	 % if at end of burnin   
  N(istep) = sum(Ntypes);
%  if (N(istep) ==0 || N(istep)>Kmax)
   if (N(istep)==0)
   fprintf(1,'Population size is %d\n',N(istep));
   N = N(1:istep-1);
   break;
  end
  freqs = Ntypes./N(istep);
  Wbar(istep) = sum(freqs.*Ws);
  Wvar(istep) = sum(freqs.*(Ws - Wbar(istep)).^2);
%  if (istep>burnin && (Wvar(istep)==0 && Wbar(istep)==Wsave))
%     fprintf(1,'Sweep is over at sweep %d\n',istep);
%   break;
%  end
  totaltypes(istep) = sum(ceil(freqs));
  shannon(istep)=0;
  pdiversity(istep)=1;
  k0(istep) = 0;
  for i=1:length(freqs)
    if freqs(i)>0,
	  shannon(istep) = shannon(istep) - freqs(i).*log(freqs(i));
          pdiversity(istep) = pdiversity(istep) - freqs(i).^2;
	  for j=1:length(freqs)
	    k0(istep) = k0(istep) + freqs(i)*freqs(j)*(nmus(i,j));
	  end
   end
  end

% reproduction of each lineage, grouped, Poisson-distributed
% mutation can independently occur with each birth
% either a beneficial or deleterious mutation can occur

  nextNtypes =zeros(size(Ntypes));
  nextW = Ws;
  nextmus = nmus;
%  K = Kmax*(1-exp(-Wfactor*Wbar(istep)));
%  I_N = 1./(1 + exp(-sN*(N(istep)-Kmax)));
  Ricker = exp(-N(istep)/Kmax);
  for iold = 1:length(Ntypes)
    if Ntypes(iold)>0
      Wtmp = Ws(iold).*Ricker;
      new = poissrnd(Ntypes(iold)*Wtmp);
      for inew = 1:new
       r = rand;
        if r<mutotal   %mutation occurs
         nextNtypes(end+1) = 1;
	 nextmus(end+1,:) = nextmus(iold,:)+1;
	 nextmus(:,end+1) = [ nextmus(end,:) 0]';
         if r/mutotal<fben     % it is beneficial
	 stmp = -sb*log(rand); 
        %fprintf(1,'ben: %f\n',stmp);
         benss = [ benss stmp];
         bentimes = [ bentimes istep];
          nextW(end+1) = Ws(iold) +stmp;
%	  nextW(end+1) = Ws(iold)*(1+sb);
         else if r/mutotal<fbd            % it is deleterious
               nextW(end+1) = Ws(iold) + sd*log(rand);
%                nextW(end+1) = Ws(iold)*(1-sd);
              else  %it is neutral
               nextW(end+1) = Ws(iold);
	      end
	 end
      else %mutation did not occur
	      nextNtypes(iold) = nextNtypes(iold)+1;
      end %mutation occurs
    end % loop for new offspring of iold
    end % non-zero Ntypes
  end % loop on iold

  % new becomes the old
    ikeep = find(nextNtypes>0);
    Ntypes = nextNtypes(ikeep);
    Ws = nextW(ikeep);
    nmus = nextmus(ikeep,ikeep);
end % loop on istep  

% compute a lot of stuff out of curiosity, most of which we won't use
Nmean = mean(N(burnin:end));
theta = 2*Nmean*mun;
thetap = 2*Nmean*(mun+mud);
kpred = theta;
Spred = log(theta + 0.5) + 0.5772;
Hpred = theta/(theta +1);
Spredp = log(thetap + 0.5) + 0.5772;
Hpredp = thetap/(thetap +1);
kpredp = theta*exp(-mud/sd);
k0mean = mean(k0(burnin:end));
k0sd = std(k0(burnin:end));
Hmean = mean(pdiversity(burnin:end));
Hsd = std(pdiversity(burnin:end));
Smean = mean(shannon(burnin:end));
Ssd = std(shannon(burnin:end));
k0cov = k0sd/k0mean;
Hcov = Hsd/Hmean;
Scov = Ssd/Smean;

Neff = pdiversity./(1-pdiversity)/2/mutotal;
ikeep =burnin+1:length(N);

L = 2;  F = 14;
figure(1)
subplot(4,1,1)
plot(N(ikeep),'k')
ylabel('N');
%hold on;
%plot([1 length(ikeep)],Nmean*[1 1],'k--');
hold off
set(gca,'fontsize',F);

subplot(4,1,2)
fitness  = Wbar(ikeep)/Wbar(ikeep(1));
fitnessvar = Wvar(ikeep)/Wbar(ikeep(1));
plot(fitness,'b','linewidth',L);
hold on;
plot(fitness+sqrt(fitnessvar),'color',[.9 .9 .9],'linewidth',3);
plot(fitness-sqrt(fitnessvar),'color',[.9 .9 .9],'linewidth',3);
if length(bentimes>0), plot(bentimes-burnin,1,'k+'); end
hold off;
ylabel('Fitness');
set(gca,'fontsize',F);

subplot(4,1,3);
plot(pdiversity(ikeep),'m','linewidth',L);
ylabel('Heterozygosity');
set(gca,'fontsize',F);

subplot(4,1,4);
semilogy(Neff(ikeep),'r','linewidth',L);
h = ylabel('N_e');
set(gca,'fontsize',F);
fprintf('min Neff = %f\n',min(Neff(ikeep)));
set(gca,'YLim',[30 3000]);
set(gca,'YTick',[100 1000]);
xlabel('time (generations)');

%subplot(5,1,5);
%greys = [.9 .95 .98];
%sbs = [.05 .1 .2];
%for i=1:length(sbs)
%   pext = exp(-N(ikeep).*sbs(i).*(1-sqrt(1-2*pdiversity(ikeep))));
%   plot(pext,'color',greys(i)*[1 1 1 ],'linewidth',L);
%hold on;
%end
%ylabel('p_e_x_t');
%hold off;

set(gcf,'color','w');
uisetfont(h);
