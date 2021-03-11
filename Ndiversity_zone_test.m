tic
nreps =25;   % number of replicate runs

Ninits = [800];  % initial population size(s)
                 % can be an array with several values in a row

mubsave = 1e-5;
sb = .05;
mud = 0.01;
sd = .005;
mun = 0.0;
nsteps = 100000;
burnin = 7500;

% first, compute N_c, (N critical), called "Nsafe" here
if sb==sd
  Nsafe = exp(mud/sd)/sqrt(6.4*mubsave*sb)
else
  i = floor(sb/sd);
  ivec = [0:1:i];
  pfixes = 2*(sb - sd*ivec);
  lambda = mud/sd;
  fi = exp(-lambda)*lambda.^ivec./factorial(ivec);
  Nsq = sd/3.2/exp(-lambda)/(mubsave*sum(pfixes.*fi.*(sb - sd*ivec)));
  Nsafe = sqrt(Nsq)
  tmp = lambda.^ivec./factorial(ivec).*(sb-sd*ivec).*(sb-sd*ivec-mud);
  Sigma = sum(tmp);
  Nsq2 = sd*exp(2*lambda)/(6.4*mubsave*Sigma);
  Nsafe2 = sqrt(Nsq2)
end

for iN = 1:length(Ninits)
 Ninit = Ninits(iN);	       
 meltdowns = 0;

 for i=1:nreps
   fprintf(1,' %d ',i);
   Kmax =  2*Ninit;
    
   clear Ntypes Ws N Wbar Wvar totaltypes
% note that Ntypes is number of each type, not total number of types
% ICs below if we want to start monomorphic
Ntypes = Ninit;
Ws = exp(Ninit/Kmax);  % then growth will be 1
Wsave = Ws;
nmus = zeros(length(Ntypes),length(Ntypes));  % number of mutations
                                              % each type is from ancestor 

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
  end	 % if at end of burnin   
  N(istep) = sum(Ntypes);
  if (N(istep)==0)
    fprintf(1,'.',N(istep));
    N = N(1:istep-1);
    meltdowns = meltdowns+1;
    break;
  end
  if (N(istep)==max(2*Ninit,100)) break;
  end
  freqs = Ntypes./N(istep);
  Wbar(istep) = sum(freqs.*Ws);
  Wvar(istep) = sum(freqs.*(Ws - Wbar(istep)).^2);
  totaltypes(istep) = sum(ceil(freqs));

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
%    new = poissrnd(Ntypes(iold)*exp(r*(1-N(istep)/(K*Ws(iold)))));
%    new = poissrnd(Ntypes(iold)*Ws(iold));
%      I_W = 1/(1 + exp(-sW*(Ws(iold)-1)));
      Wtmp = Ws(iold).*Ricker;
      new = poissrnd(Ntypes(iold)*Wtmp);
      for inew = 1:new
       r = rand;
        if r<mutotal   %mutation occurs
         nextNtypes(end+1) = 1;
	 nextmus(end+1,:) = nextmus(iold,:)+1;
	 nextmus(:,end+1) = [ nextmus(end,:) 0]';
         if r/mutotal<fben     % it is beneficial
%	  nextW(end+1) = Ws(iold) - sb*log(rand);
	  nextW(end+1) = Ws(iold)*(1+sb);
         else if r/mutotal<fbd            % it is deleterious
%               nextW(end+1) = Ws(iold) + sd*log(rand);
                nextW(end+1) = Ws(iold)*(1-sd);
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

end % loop on nreps
fprintf(1,'\n');
fprintf(1,'Meltdown prob: %f\n',meltdowns/nreps);
pmelt(iN)= meltdowns/nreps;
end % loop on iN

%Nmean = mean(N(burnin:end));

if (0)   % change 0 to 1 to produce plots
figure(1)
subplot(4,1,1)
plot(N)
ylabel('N');
hold on;
plot([1 length(N)],Nmean*[1 1],'k--');
hold off
subplot(4,1,2);
plot(totaltypes);
ylabel('total types');
subplot(4,1,3)
plot(Wbar);
hold on;
plot(Wbar+sqrt(Wvar),'color',[.9 .9 .9],'linewidth',3);
plot(Wbar-sqrt(Wvar),'color',[.9 .9 .9],'linewidth',3);
hold off;
ylabel('Wbar');
subplot(4,1,4)
hist(Ws,50);
hold off
end   % commented out plot statements

toc
