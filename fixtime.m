
s = sb;
p0 = 1/Ninit;
delx = p0/10;
xs = p0:delx:1-delx;
e2ns = exp(-2*Ninit*s);
e2nsx = exp(-2*Ninit*s*xs);

num=(1-e2nsx).*(1 - exp(-2*Ninit*s*(1-xs)));
denom = s*xs.*(1-xs).*(1-e2ns);
denoma = s*xs.*(1-xs);
term1 = trapz(xs,num./denom);
term1a = trapz(xs,num./denoma);

num2 = exp(-2*Ninit*s*p0)-e2ns;
num2a = exp(-2*s);
denom2 = (1-exp(-2*Ninit*s*p0))*s;
denom2a = (1-exp(-2*s));

x2s = delx:delx:p0;
num3 = (1-exp(-2*Ninit*s*x2s)).^2;
denom3 = x2s.*(1-x2s).*(1-e2ns).*exp(-2*Ninit*s*x2s);
denom3a = x2s.*(1-x2s).*exp(-2*Ninit*s*x2s);
term2 = num2/denom2*trapz(x2s,num3./denom3);
term2a = num2a/denom2a*trapz(x2s,num3./denom3a);

time2fix = term1 + term2
time2fix_approx = term1 + term2a
time2fix_det = 2*log(Ninit-1)/sb

clear hpred;
hpred(1) = 0;
ntmp(1) = 1;
ttf = round(time2fix);
for i=2:ttf
  hpred(i) = 1 - (1-hpred(i-1) + hpred(i-1)/ntmp(i-1))*(1-mun)^2;
  ntmp(i) = ntmp(i-1) + s*ntmp(i-1)*(1-ntmp(i-1)/Ninit);
end

for i=ttf+1:nsteps-burnin
   hpred(i) =  1 - (1-hpred(i-1) + hpred(i-1)/Ninit/(1+2*log(1+sb)))*(1-mun)^2;
end

%figure(3);
%plot([burnin burnin+ttf],[h(burnin) hpred(ttf)],'m-','linewidth',2)
plot(burnin+[1:ttf],hpred(1:ttf),'color','c','linewidth',2);
     plot(burnin+[ttf:nsteps-burnin],hpred(ttf:nsteps-burnin),'m-','linewidth',2)

tmpx = burnin+1:1:burnin+ttf-1;
tmpy = h(burnin) + [hpred(ttf)-h(burnin)]/ttf*(tmpx-burnin);
hfit = [h tmpy  hpred(ttf:nsteps-burnin)];
plot(burnin+ttf,hpred(ttf),'kp','markerfacecolor','k','markersize',14);
