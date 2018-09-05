more off; format short; format compact;
graphics_toolkit('gnuplot');
%
fname='hpclin3-time.csv';
%
system('mkdir -p tmp');
system(['cut -b -8  ',fname,' > tmp/ksptype.tmp']);
system(['cut -b 19- ',fname,' > tmp/right.tmp']);

csv=dlmread('tmp/right.tmp');
s=fileread('tmp/ksptype.tmp');

COL_RTOL=1; COL_ATOL=2; COL_DTOL=3; COL_NORM=4; COL_MAXI=5;
COL_ITER=6; COL_REAS=7; COL_TIME=8;

ksp=reshape(s,9,numel(s)/9)'; ksp=ksp(:,1:8);

f=find(csv(:,COL_REAS)<0);
csv(f,[COL_ITER,COL_TIME])=nan;

list_rtol=unique(csv(:,COL_RTOL)),
for(r=1:numel(list_rtol));
   f=find(csv(:,COL_RTOL)==list_rtol(r));
   sec(:,r)=csv(f,COL_TIME);
   its(:,r)=csv(f,COL_ITER);
end;

ok=find(~any(isnan(its')));
sec=sec(ok,:); its=its(ok,:); ksp=ksp(ok,:);
[ksp,num2str(sec),int2str(its)]

% sort by best time to rtol<1e-5
[~,s]=sort(sec(:,1)); np=numel(s);

fs=14;
%
figure(1); clf; hold on; grid off;
set (gca, "ygrid", "on");
%
title('PETSc Iterative Solver Single-core Performance','fontsize',fs);
ylabel('Time (sec)','fontsize',fs);
xlabel('Krylov Subspace Methods','fontsize',fs);
set(gca,'fontsize',8);
set(gca,'xtick',1:np);
set(gca,'xticklabel',upper(ksp(s,:)));
axis([0.5,np+0.5  , 0,200]);
%
bar(sec(s,:));
hleg=legend(num2str(list_rtol));

set(hleg,'fontsize',12);
legend('location','northwest');

if(false);
  x=0;
  for(ip=s);x+=1;%1:numel(p));
    text(x,sec(ip)+1.0,['[',strtrim(int2str(its(ip))),']'],...
	'horizontalalignment','center','verticalalignment','bottom');
  end;
end;
set(gcf,'paperposition',[0.25,0.25 , 8,6]);
print(['petsc-ksp-time','.eps'],'-depsc2','-FHelvetica');
%