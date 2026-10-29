function make_figures()
% Publication figures from frozen evidence; no optimization or new experiment.
root=fileparts(fileparts(mfilename('fullpath')));
outdir=fullfile(root,'figures');if ~isfolder(outdir),mkdir(outdir);end
s=jsondecode(fileread(fullfile(root,'matlab','fixtures','R174-matlab-inputs.json')));
n=jsondecode(fileread(fullfile(root,'matlab','fixtures','R170-matlab-inputs.json')));
sp=jsondecode(fileread(fullfile(root,'matlab','fixtures','R171-matlab-inputs.json')));
set(groot,'defaultAxesFontName','Arial','defaultTextFontName','Arial');
knots=unique([[s.safePieces.left],[s.safePieces.right],[s.failurePieces.left],[s.failurePieces.right]]);
lo=zeros(size(knots));hi=lo;
for k=1:numel(knots)
 t=knots(k);a=s.safePieces(arrayfun(@(p)p.left<=t+1e-18 && t<=p.right+1e-18,s.safePieces));
 b=s.failurePieces(arrayfun(@(p)p.left<=t+1e-18 && t<=p.right+1e-18,s.failurePieces));
 lo(k)=min([a.a]+t*[a.b]);hi(k)=min([b.a]+t*[b.b]);
end
assert(abs(max(hi-lo)-s.maximumGap)<1e-10);
f=figure('Visible','off','Color','w','Units','inches','Position',[1 1 3.45 2.7]);
ax=axes(f,'Position',[.20 .21 .76 .75]);hold(ax,'on');
fill(ax,[knots fliplr(knots)]*1000,[lo fliplr(hi)],[.85 .90 .94],'EdgeColor','none','HandleVisibility','off');
plot(ax,1000*knots,lo,'-','Color',[.04 .22 .42],'LineWidth',1.3);
plot(ax,1000*knots,hi,'--','Color',[.65 .17 .10],'LineWidth',1.3);
xlabel(ax,'Pair tolerance \delta (10^{-3})');ylabel(ax,'Critical slope G(\delta) (nm^{-1})');
xlim(ax,[0 2]);ylim(ax,[.017 .041]);grid(ax,'on');box(ax,'on');set(ax,'FontSize',8);
legend(ax,{'Certified lower bound','Witness upper bound'},'Location','northeast','FontSize',7,'Box','off');
savefigs(f,outdir,'frontier');

f=figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.05 2.7]);
tl=tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');
ax=nexttile(tl);hold(ax,'on');
fill(ax,[0 n.rhoAxis 0]*1000,[0 0 n.muAxis]*1e6,[.85 .90 .94],'EdgeColor',[.04 .22 .42],'LineWidth',1.1);
plot(ax,[5 3],[0 20],'ko','MarkerFaceColor','k','MarkerSize',4);
xlabel(ax,'Model distance \rho (10^{-3})');ylabel(ax,'Kernel L^1 bound \mu (10^{-6})');
title(ax,'(a) Normalized-kernel uncertainty','FontWeight','normal');
xlim(ax,[0 6]);ylim(ax,[0 115]);grid(ax,'on');box(ax,'on');set(ax,'FontSize',8);
ax=nexttile(tl);hold(ax,'on');
fill(ax,[0 sp.rhoAxis 0]*1000,[0 0 sp.kappaAxis]*1e6,[.85 .90 .94],'EdgeColor',[.04 .22 .42],'LineWidth',1.1);
plot(ax,[2 3],[50 20],'kd','MarkerFaceColor','k','MarkerSize',4);
xlabel(ax,'Model distance \rho (10^{-3})');ylabel(ax,'Relative SPD shape bound \kappa (ppm)');
title(ax,'(b) Shared-SPD uncertainty','FontWeight','normal');
xlim(ax,[0 6]);ylim(ax,[0 115]);grid(ax,'on');box(ax,'on');set(ax,'FontSize',8);
savefigs(f,outdir,'stability');

idx=find(arrayfun(@(w)~ismember(2,w.tests),n.witnesses));assert(isscalar(idx));w=n.witnesses(idx);
d=w.nodes(:);lambda=(380:780)';rplus=(1+d)/2;rminus=(1-d)/2;
assert(min([rplus;rminus])>=-1e-15 && max([rplus;rminus])<=1+1e-15);
assert(max(abs(diff(d)))/2<=n.g+1e-14);
P=readkernels(fullfile(root,'certificates','illumination-testing','data'));
resp=zeros(1,24);
for row=1:24
 a=P(:,1,row);b=P(:,2,row);c=P(:,3,row);z=diff(d);
 resp(row)=sum(d(1:400).*a+(d(1:400).*b+z.*a)/2+(d(1:400).*c+z.*b)/3+z.*c/4);
end
selected=reshape(3*(w.tests(:)'-1)+(1:3)',1,[]);assert(max(abs(resp(selected)))<5e-13);
resp(selected)=0;mis=max(abs(reshape(resp,3,8)),[],1);
f=figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.05 2.65]);
tl=tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');
ax=nexttile(tl);hold(ax,'on');
plot(ax,lambda,rplus,'-','Color',[.04 .22 .42],'LineWidth',1.1);
plot(ax,lambda,rminus,'--','Color',[.65 .17 .10],'LineWidth',1.1);
xlim(ax,[380 780]);ylim(ax,[0 1]);xlabel(ax,'Wavelength (nm)');ylabel(ax,'Reflectance');
title(ax,'(a) Exact four-test witness pair','FontWeight','normal');
legend(ax,{'r_+','r_-'},'Location','southeast','FontSize',7,'Box','off');set(ax,'FontSize',8);grid(ax,'on');box(ax,'on');
ax=nexttile(tl);hold(ax,'on');
b=bar(ax,1:8,mis,.65,'FaceColor','flat');
b.CData=repmat([.60 .65 .70],8,1);b.CData(2,:)=[.75 .44 .17];b.CData(6:8,:)=repmat([.04 .22 .42],3,1);
yline(ax,.02,'k--','Limit 0.02','LabelHorizontalAlignment','left','FontSize',7);
xticks(ax,1:8);xticklabels(ax,{'A','D65*','FL2','FL7','FL11','FL3','FL8','FL12'});xtickangle(ax,35);
ylabel(ax,'Maximum XYZ difference');title(ax,'(b) Test and deployment responses','FontWeight','normal');
set(ax,'FontSize',8);grid(ax,'on');box(ax,'on');
text(ax,.99,.96,'* omitted test','Units','normalized','HorizontalAlignment','right','FontSize',7);
savefigs(f,outdir,'witness');
report=struct('status','PASS','optimizerCalls',0,'frontierGap',max(hi-lo),...
 'witnessMaximumSlope',max(abs(diff(d)))/2,'witnessDeploymentMismatch',max(mis(6:8)),...
 'sourceFiles',{{'R174-matlab-inputs.json','R170-matlab-inputs.json','R171-matlab-inputs.json'}});
fid=fopen(fullfile(outdir,'generation.json'),'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));fclose(fid);disp(report);
end
function savefigs(f,folder,name)
exportgraphics(f,fullfile(folder,[name '.pdf']),'ContentType','vector');
exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',220);close(f);
end
function P=readkernels(folder)
names={'CIE_xyz_1931_2deg.csv','CIE_std_illum_A_1nm.csv','CIE_std_illum_D65.csv','CIE_illum_FLs_1nm.csv'};
tables=cell(1,4);
for k=1:4
 data=readmatrix(fullfile(folder,names{k}));[ok,ix]=ismember((380:780)',data(:,1));assert(all(ok));tables{k}=data(ix,:);
end
selection=[2 2;3 2;4 3;4 8;4 12;4 4;4 9;4 13];P=zeros(400,3,24);
for light=1:8
 E=tables{selection(light,1)}(:,selection(light,2));
 for channel=1:3
  c=tables{1}(:,channel+1);row=3*(light-1)+channel;
  P(:,:,row)=[E(1:400).*c(1:400),diff(E).*c(1:400)+E(1:400).*diff(c),diff(E).*diff(c)];
 end
 Z=sum(P(:,1,3*light-1)+P(:,2,3*light-1)/2+P(:,3,3*light-1)/3);
 assert(Z>0);P(:,:,3*light-2:3*light)=P(:,:,3*light-2:3*light)/Z;
end
end
