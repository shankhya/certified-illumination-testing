function check_170()
% Independent floating check: original CSV -> polynomial kernels -> exact-root integrals.
root=fileparts(fileparts(mfilename('fullpath')));
s=jsondecode(fileread(fullfile(root,'matlab','fixtures','R170-matlab-inputs.json')));
folder=fullfile(root,'certificates','illumination-testing','data');
names={'CIE_xyz_1931_2deg.csv','CIE_std_illum_A_1nm.csv','CIE_std_illum_D65.csv','CIE_illum_FLs_1nm.csv'};
tables=cell(1,4);
for k=1:4
    data=readmatrix(fullfile(folder,names{k}));
    [ok,ix]=ismember((380:780)',data(:,1));assert(all(ok));tables{k}=data(ix,:);
end
selection=[2 2;3 2;4 3;4 8;4 12;4 4;4 9;4 13];P=zeros(400,3,24);
for light=1:8
    E=tables{selection(light,1)}(:,selection(light,2));
    for channel=1:3
        f=tables{1}(:,channel+1);row=3*(light-1)+channel;
        P(:,:,row)=[E(1:400).*f(1:400),diff(E).*f(1:400)+E(1:400).*diff(f),diff(E).*diff(f)];
    end
    Z=sum(P(:,1,3*light-1)+P(:,2,3*light-1)/2+P(:,3,3*light-1)/3);
    assert(Z>0);P(:,:,3*light-2:3*light)=P(:,:,3*light-2:3*light)/Z;
end
maxNormError=0;U=zeros(1,9);
for k=1:9
    c=s.certificates(k);j=c.coordinate;coeff=P(:,:,15+j);
    for i=1:15,coeff=coeff-c.v(i)*P(:,:,i);end
    val=0;
    for interval=1:400,val=val+absquad(coeff(interval,:));end
    U(k)=val;err=max(abs(val-[c.Ulower,c.Uupper]));maxNormError=max(maxNormError,err);
    assert(val>=c.Ulower-5e-12 && val<=c.Uupper+5e-12);
end
axis=min([s.certificates.margin]./(2*U));assert(abs(axis-s.rhoAxis)<5e-11);
vertexMax=0;
for k=1:size(s.vertices,1)
    point=s.vertices(k,:);values=[s.certificates.base]+[s.certificates.A]*point(1)+2*U*point(2);
    assert(all(values<=s.tau+5e-12));vertexMax=max(vertexMax,max(values));
end
maxEquality=0;maxTargetError=0;
for k=1:5
    w=s.witnesses(k);d=w.nodes(:);response=zeros(1,24);
    assert(max(abs(d))<=1+1e-14 && max(abs(diff(d)))/2<=s.g+1e-14);
    for row=1:24
        a=P(:,1,row);b=P(:,2,row);c=P(:,3,row);z=diff(d);
        response(row)=sum(d(1:400).*a+(d(1:400).*b+z.*a)/2+(d(1:400).*c+z.*b)/3+z.*c/4);
    end
    selected=reshape(3*(w.tests(:)'-1)+(1:3)',1,[]);
    maxEquality=max(maxEquality,max(abs(response(selected))));
    maxTargetError=max(maxTargetError,max(abs(response(16:24)-w.allTargets(:)')));
    assert(max(abs(response(selected)))<5e-13);
    assert(max(abs(response(16:24))) - w.amplitude*s.muAxis > s.tau);
    assert(w.amplitude*s.muAxis<=s.delta);
end
assert(maxTargetError<5e-13);
covered=0;
for mask=0:30
    T=find(bitget(mask,1:5));found=false;
    for k=1:5,found=found||all(ismember(T,s.witnesses(k).tests));end
    assert(found);covered=covered+1;
end
toyU=absquad([-.5,1,0]);assert(abs(toyU-.25)<1e-15);
toyResponse=s.toy.rho/3;assert(abs(toyResponse-s.toy.response)<1e-15 && toyResponse>0);
out=struct('status','PASS','optimizerCalls',0,'rawCIETablesRead',4,'independentQuadraticRootIntegration',true,...
 'maximumResidualNormDiscrepancy',maxNormError,'modelAxisDiscrepancy',abs(axis-s.rhoAxis),...
 'verticesChecked',size(s.vertices,1),'maximumVertexUpper',vertexMax,...
 'maximumFourLightEqualityResidual',maxEquality,'maximumWitnessTargetDiscrepancy',maxTargetError,...
 'excludedProtocols',covered,'residualSubstitutionCounterexamplePassed',true);
fid=fopen(fullfile(root,'verification','matlab-170.json'),'w');assert(fid>=0);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(out,PrettyPrint=true));disp(out)
end

function val=absquad(coeff)
a=coeff(1);b=coeff(2);c=coeff(3);
if c==0
    if b==0,rs=[];else,rs=-a/b;end
else
    rs=roots([c b a]);rs=real(rs(abs(imag(rs))<1e-12));
end
cuts=sort([0;rs(rs>0 & rs<1);1]);val=0;
for k=1:numel(cuts)-1
    lo=cuts(k);hi=cuts(k+1);
    val=val+abs(a*(hi-lo)+b*(hi^2-lo^2)/2+c*(hi^3-lo^3)/3);
end
end
