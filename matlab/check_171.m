function check_171()
root=fileparts(fileparts(mfilename('fullpath')));
s=jsondecode(fileread(fullfile(root,'matlab','fixtures','R171-matlab-inputs.json')));
P=readkernels(fullfile(root,'certificates','illumination-testing','data'));
maxBlockError=0;S=zeros(1,9);
for k=1:9
    cert=s.certificates(k);j=cert.coordinate;
    val=sum(P(:,1,15+j)+P(:,2,15+j)/2+P(:,3,15+j)/3);
    for light=1:5
        coeff=zeros(400,3);
        for channel=1:3
            row=3*(light-1)+channel;coeff=coeff+cert.v(row)*P(:,:,row);
        end
        block=0;
        for interval=1:400,block=block+absquad(coeff(interval,:));end
        assert(block>=cert.blocks(light).lower-5e-12 && block<=cert.blocks(light).upper+5e-12);
        maxBlockError=max(maxBlockError,max(abs(block-[cert.blocks(light).lower,cert.blocks(light).upper])));
        val=val+block;
    end
    S(k)=val;assert(val>=cert.Slo-5e-12 && val<=cert.Shi+5e-12);
end
maxWitnessIntegralError=0;testCaps=[];targetHalfCaps=[];
for k=1:5
    w=s.witnesses(k);selected=reshape(3*(w.tests(:)'-1)+(1:3)',1,[]);
    I=zeros(1,numel(selected));
    for i=1:numel(selected),I(i)=abswitness(P(:,:,selected(i)),w.nodes(:));end
    Iq=abswitness(P(:,:,15+w.targetCoordinate),w.nodes(:));
    maxWitnessIntegralError=max([maxWitnessIntegralError,max(abs(I-w.testIntegrals(:)')),abs(Iq-w.targetIntegral)]);
    testCaps=[testCaps,s.delta./(I+s.delta)]; %#ok<AGROW>
    targetHalfCaps=[targetHalfCaps,(w.nominalTarget-s.tau)/(2*(Iq+s.tau))]; %#ok<AGROW>
end
assert(maxWitnessIntegralError<5e-13);
K=min([.5,testCaps,targetHalfCaps]);C=s.tau+s.delta*[s.certificates.N]+S;
axis=min([K,[s.certificates.margin]./C]);assert(abs(axis-s.kappaAxis)<5e-12);
maxUpper=0;
for p=1:size(s.vertices,1)
    point=s.vertices(p,:);k=point(1);rho=point(2);
    upper=([s.certificates.base]+k*(s.delta*[s.certificates.N]+S)+2*rho*[s.certificates.U])/(1-k);
    assert(all(upper<=s.tau+5e-12));maxUpper=max(maxUpper,max(upper));
    for w=s.witnesses'
        assert(max(k*w.testIntegrals/(1-k))<=s.delta+5e-14);
        assert((w.nominalTarget-k*w.targetIntegral)/(1+k)>s.tau);
    end
end
k=.1;actual=(1-(1-k))/(1+(1-k));assert(abs(actual-k/(2-k))<1e-15 && actual>k/2);
covered=0;
for mask=0:30
    T=find(bitget(mask,1:5));found=false;
    for i=1:5,found=found||all(ismember(T,s.witnesses(i).tests));end
    assert(found);covered=covered+1;
end
out=struct('status','PASS','optimizerCalls',0,'rawCIETablesRead',4,...
 'maximumBlockIntegralDiscrepancy',maxBlockError,'maximumWitnessAbsoluteIntegralDiscrepancy',maxWitnessIntegralError,...
 'kappaAxisDiscrepancy',abs(axis-s.kappaAxis),'witnessCapDiscrepancy',abs(K-s.witnessKappaCap),...
 'verticesChecked',size(s.vertices,1),'maximumVertexUpper',maxUpper,'excludedProtocols',covered,'denominatorControlPassed',true);
fid=fopen(fullfile(root,'verification','matlab-171.json'),'w');assert(fid>=0);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(out,PrettyPrint=true));disp(out)
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
        f=tables{1}(:,channel+1);row=3*(light-1)+channel;
        P(:,:,row)=[E(1:400).*f(1:400),diff(E).*f(1:400)+E(1:400).*diff(f),diff(E).*diff(f)];
    end
    Z=sum(P(:,1,3*light-1)+P(:,2,3*light-1)/2+P(:,3,3*light-1)/3);
    assert(Z>0);P(:,:,3*light-2:3*light)=P(:,:,3*light-2:3*light)/Z;
end
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
    lo=cuts(k);hi=cuts(k+1);val=val+abs(a*(hi-lo)+b*(hi^2-lo^2)/2+c*(hi^3-lo^3)/3);
end
end

function val=abswitness(coeff,nodes)
val=0;
for interval=1:400
    a=coeff(interval,1);b=coeff(interval,2);c=coeff(interval,3);
    d=nodes(interval);z=nodes(interval+1)-d;cuts=[0,1];
    if z~=0 && -d/z>0 && -d/z<1,cuts=[0,-d/z,1];end
    powers=[a*d,b*d+a*z,c*d+b*z,c*z];
    for j=1:numel(cuts)-1
        val=val+abs(sum(powers.*(cuts(j+1).^(1:4)-cuts(j).^(1:4))./(1:4)));
    end
end
end
