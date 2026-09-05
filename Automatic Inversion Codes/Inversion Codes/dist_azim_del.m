function [obc,locac,gc,mmc,NS,EW,UD]=dist_azim_del(ob,loca,g,mm,epi,maxdis,minsta,invsta,azista)
% This is a function to screen the stations based on the epicentral distance,
% the number of staions, the inter-station spacing, and inter-azimuth spacing

% Input: ob: the obverved data (1/3=EW 2/3=NS 3/3=UD ) 
%      loca: the location of stations 
%        mm: station name
%      maxdis: the maximum epicetral distance
%      minsta: the minimum number of stations  
%      invsta: inter-station spacing
%      azista: inter-azimuth spacing
%           g: green's functions

% Output: ob/loca/mm/EW/NS/UD
% -------------------------------------------------------------------------
% 中文说明
%   ob 为 N*(3S)，g 为 N*(3S)*nsub，列顺序是 [全部EW,全部NS,全部UD]；
%   loca 为 S*2，mm 为 S 行台站名，epi=[lat,lon]。
%   maxdis 是最大震中距(km)，minsta 是距离筛选后希望保留的最少台站数，
%   invsta 是台站间最小距离(km)，azista 是方位角分箱宽度(度)。
%
% 输出
%   obc、gc、locac、mmc 是三轮筛选后仍一一对应的数据。
%   EW、NS、UD 是 obc 的三个分量块，尺寸均为 N*最终台站数。
%
% 筛选顺序
%   1. 根据重复坐标删除重复台站；
%   2. 在不会使台站数降到 minsta 以下时，删除 maxdis 以外台站；
%   3. 根据两两台站距离删除过密台站；
%   4. 在每个方位角扇区内按震中距排序，再删掉距离过近的台站。
% 每轮都用三分量扩展下标同步删除 ob 和 g，避免台站元数据与波形错位。
% ob=ob1;loca=loca2;mm=mm2;maxdis=400;minsta=15;invsta=10;azista=5;
% Reove the repeated stations
[~,n]=size(ob);
% unique 返回各唯一坐标在原 loca 中的代表下标；代表下标之间的缺口被视为重复行位置。
[~,b]=unique(loca,'rows');
[ba,~]=sortrows(b);
b2=diff(ba);
c=b2>1;
nca=ba(c)+1;
% 将台站下标扩展到三个分量块，并同步删除所有关联数组。
ndca=[nca,n/3+nca,2*n/3+nca];
ob(:,ndca)=[];loca(nca,:)=[];mm(nca,:)=[];g(:,ndca,:)=[];
[m,n]=size(ob);
% 把剩余 ob 拆成 N*S 的 EW、NS、UD，后续每轮同步删除列。
oba=reshape(ob,m*3,n/3);                                 
EW=oba(1:m,:);NS=oba(m+1:2*m,:);UD=oba(2*m+1:end,:); 

% Screening the stations based on epicentral distance
% da(:,1) 为震中距；只有当前台站多于 minsta 时才尝试距离裁剪。
da=da_zh(loca,epi,1);                                           
if size(da,1)>minsta                                                          
    ndex=find(da(:,1)>maxdis);                                 
    % 仅当确有远台站且删除后仍严格多于 minsta 时执行；否则整批保留。
    if size(ndex,1)~=0&&size(da,1)-size(ndex,1)>minsta 
            mdex=[ndex,n/3+ndex,2*n/3+ndex];
            EW(:,ndex)=[];NS(:,ndex)=[];UD(:,ndex)=[];
            oba=ob;locaa=loca;mma=mm;ga=g;
            oba(:,mdex)=[];ga(:,mdex,:)=[];
            locaa(ndex,:)=[];
            mma(ndex,:)=[];
    else
        oba=ob;locaa=loca;mma=mm;ga=g;
    end
else
    oba=ob;locaa=loca;mma=mm;ga=g;
end


% screen the stations based on the inter-azimuth spacing
% Caculating the inter-station spacing
% repmat 与 kron 生成全部 S*S 台站对；distance 逐行计算球面距离。
locasta=repmat(locaa,[size(locaa,1),1]);                              
locastb=kron(locaa,ones(size(locaa,1),1));                            % Kronecker Product Operator
da_sta=distance(locasta,locastb)*6371*pi/180;                         
da_sta(all(da_sta==0,2),:)=[];
% 删除 S 个自身距离后，恢复为 (S-1)*S：每列对应一个台站与其余台站的距离。
da_sta=reshape(da_sta,size(locaa,1)-1,size(locaa,1));
% zdex 是所有小于 invsta 的近邻关系在线性展开矩阵中的位置。
zdex=find(da_sta<invsta);                                            

if size(zdex,1)>0 
    zdx1=ceil(zdex/size(da_sta,1));                                 
    zdx2=mod(zdex,size(da_sta,1));                                   
    % 把线性位置还原成候选台站对，随后修正 mod 等于 0 的列边界情况。
    zdexa=[zdx2,zdx1];                                               
    r=zdexa(:,2)-zdexa(:,1);
    zdexb=zeros(size(zdexa));                                      
    for i=1:size(zdexa,1)                                           
        if r(i)<=0
            zdexb(i,:)=[zdexa(i,1)+1,zdexa(i,2)];
        else
            zdexb(i,:)=zdexa(i,:);
        end
    end
    % If the inter-station spacing between Station A and two or more stations is less than the limit value, station A will be deleted
    % 比较配对表第二列的连续重复值，找出同时拥有多个过近邻站的台站。
    va=diff(zdexb(:,2));
    m=1;vb=zeros(1,1);
    for i=1:size(va,1)
        if i<size(va,1)&&va(i)==0&&(va(i+1)-va(i))~=0
            vb(m)=i; 
            m=m+1;
        elseif i==size(va,1)&&va(i)==0&&(va(i)-va(i-1))~=0
            vb(m)=i; 
            m=m+1;
        else
            continue
        end
    end
    if vb(1)~=0
        % idex 是本轮待删台站号，sdex 是其三分量列号。
        idex=zdexb(vb,2);idex=idex';
        [~,n]=size(oba);
        % The second update of the data
        sdex=[idex,n/3+idex,2*n/3+idex];
        EW(:,idex)=[];NS(:,idex)=[];UD(:,idex)=[];
        obb=oba;locab=locaa;mmb=mma;gb=ga;
        obb(:,sdex)=[];gb(:,sdex,:)=[];
        locab(idex,:)=[];
        mmb(idex,:)=[];
    else
        obb=oba;locab=locaa;mmb=mma;gb=ga;
    end
else
    obb=oba;locab=locaa;mmb=mma;gb=ga;
end

% screen the stations based on the inter-azimuth spacing
% 对剩余台站重新计算震中距 dax 和方位角 dat。
das=da_zh(locab,epi,1);dax=das(:,1);dat=das(:,2);
intv=azista;                                                         
inva=0:intv:360;invb=intv:intv:360+intv;
invc=[inva;invb];invc=invc';invc=invc(1:end-1,:);
% invc 每行是一个 (low,high] 方位角区间；del 累计待删台站，invd 记录空扇区。
del=[];invd=[];
for i=1:size(invc,1)
    jm=find(invc(i,1)<dat&dat<=invc(i,2));
    if size(jm,1)~=0
    % 当前扇区内按震中距由近到远排序，再检查相邻距离差。
    dis=dax(jm);
    [disa,ndex]=sortrows(dis);  
    disb=diff(disa);
    im=find(disb<invsta);    
    if size(im,1)~=0
        % 对连续过近关系每隔一个取后者，减少同一密集序列被重复删除的机会。
        im=im(1:2:end)+1;
        de=jm(ndex(im));
        del=[del;de];
    else
        continue
    end
    elseif size(jm,1)==0
        invd=[invd;invc(i,:)];   
    else
        continue
    end
end
[jdex,~]=sortrows(del);jdex=jdex';
[~,n]=size(obb);

% The second update of the data
if size(jdex,1)~=0
    % 最后一轮也把台站号展开到三个分量列，并同步裁剪输出数据。
    ydex=[jdex,n/3+jdex,2*n/3+jdex];
    EW(:,jdex)=[];NS(:,jdex)=[];UD(:,jdex)=[];
    obc=obb;locac=locab;mmc=mmb;gc=gb;
    obc(:,ydex)=[];gc(:,ydex,:)=[];
    locac(jdex,:)=[];
    mmc(jdex,:)=[];
else
    obc=obb;locac=locab;mmc=mmb;gc=gb;
end
end

