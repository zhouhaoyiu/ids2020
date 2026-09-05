function [loca0,ob0,g0,mma,idex,nst]=selec_sta(epi,loca,ob,g,mm)
% This function is to identify a partner station for each station 
% by requiring the two stations to have similar take-off angles and nearly opposite azimuths. 
% We use the two-dimensional Gaussian distribution function ( u1/u2 = take-off angle/azimuth)
% Station screening based on two-dimensional Gaussian Distribution (parameters:take-off angle and azimuth)

%  Input: 
%       epi: epicenter 
%      loca: locations of stations
%        ob: observed data 
%         g: green's functions 
%        mm： staion names
% Output: loca0/ob0/g0/mma
%       idex: the number of stations per group
%        nst: number of groups
% 
% -------------------------------------------------------------------------
% 中文说明
%   epi=[lat,lon]；loca 为 S*2；mm 为 S 行台站名。
%   ob 为 N*(3S)，g 为 N*(3S)*nsub，通道顺序必须是 [全部EW,全部NS,全部UD]。
%
% 输出
%   loca0、mma、ob0、g0 是重新分组后的台站、名称、观测和格林函数。
%   每个基准台站会配一个距离相近、方位近似相反的台站，因此输出中允许重复台站。
%   idex 最终固定为 24，即每个完整组含 12 个基准台站和 12 个配对台站。
%   nst=floor(筛选后台站数/12)；若有余数，代码另附到最后，但 nst 本身不再加 1。
%
% 数据流
%   按震中距排序 -> 删除 1000 km 外台站 -> 为每站寻找反方位配对台站
%   -> 按方位角排序 -> 从不同方位位置抽取 12 个基准台站组成一组
%   -> 将基准与配对台站的三分量数据按同一顺序重排。


%% Input：loca mm 
%----- calculate the epicentral distances, and rearrange the loca/ob/g/mm based on the epicentral distances -----
% da(:,1:2) 为震中距(km)和方位角(度)；sortrows 默认优先按第一列距离排序。
da=da_zh(loca,epi,1);                                              
[da1,dex]=sortrows(da);                                                   
n=size(loca,1);
% 把 S 个台站下标扩展到 EW、NS、UD 三个列块，并按相同距离顺序排列。
adex=[dex,n+dex,2*n+dex];
mm0=mm(1:n,:);                                                            
mm0=mm0(dex,:);                                                           % Station names are arranged according to the epicentral distance
aloca=loca(dex,:);                                                        % Station locations are arranged according to the epicentral distance
adex=adex(:);
oba=ob(:,adex);ga=g(:,adex,:);                                            % The observed data and green's functions are arranged according to the epicentral distance

% ----- remove the stations whose epicentral distances over 1000 km -----
% 台站已按距离升序，删除 aloca 的远台站后，dex 就是留下的近台站数。
aloca(da1(:,1)>1000,:)=[];                                                            
dex=size(aloca,1);
adex=[1:dex,n+1:n+dex,2*n+1:2*n+dex];
mm0=mm0(1:dex,:);
oba=oba(:,adex);ga=ga(:,adex,:); 

% ----- identify a partner station for each station based on take-off angle anda zimuth -----
% ---- we use the the two-dimensional Gaussian distribution function ( u1/u2 = take-off angle/azimuth)
% epa=[mean(aloca(:,1)),mean(aloca(:,2))];                               
da=da_zh(aloca,epi,1);                                                   
tazim=da;                                                              
segma1=std(tazim(:,1));                                                   % The standard deviation of the epicentral distance
segma2=std(tazim(:,2));                                                   % The standard deviation of the azimuth
% ndex(i) 将保存第 i 个基准台站的最佳配对台站下标。
ndex=zeros(size(tazim,1),1);
for i=1:size(tazim,1) % the number of stations
    [dex]=find_sta(tazim,i,segma1,segma2);
    ndex(i,:)=dex;
end

% ----- Pick out the relevant data of the partner stations -----
% 用同一个 ndex 同步提取配对台站的位置、名称、观测波形和格林函数。
bloca=aloca(ndex,:);                                                      % The stations locations
n=size(aloca,1);
adex=[ndex,n+ndex,2*n+ndex];
mm1=mm0(ndex,:);                                                          % The station name                                                 
adex=adex(:);
obb=oba(:,adex);gb=ga(:,adex,:);                                          % The obseved data and green's functions

% ----- The distance of the station relative to the center point is calculated and rearranged -----                                                 
[~,dex]=sortrows(da,2);
% 这里按 da 的第 2 列方位角排序；基准和配对数据都跟随同一 dex 重排。
n=size(aloca,1);
adex=[dex,n+dex,2*n+dex];
aloca=aloca(dex,:);bloca=bloca(dex,:);                                 
mm0=mm0(dex,:); mm1=mm1(dex,:);                                        
adex=adex(:); 
oba=oba(:,adex);ga=ga(:,adex,:);                                          
obb=obb(:,adex);gb=gb(:,adex,:); 
% 
% ----- Grouping based on azimuth -----                                   % Number of stations in each group
% 每组先取 12 个基准台站；配对后形成 24 个台站。12 是写死的经验设置。
idex=12;                                                              
nst=floor(size(aloca,1)/idex);                                            % The number of groups
aloc=zeros(1,2);bloc=zeros(1,2);
loca0=[];mma=[];ew=[];ns=[];ud=[];gew=[];gns=[];gud=[];
for j=1:nst 
    i=1:idex; 
    % Stations are grouped and rearranged
    % nst*(i-1)+j 相当于把方位角序列交错分组，使每组覆盖较宽的方位范围。
    aloc(i,:)=aloca(nst*(i-1)+j,:);
    bloc(i,:)=bloca(nst*(i-1)+j,:);
    cloc=[aloc;bloc];
    % 一个组中先放 12 个基准台站，再放各自的 12 个配对台站。
    loca0=[loca0;cloc];
    
    % Station names are grouped and rearranged
    amm(i,:)=mm0(nst*(i-1)+j,:);
    bmm(i,:)=mm1(nst*(i-1)+j,:);
    cmm=[amm;bmm];
    mma=[mma;cmm];
    
    % The obseved waveforms are grouped and rearranged
    % dex 是本组 12 个基准台站在方位排序数组中的下标。
    dex=nst*(i-1)+j;
    ewa=oba(:,dex);nsa=oba(:,n+dex);uda=oba(:,2*n+dex);
    ewb=obb(:,dex);nsb=obb(:,n+dex);udb=obb(:,2*n+dex);
    ewc=[ewa,ewb];nsc=[nsa,nsb];udc=[uda,udb];
    % 每个分量分别累计，函数末尾再拼回 [全部EW,全部NS,全部UD]。
    ew=[ew,ewc];                                                          
    ns=[ns,nsc];                                                         
    ud=[ud,udc];                                                          
    
    % The green's fucntions are grouped and rearranged
    gewa=ga(:,dex,:);gnsa=ga(:,n+dex,:);guda=ga(:,2*n+dex,:);
    gewb=gb(:,dex,:);gnsb=gb(:,n+dex,:);gudb=gb(:,2*n+dex,:);
    gewc=[gewa,gewb];gnsc=[gnsa,gnsb];gudc=[guda,gudb];
    gew=[gew,gewc];                                                       
    gns=[gns,gnsc];                                                       
    gud=[gud,gudc];                                                         
end

% ----- The last group/ Attach data that is not part of a group to the last group -----
if size(aloca,1)>(idex*nst)
    % 未进入完整 12 台组的余数，其基准台站和配对台站直接追加到末尾。
    loca0=[loca0;aloca(idex*nst+1:end,:);bloca(idex*nst+1:end,:)];       
    mma=[mma;mm0(idex*nst+1:end,:);mm1(idex*nst+1:end,:)];               
    dex=idex*nst+1:size(aloca,1);
    ew=[ew,oba(:,dex),obb(:,dex)];                                       
    ns=[ns,oba(:,n+dex),obb(:,n+dex)];
    ud=[ud,oba(:,2*n+dex),obb(:,2*n+dex)];
    
    gew=[gew,ga(:,dex,:),gb(:,dex,:)];                                    
    gns=[gns,ga(:,n+dex,:),gb(:,n+dex,:)];
    gud=[gud,ga(:,2*n+dex,:),gb(:,2*n+dex,:)];
end
% 对调用者而言一个完整反演组含基准和配对两部分，所以把 12 改为 24。
idex=idex*2;                                                              
% 将三个分量累积块拼成程序统一的通道顺序。
ob0=[ew,ns,ud];g0=[gew,gns,gud];
end

