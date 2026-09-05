function [obn,gn,locan,mmn,stfb]=divi_inver(epi,loca,ob,g,mm,srate,grid,gridsize,source,miu,flh,iter,Mw)
%% This function is to the station grouping, and screen the availble stations based on waveforms fittings（misfit<=0.6)
% Input:   ob: the observed data
%          g1: green's functions
%        loca: locations of stations
%          mm: station names
%        iter: maximum iteration number
% Output:obn/gn/locan/mmn
% -------------------------------------------------------------------------
% 中文说明
%   epi : [lat,lon] 震中；loca 为 S*2 台站坐标；mm 为 S 行台站名。
%   ob  : N*(3S) 观测波形，列顺序是 [全部EW,全部NS,全部UD]。
%   g   : N*(3S)*nsub，与 ob 使用相同通道顺序。
%   srate、grid、gridsize、source、miu、flh、iter、Mw 直接传给 ids_data。
%
% 输出
%   obn  : 筛选并作整数时移后的观测波形，仍按 [EW,NS,UD] 三个列块排列。
%   gn   : 与 obn、locan 对齐的格林函数。
%   locan: 保留下来的台站坐标；mmn 为对应台站名。
%   stfb : Nt*G，每列是一个台站组反演得到、时间积分归一为 1 的总源时间函数。
%
% 数据流
%   selec_sta 按距离和方位构造台站组 -> 每组去重 -> 独立 IDS 反演
%   -> cal_res 删除至少两个分量拟合不佳的台站 -> 按互相关时移观测波形
%   -> 汇总各组的 EW/NS/UD 波形和格林函数。
%
% 分组要求
%   selec_sta 当前把每组设为 24 个台站（12 个基准台站和各自配对台站）。
%   若筛选后不足一组，grp 可能为 0，末组下标会失效；调用前需保证台站数满足该条件。
% ----- station grouping  -----
% idex 是 selec_sta 返回的每组台站数；loca/ob/g/mm 已按组重新排列。
[loca,ob,g,mm,idex]=selec_sta(epi,loca,ob,g,mm);
n=size(loca,1);                                              
% 这些空数组用于逐组追加最终保留的数据。
obn=[];gn=[];locan=[];mmn=[];                                
% grp 是完整组数；循环处理前 grp-1 组，余数统一交给最后一组。
grp=floor(size(loca,1)/idex);                                             % confirm the group number,idex=the number of stations per group    
ew=[];ns=[];ud=[];gew=[];gns=[];gud=[];
stfb=[];
% ----- inversions of individual groups of stations by the IDS method -----
for i=1:grp-1
    % dex 是本组台站行号；ndex 扩展为对应 EW、NS、UD 三个通道块的列号。
    dex=idex*(i-1)+1:idex*i;
    ndex=[dex,n+dex,2*n+dex];                                             % Coolumn of EW/NS/UD components
    ob0=ob(:,ndex);g0=g(:,ndex,:);
    loca0=loca(dex,:);mm0=mm(dex,:);
        
    % ----- delete duplicates（loca0/mm0/ob0/g0） -----
    % 同坐标台站只保留 unique 返回的一条，并用相同下标同步裁剪名称、波形和格林函数。
    m=size(loca0,1);   
    [loca0,adex]=unique(loca0,'rows');                                   
    mm0=mm0(adex,:);
    vdex=[adex,m+adex,2*m+adex];
    % 转为列向量后，输出顺序仍是三个连续分量块。
    vdex=vdex(:);
    ob0=ob0(:,vdex);                                                      % Rearrange the observed data
    g0=g0(:,vdex,:);                                                      % Rearrange the green's fucntions
        
    % ----- Inversion -----
    % 此处 nsta 是波形通道数，即实际台站数的 3 倍。
    nsta=size(ob0,2);                                                  
    % ids_data 返回本组矩率、合成波形 syn 以及经过滤波并恢复振幅的观测 obr。
    [~,substf,~,syn,obr,~,~,~,~]=ids_data(ob0,g0,loca0,srate,grid,gridsize,source,nsta,miu,flh,iter,Mw); 
    % 所有子断层矩率沿第二维相加成总 STF；积分末值用于把曲线归一到单位矩。
    stf=sum(substf,2);
    stfa=cumsum(stf)/srate;
    stf=stf/stfa(end);
    stfb=[stfb,stf];
    
    % ----- Select the stations based on waveform fittings -----
    % nc 是待删台站号，nzc 是其三分量列号，dt 是每个通道的整数时移。
    [nc,nzc,dig,dt]=cal_res(obr,syn,srate,floor(1/flh(1)/4));
    % 所有数组按同一台站/通道下标同步删除，保持后续一一对应。
    dig(:,nzc)=[];
    ob0(:,nzc)=[];g0(:,nzc,:)=[];dt(:,nzc)=[];loca0(nc,:)=[];mm0(nc,:)=[];syn(:,nzc)=[];obr(:,nzc)=[];
%     ob0=ob0./dig;                                                         

    for j=1:length(dt)
        if dt(j)>0
            % 正时移删除开头、末尾补零；负时移开头补零、删除末尾。
            ob0(:,j)=[ob0(dt(j)+1:end,j);zeros(dt(j),1)];
        elseif dt(j)<0
            ob0(:,j)=[zeros(-dt(j),1);ob0(1:end+dt(j),j)];
        end
    end
        
    locan=[locan;loca0];mmn=[mmn;mm0];
    m=size(loca0,1);
    % 把本组观测重新拆成三分量，分别追加；最后才能恢复全局 [EW,NS,UD] 块顺序。
    ewa=ob0(:,1:m);nsa=ob0(:,m+1:2*m);uda=ob0(:,2*m+1:3*m);
    ew=[ew,ewa]; % the east-west components of observed data                  
    ns=[ns,nsa]; % the north-south components of observed data                                                         
    ud=[ud,uda]; % the vertical components of observed data                                                           
    
    % Rearrange the green's fucntions 
    % 格林函数按与观测相同的分量和台站顺序追加。
    gewa=g0(:,1:m,:);gnsa=g0(:,m+1:2*m,:);guda=g0(:,2*m+1:3*m,:);
    gew=[gew,gewa]; % the east-west components of green's functions
    gns=[gns,gnsa]; % the north-south components of green's functions                                                 
    gud=[gud,guda]; % the vertical components of green's functions                                   
end

% ----- the last group -----
% 最后一组从第 grp 个组首位一直取到末尾，因此同时容纳不足一整组的余数。
dex=idex*(grp-1)+1:size(loca,1);
ndex=[dex,n+dex,2*n+dex];
ob0=ob(:,ndex);g0=g(:,ndex,:);
loca0=loca(dex,:);mm0=mm(dex,:);
    
% ----- delete duplicates（loca0/mm0/ob0/g0） -----
% 与前面各组相同：按坐标去重，并同步重排三分量数据。
m=size(loca0,1);
[loca0,adex]=unique(loca0,'rows');                                      
mm0=mm0(adex,:);
vdex=[adex,m+adex,2*m+adex];
vdex=vdex(:);
ob0=ob0(:,vdex);                                                          % Rearrange the observed data
g0=g0(:,vdex,:);                                                          % Rearrange the green's fucntions
    
% ----- Inversion -----
% 对最后一组执行相同 IDS 反演，并把归一化总 STF 追加到 stfb。
nsta=size(ob0,2);
[~,substf,~,syn,obr,~,~,~,~]=ids_data(ob0,g0,loca0,srate,grid,gridsize,source,nsta,miu,flh,iter,Mw);
stf=sum(substf,2);
stfa=cumsum(stf)/srate;
stf=stf/stfa(end);
stfb=[stfb,stf];  

% -----  Select the stations based on waveform fittings -----
% 删除不合格台站及其三分量，然后按 dt 移动保留观测。
[nc,nzc,dig,dt]=cal_res(obr,syn,srate,floor(1/flh(1)/4));
dig(:,nzc)=[];
ob0(:,nzc)=[];g0(:,nzc,:)=[];dt(:,nzc)=[];loca0(nc,:)=[];mm0(nc,:)=[];
% ob0=ob0./dig;
    for j=1:length(dt)
        if dt(j)>0
            ob0(:,j)=[ob0(dt(j)+1:end,j);zeros(dt(j),1)];
        elseif dt(j)<0
            ob0(:,j)=[zeros(-dt(j),1);ob0(1:end+dt(j),j)];
        end
    end    
locan=[locan;loca0];mmn=[mmn;mm0];
m=size(loca0,1);
% 将末组也拆入三个累积块，随后一次性形成最终输出。
ewa=ob0(:,1:m);nsa=ob0(:,m+1:2*m);uda=ob0(:,2*m+1:3*m);
ew=[ew,ewa];                                                            
ns=[ns,nsa];                                                          
ud=[ud,uda];                                                          
    
% Rearrange the green's fucntions 
gewa=g0(:,1:m,:);gnsa=g0(:,m+1:2*m,:);guda=g0(:,2*m+1:3*m,:);
gew=[gew,gewa];                                                      
gns=[gns,gnsa];                                                    
gud=[gud,guda];                                                  
% 最终恢复程序统一的 [全部EW,全部NS,全部UD] 通道排列。
obn=[ew,ns,ud];gn=[gew,gns,gud];
end

