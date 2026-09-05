function [g_end,dep_out]=seekg_wang(pathdir,dda,M,rate)
%==========================================================================
% [g_end,dep_out]=seekg_wang(pathdir,dda,M,rate)
% This is a function to read the Green's functions from the database 
%--------------------------------------------------------------------------
% Input
%  pathdir: the filefolder of database
%      dda: [dep,dist,azim], azim is defined as north-east
%     depg: depth of green's function in the database
%        M: [Mxx,Mxy,Mxz,Myy,Myz,Mzz]
%     rate: the sampling rate of the output green's functions you want
%
% Output
%   g_end: output green's functions, the size of [leng,3,n_of_channels]
%  dep_out: 每个输出通道实际采用的数据库震源深度，单位 km。
%
% 中文说明
%   dda 为 N*3，每行 [目标深度(km),水平距离(km),方位角(度)]；
%   M 为六分量矩张量 [Mxx,Mxy,Mxz,Myy,Myz,Mzz]；rate 为目标采样率。
%   输出 g_end 为 L*3*N，第二维依次是数据库坐标约定下的 V、R、T 分量。
%
% 数据流
%   读取 GreenInfo.dat 的长度、采样间隔、距离网格和深度网格
%   -> 计算每个目标的 P/S 理论到时 -> 按深度分组并选最近数据库深度
%   -> getg_wang_sub0 按最近距离读取十个基本格林函数并校正到时
%   -> 用矩张量和方位角把十个基本分量线性组合成 V/R/T 三分量。
%
% 路径说明
%   当前源码用反斜杠拼接文件名，按 Windows 路径编写；在 macOS/Linux 上反斜杠
%   不是目录分隔符。这里保留原计算语句，跨平台运行时应另行确认 pathdir 的写法。
%--------------------------------------------------------------------------
%             Zhang Yong, Peking University, 2014-05-06
%==========================================================================

% 查找并打开数据库说明文件。
headfile=dir([pathdir,'\GreenInfo.dat']);
fid1=fopen([pathdir,headfile(1).name],'r');


% 跳过前 6 行文本后读取数值参数；tpara(2) 是采样间隔，tpara(3) 是记录长度。
for i=1:6;fgets(fid1);end
tpara=fscanf(fid1,'%f');

% 再跳过两行说明，读取数据库距离网格；首个数是数量或头值，所以删除。
for i=1:2;fgets(fid1);end
dist0=fscanf(fid1,'%f');
dist0(1)=[]; 

z=fgets(fid1);
while 1
    % 向后寻找深度列表标题。当前循环没有 EOF 保护，说明文件必须含该固定标题。
    if length(z)>25
        if strcmp(z(1:25),'#   list of source_depths')
        break;
        end
    end
    z=fgets(fid1);
end
%for i=1:20;fgets(fid1);end
numdep=str2num(fgets(fid1));
% 每个深度占说明文件中 11 个字符，逐行读出后统一转为数值向量 depg。
m=repmat(' ',[numdep,11]);
for i=1:numdep
    z=fgets(fid1);
    m(i,:)=z(1:11);
end
depg=str2num(m);
fclose(fid1);

leng=tpara(3);

% 数据库采样率等于采样间隔的倒数。
srate=(1/tpara(2));

if nargin<4
    % 未指定目标采样率时保持数据库原采样率。
    rate=srate;
end

%rate

dep=dda(:,1);
dist=dda(:,2);
azi=dda(:,3);

%make a transfer for azim from 'north to east' to 'south to east'
% 数据库机制组合公式使用“南起向东”的角度，故从输入“北起向东”转换。
azi=180-azi;

%--------------------------------------------------------------------------
% get the timeshift of green's function
km2deg=6371*pi/180;

% get_time 接收度，数据库和输入距离均从 km 换算为角距离。
distdeg=dist/km2deg;dist0deg=dist0/km2deg;
if max(dep)<100
    % 浅于 100 km 使用本包的 IASPEI91 走时表。
    tp=get_time(distdeg,dep,'P',1); % has the same size of dist or dep
    ts=get_time(distdeg,dep,'S',1); % has the same size of dist or dep
    tp0=get_time(dist0deg,dep,'P'); % has the size [length(dist0),length(dep)]
    ts0=get_time(dist0deg,dep,'S'); % has the size [length(dist0),length(dep)]
else
    % 深度达到 100 km 时改调 dbgrn_get_time；该函数不在本目录 39 个源码文件中。
    tp=dbgrn_get_time(distdeg,dep,'P',1); % has the same size of dist or dep
    ts=dbgrn_get_time(distdeg,dep,'S',1); % has the same size of dist or dep
    tp0=dbgrn_get_time(dist0deg,dep,'P'); % has the size [length(dist0),length(dep)]
    ts0=dbgrn_get_time(dist0deg,dep,'S'); % has the size [length(dist0),length(dep)]
end
%--------------------------------------------------------------------------
%max(dist(:))
%------------------------------------------------
[sdep,nd]=sort(dep);
% 距离、目标到时及数据库到时列都按同一深度排序，nd 保存恢复原顺序的映射。
sdist=dist(nd);
tp=tp(nd);
ts=ts(nd);
tp0=tp0(:,nd);
ts0=ts0(:,nd);

ddep=diff(sdep);
% dex 把完全相同的相邻深度划为一组，每组只打开一个深度文件。
dex=[0;find(ddep>0);length(sdep)];

% for different depth
%lengg=150*srate+leng; % srate is the sampling rate of GRN FUNCs in database 
lengg=115*rate+round(leng*rate/srate); % use rate here, (rate is the sampling rate to be used)
% 额外预留 115 s，以容纳走时校正造成的延长。

% get the vrt components of the corresponding mechanism
% six typical source mechanisms
% 将六分量矩张量分解成数据库十个基本格林函数使用的爆炸、CLVD、走滑和倾滑系数。
exp=sum(M([1,4,6]))/3;
clvd=M(6)-exp;
ss12=-M(2); %[90,90,0]; Mtp=-Mxy
ss11=(M(1)-M(4))/2; %[225,90,-180]; (Mtt-Mpp)/2=(Mxx-Myy)/2;
ds31=M(3); %[270,90,-90]; Mrt=Mxz
ds23=-M(5); %[180,90,-90]; Mpr=-Myz

sina=sin(azi*pi/180);cosa=cos(azi*pi/180);
% 一倍角用于倾滑项，二倍角用于走滑项。
sin2a=sin(2*azi*pi/180);cos2a=cos(2*azi*pi/180);

% new codes:
%--------------------------------------------------------------------------
g_end=zeros(lengg,3,length(dist));
dep_out=zeros(length(dist),1);

for i=1:length(dex)-1
    % 当前组的所有请求深度应完全相同；否则分组边界或排序出现异常。
    dep_temp=sdep(dex(i)+1);
    if sum(abs(sdep(dex(i)+1:dex(i+1))-dep_temp))>0
        error('something wrong!');
    end
    
    [~,ndep]=min(abs(depg-dep_temp));
    % 不在深度方向插值，直接选择与 dep_temp 最近的数据库深度文件。
    
    fileg=[pathdir,'\grn_d',num2str(round(depg(ndep)))];

    dex_temp=dex(i)+1:dex(i+1);
    
    dist_temp=sdist(dex_temp);
    tp_temp=tp(dex_temp);
    ts_temp=ts(dex_temp);
    t_temp=[tp_temp,ts_temp];
    t0_temp=[tp0(:,dex(i)+1),ts0(:,dex(i)+1)];% any of [dex(i)+1:dex(i+1)] is OK
    
    m=m+1;
    g_temp=getg_wang_sub0(fileg,leng,dist_temp,t_temp,t0_temp,srate,dist0,rate);
    % gdex 把当前深度组中的排序位置换回原 dda 行号。
    gdex=nd(dex(i)+1:dex(i+1));
    lgtemp=size(g_temp,1);

    for j=1:length(gdex)
        % m1、m2 是当前方位角下十个基本分量的线性组合权重。
        m1=[exp;clvd;(ss12*sin2a(gdex(j))+ss11*cos2a(gdex(j)));(ds31*cosa(gdex(j))+ds23*sina(gdex(j)))];
        m2=[(ss12*cos2a(gdex(j))-ss11*sin2a(gdex(j)));(ds31*sina(gdex(j))-ds23*cosa(gdex(j)))];
        g_end(1:lgtemp,1,gdex(j))=g_temp(:,[1,9,3,6],j)*m1;
        g_end(1:lgtemp,2,gdex(j))=g_temp(:,[2,10,4,7],j)*m1;
        g_end(1:lgtemp,3,gdex(j))=g_temp(:,[5,8],j)*m2;
        
        % 记录实际选中的离散数据库深度，便于检查深度近似误差。
        dep_out(gdex(j))=depg(ndep);
    end
end

return
%=======================================================================end


% old codes:
%--------------------------------------------------------------------------
% 下面是 return 之后的旧实现，正常调用不会执行。它先保存全部十个分量，再二次循环组合；
% 上面的新实现读取一组就立即组合，减少中间数组 g 的内存占用。
g=zeros(lengg,10,length(dist));
for i=1:length(dex)-1
    % dep( dex(i)+1:dex(i+1) ) has the same focal depth
    
    dep_temp=sdep(dex(i)+1);
    if sum(abs(sdep(dex(i)+1:dex(i+1))-dep_temp))>0
        error('something wrong!');
    end
    
    [~,ndep]=min(abs(depg-dep_temp));
    fileg=[pathdir,'\grn_d',num2str(round(depg(ndep)))];
    

    dex_temp=dex(i)+1:dex(i+1);
    sdit_temp=sdist(dex_temp);
    dist_temp=sdit_temp;
    
    tp_temp=tp(dex_temp);
    ts_temp=ts(dex_temp);
    t_temp=[tp_temp,ts_temp];
    t0_temp=[tp0(:,dex(i)+1),ts0(:,dex(i)+1)];% any of [dex(i)+1:dex(i+1)] is OK
    
    g_temp=getg_wang_sub(fileg,leng,dist_temp,t_temp,t0_temp,srate,dist0,rate);
    g(1:size(g_temp,1),:,nd(dex(i)+1:dex(i+1)))=g_temp; % equals 1+2
end

%------------------------------------------------

g_end=zeros(size(g,1),3,size(g,3));
for i=1:size(g,3)
    m1=[exp;clvd;(ss12*sin2a(i)+ss11*cos2a(i));(ds31*cosa(i)+ds23*sina(i))];
    m2=[(ss12*cos2a(i)-ss11*sin2a(i));(ds31*sina(i)-ds23*cosa(i))];
    g_end(:,1,i)=g(:,[1,9,3,6],i)*m1;
    g_end(:,2,i)=g(:,[2,10,4,7],i)*m1;
    g_end(:,3,i)=g(:,[5,8],i)*m2;
end
return

