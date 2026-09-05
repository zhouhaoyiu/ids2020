function [g,locasub,dep]=ids_getg(pathg,fault,grid,gridsize,source,loca,epi,srate,indx,index,deptha)
%==========================================================================
% g=ids_getg(pathg,fault,grid,gridsize,source,loca,epi,srate)
% This is a function to read the Green's functions from the database for
% only use of fixing mechanism (e.g., IDS method)
%--------------------------------------------------------------------------
% Input
%    pathg: path of data holder of green's functions
%    fault: [strike, dip, rake]
%     grid: numbers of sub-faults along [dip, strike] directions
% gridsize: sub-fault size in [dip, strike] directions, unit is in km
%   source: the index of the sub-fault hypocenter located on, [dip, strike]
%     loca: location of stations,[latitude,longitude]
%      epi: epicentral location, [latitude,longitude]
%    srate: sampling rate of green's functions wanted
%
% Output
%       g: output green's functions, a 3-D matrix, the size is [leng,3,nsta*nsub]
% locasub: lcoations of sub-faults,[latitude, longitude]
%     dep: depth of sub-fault, in kilometers
%
% 补充输入
%   indx,deptha,index 用于整体平移子断层深度：让第 index 个计算深度与 deptha(indx) 对齐。
%
% 输出尺寸和排列
%   nsub=prod(grid)，S=size(loca,1)。g 的第三维共有 S*nsub 个“子断层—台站”组合，
%   排列顺序是：第 1 个子断层的全部 S 个台站，随后第 2 个子断层的全部 S 个台站。
%   g(:,1:3,k) 最终依次是 E、N、U 三分量；locasub 为 nsub*2，dep 为 nsub*1。
%
% 数据流
%   走向/倾角/滑动角 -> 六分量矩张量 -> 子断层中心位置和深度
%   -> 枚举每个子断层与每个台站的距离、方位 -> 从数据库读取 V/R/T 格林函数
%   -> 按反方位角旋转为 E/N/U。
%--------------------------------------------------------------------------
%             Zhang Yong, Peking University, 2014-05-06
%==========================================================================

% get the moment tensor elements from the fault geometry
% 单位标量矩只保留机制形状，得到与数据库十个基本格林函数组合所需的六分量系数。
M=fp2mt_zh(1,fault(1),fault(2),fault(3));% 由（strike,dip,rake）计算地震矩张量6个分量

% transfer dip priority to strike priority
% 外部 grid/source 使用 [dip,strike]；get_subloca 使用 [strike,dip]，并在 source 前加段号 1。
grid=grid([2,1]);gridsize=gridsize([2,1]);source=[1,source([2,1])];

% get the location and depth of sub-faults
tic;[locasub,dep]=get_subloca(fault,grid,gridsize,source,loca,epi);toc
% 给所有子断层深度加同一个常数，使选定子断层匹配外部深度参考。
dep=dep+(deptha(indx)-dep(index)); % 新添加

% get the pairs of each sub-fault and each station, and calculate the
% distances and azimuthes
tic
% locasta 把 S 个台站坐标整体重复 nsub 次。
locasta=repmat(loca(:,1:2),[size(locasub,1),1]); 
% repmat_zh 让每个子断层坐标和深度连续重复 S 次，与 locasta 逐行配对。
locasuba=repmat_zh(locasub,size(loca,1));      
depa=repmat_zh(dep,size(loca,1));              

% dasub0 每行给出一个“台站—子断层”对的水平距离和方位角；前置深度得到 [dep,dist,azim]。
dasub0=da_zh(locasta,locasuba,1); 
dasub=[depa,dasub0];

toc

% read the Green's functions
% seekg_wang 按深度、距离查库并组合机制，返回 V/R/T 分量。
tic;[g,~]=seekg_wang(pathg,dasub,M,srate);toc;

% get the back azimuth of the subfaults relative to stations
bda=da_zh(locasuba,locasta,1);
% transfer the back azimuth to ray angle
% 加 180 度得到代码旋转公式使用的射线路径方向；sin/cos 对超过 360 度仍周期等价。
bfai=bda(:,2)+180;

%-------------------------------------------------------------------------
% convert green's functions from (V,R,T) to (E,N,U)
tic
sfai=sin(bfai*pi/180);cfai=cos(bfai*pi/180);
for i=1:size(g,3)    
%      m=[0,0,1;sfai(i),cfai(i),0;-cfai(i),sfai(i),0];
%      g(:,:,i)=g(:,:,i)*m;    
     
    v=g(:,1,i);r=g(:,2,i);t=g(:,3,i);
    % 径向 R、切向 T 按方位角投影为东向 E 和北向 N；垂向 V 直接作为 U。
    g(:,1,i)=r*sfai(i)-t*cfai(i); % E
    g(:,2,i)=r*cfai(i)+t*sfai(i); % N
    g(:,3,i)=v; % U
end
toc

