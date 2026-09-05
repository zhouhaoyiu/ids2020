function vec_out=t_shift(vec_in,sdot)
%==========================================================================
%  vec_out=t_shift(vec_in,sdot)
%  here is a function to shift the wave (vec_in) for sdot
%--------------------------------------------------------------------------
% Input
%   vec_in: the origin wave, should be a vetor or 2-dimension matrix
%   sdot: a integer scalar describing the time shift by sampling interval
%         sdot>0, pull vec_in
%         sdot<0, press vec_in
% Output
%  vec_out: shifted wave, also a vector, having the length
%           of length(vec_in)+sdot
%
% 中文说明
%   vec_in 为 N*C 波形，第一维是时间；sdot 是整数采样点数。
%   输出长度为 N+sdot：sdot>0 时用插值增加点数，sdot<0 时减少点数。
%   这并非简单在首尾补零的平移，而是把同一首末端范围重新采样到新的点数。
%
% 数据怎样变化
%   1. 根据目标点数计算旧坐标上的步长 dt=(N-1)/(N-1+sdot)。
%   2. dex 给出每个新采样点落在旧波形中的小数下标。
%   3. floor/ceil 找左右相邻旧样点，wei1/wei2 按距离作线性插值。
%
% 例：vec_in=[0;10;20]、sdot=2，输出需要 5 点；dex=[1,1.5,2,2.5,3]，
%     所以 vec_out=[0;5;10;15;20]。首末值不变，中间采样更密。
% Notice: the function 't_shift' can also work as another function of
%        'resample', it is much faster than the function 'resample' of
%         Matlab, but with a slight effect of low pass filter.
%--------------------------------------------------------------------------
%        Zhang Yong, 2012-01-31 23:55, Berlin
%==========================================================================
if sdot==0
    % 不改变点数时直接返回原数据，避免不必要的插值。
    vec_out=vec_in;
    return
end

if size(vec_in,1)==2&&sdot==-1
    % 两点压缩成一点时，线性插值的唯一结果取两点平均值。
    vec_out=mean(vec_in);
    return;
end

% 新长度为 N+sdot；dt 表示相邻新点在旧数组下标坐标中的间隔。
dt=(size(vec_in,1)-1)/((size(vec_in,1)-1)+sdot);

% dex 从旧数组第 1 点走到第 N 点，长度通常为 N+sdot。
dex=(1:dt:size(vec_in,1))';

% onevec 把同一组时间插值权重扩展到 vec_in 的所有列。
onevec=ones(1,size(vec_in,2));

% 每个小数下标两侧的整数下标决定参与线性插值的两个旧样点。
fdex=floor(dex);
cdex=ceil(dex);
vec1=vec_in(int32(fdex),:);
vec2=vec_in(int32(cdex),:);

% wei2 是距左样点的小数部分，wei1=1-wei2；两者之和始终为 1。
wei2=(dex-fdex)*onevec;
wei1=1-wei2;

% 对两侧样点作加权平均，输出每列都具有相同的新时间网格。
vec_out=vec1.*wei1+vec2.*wei2;
return
%==========================================================================
    
