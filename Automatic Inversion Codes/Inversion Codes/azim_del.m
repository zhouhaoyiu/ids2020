function [del_index]=azim_del(fai,interval)
% this is a function to reselect the stations based on inter-station/azimuth spacing
% AZIM_DEL 按方位角间隔筛掉过密的台站。
%
% 输入
%   fai      : 每个台站的方位角，单位应与 interval 一致，通常为度。
%   interval : 允许保留的最小方位角间隔。
%
% 输出
%   del_index: 应删除元素在原始 fai 中的下标，不是排序后的下标。
%
% 计算过程
%   1. 将方位角从小到大排序，同时用 n 保存它们在原数组中的位置。
%   2. 保留最小方位角作为当前基准 standard。
%   3. 后续角度若距当前基准小于 interval，就记录为待删；否则把它设为新基准。
%   4. 删除预分配向量中没有用到的尾部零值，只返回实际下标。
%
% 例：fai=[30 10 14 50], interval=10。排序后为 [10 14 30 50]，
%     14 距 10 只有 4，因此返回原数组中 14 的位置 3。
[m,n]=sort(fai);
% j 统计已经找到的待删台站数；先按最大可能长度分配，避免循环中反复扩容。
j=0;
del_index=zeros(length(fai),1);
standard=m(1);
for i=2:length(m)
    if m(i)-standard<interval
        %standard=m(i-1);
        % 该方位角太靠近最近一次保留的方位角，记录其原始下标。
        j=j+1;
        del_index(j)=n(i);
    else
        % 间隔足够，保留该台站，并以它作为下一轮比较基准。
        standard=m(i);
    end
end
% 预分配空间中 j 之后的元素仍为 0；MATLAB 删除语法把它们裁掉。
del_index(j+1:end)=[];
return
