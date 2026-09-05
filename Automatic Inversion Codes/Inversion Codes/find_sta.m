function [ndex]=find_sta(tazim,ndex,segma1,segma2)
%% identify a partner station by requiring the two stations to have similar take-off angles and nearly opposite azimuths.
% Input:     tazim: The epicentral distance and azimuth of the stations
%             ndex: The line number point to the tazim 
%              X/Y: X,the epicentral distance;Y,the azimuth
%    segma1/segma2: The standard deviations of X and Y

% Output:     ndex: The line number point to the tazim, this is the partner
% station
% -------------------------------------------------------------------------
% 中文说明
%   tazim 每行对应一个台站：[距离, 方位角]。
%   输入 ndex 指向基准台站；输出 ndex 改为评分最高的配对台站下标。
%   segma1、segma2 分别控制允许的距离差和方位角差，数值越大，限制越宽。
%
% 计算过程
%   1. 目标距离 mu1 取基准台站距离。
%   2. 目标方位角 mu2 取基准方位角的反方向，即加或减 180 度。
%   3. 对所有台站计算二维高斯权重：距离和目标越近、方位角和反方向越近，权重越大。
%   4. 返回最大权重所在行。
%
% 输出原因
%   这样选出的台站与基准台站距离相近、方位大致相反，可形成近似对称的台站对。
%   当前代码直接计算角度差，没有把 0/360 度视为相邻；靠近该边界时应留意结果。

% 二维高斯分布的归一化系数。它对所有候选台站相同，不影响最大值的位置。
segma=1/(2*pi*segma1*segma2);
mu1=tazim(ndex,1);                                                       
if tazim(ndex,2)<=180
    % 方位角位于前半圆时，加 180 度得到反方向。
    mu2=tazim(ndex,2)+180;                                               
else
    % 方位角位于后半圆时，减 180 度避免超过 360 度。
    mu2=tazim(ndex,2)-180;
end
X=tazim(:,1); % The epicentral distance
Y=tazim(:,2); % The azimuth
% para1、para2 是距离差和方位角差各自对指数的贡献。
para1=((X-mu1).^2/segma1^2)/2;para2=((Y-mu2).^2/segma2^2)/2;
% 差异越小，负指数越接近 0，funa 越大。
funa=segma*exp(-para1-para2);
[~,ndex]=max(funa);
% 若出现并列最大值，max 可能返回首个位置；这里明确取第一个。
ndex=ndex(1);
end
