function acc=cor_baseline(acc,len1)
% This function is to remove the baseline shifts
% COR_BASELINE 通过位移趋势拟合校正一条加速度记录的基线漂移。
%
% 输入
%   acc : N*1 加速度序列。代码按单列数据设计。
%   len1: 用于第一次二次趋势拟合的前段样点数，必须位于 1 到 N 之间。
%
% 输出
%   acc : 校正后加速度。通常为 N*1；若 dex 落在 (t1,t2) 内，当前分支只返回 dex:N 段。
%
% 数据怎样向前传递
%   acc0 --累加--> vel0 --累加--> dis0
%   dis0 --去除前段外推的二次趋势--> dis1
%   dis1 --一阶差分--> vel1 --一阶差分--> acc1
%   acc1 的累计平方能量确定 30%%、70%% 能量时刻 t1、t2；
%   再由 t2 之后的位移拟合残余速度基线，最后差分回 acc。
%
% 为什么用 cumsum 和 diff
%   基线中的很小常量或斜率在积分后会放大成明显的速度、位移漂移，较容易拟合；
%   校正位移或速度后再差分，可以回到加速度域。代码没有显式乘采样间隔，
%   因而中间 vel/dis 是“每采样点”的相对量，而不是带严格物理单位的积分结果。
% 保存原始输入，后面的计算都在副本上进行。
acc0=acc;
% 第一次累加得到相对速度，第二次累加得到相对位移。
vel0=cumsum(acc0);
dis0=cumsum(vel0);

% 用最前 len1 点拟合二次 Legendre 趋势，并外推到整条记录长度。
syndis0=polyfit_zh(dis0(1:len1),2,size(acc,1)-len1,2);

% 从累计位移中减去第一条基线趋势，得到初步校正位移。
dis1=dis0-syndis0;

% 对位移作一阶差分恢复速度，再对速度差分恢复加速度。
vel1=[dis0(1);diff(dis1)];
acc1=[vel0(1);diff(vel1)];

% 累计平方值描述校正后信号能量增长；除以末值后范围变为 0 到 1。
energy=cumsum(acc1.^2);
energy=energy/energy(end);
% t1、t2 分别取最接近累计能量 30%% 和 70%% 的采样点。
[~,t1]=min(abs(energy-0.30));
[~,t2]=min(abs(energy-0.70));

% 对 t2 之后的初步校正位移再作二次拟合，ss 是尾段拟合位移。
[x,ss]=fit_dex(dis1(t2:end),2);

% t2 之前补 0，使尾段拟合恢复成与整条记录等长的趋势；一阶差分得到速度改正量。
syndis1=[zeros(t2-1,1);ss(:)];
covel1=[0;diff(syndis1)];
% dex 是速度改正量绝对值最小的位置，可视为该趋势接近零的连接点。
[~,dex]=min(abs(covel1));
if dex>t1&&dex<t2
    % 连接点落在主要能量区间内：连接点以前被舍去，之后直接减去拟合的速度漂移。
    vel2=vel1(dex:end)-covel1(dex:end);
else
    % 连接点不合适时，构造分段速度：前段保留，中段线性过渡，后段减去拟合漂移。
    vel2=zeros(size(vel1));
    
    vel2(1:t1)=vel1(1:t1);

    vel2(t1:t2+1)=vel1(t1:t2+1)-(0:(covel1(t2+1))/(t2-t1+1):covel1(t2+1))';
    vel2(t2+1:end)=vel1(t2+1:end)-covel1(t2+1:end);
end

% 最后对校正速度作差分，得到校正加速度。
% 注意：若进入 dex 分支，vel2 从 dex 开始，输出会短于输入；这是当前源码的实际行为。
acc=[vel2(1);diff(vel2)];
