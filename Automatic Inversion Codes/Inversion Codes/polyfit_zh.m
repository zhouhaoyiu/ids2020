function [wavepoly,L,A]=polyfit_zh(wave,numth,leng_add,option,dtime)
% fit the wave with legendre functions
% POLYFIT_ZH 用 Legendre 基函数拟合波形，并可向后外推。
%
% 输入
%   wave    : N*C 波形矩阵，第一维为采样点，每列独立拟合。
%   numth   : Legendre 多项式最高阶数，输出共有 numth+1 个基函数。
%   leng_add: 需要在原 N 点之后增加的外推点数。
%   option  : 可选。option=1 使用离散积分形式求系数；其他值用最小二乘。
%   dtime   : 可选的时间坐标，长度应为 N+leng_add。
%
% 输出
%   wavepoly: (N+leng_add)*C 的拟合及外推波形。
%   L       : (N+leng_add)*(numth+1) 的 Legendre 设计矩阵，每列对应一个阶次。
%   A       : (numth+1)*C 系数矩阵，满足 wavepoly=L*A。
%
% 调用规则
%   3 个参数：自动 option=2，建立等间隔 [-1,1] 坐标，最小二乘拟合。
%   4 个参数：建立等间隔坐标；只有这种调用可令 option=1。
%   5 个参数：使用传入 dtime，但代码会把 option 重设为 2。

if nargin==3
    % 缺省采用最小二乘，适合普通波形拟合。
    option=2;
end

% sw(1)=原采样点数 N，sw(2)=待拟合列数 C。
sw=size(wave);

% 每阶多项式占一列；行数包括原始点和待外推点。
L=zeros(sw(1)+leng_add,numth+1);

if nargin==4
    % 把完整拟合区间均匀映射到 [-1,1]；dt 是相邻归一化坐标间隔。
    dt=2/(sw(1)+leng_add-1);
    dtime=-1:dt:1;
else
    % 五参数方式保留自定义时间坐标的相对位置，并固定采用最小二乘。
    option=2;
    dtime=(dtime-dtime(1))/(dtime(end)-dtime(1))*2-1;
    
end

for i=0:numth
    % legendre(...,'sch') 返回 Schmidt 半归一化连带 Legendre 函数；这里只取 m=0 的首行。
    y=legendre(i,dtime,'sch');
    L(:,i+1)=y(1,:)';
end

if option==1
    % 离散积分方式当前只建立一列系数，适用于单列 wave。
    A=zeros(numth+1,1);
    
    % 梯形积分的两个端点权重为一半，因此同时把 wave 和 L 的首尾行除以 2。
    wave([1,end],:)=wave([1,end],:)/2;
    L([1,end],:)=L([1,end],:)/2;
    for i=1:numth+1
        waveL=wave.*L(:,i);
        % 用离散正交积分近似得到第 i-1 阶系数。
        A(i)=dt*(2*i+1)/2*sum(waveL);
    end
    wavepoly=L*A;
    
else
    % 最小二乘分支为每一列 wave 分别求一组系数。
    wavepoly=zeros(sw(1)+leng_add,sw(2));
    A=zeros(numth+1,sw(2));
    for i=1:sw(2)
        %A=cgls_zh0(L,wave,5e2);
        % 只用 L 的前 N 行拟合观测点；完整 L*A 自然给出后面的 leng_add 个外推点。
        A(:,i)=L(1:sw(1),:)\wave(:,i);
        %%A=inv(L'*L)*L'*wave;
        wavepoly(:,i)=L*A(:,i);
    end
end

