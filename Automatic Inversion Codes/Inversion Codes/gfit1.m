function [cxy]=gfit1(x,y)
% GFIT1 计算 x 与 y 对应列之间的归一化内积。
%
% 输入
%   x, y: 尺寸相同的列向量或矩阵。
%
% 输出
%   cxy : 单列输入得到标量；多列输入得到 1*size(x,2) 行向量。
%         每个值通常位于 [-1,1]，1 表示波形同向且形状完全成比例。
%
% 与 gfit0 的关系
%   两者使用同一公式，都能处理矩阵。gfit1 额外区分单列和多列，
%   并对多列中的零能量列逐列置 0，避免某一列使全部结果提前返回。
%   该公式不减均值，属于归一化内积，不是去均值的 Pearson 相关系数。

% s(2) 是列数，用它判断当前是单条波形还是多条波形。
s=size(x);
if s(2)==1
    % 单列时用矩阵乘法直接得到三个标量：互内积和各自能量。
    xy=x'*y;
    xx=x'*x;
    yy=y'*y;
    if xx==0||yy==0
        % 任一输入全零都会使分母为 0，因此约定相似度为 0。
        cxy=0;
        return;
    else
        cxy=xy./sqrt(xx.*yy);
    end
    return
else
    % 多列时沿第一维求和，一次得到每一列的互内积和能量。
    xy=sum(x.*y);
    xx=sum(x.*x);
    yy=sum(y.*y);
    
    % 先预置为 0；只有分母非零的列才执行除法。
    cxy=zeros(size(xy));
    sxy=sqrt(xx.*yy);
    cxy(sxy==0)=0;
    
    nz=sxy~=0;
    % nz 只选择可安全相除的列，结果不会出现由 0/0 产生的 NaN。
    cxy(nz)=xy(nz)./sxy(nz);
end
