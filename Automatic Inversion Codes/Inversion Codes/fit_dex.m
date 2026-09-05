function [x,syn]=fit_dex(dot,num,flag)
% FIT_DEX 用多项式拟合一段序列，并返回系数和拟合曲线。
%
% 输入
%   dot : 待拟合序列，通常为列向量。
%   num : 多项式次数。
%   flag: 可选占位参数。程序不读取它的值；只要提供第三个参数，就进入第二种拟合写法。
%
% 输出
%   x   : 从常数项到高次项排列的系数。
%   syn : 在原采样点上的拟合值，尺寸与 dot 相同。
%
% 两个调用方式
%   fit_dex(dot,num)      自行建立幂函数矩阵 G，用最小二乘 G\dot 求系数。
%   fit_dex(dot,num,flag) 先减去首点，再调用 MATLAB polyfit，最后把首点加回。
%   两种分支的拟合约束不同，不能仅把 flag 当作普通数值选项理解。
if nargin==2
    % 后面 j-1 从 0 到原 num，因此先加 1 表示所需的系数个数。
    num=num+1;

    % i 是采样点编号，j 是幂次列编号；G 每列依次为 i^0、i^1、...。
    [j,i]=meshgrid(1:num,1:length(dot));
    G=i.^(j-1);
    % 反斜杠求最小二乘解，x 此时为列向量；syn 是 G*x 生成的拟合序列。
    x=G\dot; 
    syn=G*x;
    % 转成行向量，统一本函数的系数输出方向。
    x=x';
else
    % 保存首点并把整条序列减去它，使 polyfit 拟合相对首点的变化量。
    firdot=dot(1);
    dot=dot-firdot;
    dex=(1:length(dot))';
    % polyfit 返回从最高次项到常数项排列的 num+1 个系数。
    [x] = polyfit(dex,dot,num);
    
    % 按 polyfit 系数逐项重建拟合曲线，随后恢复先前减去的首点。
    syn=zeros(size(dot));
    for i=1:length(x)
        syn=syn+x(i)*dex.^(num-i+1);
    end
    syn=syn+firdot;
    % 反转系数，使输出 x 改为常数项、一次项、...、最高次项的顺序。
    x=x(end:-1:1);
end
