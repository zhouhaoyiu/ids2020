%% This function is to creat a sparse matrix
% 
% i=1,nx=[1,2];numx=[2,1]/3;
% i=grid(1)：e.g. i=12;nx=[11,12],numx=[1,2]/3
% i=other，e.g. i=2,nx=[1,2,3],numx=[1,1,1]/3
% mm=[1 1 2 2 2 3 3 3 ...... 12 12];
% nn=[1 2 1 2 3 2 3 4 ......10 11 12 11 12]
% num=[2 1 1 1 1 ...... 1 1 1 2]/3
% Result
%    1      2     3      4     5  ......    11     12
% 1 0.67  0.33
% 2       0.33  0.33  0.33
% 3             0.33  0.33   0.33
% 4                   0.33   0.33 0.33
% 5                          0.33 0.33 0.33 
% ...                              ... ...
% ...                              ... ...
% ...                              ... ...  0.33
% 11                                        0.33  0.33
% 12                                        0.33  0.67

function SM=smmat_gfz(grid,se)
% SMMAT_GFZ 建立一维三点滑动平均的稀疏矩阵。
%
% 输入
%   grid: 程序只使用 grid(1)，记为 N；其余元素不会参与计算。
%   se  : 可选二元素边界标志 [first,last]，默认 [1,1]。
%         标志为 1 时，边界自身权重为 2/3、相邻点为 1/3；
%         标志为 0 时，两项均为 1/3，相当于把缺失的外侧点当作 0。
%
% 输出
%   SM  : N*N 稀疏矩阵。SM*v 得到 v 的一次三点平滑结果。
%
% 例：N=3、se=[1,1] 时，SM 约为
%   [2/3 1/3 0; 1/3 1/3 1/3; 0 1/3 2/3]。
% 在 IDS_2020 中，分别对倾向和走向构造该矩阵，实现二维空间平滑。
if nargin==1
    % 默认保持边界行权重和为 1。
    se=[1,1];
end

% 每行最多有 3 个非零项，预分配行下标 mm、列下标 nn 和权重 num。
mm=zeros(grid(1)*3,1);nn=zeros(grid(1)*3,1);num=zeros(grid(1)*3,1);
dnum=0;
for i=1:grid(1)
    if i==1
        % 第一行只有自身和右邻点。
        if se(1)==0
            nx=[i,i+1];numx=[1,1]/3;%==sm(,3)
        else
        nx=[i,i+1];numx=[2,1]/3;%==sm(,3)
        end
    elseif i==grid(1)
        % 最后一行只有左邻点和自身。
        if se(2)==0
            nx=[i-1,i];numx=[1,1]/3;%==sm(,3)
        else
        nx=[i-1,i];numx=[1,2]/3;%==sm(,3)
        end
    else
        % 内部行对左邻、自身、右邻各取 1/3。
        nx=[i-1,i,i+1];numx=[1,1,1]/3;%==sm(,3)
    end
    for j=1:numel(nx)
        % 每个非零权重写成 sparse 所需的 (行,列,值) 三元组。
        dnum=dnum+1;
        mm(dnum)=i;
        nn(dnum)=nx(j);
        num(dnum)=numx(j);
    end
end
% 裁掉预分配但未使用的尾部元素，再一次性建立稀疏矩阵。
mm(dnum+1:end)=[];nn(dnum+1:end)=[];num(dnum+1:end)=[];
SM=sparse(mm,nn,num);

