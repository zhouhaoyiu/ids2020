function [patch,sumslip]=slippatch_new(slip,num)
% To find the first 'num' slip patches, the old codes please see 'slippatch'
% SLIPPATCH_NEW 把二维滑移分布按主要峰值逐块拆分。
%
% 输入
%   slip: ndip*nstrike 非负滑移矩阵。
%   num : 最多提取的滑移块数；省略时上限取 numel(slip)。
%
% 输出
%   patch  : ndip*nstrike*K，第 i 层是第 i 次提取的滑移块；实际 K 可小于 num。
%   sumslip: K*1，每个滑移块全部元素之和。
%
% 当前实际执行的过程
%   1. 用接近 1 的相对阈值调用 clearpatch，保留当前最大峰所属的主要区域。
%   2. 把该区域写入 patch(:,:,i)，并从剩余 slip 中减掉。
%   3. 对剩余滑移重复，直到达到 num 或剩余矩阵全零。
%
% 注意：第 20 行 return 之后是旧版 slippatch 算法，公共函数运行时不会到达；
% 本文件保留它仅供追溯。下面两个局部函数只服务于这段旧算法。
if nargin==1
    % 未限制块数时，理论上最多每个网格单元形成一块。
    num=numel(slip);
end
% 阈值略低于全局最大值，避免浮点舍入把最大值本身判为“低于阈值”。
rat=1-1e-10;
% 先按最大块数分配，提前结束时再裁掉空层。
patch=zeros([size(slip),num]);
sumslip=zeros(num,1);
for i=1:num
    if max(slip(:))==0
        % 没有剩余正滑移时，删除未使用的输出层并结束。
        patch(:,:,i:end)=[];
        sumslip(i:end)=[];
        return
    end
    % clearpatch 从当前剩余滑移中分离本轮主要区域。
    slip0=clearpatch(slip,rat);
    sumslip(i)=sum(slip0(:));
    patch(:,:,i)=slip0;
    % 已提取区域从残量中扣除，下一轮寻找下一个滑移块。
    slip=slip-slip0;
end
return
%==================
% below is the codes of 'slippatch'
% 以下为不可到达的旧实现：先建立每个网格点的四邻域，再从全局最大值沿严格下降方向扩展。

if nargin==1
    num=numel(slip);
end

grid=size(slip);

% 四个方向依次为左、右、上、下的行列偏移。
s1vec=[0,0,-1,1];
s2vec=[-1,1,0,0];
dex=cell(numel(slip),1);
for i=1:grid(1)
    for j=1:grid(2)
        % dex{线性下标} 保存该点位于网格内的四邻点线性下标。
        [dex{(j-1)*grid(1)+i}]=vecslip(i+s1vec,j+s2vec,grid);
    end
end

numslip=1:numel(slip);
s12=numslip./grid(1);
s1=mod(numslip,grid(1));s1(s1==0)=grid(1);
s2=ceil(s12);

patch=zeros([grid,num]);
pdex=cell(num,1);
mslip=zeros(num,1);
for i=1:num
    if sum(abs(slip(:)))==0
        patch(:,:,i:end)=[];
        pdex(i:end)=[];
        mslip(i:end)=[];
        return
    end
    [slip0,pdex0]=sub_slippatch(slip,s1,s2,grid,dex);
    % 旧算法也在每轮从残量中扣除已找到的区域。
    slip=slip-slip0;
    
    mslip(i)=max(slip0(:));
    patch(:,:,i)=slip0;
    pdex{i}=pdex0;
end

function [slip0,pdex]=sub_slippatch(slip,s10,s20,grid,dex)
% SUB_SLIPPATCH（旧版局部函数）从全局最大值出发，沿数值严格变小的四邻点扩展。

pdex=[];

slip0=zeros(grid);
[m,n]=max(slip(:));
% n 是本轮全局最大值位置，先把它放入滑移块。
pdex=[pdex,n];

dex1=dex{n};
vec=slip(dex1);
dexnew=dex1(vec<m);
% 只接收比当前最大值低的邻点，因此搜索从峰顶向外下降。

pdex=[pdex,dexnew];


dexnew0=dexnew;
dexiter=dexnew;


while ~isempty(dexiter)    
    dexiter=[];
    for i=1:length(dexnew0)

        dex0=dexnew0(i);

        dex1=dex{dex0};
        vec=slip(dex1);
        
        dexget=dex1(vec<slip(dex0));
        % 每一层继续寻找比当前点更低的邻点，直到没有新位置。
        dexiter=[dexiter,dexget];
    end
    dexnew0=dexiter;
    pdex=[pdex,dexnew0];
end
slip0(pdex)=slip(pdex);


function [dex]=vecslip(s1,s2,grid)
% VECSLIP（旧版局部函数）把候选四邻点裁到网格内，并转换为 MATLAB 线性下标。
xx=s1<1|s1>grid(1)|s2<1|s2>grid(2);

% MATLAB 列优先线性下标公式：(列-1)*行数+行。
dex0=(s2-1)*grid(1)+s1;
dex=dex0(~xx);

return
%========================== END ================================


