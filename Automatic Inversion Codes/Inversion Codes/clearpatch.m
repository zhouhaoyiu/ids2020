% This fucntion is to find the rupture area, the number of rupture area can
% be user-defined.
function slip=clearpatch(slip,rat)
% CLEARPATCH 删除低于给定幅值门槛的孤立局部峰值。
%
% 输入
%   slip: 二维滑移量矩阵，每个元素对应一个子断层。
%   rat : 相对阈值，实际门槛为 max(slip(:))*rat。
%
% 输出
%   slip: 尺寸不变。满足“低于门槛且大于上、下、左、右四邻点”的元素被置 0。
%
% 为什么反复检查
%   某个峰值清零后，邻近元素的比较关系可能改变；while 循环一直执行到一整轮
%   没有元素再被清零。该函数只检查四邻域，不比较对角邻点。
ss=size(slip);
% 外围补一圈 0，原矩阵边缘也能用统一的四邻域下标检查。
slip0=zeros(ss+2);
slip0(2:end-1,2:end-1)=slip;  

mslip=max(slip(:));           % the maximum slip
% 低于 level 才属于可删除的小幅值候选点。
level=mslip*rat;

% num=1 表示上一轮发生过删除；进入下一轮继续检查。
num=1;
while num==1
    num=0;
    for i=1:ss(1)
        % slip0 比原 slip 多一圈边界，所以原第 i 行对应 i+1。
        i1=i+1;
        for j=1:ss(2)
            j1=j+1;
            % 当前点既低于幅值门槛，又严格高于四个邻点时，它是孤立的小峰值。
            if slip0(i1,j1)<level&&slip0(i1,j1)>slip0(i,j1)&&slip0(i1,j1)...
                    >slip0(i1,j)&&slip0(i1,j1)>slip0(i1+1,j1)&&...
                    slip0(i1,j1)>slip0(i1,j1+1)                            % If a subfault is larger than the slip of all four surrounding subfaults, it is set to zero
                slip0(i1,j1)=0;
                % 标记本轮发生了变化，完成扫描后还要再检查一轮。
                num=1;
            end
        end
    end
end

% 去掉计算时增加的零边框，恢复输入尺寸。
slip=slip0(2:end-1,2:end-1);
