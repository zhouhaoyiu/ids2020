function subdist=faultdist(grid,gridsize,source,flag)
%==========================================================================
%  subdist=faultdist(grid,gridsize,source)
%  here is a function to get the distances of each subfault to the
%  hypocenter
%--------------------------------------------------------------------------
% Input
%       grid: the number of subfault in [dip,strike] direction
%   gridsize: the size of the subfault, in [dip,strike] direction
%     source: the subfault hypocenter located in
% Output
%   subdist: the distances of each subfault to the hypocenter
%
% 中文说明
%   grid=[ndip,nstrike]，gridsize=[dip_km,strike_km]，source=[row,column]。
%   输出 subdist 为 ndip*nstrike，每个元素单位为 km。
%
% 两种计算方式
%   faultdist(grid,gridsize,source)
%     把震源放在 source 指定子断层的中心，计算所有子断层中心到该中心的直线距离。
%   faultdist(grid,gridsize,source,flag)
%     flag 的具体数值不参与运算；只要传入第 4 个参数，每个方向就用
%     “目标子断层中心距 source 中心的距离 - 半个网格尺寸”，再把负值截为 0。
%     这可理解为目标中心到震源子断层矩形区域的最短距离；相邻单元在该方向为半格。
%
% 例：grid=[2,3]、gridsize=[5,10]、source=[1,2]，无 flag 时
%   source 位置距离为 0；同行相邻列距离为 10 km；下一行同列距离为 5 km。
%--------------------------------------------------------------------------
%          Zhang Yong, 2012-03-14 09:46, GFZ, Potsdam
%==========================================================================
%[ii,jj]=meshgrid(1:grid(2),1:grid(1));

if nargin==3
    % ii 保存列号、jj 保存行号；两者尺寸均为 grid。
    ii=((1:grid(2))'*ones(1,grid(1)))';jj=(1:grid(1))'*ones(1,grid(2));
    % 行、列索引差分别乘对应子断层尺寸，再由勾股关系合成距离。
    subdist=sqrt((abs(jj-source(1))*gridsize(1)).^2+(abs(ii-source(2))*gridsize(2)).^2);
end

if nargin>3
    % 第四参数分支逐个计算目标子断层中心到震源子断层区域的最短距离。
    subdist=zeros(grid);
    for i=1:grid(1)
        for j=1:grid(2)
            di=(abs(i-source(1))-0.5).*gridsize(1);
            dj=(abs(j-source(2))-0.5).*gridsize(2);
            % 负值表示两个子断层在该方向相接或重叠，最短间隔应记为 0。
            di=max(di,0);dj=max(dj,0);
            dd=sqrt(di.^2+dj.^2);
            subdist(i,j)=dd;
        end
    end
end
return
