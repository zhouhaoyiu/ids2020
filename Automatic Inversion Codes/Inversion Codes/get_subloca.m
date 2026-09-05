function [loca1,dep]=get_subloca(fault,grid,sizegrid,source,sta_loca,epi)
%==========================================================================
%  delt=ruptime_nf(fault,grid,sizegrid,source,sta_loca,epi,phase);
%  This is a function to calculate relative times between fault consists of
%   several segments and stations.
% -------------------------------------------------------------------------
% 
%  sgement 1
%  ----------
%            \
%             \ segment 2
%              \
%               ----------------- .......
%                  segment 3 .......
%
%   Input:
%          fault: [strike,dip]. both strike and dip are vectors with size
%                  of [num of segments,1]
%          grid:  [Nstrike, Ndip]. both strike and dip are vectors with size
%                  of [num of segments,1]
%      sizegrid:  describe the size of grid. the units are km
%        source:  a row vector, and has 3 elements. source(1) is the index
%                of the segments source locats in. source(2) is the index
%                of grid for strike direction and source(3) is the index
%                of grid for dip direction in the segment--source(1).
%      sta_loca: locations of stations. [Lat,Long], bot are vectors
%           epi: [lat,long]. the epicenter
%         phase: only can be 'P' or 'S'!
%
%   Output:
%        loca1: nsub*2，各子断层中心的 [纬度,经度]。
%          dep: nsub*1，各子断层中心深度，单位 km。
%
% 中文说明
%   fault 每行是一段断层的 [strike,dip,...]；grid 每行为 [Nstrike,Ndip]；
%   sizegrid 每行为 [strike_km,dip_km]；source=[segment,strike_index,dip_index]。
%   epi=[lat,lon] 给出震中。sta_loca 在当前执行路径中没有被使用。
%
% 数据怎样生成
%   先从包含震源的断层段出发，用 rup_locanew 求该段四角坐标和深度；
%   再用 interm 在四角之间插值出全部子断层中心。随后向段号减小和增大的两侧推进，
%   相邻段共用边界，程序临时多建一列，插值后删去重叠列。
%   最终顺序为逐段排列，每段内部按 MATLAB 列优先顺序展开。
%
% 当前限制
%   grid(1,2)==1 的特殊分支调用本目录中不存在的 ruptime_nfnew，并引用未定义变量 phase，
%   而且没有给 dep 赋值；单个倾向网格的输入会在该分支失败。常规多网格分支不经过它。
%
%  This is a program to calculate relative times for several segments. To
%  kown more details about this, you can see ruptime which can calculate
%  relative times for 3-d fault. And for a simple plane fault, This 
%  function has the same results with ruptime.
%  
%  Two function used in this function is important: rup_locanew and get_time
% -------------------------------------------------------------------------
%                                                  Zhang Yong 
%                                2007/04/16/22:46  Peking University
%==========================================================================

if grid(1,2)==1
    % 旧的单倾向网格兼容路径：临时复制成两行，再删除隔行结果。
    grid(:,2)=2;
    [delt,loca1]=ruptime_nfnew(fault,grid,sizegrid,source,sta_loca,epi,phase);
    delt(2:2:end,:)=[];
    return
end

ss=size(fault);
% 每段子断层数为 Nstrike*Ndip，求和得到总输出行数。
nsub=sum(prod(grid'));

%find the location of the source, in which segment
index_s=source(1);

% 预分配全部子断层的位置和深度；sumgrid 给出每段在输出中的起止偏移。
loca=zeros(nsub,2); 
dep=zeros(nsub,1);  
sumgrid=[0;cumsum(prod(grid')')];
% find location of the segment from index_s to 1:
for i=index_s:-1:1
    if i==index_s
        % 源段以真实震中和震源网格下标确定四角几何。
        [loca0,dep0]=rup_locanew(fault(i,:),grid(i,:),sizegrid(i,:),[epi,source(2:3)]);
        locanext=loca0;% for next
        locasource=loca0;% for along the strike direction
        lat=reshape(loca0(:,1),[2,2])';
        long=reshape(loca0(:,2),[2,2])';
        % 2*2 四角数组分别加密到 Ndip*Nstrike，再按列展开写入该段输出区间。
        lat=interm(lat,[grid(i,2)-1,grid(i,1)-1]);
        long=interm(long,[grid(i,2)-1,grid(i,1)-1]);
      
        loca(sumgrid(i)+1:sumgrid(i+1),1)=lat(:);
        loca(sumgrid(i)+1:sumgrid(i+1),2)=long(:);
        dep0=reshape(dep0,[2,2])';
        % 深度使用与经纬度相同的网格插值和展开顺序。
        dep0=interm(dep0,[grid(i,2)-1,grid(i,1)-1]);
        dep(sumgrid(i)+1:sumgrid(i+1))=dep0(:);
    else
        % 源段之前的段以后一段的第一个上边界角点作为连接点。
        epi_o=locanext(1,:);%[lat,long]
        % 走向临时增加一个网格，以便生成与相邻段共用的连接边界。
        [loca0,dep0]=rup_locanew(fault(i,:),grid(i,:)+[1,0],sizegrid(i,:),[epi_o,grid(i,1)+1,1]);%add a vector
        locanext=loca0;% for next
        lat=reshape(loca0(:,1),[2,2])';
        long=reshape(loca0(:,2),[2,2])';
        lat=interm(lat,[grid(i,2)-1,grid(i,1)]); %add a vector
        % 删除临时增加的最后一列，避免相邻断层段重复保存边界。
        lat(:,end)=[]; % delete the added vector: the last one
        long=interm(long,[grid(i,2)-1,grid(i,1)]);%add a vector
        long(:,end)=[]; % delete the added vector: the last one
        loca(sumgrid(i)+1:sumgrid(i+1),1)=lat(:);
        loca(sumgrid(i)+1:sumgrid(i+1),2)=long(:);
        dep0=reshape(dep0,[2,2])';
        dep0=interm(dep0,[grid(i,2)-1,grid(i,1)-1]);
        dep(sumgrid(i)+1:sumgrid(i+1))=dep0(:);
    end
end
% find location of the segment from index_s+1 to ss(1): ss(1)---number of
% segments
if index_s<ss(1)
    for i=index_s+1:1:ss(1)
        if i==index_s+1
            % 紧邻源段的后一段从源段第二个上边界角点开始。
            epi_o=locasource(2,:);
        else
            % 更后面的段继续使用上一段的第二个连接角点。
            epi_o=locanext(2,:);%[lat,long]
        end
        [loca0,dep0]=rup_locanew(fault(i,:),grid(i,:)+[1,0],sizegrid(i,:),[epi_o,1,1]);%add a vector
        locanext=loca0;% for next
        lat=reshape(loca0(:,1),[2,2])'; 
        long=reshape(loca0(:,2),[2,2])';
        lat=interm(lat,[grid(i,2)-1,grid(i,1)]); %add a vector
        % 后续段删除临时增加的第一列，把共享边界只留在前一段。
        lat(:,1)=[]; % delete the added vector: the first one
        long=interm(long,[grid(i,2)-1,grid(i,1)]);%add a vector
        long(:,1)=[]; % delete the added vector: the first one
        loca(sumgrid(i)+1:sumgrid(i+1),1)=lat(:);
        loca(sumgrid(i)+1:sumgrid(i+1),2)=long(:);
        dep0=reshape(dep0,[2,2])';
        dep0=interm(dep0,[grid(i,2)-1,grid(i,1)-1]);        
        dep(sumgrid(i)+1:sumgrid(i+1))=dep0(:);
    end
end

% 对外使用 loca1 这个名称返回完整位置矩阵。
loca1=loca;

