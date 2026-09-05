function [loca,dep,epiloca,epidep]=rup_locanew(fault,grid,sizegrid,epi)
%==========================================================================
%[loca,dep,epiloca,epidep]=rup_loca0(fault,grid,sizegrid,epi);
% rup_loca0 is used in ruptime. it is used for calculating the positions of 
% 4 corners for a plane fault.
% 
%         strike:
%        -----------------------> 
%       /   1                 2
%      / 
%     /   3                 4           1,2,3,4 is the 4 corners 
%    dip
%
% as an example:
%   clear all
%   fault=[90,90,0];grid=[3,3];sizegrid=[100,100];epi=[36,90,2,2];
%--------------------------------------------------------------------------
%  Input:
%       fault: [stike,dip rake]
%        grid: also is the 'sizemat' in other functions.
%              num in [strike,dip,rake] for 3-d model
%              and num in [strike,dip] for 2-d model
%    sizegrid: the size in [strike,dip,rake] for 3-d model
%              and size in [strike,dip] for 2-d model. Units: km
%         epi: the first 2 elements are location of the epicenter: lat and
%              long, the next others are source position in the grids for
%              both 2-d and 3-d models
%
%  Output:
%       loca: a matrix has a size of [4,2]. first column is lat, second is
%             long. from top to bottom, are the locations of for corners:
%             1,2,3,4, respectively
%        dep: depthes for corners:1,2,3,4, respectively
%    epiloca: same as epi(1:2)
%     epidep: depth of the hypocenter
%
% 中文说明
%   fault=[strike,dip,rake]，其中本函数只使用 strike 和 dip。
%   grid=[Nstrike,Ndip]，sizegrid=[strike_km,dip_km]。
%   epi=[lat,lon,source_strike_index,source_dip_index]。
%
% 输出怎样形成
%   1. 在断层面内计算震源到四个角点的距离 f 和夹角 fai。
%   2. 用倾角 dip 把断层面距离、方向投影到水平面。
%   3. 根据走向 strike 把四个水平向量换成纬度、经度增量，得到 loca(1:4,:)。
%   4. 用沿倾向距离乘 sin(dip) 得到四角深度 dep 和震源深度 epidep。
%
% 坐标近似
%   水平距离用 111.2 km/度换算纬度，经度再除以 cos(epi_lat)。
%   这是局部平面近似，断层范围很大或靠近极区时误差会增大。
% -------------------------------------------------------------------------
%                                  Zhang Yong, Chen Yun-Tai, Xu Li-Sheng
%                                        2006/04,Peking University 
%                              repaired 2007/01/30/01:00,Peking University 
% =========================================================================

if epi(3)<0||epi(3)>grid(1)||epi(4)<0||epi(4)>grid(2)
    % 震源网格索引超出断层面时，四角几何没有定义，立即停止。
    error('position of hypocenter should be correctly on the fault plane!')
end

% strike 保持度；dip 转为弧度供 sin、cos、atan 使用。
strike=fault(1);
dip=fault(2).*pi./180;
fai=zeros(4,1);
f=zeros(4,1);
sita=zeros(4,1);
loca=zeros(4,2);

%calculate the angles:
if epi(3)==1
    % 震源位于第一条走向网格时，通往相应两个角点的断层面夹角取 pi/2，避免除以 0。
    fai(1)=pi./2;
    fai(3)=pi./2;
else
    % atan(倾向距离/走向距离) 得到断层面内指向角点的夹角。
    fai(1)=atan(abs((sizegrid(2).*(epi(4)-1))./(sizegrid(1).*(epi(3)-1))));
    fai(3)=atan(abs((sizegrid(2).*(grid(2)-epi(4)))./(sizegrid(1).*(epi(3)-1))));
end
if epi(3)==grid(1)
    % 震源位于最后一条走向网格时，对另一侧两个角点同样处理除零情况。
    fai(2)=pi./2;
    fai(4)=pi./2;
else
    fai(2)=atan(abs((sizegrid(2).*(epi(4)-1))./(sizegrid(1).*(grid(1)-epi(3)))));
    fai(4)=atan(abs((sizegrid(2).*(grid(2)-epi(4)))./(sizegrid(1).*(grid(1)-epi(3)))));
end

%calculate the distances for the four conners to the epicenter
% 勾股关系计算断层面内震源至四角的距离，单位 km。
f(1)=sqrt((sizegrid(1).*(1-epi(3))).^2+(sizegrid(2).*(1-epi(4))).^2);
f(2)=sqrt((sizegrid(1).*(grid(1)-epi(3))).^2+(sizegrid(2).*(1-epi(4))).^2);
f(3)=sqrt((sizegrid(1).*(1-epi(3))).^2+(sizegrid(2).*(grid(2)-epi(4))).^2);
f(4)=sqrt((sizegrid(1).*(grid(1)-epi(3))).^2+(sizegrid(2).*(grid(2)-epi(4))).^2);

%project the distances to horizontal plane (with depth of 1 or 2)
% 扣除倾向造成的垂向分量，得到每条角点向量的水平投影长度。
f=f.*cos(asin(sin(dip).*sin(fai)));

%project the angles to horizontal plane (with depth of 1 or 2)
% 倾斜断层投影后，断层面内角 fai 变为水平面内角。
fai=atan(cos(dip).*tan(fai));

%the horizontal angles for the orientation of each distance vector
% 将相对断层走向的角度换成地理坐标中使用的水平向量方向 sita。
sita(1)=(strike+fai(1)*180/pi-270).*pi./180;
sita(2)=(strike-fai(2)*180/pi-270).*pi./180;
sita(3)=(strike-fai(3)*180/pi-270).*pi./180;
sita(4)=(strike+fai(4)*180/pi-270).*pi./180;

%get the positions for the 4 corner projected into the horizontal planes
% 纬度增量约为 north_km/111.2；经度每度长度再乘 cos(latitude)，所以经度增量需除该项。
loca(1,:)=f(1).*[-sin(sita(1)),cos(sita(1))./cosd(epi(1))]./111.2+epi(1:2);
loca(2,:)=f(2).*[sin(sita(2)),-cos(sita(2))./cosd(epi(1))]./111.2+epi(1:2);
loca(3,:)=f(3).*[-sin(sita(3)),cos(sita(3))./cosd(epi(1))]./111.2+epi(1:2);
loca(4,:)=f(4).*[sin(sita(4)),-cos(sita(4))./cosd(epi(1))]./111.2+epi(1:2);

%depth of 4 corners
% 上边界角点位于半个倾向网格深度，下边界位于 (Ndip-0.5) 个网格深度。
dep1=sizegrid(2)./2;
dep2=sizegrid(2)./2;
dep3=(grid(2)-1/2).*sizegrid(2);
dep4=(grid(2)-1/2).*sizegrid(2);
dep=[dep1;dep2;dep3;dep4].*sin(dip);

%location and delpth of the hypocenter 
% 震中坐标直接沿用输入；震源深度由震源所在倾向网格中心计算。
epiloca=epi(1:2);
epidep=(epi(4)-1/2).*sizegrid(2).*sin(dip);
return
%==================================end=====================================

