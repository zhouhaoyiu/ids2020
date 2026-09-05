function da=da_zh(loca,epi,flag)
%==========================================================================
%  da=da_zh(loca,epi,flag)
%  to calculate the distance and azimuth
%--------------------------------------------------------------------------
%  Input
%    loca: [lat,long], the locations of stations
%     epi: epicenter
%   if nargin<3, calculate with ellipse model, else with spherical model
% Output
%      da: distance and azimuth,[distance, azimuth]
%
% 中文说明
%   loca 为 N*2 [纬度,经度]，epi 可为 1*2，也可与 loca 逐行配对。
%   输出第一列是距离(km)，第二列是从 epi 指向 loca 的方位角(度、北起顺时针)。
%   第三个参数决定计算方式：
%     不提供 flag：调用 Mapping Toolbox 的 WGS84 大圆距离；
%     flag==1：先按球面求角距离，再乘平均地球半径换算 km；
%     其他 flag：使用本文件展开的球面公式，并额外输出反方位角第三列。
%--------------------------------------------------------------------------
%       Zhang Yong, 2012-02-09 12:49, GFZ, Potsdam
%==========================================================================
if nargin<3
    % Supply WGS84 earth ellipsoid axis lengths in kilometers:
    a = 6378.137; % definitionally
    b = 6356.75231424518; % computed from WGS84 earth flattening coefficient definition
    c = sqrt(a*a-b*b);
    % e 是 WGS84 第一偏心率；[a,e] 把 distance 的长度单位设为 km。
    e=c./a;
    
    [dist,azi]=distance('gc',epi,loca,[a,e]);
    %[d,a]=distance(epi,loca,[a,e]);
    da=[dist,azi];
elseif flag==1
    % distance 默认返回角距离(度)，乘 6371*pi/180 换算为球面弧长(km)。
    [dist,azi]=distance(epi,loca);

    da=[dist*6371*pi/180,azi];
else
    % 手写球面分支先把经纬度从度转换为弧度。
    k=pi./180;
    slon=loca(:,2).*k;
    elon=epi(:,2).*k;
    % slat=atan(0.9933.*tan(statlocat(:,1).*k));
    % elat=atan(0.9933.*tan(epilocat(:,1).*k));
    slat=atan(tan(loca(:,1).*k));
    elat=atan(tan(epi(:,1).*k));
    
    %Calculation of distance in rad
    % 球面余弦定理给出两点的大圆圆心角 drad。
    drad=acos(sin(slat).*sin(elat)+cos(slat).*cos(elat).*cos(slon-elon));
    
    %calculation of the azimuth
    saz=acos((sin(elat)-sin(slat).*cos(drad))./(cos(slat).*sin(drad)));
    eaz=acos((sin(slat)-sin(elat).*cos(drad))./(cos(elat).*sin(drad)));
    
    % conversion rad-deg
    delta=drad./k;
    
    epicaz=real(eaz./k);
    stataz=real(saz./k);
    %distance in km
    dist=delta*6371*pi/180;
    
    % make sure the azimuth is measured from N clockwise
    % acos 只能给出 0 到 180 度，利用经度差符号把方位补到 0 到 360 度。
    x=sin(slon-elon);
    y=sin(elon-slon);
    xx=(x>0);
    stataz(xx)=360-stataz(xx);
    yy=(y>0);
    epicaz(yy)=360-epicaz(yy);
    %Change variables
    temp=epicaz;
    epicaz=stataz;
    stataz=temp;
    % 手写分支返回 [距离,正方位角,反方位角] 三列。
    da=[dist,stataz,epicaz];
end
return
