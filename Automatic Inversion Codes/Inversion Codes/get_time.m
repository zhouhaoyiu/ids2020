function [t]=get_time(dist,dep,phase,endnar)
% ========================================================================
%function [t]=get_time(dist,dep,phase,endnar);
%
%GET_TIME is to get travel times of P or S from the given time table,which
%may be created by other applications.In this case, the time table is based
%on the IASPEI91 model.
%
% input:
%        dist: epicentral distance(s) in degree, which may be a vector with
%              the same step. Ranges 0 to 99deg
%         dep: depth of the hypocenter in kilometers, which can be a vector
%              with the same step. Ranges 0 to 99km
%       phase: phase name of interests, the option for which is 'P' or 'S'
%      endnar: a parameter which controls the input style. If ANY value is
%              given to this argument. The dist or dep above does NOT have 
%              to have the same step.  
% output:
%          t: the travel time for the given phase,epicentral distances and
%             focal depths. If nargin<4, t is a matrix with rows of length
%             (dist) and columns of length(dep); If nargin==4, the DIST and
%             the DEP should have the same size, and t will have the same
%             size as the inputs.
%
% 中文说明
%   本函数从当前 MATLAB 路径中的走时表文件读取 P 或 S 波走时。
%   表格距离间隔为 0.01 度、深度间隔为 1 km；非整数网格点用线性或双线性插值。
%
% 两种输出方式
%   get_time(dist,dep,phase)
%     把 dist 中每个距离与 dep 中每个深度全部组合，输出 length(dist)*length(dep)。
%   get_time(dist,dep,phase,endnar)
%     endnar 的值不参与计算；它只触发逐元素配对模式。dist 与 dep 必须同尺寸，
%     输出 t 也保持该尺寸，t(i,j) 对应 dist(i,j)、dep(i,j) 这一对。
%
% 例：dist=[1,2]、dep=[5,10]，三参数模式计算 4 个距离—深度组合；
%     四参数模式只计算 (1,5) 和 (2,10) 两对，输出仍为 1*2。
%
% 输出 t 的单位由 tp/ts 数据表决定，本套表按秒使用。dist、dep 超出表格范围会下标越界。
%--------------------------------------------------------------------------
%           Zhang Yong and Xu Lisheng, 2006/09/08, Bejing
%==========================================================================

% Load the Time_table
if nargin==2 %Default phase is P
    % 只给距离和深度时默认查询 P 波。
    phase='P';
end
if upper(phase)=='P' %Load P time table
    % tp、tp0 是两套 P 波走时；0 被当作缺测值，合并时取两者中较早的有效到时。
    load tp
    load tp0
    tp(tp==0)=NaN;
    tp0(tp0==0)=NaN;
    tp(1:length(tp0),:)=min(tp(1:length(tp0),:),tp0);
    % 合并后仍缺测的位置恢复为 0，供后续插值读取。
    tp(isnan(tp))=0;
    t0=tp;
elseif upper(phase)=='S'%Load S time table
    % S 波表按与 P 波相同的规则合并 ts 和 ts0。
    load ts
    load ts0
    ts(ts==0)=NaN;
    ts0(ts0==0)=NaN;
    ts(1:length(ts0),:)=min(ts(1:length(ts0),:),ts0);
    ts(isnan(ts))=0;
    t0=ts;
else
    % phase 不是单字符 P/p/S/s 时终止，避免误用错误走时表。
    error('The phase can not be found!');
end

if nargin<4 %Get times for the given distance and depth range
    % 组合模式：预分配“距离数*深度数”的二维输出。
    sdist=length(dist);
    sdep=length(dep);
    t=zeros(sdist,sdep);
    dist_mod=mod(dist,0.01)*100;
    % dep_mod 是深度的小数部分，作为上下两个整数深度之间的插值权重。
    dep_mod=mod(dep,1);
    for i=1:sdist
        for j=1:sdep
            if dist_mod(i)==0
                % 距离正好落在 0.01 度网格上，只需沿深度方向作一次线性插值。
                t(i,j)=t0(round(dist(i)*100)+1,fix(dep(j))+1)...
                    *(1-dep_mod(j))...
                    +t0(round(dist(i)*100)+1,fix(dep(j))+2)*(dep_mod(j));
            else
                % 距离不在网格上：先在两个相邻距离行内插，再在两个相邻深度列之间插值。
                t(i,j)=((t0(fix(dist(i)*100)+2,fix(dep(j))+1)...
                    -t0(fix(dist(i)*100)+1,fix(dep(j))+1))*dist_mod(i)...
                    +t0(fix(dist(i)*100)+1,fix(dep(j))+1))*(1-dep_mod(j))...
                    +((t0(fix(dist(i)*100)+2,fix(dep(j))+2)-t0(fix(dist(i)...
                    *100)+1,fix(dep(j))+2))*dist_mod(i)...
                    +t0(fix(dist(i)*100)+1,fix(dep(j))+2))*(dep_mod(j));
            end
            %t(i,j)=t0(round(dist(i)*100)+1,round(dep(j))+1);
        end
    end
    return
else %Get times for the given pairs of distances and depths
    % 配对模式：每个 dist(i,j) 只与同位置 dep(i,j) 组合。
    sd=size(dist);
    t=zeros(sd);
    dist_mod=mod(dist,0.01)*100;
    dep_mod=mod(dep,1);
    for i=1:sd(1)
        for j=1:sd(2)
            if dist_mod(i,j)==0
                % 距离恰在表格网格上，只对这一对数据作深度线性插值。
                t(i,j)=t0(round(dist(i,j)*100)+1,fix(dep(i,j))+1)...
                    *(1-dep_mod(i,j))...
                    +t0(round(dist(i,j)*100)+1,fix(dep(i,j))+2)...
                    *(dep_mod(i,j));
            else
                % 对当前距离—深度对执行双线性插值，写回相同位置的 t(i,j)。
                t(i,j)=((t0(fix(dist(i,j)*100)+2,fix(dep(i,j))+1)...
                    -t0(fix(dist(i,j)*100)+1,fix(dep(i,j))+1))*dist_mod(i,j)...
                    +t0(fix(dist(i,j)*100)+1,fix(dep(i,j))+1))*(1-dep_mod(i,j))...
                    +((t0(fix(dist(i,j)*100)+2,fix(dep(i,j))+2)...
                    -t0(fix(dist(i,j)*100)+1,fix(dep(i,j))+2))*dist_mod(i,j)...
                    +t0(fix(dist(i,j)*100)+1,fix(dep(i,j))+2))*(dep_mod(i,j));
            end
            %t(i,j)=t0(round(dist(i,j)*100)+1,round(dep(i,j))+1);
        end
    end
end
return
%=====================================End==================================

