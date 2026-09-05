function [grid,source]=subfaults_making(slipb,grid,source)
% This function is to adjust the size of fault plane
% SUBFAULTS_MAKING 根据归一化滑移分布决定下一轮反演的断层网格大小和震源位置。
%
% 输入
%   slipb : grid(1)*grid(2) 相对滑移矩阵；代码把等于 0 的边缘视为空白。
%   grid  : [ndip,nstrike] 当前倾向、走向子断层数。
%   source: [dip_index,strike_index] 当前震源网格下标。
%
% 输出
%   grid  : 下一轮建议的网格数。
%   source: 网格移动或扩展后，震源在新网格中的下标。
%
% 判断方法
%   对行和列分别从两端累计滑移。累计和持续为 0 的长度，就是该边缘连续空白网格数。
%   空白超过 2 格时收缩或平移；滑移碰到边缘时，根据边缘最大滑移/0.2 向外增加网格。
%   公式中的 0.2、4、7、8、9 是原程序固定的经验余量，作用是让震源和有效滑移
%   与新边界之间保留若干网格；它们不是由输入数据自动估计的参数。

% The adjustment of teh fault plane in the strike direction
% 第一部分实际处理倾向（矩阵行）方向。
% da 从顶部累计，db 从底部累计；num1、num2 是顶部和底部连续空白行数。
da=cumsum(sum(slipb,2));                           
db=cumsum(flipud(sum(slipb,2)));                    % Flip it 90 degrees up and down, and sum it up
num1=sum(da(:)==0);num2=sum(db(:)==0);

if num1>2&&num2<=2;                                 % Zero appears at the top of the dip direction and the fault moves downward
    % 顶部空白较多、底部滑移靠边：裁顶部，并按底边滑移强度决定是否向下扩展。
    num=sum(da(:)==0);
    if source(1)-num<=0&&db(1)~=0;
        grid(1)=grid(1)-source(1)+ceil(max(slipb(end,:))/0.2)+8;
        source(1)=5;                                % Make sure there are two rows above the epicenter      
    elseif source(1)-num<=0&&db(1)==0;
        grid(1)=grid(1)-source(1)+8;
        source(1)=5; % 
    elseif source(1)-num>0&&db(1)~=0;
        grid(1)=grid(1)-num+ceil(max(slipb(end,:))/0.2)+7;
        source(1)=source(1)-num+4;
    elseif source(1)-num>0&&db(1)==0;
        grid(1)=grid(1)-num+7;
        source(1)=source(1)-num+4;
    end
elseif num1<=2&&num2>2;                             % Zero appears at the bottom of the dip direction and the fault moves upward
    % 底部空白较多：裁底部；若顶部有滑移贴边，则同时向上扩展并移动 source(1)。
    num=sum(db(:)==0);                              
    if grid(1)-num>=source(1)&&da(1)==0;             % Determine if the fault plane has reached the surface and make sure there are two rows below the source
       grid(1)=grid(1)-num+4;     
    elseif grid(1)-num>=source(1)&&da(1)~=0;       
        grid(1)=grid(1)-num+ceil(max(slipb(1,:))/0.2)+8;
        source(1)=source(1)+ceil(max(slipb(1,:))/0.2)+4;                       
    elseif grid(1)-num<source(1)&&da(1)~=0
        grid(1)=source(1)+ceil(max(slipb(1,:))/0.2)+8;
        source(1)=grid(1)-5;
    elseif grid(1)-num<source(1)&&da(1)==0 ;         
        grid(1)=source(1)+7;
        source(1)=grid(1)-2;
    end
elseif num1<=2&&num2<=2; % 
    % 上、下两边都没有足够空白：断层在触边方向扩展，并保留固定缓冲行。
    if da(1)~=0&&db(1)~=0;
        grid(1)=grid(1)+ceil(max(slipb(1,:))/0.2)+ceil(max(slipb(end,:))/0.2)+7;
        source(1)=source(1)+ceil(max(slipb(1,:))/0.2)+4; 
    elseif da(1)~=0&&db(1)==0;
        grid(1)=grid(1)+ceil(max(slipb(1,:))/0.2)+7;
        source(1)=source(1)+ceil(max(slipb(1,:))/0.2)+4;
    elseif da(1)==0&&db(1)~=0;
        grid(1)=grid(1)+ceil(max(slipb(end,:))/0.2)+7;
        source(1)=source(1)+4;
    elseif da(1)==0&&db(1)==0;
        grid(1)=grid(1)+7;
        source(1)=source(1)+4;
    end
elseif num1>2&&num2>2;
    % 两边都有较宽空白：主要执行裁剪，同时确保震源仍落在网格内部并留缓冲。
    numa=sum(da(:)==0);numb=sum(db(:)==0);
    if numa>=source(1)&&source(1)<grid(1)-numb
        grid(1)=grid(1)-source(1)-numb+7;
        source(1)=4;
    elseif grid(1)-numb<=source(1)
        grid(1)=source(1)-numa+7;
        source(1)=grid(1)-4;
    elseif grid(1)-numb>source(1)
        grid(1)=grid(1)-numa-numb+9;
        source(1)=source(1)-numa+4;
    end
end

% The adjustment of teh fault plane in the strike direction
% 第二部分处理走向（矩阵列）方向，逻辑与倾向方向相同。
% sa 从左侧累计，sb 从右侧累计；num3、num4 是两侧连续空白列数。
sa=cumsum(sum(slipb,1));                            
sb=cumsum(fliplr(sum(slipb,1)));                    % Flip it 90 degrees left and right, and sum it up
num3=sum(sa(:)==0);num4=sum(sb(:)==0);

% Zero appears at the head end of the strike direction, and the fault plane moves along the strike direction
if num3>2&&num4<=2;                                 
    % 左侧空白较多、右侧可能触边：裁左侧，并按右边缘滑移决定向右扩展量。
    num=sum(sa(:)==0); 
    if source(2)-num<=0&&sb(1)~=0;
        grid(2)=grid(2)-source(2)+ceil(max(slipb(:,end))/0.2)+8;
        source(2)=5;
    elseif source(2)-num<=0&&sb(1)==0;
        grid(2)=grid(2)-source(2)+7;
        source(2)=5;
    elseif source(2)-num>0&&sb(1)~=0;
        grid(2)=grid(2)-num+ceil(max(slipb(:,end))/0.2)+7;
        source(2)=source(2)-num+4;
    elseif source(2)-num>0&&sb(1)==0;
        grid(2)=grid(2)-num+7;
        source(2)=source(2)-num+4; 
    end
    
% Zero appears at the end of the strike direction, and the fault plane moves in the opposite direction
elseif num3<=2&&num4>2;                             
    % 右侧空白较多：裁右侧；若左边缘有滑移，则向左扩展并平移震源列号。
    num=sum(sb(:)==0);
    if grid(2)-num>=source(2)&&sa(1)~=0;
       grid(2)=grid(2)-num+ceil(max(slipb(:,1))/0.2)+7;
       source(2)=source(2)+ceil(max(slipb(:,1))/0.2)+3; 
    elseif grid(2)-num>=source(2)&&sa(1)==0; 
       grid(2)=grid(2)-num+7; 
       source(2)=source(2)+3;   
    elseif grid(2)-num<source(2)&&sa(1)~=0;
        grid(2)=source(2)+ceil(max(slipb(:,1))/0.2)+7;
        source(2)=grid(2)-4;
    elseif grid(2)-num<source(2)&&sa(1)==0;
        grid(2)=source(2)+7;
        source(2)=grid(2)-4;        
    end

elseif num3<=2&&num4<=2; 
    % 左右均触边或仅有少量空白：按两侧边缘滑移幅值扩展走向网格。
    if sa(1)~=0&&sb(1)~=0;
     grid(2)=grid(2)+ceil(max(slipb(:,1))/0.2)+ceil(max(slipb(:,end))/0.2)+7;
     source(2)=source(2)+ceil(max(slipb(:,1))/0.2)+4;
    elseif sa(1)~=0&&sb(1)==0;
     grid(2)=grid(2)+ceil(max(slipb(:,1))/0.2)+7;
     source(2)=source(2)+ceil(max(slipb(:,1))/0.2)+4;
    elseif sa(1)==0&&sb(1)~=0;
     grid(2)=grid(2)+ceil(max(slipb(:,end))/0.2)+7;
     source(2)=source(2)+4;
    elseif sa(1)==0&&sb(1)==0;
     grid(2)=grid(2)+7;
     source(2)=source(2)+4;
    end

elseif num3>2&&num4>2;
    % 左右都有足够空白：裁掉空白列，再把震源换算到新列号。
    numa=sum(sa(:)==0);
    numb=sum(sb(:)==0);
    if source(2)>=numa&&grid(2)-source(2)>=numb 
        grid(2)=grid(2)-numa-numb+8;
        source(2)=source(2)-numa+4;    
    elseif numa>source(2)
        grid(2)=grid(2)-source(2)+9-numb;
        source(2)=5;
    elseif numb>grid(2)-source(2)
        grid(2)=grid(2)+9-numa;
        source(2)=grid(2)-5;
    end
end
end

