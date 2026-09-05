function [dex,grid_new,source]=adjust_grid(slip,grid,source,factor)
% This function is to adjust the fault plane
% ADJUST_GRID 根据已有滑移分布裁掉断层面四周的低滑移空白区。
%
% 输入
%   slip  : prod(grid) 个滑移量，可为向量；MATLAB 按列顺序重排为 grid。
%   grid  : [ndip,nstrike]，原断层面沿倾向、走向的子断层数。
%   source: [row,column]，震源子断层在原网格中的下标；函数会更新它。
%   factor: 相对阈值，slip/max(slip) <= factor 的位置被视为空白。
%
% 输出
%   dex     : 新网格各子断层在原 slip 向量中的线性下标，顺序与 MATLAB (: ) 一致。
%   grid_new: 裁剪后的 [ndip,nstrike]。
%   source  : 震源在新网格中的 [row,column]。
%
% 为什么保留边界
%   程序尽量裁掉连续空白行、列，同时在有效滑移或震源周围保留约两格缓冲，
%   避免下一轮反演把断层面裁得贴住震源。输入 slip 应至少含一个正值；
%   全零输入会在除以 max(slip(:)) 时产生 NaN。
% 把一维滑移量恢复成断层面；MATLAB 按列填充，即倾向索引变化最快。
slipa=reshape(slip,grid);
% 用全局最大滑移归一化到相对幅值，再释放临时矩阵。
slipb=slipa/max(slip(:));clear slipa
% 小于或等于 factor 的位置归零，后面只关心边缘连续零区。
slipb(slipb<=factor)=0; 
% 正向/反向累计行和、列和；累计值仍为 0 的个数就是对应边缘的连续空白数。
da=cumsum(sum(slipb,2));          num1=sum(da(:)==0);     % The number of top zeros                            
db=cumsum(flipud(sum(slipb,2)));  num2=sum(db(:)==0);     % Flip it up and down, the number of zeros at the bottom
sa=cumsum(sum(slipb,1));          num3=sum(sa(:)==0);     % The number of zeros on the left-hand side                   
sb=cumsum(fliplr(sum(slipb,1)));  num4=sum(sb(:)==0);     % Flip left and right, the number of zeros on the right
% The top
if num1<=2||source(1)<=3
    % 顶部空白不多或震源已靠近顶部：从原第 1 行开始，不移动震源行号。
    dex_top=1;source(1)=source(1);
elseif num1<source(1)&&num1>2 
    % 可裁掉顶部空白，并把旧 source 行号换算到新网格。
    dex_top=num1;source(1)=source(1)-num1+1;
elseif num1>=source(1)&&num1>2&&source(1)>3
    % 空白区跨过震源时，只裁到震源上方两行，确保震源仍在网格内。
    dex_top=source(1)-2;source(1)=source(1)-1;
end

% The bottom
if num2<=2||grid(1)-source(1)<3
    % 底部空白不多或震源靠近底部：保留到原网格最后一行。
    dex_bottom=grid(1);source(1)=source(1);
elseif num2>=grid(1)-source(1) 
    dex_bottom=source(1)+2;source(1)=source(1);
elseif num2>2&&num2<grid(1)-source(1)+1 
    dex_bottom=grid(1)-num2+2;source(1)=source(1);
end

% The left
if num3<=2||source(2)<=3
    % 左侧逻辑与顶部相同，但处理的是走向列和 source(2)。
    dex_left=1;source(2)=source(2);
elseif num3<source(2)&&num3>2
    dex_left=num3;source(2)=source(2)-num3+1;
elseif num3>=source(2)&&num3>2&&source(2)>3
    dex_left=source(2)-2;source(2)=source(2)-1;
end
    
% The right
if num4<=2||grid(2)-source(2)<3
    % 右侧逻辑与底部相同，确定最后保留的原网格列号。
    dex_right=grid(2);source(2)=source(2);
elseif num4>=grid(2)-source(2) 
    dex_right=source(2)+2;source(2)=source(2);
elseif num4>2&&num4<grid(2)-source(2)+1 
    dex_right=grid(2)-num4+2;source(2)=source(2);
end

% x、y 是保留的原列号和原行号；meshgrid 枚举新矩形中的全部子断层。
x=dex_left:dex_right;y=dex_top:dex_bottom;
[X,Y] = meshgrid(x,y);grid_new=size(X);
dex=[X(:) Y(:)];
% 二维 [列,行] 下标换成原 grid 的列优先线性下标，调用者可直接 slip(dex)。
dex=(dex(:,1)-1)*grid(1)+dex(:,2);
    
