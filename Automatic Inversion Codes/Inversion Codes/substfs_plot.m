function [dex,dey]=substfs_plot(epi,grida,gridsize,source,substf,srate)
% This is a function to plot the sub-fault source time functions 
% SUBSTFS_PLOT 把每个子断层的源时间函数画在断层面对应网格内。
%
% 输入
%   epi     : 震中 [lat,lon]；程序只使用纬度 epi(1) 修正绘图纵横比。
%   grida   : [ndip,nstrike] 子断层网格数。
%   gridsize: [dip_km,strike_km] 子断层尺寸。
%   source  : [dip_index,strike_index] 震源网格下标。
%   substf  : Nt*prod(grida)，每列是一条子断层源时间函数。
%   srate   : 采样率，用于给每条曲线末尾留出约 1 s 的横向空白。
%
% 输出
%   dex: 走向坐标显示范围 [xmin,xmax]，单位 km。
%   dey: 倾向坐标显示范围 [0,ndip*dip_km]，单位 km。
%   函数还会新建图窗并直接绘图。
%
% 绘图逻辑
%   每条 STF 先除以全局最大值，使振幅落在相对尺度；随后将时间映射到网格宽度，
%   将振幅映射为单元格内向上的填充高度。这样可同时看到破裂起始、持续时间和相对幅值。
% 
% close all
% 末尾补 3 个零，使填充多边形回到底边，避免曲线末端悬空。
substf=[substf;zeros(3,prod(grida))];   
% 所有子断层共用同一幅值尺度；若 substf 全零，这一步会产生 NaN。
substf=substf/max(substf(:));
% 第二维、第三维分别对应倾向行和走向列。
b=reshape(substf,size(substf,1),grida(1),grida(2));                       
if mod(grida(2),2)==0
    % 偶数走向网格没有正中的列，坐标在震源附近作半格偏移。
    xz=(0:grida(2))-source(2)+1;xz=repmat(xz*gridsize(2),[grida(1),1]);xz=xz-gridsize(1)/2;
else
    % 奇数走向网格以震源列为中心，再整体减去半个走向网格宽。
    xz=(1:grida(2)+1)-source(2);xz=repmat(xz*gridsize(2),[grida(1),1])-gridsize(2)/2;
end
dex=[min(xz(:)),max(xz(:))];                                         
dey=[0,grida(1)*gridsize(1)];                                          

% ------ Plot the epi----
figure
% 星形标记位于走向 0 km、震源所在倾向网格中心。
plot(0,source(1)*gridsize(1)-gridsize(1)/2,'p','LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor','y','MarkerSize',20);% 画震源位置
set(gca,'linewidth',2);                                              
axis([dex,0,grida(1)*gridsize(1)]);                               
set(gca,'xtick',dex(1):gridsize(1):dex(2),'FontSize',10)              
set(gca,'ytick',dey(1):gridsize(1):dey(2),'FontSize',10);                
set(gca,'YDir','reverse')                                               
set(gca,'dataaspectratio',[1,cosd(epi(1)),1])                         
grid on
set(gca,'GridLineStyle','-');                                        
xlabel('Distance along strike (km)','LineWidth',10,'FontSize',10);
ylabel('Distance down dip (km)','LineWidth',10,'FontSize',10);
hold on

% ----- Plot the source time functions ----
% st 收集全部填充多边形坐标；当前函数不输出它，只用于循环内累积。
st=[];
dp=gridsize(1);                                                        
str=gridsize(2);                                                       
for j=1:grida(2)                                                       
    % stf 为第 j 个走向列上的全部倾向子断层，尺寸 Nt*ndip。
    stf=b(:,:,j);
    for i=1:grida(1)                                                     
        y=i*dp-stf(:,i)*dp;
        % 把 Nt 个采样点压入一个走向网格宽，并额外按 srate 留出 1 s 对应的空白比例。
        nt=str/(size(substf,1)+srate*1);                               
        x=min(xz(:))+str*(j-1):nt:min(xz(:))+str+str*(j-1);             
        x=x';
        x(size(y,1)+1:end,:)=[];
        sa=[x,y];sa=[sa;sa(1,:)];
        st=[st;sa];
        plot(x,y);
        fill(x,y,'c');
        hold on 
    end
end
end

