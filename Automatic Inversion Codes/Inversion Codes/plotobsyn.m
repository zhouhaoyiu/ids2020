function [fit]=plotobsyn(epi,locan,obr,syn,mmn,srate,rc)
%==========================================================================
% epi: The epicenter [lat,lon]
% locan:The locations of stations [lat,lon]
% obr: The observed data，[leng,nsta*3],arrangement：EW..EW..EW  NS..NS..NS  UD..UD..UD
% syn: The synthetic data，[leng,nsta*3],arrangement：EW..EW..EW  NS..NS..NS  UD..UD..UD
% rc:  Required ranks and columns to plot the waveforms [row,column]
% mmn: The name of stations
% srate: Sampling rate /sps
% =========================================================================
% 中文说明
% 输入
%   obr、syn 均为 N*(3S)，列顺序是 [全部EW,全部NS,全部UD]；locan、mmn 各有 S 行。
%   rc=[每列图中的台站数,图的列数]。省略 rc 时，rowcolumn 根据 S 自动排版。
%
% 输出
%   fit: 1*3*S，每个台站三个分量的观测—合成归一化内积。
%   函数还会新建图窗：黑线是观测，红线是合成；每个分量显示前会按两者共同峰值归一化。
%
% 数据流
%   按震中距给台站排序 -> 把 [EW块,NS块,UD块] 改成每台站连续三个分量
%   -> 形成 N*3*S 三维数组 -> 计算拟合度 -> 分栏、归一化、加垂直偏移后绘图。
% rc
% rat
if nargin==6
    % 六参数调用没有 rc，由台站数自动选择不超过 8 行的布局。
    [rc]=rowcolumn(locan);   
end
    
nsta=size(obr,2);
[~,n]=size(obr);
% Rearrangement of the waveforms;EW..EW..EW NS..NS..NS UD..UD..UD ==> EW/NS/UD...EW/NS/UD
% Arrangement based on the epicenter distance
da=da_zh(locan,epi,1);
da1=da(:,1);
% ndex 按震中距从近到远排列，名称和三分量波形都跟随它重排。
[~,ndex]=sortrows(da1);
mmn=mmn(ndex,:);
EW=obr(:,1:n/3);NS=obr(:,n/3+1:2*n/3);UD=obr(:,2*n/3+1:n);
EW=EW(:,ndex);NS=NS(:,ndex);UD=UD(:,ndex);
ob=[EW;NS;UD];clear EW NS UD

% 先纵向堆成 3N*S，再 reshape 为 N*(3S)，得到每个台站的 EW/NS/UD 连续列。
[m,n]=size(ob);
ob=reshape(ob,m/3,n*3);
[~,n]=size(syn);
EW=syn(:,1:n/3);NS=syn(:,n/3+1:2*n/3);UD=syn(:,2*n/3+1:n);
EW=EW(:,ndex);NS=NS(:,ndex);UD=UD(:,ndex);
syn=[EW;NS;UD];clear EW NS UD
[m,n]=size(syn);
syn=reshape(syn,m/3,n*3);

% 最终转为 N*3*S：第二维是分量，第三维是按距离排序的台站。
ob=reshape(ob,size(ob,1),3,nsta/3);
syn=reshape(syn,size(syn,1),3,nsta/3);
t=0:1/srate:size(ob,1)/srate;                                          
% 上式会多生成末端一个点，裁掉后时间轴长度与 N 完全一致。
t(size(ob,1)+1:end)=[];
                                                         
rat=0.05;                                                            
wid0=1/(rc(2)*(1+rat)+rat);                                               % The width of each subgraph
wid=wid0*rat;                                                             % The interval of each subgraph

% rat=0.05 控制子图列间留白；wid0 是每列宽度，wid 是实际间隔。
rat_sta=1.0;
nsta=size(ob,3);
rc1=rc(1);
ylim=[rc1*6+(rc1-1)*rat_sta,0];
yzero=zeros(3,rc1);
% 为每一行台站的三个分量生成基线位置，分量间隔 2，台站间另加 rat_sta。
for j=1:rc1
    yzero(:,j)=ylim(1)-[1;3;5]-(j-1)*rat_sta-6*j;
end
yzero=yzero(:);
    
% fit 在原始振幅三维数组上计算；随后才为显示创建白色图窗。
fit=gfit1(ob,syn);                                                        % The correlation coefficient between the observed and synthetic waveforms
figure('color','w')                                                      
% 
nt=min(t);mt=max(t);
for i=1:rc(2)                                                             % The number of columns
    % 每次循环建立一个 axes，容纳最多 rc1 个台站。
    axes('position',[(i-1)*wid0+wid*i,0.08,wid0,0.9]);                    % Coordinate setting

    dex0=(i-1)*rc1+1;
    if dex0>nsta;
        % 所有台站已画完时关闭多余坐标轴并结束。
        axis off
        return;
    end
    
    ob0=ob(:,:,(i-1)*rc1+1:min(i*rc1,size(ob,3)));                        
    syn0=syn(:,:,(i-1)*rc1+1:min(i*rc1,size(ob,3)));                   
    mm0=mmn((i-1)*rc1+1:min(i*rc1,size(ob,3)),:);                  
     
    % 展平成 N*(3*本栏台站数)，让 vecnor 对每个分量列独立归一化。
    [ob0,syn0,maxwave]=vecnor(ob0(:,:),syn0(:,:),1);                      % Normalization of the waveforms, maxwave:amplitude peak
    % 每列加对应 yzero，使多条波形在同一坐标轴中上下错开而不重叠。
    for j=1:size(ob0,2)
        ob0(:,j)=ob0(:,j)+yzero(j);
        syn0(:,j)=syn0(:,j)+yzero(j);
    end    
    
     plot(t,ob0,'k','linewidth',1);%1.0
     hold on
     plot(t,syn0,'r','linewidth',1);%0.5
     fsize=7;                                                             % Font size
     for j=1:size(mm0,1)
        text(nt,yzero((j-1)*3+1)+0.1,mm0(j,:),'fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','bottom','FontName','Times New Roman');
        if i==1&&j==1
            text(mt-30/srate,yzero((j-1)*3+1)+0.5,'EW','fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','top','FontName','Times New Roman');
            text(mt-30/srate,yzero((j-1)*3+1)-1.5,'NS','fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','top','FontName','Times New Roman');
            text(mt-30/srate,yzero((j-1)*3+1)-3.5,'UD','fontsize',fsize,'HorizontalAlignment','left','VerticalAlignment','top','FontName','Times New Roman');
        end
    end
       
    axis tight
    set(gca,'ycolor','w','ytick',[],'tickdir','out','xminortick','on',...
        'ylim',[-7,max(yzero)]+1,'ticklength',[0.02,0.02]);%,'grid','on')
    box off
    xlabel('Time(s)')
end

%==========================================================================
function [mout1,mout2,mm]=vecnor(min1,min2,rat)
%% % Normalization of the waveforms（ob/max(ob,syn);syn/max(ob,syn)）
% VECNOR 对每个通道使用观测与合成两者中较大的绝对峰值作共同归一化。
if nargin==2;
    rat=1;
end
sm=size(min1);
mm1=max(abs(min1));mm2=max(abs(min2));
% mm(i) 是第 i 个通道共同的振幅标尺；两条波形除以同一值才能保留相对振幅。
mm=max(mm1,mm2);
for i=1:sm(2)
    if mm(i)==0
        continue;
    end
    min1(:,i)=min1(:,i)/mm(i)*rat;%-i*2;
    min2(:,i)=min2(:,i)/mm(i)*rat;%-i*2;    
end
mout1=min1;mout2=min2;

function [rc]=rowcolumn(locan)
%% Determining the arrangement of the waveforms
% ROWCOLUMN 令图形尽量接近方形，同时限制每列最多放 8 个台站。
n=size(locan,1);
dex=sqrt(n);
if dex<=8
    rca=ceil(dex);
    rcb=ceil(n/rca);
    rc=[rca rcb];
else
    rca=8;
    rcb=ceil(n/rca);
    rc=[rca rcb];
end

