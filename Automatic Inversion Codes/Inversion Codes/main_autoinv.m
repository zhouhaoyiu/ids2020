% This is the automatic finite-fault inversion system, including automatic station screening, 
% finite-fault inversion, and updating of the fault plane.
% Data: near-field full waveforms, e.g. strong motion data, high-rate GPS data
% Green function: QSSP (Wang et al, 2017, GJI)
% If any qustion， please email us （zhang-yong@pku.edu.cn;zhengxujun@pku.edu.cn）
% Here we carry an inversion in application to the earthquake (Mw7.7) occurred at 15:15 on 11 March 2011 JST (See Figure 2).
% You should prepare for the "earthquake_Info.txt" in advance
% -------- 07/2020 --------
% MAIN_AUTOINV 是本包的顶层脚本，按三个阶段完成一次自动有限断层反演。
%
% 运行前输入
%   1. folder：事件目录的上级目录。目录中需有一个地震信息 txt 和一个事件子目录。
%   2. pathg：QSSP 格林函数数据库目录，需含 GreenInfo.dat、grn_d* 及走时表依赖。
%   3. inEa：earthquake_Info.txt 数值表中的事件行号。
%   4. inMt：选择第一或第二节面。
%
% 主要磁盘输入
%   earthquake_Info.txt：发震时刻及每个事件 11 个数值字段。
%   事件子目录/*.dat：前若干文件为三分量记录，最后一个 dat 文件作为台站信息表读取。
%   earth_offshore.mat：变量 earth，每行 [深度,Vp,Vs,密度]。
%
% 主要结果
%   工作区留下 grid、source、locasub、locadep、slip、substf、syn、res、Mw、Moment 等变量；
%   Slip_Info 另写 Dname_Rupture_Info.txt；脚本还生成台站、波形、滑移和 STF 图。
%
% 三阶段数据流
%   第一阶段：事件与台站文件 -> 重采样观测 ob -> 初始断层 grid/source -> 格林函数 g
%   第二阶段：台站分组试反演 -> 波形拟合筛选 -> 正式 IDS 反演 -> 迭代调整断层边界
%   第三阶段：裁剪最终断层 -> 绘图 -> 计算 Mw/Moment -> 写破裂模型文本。
%
% 当前示例边界
%   路径使用反斜杠，按 Windows 编写；数据读取逻辑假定一次只处理一个事件子目录、
%   一个地震信息 txt。depth 在第 36 行读成所有事件的整列，因此多事件配置需特别核对。
% 第一次 clear all 清空工作区变量、函数持久状态；下一行重复执行，效果相同。
clear all
clear all
% 关闭已有图窗，避免新结果与旧图混在一起。
close all

%% Input variables 
% Path of the database. Note that it is not 'F:\Automatic Inversion Codes\strong_motion\20110311151500'
% folder 指向包含事件子目录的根目录；pathg 指向格林函数数据库，不是某个事件目录。
folder='E:\Automatic Inversion Codes\strong_motion\'; 
pathg='E:\green_Tohoku\green_func\'; % Path of Green's functions.
inMt=1; % Selection of focal machanism (1= Nodle plane I; 2=Nodal plane II)
inEa=1; % List of earthqaukes. Here we just give an example
% --------------------- END -----------------------

%% ----------------------Starting-------------------------------
% ====== Part 1: Data Processing and Green's functions ======
% Information of the earthquakes, you can invert any numer of eathquakes. 
% Here we just give an example 
% -------------------------------------------------------------------------
file=dir([folder,'*.txt']); %（see earthqauke_Info.txt）
% 当前写法直接使用 file.name，要求 folder 中匹配到唯一 txt 文件。
fid=fopen([folder,'\',file.name],'r');
% 第一行跳过；第二行作为发震时刻文本，并在末尾加 JST。
fgets(fid);
for k=1:1;epitime=fgets(fid);end;
Event=strcat(epitime,32,'JST'); clear epitime % The origion time，here use the JST
for k=1:1;ch=fgets(fid);end;clear ch; % The k should be changed with the number of earthquake
ear_all=fscanf(fid,'%f');
% 每个事件固定 11 个数，先按 11 行重排再转置，使每行对应一个事件。
ear_all=reshape(ear_all,11,size(ear_all,1)/11);ear_all=ear_all';
fclose(fid);
epi=[ear_all(inEa,2) ear_all(inEa,1)];  % Epicenter, [latitude,longitude]
% depth 当前保存所有事件的深度列，而不是只取 depth(inEa)；示例只有一个事件时二者相同。
depth=ear_all(:,3); % Depth                                              
mag=ear_all(inEa,10); %  Moment magnitude
if inMt==1
    % 第一节面取第 4:6 列 [strike,dip,rake]。
    fault=ear_all(inEa,4:6); % Nodle plane I
elseif inMt==2
    % 第二节面取第 7:9 列。
    fault=ear_all(inEa,7:9); % Nodle plane II
end
clear inMt;

% ------ Information of stations ------
% Point to path of stations
file_sm=dir(folder); % The strong motion records which download from NIED
% 只保留子目录名并排除 . 和 ..；示例假定最终仅剩一个事件目录。
isub=[file_sm(:).isdir];
nameFolds={file_sm(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
Dname=char(nameFolds);clear nameFolders;
% datafolder 是事件记录目录；多个 nameFolds 会形成字符矩阵，后续路径拼接不再适用。
datafolder=strcat(folder,Dname,'\');                                    
% 事件目录中最后一个 dat 被当作台站信息文件，其余 dat 被当作波形文件。
file=dir([datafolder,'*.dat']);
fid=fopen([datafolder,'\',file(end).name],'r');
StaInfo=textscan(fid,'%s %f %f %f %f','HeaderLines',13);
fclose(fid);
c1=cell2mat(StaInfo(2));c2=cell2mat(StaInfo(3));c3=cell2mat(StaInfo(4));
c3=double(c3);
locat=[c1 c2 c3];clear c1 c2 c3; % the station locations and their relative times to the origin time
% ob 初始为 1*3*S，后面赋入更长记录时 MATLAB 自动扩展第一维；第三维对应台站。
ob=zeros(1,3,length(file)-1); % the observed seismograms
loca=locat(:,1:2); % Locations of station 
dt=locat(:,3); % Relative times
% dt 随后直接用于数组下标，因而必须是非负整数采样点，而不是任意秒数。
clear locat;

% ------ Set the sampling rate ----
rate0=4; % the sampling rate of the original data (sps)
if mag<7.5
    % 较小事件保留 2 sps，较大事件降到 1 sps，以控制长时窗反演规模。
    srate=2; % sampling rate we wanted (sps)
elseif mag>=7.5
    srate=1; % sampling rate we wanted (sps)
end
dtime=0; % the time difference between the data beginning and earthquake origin time

% ------- Load the observed data -------
for i=1:length(file)-1; % the last file is the station information, not include here
    % data0 应为 Ni*3，每列一个分量；resample 同时改变采样点数并作抗混叠滤波。
    data0=load([datafolder,file(i).name]); % load the data
    data0=resample(data0,srate,rate0); % change the sampling data
    ob(1+dt(i):size(data0,1)+dt(i),:,i)=data0; % preserve the data in the data matrix 
end
clear data0;
ob(1:round(dtime*srate),:,:)=[]; % Fix the data beginnings at the earthquake origin time
% --- get the station names ---
mm=' ';
for i=1:length(file)-1
    % 删除 .dat 扩展名，把文件基本名逐行写入字符矩阵 mm。
    mm0=file(i).name;
    mm(i,1:length(mm0)-4)=mm0(1:end-4);
end
clear mm0 dtime file;

% ------ The length and width of the fault are estmated by the empirical ----- 
% scaling law of Wells and Coppersmith (1994,BSSA）
% Type: overthrust/strike-slip/nornal fault
% earthquake_Info.txt 第 11 列给断层类型；下式把 Mw 换算为初始破裂长度 RLD 和宽度 RW(km)。
if ear_all(inEa,11)==1 % The overthrust fault
    RLD=10^(-2.42+0.58*mag); 
    RW=10^(-1.61+0.41*mag); 
elseif ear_all(inEa,11)==2 % The strike-slip fault 
    RLD=10^(-2.57+0.62*mag); 
    RW=10^(-0.76+0.27*mag); 
elseif ear_all(inEa,11)==3 % The nornal fault 
    RLD=10^(-1.88+0.5*mag); 
    RW=10^(-1.14+0.35*mag); 
end
clear earth_all;

% ----- Subfault sizes in [dip,strike] directions and frequency bands
% based on earthquake magnitudesSub-fault size 
% --- Subfault sizes, Mw<6.5 & 6.5≤Mw<7.5 & 7.5≤Mw<8.5 & Mw>8.5 ---
% 不同震级采用 2、5、10、20 km 网格；grid 被构造成奇数，使 source 位于初始网格中央。
if mag<6.5
    gridsize=[2,2]; 
    grid=[ceil(ceil(RW/2)/2)*2+1,ceil(ceil(RLD/2)/2)*2+1];
    source=[ceil(ceil(RW/2)/2)+1,ceil(ceil(RLD/2)/2)+1];
elseif mag>=6.5&&mag<7.5
    gridsize=[5,5];
    grid=[ceil(ceil(RW/5)/2)*2+1,ceil(ceil(RLD/5)/2)*2+1];
    source=[ceil(ceil(RW/5)/2)+1,ceil(ceil(RLD/5)/2)+1];
elseif mag>=7.5&&mag<8.5
    gridsize=[10,10];
    grid=[ceil(ceil(RW/10)/2)*2+1,ceil(ceil(RLD/10)/2)*2+1];
    source=[ceil(ceil(RW/10)/2)+1,ceil(ceil(RLD/10)/2)+1];
else
    gridsize=[20,20];
    % 本分支 grid(2) 使用 ceil(ceil(RLD/2)/20)，与其他分支的 ceil(RLD/gridsize) 结构不同；
    % 这里记录原式，不判断它是否为作者有意设置。
    grid=[ceil(ceil(RW/20)/2)*2+1,ceil(ceil(RLD/2)/20)*2+1];
    source=[ceil(ceil(RW/20)/2)+1,ceil(ceil(RLD/20)/2)+1];
end
% --- Frequency bands ---
% 震级越大，采用的最高频率越低；flh 将作为三阶 Butterworth 带通范围。
if mag>=8.5
    flh=[0.02,0.05]; 
elseif mag>=7.5&&mag(inEa)<8.5
    % mag 已是标量；当 inEa>1 时 mag(inEa) 会越界。当前单事件示例 inEa=1 可执行。
    flh=[0.02,0.1];
elseif mag>=6.5&&mag<7.5
    flh=[0.02,0.2];
elseif mag<6.5
    flh=[0.02,0.5];
end

% ----- Remove the stations whose epicentral distances over three times the fault length -----
% 距离门槛为最近台站距离加三倍经验破裂长度；三分量记录、坐标和名称同步删除。
da=da_zh(loca,epi,1);dis_epi_sta=da(:,1);clear da;
[dex_max]=find(dis_epi_sta>3*RLD+min(dis_epi_sta));
ob(:,:,dex_max)=[];loca(dex_max,:)=[];mm(dex_max,:)=[];
clear dex_max dis_epi_sta RLD RW;

% ------ Median screening -----
% 对每个台站先把三分量相加，再取时间绝对峰值；峰值高于中位数 3 倍或低于 0.1 倍的台站被删。
max_obs=max(abs(sum(ob,2)));max_obs=max_obs(:);
med_obs=median(max_obs);
[dex_max]=find(max_obs>3*med_obs|max_obs<0.1*med_obs);
ob(:,:,dex_max)=[];loca(dex_max,:)=[];mm(dex_max,:)=[];
clear max_obs med_obs dex_max;

% ------ Plot the locations of stations and epicenter -----
% 青色三角是初筛后的台站，白色五角星是震中。
figure
plot(loca(:,2),loca(:,1),'k^','markerfacecolor','c');hold on
plot(epi(2),epi(1),'kp','markerfacecolor','w','markersize',12);

% ----- Make the depth of designed source fix to the real source depth ----
% get_subloca 使用 [strike,dip] 顺序，所以先把外部 [dip,strike] 的 grid 和 source 对调。
[~,dep]=get_subloca(fault,grid([2,1]),gridsize,[1,source([2,1])],loca,epi);
% MATLAB 列优先线性下标把 [source_row,source_col] 变为一维子断层号。
index=(source(2)-1)*grid(1)+source(1); % The location of epicenter in subfaults
% 整体平移所有子断层深度，使 source 对应深度等于输入 depth。
dep=dep+(depth-dep(index));
if dep(1)<0
    % 若顶部若干完整子断层行高出地表，把 source 向上移动相同行数。
    num=sum(dep(:)<0)/grid(2);
    source(1)=source(1)-num;
end
index=(source(2)-1)*grid(1)+source(1);

% ------ Below is to get the Green's functions -----
% Get the green's functions， N/E/U is the positive
% ids_getg 输出 L*3*(S*nsub)，第三维依次列出每个子断层对应的全部台站。
[g,locasub,locadep]=ids_getg(pathg,fault,grid,gridsize,source,loca,epi,srate,inEa,index,depth); clear index;

% Make sure the green's fucntions have the same length as the observations
g(size(ob,1)+1:end,:,:)=[];
nsta=size(ob,3); % number of components
% 此处 nsta 实际是台站数 S；nsub 是子断层总数。
nsub=prod(grid); % number of sub-faults
% 分别取 E/N/U，并把组合维恢复为 S*nsub，得到 N*S*nsub。
g1=g(:,1,:);g1=g1(:,:);g1=reshape(g1,[size(g,1),nsta,nsub]);
g2=g(:,2,:);g2=g2(:,:);g2=reshape(g2,[size(g,1),nsta,nsub]);
g3=g(:,3,:);g3=g3(:,:);g3=reshape(g3,[size(g,1),nsta,nsub]);
g=cat(2,g1,g2,g3);
% 沿第二维拼接后 g 为 N*(3S)*nsub，通道顺序 [全部E,全部N,全部U]。
clear g1 g2 g3;
% QSSP 输出在此按速度处理，一次时间积分后变为位移格林函数。
g=cumsum(g)/srate; % Integrate velocities to displacements
 
% ---- Change the order of components, the observation and green's function must have the same orders -----
% 原 ob 为 N*3*S；permute 后成 N*S*3，再压平为 N*(3S)，顺序与 g 相同。
ob=permute(ob,[1,3,2]);ob=ob(:,:); % Make ob=[EW..EW..EW/NS..NS..NS/UD..UD..UD]
% Remove the baseline drift. this just a simple operation
for i=1:size(ob,2)
    % 每个通道独立校正基线；len1 取整条记录长度。
    acc=cor_baseline(ob(:,i),size(ob,1));
    ob(:,i)=acc;
end
clear acc;
% Integrate observed accelerations to velocities
% 注释写“to velocities”，但代码连续积分两次：加速度 -> 速度 -> 位移，与位移格林函数匹配。
ob=cumsum(ob)/srate;
ob=cumsum(ob)/srate;
% 通道数从 S 更新为 3S；台站名重复三份以对应三个分量块。
nsta=nsta*3;mm=[mm;mm;mm];
% Plot the observed waveforms
[fit]=plotobsyn(epi,loca,ob(:,:),ob(:,:),mm,srate);

% ------ Below is to correct the green's functions by shear modulus ------
% Get the velocity model
% the crust is the velocity model used to calculate the green's fucntions,it has four colomns ... 
% Depth,Vp,Vs,Density, here we use only the Vs and Density to get the shear modulus  Velocity Structure
load earth_offshore.mat;
% Calculate the shear modulus
% histc 给每个子断层深度分配速度模型层号 y。
[~,y]=histc(locadep(:),earth(:,1));                         
% 剪切模量 μ=ρVs^2；密度和速度单位通过 1e9 转为 Pa。
miu=earth(:,4).*earth(:,3).^2.*1e9;
miu=miu(y); % miu is the shear modulus of all sub-faults
for i=1:nsub
    % for sub-fault, correct its amplitude to make sure that the smoothing
    % is made for the slips, not for the moments 
    g(:,:,i)=g(:,:,i)*miu(i)/3e10;
end
% 这一比例把不同深度的矩响应调整为更接近滑移参数化的统一参考模量 3e10 Pa。
clear x y;

% ----- Sceening the station based on epicentral distance, inter-station spacing, and inter-azimuth spacing -----
% 此处 maxdis=400 km，而 minsta、invsta、azista 都传 0；主要执行重复站和距离筛选。
[ob,loca,g,mm,~,~,~]=dist_azim_del(ob,loca,g,mm,epi,400,0,0,0);

% Select the length of waveforms, because some waveforms are too long
% 每个通道先除以自身绝对峰值，再沿通道求和并累计；达到总累计量 90%% 的点作为统一截断点。
oba=abs(ob);[idx]=max(oba);idx=repmat(idx,size(oba,1),1);oba=oba./idx;
cumene=(sum(oba,2));
cumene=cumsum(cumene);
cumene=cumene/cumene(end); 
[~,ndex]=min(abs(cumene-0.9)); % select a threshold 
clear oba idx cumene;
ob=ob(1:ndex,:);
% 观测和格林函数使用同一 ndex 裁剪，保持卷积时间长度一致。
g=g(1:ndex,:,:);
[fit]=plotobsyn(epi,loca,ob(:,:),ob(:,:),mm,srate);
clear ndex fit;
close Figure 2 Figure 3;

%% --------------------------------------------------------------
% ==== Part 2: Below is to invert the strong motion records automatically ====
% ------------------------------------------------------------------------
% ----- One: Station Screening based on waveform fits (misfit<=0.6) -----
if size(loca,1)>=12 % if the number of stations over 12, we carry out the station grouping
    % 至少 12 台时先用每组 5 次 IDS 快速反演评估波形拟合，不作为最终 30 次迭代结果。
    iter=5;
    [obn,gn,locan,mmn,stfb]=divi_inver(epi,loca,ob,g,mm,srate,grid,gridsize,source,miu,flh,iter,mag);
    % ----- Rmove the duplicates（locan/mmn/obn/gn） -----
    % 配对分组允许同一台站多次出现；这里按坐标 unique，并同步保留三分量通道。
    m=size(locan,1);
    [locan,adex]=unique(locan,'rows');  % Removed the duplicated stations
    mmn=mmn(adex,:); % The station names 
    vdex=[adex,m+adex,2*m+adex];
    vdex=vdex(:);
    obn=obn(:,vdex); % The rearranged obseved data
    gn=gn(:,vdex,:); % The rearranged green's functions  
    clear vdex adex m stfb;
else
    % 台站不足 12 时跳过分组试反演，全部初筛数据直接进入正式反演。
    obn=ob;gn=g;locan=loca;mmn=mm;
end
[~,n]=size(obn);
clear loca ob g mm

% If the number of stations over 40, we will reselect the stations with a larger inter-station spacing
% dist 是所有台站相对第一个台站的距离。azim_del 在这里实际按该距离的一维间隔筛选，
% 并非按地理方位角筛选；先试 15 km，若剩余不足 20 台则改用 5 km。
da=da_zh(locan,locan(1,:),1);dist=da(:,1);
if size(locan,1)>=40
    [del_index]=azim_del(dist,15);del_ob=[del_index,n/3+del_index,2*n/3+del_index];
    if size(locan,1)-size(del_index,1)<20
        [del_index]=azim_del(dist,5);del_ob=[del_index,n/3+del_index,2*n/3+del_index];
    end
    locan(del_index,:)=[];mmn(del_index,:)=[];gn(:,del_ob,:)=[];obn(:,del_ob)=[];
end
ob=obn;g=gn;loca=locan;mm=mmn; clear obn gn locan del_ob del_index dist;
% 以上赋值把筛选结果设为正式反演的新输入；ob/g/loca/mm 仍保持逐台站对应。

% ------ After finishing the station screnning, we make a final inversion ------
% Updating the fault plane automatically
for j=1:10 % the maximum updating number
    % 每轮都在当前断层网格上进行最多 30 次 IDS 迭代。
    nsta=size(ob,2);
    iter=30; % the maximum iteration number
    % ----- Two: Inversion -----
    [~,substf,substfa,syn,obr,res,~,loca,time_2]=ids_data(ob,g,loca,srate,grid,gridsize,source,nsta,miu,flh,iter,mag);
    clear time_2;
    
    % Plot the source time function per updating
    % substf 沿子断层求和得到本轮总源时间函数；这里只显示前 200 个采样点。
    figure;
    stf=sum(substf,2);plot(stf(1:200)); % Source time function
    
    % slip distribution per updating
    % 时间积分形式的离散求和按原程序除以 3e16 和子断层面积，得到每块滑移并恢复二维网格。
    slip=sum(substf(:,:))/3e16/prod(gridsize);slip=reshape(slip,grid);
    figure; 
    xz=(1:grid(2))-source(2);xz=repmat(xz*gridsize(2),[grid(1),1]);  
    yz=(1:grid(1))';yz=repmat(yz*gridsize(1)-gridsize(1)/2,[1,grid(2)]);
    pcolor(xz,yz,slip);
    axis ij equal tight 
    xlabel('Distance along strike (km)','FontSize',14); ylabel('Distance down dip (km)','FontSize',14);
    h=colorbar; set(get(h,'Title'),'string','m','FontSize',14);
    
    % Calculate the moment magnitude 
    % 每块 M0=slip*miu*area；汇总后由 m2m 的 M0>=20 分支换算 Mw。
    disp(['Mw ',num2str(m2m(sum(slip(:).*miu(:)*prod(gridsize)*1e6)))])   
    
 % ------ Three: The automatic adjustment of the fault sizes -----
    % 以本轮最大滑移归一化，低于或等于 20%% 的单元置 0，只用主要破裂区判断边界。
    slipb=slip/max(slip(:));
    slipb(slipb<=0.2)=0; % Set the slip values below 20% of the maximum slip value to zero
    da=cumsum(sum(slipb,2));          num1=sum(da(:)==0);    % The number of top zeros                        
    db=cumsum(flipud(sum(slipb,2)));  num2=sum(db(:)==0);    % Flip it up and down, the number of zeros at the bottom
    sa=cumsum(sum(slipb,1));          num3=sum(sa(:)==0);    % The number of zeros on the left-hand side                 
    sb=cumsum(fliplr(sum(slipb,1)));  num4=sum(sb(:)==0);    % Flip left and right, the number of zeros on the right
    clear da db sa sb
    % 若顶部已经接近地表或至少有一行空白，并且其余三边各有空白，当前断层已包住主要滑移区。
    if locadep(1)<gridsize(1)||num1>=1
        if num2>=1&&num3>=1&&num4>=1 
            break
        end
    end
    if j==10 % the maximum updating number
        % 达到十轮上限时保留本轮反演结果，不再建立新断层网格。
        break
    end
    % Adjust the fault sizes, [grid/spurce]
    % subfaults_making 根据四边空白和触边幅值给出下一轮 grida/sourcea。
    [grida,sourcea]=subfaults_making(slipb,grid,source); 
    
    % Determine the depth of the sub-faults, and if it reaches above the surface, delete the sub-faults
    [~,dep]=get_subloca(fault,grida([2,1]),gridsize,[1,sourcea([2,1])],loca,epi);
    % 先按新网格计算深度，再整体平移，使新 source 仍落在观测震源深度 depth(inEa)。
    index=(sourcea(2)-1)*grida(1)+sourcea(1);
    dep=dep+(depth(inEa)-dep(index));
    
    % Renew the [grid/source]
    if dep(1)<0
        % 删除所有完整位于地表以上的倾向行，并同步上移震源行号。
        num=sum(dep(:)<0)/grida(2);
        sourcea(1)=sourcea(1)-num;
        grida(1)=grida(1)-num;
    end
    grid=grida;source=sourcea;
    % 从这一行开始，grid/source 已更新；旧 g 与旧子断层数不再匹配，必须重新读取。
    clear gn substf 
    index=(source(2)-1)*grid(1)+source(1);
    clear grida sourcea;
    
% ----- Reprepare the green's functions -----
   % 对新网格重复第一阶段的流程：查库 -> 裁到观测长度 -> 重排 E/N/U -> 计算 μ -> 位移积分。
   [g,locasub,locadep]=ids_getg(pathg,fault,grid,gridsize,source,loca,epi,srate,inEa,index,depth);
   g(size(ob,1)+1:end,:,:)=[];
   nsub=prod(grid); % number of sub-faults
   nst=size(loca,1);
   g1=g(:,1,:);g1=g1(:,:);g1=reshape(g1,[size(g,1),nst,nsub]);
   g2=g(:,2,:);g2=g2(:,:);g2=reshape(g2,[size(g,1),nst,nsub]);
   g3=g(:,3,:);g3=g3(:,:);g3=reshape(g3,[size(g,1),nst,nsub]);
   g=cat(2,g1,g2,g3);clear g1 g2 g3;
   % 新 locadep 重新映射到速度层，每个新子断层都得到对应剪切模量。
   [x,y]=histc(locadep(:),earth(:,1));
   miu=earth(:,4).*earth(:,3).^2.*1e9;
   miu=miu(y); % miu is the shear modulus of all sub-faults
   for i=1:nsub
       g(:,:,i)=g(:,:,i)*miu(i)/3e10;
   end
   g=cumsum(g)/srate; % integrate velocities to displacements
   clear x y
end
clear num num1 num2 num3 num4 slipb;
% ----- Last adjustment of fault plane -----
% substfa(:,:,end) 是最后一轮 IDS 的累计矩率结果；先换算一维滑移供 adjust_grid 判断。
slip=sum(substfa(:,:,end))/3e16/prod(gridsize);
% 最终只裁掉 20%% 阈值外的多余边缘，不再读取新格林函数；dex 给出旧网格中保留的子断层号。
[dex,grid,source]=adjust_grid(slip,grid,source,0.2); % 
% 所有带子断层维的数据都用同一个 dex 裁剪，确保新 grid 与结果一致。
substf=substfa(:,dex,end);locadep=locadep(dex,:);
locasub=locasub(dex,:);miu=miu(dex);g=g(:,:,dex);
slip=sum(substf)/3e16/prod(gridsize);
slip=reshape(slip,grid);
% ids_data 返回的 obr 是本轮滤波并恢复振幅的观测，作为最终绘图观测。
ob=obr;clear obr;

%% --------------------------------------------------------------
% ==============Part 3: Plot the results ===============
% ------------------------------------------------------------------------
close Figure 2 Figure 4
% ----- Plot the slip distribution -----
% locasub 的经纬度分别恢复成 grid 矩阵，pcolor 在地图坐标上显示滑移。
figure; 
xx=reshape(locasub(:,2),grid);yy=reshape(locasub(:,1),grid);
plot(loca(:,2),loca(:,1),'k^','markerfacecolor','c');hold on
pcolor(xx,yy,slip); %shading interp
set(gca,'dataaspectratio',[1,cosd(epi(1)),1])
axis tight
plot(epi(2),epi(1),'kp','markerfacecolor','w','markersize',12);
colorbar

% ----- Comparisons of the observed (black lines) and synthetic (red lines) strong-motion seismograms -----
% plotobsyn 返回逐分量拟合度 fit，并同时生成按震中距排序的波形比较图。
[fit]=plotobsyn(epi,loca,ob,syn,mmn,srate); clear fit;

% ----- Plot the subfaults source time functions -----
% 总 STF 积分并归一化；累计矩达到 99%% 的时刻 durtime 作为绘图和输出截止点。
stf=sum(substf,2);
stfa=cumsum(stf)/srate;
stfa=stfa/stfa(end); 
[~,durtime]=min(abs(stfa-0.99)); % Just a reference
[~,dey]=substfs_plot(epi,grid,gridsize,source,substf(1:durtime,:),srate);
clear stfa dex dey;

% ----- Save the slip models as a txt -----
% 总地震矩为各子断层 μ*面积*滑移之和，面积由 km^2 乘 1e6 转为 m^2。
Mw=m2m(sum(slip(:).*miu(:)*prod(gridsize)*1e6)); % Magnitude
Moment=sum(slip(:).*miu(:)*prod(gridsize)*1e6); % Moment
[slip_model]=Slip_Info(epi,depth,grid,gridsize,source,locasub,locadep,fault,slip,substf(1:durtime,:),Moment,Mw,Dname,earth,Event,datafolder);
% slip_model 同时留在工作区，文本文件写到 datafolder。
% The save path can be changed based on your needs.
% The all source parameters can be seen in the Workspace of Matlab
% ------------- END ---------------

