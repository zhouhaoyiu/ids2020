function [g_new]=getg_wang_sub0(fileg,leng,dist,t,t0,srate,dist0,rate)
%==========================================================================
%
%--------------------------------------------------------------------------
%  Input
%   fileg: the file of the green's function
%    leng: the length of the green's function
%    dist: epicentral distance, a vector
%       t: [t_p,t_s], the arrival times of the epicentral distances
%      t0: [t0_p,t0_s], the arrival times of the epicentral distances
%          contained in the database of green's functions
%   srate: sampling rate
%   dist0: the epicentral distances contained in the database of green's functions
%  Outpit 
%    g_new: obtained green's functions
%
% 中文说明
%   fileg 是某一震源深度的二进制格林函数文件；每个距离记录含 3 个头字段和
%   10 个基本机制分量，每个分量原始长度为 leng+2。
%   dist 为目标距离向量，dist0 为数据库已有距离网格，二者单位必须一致。
%   t 为目标距离的 [P,S] 到时，t0 为数据库距离网格的 [P,S] 到时，单位为秒。
%   srate 是数据库采样率，rate 是需要输出的采样率。
%
% 输出
%   g_new: L*10*length(dist)。第三维顺序与输入 dist 完全一致；L 会随最大时移增加。
%
% 数据流
%   每个目标距离 -> 找最近数据库距离 -> 相同距离只读一次二进制记录
%   -> 去掉头字段和每列首尾保护点 -> 必要时重采样到 rate
%   -> 按文件内起始时间对齐震源时刻 -> 分别调整 P 前段和 P-S 段长度
%   -> 恢复输入 dist 的重复项和原顺序。
%--------------------------------------------------------------------------
%    Zhang Yong, 2012-02-06 14:04, GFZ, Potsdam
%==========================================================================

% ensure the epicentral distance is monotonely inceased: it is not necessary
% because we 'sort' and 'unique' the index of the dist next
%[dist,dex]=sort(dist); 

% find the distance index
% distnum(i) 是离 dist(i) 最近的数据库距离下标；这里采用最近邻查库，不作距离振幅插值。
distnum=zeros(size(dist));
for i=1:length(dist)
    [~,nd]=min(abs(dist(i)-dist0));
    distnum(i)=nd;
end

% get the arrival times for the green's function in database
% 把完整数据库到时表裁成每个目标距离实际选中的最近邻到时。
t0=t0(distnum,:);

% make unique and sort for the distance
% udnum 是需要真正读取的唯一距离记录；adex 将唯一记录映射回原 dist 顺序。
[udnum,~,adex]=unique(distnum);


% the shift bytes for each epicentral distance
% 单精度 float 占 4 字节：10*(leng+2) 个波形值，再加 3 个头字段。
shift0=(leng+2)*10*4+4*3;
% the length of data read each time
leng0=(leng+2)*10+3;

ts=zeros(length(udnum),1);
% Suitable for lower versions of Matlab
% 先按 rate/srate 估计重采样后的长度；每个唯一距离保存 10 个基本格林分量。
g_temp=zeros(round(leng*rate/srate),10,length(udnum));
% Suitable for higher versions of Matlab
% g_temp=zeros(ceil(leng*rate/srate),10,length(udnum));


fid=fopen(fileg,'r');
for i=1:length(udnum)
    % move the file pointer and read data
    % 文件指针从开头移动到第 udnum(i) 条定长记录，再读这一整条记录。
    fseek(fid,shift0*(udnum(i)-1),-1);
    z=fread(fid,leng0,'float');
    
    % the shift value
    % 第二个头字段是该记录相对震源时刻的起始偏移，单位按秒处理。
    ts(i)=z(2);
    %z(1:3)=[];
    z=z(4:end);
    
    % the data
    % 剩余数据恢复成 (leng+2)*10，并删除每个分量首尾各一个保护点。
    z=reshape(z,[leng+2,10]);
    %z([1,end],:)=[];
    z=z(2:end-1,:);
    
    % make resample for the green's functions and the sampling rate is
    % 'rate' now, not 'srate' again
    if rate~=srate
        % MATLAB resample(z,rate,srate) 把采样率从 srate 变到 rate，并带抗混叠滤波。
        z=resample(z,rate,srate);
    end
    
    g_temp(:,:,i)=z; % 1
%    g_temp(:,:,adex(i))=z; % 1+2+3
end
fclose(fid);

% adex 把只读一次的唯一距离结果复制、重排成与原始 dist 一一对应的第三维。
g_temp=g_temp(:,:,adex); %2
ts=ts(adex);
%--------------------------------------------------------------------------

% shift the data by fixing their beginnings at the origin time
% 文件头的秒偏移乘输出采样率，得到整数采样点偏移。
ts=round(ts*rate);
% 预留最大正偏移所需长度；未写入位置保持 0。
g=zeros(max(ts)+size(g_temp,1),size(g_temp,2),size(g_temp,3));
for i=1:size(g_temp,3)
    % get the synthetic seismograms by shift with ts
    if ts(i)<0
        % 记录早于震源时刻时，裁掉震源前的 -ts(i) 个点。
        zg=g_temp(1-ts(i):end,:,i);
        g(1:size(zg,1),:,i)=zg;
    else
        % 记录晚于震源时刻时，在前面留 ts(i) 个零点。
%         zg=[zeros(ts(i),size(g_temp,2));g_temp(:,:,i)];
%         g(1:size(zg,1),:,i)=zg;
        g(ts(i)+1:ts(i)+size(g_temp,1),:,i)=g_temp(:,:,i);
    end
end

clear g_temp;

% make shifting for the P wave and S wave, respectively
% 所有到时从秒换成输出采样率下的整数点号。
t=round(t*rate);
t0=round(t0*rate);
tps=round(t(:,1)-t0(:,1));
tss=round(t(:,2)-t0(:,2));
% tps、tss 是目标距离到时相对最近数据库距离到时的 P、S 点数改变量。
% size(g,1)+max(tss)
% size(g)
g_new=zeros(size(g,1)+max(tss),size(g,2),size(g,3));

for i=1:size(g,3)
    % 把序列分为 P 到时以前、P 到 S 之间、S 以后三段。
    gp=g(1:t(i,1),:,i); % t0?
    gps=g(t(i,1)+1:t(i,2),:,i); % t0?
    %[size(gp) tps(i)]
    gp_new=t_shift(gp,tps(i)); % for waves before P arrival times
    % 第一段长度改变 tps，使 P 到时对齐；第二段再改变 tss-tps，使 S 到时也对齐。
    gps_new=t_shift(gps,tss(i)-tps(i)); % for waves between P ans S arrival times
    
    %g_new(1:size(g,1)+tss(i),:,i)=[gp_new;gps_new;g(t(i,2)+1:end,:,i)];
    g_new(1:size(gp_new,1),:,i)=gp_new;
    g_new(1+size(gp_new,1):size(gp_new,1)+size(gps_new,1),:,i)=gps_new;
    
%     size(g(t(i,2)+1:end,:,i))
%     size(g_new(1+size(gp_new,1)+size(gps_new,1):size(g,1)+tss(i),:,i))
    
    g00=g(t(i,2)+1:end,:,i);
    % S 到时后的第三段不再重采样，直接接到两个已调整分段之后。
    g_new(size(gp_new,1)+size(gps_new,1)+(1:size(g00,1)),:,i)=g00;
%    g_new(1+size(gp_new,1)+size(gps_new,1):size(g,1)+tss(i),:,i)=g(t(i,2)+1:end,:,i);
end
return
%-----------------------------------------------------------------------old
