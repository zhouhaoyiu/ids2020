function [slip,substf,substfa,syn,ob,res,Mw,loca,time_2]=ids_data(ob,g,loca,srate,grid,gridsize,source,nsta,miu,flh,iter,Mw)
% This function is to invert the waveforms using the IDS method
% IDS_DATA 完成一次 IDS 反演前后的数据准备，并把最终矩速率换算为滑移量。
%
% 输入
%   ob      : N*nsta 观测波形，第一维为时间，第二维为波形通道。
%   g       : N*nsta*nsub 格林函数，第三维对应子断层。
%   loca    : 台站位置；本函数不修改它，只原样返回。
%   srate   : 采样率，单位 samples/s。
%   grid    : [ndip,nstrike] 子断层网格数，nsub=prod(grid)。
%   gridsize: [dip_km,strike_km] 单个子断层尺寸。
%   source  : [dip_index,strike_index] 震源所在子断层。
%   nsta    : 波形通道数，应等于 size(ob,2) 和 size(g,2)。
%   miu     : 各子断层剪切模量，元素数应为 nsub。
%   flh     : Butterworth 滤波截止频率；标量表示低通，二元素向量表示带通，单位 Hz。
%   iter    : IDS 最大迭代次数。
%   Mw      : 参考矩震级，传给 IDS_2020 作矩震级约束，并原样作为本函数输出之一。
%
% 输出
%   slip   : ndip*nstrike 最终滑移分布，单位换算沿用原程序公式。
%   substf : fftlen*nsub 最终子断层矩率函数。
%   substfa: fftlen*nsub*K，各次迭代累计得到的子断层矩率函数；提前收敛时 K<iter。
%   syn,ob : 恢复原振幅后的合成波形和滤波观测波形，尺寸均为 N*nsta。
%   res     : IDS 每次迭代的相对残差。
%   Mw,loca: 输入值原样传出。
%   time_2 : IDS 记录的各迭代用时。
%
% 数据流
%   原波形和格林函数 -> 同一滤波器 -> 按每个观测通道能量归一化 -> IDS_2020
%   -> 恢复通道振幅 -> 对最终矩率沿时间求和 -> 除以面积及单位系数得到 slip。

%% filter the data and green's fucntions 
% low and high cut-off frequencis, in Hz, can be ajusted
% butter 的归一化截止频率为 flh/(srate/2)，所以写成 flh*2/srate。
[bb,aa]=butter(3,flh*2/srate); % butterworth filter
% 两者使用相同因果滤波器，保持观测与格林函数的频带一致。
ob=filter(bb,aa,ob); % filter the data
g=filter(bb,aa,g); % filter the green's functions

%% normalizing the data to make sure they have the equal weights in inversion
% Wob(i) 是第 i 个观测通道的 L2 范数，即 sqrt(sum(ob(:,i).^2))。
Wob=sqrt(sum(ob.^2));
for i=1:nsta
    % 同时除观测和对应格林函数，既平衡各通道权重，又不改变二者的相对振幅关系。
    ob(:,i)=ob(:,i)./Wob(i);
    g(:,i,:)=g(:,i,:)./Wob(i);
end

%% some inversion parameters for IDS inversion, usually needn't change
iter_ids=iter; % maximum iteration number
vlim=[1,6];   % minimum and maximum rupture velocity, in km/s
% 允许的源时间长度先取整条观测记录长度。
lenrup=size(ob,1); % rupture length

%% the following five parameters is meaningless 
fh=0.5; % highest frequency of filter for temporal smooth
weit=ones(nsta,1);% all stations are equally weighted in inversions.
% 先把 N 提升到最近的 2 次幂，再乘 srate；这是原程序传给 IDS 的频域长度。
fftlen=round(2^nextpow2(size(ob,1))*srate); % length of fast fourier transforms
% 缩放时使用全部观测元素，线性索引范围为 1 到 numel(ob)。
obdex=1:numel(ob); % the index of observed data used in scalings

%% carry out the inversion with IDS method
% IDS_2020 返回归一化数据尺度下的矩率函数、合成波形和迭代历史。
[substf,syn,res,~,substfa,~,time_2,fg]=IDS_2020(ob,g,iter_ids,vlim,lenrup,srate,grid,...
    gridsize,source,fh,weit,fftlen,obdex,Mw,miu);

%% recover the amplitude of the observed and synthetic data
for i=1:nsta
    % 乘回各通道原来的 L2 范数，撤销反演前的归一化。
    ob(:,i)=ob(:,i)*Wob(i);
    syn(:,i)=syn(:,i)*Wob(i);
end

%% plot the result: fault slip distribution
% slip: finial fault slip
% 取最后一次累计矩率，沿时间求和；3e16 和子断层面积完成原程序采用的单位换算。
slip=sum(substfa(:,:,end))/3e16/prod(gridsize);
% 一维子断层顺序恢复为 [ndip,nstrike] 断层面。
slip=reshape(slip,grid);
% xz=[1:grid(2)]-source(2);xz=repmat(xz*gridsize(2),[grid(1),1]);
% yz=[1:grid(1)]';yz=repmat(yz*gridsize(1)-gridsize(1)/2,[1,grid(2)]);
% Mw=num2str(m2m(sum(slip(:).*miu(:)*prod(gridsize)*1e6)));
% disp(['Mw ',num2str(m2m(sum(slip(:).*miu(:)*prod(gridsize)*1e6)))])

end

