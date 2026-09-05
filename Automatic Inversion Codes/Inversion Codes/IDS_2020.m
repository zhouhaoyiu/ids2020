function [SUBSTF,syn,resa,syna,substfa,smnum,time_2,fg]=IDS_2020(ob,g,iter,vlim,lenrup,srate,...
    grid,gridsize,source,fh,weit,fftlen,obdex,Mw,miu)   
%==========================================================================
% This is a function of IDS method to solve the rupture process   // 破裂过程
%--------------------------------------------------------------------------
% Input
%       ob: observed seismograms,[leng,nsta]  ob是观测波形，[leng,nsta]是观测波形的长度和台站数
%        g: green's functions with given mechanisn,[leng,nsta,nsub] g是给定机制的格林函数，[leng,nsta,nsub]是格林函数的长度、台站数和子断层数
%     iter: maximum iteration number   iter是最大迭代次数 这个迭代次数是指在IDS方法中进行的迭代次数，用于逐步优化破裂过程的解
%     vlim: [minimum, maximum] of rupture velocity vlim是破裂速度的最小值和最大值，用于约束破裂过程中的速度范围
%   lenrup: permitted source-time-function duration, in samples. // 允许的源时间函数持续长度，按采样点计
%    srate: sampling rate of the data and green's functions // srate是数据和格林函数的采样率，用于确定时间分辨率
%     grid: numbers of sub-faults along [dip, strike] directions  // grid是沿[倾角，走向]方向的子断层数量，用于定义破裂过程的空间分辨率
% gridsize: sub-fault size in [dip, strike] directions, unit is in km  // gridsize是沿[倾角，走向]方向的子断层大小，单位为公里，用于定义每个子断层的空间范围
%   source: the index of the sub-fault hypocenter located on, [dip, strike] // source是子断层震源所在的索引，[倾角，走向]，用于确定破裂过程的起始位置
%       fh: highest frequency of filter for temporal smooth // fh是用于时间平滑的滤波器的最高频率，用于控制破裂过程的时间分辨率
%     weit: weight of each seismograms  // weit是每个地震波形的权重，用于在反演过程中对不同波形进行加权
%   fftlen: length of fast fourier transforms  // fftlen是快速傅里叶变换的长度，用于计算频域表示
%    obdex: the index of observed data used in scalings  // obdex是用于缩放的观测数据索引，用于选择用于计算缩放因子的观测数据子集
% 
% NOTICE: for use of post-earthquake studies, we needn't care about the input parameters
% after 'fh'. // 原注释表示常规震后应用通常沿用 fh 之后参数的默认设置；这些参数在代码中仍会参与计算，不能理解为没有影响
%
% Output
%  SUBSTF: output results. has the size of [fftlen, nsub], each column is 
%          the momen rate of one sub-fault // SUBSTF是输出结果，大小为[fftlen, nsub]，每一列是一个子断层的矩速率 example: SUBSTF(:,1) 是第一个子断层的矩速率，SUBSTF(:,2) 是第二个子断层的矩速率，以此类推 fftlen表示快速傅里叶变换的长度，nsub表示子断层的数量
%     syn: synthetic waves, has the same size as 'ob' // syn是合成波形，大小与'ob'相同 example: syn(:,1) 是第一个台站的合成波形，syn(:,2) 是第二个台站的合成波形，以此类推
%    resa: relative misfit curves of the iterations // resa是迭代的相对拟合曲线 example: resa(1) 是第一次迭代的相对拟合值，resa(2) 是第二次迭代的相对拟合值，以此类推
%    syna: all synthetic waves obtained in each iteration // syna是每次迭代中获得的所有合成波形 example: syna(:,:,1) 是第一次迭代的合成波形，syna(:,:,2) 是第二次迭代的合成波形，以此类推
% substfa: all cumulative source-time functions obtained in each iteration. // substfa(:,:,k) 是第 k 次迭代结束后的 SUBSTF，不是 syn
%   smnum: temporal rows by iterations. // smnum(i,k) 是第 k 次迭代中第 i 个时间采样行被空间平滑的次数
%  time_2: elapsed seconds recorded for completed non-converged iterations. // 各轮耗时；触发提前收敛的那一轮不会写入
%       fg: FFT of the input Green's functions, size [fftlen,nsta,nsub]. // 输入格林函数的频谱
%
% 整体数据流
%   ob(N*nsta) 与 g(N*nsta*nsub)
%   -> 由观测累计能量估计最长有效源时长 leng_source
%   -> 由断层距离和破裂速度生成允许时间窗 wid(fftlen*nsub)
%   -> 对残差波形逐子断层反卷积、叠加，得到本轮增量 tri_sub
%   -> 灵敏度加权、空间平滑、时间窗裁剪、滑移块筛选
%   -> 合成新增波形 ob0，并累加到 syn 和 SUBSTF
%   -> 用 obori-syn 形成下一轮残差，直到改善量小于 1e-3 或达到 iter。
%--------------------------------------------------------------------------
%             Zhang Yong, Peking University, 2014-05-07 
%==========================================================================



nsta=size(g,2); % the number of stations  // 这里实际是波形通道数；三分量数据通常为台站数的 3 倍
nsub=size(g,3); % the number of sub-faults  // 从格林函数g中获取子断层数量，g的第三维表示子断层数量 
lensub=size(ob,1); % the length of the observed waveforms // 从观测波形ob中获取波形长度，ob的第一维表示波形长度

%% Roughly estimate the maximum source duration time: // 粗略估计最大震源持续时间
mampli=0.95; % 累计归一化能量阈值为 95%，不是最大振幅的 95%

% encuv(ob,3) 先把 3 变为半窗长 round(3/2)=2，再对每个时间点前后共 5 点求平方和。
% 例如单列 ob=[1;2;3]，首尾补 2 个零后，各点局部能量为
% [0^2+0^2+1^2+2^2+3^2; 0^2+1^2+2^2+3^2+0^2; 1^2+2^2+3^2+0^2+0^2]
% = [14;14;14]。参数 3 控制窗口宽度，代码中的能量幂次始终是 2。
ene=encuv(ob,3); 
% cumsum 默认沿第一维累加，把局部能量变成每个通道随时间增长的累计能量。
ene=cumsum(ene);
% max_ene 是每个通道末段的最大累计能量；静默通道最大值为 0，其倒数被安全改为 0。
max_ene=max(ene);max_ene_d1=1./max_ene;max_ene_d1(max_ene==0)=0;
% sparse(i,i,max_ene_d1) 建立对角矩阵；右乘等价于每一列除以自己的最大累计能量。
ene=ene*sparse(1:size(ene,2),1:size(ene,2),max_ene_d1);

% 所有通道的归一化累计曲线逐行相加，再用末值归一到 0~1。
cumene=(sum(ene,2));
cumene=cumene/cumene(end); 
% ndex 取最接近 95% 累计能量的采样点，而不是第一个严格超过 95% 的点。
[~,ndex]=min(abs(cumene-mampli)); 
leng_source=ndex;
% 最长有效源时长不能超过实际观测行数 lensub。
leng_source=min(leng_source,lensub);% The maximum source duration time
disp(['Source duration (s): ',num2str(leng_source/srate)])

%% prepare for the limitation of rupture velocity and source duration time
% The possible time which the subfaults start to rupture
% wid(t,j) 为 1 表示第 j 个子断层允许在第 t 个采样点释放矩率，为 0 表示禁止。
wid=substfwid(grid,gridsize,source,vlim,srate,fftlen,lenrup,nsub,leng_source);
% 非零元素数就是本轮反演允许求解的“时间点—子断层”参数总数。
disp(['Number of unknown parameters: ',num2str(length(find(wid>0)))]);

% All waveform components were equally weighted in inversions
% expob 与 ob 同尺寸且全为 1；ob.*expob 当前不改变数值，保留了未来逐样点加权的接口。
expob=ones(size(ob,1),1);
expob=repmat(expob(:),[1,nsta]);  

%% prepare for the deconvolution speed up; cfgDfg02=conj(fg)/max(fg.^2,max(fg.^2)*alpha^2) 
alpha=0.1; % Water level parameter
% fg 保存原格林函数的 fftlen 点复频谱，后续用于正演合成。
fg=fft(g,fftlen);
% predecon 返回加水位稳定后的逆滤波器，而非原始频谱。
[cfgDfg02]=predecon(g,fftlen,nsta,nsub,alpha); % =conj(fg)/max(fg.^2,max(fg.^2)*alpha^2) 

% 同时用 0.1 和 0.6 两个水位反卷积；后面只保留两种水位都给出正值的位置。
levels=[0.1,0.6];
cfgDfg02=repmat(cfgDfg02,[1,1,1,numel(levels)]);

for i=1:numel(levels);
    cfgDfg02(:,:,:,i)=predecon(g,fftlen,nsta,nsub,levels(i)); % Different Water level factors 
end

%% Prepare and start the iteration:
% 反卷积逆滤波器和原始频谱 fg 已准备好，时间域 g 可释放以节省内存。
clear g;
% 以下数组保存累计解和逐轮历史；第三维或第二维的 k 都对应迭代号。
SUBSTF=zeros(fftlen,nsub); % a matrix to memory the subfaults STFs
resa=zeros(iter,1); % relative misfit curves of the iterations
syna=zeros([size(ob),iter]);% all synthetic waves obtained in each iteration
substfa=zeros(fftlen,nsub,iter); % all results obtained in each iteration, substfa(:,:,end) is 'syn'
smnum=zeros(size(ob,1),iter); % 每行对应一个时间采样点，每列对应一次迭代

% 初始合成波形为标量 0；第一次 syn+ob0 时由 MATLAB 标量扩展成观测同尺寸矩阵。
syn=0; % set the initial synthetical waveforms
obori=ob; % the observed data
% res_ori 是原始观测总平方能量，作为各轮相对残差的固定分母。
res_ori=ob(:)'*ob(:); % the energy of the observed waveforms

% SM1 左乘平滑倾向，SM2 右乘平滑走向；D 用于量化二维滑移粗糙度。
SM1=smmat_gfz(grid,[1,1]);SM2=smmat_gfz(grid([2,1]),[1,1])'; % creat a sparse matrix
D=sm_space_inter(grid,1); % a sparse matrix which is used for smoothing the sliprate on the fault

% Iteration
for k=1:iter;
    tic;
    % Stacking, retrieve the ASTFs, only main lobes(maximum area) were chosen 
    % 输入 ob 已是上一轮剩余残差；tri_sub 为本轮每个子断层新提取的主正波包，尺寸 fftlen*nsub。
    [~,tri_sub]=wldecon2016(ob.*expob,cfgDfg02,wid,fftlen,nsub,nsta,weit,fh,srate);

    % Scaling 
    % get_res 分别让每个子断层的候选矩率最佳拟合当前残差，并返回单块残差能量和相关性。
    [res_sub,tri_sub,fit]=get_res(ob,fg,tri_sub,fftlen,nsta,nsub,obori);

    % Regularization
    % res 是“只用该子断层拟合后的误差/原始观测能量”；大于等于 1 的值截为 1。
    res=res_sub(:)'/res_ori;
    res(res>=1)=1;
    vr=1-res;
    % vr 同时要求该子断层能降低残差且与原始观测正相关；平方扩大高灵敏度块的相对权重。
    vr=(vr(:).^1.*fit(:).^1).^2; % sensitivity factor
    
    if max(vr)==0
        vr(:)=0;
    else
        % 除以最大值后 vr 范围为 0~1，最大贡献子断层权重为 1。
        vr=vr/max(vr);
    end
    % 把 1*nsub 权重复制到全部时间点，逐列缩放本轮矩率增量。
    tri_sub=tri_sub.*(ones(fftlen,1)*vr(:)');
    
    %----------------------------------------------------------------------
    % for substf's spatial smoothing    
    if k==1
        % 首轮沿时间求和得到滑移形状，并把其归一化粗糙度保存为后续各轮允许的上限 rough1。
        slip1=reshape(sum(tri_sub,1),grid);
        sm1=D*slip1(:);
        rough1=norm(sm1)./norm(slip1(:));
        if norm(slip1(:))==0;rough1=0;end
    end     
    % sumtri(t) 为该时刻所有子断层矩率之和；等于 0 的时间行无需空间平滑。
    sumtri=sum(tri_sub,2);
    
    %---------------------------------------
    % see if slip is smoothed enough
    if k>1
        % 后续轮先计算当前新增滑移的粗糙度 rough。
        slip=reshape(sum(tri_sub,1),grid);
        % calculate the roughness of the current slip distribution
        sm=D*slip(:);
        rough=norm(sm)./norm(slip(:));
        if norm(slip(:))==0;rough=0;end
        
        while rough>rough1
            % 每执行一轮 while，就对所有非零时间切片各做一次二维邻域平均。
            for i=1:size(ob,1)
                if sumtri(i)==0;
                    continue;
                end
                sliprate=reshape(tri_sub(i,:),grid);
                % 分别沿倾向和走向平滑，再取两者平均，最后恢复成 tri_sub 的一行。
                sliprate1=SM1*sliprate;sliprate2=sliprate*SM2;
                sliprate=(sliprate1+sliprate2)/2;
                tri_sub(i,:)=sliprate(:)';
                smnum(i,k)=smnum(i,k)+1;
            end
            slip=reshape(sum(tri_sub,1),grid);
            sm=D*slip(:);
            rough=norm(sm)./norm(slip(:));
            % 全零滑移的范数分母为 0，程序把粗糙度约定为 0 以退出循环。
            if norm(slip(:))==0;rough=0;end
        end
    end
    %---------------------------------------
    % 空间平滑可能把值扩散到不允许时窗，再乘 wid 把这些值重新清零。
    tri_sub=tri_sub.*wid;
    %----------------------------------------------------------------------
    % Get the synthetic waveforms
    % getsyn 在频域卷积所有子断层贡献；逆变换长度为 fftlen，再裁为观测的 N 行。
    ob0=getsyn(fg,tri_sub);ob0(size(ob,1)+1:end,:)=[];
    % Recaling the syn and substf should be scaled every iteration
    fac0=bestfac(ob(obdex),ob0(obdex));
    % 负比例意味着反相贡献，当前非负矩率约束将它截为 0。
    fac0=max(fac0,0);
    % 这里只先缩放 tri_sub；后面的滑移块检查会重新用它计算 ob0。
    tri_sub=tri_sub*fac0; 
    
    %% Find the substfs which is less valuable for the waveforms, and then, set them to zero. 
    %----------------------------------------------------------------------
    % check the synthetic waves of each slip patch
    % 时间求和后恢复二维滑移，把空间上分离的主要区域拆成 patch(:,:,iii)。
    slip0=reshape(sum(tri_sub),grid);
    [patch,sumslip]=slippatch_new(slip0,nsub);
    numpatch=size(patch,3);
    
    slippatch=zeros(numpatch,1);
    misfitpatch=slippatch;
    for iii=1:numpatch
        % tri_sub0 只保留当前滑移块覆盖的子断层，其他列清零。
        patch0=patch(:,:,iii);
        tri_sub0=tri_sub;
        tri_sub0(:,patch0==0)=0;
        ob0=getsyn(fg,tri_sub0);
        ob0(size(ob,1)+1:end,:)=[];

        % slippatch 记录该块总滑移；misfitpatch 是只用该块时相对当前残差的误差能量。
        slippatch(iii)=sum(patch0(:));
        misfitpatch(iii)=sum((ob(obdex)-ob0(obdex)).^2)./sum(ob(obdex).^2);
    end
    
    % Clear the slip patch which cannot fit the data
    % 相对误差仍大于 0.9999 的块几乎没有降低当前残差，因此整块删除。
    xx=misfitpatch>0.9999;
    patch(:,:,xx)=0;
    slip0=sum(patch,3);
    % 所有保留块都为 0 的子断层，其整条新增矩率列清零。
    tri_sub(:,slip0==0)=0;
    
    % Recaling the syn and substf should be scaled every iteration
    % 删除无效块后重新正演，并再次求非负最佳比例；这次同时缩放 ob0 和 tri_sub。
    ob0=getsyn(fg,tri_sub);ob0(size(ob,1)+1:end,:)=[];
    fac0=bestfac(ob(obdex),ob0(obdex));
    fac0=max(fac0,0);%fac0=min(fac0,1);
    ob0=ob0*fac0;
    tri_sub=tri_sub*fac0;
    %----------------------------------------------------------------------
    % 本轮合成贡献和矩率增量加入此前累计结果。
    syn=syn+ob0;
    SUBSTF=SUBSTF+tri_sub;
    %---------------------------------------------
    % % codes below is to normalize the synthetic waves, it will result in a
    % % more accurate moment magnitude, but it makes the iteration converge slowly
    
    % Recaling the syn and substf should be scaled every iteration
    % faca 用全部累计合成 syn 最佳拟合原始观测 obori；每轮都会整体重标定历史累计结果。
    faca=bestfac(obori(obdex),syn(obdex));
    faca=max(faca,0);
    syn=syn*faca;
    SUBSTF=SUBSTF*faca; 
   
    % 每个子断层的累计矩率沿时间求和并换算为滑移，再与 μ*面积相乘得到总 M0。
    slip=sum(SUBSTF/3e16/prod(gridsize));
    Mw0=m2m(sum(slip(:).*miu(:)*prod(gridsize)*1e6));
    if Mw0-Mw>0.2
        % 若反演 Mw 超过输入 Mw 0.2，仅保留第一轮历史结果中曾为正的“时间—子断层”位置。
        % 在第 1 轮写入 substfa 之前触发此条件时，substfa(:,:,1) 仍是预分配的全零矩阵。
        substf0=substfa(:,:,1);
        substf0((substf0>0))=1;
        SUBSTF=SUBSTF.*substf0;
    end       
    %---------------------------------------------
    % 下一轮输入残差始终由固定原始观测减当前累计合成得到，避免逐轮减法积累舍入误差。
    ob=obori-syn; % after one iteration�� the syn should be removed from the ob
    
    if nargout>3
        % 只有调用者请求第 4 个及之后输出时才保存逐轮历史，以节省赋值工作。
        syna(:,:,k)=syn;
        if nargout>4
            substfa(:,:,k)=SUBSTF;
        end
    end
    
    % 本轮相对残差使用原始观测能量作分母；resa=0 为完全拟合，resa=1 与零合成同级。
    resa(k)=(ob(:)'*ob(:))/res_ori; 
 
    if k>1
        % 改善量小于 0.001 即停止；若残差变差，该差值为负，也会满足当前条件。
        if resa(k-1)-resa(k)<1e-3; % convergence criterion
            % 停止前再用一个全局最小二乘比例缩放累计合成和累计矩率。
            % 此处 fac 没有像前面一样截为非负值，这是原程序行为。
            fac=bestfac(obori(obdex),syn(obdex));           
            syn=syn*fac;
            if nargout>3
                syna=syna*fac;
                if nargout>4
                    substfa=substfa*fac;
                end
            end
            SUBSTF=SUBSTF*fac;
           
            % 删除尚未执行的迭代位置，使输出长度与实际迭代次数一致。
            resa(k+1:end)=[];
            
            if nargout>3
                syna(:,:,k+1:end)=[];
                if nargout>4
                    substfa(:,:,k+1:end)=[];smnum(:,k+1:end)=[];
                end
            end
            
            tt=toc;
            disp(['The ',num2str(k),'-th iteration, ',num2str(tt,'%8.5f'),' second.']);
            % 提前 return 发生在 time_2(k)=tt 之前，所以收敛这一轮耗时只显示、不写入 time_2。
            return
        end
    end
    
    tt=toc;   
    disp(['The ',num2str(k),'-th iteration, ',num2str(tt,'%8.5f'),' second.']);
    % 未触发提前收敛的轮次把耗时写入 time_2(k)。
    time_2(k)=tt;
end

%==========================================================================
function [fg]=predecon(g,fftlen,nsta,nsub,alpha)
%--------------------------------------------------------------------------
% prepare for the water level deconvolution
%
% Input
%       g: green's functions, has the size of [leng,nsta,nsub]
%  fftlen: length of fft
%    nsta: number of seismograms used for deconvolution
%    nsub: number of subfaults
%   plpha: water level
%
% Output
%     fg: water-level inverse spectrum, not the original spectrum
%   cfgFfg02: prepared results, equals conj(fg)/max(fg.^2,max(fg.^2)*alpha^2)
%            also can have this form: conj(fg)/(fg.^2*(1+alpha^2))
%--------------------------------------------------------------------------
% 中文说明
%   输入 g 为 Ng*nsta*nsub，fftlen 是变换长度，alpha 是水位比例。
%   nsta 只为保持接口完整，本局部函数没有直接读取它。
%   输出 fg 与 FFT(g,fftlen) 同尺寸，但每个元素已变成稳定逆滤波器：
%       conj(G) / max(|G|^2, alpha^2*该通道最大|G|^2)。
%   水位下限防止格林函数频谱很小时除法把噪声无限放大。
% 把每个通道的标量水位复制到 fftlen 个频率点。
onefftlen=ones(fftlen,1);
% 原始格林函数沿时间轴转到频域；随后清除本地时间域副本。
fg=fft(g,fftlen);clear g;
for i=1:nsub
    % fg20=|G|^2，尺寸 fftlen*nsta。
    fg20=real(fg(:,:,i)).^2+imag(fg(:,:,i)).^2;
    % 每个通道以自身最大谱功率的 alpha^2 倍作为水位。
    fac=(alpha.^2)*max(fg20);
    fac=onefftlen*fac;
    fg20=max(fg20,fac);
    
    % avoid the holes of green's functions
    % 若某通道整个频谱为 0，水位仍为 0；先记住这些位置，除法后再置 0。
    dexfg20=(fg20==0);
    
    % 共轭谱除以稳定后的功率谱，得到反卷积所需逆滤波器。
    cfg0=conj(fg(:,:,i))./fg20;
    cfg0(dexfg20)=0;
    fg(:,:,i)=cfg0;
end

%==========================================================================
function [stf_sub,stf_pk]=wldecon2016(ob,cfgDfg02,wid,fftlen,nsub,nsta,weit,fh,srate,stf_sub_vec)
%--------------------------------------------------------------------------
% sub-function to retrieve the ASTFs
%
% Input
%   g: has the size of [leng,nsta,nsub], and stf has the size of [leng,nsub,nsta]
%   cfgDfg02: has the size of [fftlen,nsta,nsub]
%   wid: use substf.*wid to constrain the rupture velocity and rupture
%      duration
%   fftlen: length of fft
%   nsub: number of subfaults
%   nsta: number of stations
%   weit: weights of stations while it has the size [nsta,1], else the weights
%       of stations for each subfault while it has the size of [nsta,nsub]
%   fh: highest frequency of filter for temporal smooth (if is given)
%   srate: needed by filter for temporal smooth (if is given)
%
% Output
%   stf_sub: stacked stfs of subfaults
%   stf_pk: only main lobes were chosen from the stf_sub
%  
%--------------------------------------------------------------------------
% 中文说明
%   ob 为 N*nsta 当前残差，cfgDfg02 为 fftlen*nsta*nsub*nlevel 水位逆滤波器，
%   wid 为 fftlen*nsub 允许时间窗。weit 是 nsta*1 通道权重，或 nsta*nsub 的逐块权重。
%   标准调用省略 stf_sub_vec，共 9 个输入；十参数分支在原注释中标为 Does not work。
%
% 输出
%   stf_sub: fftlen*nsub，多水位一致性筛选、时间窗和可选滤波后的叠加表观 STF。
%   stf_pk : fftlen*nsub，从 stf_sub 每列只保留面积最大的正波包。
%
% 计算原因
%   每个台站先分别反卷积，再按其反卷积谱的 L1 大小归一化后叠加，避免大振幅台站独占结果。
%   两个水位都为正的时间点才保留，可减少对单一水位设置敏感的伪峰。

% 当前残差沿时间轴做 fftlen 点 FFT，一次供全部子断层和水位使用。
fob=fft(ob,fftlen);
% 第四维长度是并行比较的水位个数，本程序为 2。
numlevel=size(cfgDfg02,4);

% stfsuba(:,:,num) 保存第 num 个水位得到的全部子断层 STF。
stfsuba=zeros(fftlen,nsub,numlevel);
for num=1:numlevel
    if nargin==9
        % % full waveform inversion
        % fspec(:,i) 将保存第 i 个子断层跨台站加权叠加后的反卷积频谱。
        fspec=zeros(fftlen,nsub);
        %norvec=zeros(nsta,nsub);
        for i=1:nsub
            % 观测频谱逐元素乘该子断层的逆格林函数，得到每个台站的表观源谱。
            fspec0=fob.*cfgDfg02(:,:,i,num);         
            % 沿频率求复谱幅值之和，得到每个台站的 L1 标尺。
            sfabs=sum(sqrt(real(fspec0).^2+imag(fspec0).^2),1)'; % L1 normalize, determine the energy of each station
            norvec=1./sfabs;
            norvec(sfabs==0)=0;
            
            % % add weit, here weit contains only the station weight of spatial
            % % distribution
            if size(weit,2)==1
                % 单列 weit 对所有子断层使用同一组通道权重。
                norvec=norvec.*weit;
            else
                % 多列 weit 允许每个子断层使用不同通道权重。
                norvec=norvec.*weit(:,i);
            end
            %-------------------------------------------------------------- 
            fspec(:,i)=fspec0*norvec;% energy normalization for each station,make all stations have euqal weight.
        end
        % 逆 FFT 回到时间域，并除以通道数得到平均叠加 STF；real 去掉数值圆整虚部。
        stf_sub=real(ifft(fspec))/nsta;  % The formula��6��of Zhang et al.(2014, JGR)
    elseif nargin>9 % Does not work
        % 旧十参数分支保留每个台站和子断层的时域结果，再由 stf_sub_vec 作逐点选择。
        % while the dt_ts_cut is specified, size of [nsta,nsub]
        fspec=zeros(fftlen,nsta,nsub);
        %norvec=zeros(nsta,nsub);
        for i=1:nsub
            fspec0=fob.*cfgDfg02(:,:,i,num);
            sfabs=sum(sqrt(real(fspec0).^2+imag(fspec0).^2),1)'; % L1 normalize
            
            norvec=1./sfabs;
            norvec(sfabs==0)=0;
            fspec(:,:,i)=fspec0*sparse(1:nsta,1:nsta,norvec);
        end
        
        stf_sub=real(ifft(fspec));
        stf_sub=stf_sub.*stf_sub_vec;
        % 沿台站维求和，再恢复为 fftlen*nsub 并取平均。
        stf_sub=sum(stf_sub,2);
        stf_sub=stf_sub(:,:)/nsta;   
    end
    stfsuba(:,:,num)=stf_sub;
end

stf_sub=stfsuba(:,:,1);
% dexGT0 初始全真，逐水位相与后，只在所有水位结果均严格大于 0 的位置保持真。
dexGT0=true(fftlen,nsub);
for num=1:numlevel
    dexGT0=dexGT0&logical(stfsuba(:,:,num)>0); % Use th last water level factor (allph=0.6)
end

% 保留第一个水位的幅值，但用全部水位的正值交集决定保留位置。
stf_sub(~dexGT0)=0;
%--------------------------------------------------------------------
% 先应用破裂速度/持续时间窗，再把残留负值截为 0。
stf_sub=stf_sub.*wid;
stf_sub(stf_sub<0)=0;

% filter smooth
if srate==2
% 仅当采样率恰为 2 sps 时执行零相位 Butterworth 平滑；其他采样率直接跳过。
[bb,aa]=butter(3,fh*2/srate);
stf_sub=filtfilt(bb,aa,stf_sub);
% 滤波可能在允许窗外产生振铃，故再次乘 wid 并截掉负值。
stf_sub=stf_sub.*wid;
stf_sub(stf_sub<0)=0;
end
%----------------------------------
if nargout==2
    % 第二输出只保留每列累计面积最大的正波包。
    stf_pk=peakstf(stf_sub,1); % area maximum
end
%==========================================================================
%==========================================================================
function [res_sub,tri_sub,fit,synsub]=get_res(ob,fg,tri_sub,fftlen,nsta,nsub,obori,obdex)
% GET_RES 逐个子断层缩放候选 STF，并衡量它单独解释当前残差的能力。
%
% 输入
%   ob      : N*nsta 当前迭代残差。
%   fg      : fftlen*nsta*nsub 原格林函数频谱。
%   tri_sub : fftlen*nsub 候选 STF，每列对应一个子断层；函数会就地缩放各列。
%   obori   : N*nsta 原始观测，用于计算相关性 fit。
%   obdex   : 可选缩放下标；原注释标记该八参数分支 does not work。
%
% 输出
%   res_sub : nsub*1，只用每个子断层时的残差平方和。
%   tri_sub : 每列乘过自身非负最佳比例后的 STF。
%   fit     : nsub*1，该子断层合成波形与原始观测的归一化内积，负值最终截为 0。
%   synsub  : 可选，numel(ob)*nsub；每列是一个子断层合成波形按列展开后的结果。
%
% 这些输出随后组成灵敏度权重：残差降低得多且与原观测相关性高的子断层权重更大。
% 所有候选 STF 一次变换到频域，后面逐列与对应格林函数相乘。
f_trisub=fft(tri_sub,fftlen);

% 预分配每个子断层的误差、相关性和可选合成波形。
res_sub=zeros(nsub,1);
mtri=max(tri_sub);
fit=zeros(nsub,1);
if nargout>3
    synsub=zeros(numel(ob),nsub);
end
for i=1:nsub
    if mtri(i)==0
        % 该子断层没有候选矩率时，合成贡献为 0，所以误差就是当前 ob 本身。
        err=ob(:);
    else
        fg0=fg(:,:,i); 
        ftri0=f_trisub(:,i);
        
        fob0=zeros(fftlen,nsta);
        for j=1:nsta
            % 频域卷积：第 i 个子断层 STF 乘其到第 j 个通道的格林函数。
            fob0(:,j)=fg0(:,j).*ftri0;
        end
        
        ob0=real(ifft(fob0));    
        % 逆 FFT 长度为 fftlen，裁掉观测 N 点以后的尾部。
        ob0(size(ob,1)+1:end,:)=[];  
        
        if nargin<8
            % 标准路径使用全部当前残差求最小二乘比例 Aj。
            fac=bestfac(ob(:),ob0(:)); % Aj 
        else
            fac=bestfac(ob(obdex),ob0(obdex));% does not work
        end
        
        fac=max(fac,0);
        % 合成波形和对应 STF 必须乘相同比例，保持正演关系不变。
        ob0=ob0*fac;
        tri_sub(:,i)=tri_sub(:,i)*fac; 
        
        if nargout>3
            % 每个 N*nsta 合成矩阵按 MATLAB 列顺序拉直，写成 synsub 的一列。
            synsub(:,i)=ob0(:);
        end
        
        err=ob(:)-ob0(:);
        % fit 对比固定原始观测 obori，而残差能量对比当前迭代 ob。
        fit(i)=gfit1(obori(:),ob0(:)); % correlation coefficients
    end
    res_sub(i)=err'*err;
end
% 反相关子断层不参与后续正灵敏度加权。
fit(fit<0)=0;
%==========================================================================

%==========================================================================
function syn=getsyn(fg,substf)
% this is a funcion to get the synthetic seismograms for sub-fault source
% time functions
% GETSYN 将所有子断层源时间函数与对应格林函数卷积并相加。
%
% 输入
%   fg    : fftlen*nsta*nsub 格林函数频谱。
%   substf: fftlen*nsub 子断层源时间函数。
%
% 输出
%   syn   : fftlen*nsta 合成波形。调用者按需要再裁到观测长度。
%
% 频域中卷积变为逐元素乘法；同一子断层 STF 要作用于全部 nsta 通道，
% 所以先把该列频谱横向复制，再将所有子断层贡献累加。

% 每个子断层 STF 沿时间轴转到频域。
f_substf=fft(substf);
nsub=size(substf,2);
nsta=size(fg,2);

fsyn=zeros(size(fg,1),nsta);
% msub 用于快速跳过整列和为 0 的子断层；在本流程 substf 非负，因此和为 0 即整列为 0。
msub=sum(substf);
onesta=ones(1,nsta);
for i=1:nsub
    if msub(i)==0
        continue;
    end
    fg0=fg(:,:,i);
    fsub0=f_substf(:,i);
%     for j=1:nsta
%         fsyn(:,j)=fsyn(:,j)+fg0(:,j).*fsub0;
%     end
    % fsub0*onesta 从 fftlen*1 扩展为 fftlen*nsta，与 fg0 尺寸一致。
    fsub0=fsub0*onesta;
    % 当前子断层的全部通道频谱贡献加入总合成谱。
    fsyn=fsyn+fg0.*fsub0;
end
% 逆 FFT 回到时间域，理论结果为实数，real 去除数值误差产生的极小虚部。
syn=real(ifft(fsyn));

%==========================================================================
function fac=bestfac(ob,syn) 
% BESTFAC 求使 ||ob-fac*syn||^2 最小的无截距比例系数。
%
% 输入 ob、syn 尺寸相同。单列时输出标量；多列时沿第一维逐列输出 1*C。
% 公式为 fac=(ob'*syn)/(syn'*syn)。若某列 syn 能量为 0，该列 fac 设为 0。
% 本函数本身允许负比例；调用者是否截为非负取决于各自上下文。
if size(ob,2)==1
    % 单列向量用矩阵内积得到合成能量。
    msyn=syn'*syn;
    if msyn==0
        fac=0;
        return
    else
        % 正规方程的一参数闭式解。
        fac=ob'*syn/msyn;
    end    
else
    % 多列时分别求每列分母和分子，再逐元素相除。
    syn2=sum(syn.*syn);
    obsyn=sum(ob.*syn);
    fac=obsyn./syn2;
    fac(syn2==0)=0;
end
%==========================================================================

%==========================================================================
function wid=substfwid(grid,gridsize,source,vlim,srate,fftlen,lenrup,nsub,lengob)
% SUBSTFWID 为每个子断层建立允许释放矩率的时间窗。
%
% 输入
%   grid/gridsize/source: 断层网格、网格尺寸(km)和震源网格下标。
%   vlim  : 破裂速度设置(km/s)。主流程使用 [vmin,vmax] 两元素形式。
%   srate : 采样率；fftlen 为输出时间行数。
%   lenrup: 两元素 vlim 分支中作为最短允许窗长度，单位按采样点使用。
%   nsub  : 子断层数，应等于 prod(grid)。
%   lengob: 全部子断层统一允许的最晚时间点；主流程传入 95%% 累计能量点 leng_source。
%
% 输出
%   wid: fftlen*nsub 的 0/1 矩阵。第 j 列的 1 区间由震源距和破裂速度上下限决定，
%        反卷积结果与它逐元素相乘后，过早、过晚的子断层矩率被强制清零。
%
% 两元素速度的物理含义
%   earliest = distance/vmax，最快破裂速度给出最早起裂时刻；
%   latest   = distance/vmin，最慢破裂速度给出最晚时刻。
%   代码还保证 latest 至少比 earliest 晚 lenrup 个采样点，然后再由 lengob 统一截断。
% 先计算每个子断层中心到震源子断层中心的断层面距离，输出尺寸为 grid。
subdist=faultdist(grid,gridsize,source); % the distances between each sub-fault to epicenter
if length(vlim)==3 % Does not work
    % 三速度分支原注释标为不可用；它构造“零—正弦渐入—平台—零”的软时间窗。
    tr1=subdist/vlim(3)*srate;tr2=subdist/vlim(2)*srate;tr3=subdist/vlim(1)*srate;
    mintr3=tr3<lenrup*srate+tr2;
    tr3(mintr3)=lenrup*srate+tr2(mintr3);
    tr=round([tr1(:),tr2(:),tr3(:)]);
    % 超过 FFT 长度的界限截在 fftlen，防止下标越界。
    tr(tr>fftlen)=fftlen;
    
    wid=zeros(fftlen,nsub);
    for i=1:nsub
        % constraints of the rupture velocities
        t1=tr(i,1);t2=tr(i,2);t3=tr(i,3);
        if t1==0&&t2==0
            wid(:,i)=[ones(t3-t2,1);zeros(fftlen-t3,1)];
        else
            if t1==t2
                % 渐入区长度为 0 时直接从 1 开始。
                zw=[zeros(t1-1,1);ones(t3-t2+1,1);zeros(fftlen-t3,1)];
            else
                % t1 到 t2 用 sin 从 0 平滑增至 1，t2 到 t3 保持允许。
                zw=[zeros(t1-1,1);sin(0:pi/2/(t2-t1):pi/2)';ones(t3-t2,1);zeros(fftlen-t3,1)];
            end
            if numel(zw)==fftlen
                wid(:,i)=zw;
            else
                wid(:,i)=zw(1:fftlen);
            end %[zeros(t1-1,1);sin(0:pi/2/(t2-t1):pi/2)';ones(t3-t2,1);zeros(fftlen-t3,1)];
        end
    end
    wid(lengob+1:end,:)=0;
elseif length(vlim)==2
    % 主流程分支：vlim(2)=vmax 给最早时刻 tr1，vlim(1)=vmin 给最晚时刻 tr2。
    tr1=subdist/vlim(2)*srate;tr2=subdist/vlim(1)*srate;% the maximun and minmum travel time ��*srate�� from epicenter to each sub-fault
    % 若速度范围给出的窗太短，把最晚时刻延长到 earliest+lenrup。
    mintr2=tr2<lenrup+tr1; 
    tr2(mintr2)=lenrup+tr1(mintr2);  
    
    tr=round([tr1(:),tr2(:)]);
    
    wid=zeros(fftlen,nsub);
    for i=1:nsub
        % constraints of the rupture velocities
        t1=tr(i,1);t2=tr(i,2);
        % MATLAB 时间下标从 1 开始；前 t1+1 点为 0，随后 t2-t1 点允许为 1。
        zw=[zeros(t1+1,1);ones(t2-t1,1);zeros(fftlen-t2-1,1)];  
        if numel(zw)<=fftlen
            wid(1:length(zw),i)=zw;
        else
            wid(:,i)=zw(1:fftlen);
        end%[zeros(t1-1,1);sin(0:pi/2/(t2-t1):pi/2)';ones(t3-t2,1);zeros(fftlen-t3,1)];
    end
    % 不论各子断层自己的速度窗多长，lengob 之后一律禁止释放矩率。
    wid(lengob+1:end,:)=0;
elseif length(vlim)==1 % Does not work
    % 单速度分支原注释标为不可用：起点由 distance/vlim 决定，随后开放 lenrup*srate 点。
    tr=subdist/vlim*srate;
    
    wid=zeros(fftlen,nsub);
    for i=1:nsub
        % constraints of the rupture velocities
        t0=round(tr(i));
        zw=[zeros(t0+1,1);ones(lenrup*srate,1)];
        if numel(zw)<=fftlen
            wid(1:length(zw),i)=zw;
        else
            wid(:,i)=zw(1:fftlen);
        end%[zeros(t1-1,1);sin(0:pi/2/(t2-t1):pi/2)';ones(t3-t2,1);zeros(fftlen-t3,1)];
    end
    wid(lengob+1:end,:)=0;
end
