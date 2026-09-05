function [obs,dt,mco]=shift_obssyn(obs,syn,srate,dtmax)
% SHIFT_OBSSYN 用互相关估计观测波形相对合成波形的整数采样点时移。
%
% 输入
%   obs   : N*C 或 N*S*K 观测波形，第一维是时间。
%   syn   : 与 obs 尺寸一致的合成波形。
%   srate : 采样率，单位 samples/s。
%   dtmax : 允许的最大时移，单位 s；超出范围的估计会被改为 0。
%
% 输出
%   obs: 按估计时移移动后的观测波形，尺寸不变；空出的样点用 0 填充。
%   dt : 每列或每个台站的时移，单位是采样点，不是秒。换算成秒需用 dt/srate。
%   mco: 移动后 obs 与 syn 的归一化内积。
%
% 计算过程
%   1. 对时间轴补到 2 的幂并作 FFT。
%   2. ifft(fft(obs).*conj(fft(syn))) 得到循环互相关，峰值位置给出整数时移。
%   3. 把 FFT 后半段的峰值下标换成负时移。
%   4. 超过 dtmax 的时移归零，再用裁剪和补零移动观测波形。
%
% 注意：源码虽然包含 nargin==2 分支，但后面仍使用 dtmax；实际调用应提供四个参数。
so=size(obs);

% 取不小于 N 的最小 2 次幂，通常能提高 FFT 速度。
fftlen=2^nextpow2(so(1));
dim=length(so);

% fft 默认沿第一维计算，因此每条波形独立变换到频域。
fob=fft(obs,fftlen);
fsyn=fft(syn,fftlen);

% 频域相乘再逆变换，得到各波形的循环互相关序列。
R=ifft(fob.*conj(fsyn));

if dim==3
    % 三维输入把第三维分量的相关值相加，让同一台站的分量共用一个时移。
    R=sum(R,3);
end

% max 沿时间轴找相关峰；MATLAB 下标从 1 开始，所以减 1 得到零起点滞后量。
[~,nr]=max(R);
nr=nr-1;
% 循环相关中大于 fftlen/2 的位置代表负时移，减去 fftlen 后恢复符号。
nr(nr>fftlen/2)=nr(nr>fftlen/2)-fftlen;
dt=nr;

if nargin==2
    % 该默认值只补 srate；由于 dtmax 没有默认值，两参数调用随后仍会报未定义变量。
    srate=1;
end

% 把秒阈值换成采样点阈值，拒绝过大的相关峰时移。
xz=(abs(dt)>dtmax*srate);
dt(xz)=0;%mco(xz)=0;

for i=1:length(dt)
    if dim==2
        if dt(i)>0
            % 正时移：删除观测波形开头 dt 点，并在末尾补 0，相当于把波形向前移。
            obs(:,i)=[obs(dt(i)+1:end,i);zeros(dt(i),1)];
        elseif dt(i)<0
            % 负时移：在开头补 0，并裁掉末尾相同点数，相当于把波形向后移。
            obs(:,i)=[zeros(-dt(i),1);obs(1:end+dt(i),i)];
        end
        % 对全部列计算移动后的归一化内积；循环结束时 mco 保存最终结果。
        mco=gfit0(obs,syn);
    elseif dim==3
        if dt(i)>0
            % 三维输入沿第一维移动第 i 个台站的全部分量。
            obs(:,i,:)=[obs(dt(i)+1:end,i,:);zeros(dt(i),1,so(3))];
        elseif dt(i)<0
            obs(:,i,:)=[zeros(-dt(i),1,so(3));obs(1:end+dt(i),i,:)];
        end
        mco=gfit0(obs(:,:),syn(:,:));
    end
end

