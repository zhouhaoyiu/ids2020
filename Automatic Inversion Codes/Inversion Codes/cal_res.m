function [nzb,nzc,dig,dt]=cal_res(obr,syn,srate,dtmax)

% The function is to screen the stations based on energy and cross-correlation;
% A waveform station is retained if the average misfit of its three components is less than 0.6.
%% The theory
%    Error signal： v(n)=s1(n)-A*s2(n)  /s1 is the oberved data； s2 is the synthetic data 
% Judgment method： Ev=∑v2(n)=∑[s1(n)-A*s2(n)]2  /Minimum
%                   A=∑s1(n)* s2(n)/E2    /Ev is the error energy，E1 and E2 are the energy of s1 and s2
% ------------------------
%  Input: obr: the observed data
%         syn: the synthetic data
% Output: nzb/nzc：the column number is to be removed
%  
% 中文说明
%   输入 obr、syn 均为 N*(3*S)，列顺序必须是 [全部分量A, 全部分量B, 全部分量C]；
%   S 是台站数，所以总列数必须能被 3 整除。srate 为采样率，dtmax 为最大允许时移(s)。
%
% 输出
%   nzb: 应删除的台站下标。当前判据是该台站至少两个分量的相对残差大于 0.6。
%   nzc: 上述台站对应到 obr/syn 三个分量块中的全部列下标。
%   dig: N*(3*S) 缩放系数矩阵，与合成波形逐元素相乘可得到分量独立缩放结果。
%   dt : shift_obssyn 求得的整数采样点时移。
%
% 计算过程
%   先对齐观测和合成波形；再把三分量拆开。每条合成波形用
%   A=sum(ob.*syn)/sum(syn.^2) 作最小二乘振幅缩放，随后计算相对残差能量。
%   最后的重复下标统计找出至少有两个分量超过阈值的台站。
% 先按互相关移动 obr；返回的 obr 已对齐，dt 记录移动量。
[obr,dt,~]=shift_obssyn(obr,syn,srate,dtmax);


[~,n]=size(obr);                                                        
% 按三个连续列块拆出观测三分量，每块为 N*S。
oba=obr(:,1:n/3);
obb=obr(:,n/3+1:2*n/3);
obc=obr(:,2*n/3+1:end);

% 用同样列规则拆出合成三分量，确保每个台站、分量一一对应。
syna=syn(:,1:n/3);                                                
synb=syn(:,n/3+1:2*n/3);
sync=syn(:,2*n/3+1:end);

% Energy normalization
% A1/A2/A3 是每条合成波形拟合对应观测波形的最小二乘振幅系数。
A1=sum(oba.*syna)./sum(syna.^2);Aa=repmat(A1,size(obr,1),1);
A2=sum(obb.*synb)./sum(synb.^2);Ab=repmat(A2,size(obr,1),1);
A3=sum(obc.*sync)./sum(sync.^2);Ac=repmat(A3,size(obr,1),1);
% 缩放后的三分量重新拼接为 synr；dig 保存相同尺寸的缩放系数。
synd=syna.*Aa;syne=synb.*Ab;synf=sync.*Ac;synr=[synd syne synf];
dig=[Aa,Ab,Ac];

% 每个分量按“误差平方和/观测平方和”计算相对残差；越小表示拟合越好。
resa=sum((oba-synd).^2)./sum(oba.^2);
resb=sum((obb-syne).^2)./sum(obb.^2);
resc=sum((obc-synf).^2)./sum(obc.^2);

mis=0.6; % the threshold value of misfit
% na、nb、nc 分别保存三个分量中残差超过 0.6 的台站编号。
[~,na]=find(resa>mis);
[~,nb]=find(resb>mis);
[~,nc]=find(resc>mis);

% 合并并排序后，相同台站号会相邻；重复出现表示多个分量同时不合格。
nz=[na,nb,nc];nz=sort(nz);nza=diff(nz);
% 第一次提取相邻重复值；若一个台站出现三次，再删除重复项，只保留一个台站号。
a= nza==0;nzb=nz(a);nzc=diff(nzb);
aa= nzc==0;nzb(aa)=[];
% 将台站号展开成原始三分量列号，供调用者一次删除整台站数据。
nzc=[nzb,n/3+nzb,2*n/3+nzb];

end

