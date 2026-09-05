function stf_out=peakstf(stf_in,isarea)
% This fuction is to find the maximum positive wavelet of the stacked apparent source time functions
% PEAKSTF 从每列表观源时间函数中保留一个主要正波包，其余位置置 0。
%
% 输入
%   stf_in: Nt*nsub，每列是一条源时间函数。
%   isarea: 可选；不提供时选择“包含最高峰的连续正波包”，提供任意值时选择
%           “由零值分隔、累计面积最大的正波包”。参数的具体数值不参与判断。
%
% 输出
%   stf_out: 与 stf_in 尺寸相同，只保留每列选中的一个正波包。
%
% 为什么先去掉负值
%   IDS 此处只提取正的矩率贡献；所有负样点先设为 0，同时零值也成为波包边界。
% 把负值截为 0，输出只可能包含非负样点。
stf_in(stf_in<0)=0;
nsub=size(stf_in,2);
stf_out=zeros(size(stf_in));

if nargin==1 % seek the maximum value 
    % mstf 是每列最大值，nstf 是该最大值首次出现的时间下标。
    [mstf,nstf]=max(stf_in);
    
    ls=size(stf_in);
    for i=1:nsub
        if mstf(i)==0
            continue;
        end
        stf0=stf_in(:,i);
        
        % 从峰值向前和向后搜索第一个 0，确定包含最高峰的连续区间。
        s1=stf0(nstf(i):-1:1);
        s2=stf0(nstf(i):end);
        
        d1=find(s1==0);
        d2=find(s2==0);
        
        if isempty(d1)
            % 峰前没有 0 时，从整列第 1 点开始保留。
            dex1=1;
        else
            dex1=nstf(i)-d1(1)+1;
        end
        
        if isempty(d2)
            % 当前代码把 size(stf_in) 整体赋给 dex2；多列且峰后无 0 时可能触发非标量下标错误。
            dex2=ls;
        else
            dex2=nstf(i)+d2(1)-1;
        end
        
        stf_out(dex1:dex2,i)=stf0(dex1:dex2);
    end
else % seek the maximum area
    % 累计和让任意两个边界之间的面积可用一次相减得到。
    leng=size(stf_in,1);
    cstf=cumsum(stf_in);  
    
    mstf=max(stf_in);
    for i=1:nsub
        if mstf(i)==0
            continue;
        end  

%% Xujun Zheng
         stf0=cstf(:,i);
         % 找全部零值下标；相邻零值之间没有正面积，只有不连续零值之间形成候选波包。
         idx0=find(stf_in(:,i)==0);
         idx=diff(idx0);
         idx_1=find(idx~=1);
         idx_2=idx_1+1;
         d1=idx0(idx_1);
         d2=idx0(idx_2);
         area=stf0(d2)-stf0(d1);
         % 选择累计幅值最大的候选区间，并把原始波包复制到输出。
         [~,na]=max(area);
         stf_out(d1(na):d2(na),i)=stf_in(d1(na):d2(na),i);
         
    end
end

