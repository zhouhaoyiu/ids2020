function  [slip_model]=Slip_Info(epi,depth,grid,gridsize,source,locasub,locadep,fault,slip,substf,Moment,Mw,Dname,crust,Event,outpagth)
% This function is to sace the rupture models
% SLIP_INFO 把震源参数、速度结构、滑移和子断层源时间函数写入文本结果文件。
%
% 输入
%   epi=[lat,lon]，depth 为震源深度(km)，fault=[strike,dip,rake]。
%   grid=[ndip,nstrike]，gridsize=[dip_km,strike_km]，source 为震源网格下标。
%   locasub 为 nsub*2 子断层中心坐标，locadep 为 nsub*1 深度，slip 为 nsub 个滑移量。
%   substf 为 Nt*nsub 子断层源时间函数；Moment、Mw 为总地震矩和矩震级。
%   Dname 是文件名前缀，crust 每行通常为 [深度,Vp,Vs,密度]，Event 是发震时刻文本。
%   outpagth 是输出目录字符串，末尾需自行包含文件分隔符。
%
% 输出
%   slip_model: nsub*(9+Nt+1) 数值矩阵。前 9 列依次为纬度、经度、深度、
%   倾向尺寸、走向尺寸、走向、倾角、滑动角、滑移；后面是缩放后的 STF。
%   同一内容及元数据写入 [outpagth,Dname,'_Rupture_Info.txt']。
%
% 写文件过程
%   STF 末尾补零并按 10 的整次幂缩放 -> 组合每个子断层的一行数据
%   -> 写事件和断层头信息 -> 写速度结构 -> 写逐子断层破裂模型 -> 写列说明。
if mod(grid(2),2)==0
% 该分支只重新构造同样的二元素 source，不改变数值。
source=[source(1) source(2)];
end

% ----- File name -----
% dex1=find(event==':');dex2=find(event=='(');
% name=event(dex1(1)+11:dex2(1)-13);
name=Dname;
% strcat 生成最终文件名，例如 Dname='event1' 得到 event1_Rupture_Info.txt。
rupname=strcat(name,'_Rupture_Info.txt');
% 
% -----  The sub-fault source time functions -----
% 每条 STF 末尾加一个 0，使输出文件中的时间函数明确回到零；转置后每个子断层占一行。
subtimfun=[substf; zeros(1,size(locasub,1))];
subtimfun=subtimfun';

% Extract the magnitude of the sub-fault source time function
% 用全局最大 STF 的十进制数量级 ma 缩放所有值，表头同时记录 e+ma Nm/s。
ma=max(subtimfun(:));ma=floor(log10(ma));
subtimfun=subtimfun/10^ma;
ma=strcat('[','e+',num2str(ma),'Nm/s',']');

% ----- slip models -----
% 前 8 列几何量按子断层逐行复制，随后接 slip 和完整 STF。
slip_model=[locasub locadep gridsize(1)*ones(size(locadep,1),1) gridsize(2)...
    *ones(size(locadep,1),1) ones(size(locadep,1),1)*fault(1) ones(size...
    (locadep,1),1)*fault(2) ones(size(locadep,1),1)*fault(3) slip(:) subtimfun];

% Seismic moment and moment magnitude
% 两个总量以 2 位有效数字转成文本，仅用于写文件，不改变 slip_model。
Moment=num2str(Moment,2);Mw=num2str(Mw,2);

% ----- Save the rupture model ----
% 以文本写模式创建或覆盖结果文件；fid 是后续 fprintf 使用的文件句柄。
fid=fopen([outpagth,rupname],'wt');
fprintf(fid,'%s\n',['# Origin_time: ',Event]);
% fprintf(fid,'%s\n',strcat('# Version Pubdate:',32,datestr(now,31)));       % Add current time
fprintf(fid,'%s\n','# Reference:  Zhang Y,Wang R,Zschau J,et al.2014a.Automatic imaging of earthquake rupture processes by iterative deconvolution');
fprintf(fid,'%s\n','              and stacking of high-rate GPS and strong motion seismograms.J.Geophys.Res.Solid Earth,119(7),5633-5650.');
fprintf(fid,'%s\n','# Earthquake hypocenter: Latitude[deg] Longitude[deg] Depth[km]');
fprintf(fid,'%0.3f ',epi);
fprintf(fid,'%d\t',depth);
fprintf(fid,'\n%s\n','# Focal mechanism: Strike[deg] Dip[deg] Rake[deg]');
fprintf(fid,'%d\t',fault);
fprintf(fid,'\n%s\n','# Number of sub_faults in [dip,strike] directions');
fprintf(fid,'%d ',grid);
fprintf(fid,'\n%s\n','# Subfault size in [dip,strike] directions');
fprintf(fid,'%d ',gridsize);
fprintf(fid,'\n%s\n','# The index of the hypocentre in [dip,strike] directions');
fprintf(fid,'%d ',source);
fprintf(fid,'\n%s\n','## Velocity_density structure model used for the inversion');
fprintf(fid,'%s\n','# Depth[km]  Vp(km/s)  Vs(km/s)  Density(g/cm^3)');
[m_1,n_1]=size(crust);
% 速度结构逐行输出，列之间用制表符，最后一列后换行。
for i_1=1:m_1
    for j_1=1:n_1
        if j_1==n_1
            fprintf(fid,'%8.4f\n',crust(i_1,j_1));
        else
            fprintf(fid,'%8.4f\t',crust(i_1,j_1));
        end
    end
end
fprintf(fid,'\n%s\n','# Seismic_moment[Nm]  Moment_magnitude[Mw]');
fprintf(fid,'%s ',strcat(Moment,32,32,Mw));
fprintf(fid,'\n%s\n','########   Rupture Model   ########');
fprintf(fid,'%s\n',strcat('# Lat[deg]  Lon[deg]  Depth[km]  Size_dip[km]  Size_strike[km]  Strike[deg]  Dip[deg]  Rake[deg]  Slip[m] Subfault_source_time_functions',ma));
[m_1,n_1]=size(slip_model);
% 破裂模型按列类型选择格式：坐标/深度、整数几何参数、滑移和 STF 使用不同精度。
for i_1=1:m_1
    for j_1=1:n_1
        if j_1==1||j_1==2||j_1==3
            fprintf(fid,'%10.4f\t',slip_model(i_1,j_1));
        elseif j_1==4||j_1==5||j_1==6||j_1==7||j_1==8
            fprintf(fid,'%5d\t',slip_model(i_1,j_1));
        elseif j_1==9
            fprintf(fid,'%10.4f\t',slip_model(i_1,j_1));
        elseif j_1==n_1
            fprintf(fid,'%10.4f\n',slip_model(i_1,j_1));
        else
            fprintf(fid,'%10.7f\t',slip_model(i_1,j_1));
        end
    end
end
% Annotation
% 文件末尾补充各几何列含义，最后关闭文件句柄并刷新缓冲区。
fprintf(fid,'\n%s\n','# Lat[deg]：latitude of center of subfault');
fprintf(fid,'%s\n','# Lon[deg]：longitude of center of subfault');
fprintf(fid,'%s\n','# Depth[km]：depth of center of subfault');
fprintf(fid,'%s\n','# Size_dip[km]：sub-fault size in dip direction');
fprintf(fid,'%s\n','# Size_strike[km]：sub-fault size in strike direction');
fclose(fid);
end

