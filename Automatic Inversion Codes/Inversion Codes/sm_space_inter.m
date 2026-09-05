function x=sm_space_inter(grid,ndot)
% SM_SPACE_INTER 建立断层面空间平滑的稀疏差分矩阵。
%
% 输入
%   grid: 二维时为 [ndip,nstrike]；三维时为 [nvertical,ndip,nstrike]。
%   ndot: 每个子断层包含的参数数目；IDS_2020 传入 1，表示每块只有一个总滑移量。
%
% 输出
%   x: (prod(grid)*ndot) 方阵。对参数向量 v，x*v 的每个元素等于
%      “同一参数在相邻子断层上的平均值 - 当前子断层值”。
%      因此 norm(x*v) 越小，空间分布越平滑。
%
% 参数排列
%   输出矩阵假定 v=[所有子断层的参数1;所有子断层的参数2;...]。
%   sm_space_dex 先按“每个子断层连续放 ndot 个参数”生成下标，
%   con2con 再把它换成上述按参数分块的顺序。
sizemat=grid;lensub=ndot;
% 得到稀疏矩阵每个非零项的原始行下标 mm、列下标 nn 和权重 num。
[mm,nn,num]=sm_space_dex(sizemat,lensub); 

nsub=prod(sizemat);
% 同一种换序同时作用于行、列，矩阵表达的邻接关系不变。
mm0=con2con(mm,nsub,lensub);
nn0=con2con(nn,nsub,lensub);

x=sparse(mm0,nn0,num);


function mm0=con2con(mm,nsub,lensub)
% CON2CON 把“子断层优先”线性下标转换为“参数优先”线性下标。
% 例：2 个子断层、每块 3 个参数时，旧顺序 [s1p1,s1p2,s1p3,s2p1,s2p2,s2p3]
% 变成 [s1p1,s2p1,s1p2,s2p2,s1p3,s2p3]。
% sub_dex 是旧下标对应的子断层号。
sub_dex=ceil(mm/lensub); 
% dian_dex 是该子断层内的参数号；整除时 mod 为 0，需改回 lensub。
dian_dex=mod(mm,lensub);

dian_dex(dian_dex==0)=lensub;

% 新下标先跨过前面完整的参数块，再加当前子断层号。
mm0=nsub*(dian_dex-1)+sub_dex;

function [mm,nn,num]=sm_space_dex(sizemat,lensub)
%   clear
%   sizemat=[10,20];lensub=40;flag=1;
% =========================================================================
%  x=con_sm_mat_invt_multi0(sizemat,lensub,flag)
%zha
%  For smoothing the sliprate on the fault, x can be used in waveform
%  inversion of the rupture process, where lambda*x is a part of the
%  green's matrix, and lambda is a smoothing factor---whose value can be
%  determined by seeking the minimum of ABIC (Akaike's Bayesian Information
%  Criterion). (Akaike,1980)
%
%  sizemat is the nums of the fault's grids.sizemat(1) is for dip direction
%  and sizemat(2) for strike direction when dim==2. while dim==3,sizemat(1)
%  is for vertical direction, sizemat(2) for dip direction, and sizemat(3)
%  for strike direction.
%
%  lensub is the length of STFs of the sub-faults.
%
%  This program also can be used in moment tensor reversion study, where
%  lensub is the number of free parameters of the moment tensor for each
%  subfault.Usually lensub=6;
%
%  flag==1, edge elements donot change
%     else, edge elements are constrained to be zeros
%
% 本文件中的实际调用只有 sizemat、lensub 两个参数，所以 nargin==2 总是令 flag=1。
% 这表示边缘点用现有邻点的平均值，不把断层面外部强制设为零；两个 else 分支均为空或注释代码。
%
%  See also con_sm_mat_invt, con_sm_mat, sm_mat, sm_mat1,
%  con_sm_mat_invt_multi.
%--------------------------------------------------------------------------
%                          Zhang Yong , Chen Yuntai and Xu Lisheng
%                             2006/09/25/   IGCEA (dim==2)
%                             2006/10/24/   Peking University (dim==3)
%==========================================================================
sb=sizemat; 
dim=length(sb);
if nargin==2;flag=1;end
%x=sparse(prod(sb)*lensub,prod(sb)*lensub);

if dim==2
    % 二维每个点最多涉及自身和四邻点，所以最多预分配 5*nsub*lensub 个非零项。
    mm=zeros(prod(sizemat)*5*lensub,1);
    nn=zeros(prod(sizemat)*5*lensub,1);
    num=zeros(prod(sizemat)*5*lensub,1);
    iter=0;
    
    nx=zeros(5,1);
    for i=1:sb(1) % Dip
        for j=1:sb(2) % Strike
            %-------------------------------------
            %          Strike: --------->        |
            %        Dip       2                 |
            %         |    1   3   5             |
            %         |        4                 |
            %        \|/                         |
            %-------------------------------------

            nx(1)=lensub*( (j-2)*sb(1)+ i -1 );
            nx(2)=lensub*( (j-1)*sb(1)+ i -1-1 );
            nx(3)=lensub*( (j-1)*sb(1)+ i -1 );
            nx(4)=lensub*( (j-1)*sb(1)+ i -1+1 );
            nx(5)=lensub*( (j  )*sb(1)+ i -1 );
            % nx(3) 是当前点参数块起始偏移；1、2、4、5 分别是左、上、下、右邻点。
            % find the positions of the four elements
            if flag==1
                %-------------------------------------------------
                % 4 edges:
                if i==1
                    nx(2)=NaN;
                end
                if i==sb(1)
                    nx(4)=NaN;
                end
                if j==1
                    nx(1)=NaN;
                end
                if j==sb(2)
                    nx(5)=NaN;
                end
                %--------------------------------------------------
                % NaN 和负下标不会通过筛选；xn 只保留当前点及网格内邻点。
                xn=find(nx>-1e0); 
                for k=1:lensub
                    for zz=1:length(xn)
                        iter=iter+1;
                        if xn(zz)==3 
                            % 当前点系数为 -1。
                            mm(iter)=nx(3)+k; 
                            nn(iter)=nx(3)+k; 
                            num(iter)=-1; 
                        else
                            % 每个有效邻点平分 +1，总邻点权重之和为 +1。
                            mm(iter)=nx(3)+k;
                            nn(iter)=nx(xn(zz))+k;
                            num(iter)=1./(length(xn)-1);
                        end
                    end
                end
            else
            end % if or not : set edge elements to be zeros
        end % for i=1:sb(2)
    end % for i=1:sb(1)
elseif dim==3
    % 三维每个点最多涉及自身和六个面邻点，因此每行最多 7 个非零项。
    mm=zeros(prod(sizemat)*7*lensub,1);
    nn=zeros(prod(sizemat)*7*lensub,1);
    num=zeros(prod(sizemat)*7*lensub,1);
    iter=0;
    
    nx=zeros(7,1);
    for i=1:sb(1)
        for j=1:sb(2)
            for k=1:sb(3)
                % ---------------------------------------------------------
                %                           2      3
                %                           |    /
                %                           | /
                %                  1 ------ 4 ------- 7
                %                         / |
                %                      /    |
                %                    5      6
                %
                %                ___________________
                %               /                  /|
                %           i  /                  / |
                %             /__________________/  |
                %            |       k           |  |
                %          j |                   |  /
                %            |                   | /
                %            |___________________|/Left: the 3d fault model
                %
                %     3,5 ---- i th direction , the vertical direction
                %     2,6 ---- j th direction , also the dip direction
                %     1,7 ---- k th direction , also the strike direction
                %     nx(1) < nx(2) < nx(3) < nx(4) < nx(5) < nx(6) < nx(7)
                % ---------------------------------------------------------
                nx(1)=lensub*( (k-2)*sb(1)*sb(2) +(j-1)*sb(1) + i-1);
                nx(2)=lensub*( (k-1)*sb(1)*sb(2) +(j-2)*sb(1) + i-1);
                nx(3)=lensub*( (k-1)*sb(1)*sb(2) +(j-1)*sb(1) + i-1-1);
                nx(4)=lensub*( (k-1)*sb(1)*sb(2) +(j-1)*sb(1) + i-1);
                nx(5)=lensub*( (k-1)*sb(1)*sb(2) +(j-1)*sb(1) + i-1+1);
                nx(6)=lensub*( (k-1)*sb(1)*sb(2) +(j  )*sb(1) + i-1);
                nx(7)=lensub*( (k  )*sb(1)*sb(2) +(j-1)*sb(1) + i-1);
                % nx(4) 为当前体素，其他六项对应三个坐标方向的前后邻点。
                if flag==1
                    %---------------------------------------------------
                    % smooth edge elements:
                    % 6 sides:
                    if i==1
                        nx(3)=NaN;
                    end
                    if i==sb(1)
                        nx(5)=NaN;
                    end
                    if j==1
                        nx(2)=NaN;
                    end
                    if j==sb(2)
                        nx(6)=NaN;
                    end
                    if k==1
                        nx(1)=NaN;
                    end
                    if k==sb(3)
                        nx(7)=NaN;
                    end
                    %-------------------------------------------------
                    xn=find(nx>-1e5);
                    for l=1:lensub
                        %                         x(nx(4)+l,nx(xn)+l)=1./(length(xn)-1);
                        %                         x(nx(4)+l,nx(4)+l)=-1;
                        for zz=1:length(xn)
                            if xn(zz)==4 
                                % 当前体素权重为 -1。
                                iter=iter+1;
                                mm(iter)=nx(4)+l;
                                nn(iter)=nx(4)+l;
                                num(iter)=-1;
                            else
                                % 所有网格内面邻点均分 +1 权重。
                                iter=iter+1;
                                mm(iter)=nx(4)+l;
                                nn(iter)=nx(xn(zz))+l;
                                num(iter)=1./(length(xn)-1);
                            end
                        end
                    end
                else
                    % set edges to be zeros
%                     for l=1:lensub
%                         if i==1||i==sb(1)||j==1||j==sb(2)||k==1||k==sb(3)
%                             % all sides: contain all edges and all corners
%                             x(nx(4)+l,nx(4)+l)=1e3;
%                         else
%                             x(nx(4)+l,nx(4)+l)=-1;
%                             x(nx(4)+l,nx([1,2,3,5,6,7])+l)=1/6;
%                         end
%                     end
                end % for flag
            end % for k=1:sb(3)
        end % for j=1:sb(2)
    end % for i=1:sb(1)
end % for dim

% 删除预分配中未写入的尾部零下标，避免 sparse 把它们解释成非法第 0 行/列。
mm(iter+1:end)=[];
nn(iter+1:end)=[];
num(iter+1:end)=[];

%================================end=======================================
%
