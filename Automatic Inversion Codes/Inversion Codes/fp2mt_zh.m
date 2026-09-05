function moment=fp2mt_zh(m,x,y,z)
%moment=fp2mt_zh(m,x,y,z)
%m: scalar seismic moment
%x: strike
%y: dip
%z: rake
%moment: seismic moment tensor
%
% 中文说明
%   输入 m 为标量地震矩，x、y、z 分别为走向、倾角和滑动角，角度单位为度。
%   四个输入可以是标量或等长向量；函数先统一转成行向量。
%   输出 moment 为 6*N，六行依次为 [Mxx,Mxy,Mxz,Myy,Myz,Mzz]，
%   N=length(x)，每一列对应一组震源参数。
%
% 计算过程
%   1. 把走向、倾角、滑动角从度转换为弧度，因为 MATLAB 的 sin、cos 接收弧度。
%   2. 预先计算一倍角和二倍角三角函数，供六个张量分量重复使用。
%   3. 按双力偶震源的走向—倾角—滑动角关系，逐列生成一个矩张量。
%   所有乘法使用 .*，因此可一次处理多组震源参数。
%
% 例：fp2mt_zh(1,0,90,0) 返回一列六分量；m 改为 2 时六个分量也都乘 2。
% 注意：坐标轴方向和六分量顺序必须与格林函数数据库的约定一致。
% 统一为行向量，保证每一列输出对应一组参数。
m=m(:)';
x=x(:)';
y=y(:)';
z=z(:)';
% 先分配 6*N 输出矩阵。
moment=zeros(6,length(x));

% 度转弧度：180 度对应 pi 弧度。
x=x.*pi./180;
y=y.*pi./180;
z=z.*pi./180;

% 缓存三角函数，避免在后面六条公式中反复求值。
sinx=sin(x);cosx=cos(x);sin2x=2.*sinx.*cosx;cos2x=cosx.*cosx-sinx.*sinx;
siny=sin(y);cosy=cos(y);sin2y=2.*siny.*cosy;cos2y=cosy.*cosy-siny.*siny;
sinz=sin(z);cosz=cos(z);
% 以下六式输出对称矩张量的六个独立元素；每列对应一组 x、y、z。
moment(1,:)=-m.*(cosz.*siny.*sin2x+sinz.*sin2y.*(sinx.*sinx));
moment(2,:)=m.*(cosz.*siny.*cos2x+sinz.*sin2y.*sin2x./2);
%moment(2,1)=moment(1,2);
moment(3,:)=-m.*(cosz.*cosy.*cosx+sinz.*cos2y.*sinx);
%moment(3,1)=moment(1,3);
moment(4,:)=m.*(cosz.*siny.*sin2x-sinz.*sin2y.*cosx.*cosx);
moment(5,:)=m.*(-cosz.*cosy.*sinx+sinz.*cos2y.*cosx);
%moment(3,2)=moment(2,3);
moment(6,:)=m.*sinz.*sin2y;
return
% ================================end======================================
