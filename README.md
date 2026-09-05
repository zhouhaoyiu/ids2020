# IDS_2020 自动有限断层反演

本仓库保存 IDS_2020 自动有限断层反演方法的公开说明和中文注释版程序，便于查阅、学习与复现。

## 仓库内容

- `IDS_2020_algorithm_guide.pdf`：中文算法与 MATLAB 程序说明。
- `Automatic Inversion Codes/Inversion Codes/`：含中文学习注释的 IDS_2020 MATLAB 反演程序及配套 MAT 文件。
- `Automatic Inversion Codes/green/`：格林函数计算输入示例和 Windows 程序 `spgrn2013.exe`。
- `Automatic Inversion Codes/Readme.txt`：原始英文使用说明。
- `Automatic Inversion Codes/Figure 2.pdf`：原始发布包附带的示例图。

程序来源为北京大学地球物理课题组的[程序下载页](https://pku-geophysics-source.group/htmls/codes.html)及其公开的 [Automatic Inversion Codes 发布包](https://pku-geophysics-source.group/release/Automatic%20Inversion%20Codes.rar)。本仓库未将这些程序标注为仓库维护者原创。

仓库中的 39 个 `.m` 文件在官网发布版基础上增加了中文学习注释，并修复了注释中的乱码标点。2026 年 9 月 6 日以官网 RAR 文件（SHA-256：`26f779b837dea7790c1f4612d1756d0d896825e950b393d3106554d8c7ff42d5`）核对：逐文件去除注释、空白和换行差异后，可执行内容一致。其余 9 个随包文件与该发布包逐字节一致。

## 数据说明

原始发布包中的 `strong_motion/` 包含 NIED K-NET、KiK-net 强震记录。NIED 当前条款禁止再分发，因此本仓库不包含这些波形文件。请在 [NIED K-NET、KiK-net](https://www.kyoshin.bosai.go.jp/en/) 注册并按其条款获取数据；使用时应引用数据 DOI：[`10.17598/NIED.0004`](https://doi.org/10.17598/NIED.0004)。

## 运行

原始说明要求 MATLAB Signal Processing Toolbox 和 Mapping Toolbox。

1. 按 `Automatic Inversion Codes/Readme.txt` 准备强震数据和格林函数数据库。
2. 修改 `main_autoinv.m` 中的 `folder`、`pathg`，以及 `green/offshore.inp` 第 72、139 行的 Windows 示例路径。
3. 将 `Automatic Inversion Codes/Inversion Codes/` 加入 MATLAB 路径，运行 `main_autoinv.m`。

## 许可与引用

原发布页将程序列为“开源代码”，但下载包未附具体许可证。本仓库不额外授予许可；程序、图件、可执行文件和数据的权利归各自权利人。使用前请核对原发布方要求，并引用相关论文：

- Zheng, X. et al. (2020). *Automatic Inversions of Strong-Motion Records for Finite-Fault Models of Significant Earthquakes in and around Japan*. Journal of Geophysical Research: Solid Earth.
- Zhang, Y. et al. (2014). *Automatic imaging of earthquake rupture processes by iterative deconvolution and stacking of high-rate GPS and strong motion seismograms*. Journal of Geophysical Research: Solid Earth, 119(7), 5633-5650.
- Wang, R. et al. (2017). *Complete synthetic seismograms based on a spherical self-gravitating Earth model with an atmosphere-ocean-mantle-core structure*. Geophysical Journal International, 210(3), 1739-1764.
