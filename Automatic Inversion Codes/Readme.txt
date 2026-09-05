-------------------------------------------------------------------------------------------------------------------------
Part I: Preparation of database of green's functions

In the folder "Green", the procedure "spgrn2013.exe" is a windows procedure to calculate the green's functions and to build the database (Wang et al., 2017).
Before calculation, please set the parameters by editing the input file "offshore.inp". 
The original "offshore.inp" gives an example，one can reset the path (in lines 72 and 139) and run "spgrn2013.exe".
In running "spgrn2013.exe"，just double-click the spgrn2013.exe，and type "offshore.inp".

-------------------------------------------------------------------------------------------------------------------------
Part II: Read the data

Here we use the strong-motion data in folder "strong_motion" which have been downloaded from the NIED website (https://www.kyoshin.bosai.go.jp/kyoshin/search/index_en.html)
The data downloaded is a compressed package，just unzip it for use.
Before reading the data, you should prepare the "earthquake_Info.txt". It contains origin time, hypocantral location, best double couple solution of the earthquake and the type of earthquake.
Then, edit the path of data and run the "Read_SM_NIED.m" on Matlab. Some data files will be organized to the next inversion.

-------------------------------------------------------------------------------------------------------------------------
Part III: Automatic finite-fault inversion

After finishing the Parts I and II，we can carry out the automatic finit-fault inversion. All the codes are contained in folder "Inversion Codes".
By adding the "Inversion Codes" to Matlab path, we can perform the inversion by simply running "main_autoinv.m". Nothing else is required since the inversion has been automated.

The procedure "main_autoinv.m" includes three parts:
1. Processing of the data and Green's functions
2. Automatic Finite-fault Inversion (Station Screening/Inversions by the IDS method/Automatic Updating of the Fault Plane) (Zhang et al., 2014)
3. Plot the slip models

Here we take the 11 March 2011 Mw7.7 earthquake as an example (See Figure 2 in Zheng et al., 2020).

-------------------------------------------------------------------------------------------------------------------------
Note: This procedure requires supports from Singnal Processing Toolbox and Mapping Toolbox of Matlab, and some built-in functions may change with the software versions.

           If you have any problems, please contact us （zhang-yong@pku.edu.cn; zhengxujun@pku.edu.cn）”

-------------------------------------------------------------------------------------------------------------------------
Major references
Zheng, X. et al. 2020. Automatic Inversions of Strong-Motion Records for Finite-Fault Models of Significant Earthquakes in and around Japan.J. Geophys. Res. Solid Earth.
Zhang, Y. et al. 2014. Automatic imaging of earthquake rupture processes by iterative deconvolution and stacking of high-rate GPS and strong motion seismograms. J. Geophys. Res. Solid Earth, 119(7), 5633-5650.
Wang, R. et al. 2017. Complete synthetic seismograms based on a spherical self-gravitating Earth model with an atmosphere-ocean-mantle-core structure. Geophys. J. Int., 210(3): 1739-1764.