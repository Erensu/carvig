# carvig
Integrated navigation for INS/GNSS/vision for vehicle positioning

# Introduction
CARVIG is an INS/GNSS integrated navigation algorithm library developed based on RTKLIB and written in C. It is designed for vehicle applications. It is currently being extended with INS/GNSS fusion algorithms that incorporate visual information and LiDAR, and these updates will be added gradually in the future.

CARVIG mainly includes the following functions:
1. INS/GNSS/VO loose-coupling algorithms;
2. SPP, PPP, DGPS, and tightly coupled RTK algorithms;
3. Odometer-assisted navigation;
4. Magnetometer-assisted navigation;
5. Motion constraints such as NHC, ZUPT, and ZARU;
6. Doppler-aided INS/GNSS;
7. Dual-antenna heading assistance;
8. Static and dynamic alignment initialization;
9. Forward and backward INS mechanization;
10. Forward and backward INS/GNSS integrated filtering;
11. Preliminary support for visual-information-assisted positioning and attitude determination;
12. RTS forward/backward smoothing;
13. RTK/INS/GNSS/VO post-processing;
14. Generation of simulated vehicle data for INS/GNSS (position)/VO;
15. Dynamic display of vehicle trajectories.

# Configuration
Before processing data, CARVIG requires editing the relevant configuration files. The format of the configuration file can be referenced from the RTKLIB documentation.

# Processing
By default, CARVIG processing results are stored in the directory of the executable file carvignav under the bin folder.

# Build
CARVIG currently supports compilation only in a Linux environment. The process is the same as compiling a CMake project.

# Run
Qt must be installed for plot functionality.

Run command:
carvignav -o ../example/conf/navlib.conf -m 52716

# References
[1] http://www.rtklib.com

[2] P. D. Groves, Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems

[3] Wuhan MAP Navigation Technology Co., Ltd., http://www.whmpst.com/cn/

[4] http://www.cvlibs.net/

# Contact
Address: Satellite Navigation and Positioning Technology Research Center, Wuhan University, Su Jinglan

QQ: 1971129844

Email: sujinglan@whu.edu.cn

# Note
If you need the professional version of CARVIG, please contact the author.
