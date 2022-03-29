# BaVI-pose-tracking
Particle Filter and least squares based tracking of climbers with IMUs and route maps.

This repository contains the code for the publication: 

Event-Domain Knowledge in Inertial Sensor Based State Estimation of Human Motion by
Tom L. Koller, Tim Laue and Udo Frese (all University of Bremen)
Currently available as preprint only:
https://easychair.org/publications/preprint/rhMf

Please read the publication for an explanation of the algorithm in this repository. 


This implementation uses the orientation estimator of: https://github.com/dlaidig/qmt
which was published under MIT license by Thomas Seel and Stefan Ruppin

Please contact me if you are interested in using this code or reproducing my results but have trouble to run the code.

# Dependencies
(Please look up the dependencies of the repositories)
1.  ADEKF https://github.com/TomLKoller/ADEKF
1.  ADEKF VIZ https://github.com/TomLKoller/ADEKF_VIZ
1.  libconfig++ https://github.com/hyperrealm/libconfig
1.  glog https://github.com/google/glog
3.  JKQTPlotter https://github.com/jkriege2/JKQtPlotter
4.  QT5  https://doc.qt.io/qt-5/qt5-intro.html
5.  VTK https://vtk.org/
6.  Boost https://www.boost.org/
7.  Eigen3 https://eigen.tuxfamily.org/index.php?title=Main_Page
8.  Json reader https://github.com/nlohmann/json (As submodule of the repository)

You will need to adapt the CMakeLists.txt to point to your installation locations
