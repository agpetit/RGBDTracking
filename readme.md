# RGBDTracking plugin

<p align="center"><img src="doc/img/realsense.png" width="70%" /><br><br></p>


## Introduction
**RGBDTracking** is a [SOFA](www.sofa-framework.com) plugin to register and track deformable objects using an RGB-D sensor in real-time.

The framework relies on a prior visual segmentation of the object in the image. The segmented point cloud is registered first in a rigid manner and then by
non-rigidly fitting the mesh, based on the Finite Element Method to model elasticity, and on vision-based external forces exerted on the mesh.

The developped methods are described in the papers [RAS 2017](http://wpage.unina.it/antoine.petit/RAS_2016_apetit.pdf) and [IROS 2015](http://wpage.unina.it/antoine.petit/deformable_object_tracking_vf.pdf)

## Installation

Using Ubuntu up to 16-04 LTS as OS is preferable to benefit form all the features of the plugin (CUDA based segmentation).

The plugin is developed with the SOFA platform, we refer to the [SOFA documentation](https://www.sofa-framework.org/documentation) for its up-to-date (master branch) installation (on Linux here).


Use gcc-5 preferably, apart from the optional custom CUDA installation which requires gcc-4.8 (see below how you can install CUDA and switch between compilers)

### Step 1: dependencies
The plugin itself has the following dependencies:

#### OpenCV
Installation from [source] (https://docs.opencv.org/trunk/d7/d9f/tutorial_linux_install.html)
OpenCV 3.2.0 is recommended.

#### PCL 1.8.1

Required dependencies for PCL:

```
sudo apt-get install doxygen
sudo apt-get install mpi-default-dev openmpi-bin openmpi-common
sudo apt-get install libflann1.8 libflann-dev
sudo apt-get install libeigen3-dev
sudo apt-get install libboost-all-dev
sudo apt-get install libvtk6.2-qt4 libvtk6.2 libvtk6-dev
sudo apt-get install libqhull*
sudo apt-get install libusb-dev
sudo apt-get install libgtest-dev
sudo apt-get install git-core freeglut3-dev pkg-config
sudo apt-get install build-essential libxmu-dev libxi-dev
sudo apt-get install libusb-1.0-0-dev graphviz mono-complete
sudo apt-get install libproj-dev
sudo apt-get install qt-sdk openjdk-8-jdk openjdk-8-jre
sudo apt-get install phonon-backend-gstreamer
sudo apt-get install phonon-backend-vlc
```

Installation

```
Sources for PCL 1.8.1 can be found here https://github.com/PointCloudLibrary/pcl/releases/tag/pcl-1.8.1
Unzip and then:
cd pcl-pcl-1.8.1
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=None -DBUILD_GPU=ON -DBUILD_apps=ON -DBUILD_examples=ON ..
make
sudo make install
```

#### Freeimage

```
sudo apt-get install libfreeimage3 libfreeimage-dev
```

#### ViSP

```
sudo apt-get install libvisp-dev
```

#### Cuda (optional, if you have an Ubuntu > 16.04 LTS or no available NVIDIA Graphic card, directly go to step 2)

Cuda (<= 7.0) is required to enable some features for the CUDA based segmentation (a CUDA free CPU implementation of the segmentation method, based on OpenCV, is also availble)

For Ubuntu 16.04 and Cuda 7.0:


First uninstall current Nvidia driver version

```
sudo apt-get --purge remove nvidia-*
```
No need to create an xorg.conf file. If you have one, remove it (assuming you have a fresh OS install).
```
sudo rm /etc/X11/xorg.conf
```

Blacklist the "nouveau" driver

```
echo -e "blacklist nouveau\nblacklist lbm-nouveau\noptions nouveau modeset=0\nalias nouveau off\nalias lbm-nouveau off\n" | sudo tee /etc/modprobe.d/blacklist-nouveau.conf
echo options nouveau modeset=0 | sudo tee -a /etc/modprobe.d/nouveau-kms.conf
sudo update-initramfs -u
```

Install Nvidia driver 384.90
```
sudo apt-get install nvidia-384-dev
```

Cuda 7.0 supports up to gcc-4.8 compiler.  Install it, assuming you have gcc-5 already installed and used to install the previously mentioned libraries:
```
sudo apt-get install gcc-4.8 g++-4.8

sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 10
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-5 10

sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 20
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.8 20
```
You can then manually switch with gcc-5

```
sudo update-alternatives --config gcc
sudo update-alternatives --config g++
```

Download CUDA toolkit 7.0 (officially supported by Ubuntu 14.04) [here](https://developer.nvidia.com/cuda-toolkit-70) as a runfile

Go to tty (Ctrl + Alt + F5) and do:
```
sudo service lightdm stop
```

In your Downloads folder do:
```
chmod +x cuda_7.0.28_linux.run
sudo ./cuda_7.0.28_linux.run --no-opengl-lib
```
Answers to cuda installation requests:
```
Do you accept the previously read EULA? (accept/decline/quit): accept
You are attempting to install on an unsupported configuration. Do you wish to continue? ((y)es/(n)o) [ default is no ]: yes
Install NVIDIA Accelerated Graphics Driver for Linux-x86_64 352.39? ((y)es/(n)o/(q)uit): no
Install the CUDA 7.0 Toolkit? ((y)es/(n)o/(q)uit): yes
Enter Toolkit Location [ default is /usr/local/cuda-7.0 ]: press enter or specify a path
Do you want to install a symbolic link at /usr/local/cuda? ((y)es/(n)o/(q)uit): y
Install the CUDA 7.0 Samples? ((y)es/(n)o/(q)uit): y
Enter CUDA Samples Location [ default is /home/user ]: presse enter or specify a path
```
Reboot:

```
sudo reboot
```

Set environment path variables in .bashrc:

```
export PATH=/usr/local/cuda-7.0/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-7.0/lib64:$LD_LIBRARY_PATH
```

To check your installation, compile the CUDA Samples

```
cd home/user/NVIDIA_CUDA-7.0_Samples
make
cd bin/x86_64/linux/release/
./deviceQuery
```
The configurations and versions of your NVIDIA graphic card and drivers and CUDA toolkit should then appear.

Switch back to gcc-5

### Step 2: Compilation of the plugin

Edit the SOFA plugin CMakeLists to add the plugin `gedit /home/.../sofa/master/src/applications/plugins/CMakelists.txt`
and add the line `sofa_add_plugin(RGBDTracking RGBDTracking)`

In your build directory
```
cd /home/.../sofa/master/build
ccmake ../src
```
In the cmake GUI, activate the plugin 'RGBDTracking' and the plugin 'image', set 'OpenCV_DIR' to the build directory of your compiled OpenCV 3.2,
and if CUDA 7.0 is installed, set 'CUDA_HOST_COMPILER' to /usr/bin/g++-4.8

Configure ang generate your cmake and run `make`.
