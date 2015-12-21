# Advanced-Computer-Graphics
___

## 00_Mitsuba Rendering Scene

#### install mitsuba (Ubuntu 15.04) :

+ 'sudo apt-get install libboost-all-dev'
+ download source-code (all platforms) from https://www.mitsuba-renderer.org/download.html
+ change into build directory of source-code (in terminal)
+ 'cmake ..'
+ 'make install'
+ add binaries to path variable (in ~/.profiles add 'export PATH="$PATH:$HOME/your installation path/build/binaries"' to end of file)
+ relog

#### render scene :

+ cd 'repository directory/Renderscenes/00_Mitsuba'
+ 'source MitsubaPaths.source'
+ 'mtsgui'

___

## 01_Radiosity

#### render scene :

+ make
+ executing radiosity will give 2 output images

___

## 02_PathTracing
![alt tag](https://raw.githubusercontent.com/MasterMoritz/Advanced-Computer-Graphics/master/Images/no_focus.jpg)
![alt tag](https://raw.githubusercontent.com/MasterMoritz/Advanced-Computer-Graphics/master/Images/glass_focus.jpg)
![alt tag](https://raw.githubusercontent.com/MasterMoritz/Advanced-Computer-Graphics/master/Images/metal_focus.jpg) ![alt tag](https://raw.githubusercontent.com/MasterMoritz/Advanced-Computer-Graphics/master/Images/dahoe.jpg)
![alt tag](https://raw.githubusercontent.com/MasterMoritz/Advanced-Computer-Graphics/master/Images/triangle_geometry.jpg)

### render scene :

+ make
+ ./PathTracing [samples per subpixel|int] [samples on lens|int] [focal distance|cm] [aperture|f-stops]
+ e.g.
  + ./PathTracing 40
  + ./PathTracing 20 10
  + ./PathTracing 1 1 248.6 2.8
+ note: to disable DOF, simply use a high f-stop value (e.g. 2000)
