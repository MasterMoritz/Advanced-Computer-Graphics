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
![Alt text](images/dahoe.jpg?raw=true "happy little accident")

### render scene :

+ make
+ ./PathTracing [samples per subpixel|int] [samples on lens|int] [focal distance|cm] [aperture|f-stops]
+ e.g.
  + ./PathTracing 40
  + ./PathTracing 20 10
  + ./PathTracing 1 1 248.6 2.8
