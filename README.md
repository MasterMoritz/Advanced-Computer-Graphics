# Advanced-Computer-Graphics
___

#2 00_Mitsuba Rendering Scene

#4 install mitsuba (Ubuntu 15.04) :

+ 'sudo apt-get install libboost-all-dev'
+ download source-code (all platforms) from https://www.mitsuba-renderer.org/download.html
+ change into build directory of source-code (in terminal)
+ 'cmake ..'
+ 'make install'
+ add binaries to path variable (in ~/.profiles add 'export PATH="$PATH:$HOME/<your installation path>/build/binaries')
+ relog

#4 render scene :

+ cd '<repository directory>/Renderscenes/00_Mitsuba'
+ source MitsubaPaths.source
+ mtsgui
___
