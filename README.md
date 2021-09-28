# inplanefkashmethod

This is the program for calculating the the thermal conductivity using in plane analytical model.

## Build 

To build this project you need a compiler that support C++17 standard. To build simply run the following commands in bash:
	
	$ git clone https://github.com/AlexKiryushkin/inplaneflashmethod
	$ cd inplaneflashmethod
	$ mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -A Win32 ..
	$ cmake --build ./ --config Release


## Run

To run the program find the executable. By default the program reads parameter values from parameters.txt file. This file is located in the repository directory and is copied to executable folder as a post-build event. To change parameters, rewrite parameters.txt in the repository directory and it will be copied after the next build. To run the program you need to specify a path to file with experimental data like this:
	
	$ ./ipfm.exe --input-file=path_to_file

where `path_to_file` is path to the current expermintal data. If you provide relative path, path will be computed with respect to directory of executable. Note: if path contains whitespaces, you should wrap it around with quotes:
	
	$ ./ipfm.exe --input-file="C:/Program Files/data.csv"

Also note that there should be no whitespaces before and after '=' sign. Otherwise parsing will fail.

The format of output file is csv-like. See example:
```
Untitled_Time;Untitled
0;0.866876826
1E-06;0.866644457
2E-06;0.866417613
3E-06;0.86619624
```
Sample files are located in data folder of the repository.

## Algorithm

The program calculates the values according to the formulas from Gembarovich, J., et al., "In-plane thermal diffusivity measurement of highly thermal conductive thin films by the flash method", 2018.