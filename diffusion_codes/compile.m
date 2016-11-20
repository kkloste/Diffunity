function compile
mypath = fileparts(mfilename('fullpath'));
curdir = pwd();
try
    cd(mypath);
    if ismac
       	%mex -O CXXFLAGS="\$CXXFLAGS -Wall -I. -std=c++11" -largeArrayDims hkgrow_mex.cpp
        mex -O CXXFLAGS="\$CXXFLAGS -Wall -I. -std=c++11" -largeArrayDims gendiff_mex.cpp
        %mex -O CXXFLAGS="\$CXXFLAGS -Wall -I. -std=c++11" -largeArrayDims pprgrow_mex.cc
       	mex -O CXXFLAGS="\$CXXFLAGS -Wall -I. -std=c++11" -largeArrayDims reg_power_mex.cpp
        mex -O CXXFLAGS="\$CXXFLAGS -Wall -I. -std=c++11" -largeArrayDims ../util/sweepcut_mex.cpp
        mex -O CXXFLAGS="\$CXXFLAGS -Wall -I. -std=c++11" -largeArrayDims ../util/sparse_degpower_mex.cpp
    else
        %mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -I. -largeArrayDims hkgrow_mex.cpp
       	mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -I. -largeArrayDims gendiff_mex.cpp
        %mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -I. -largeArrayDims pprgrow_mex.cc
        mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -I. -largeArrayDims reg_power_mex.cpp
        mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -I. -largeArrayDims ../util/sweepcut_mex.cpp
        mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -I. -largeArrayDims ../util/sparse_degpower_mex.cpp
    end
    cd(curdir);
catch me
    cd(curdir)
    rethrow(me)
end
