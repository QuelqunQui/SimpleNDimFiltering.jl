# SimpleNDimFiltering

[![Build Status](https://travis-ci.org/QuelqunQui/SimpleNDimFiltering.jl.svg?branch=master)](https://travis-ci.org/QuelqunQui/SimpleNDimFiltering.jl)

[![Coverage Status](https://coveralls.io/repos/QuelqunQui/SimpleNDimFiltering.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/QuelqunQui/SimpleNDimFiltering.jl?branch=master)

[![codecov.io](http://codecov.io/github/QuelqunQui/SimpleNDimFiltering.jl/coverage.svg?branch=master)](http://codecov.io/github/QuelqunQui/SimpleNDimFiltering.jl?branch=master)

## ColArray(::Int64,::Int64,::Any)
The purpose of this function is to vectorize ::Colon, and replace one with something else, like a range. To allow to work on a n-1 dimension Array from a n-dimension Array.  Arguments are the wanted length, the position of the non Colon element and the value put in that position.
```julia
X=ones(4,4,4) # you want one slice of the cube of ones
NDims=length(size(X))
X[1,:,:]==X[ColArray(NDims,1,1)...]
X[:,4,:]==X[ColArray(NDims,2,4)...]
X[:,:,1]==X[ColArray(NDims,3,1)...]
X[:,:,1:2]==X[ColArray(NDims,3,1:2)]
```
Which is usefull to standirdize some function towards n-dimension.

## GettingIndex(::Array{Any})

Purpose is to get the indexes of the elements of an n-dim Array in a Matrix (n x #elmts) with a centered origin.
```julia
X=ones(6,4,8)
GettingIndex(X)==[1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3;
                  1 1 1 2 2 2 1 1 1 2 2 2 1 1 1 2 2 2 1 1 1 2 2 2;
                  1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4]
```
## NewLowFilter(::String,::Array,[::Array])

Creating n-dimension filters of the same size than a provided array. Arguments are the type of filter asked and a provided Array, if wanted somme parameters can be added, otherwize those are set to 1.

### Implemented types:

* "Rect",with parameters the size of the n-parallelipipedal filter composed of 1 and 0.
* "Rond", same as previous for a n-ellipsoïdal
* "Triangle", Lineaire filter, parameters can be passed to change the value at the edge of the filter (otherwise = 0)
* "Hanning", n-dim Hanning,
* "Blackman", n-dim Blackman
* "Pcos", Cosinus sum filter with given parameter : ```math p_i -0.5cos(pi*(x-1)/N)+(0.5-p_i)cos(2pi(x-1)/N)```
$$p_i -0.5cos(pi*(x-1)/N)+(0.5-p_i)cos(2pi(x-1)/N)$$ where N is the size of the i dimension in which this is being applied.

## MeanFilter(::String, ::Array, ::Array, [::Array])

Creating filter based on the average of the elements around one elements.
Arguments are the type of mean wanted, the array on which the filter as to be applied and the size of the filter, an additional parameter can be pas in some cases, its set to 1 by default.

### Implemented types

* "Rect", all the elements in the range asked have the same weight on the mean
* "Cross", diagonal elements (elements that do not have a coordonate in common with the considered element) are not taken into account in the average.
* "Triangle", all elements in the range asked are taken but are weighted so that the further they are from the centre elements the less they weight.

### Examples
```julia
X=[0 0 0 0 0;
   0 0 0 0 0;
   0 0 1 0 0;
   0 0 0 0 0;
   0 0 0 0 0]
MeanFilter("Rect",X,[3,3])==[0  0   0   0  0;
                             0 1/9 1/9 1/9 0;
                             0 1/9 1/9 1/9 0;
                             0 1/9 1/9 1/9 0;
                             0  0   0   0  0]

MeanFilter("Cross",X,[3,3])==[0  0   0   0  0;
                              0  0  1/6  0  0;
                              0 1/6 1/3 1/6 0;
                              0  0  1/6  0  0;
                              0  0   0   0  0]

MeanFilter("Triangle",X,[3,3])==[0   0   0   0   0;
                                 0 1/12 1/8 1/12 0;
                                 0 1/8  1/6 1/8  0;
                                 0 1/12 1/8 1/12 0;
                                 0   0   0   0   0]
```
## FilterS1(::Array,::Array,[::Int64])

Arguments are a 2-D array to filter, a 2-D array containing associated noise and optionnaly a parameter to determine the amount of elements eliminated by a last step described below.

FilterS1 is more of an example in hich a cascade of filters are applied using the function described above.  Another type of filter is added wich is a filteing on value in this suite of filters all the elements that are still smaller than 1/20th of the average over the whole array are set to 0.
