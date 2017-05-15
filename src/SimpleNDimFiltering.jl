__precompile__()

module SplNDimsFiltering
############################################################################
# Col Array ceates a vector of type Colon where one of the elements is
# replaced by another value typically a Range type but not mandatory

function ColArray(L::Int64,Pos::Int64,Val)
  Arr=[]
  for j=1:L
    if j==Pos
      V=Val
      push!(Arr,V)
    else
      push!(Arr,[:][1])
    end
  end
return Arr
end
##############################################################################
# Getting a matrix 2xn of the indexes of the point of an n-dimension Array

function GettingIndex(OneImg)
  NDim=length(size(OneImg))
  WindowPc=ones((convert(Vector{Int64},[size(OneImg)...]/2)...))

  VWinPc=[WindowPc...]
  IndexVWP=zeros(Int64,NDim,length(VWinPc))
  SPrec=1
  NRepVec=length(VWinPc)

  for i=1:NDim
    S=size(WindowPc,i)

    NRepVal=length(VWinPc)/NRepVec
    SNew=SPrec*S
    NRepVec=length(VWinPc)/SNew
    #println(NRepVal)
    #println("S=$S,SPrec=$SPrec,Vec=$NRepVec")
    for j=1:NRepVec
      for q=1:SNew
        k=floor((q+NRepVal-1)/NRepVal)
        IndexVWP[i,convert(Int64,q+(j-1)*S)]=k
      end
    end
    SPrec=SNew
  end
  return(IndexVWP)
end

###########################################################################
# n-dimension low pass FilterS1

function NewLowFilter(FType::String,OneImg::Array,VParam=[1]::Array)
  NDim=length(size(OneImg)) # Symetry on 2^NDims
  ImpDim=[]
  if sum(!isinteger.([size(OneImg)...]/2))>0
    ImpDim=find(!isinteger.([size(OneImg)...]/2))
    for dd=ImpDim
      V=1:size(OneImg,dd)-1
      Col=ColArray(NDim,dd,V)
      OneImg=OneImg[Col...]
    end
  end
  WindowPc=ones((convert(Vector{Int64},[size(OneImg)...]/2)...))
  # Creation of the Piece of Window
  if sum(FType.==["Rect","Rond"])>=1
    Sharpe=true
  else
    Sharpe=false
  end
  IndexVWP=GettingIndex(OneImg)
  LI=size(IndexVWP,2)
  WinV=zeros(LI)
  if Sharpe
    if FType=="Rect"
      if length(VParam)!=NDim
        VParam=VParam[1]/2*ones(NDim) # Same length for all dimentions
      end
      I2=zeros(size(IndexVWP))
      for i=1:NDim
        MaxL=maximum(IndexVWP[i,:])
        I2[i,:]=IndexVWP[i,:].>=(MaxL-VParam[i]/2)
      end
      WinV=(sum(I2,1).==NDim)+0
    elseif FType=="Rond"
      if length(VParam)!=NDim
        VParam=VParam[1]*ones(NDim) # Same radius for all dimentions
      end
      Xyz=zeros(NDim)
      for pp=1:LI
        for i=1:NDim
          Xyz[i]=IndexVWP[i,pp]-maximum(IndexVWP[i,:])
        end
        if sum((Xyz./VParam).^2)<1
          WinV[pp]=1
        end
      end
    end
  else
    if FType=="Triangle"
      if length(VParam)!=NDim
        VParam=zeros(NDim,1)
      end
      WinV=prod((IndexVWP-1).+VParam,1) #Allowing different slope per direction
    elseif FType=="Hanning"|| FType=="Pcos" || FType=="Blackman"
      # Some usual Filter
      if FType=="Hanning"
        VParam=0.5*ones(NDim)
      elseif FType=="Blackman"
        VParam=0.42*ones(NDim)
      elseif FType=="Pcos"
        if length(VParam)!=NDim
          VParam=VParam[1]/2*ones(NDim) # Same length for all dimentions
        end
      end
      I2=zeros(size(IndexVWP))
      for i=1:NDim
        N=(size(OneImg)[i]/2-1)
        I2[i,:]=VParam[i]-0.5*cos.(pi*(IndexVWP[i,:]-1)/N)+(0.5-VParam[i])*cos.(2*pi*(IndexVWP[i,:]-1)/N)
      end
      WinV=prod(I2,1)
    end
  end
  # Reshape
  WindowPc=reshape(WinV,convert(Vector{Int64},[size(OneImg)...]/2)...)
  #--------------------------------------------------------------------------
  # Piece of the window has been made
  W=WindowPc
  RParam=ones(size(VParam))
  for i=1:NDim
    Wfp=flipdim(W,i)
    # All dimensions where even
    if !isempty(findin(ImpDim,i).>0) && size(OneImg,i)>1
      WMid=W[ColArray(NDim,i,size(W,i))...]
      if NDim==2 && i==1
        WMid=transpose(WMid) # due to reading in matrix always yielding Array{Ax1} and so no row-vector
      end
      W=cat([i],W,WMid)
      W=cat([i],W,Wfp)
      if isinteger(VParam[i])
        RParam[i]=VParam[i]+1
      else
        RParam[i]=VParam[i]
      end
    else
      W=cat([i],W,Wfp)
    end
  end
  if FType=="Triangle" && sum(ImpDim.>0)==NDim
    W[convert(Vector{Int64},[size(OneImg)...]/2+1)...]=W[convert(Vector{Int64},[size(OneImg)...]/2+1)...]+1
  end
  return W
end

############################################################################

# Mean based FilterS1

function MeanFilter(FType::String,OneImg::Array,MSize::Array,VParam=[1]::Array)
  # -------------------------------- Head --------------------------------------
  NDim=length(size(OneImg))
  for mm=1:NDim
    if isinteger(MSize[mm]/2)
      MSize[mm]=MSize[mm]+1 # Only non even mean have meanings here
    end
  end
  ZBorder=convert(Vector{Int64},(MSize-1)/2) # number of Zeroes needed to counter border effect
  BImg=zeros(([size(OneImg)...]+2*ZBorder)...)
  FirstElmt=ZBorder+ones(size(ZBorder))
  VRange=[]
  for rr=1:NDim
   R=ZBorder[rr]+1:ZBorder[rr]+size(OneImg,rr)
   push!(VRange,R)
  end
  BImg[VRange...]=OneImg
  # ------------------------ Switch Small Window --------------------------------
  if FType=="Rect"||FType=="Cross"
    SWin=ones(MSize...)# Small Window
  elseif FType=="Triangle"
    SWin=NewLowFilter(FType,ones(MSize...))+1 #No zero elements are wanted
  elseif FType==""
  end
  # ========================== Making the NewImg ==============================
  NewImg=zeros(size(OneImg))
  BNew=zeros(size(BImg))
  for i=1:NDim
    for pp=VRange[i]
      ColV=ColArray(NDim,i,pp)
      for kk=-(MSize[i]-1)/2:(MSize[i]-1)/2
        if FType=="Cross"
          PosSWin=(MSize-1)/2+1# Centre every direction
          kk=convert(Int64,kk)
          PosSWin[i]=(MSize[i]-1)/2+1-kk
          FacSWin=SWin[convert(Vector{Int64},PosSWin)...]
          BNew[ColV...]=BNew[ColV...]+BImg[ColArray(NDim,i,pp+kk)...]*FacSWin
        else
          for xx=-((size(SWin,i)-1)/2):((size(SWin,i)-1)/2)
            xx=convert(Int64,xx)
            PosSWin=(MSize-1)/2+1
            kk=convert(Int64,kk)
            PosSWin[i]=(MSize[i]-1)/2+1-kk
            FacSWin=SWin[convert(Vector{Int64},PosSWin)...]
            ColV2=ColArray(NDim,i,pp+kk)
            BVec=BImg[ColV2...]
            if xx>=1
              for zz=1:xx
                LastE=pop!(BVec)
                BVec=unshift!(BVec,LastE)
              end
            end
            if xx<=-1
              for zz=xx:-1
                FirstE=shift!(BVec)
                BVec=push!(BVec,FirstE)
              end
            end
            BNew[ColV...]=BNew[ColV...]+BVec*FacSWin
          end
        end
      end
    end
  end
  NewImg=BNew[VRange...]
  NewImg=NewImg*sum(abs(OneImg))/sum(abs(NewImg)) #keeping the Energy
return NewImg
end
############################################################################

# A serie of Filters design for a LIF combustion application,
# more of an example than anything else

function FilterS1(Img::Array,Noise::Array,X4Param=20::Int64)
  X=Img-Noise
  F1=NewLowFilter("Rect",X,[1024,1340]);
  Fi=F1
  #matshow(X)
  X1=MeanFilter("Cross",X,[7,1]) # 7 = lambda sinus (1/freq spaciale)
  X2=real(ifft(real(fft(X1)).*F1+imag(fft(X1))*im.*Fi))
  #matshow(X2)
  X3=MeanFilter("Triangle",X2,[5,5])
  #matshow(X3)

  X4=X3
  X4[X4.<mean(X4)/X4Param]=0
  #matshow(X4)
  #title("Suite Filtre 1")
return X3
end
#############################################################################
end # module
