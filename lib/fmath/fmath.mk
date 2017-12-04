##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=fmath
ConfigurationName      :=Debug
WorkspacePath          :=/home/flor/projects/github/radproc_linux/ionTori
ProjectPath            :=/home/flor/projects/github/radproc_linux/lib/fmath
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=flor
Date                   :=04/12/17
CodeLitePath           :=/home/flor/.codelite
LinkerName             :=/usr/bin/g++
SharedObjectLinkerName :=/usr/bin/g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.i
DebugSwitch            :=-g 
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/lib$(ProjectName).a
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="fmath.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++ -std=c++11
CC       := /usr/bin/gcc
CXXFLAGS :=  -g $(Preprocessors)
CFLAGS   :=  -g $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/elimiGaussiana.cpp$(ObjectSuffix) $(IntermediateDirectory)/gaussianScaledPivoting.cpp$(ObjectSuffix) $(IntermediateDirectory)/interpolation.cpp$(ObjectSuffix) $(IntermediateDirectory)/matrixInit.cpp$(ObjectSuffix) $(IntermediateDirectory)/RungeKutta.cpp$(ObjectSuffix) $(IntermediateDirectory)/bisection.cpp$(ObjectSuffix) $(IntermediateDirectory)/mathFunctions.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(IntermediateDirectory) $(OutputFile)

$(OutputFile): $(Objects)
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(AR) $(ArchiveOutputSwitch)$(OutputFile) @$(ObjectsFileList) $(ArLibs)
	@$(MakeDirCommand) "/home/flor/projects/github/radproc_linux/ionTori/.build-debug"
	@echo rebuilt > "/home/flor/projects/github/radproc_linux/ionTori/.build-debug/fmath"

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


./Debug:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/elimiGaussiana.cpp$(ObjectSuffix): elimiGaussiana.cpp $(IntermediateDirectory)/elimiGaussiana.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/flor/projects/github/radproc_linux/lib/fmath/elimiGaussiana.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/elimiGaussiana.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/elimiGaussiana.cpp$(DependSuffix): elimiGaussiana.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/elimiGaussiana.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/elimiGaussiana.cpp$(DependSuffix) -MM elimiGaussiana.cpp

$(IntermediateDirectory)/elimiGaussiana.cpp$(PreprocessSuffix): elimiGaussiana.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/elimiGaussiana.cpp$(PreprocessSuffix) elimiGaussiana.cpp

$(IntermediateDirectory)/gaussianScaledPivoting.cpp$(ObjectSuffix): gaussianScaledPivoting.cpp $(IntermediateDirectory)/gaussianScaledPivoting.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/flor/projects/github/radproc_linux/lib/fmath/gaussianScaledPivoting.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/gaussianScaledPivoting.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/gaussianScaledPivoting.cpp$(DependSuffix): gaussianScaledPivoting.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/gaussianScaledPivoting.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/gaussianScaledPivoting.cpp$(DependSuffix) -MM gaussianScaledPivoting.cpp

$(IntermediateDirectory)/gaussianScaledPivoting.cpp$(PreprocessSuffix): gaussianScaledPivoting.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/gaussianScaledPivoting.cpp$(PreprocessSuffix) gaussianScaledPivoting.cpp

$(IntermediateDirectory)/interpolation.cpp$(ObjectSuffix): interpolation.cpp $(IntermediateDirectory)/interpolation.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/flor/projects/github/radproc_linux/lib/fmath/interpolation.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/interpolation.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/interpolation.cpp$(DependSuffix): interpolation.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/interpolation.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/interpolation.cpp$(DependSuffix) -MM interpolation.cpp

$(IntermediateDirectory)/interpolation.cpp$(PreprocessSuffix): interpolation.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/interpolation.cpp$(PreprocessSuffix) interpolation.cpp

$(IntermediateDirectory)/matrixInit.cpp$(ObjectSuffix): matrixInit.cpp $(IntermediateDirectory)/matrixInit.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/flor/projects/github/radproc_linux/lib/fmath/matrixInit.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/matrixInit.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/matrixInit.cpp$(DependSuffix): matrixInit.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/matrixInit.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/matrixInit.cpp$(DependSuffix) -MM matrixInit.cpp

$(IntermediateDirectory)/matrixInit.cpp$(PreprocessSuffix): matrixInit.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/matrixInit.cpp$(PreprocessSuffix) matrixInit.cpp

$(IntermediateDirectory)/RungeKutta.cpp$(ObjectSuffix): RungeKutta.cpp $(IntermediateDirectory)/RungeKutta.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/flor/projects/github/radproc_linux/lib/fmath/RungeKutta.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/RungeKutta.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/RungeKutta.cpp$(DependSuffix): RungeKutta.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/RungeKutta.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/RungeKutta.cpp$(DependSuffix) -MM RungeKutta.cpp

$(IntermediateDirectory)/RungeKutta.cpp$(PreprocessSuffix): RungeKutta.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/RungeKutta.cpp$(PreprocessSuffix) RungeKutta.cpp

$(IntermediateDirectory)/bisection.cpp$(ObjectSuffix): bisection.cpp $(IntermediateDirectory)/bisection.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/flor/projects/github/radproc_linux/lib/fmath/bisection.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/bisection.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/bisection.cpp$(DependSuffix): bisection.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/bisection.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/bisection.cpp$(DependSuffix) -MM bisection.cpp

$(IntermediateDirectory)/bisection.cpp$(PreprocessSuffix): bisection.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/bisection.cpp$(PreprocessSuffix) bisection.cpp

$(IntermediateDirectory)/mathFunctions.cpp$(ObjectSuffix): mathFunctions.cpp $(IntermediateDirectory)/mathFunctions.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/flor/projects/github/radproc_linux/lib/fmath/mathFunctions.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/mathFunctions.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/mathFunctions.cpp$(DependSuffix): mathFunctions.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/mathFunctions.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/mathFunctions.cpp$(DependSuffix) -MM mathFunctions.cpp

$(IntermediateDirectory)/mathFunctions.cpp$(PreprocessSuffix): mathFunctions.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/mathFunctions.cpp$(PreprocessSuffix) mathFunctions.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


