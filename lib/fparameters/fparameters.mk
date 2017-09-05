##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=fparameters
ConfigurationName      :=Debug
WorkspacePath          := "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori"
ProjectPath            := "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/lib/fparameters"
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Eduardo Mario Gutiérrez
Date                   :=30/08/17
CodeLitePath           :="/home/eduardog/.codelite"
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
ObjectsFileList        :="fparameters.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). $(IncludeSwitch).. $(IncludeSwitch)../../external/include 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). $(LibraryPathSwitch)../../lib/ $(LibraryPathSwitch)../../external/include/ 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++  -std=c++11
CC       := /usr/bin/gcc
CXXFLAGS :=  -g $(Preprocessors)
CFLAGS   :=  -g $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/parameters.cpp$(ObjectSuffix) $(IntermediateDirectory)/ParamSpace.cpp$(ObjectSuffix) $(IntermediateDirectory)/ParamSpaceValues.cpp$(ObjectSuffix) $(IntermediateDirectory)/SpaceCoord.cpp$(ObjectSuffix) $(IntermediateDirectory)/SpaceIterator.cpp$(ObjectSuffix) $(IntermediateDirectory)/Dimension.cpp$(ObjectSuffix) $(IntermediateDirectory)/DimensionIterator.cpp$(ObjectSuffix) 



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
	@$(MakeDirCommand) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori/.build-debug"
	@echo rebuilt > "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori/.build-debug/fparameters"

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


./Debug:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/parameters.cpp$(ObjectSuffix): parameters.cpp $(IntermediateDirectory)/parameters.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/lib/fparameters/parameters.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/parameters.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/parameters.cpp$(DependSuffix): parameters.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/parameters.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/parameters.cpp$(DependSuffix) -MM "parameters.cpp"

$(IntermediateDirectory)/parameters.cpp$(PreprocessSuffix): parameters.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/parameters.cpp$(PreprocessSuffix) "parameters.cpp"

$(IntermediateDirectory)/ParamSpace.cpp$(ObjectSuffix): ParamSpace.cpp $(IntermediateDirectory)/ParamSpace.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/lib/fparameters/ParamSpace.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/ParamSpace.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/ParamSpace.cpp$(DependSuffix): ParamSpace.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/ParamSpace.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/ParamSpace.cpp$(DependSuffix) -MM "ParamSpace.cpp"

$(IntermediateDirectory)/ParamSpace.cpp$(PreprocessSuffix): ParamSpace.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/ParamSpace.cpp$(PreprocessSuffix) "ParamSpace.cpp"

$(IntermediateDirectory)/ParamSpaceValues.cpp$(ObjectSuffix): ParamSpaceValues.cpp $(IntermediateDirectory)/ParamSpaceValues.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/lib/fparameters/ParamSpaceValues.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/ParamSpaceValues.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/ParamSpaceValues.cpp$(DependSuffix): ParamSpaceValues.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/ParamSpaceValues.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/ParamSpaceValues.cpp$(DependSuffix) -MM "ParamSpaceValues.cpp"

$(IntermediateDirectory)/ParamSpaceValues.cpp$(PreprocessSuffix): ParamSpaceValues.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/ParamSpaceValues.cpp$(PreprocessSuffix) "ParamSpaceValues.cpp"

$(IntermediateDirectory)/SpaceCoord.cpp$(ObjectSuffix): SpaceCoord.cpp $(IntermediateDirectory)/SpaceCoord.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/lib/fparameters/SpaceCoord.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/SpaceCoord.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/SpaceCoord.cpp$(DependSuffix): SpaceCoord.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/SpaceCoord.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/SpaceCoord.cpp$(DependSuffix) -MM "SpaceCoord.cpp"

$(IntermediateDirectory)/SpaceCoord.cpp$(PreprocessSuffix): SpaceCoord.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/SpaceCoord.cpp$(PreprocessSuffix) "SpaceCoord.cpp"

$(IntermediateDirectory)/SpaceIterator.cpp$(ObjectSuffix): SpaceIterator.cpp $(IntermediateDirectory)/SpaceIterator.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/lib/fparameters/SpaceIterator.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/SpaceIterator.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/SpaceIterator.cpp$(DependSuffix): SpaceIterator.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/SpaceIterator.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/SpaceIterator.cpp$(DependSuffix) -MM "SpaceIterator.cpp"

$(IntermediateDirectory)/SpaceIterator.cpp$(PreprocessSuffix): SpaceIterator.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/SpaceIterator.cpp$(PreprocessSuffix) "SpaceIterator.cpp"

$(IntermediateDirectory)/Dimension.cpp$(ObjectSuffix): Dimension.cpp $(IntermediateDirectory)/Dimension.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/lib/fparameters/Dimension.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Dimension.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Dimension.cpp$(DependSuffix): Dimension.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/Dimension.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/Dimension.cpp$(DependSuffix) -MM "Dimension.cpp"

$(IntermediateDirectory)/Dimension.cpp$(PreprocessSuffix): Dimension.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Dimension.cpp$(PreprocessSuffix) "Dimension.cpp"

$(IntermediateDirectory)/DimensionIterator.cpp$(ObjectSuffix): DimensionIterator.cpp $(IntermediateDirectory)/DimensionIterator.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/lib/fparameters/DimensionIterator.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/DimensionIterator.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/DimensionIterator.cpp$(DependSuffix): DimensionIterator.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/DimensionIterator.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/DimensionIterator.cpp$(DependSuffix) -MM "DimensionIterator.cpp"

$(IntermediateDirectory)/DimensionIterator.cpp$(PreprocessSuffix): DimensionIterator.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/DimensionIterator.cpp$(PreprocessSuffix) "DimensionIterator.cpp"


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


