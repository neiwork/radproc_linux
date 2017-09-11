##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=fluminosities
ConfigurationName      :=Debug
WorkspacePath          := "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori"
ProjectPath            := "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/lib/fluminosities"
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Eduardo Mario Gutiérrez
Date                   :=10/09/17
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
ObjectsFileList        :="fluminosities.txt"
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
Objects0=$(IntermediateDirectory)/luminosityIC.cpp$(ObjectSuffix) $(IntermediateDirectory)/luminositySynchrotron.cpp$(ObjectSuffix) $(IntermediateDirectory)/opticalDepthSSA.cpp$(ObjectSuffix) 



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
	@echo rebuilt > "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori/.build-debug/fluminosities"

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


./Debug:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/luminosityIC.cpp$(ObjectSuffix): luminosityIC.cpp $(IntermediateDirectory)/luminosityIC.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/lib/fluminosities/luminosityIC.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/luminosityIC.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/luminosityIC.cpp$(DependSuffix): luminosityIC.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/luminosityIC.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/luminosityIC.cpp$(DependSuffix) -MM "luminosityIC.cpp"

$(IntermediateDirectory)/luminosityIC.cpp$(PreprocessSuffix): luminosityIC.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/luminosityIC.cpp$(PreprocessSuffix) "luminosityIC.cpp"

$(IntermediateDirectory)/luminositySynchrotron.cpp$(ObjectSuffix): luminositySynchrotron.cpp $(IntermediateDirectory)/luminositySynchrotron.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/lib/fluminosities/luminositySynchrotron.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/luminositySynchrotron.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/luminositySynchrotron.cpp$(DependSuffix): luminositySynchrotron.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/luminositySynchrotron.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/luminositySynchrotron.cpp$(DependSuffix) -MM "luminositySynchrotron.cpp"

$(IntermediateDirectory)/luminositySynchrotron.cpp$(PreprocessSuffix): luminositySynchrotron.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/luminositySynchrotron.cpp$(PreprocessSuffix) "luminositySynchrotron.cpp"

$(IntermediateDirectory)/opticalDepthSSA.cpp$(ObjectSuffix): opticalDepthSSA.cpp $(IntermediateDirectory)/opticalDepthSSA.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/lib/fluminosities/opticalDepthSSA.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/opticalDepthSSA.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/opticalDepthSSA.cpp$(DependSuffix): opticalDepthSSA.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/opticalDepthSSA.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/opticalDepthSSA.cpp$(DependSuffix) -MM "opticalDepthSSA.cpp"

$(IntermediateDirectory)/opticalDepthSSA.cpp$(PreprocessSuffix): opticalDepthSSA.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/opticalDepthSSA.cpp$(PreprocessSuffix) "opticalDepthSSA.cpp"


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


