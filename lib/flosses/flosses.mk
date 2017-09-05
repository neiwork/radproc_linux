##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=flosses
ConfigurationName      :=Debug
WorkspacePath          :=/home/flor/projects/github/radproc_linux/ionTori
ProjectPath            :=/home/flor/projects/github/radproc_linux/lib/flosses
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=flor
Date                   :=05/09/17
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
ObjectsFileList        :="flosses.txt"
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
Objects0=$(IntermediateDirectory)/crossSectionInel.cpp$(ObjectSuffix) $(IntermediateDirectory)/lossesBrem.cpp$(ObjectSuffix) $(IntermediateDirectory)/lossesIC.cpp$(ObjectSuffix) $(IntermediateDirectory)/lossesPhotoHadronic.cpp$(ObjectSuffix) $(IntermediateDirectory)/lossesSyn.cpp$(ObjectSuffix) $(IntermediateDirectory)/nonThermalLosses.cpp$(ObjectSuffix) $(IntermediateDirectory)/lossesHadronics.cpp$(ObjectSuffix) 



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
	@echo rebuilt > "/home/flor/projects/github/radproc_linux/ionTori/.build-debug/flosses"

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


./Debug:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/crossSectionInel.cpp$(ObjectSuffix): crossSectionInel.cpp $(IntermediateDirectory)/crossSectionInel.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/flor/projects/github/radproc_linux/lib/flosses/crossSectionInel.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/crossSectionInel.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/crossSectionInel.cpp$(DependSuffix): crossSectionInel.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/crossSectionInel.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/crossSectionInel.cpp$(DependSuffix) -MM crossSectionInel.cpp

$(IntermediateDirectory)/crossSectionInel.cpp$(PreprocessSuffix): crossSectionInel.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/crossSectionInel.cpp$(PreprocessSuffix) crossSectionInel.cpp

$(IntermediateDirectory)/lossesBrem.cpp$(ObjectSuffix): lossesBrem.cpp $(IntermediateDirectory)/lossesBrem.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/flor/projects/github/radproc_linux/lib/flosses/lossesBrem.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/lossesBrem.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/lossesBrem.cpp$(DependSuffix): lossesBrem.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/lossesBrem.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/lossesBrem.cpp$(DependSuffix) -MM lossesBrem.cpp

$(IntermediateDirectory)/lossesBrem.cpp$(PreprocessSuffix): lossesBrem.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/lossesBrem.cpp$(PreprocessSuffix) lossesBrem.cpp

$(IntermediateDirectory)/lossesIC.cpp$(ObjectSuffix): lossesIC.cpp $(IntermediateDirectory)/lossesIC.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/flor/projects/github/radproc_linux/lib/flosses/lossesIC.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/lossesIC.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/lossesIC.cpp$(DependSuffix): lossesIC.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/lossesIC.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/lossesIC.cpp$(DependSuffix) -MM lossesIC.cpp

$(IntermediateDirectory)/lossesIC.cpp$(PreprocessSuffix): lossesIC.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/lossesIC.cpp$(PreprocessSuffix) lossesIC.cpp

$(IntermediateDirectory)/lossesPhotoHadronic.cpp$(ObjectSuffix): lossesPhotoHadronic.cpp $(IntermediateDirectory)/lossesPhotoHadronic.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/flor/projects/github/radproc_linux/lib/flosses/lossesPhotoHadronic.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/lossesPhotoHadronic.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/lossesPhotoHadronic.cpp$(DependSuffix): lossesPhotoHadronic.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/lossesPhotoHadronic.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/lossesPhotoHadronic.cpp$(DependSuffix) -MM lossesPhotoHadronic.cpp

$(IntermediateDirectory)/lossesPhotoHadronic.cpp$(PreprocessSuffix): lossesPhotoHadronic.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/lossesPhotoHadronic.cpp$(PreprocessSuffix) lossesPhotoHadronic.cpp

$(IntermediateDirectory)/lossesSyn.cpp$(ObjectSuffix): lossesSyn.cpp $(IntermediateDirectory)/lossesSyn.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/flor/projects/github/radproc_linux/lib/flosses/lossesSyn.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/lossesSyn.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/lossesSyn.cpp$(DependSuffix): lossesSyn.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/lossesSyn.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/lossesSyn.cpp$(DependSuffix) -MM lossesSyn.cpp

$(IntermediateDirectory)/lossesSyn.cpp$(PreprocessSuffix): lossesSyn.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/lossesSyn.cpp$(PreprocessSuffix) lossesSyn.cpp

$(IntermediateDirectory)/nonThermalLosses.cpp$(ObjectSuffix): nonThermalLosses.cpp $(IntermediateDirectory)/nonThermalLosses.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/flor/projects/github/radproc_linux/lib/flosses/nonThermalLosses.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/nonThermalLosses.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/nonThermalLosses.cpp$(DependSuffix): nonThermalLosses.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/nonThermalLosses.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/nonThermalLosses.cpp$(DependSuffix) -MM nonThermalLosses.cpp

$(IntermediateDirectory)/nonThermalLosses.cpp$(PreprocessSuffix): nonThermalLosses.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/nonThermalLosses.cpp$(PreprocessSuffix) nonThermalLosses.cpp

$(IntermediateDirectory)/lossesHadronics.cpp$(ObjectSuffix): lossesHadronics.cpp $(IntermediateDirectory)/lossesHadronics.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/flor/projects/github/radproc_linux/lib/flosses/lossesHadronics.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/lossesHadronics.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/lossesHadronics.cpp$(DependSuffix): lossesHadronics.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/lossesHadronics.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/lossesHadronics.cpp$(DependSuffix) -MM lossesHadronics.cpp

$(IntermediateDirectory)/lossesHadronics.cpp$(PreprocessSuffix): lossesHadronics.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/lossesHadronics.cpp$(PreprocessSuffix) lossesHadronics.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


