##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=ionTori
ConfigurationName      :=Debug
WorkspacePath          := "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori"
ProjectPath            := "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori"
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Eduardo Mario Gutiérrez
Date                   :=09/09/17
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
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
Preprocessors          :=$(PreprocessorSwitch)BOOST_SYSTEM_NO_DEPRECATED $(PreprocessorSwitch)BOOST_NO_CXX11_SCOPED_ENUMS 
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="ionTori.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). $(IncludeSwitch)../lib $(IncludeSwitch)../external/include 
IncludePCH             := 
RcIncludePath          := 
Libs                   := $(LibrarySwitch)boost_filesystem $(LibrarySwitch)boost_system $(LibrarySwitch)fmath $(LibrarySwitch)fparameters $(LibrarySwitch)fparticle $(LibrarySwitch)flosses $(LibrarySwitch)fluminosities $(LibrarySwitch)inout 
ArLibs                 :=  "boost_filesystem" "boost_system" "fmath" "fparameters" "fparticle" "flosses" "fluminosities" "inout" 
LibPath                := $(LibraryPathSwitch). $(LibraryPathSwitch)../external/lib $(LibraryPathSwitch)../lib/fmath/Debug $(LibraryPathSwitch)../lib/fparameters/Debug $(LibraryPathSwitch)../lib/fparticle/Debug $(LibraryPathSwitch)../lib/flosses/Debug $(LibraryPathSwitch)../lib/fluminosities/Debug $(LibraryPathSwitch)../lib/inout/Debug 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++  -std=c++11
CC       := /usr/bin/gcc
CXXFLAGS :=  -g -O0 -Wall $(Preprocessors)
CFLAGS   :=  -g -O0 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IntermediateDirectory)/modelParameters.cpp$(ObjectSuffix) $(IntermediateDirectory)/radiativeLosses.cpp$(ObjectSuffix) $(IntermediateDirectory)/State.cpp$(ObjectSuffix) $(IntermediateDirectory)/write.cpp$(ObjectSuffix) $(IntermediateDirectory)/functions.cpp$(ObjectSuffix) $(IntermediateDirectory)/toroParam.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d ".build-debug/fmath" ".build-debug/fparameters" ".build-debug/fparticle" ".build-debug/inout" ".build-debug/flosses" ".build-debug/fluminosities" $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

".build-debug/fmath":
	@$(MakeDirCommand) ".build-debug"
	@echo stam > ".build-debug/fmath"


".build-debug/fparameters":
	@$(MakeDirCommand) ".build-debug"
	@echo stam > ".build-debug/fparameters"


".build-debug/fparticle":
	@$(MakeDirCommand) ".build-debug"
	@echo stam > ".build-debug/fparticle"


".build-debug/inout":
	@$(MakeDirCommand) ".build-debug"
	@echo stam > ".build-debug/inout"


".build-debug/flosses":
	@$(MakeDirCommand) ".build-debug"
	@echo stam > ".build-debug/flosses"


".build-debug/fluminosities":
	@$(MakeDirCommand) ".build-debug"
	@echo stam > ".build-debug/fluminosities"




MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


$(IntermediateDirectory)/.d:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/main.cpp$(ObjectSuffix): main.cpp $(IntermediateDirectory)/main.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori/main.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/main.cpp$(DependSuffix): main.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/main.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/main.cpp$(DependSuffix) -MM "main.cpp"

$(IntermediateDirectory)/main.cpp$(PreprocessSuffix): main.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/main.cpp$(PreprocessSuffix) "main.cpp"

$(IntermediateDirectory)/modelParameters.cpp$(ObjectSuffix): modelParameters.cpp $(IntermediateDirectory)/modelParameters.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori/modelParameters.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/modelParameters.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/modelParameters.cpp$(DependSuffix): modelParameters.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/modelParameters.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/modelParameters.cpp$(DependSuffix) -MM "modelParameters.cpp"

$(IntermediateDirectory)/modelParameters.cpp$(PreprocessSuffix): modelParameters.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/modelParameters.cpp$(PreprocessSuffix) "modelParameters.cpp"

$(IntermediateDirectory)/radiativeLosses.cpp$(ObjectSuffix): radiativeLosses.cpp $(IntermediateDirectory)/radiativeLosses.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori/radiativeLosses.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/radiativeLosses.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/radiativeLosses.cpp$(DependSuffix): radiativeLosses.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/radiativeLosses.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/radiativeLosses.cpp$(DependSuffix) -MM "radiativeLosses.cpp"

$(IntermediateDirectory)/radiativeLosses.cpp$(PreprocessSuffix): radiativeLosses.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/radiativeLosses.cpp$(PreprocessSuffix) "radiativeLosses.cpp"

$(IntermediateDirectory)/State.cpp$(ObjectSuffix): State.cpp $(IntermediateDirectory)/State.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori/State.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/State.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/State.cpp$(DependSuffix): State.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/State.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/State.cpp$(DependSuffix) -MM "State.cpp"

$(IntermediateDirectory)/State.cpp$(PreprocessSuffix): State.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/State.cpp$(PreprocessSuffix) "State.cpp"

$(IntermediateDirectory)/write.cpp$(ObjectSuffix): write.cpp $(IntermediateDirectory)/write.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori/write.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/write.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/write.cpp$(DependSuffix): write.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/write.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/write.cpp$(DependSuffix) -MM "write.cpp"

$(IntermediateDirectory)/write.cpp$(PreprocessSuffix): write.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/write.cpp$(PreprocessSuffix) "write.cpp"

$(IntermediateDirectory)/functions.cpp$(ObjectSuffix): functions.cpp $(IntermediateDirectory)/functions.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori/functions.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/functions.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/functions.cpp$(DependSuffix): functions.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/functions.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/functions.cpp$(DependSuffix) -MM "functions.cpp"

$(IntermediateDirectory)/functions.cpp$(PreprocessSuffix): functions.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/functions.cpp$(PreprocessSuffix) "functions.cpp"

$(IntermediateDirectory)/toroParam.cpp$(ObjectSuffix): toroParam.cpp $(IntermediateDirectory)/toroParam.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori/toroParam.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/toroParam.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/toroParam.cpp$(DependSuffix): toroParam.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/toroParam.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/toroParam.cpp$(DependSuffix) -MM "toroParam.cpp"

$(IntermediateDirectory)/toroParam.cpp$(PreprocessSuffix): toroParam.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/toroParam.cpp$(PreprocessSuffix) "toroParam.cpp"


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


