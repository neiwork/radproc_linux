##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=fparticle
ConfigurationName      :=Debug
WorkspacePath          := "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori"
ProjectPath            := "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/lib/fparticle"
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
ObjectsFileList        :="fparticle.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). $(IncludeSwitch).. $(IncludeSwitch)../../external/include 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). $(LibraryPathSwitch)../../lib/ 

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
Objects0=$(IntermediateDirectory)/Particle.cpp$(ObjectSuffix) 



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
	@echo rebuilt > "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/ionTori/.build-debug/fparticle"

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


./Debug:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/Particle.cpp$(ObjectSuffix): Particle.cpp $(IntermediateDirectory)/Particle.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/lib/fparticle/Particle.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Particle.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Particle.cpp$(DependSuffix): Particle.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/Particle.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/Particle.cpp$(DependSuffix) -MM "Particle.cpp"

$(IntermediateDirectory)/Particle.cpp$(PreprocessSuffix): Particle.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Particle.cpp$(PreprocessSuffix) "Particle.cpp"


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


