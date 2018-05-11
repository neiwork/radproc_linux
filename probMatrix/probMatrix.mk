##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Release
ProjectName            :=probMatrix
ConfigurationName      :=Release
WorkspacePath          :=/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/probMatrix
ProjectPath            :=/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/probMatrix
IntermediateDirectory  :=./Release
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Eduardo Mario Gutiérrez
Date                   :=11/05/18
CodeLitePath           :=/home/eduardog/.codelite
LinkerName             :=/usr/bin/g++-4.8
SharedObjectLinkerName :=/usr/bin/g++-4.8 -shared -fPIC
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
Preprocessors          :=$(PreprocessorSwitch)NDEBUG 
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="probMatrix.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  -lgsl -lgslcblas -lm
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch)/usr/include 
IncludePCH             := 
RcIncludePath          := 
Libs                   := $(LibrarySwitch)gsl 
ArLibs                 :=  "gsl" 
LibPath                := $(LibraryPathSwitch). $(LibraryPathSwitch)/home/eduardog/gsl/lib 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++-4.8
CC       := /usr/bin/gcc-4.8
CXXFLAGS :=  -O2 -Wall $(Preprocessors)
CFLAGS   :=  -O1 -O -O3 -O2 -std=c99 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/main.c$(ObjectSuffix) $(IntermediateDirectory)/nrutil.c$(ObjectSuffix) $(IntermediateDirectory)/metric.c$(ObjectSuffix) $(IntermediateDirectory)/torusParameters.c$(ObjectSuffix) $(IntermediateDirectory)/auxFunctions.c$(ObjectSuffix) $(IntermediateDirectory)/torusSampling.c$(ObjectSuffix) $(IntermediateDirectory)/bisection.c$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@test -d ./Release || $(MakeDirCommand) ./Release


$(IntermediateDirectory)/.d:
	@test -d ./Release || $(MakeDirCommand) ./Release

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/main.c$(ObjectSuffix): main.c $(IntermediateDirectory)/main.c$(DependSuffix)
	$(CC) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/probMatrix/main.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/main.c$(DependSuffix): main.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/main.c$(ObjectSuffix) -MF$(IntermediateDirectory)/main.c$(DependSuffix) -MM main.c

$(IntermediateDirectory)/main.c$(PreprocessSuffix): main.c
	$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/main.c$(PreprocessSuffix) main.c

$(IntermediateDirectory)/nrutil.c$(ObjectSuffix): nrutil.c $(IntermediateDirectory)/nrutil.c$(DependSuffix)
	$(CC) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/probMatrix/nrutil.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/nrutil.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/nrutil.c$(DependSuffix): nrutil.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/nrutil.c$(ObjectSuffix) -MF$(IntermediateDirectory)/nrutil.c$(DependSuffix) -MM nrutil.c

$(IntermediateDirectory)/nrutil.c$(PreprocessSuffix): nrutil.c
	$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/nrutil.c$(PreprocessSuffix) nrutil.c

$(IntermediateDirectory)/metric.c$(ObjectSuffix): metric.c $(IntermediateDirectory)/metric.c$(DependSuffix)
	$(CC) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/probMatrix/metric.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/metric.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/metric.c$(DependSuffix): metric.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/metric.c$(ObjectSuffix) -MF$(IntermediateDirectory)/metric.c$(DependSuffix) -MM metric.c

$(IntermediateDirectory)/metric.c$(PreprocessSuffix): metric.c
	$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/metric.c$(PreprocessSuffix) metric.c

$(IntermediateDirectory)/torusParameters.c$(ObjectSuffix): torusParameters.c $(IntermediateDirectory)/torusParameters.c$(DependSuffix)
	$(CC) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/probMatrix/torusParameters.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/torusParameters.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/torusParameters.c$(DependSuffix): torusParameters.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/torusParameters.c$(ObjectSuffix) -MF$(IntermediateDirectory)/torusParameters.c$(DependSuffix) -MM torusParameters.c

$(IntermediateDirectory)/torusParameters.c$(PreprocessSuffix): torusParameters.c
	$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/torusParameters.c$(PreprocessSuffix) torusParameters.c

$(IntermediateDirectory)/auxFunctions.c$(ObjectSuffix): auxFunctions.c $(IntermediateDirectory)/auxFunctions.c$(DependSuffix)
	$(CC) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/probMatrix/auxFunctions.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/auxFunctions.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/auxFunctions.c$(DependSuffix): auxFunctions.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/auxFunctions.c$(ObjectSuffix) -MF$(IntermediateDirectory)/auxFunctions.c$(DependSuffix) -MM auxFunctions.c

$(IntermediateDirectory)/auxFunctions.c$(PreprocessSuffix): auxFunctions.c
	$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/auxFunctions.c$(PreprocessSuffix) auxFunctions.c

$(IntermediateDirectory)/torusSampling.c$(ObjectSuffix): torusSampling.c $(IntermediateDirectory)/torusSampling.c$(DependSuffix)
	$(CC) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/probMatrix/torusSampling.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/torusSampling.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/torusSampling.c$(DependSuffix): torusSampling.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/torusSampling.c$(ObjectSuffix) -MF$(IntermediateDirectory)/torusSampling.c$(DependSuffix) -MM torusSampling.c

$(IntermediateDirectory)/torusSampling.c$(PreprocessSuffix): torusSampling.c
	$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/torusSampling.c$(PreprocessSuffix) torusSampling.c

$(IntermediateDirectory)/bisection.c$(ObjectSuffix): bisection.c $(IntermediateDirectory)/bisection.c$(DependSuffix)
	$(CC) $(SourceSwitch) "/home/eduardog/FACULTAD/DOCTORADO/Tori/neiwork/radproc_linux/probMatrix/bisection.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/bisection.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/bisection.c$(DependSuffix): bisection.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/bisection.c$(ObjectSuffix) -MF$(IntermediateDirectory)/bisection.c$(DependSuffix) -MM bisection.c

$(IntermediateDirectory)/bisection.c$(PreprocessSuffix): bisection.c
	$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/bisection.c$(PreprocessSuffix) bisection.c


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Release/


