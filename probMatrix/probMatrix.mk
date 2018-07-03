##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=probMatrix
ConfigurationName      :=Debug
WorkspacePath          :=/home/luciano/Doctorado/Accretion/radproc_linux/probMatrix
ProjectPath            :=/home/luciano/Doctorado/Accretion/radproc_linux/probMatrix
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Luciano
Date                   :=29/06/18
CodeLitePath           :=/home/luciano/.codelite
LinkerName             :=/usr/bin/x86_64-linux-gnu-g++
SharedObjectLinkerName :=/usr/bin/x86_64-linux-gnu-g++ -shared -fPIC
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
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="probMatrix.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  -lgsl -lgslcblas -lm
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch)/home/luciano/gsl/include 
IncludePCH             := 
RcIncludePath          := 
Libs                   := $(LibrarySwitch)gsl $(LibrarySwitch)gsl 
ArLibs                 :=  "gsl" "gsl" 
LibPath                := $(LibraryPathSwitch)/home/luciano/gsl/lib $(LibraryPathSwitch)/home/luciano/gsl/lib 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/x86_64-linux-gnu-ar rcu
CXX      := /usr/bin/x86_64-linux-gnu-g++
CC       := /usr/bin/x86_64-linux-gnu-gcc
CXXFLAGS :=  -g -O0 -std=c++11 -Wall $(Preprocessors)
CFLAGS   := -fPIC -g -O0 -std=c99 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/x86_64-linux-gnu-as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/torusSampling2.c$(ObjectSuffix) $(IntermediateDirectory)/interpol.c$(ObjectSuffix) $(IntermediateDirectory)/functions.c$(ObjectSuffix) $(IntermediateDirectory)/nrutil.c$(ObjectSuffix) $(IntermediateDirectory)/main.c$(ObjectSuffix) 



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
	@test -d ./Debug || $(MakeDirCommand) ./Debug


$(IntermediateDirectory)/.d:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/torusSampling2.c$(ObjectSuffix): torusSampling2.c $(IntermediateDirectory)/torusSampling2.c$(DependSuffix)
	$(CC) $(SourceSwitch) "/home/luciano/Doctorado/Accretion/radproc_linux/probMatrix/torusSampling2.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/torusSampling2.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/torusSampling2.c$(DependSuffix): torusSampling2.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/torusSampling2.c$(ObjectSuffix) -MF$(IntermediateDirectory)/torusSampling2.c$(DependSuffix) -MM torusSampling2.c

$(IntermediateDirectory)/torusSampling2.c$(PreprocessSuffix): torusSampling2.c
	$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/torusSampling2.c$(PreprocessSuffix) torusSampling2.c

$(IntermediateDirectory)/interpol.c$(ObjectSuffix): interpol.c $(IntermediateDirectory)/interpol.c$(DependSuffix)
	$(CC) $(SourceSwitch) "/home/luciano/Doctorado/Accretion/radproc_linux/probMatrix/interpol.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/interpol.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/interpol.c$(DependSuffix): interpol.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/interpol.c$(ObjectSuffix) -MF$(IntermediateDirectory)/interpol.c$(DependSuffix) -MM interpol.c

$(IntermediateDirectory)/interpol.c$(PreprocessSuffix): interpol.c
	$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/interpol.c$(PreprocessSuffix) interpol.c

$(IntermediateDirectory)/functions.c$(ObjectSuffix): functions.c $(IntermediateDirectory)/functions.c$(DependSuffix)
	$(CC) $(SourceSwitch) "/home/luciano/Doctorado/Accretion/radproc_linux/probMatrix/functions.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/functions.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/functions.c$(DependSuffix): functions.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/functions.c$(ObjectSuffix) -MF$(IntermediateDirectory)/functions.c$(DependSuffix) -MM functions.c

$(IntermediateDirectory)/functions.c$(PreprocessSuffix): functions.c
	$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/functions.c$(PreprocessSuffix) functions.c

$(IntermediateDirectory)/nrutil.c$(ObjectSuffix): nrutil.c $(IntermediateDirectory)/nrutil.c$(DependSuffix)
	$(CC) $(SourceSwitch) "/home/luciano/Doctorado/Accretion/radproc_linux/probMatrix/nrutil.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/nrutil.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/nrutil.c$(DependSuffix): nrutil.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/nrutil.c$(ObjectSuffix) -MF$(IntermediateDirectory)/nrutil.c$(DependSuffix) -MM nrutil.c

$(IntermediateDirectory)/nrutil.c$(PreprocessSuffix): nrutil.c
	$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/nrutil.c$(PreprocessSuffix) nrutil.c

$(IntermediateDirectory)/main.c$(ObjectSuffix): main.c $(IntermediateDirectory)/main.c$(DependSuffix)
	$(CC) $(SourceSwitch) "/home/luciano/Doctorado/Accretion/radproc_linux/probMatrix/main.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/main.c$(DependSuffix): main.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/main.c$(ObjectSuffix) -MF$(IntermediateDirectory)/main.c$(DependSuffix) -MM main.c

$(IntermediateDirectory)/main.c$(PreprocessSuffix): main.c
	$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/main.c$(PreprocessSuffix) main.c


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


