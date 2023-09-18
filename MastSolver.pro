TEMPLATE = app
#CONFIG += console c++17 il_std ilo_window stl
CONFIG += console c++17
CONFIG -= app_bundle
QT += gui core charts
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
DEFINES += USE_TIMER USE_GENERATORWINDOW
DEFINES += ILOUSESTL
win32{
DEFINES += ILO_WINDOWS
}
DEFINES += IL_STD
unix{
INCLUDEPATH += /opt/ibm/ILOG/CPLEX_Studio128/cplex/include
INCLUDEPATH += /opt/ibm/ILOG/CPLEX_Studio128/concert/include
}

SOURCES += \
    GeneratorBase.cpp \
    GeneratorTerminal.cpp \
    GeneratorUtilities.cpp \
#    GeneratorWindow.cpp \
#    MastTimer.cpp \
#    Scenario.cpp \
#    SolverConstructor.cpp \
    Utilities.cpp \
    main.cpp \
    Solver.cpp \
    ReaderHelper.cpp \
    WriteHelper.cpp \
    SolverBase.cpp \
    TripsInfo.cpp\
    Hybrid.cpp \
    Heuristic.cpp

DISTFILES += \
    README.md

HEADERS += \
    GeneratorBase.h \
    GeneratorTerminal.h \
    GeneratorUtilities.h \
#    GeneratorWindow.h \
#    MastTimer.h \
#    Scenario.h \
    Solver.h \
    ReaderHelper.h \
#    SolverConstructor.h \
    Utilities.h \
    WriteHelper.h \
    SolverBase.h \
    TripsInfo.h\
    Macro.h\
    Hybrid.h \
    Heuristic.h
unix:!macx: LIBS += -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/lib/x86-64_linux/static_pic/ -lilocplex
unix:!macx: PRE_TARGETDEPS += /opt/ibm/ILOG/CPLEX_Studio128/cplex/lib/x86-64_linux/static_pic/libilocplex.a

unix:!macx: LIBS += -L/opt/ibm/ILOG/CPLEX_Studio128/concert/lib/x86-64_linux/static_pic/ -lconcert
unix:!macx: PRE_TARGETDEPS += /opt/ibm/ILOG/CPLEX_Studio128/concert/lib/x86-64_linux/static_pic/libconcert.a

unix:!macx: LIBS += -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/ -lcplex1280
unix:!macx: LIBS += -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/ -lcplex1280mpitransport
unix:!macx: LIBS += -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/ -lcplex1280mpiworker
unix:!macx: LIBS += -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/ -lcplex1280processtransport
unix:!macx: LIBS += -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/ -lcplex1280processworker
unix:!macx: LIBS += -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/ -lcplex1280remote
unix:!macx: LIBS += -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/ -lcplex1280remotejni
unix:!macx: LIBS += -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/ -lcplex1280tcpiptransport
unix:!macx: LIBS += -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/ -lcplex1280tcpipworker

win32{

CONFIG(release, debug|release){
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/ -lilocplex
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/ -lcplex1280
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/ -lcplex1280processtransport
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/ -lcplex1280processworker
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/ -lcplex1280remote
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/ -lcplex1280tcpiptransport
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/ -lcplex1280tcpipworker
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/ -linteractive
        win32:LIBS += -L$$(CPLEX_STUDIO_DIR128)/concert/lib/x64_windows_vs2017/stat_mda/ -lconcert
        }
CONFIG(debug, debug|release){
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/ -lilocplex
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/ -lcplex1280
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/ -lcplex1280processtransport
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/ -lcplex1280processworker
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/ -lcplex1280remote
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/ -lcplex1280tcpiptransport
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/ -lcplex1280tcpipworker
        win32:LIBS += -L$$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/ -linteractive
        win32:LIBS += -L$$(CPLEX_STUDIO_DIR128)/concert/lib/x64_windows_vs2017/stat_mdd/ -lconcert
        }
}

win32{
INCLUDEPATH += $$(CPLEX_DIR)/include
DEPENDPATH += $$(CPLEX_DIR)/include
INCLUDEPATH += $$(CPLEX_STUDIO_DIR128)/concert/include
DEPENDPATH += $$(CPLEX_STUDIO_DIR128)/concert/include
}

#win32{
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/ilocplex.lib
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/cplex1280.lib
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/cplex1280processtransport.lib
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/cplex1280processworker.lib
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/cplex1280remote.lib
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/cplex1280tcpiptransport.lib
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/cplex1280tcpipworker.lib
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/interactive.lib
#        PRE_TARGETDEPS += $$(CPLEX_STUDIO_DIR128)/concert/lib/x64_windows_vs2017/stat_mda/concert.lib

#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/ilocplex.lib
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/cplex1280.lib
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/cplex1280processtransport.lib
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/cplex1280processworker.lib
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/cplex1280remote.lib
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/cplex1280tcpiptransport.lib
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/cplex1280tcpipworker.lib
#        PRE_TARGETDEPS += $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mdd/interactive.lib
#        PRE_TARGETDEPS += $$(CPLEX_STUDIO_DIR128)/concert/lib/x64_windows_vs2017/stat_mdd/concert.lib
#}

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
#win32: target.path = $$(CPLEX_DIR)/lib/x64_windows_vs2017/stat_mda/
!isEmpty(target.path): INSTALLS += target

FORMS += \
    GeneratorWindow.ui



