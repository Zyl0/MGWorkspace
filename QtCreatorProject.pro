QT += core gui opengl
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
greaterThan(QT_MAJOR_VERSION, 5): QT += openglwidgets

CONFIG += c++11

INCLUDEPATH += AppTinyMesh/Include
INCLUDEPATH += $$(GLEW_DIR)
INCLUDEPATH += $$(OUT_PWD)

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp

VPATH += AppTinyMesh

SOURCES += \
    AppTinyMesh/Source/bezier.cpp \
    AppTinyMesh/Source/box.cpp \
    AppTinyMesh/Source/deformations.cpp \
    AppTinyMesh/Source/distancefieldhierarchy.cpp \
    AppTinyMesh/Source/evector.cpp \
    AppTinyMesh/Source/implicits.cpp \
    AppTinyMesh/Source/main.cpp \
    AppTinyMesh/Source/camera.cpp \
    AppTinyMesh/Source/mat.cpp \
    AppTinyMesh/Source/mesh.cpp \
    AppTinyMesh/Source/meshcolor.cpp \
    AppTinyMesh/Source/mesh-widget.cpp \
    AppTinyMesh/Source/primitives.cpp \
    AppTinyMesh/Source/qtemainwindow.cpp \
    AppTinyMesh/Source/ray.cpp \
    AppTinyMesh/Source/shader-api.cpp \
    AppTinyMesh/Source/triangle.cpp \
    AppTinyMesh/Source/vec.cpp

HEADERS += \
    AppTinyMesh/Include/bezier.h \
    AppTinyMesh/Include/box.h \
    AppTinyMesh/Include/camera.h \
    AppTinyMesh/Include/color.h \
    AppTinyMesh/Include/deformations.h \
    AppTinyMesh/Include/distancefieldhierarchy.h \
    AppTinyMesh/Include/implicits.h \
    AppTinyMesh/Include/mat.h \
    AppTinyMesh/Include/mathematics.h \
    AppTinyMesh/Include/mesh.h \
    AppTinyMesh/Include/meshcolor.h \
    AppTinyMesh/Include/primitives.h \
    AppTinyMesh/Include/qte.h \
    AppTinyMesh/Include/ray.h \
    AppTinyMesh/Include/realtime.h \
    AppTinyMesh/Include/shader-api.h \
    AppTinyMesh/Include/vec.h

FORMS += \
    AppTinyMesh/UI/interface.ui

win32 {
    LIBS += -L$$(GLEW_DIR) -lglew32
    LIBS += -lopengl32 -lglu32
}
unix:!macx {
    LIBS += -lGLEW -lGL -lGLU
}
macx {
    LIBS += -lGLEW -lGL -lGLU
}


# Copy shader files
# $$shell_quote puts quote around the path, to make it work if it contains space or other special characters.
copydata.commands = $(COPY_DIR) $$shell_quote($$PWD/AppTinyMesh/Shaders) $$shell_quote($$OUT_PWD/Shaders)
first.depends = $(first) copydata
export(first.depends)
export(copydata.commands)
QMAKE_EXTRA_TARGETS += first copydata

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
