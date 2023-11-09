#ifndef __Qte__
#define __Qte__

#include <QtWidgets/qmainwindow.h>
#include "realtime.h"
#include "meshcolor.h"
#include "implicits.h"

QT_BEGIN_NAMESPACE
	namespace Ui { class Assets; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
  Q_OBJECT
private:
  Ui::Assets* uiw;           //!< Interface

  MeshWidget* meshWidget;   //!< Viewer
  MeshColor meshColor;		//!< Mesh.

public:
  MainWindow();
  ~MainWindow();
  void CreateActions();
  void UpdateGeometry();

public slots:
  void editingSceneLeft(const Ray&);
  void editingSceneRight(const Ray&);
  void editingErosion(const Ray&);
  void BoxMeshExample();
  void SphereImplicitExample();
  void SpherePrimitive();
  void BoxPrimitive();
  void CapsulePrimitive();
  void TorePrimitive();
  void Shape1();
  void rollbackTree();
  void TransformTree();
  void Refresh();
  void SaveMesh();
  void MeshifyImplicit(AnalyticScalarField* implicit);
  void ResetCamera();
  void UpdateMaterial();
};

#endif
