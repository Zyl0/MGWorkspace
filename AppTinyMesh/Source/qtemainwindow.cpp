#include "qte.h"
#include "implicits.h"
#include "ui_interface.h"

#include "distancefieldhierarchy.h"
#include "primitives.h"

#include "bezier.h"
#include "deformations.h"

#include "mat.h"

#include <vector>
#include <list>
#include <random>
#include <chrono>
#include <iostream>

AnalyticScalarField *implicitTree = nullptr;

inline double getBrushSphereSize(Ui::Assets *uiw)
{
  return uiw->PaintBrushSphereSize->text().toDouble();
}

inline double getBrushAngle(Ui::Assets *uiw)
{
  return uiw->PaintBrushAngle->text().toDouble();
}

inline double getBrushCount(Ui::Assets *uiw)
{
  return uiw->PaintBrushCount->text().toDouble();
}

inline double getBlendRadius(Ui::Assets *uiw)
{
  return uiw->SmoothRadius->text().toDouble();
}

inline Transform getTransform(Ui::Assets *uiw)
{
  return Scale(
             uiw->TScaleX->text().toDouble(),
             uiw->TScaleY->text().toDouble(),
             uiw->TScaleZ->text().toDouble()) *
         RotationX(uiw->TRotateX->text().toDouble()) * RotationY(uiw->TRotateY->text().toDouble()) * RotationZ(uiw->TRotateZ->text().toDouble()) * Translation(uiw->TTranslateX->text().toDouble(), uiw->TTranslateY->text().toDouble(), uiw->TTranslateZ->text().toDouble());
}

inline AnalyticScalarField *getSphere(Ui::Assets *uiw)
{
  return new SDFSphere(
      uiw->SphereRadiusInput->text().toDouble(),
      Vector(
          uiw->SphereOriginXInput->text().toDouble(),
          uiw->SphereOriginYInput->text().toDouble(),
          uiw->SphereOriginZInput->text().toDouble()));
}

inline AnalyticScalarField *getBox(Ui::Assets *uiw)
{
  return new SDFBox(
      0.0,
      Vector(
          uiw->BoxAXInput->text().toDouble(),
          uiw->BoxAYInput->text().toDouble(),
          uiw->BoxAZInput->text().toDouble()),
      Vector(
          uiw->BoxBXInput->text().toDouble(),
          uiw->BoxBYInput->text().toDouble(),
          uiw->BoxBZInput->text().toDouble()));
}

inline AnalyticScalarField *getCapsule(Ui::Assets *uiw)
{
  return new SDFCapsule(
      uiw->CapsuleRadiusInput->text().toDouble(),
      Vector(
          uiw->CapsuleAXInput->text().toDouble(),
          uiw->CapsuleAYInput->text().toDouble(),
          uiw->CapsuleAZInput->text().toDouble()),
      Vector(
          uiw->CapsuleBXInput->text().toDouble(),
          uiw->CapsuleBYInput->text().toDouble(),
          uiw->CapsuleBZInput->text().toDouble()));
}

inline AnalyticScalarField *getTore(Ui::Assets *uiw)
{
  return new SDFTore(
      uiw->ToreInnerRadiusInput->text().toDouble(),
      uiw->ToreRadiusInput->text().toDouble(),
      Vector(
          uiw->ToreOriginXInput->text().toDouble(),
          uiw->ToreOriginYInput->text().toDouble(),
          uiw->ToreOriginZInput->text().toDouble()));
}

void deleteTree()
{
  if (implicitTree != nullptr)
  {
    delete implicitTree;
    implicitTree = nullptr;
  }
}

void appendTree(Ui::Assets *uiw, AnalyticScalarField *node)
{
  if (implicitTree == nullptr)
  {
    implicitTree = node;
    return;
  }

  if (uiw->Union->isChecked())
  {
    implicitTree = new HDFUnion(
        implicitTree,
        node);

    return;
  }

  if (uiw->Blend->isChecked())
  {
    implicitTree = new HDFBlend(
        implicitTree,
        node,
        getBlendRadius(uiw));

    return;
  }

  if (uiw->Diff->isChecked())
  {
    implicitTree = new HDFDiff(
        implicitTree,
        node);

    return;
  }

  if (uiw->Intersect->isChecked())
  {
    implicitTree = new HDFIntersection(
        implicitTree,
        node);

    return;
  }

  if (uiw->SmoothUnion->isChecked())
  {
    implicitTree = new HDFSmouthUnion(
        implicitTree,
        node,
        getBlendRadius(uiw));

    return;
  }
}

MainWindow::MainWindow() : QMainWindow(), uiw(new Ui::Assets)
{
  // Chargement de l'interface
  uiw->setupUi(this);

  // Chargement du GLWidget
  meshWidget = new MeshWidget;
  QGridLayout *GLlayout = new QGridLayout;
  GLlayout->addWidget(meshWidget, 0, 0);
  GLlayout->setContentsMargins(0, 0, 0, 0);
  uiw->widget_GL->setLayout(GLlayout);

  // Creation des connect
  CreateActions();

  meshWidget->SetCamera(Camera(Vector(10, 0, 0), Vector(0.0, 0.0, 0.0)));
}

MainWindow::~MainWindow()
{
  delete meshWidget;
}

void MainWindow::CreateActions()
{
  // Buttons
  connect(uiw->boxMesh, SIGNAL(clicked()), this, SLOT(BoxMeshExample()));
  connect(uiw->sphereImplicit, SIGNAL(clicked()), this, SLOT(SphereImplicitExample()));
  connect(uiw->resetcameraButton, SIGNAL(clicked()), this, SLOT(ResetCamera()));
  connect(uiw->wireframe, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
  connect(uiw->radioShadingButton_1, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
  connect(uiw->radioShadingButton_2, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
  connect(uiw->BezierMesh, SIGNAL(clicked()), this, SLOT(BezierSurfacePreview()));
  connect(uiw->RevolutionMesh, SIGNAL(clicked()), this, SLOT(RevolutionSurfacePreview()));
  connect(uiw->TwistSDF, SIGNAL(clicked()), this, SLOT(TwistSDF()));
  connect(uiw->SpawnPrimitiveSphere, SIGNAL(clicked()), this, SLOT(SpherePrimitive()));
  connect(uiw->BoxPrimitiveSphere, SIGNAL(clicked()), this, SLOT(BoxPrimitive()));
  connect(uiw->CapsulePrimitiveSphere, SIGNAL(clicked()), this, SLOT(CapsulePrimitive()));
  connect(uiw->TorePrimitiveSphere, SIGNAL(clicked()), this, SLOT(TorePrimitive()));
  connect(uiw->SpawnShape1, SIGNAL(clicked()), this, SLOT(Shape1()));
  connect(uiw->rb, SIGNAL(clicked()), this, SLOT(rollbackTree()));
  connect(uiw->SDFTransform, SIGNAL(clicked()), this, SLOT(TransformTree()));
  connect(uiw->refreshMesh, SIGNAL(clicked()), this, SLOT(Refresh()));
  connect(uiw->saveMesh, SIGNAL(clicked()), this, SLOT(SaveMesh()));

  // Widget edition
  connect(meshWidget, SIGNAL(_signalEditSceneLeft(const Ray &)), this, SLOT(editingSceneLeft(const Ray &)));
  connect(meshWidget, SIGNAL(_signalEditSceneRight(const Ray &)), this, SLOT(editingSceneRight(const Ray &)));
  connect(meshWidget, SIGNAL(_signalErosion(const Ray &)), this, SLOT(editingErosion(const Ray &)));

  // Marching cube
  uiw->MCResInput->setText(QString::number(64));
  uiw->MCAXInput->setText(QString::number(-2.f));
  uiw->MCAYInput->setText(QString::number(-2.f));
  uiw->MCAZInput->setText(QString::number(-2.f));
  uiw->MCBXInput->setText(QString::number(2.f));
  uiw->MCBYInput->setText(QString::number(2.f));
  uiw->MCBZInput->setText(QString::number(4.f));

  // Sphere
  uiw->SphereRadiusInput->setText(QString::number(1.f));
  uiw->SphereOriginXInput->setText(QString::number(0.f));
  uiw->SphereOriginYInput->setText(QString::number(0.f));
  uiw->SphereOriginZInput->setText(QString::number(0.f));

  // Capsule
  uiw->CapsuleRadiusInput->setText(QString::number(1.f));
  uiw->CapsuleAXInput->setText(QString::number(0.f));
  uiw->CapsuleAYInput->setText(QString::number(0.f));
  uiw->CapsuleAZInput->setText(QString::number(0.f));
  uiw->CapsuleBXInput->setText(QString::number(0.f));
  uiw->CapsuleBYInput->setText(QString::number(0.f));
  uiw->CapsuleBZInput->setText(QString::number(1.f));

  // Box
  uiw->BoxAXInput->setText(QString::number(-1.f));
  uiw->BoxAYInput->setText(QString::number(-1.f));
  uiw->BoxAZInput->setText(QString::number(0.f));
  uiw->BoxBXInput->setText(QString::number(1.f));
  uiw->BoxBYInput->setText(QString::number(1.f));
  uiw->BoxBZInput->setText(QString::number(2.f));

  // Tore
  uiw->ToreRadiusInput->setText(QString::number(1.f));
  uiw->ToreInnerRadiusInput->setText(QString::number(0.25f));
  uiw->ToreOriginXInput->setText(QString::number(0.f));
  uiw->ToreOriginYInput->setText(QString::number(0.f));
  uiw->ToreOriginZInput->setText(QString::number(1.f));

  // Transform
  uiw->TTranslateX->setText(QString::number(0.0f));
  uiw->TTranslateY->setText(QString::number(0.0f));
  uiw->TTranslateZ->setText(QString::number(0.0f));
  uiw->TRotateX->setText(QString::number(0.0f));
  uiw->TRotateY->setText(QString::number(0.0f));
  uiw->TRotateZ->setText(QString::number(0.0f));
  uiw->TScaleX->setText(QString::number(1.0f));
  uiw->TScaleY->setText(QString::number(1.0f));
  uiw->TScaleZ->setText(QString::number(1.0f));

  // Append
  uiw->Union->setChecked(true);
}

void MainWindow::editingSceneLeft(const Ray &ray)
{
}

void MainWindow::editingSceneRight(const Ray &ray)
{
}

void MainWindow::editingErosion(const Ray &ray)
{
  if (implicitTree == nullptr)
    return;

  // Ray marching settings
  const double maxTrace = 100.;
  const double deltaTrace = 0.05;

  // Paint UI param
  const unsigned int count = getBrushCount(uiw);
  const double angle = getBrushAngle(uiw);
  const double sphereSize = getBrushSphereSize(uiw);

  const bool useRandom = (count != 0);
  std::vector<Vector> spheres;

  // random generator
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<double> d{0.0, 1.0};

  Box bounds = Box(Vector(uiw->MCAXInput->text().toDouble(), uiw->MCAYInput->text().toDouble(), uiw->MCAZInput->text().toDouble()), Vector(uiw->MCBXInput->text().toDouble(), uiw->MCBYInput->text().toDouble(), uiw->MCBZInput->text().toDouble()));

  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<float>> sp;
  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<float>> ep;

  sp = std::chrono::steady_clock::now();
  for (unsigned int i = 0; i < count; i++)
  {
    Vector p = ray.Origin();
    Transform rotation = RotationX(d(gen) * angle) * RotationY(d(gen) * angle) * RotationZ(d(gen) * angle);
    Vector direction = rotation(ray.Direction());

    // Ray marching
    double traceDist = 0.;
    while ((traceDist += deltaTrace) < maxTrace)
    {
      p = p + (deltaTrace * direction);
      if (!bounds.Inside(p)) continue;
      
      double v = implicitTree->Value(p);
      if (v < 0.)
      {
        spheres.push_back(p);
        break;
      }
    }
  }
  ep = std::chrono::steady_clock::now();

  std::cout << "Ray casting (by ray marching) of " << count << " rays took " << (ep - sp).count()*1000. << " ms" << std::endl;

  if(spheres.size() == 0)
    return;
  
  if(spheres.size() == 1)
    implicitTree = new HDFDiff
      (
        implicitTree,
        new SDFSphere(sphereSize, spheres[0])
      );
  else
  {
    HierarchalDistanceField *modifier = nullptr;

    std::list<HierarchalDistanceField*> nodes;

    for (size_t i = 0; i < spheres.size(); i+=2)
    {
      nodes.push_back(new HDFUnion
      (
        new SDFSphere(sphereSize, spheres[i]),
        new SDFSphere(sphereSize, spheres[i+1])
      ));
    }
    
    while (nodes.size() > 1)
    {
      size_t sphereCount = nodes.size();

      for (size_t i = 0; i < sphereCount; i+=2)
      {
        HierarchalDistanceField *A = nodes.front();
        nodes.pop_front();
        HierarchalDistanceField *B = nodes.front();
        nodes.pop_front();
        nodes.push_back(new HDFUnion(A, B));
      }
    }

    appendTree(uiw, nodes.front());
  }

  MeshifyImplicit(implicitTree);
}

void MainWindow::BoxMeshExample()
{
  Mesh boxMesh = Mesh(Box(1.0));

  std::vector<Color> cols;
  cols.resize(boxMesh.Vertexes());
  for (size_t i = 0; i < cols.size(); i++)
    cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

  meshColor = MeshColor(boxMesh, cols, boxMesh.VertexIndexes());
  UpdateGeometry();
  deleteTree();
}

void MainWindow::BezierSurfacePreview()
{
  //  size_t n = 3, m = 3;
  //  std::vector<Vector> controlPoints =
  //  {
  //    Vector(-1,-1, 1), Vector( 0,-1,-1), Vector( 1,-1, 1),
  //    Vector(-1, 0,-1), Vector( 0, 0, 5), Vector( 1, 0,-1),
  //    Vector(-1, 1, 1), Vector( 0, 1,-1), Vector( 1, 1, 1)
  //  };
  size_t n = 4, m = 3;
  std::vector<Vector> controlPoints = 
  {
    Vector(-1,-1,-1), Vector( 0,-1,-1), Vector( 1,-1, 1), Vector(-1,-1, 1),
    Vector(-1, 0,-1), Vector( 0, 0, 0), Vector( 1, 0, 0), Vector(-1, 0, 0.5),
    Vector(-1, 1, 1), Vector( 0, 1,-2), Vector( 1, 1, 0), Vector( 1, 1, 0)
  };
  BezierSurface b(n, m, controlPoints);

  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<float>> sp;
  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<float>> ep;
  
  sp = std::chrono::steady_clock::now();
  Mesh BezierMesh = b.Polygonize(100,100);
  ep = std::chrono::steady_clock::now();
  std::cout << "Meshify (by sampling bezier surface) took " << (ep - sp).count()*1000. << " ms" << std::endl;

  BezierMesh.SaveObj("./BezierSurface2.obj", "bPlan");

  std::vector<Color> cols;
  cols.resize(BezierMesh.Vertexes());
  for (size_t i = 0; i < cols.size(); i++)
    cols[i] = Color(0.8, 0.8, 0.8);

  meshColor = MeshColor(BezierMesh, cols, BezierMesh.VertexIndexes());
  UpdateGeometry();
}

void MainWindow::TwistSDF()
{
  //if (!implicitTree)
  //  return;
  //std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<float>> sp;
  //std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<float>> ep;

  //Mesh implicitMesh;
  //sp = std::chrono::steady_clock::now();
  //implicitTree->Polygonize(
  //    uiw->MCResInput->text().toDouble(),
  //    implicitMesh,
  //    Box(Vector(uiw->MCAXInput->text().toDouble(), uiw->MCAYInput->text().toDouble(), uiw->MCAZInput->text().toDouble()), Vector(uiw->MCBXInput->text().toDouble(), uiw->MCBYInput->text().toDouble(), uiw->MCBZInput->text().toDouble()))
  //    );
  //ep = std::chrono::steady_clock::now();
  //std::cout << "Meshify (by marching cube) of " << uiw->MCResInput->text().toDouble() << " box count took " << (ep - sp).count()*1000. << " ms" << std::endl;


    size_t n = 3, m = 3;
    std::vector<Vector> controlPoints =
    {
      Vector(-1,-1, 1), Vector( 0,-1,-1), Vector( 1,-1, 1),
      Vector(-1, 0,-1), Vector( 0, 0, 5), Vector( 1, 0,-1),
      Vector(-1, 1, 1), Vector( 0, 1,-1), Vector( 1, 1, 1)
    };
    BezierSurface b(n, m, controlPoints);

    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<float>> sp;
    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<float>> ep;

    sp = std::chrono::steady_clock::now();
    Mesh BezierMesh = b.Polygonize(100,100);
    ep = std::chrono::steady_clock::now();

  Twist twist = Twist(8, Vector(0,0,1));
  Mesh mesh = twist.WarpMesh(BezierMesh);

  mesh.SaveObj("./BezierTwistF8.obj", "rSurface");

  std::vector<Color> cols;
  cols.resize(mesh.Vertexes());
  for (size_t i = 0; i < cols.size(); i++)
    cols[i] = Color(0.8, 0.8, 0.8);

  meshColor = MeshColor(mesh, cols, mesh.VertexIndexes());
  UpdateGeometry();
}

void MainWindow::RevolutionSurfacePreview()
{
  //std::vector<Vector> controlPoints =
  //{
  //  Vector(0, 1.5, -1), Vector( 0, 1.25, 0), Vector( 0, 0.125, 1), Vector( 0, 0.125, 1.25), Vector( 0, 0.5, 1.5)
  //};

  std::vector<Vector> controlPoints =
  {
    Vector(0, 0.01, -1), Vector(0, 1, -1), Vector( 0, 1.25, -1), Vector( 0, 1.5, -0.5), Vector( 0, 0.125, 1.25), Vector( 0, 0.125, 1.25), Vector( 0, 0.125, 1.25), Vector( 0, 0.25, 1.5), Vector( 0, 1.5, 1.7), Vector( 0, 0, 2)
  };

  Revolution revolutionSurvace(Vector(0,0,-1), Vector(0,0,1.5),  controlPoints);

  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<float>> sp;
  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<float>> ep;
  
  sp = std::chrono::steady_clock::now();
  Mesh revolutionMesh = revolutionSurvace.Polygonize(100,100);
  ep = std::chrono::steady_clock::now();
  std::cout << "Meshify (by sampling revolution surface) took " << (ep - sp).count()*1000. << " ms" << std::endl;

  revolutionMesh.SaveObj("./RevolutionSurface2.obj", "rSurface");

  std::vector<Color> cols;
  cols.resize(revolutionMesh.Vertexes());
  for (size_t i = 0; i < cols.size(); i++)
    cols[i] = Color(0.8, 0.8, 0.8);

  meshColor = MeshColor(revolutionMesh, cols, revolutionMesh.VertexIndexes());
  UpdateGeometry();
}

void MainWindow::SphereImplicitExample()
{
  AnalyticScalarField *implicit =
      new HDFTransform(
          new SDFTore(0.25, 1, Vector(0, 0, 0)),
          RotationX(90));

  MeshifyImplicit(implicit);

  delete implicit;

  deleteTree();
}

void MainWindow::SpherePrimitive()
{
  AnalyticScalarField *implicit = getSphere(uiw);

  appendTree(uiw, implicit);

  MeshifyImplicit(implicitTree);
}

void MainWindow::BoxPrimitive()
{
  AnalyticScalarField *implicit = getBox(uiw);

  appendTree(uiw, implicit);

  MeshifyImplicit(implicitTree);
}

void MainWindow::CapsulePrimitive()
{
  AnalyticScalarField *implicit = getCapsule(uiw);

  appendTree(uiw, implicit);

  MeshifyImplicit(implicitTree);
}

void MainWindow::TorePrimitive()
{
  AnalyticScalarField *implicit = getTore(uiw);

  appendTree(uiw, implicit);

  MeshifyImplicit(implicitTree);
}

void MainWindow::Shape1()
{
  AnalyticScalarField *implicit =
      new HDFTransform(
          new HDFUnion(
              new HDFDiff(
                  new HDFUnion(
                      new HDFUnion(
                          new HDFDiff(
                              new HDFDiff(
                                  new SDFBox(0, Vector(-1, -1, 0), Vector(1, 1, 1)),
                                  new SDFBox(0, Vector(-1, -1, 0.35), Vector(1, 1, 0.75))),
                              new SDFBox(0, Vector(-0.9, -0.9, 0), Vector(0.9, 0.9, 0.1))),
                          new SDFBox(0, Vector(-.8, -.8, 0.25), Vector(.8, .8, 0.85))),
                      new HDFUnion(
                          new SDFSphere(1, Vector(0, 0, 1)),
                          new HDFUnion(
                              new HDFTransform(
                                  new HDFUnion(
                                      new SDFTore(0.05, 1, Vector(0, 0.15, 0)),
                                      new SDFTore(0.05, 1, Vector(0, -0.15, 0))),
                                  Translation(0, 0, 1) * RotationX(45)),
                              new HDFTransform(
                                  new SDFTore(0.1, 1, Vector(0, 0, 0)),
                                  Translation(0, 0, 1) * RotationX(-45))))),
                  new SDFSphere(1.25, Vector(0, 1, 2))),
              new SDFSphere(0.75, Vector(0, 0, 1))),
          Translation(0, 0, -1));

  appendTree(uiw, implicit);

  MeshifyImplicit(implicitTree);
}

void MainWindow::rollbackTree()
{
  HierarchalDistanceField *tree = dynamic_cast<HierarchalDistanceField *>(implicitTree);
  if (!tree)
    return;

  AnalyticScalarField *son = tree->PopLeftSon();
  deleteTree();

  appendTree(uiw, son);

  MeshifyImplicit(implicitTree);
}

void MainWindow::TransformTree()
{
  if (implicitTree == nullptr)
    return;

  Transform t = getTransform(uiw);
  implicitTree = new HDFTransform(implicitTree, t);
  MeshifyImplicit(implicitTree);
}

void MainWindow::Refresh()
{
  if (implicitTree == nullptr)
    return;
  MeshifyImplicit(implicitTree);
}

void MainWindow::SaveMesh()
{
  if (implicitTree == nullptr)
    return;

  Mesh implicitMesh;
  implicitTree->Polygonize(
      uiw->MCResInput->text().toDouble(),
      implicitMesh, Box(Vector(uiw->MCAXInput->text().toDouble(), uiw->MCAYInput->text().toDouble(), uiw->MCAZInput->text().toDouble()), Vector(uiw->MCBXInput->text().toDouble(), uiw->MCBYInput->text().toDouble(), uiw->MCBZInput->text().toDouble())));


  implicitMesh.SaveObj("./DistanceFieldObject.obj", "implicit");
}

void MainWindow::MeshifyImplicit(AnalyticScalarField *implicit)
{
  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<float>> sp;
  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<float>> ep;

  Mesh implicitMesh;
  sp = std::chrono::steady_clock::now();
  implicit->Polygonize(
      uiw->MCResInput->text().toDouble(),
      implicitMesh, 
      Box(Vector(uiw->MCAXInput->text().toDouble(), uiw->MCAYInput->text().toDouble(), uiw->MCAZInput->text().toDouble()), Vector(uiw->MCBXInput->text().toDouble(), uiw->MCBYInput->text().toDouble(), uiw->MCBZInput->text().toDouble()))
      );
  ep = std::chrono::steady_clock::now();
  std::cout << "Meshify (by marching cube) of " << uiw->MCResInput->text().toDouble() << " box count took " << (ep - sp).count()*1000. << " ms" << std::endl;

  //implicitMesh.SaveObj("./", "implicit_export");

  std::vector<Color> cols;
  cols.resize(implicitMesh.Vertexes());
  for (size_t i = 0; i < cols.size(); i++)
    cols[i] = Color(0.8, 0.8, 0.8);

  meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
  UpdateGeometry();
}

void MainWindow::UpdateGeometry()
{
  meshWidget->ClearAll();
  meshWidget->AddMesh("BoxMesh", meshColor);

  uiw->lineEdit->setText(QString::number(meshColor.Vertexes()));
  uiw->lineEdit_2->setText(QString::number(meshColor.Triangles()));

  UpdateMaterial();
}

void MainWindow::UpdateMaterial()
{
  meshWidget->UseWireframeGlobal(uiw->wireframe->isChecked());

  if (uiw->radioShadingButton_1->isChecked())
    meshWidget->SetMaterialGlobal(MeshMaterial::Normal);
  else
    meshWidget->SetMaterialGlobal(MeshMaterial::Color);
}

void MainWindow::ResetCamera()
{
  meshWidget->SetCamera(Camera(Vector(-10.0), Vector(0.0)));
}
