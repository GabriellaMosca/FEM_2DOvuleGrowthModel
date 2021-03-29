#include "FemMembranes.hpp"
namespace mdx
{

  bool FemMembranes::initialize(QWidget *parent)
  {
    mesh = currentMesh();
    if(!mesh)
      throw(QString("FemMembranes::initialize No current mesh"));

    ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("FemMembranes::initialize No cell complex"));

    cs = &mesh->ccStructure(ccName);

    // Get the solver process
    if(!getProcess(parm("Solver Process"), solverProcess))
      throw(QString("FemMembranes::initialize Unable to make solver process: %1").arg(parm("Solver Process")));
    solverProcess->initialize(parent);

    // Get the growth process
    if(!getProcess(parm("Growth Process"), growthProcess))
      throw(QString("FemMembranes::initialize Unable to make growth process: %1").arg(parm("Growth Process")));
    growthProcess->initialize(parent);

    // Get the subdivide process
    if(!getProcess(parm("Subdivide Process"), subdivideProcess))
      throw(QString("FemMembranes::initialize Unable to make subdivide process: %1").arg(parm("Subdivide Process")));
    subdivideProcess->initialize(parent);

    savingTimes = parm("Simulation time to save mesh");
    listSavingTimes = savingTimes.split(',');
    
    if(!getProcess(parm("Save Mesh Process Name"), meshSave))
        throw(QString("FemMembranes::initialize Unable to make save mesh process : %1").arg(parm("Save Mesh Process Name")));
    meshSave->initialize(parent);

    return true;
  }

  bool FemMembranes::step()
  {
    if(!solverProcess)
      throw(QString("FemMembranes::run Solver process pointer invalid"));
    if(!growthProcess)
      throw(QString("FemMembranes::run Growth process pointer invalid"));
    if(!subdivideProcess)
      throw(QString("FemMembranes::run Subdivide process pointer invalid"));

    // If converged subdivide then grow
    if(!solverProcess->step()) {

      //take a screenshot
      static int screenShotCount = 0;
      QString fileName = QString("OvuleFEMGrowth-%1.JPG").arg(screenShotCount++, 4, 10, QChar('0'));
      takeSnapshot(fileName);
      for (int i=0; i<listSavingTimes.size(); ++i)
      {
        if(listSavingTimes[i].toInt() == screenShotCount)
        {
           meshName = parm("Name of mesh");
           QString fileNameMesh = QString("%1-%2.mdxm").arg(meshName).arg(screenShotCount, 4, 10, QChar('0'));
           meshSave->run(mesh, fileNameMesh, true);
           break;
        }

      }
      growthProcess->run();
      double growthDt = growthProcess->parm("Growth Dt").toDouble();
      growthTime += growthDt;
      if(stringToBool(solverProcess->parm("Print Stats")))
        mdxInfo << QString("Growth Step Time %1, step %2").arg(growthTime).arg(growthDt) << endl;

      // If we subdivide, we need to re-initialize the solver
      if(subdivideProcess->run()) {
        mesh->updateAll(ccName);
        // Reinitialize solver with new graph
        solverProcess->initSolver(cs);
      }
    }
    mesh->updatePositions(ccName);

    return true;
  }

  bool FemMembranes::rewind(QWidget *parent)
  {
    // To rewind, we'll reload the mesh
    Mesh *mesh = currentMesh();
    if(!mesh or mesh->file().isEmpty())
      throw(QString("No current mesh, cannot rewind"));
    MeshLoad meshLoad(*this);
    meshLoad.setParm("File Name", mesh->file());
    growthTime = 0;
    return meshLoad.run();
  }

  bool FemMembranes::finalize(QWidget *parent)
  {
    if(!solverProcess)
      throw(QString("FemMembranes::run Solver process pointer invalid"));

    bool result = solverProcess->finalize(parent);

    // Cleanup
    mesh = 0;
    solverProcess = 0;
    growthProcess = 0;
    subdivideProcess = 0;

    return result;
  }

  
  bool FemMembraneBisect::run()
  {
    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("FemMembraneSubdivide::run No current mesh"));
    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("FemMembraneSubdivide::run Invalid cell complex"));

    fem::ElasticTriangle3Attr &elementAttr = mesh->attributes().attrMap<CCIndex, fem::ElasticTriangle3>(parm("Element Attribute"));
    fem::TransIsoMaterialAttr &materialAttr = mesh->attributes().attrMap<CCIndex, fem::TransIsoMaterial>(parm("Material Attribute"));
    fem::PressureAttr &pressureAttr = mesh->attributes().attrMap<CCIndex, fem::Pressure>(parm("Pressure Attribute"));
    fem::PressureEdgeAttr &pressureEdgeAttr = mesh->attributes().attrMap<CCIndex, fem::PressureEdge>(parm("Pressure Edge Attribute"));
 
    fem::GrowthAttr &growthAttr = mesh->attributes().attrMap<CCIndex, fem::Growth>(parm("Growth Attribute"));
    fem::DirichletAttr &dirichletAttr = mesh->attributes().attrMap<CCIndex, fem::Dirichlet>(parm("Dirichlet Attribute"));

    fem::Triangle3Attr<Point1d> &elementMorphoAnisoAttr = mesh->attributes().attrMap<CCIndex, fem::Triangle3<Point1d>>(parm("Morphogen Anisotropy Element Attribute"));
    fem::Triangle3Attr<Point1d> &elementMorphoGrowthAttr = mesh->attributes().attrMap<CCIndex, fem::Triangle3<Point1d>>(parm("Morphogen Growth Element Attribute"));
    
    fem::MorphogenDataAttr &morphogenE2DataAttr = mesh->attributes().attrMap<CCIndex, fem::MorphogenData>(parm("Morphogen Anisotropy E2 Data"));
    fem::MorphogenDataAttr &morphogenKParDataAttr = mesh->attributes().attrMap<CCIndex, fem::MorphogenData>(parm("Morphogen Anisotropy KPar Data"));
    fem::MorphogenDataAttr &morphogenGrowth1DataAttr = mesh->attributes().attrMap<CCIndex, fem::MorphogenData>(parm("Morphogen Growth Data 1"));
    fem::MorphogenDataAttr &morphogenGrowth2DataAttr = mesh->attributes().attrMap<CCIndex, fem::MorphogenData>(parm("Morphogen Growth Data 2"));
     
    fem::MorphogenDirichletAttr &moprhogenDirichletE2Attr = mesh->attributes().attrMap<CCIndex, fem::MorphogenDirichlet>(parm("Morphogen Anisotropy E2 Dirichlet"));
    fem::MorphogenDirichletAttr &moprhogenDirichletKParAttr = mesh->attributes().attrMap<CCIndex, fem::MorphogenDirichlet>(parm("Morphogen Anisotropy KPar Dirichlet"));
    fem::MorphogenDirichletAttr &moprhogenDirichletGrowth1Attr = mesh->attributes().attrMap<CCIndex, fem::MorphogenDirichlet>(parm("Morphogen Growth Dirichlet 1"));
    fem::MorphogenDirichletAttr &moprhogenDirichletGrowth2Attr = mesh->attributes().attrMap<CCIndex, fem::MorphogenDirichlet>(parm("Morphogen Growth Dirichlet 2"));


    FemMembraneSubdivide sDiv(*mesh, elementAttr, elementMorphoAnisoAttr, elementMorphoGrowthAttr, materialAttr, pressureAttr, pressureEdgeAttr, dirichletAttr, growthAttr, morphogenE2DataAttr, morphogenKParDataAttr, 
                           morphogenGrowth1DataAttr, morphogenGrowth2DataAttr, moprhogenDirichletE2Attr, moprhogenDirichletKParAttr, moprhogenDirichletGrowth1Attr, moprhogenDirichletGrowth2Attr);
    mesh->updateAll(ccName);

    return run(*mesh, mesh->ccStructure(ccName), parm("Max Area").toDouble(), &sDiv);
  }



  ////////////////////////////////////////////////////////

  REGISTER_PROCESS(FemMembranes);
  REGISTER_PROCESS(FemMembraneSolver);
  REGISTER_PROCESS(FemMembraneRefCfg);
  REGISTER_PROCESS(FemMembraneStressStrain);
  REGISTER_PROCESS(FemMembraneDerivs);
  REGISTER_PROCESS(FemMembraneMaterial);
  REGISTER_PROCESS(FemMembraneCellFacesMaterial);
  REGISTER_PROCESS(FemMembraneMaterialMorphogenFaces);
  REGISTER_PROCESS(FemMembraneAnisoDir);
  REGISTER_PROCESS(FemMembraneAnisoDirMorphogens);
  REGISTER_PROCESS(FemMembranePressure);
  REGISTER_PROCESS(FemMembranePressureDerivs);
  REGISTER_PROCESS(FemMembraneEdgePressure);
  REGISTER_PROCESS(FemMembranePressureEdgeDerivs);
 
  REGISTER_PROCESS(FemMembraneSetGrowth);
  REGISTER_PROCESS(FemMembraneCellFacesGrowth);
  REGISTER_PROCESS(FemMembraneSetDirichlet);
  REGISTER_PROCESS(FemMembraneDirichletDerivs);
  REGISTER_PROCESS(FemMembraneGrowth);
  REGISTER_PROCESS(FemMembraneSetGrowthMorphogens);

  REGISTER_PROCESS(FemMembraneBisect);
  REGISTER_PROCESS(FemMembraneDiffusionAnisotropyDerivs);
  REGISTER_PROCESS(FemMembraneDiffusionAnisotropyDirichletDerivs);
  REGISTER_PROCESS(FemMembraneDiffusionAnisotropySetDirichlet);
  REGISTER_PROCESS(FemMembraneDiffusionAnisotropyVisualize);
  REGISTER_PROCESS(FemMembraneDiffusionAnisotropySolver);
  REGISTER_PROCESS(FemMembraneCreateAnisotropyMorphogenElement);
  REGISTER_PROCESS(FemMembraneDiffusionGrowthDerivs);
  REGISTER_PROCESS(FemMembraneDiffusionGrowthDirichletDerivs);
  REGISTER_PROCESS(FemMembraneDiffusionGrowthSetDirichlet);
  REGISTER_PROCESS(FemMembraneDiffusionGrowthVisualize);
  REGISTER_PROCESS(FemMembraneDiffusionGrowthSolver);
  REGISTER_PROCESS(FemMembraneCreateGrowthMorphogenElement);


  REGISTER_PROCESS(FemMembraneVisMaterial);
  REGISTER_PROCESS(FemMembraneVisGrowth);
  REGISTER_PROCESS(FemMembraneComputeCurrentDirections);
  REGISTER_PROCESS(FemMembraneVisDirections);
  //REGISTER_PROCESS(FemMembraneRender);
  REGISTER_PROCESS(FemAnisotropyPropagationFailure);

}
