#ifndef FEM_MEMBRANES_HPP
#define FEM_MEMBRANES_HPP

#include <MDXProcessFem.hpp>
#include <MDXProcessFemMorphogen.hpp>
#include <MeshProcessSystem.hpp>
#include <MeshProcessStructure.hpp>
#include <MDXProcessTissue.hpp>
#include <MDXProcessCellDivide.hpp>
#include <Attributes.hpp>
#include <MeshProcessSystemRender.hpp>
#include <MeshUtils.hpp>
#include <MeshProcessSelection.hpp>
#include <MeshProcessMeasures3D.hpp>
#include <Triangle3.hpp>

#include <Solver.hpp>
namespace mdx
{
  class FemMembranes;
  class FemMembraneSolver;
  class FemMembraneGrowth;
  class FemMembraneBisect;

  // Main model class
  class FemMembranes : public Process
  {
  public:
    FemMembranes(const Process &proc) : Process(proc) 
    {
      setName("Model/CCF/01 FEM Membranes");
      setDesc("FEM simulation with growth and subdivision");

      addParm("Solver Process", "Name of solver process", "Model/CCF/02 FEM Solver");
      addParm("Growth Process", "Name of growth process", "Model/CCF/13 Growth");
      addParm("Subdivide Process", "Name of subdivision process", "Model/CCF/14 Subdivide");
      addParm("Simulation time to save mesh", "At which growth time points save the mesh -- need to be in crescent order", "49, 46, 32" );
      addParm("Save Mesh Process Name", "Name of the Process to save the mesh", "Mesh/System/Save" );
      addParm("Name of mesh", "Name of the mesh to be saved", "OvuleFEMGrowth");
    }
    bool initialize(QWidget *parent);
    bool step();
    bool rewind(QWidget *parent);
    bool finalize(QWidget *parent);
    
    double growthTime = 0;
    
  private:
    Mesh *mesh = 0;
    QString ccName;
    CCStructure *cs = 0;

    FemMembraneSolver *solverProcess = 0;
    FemMembraneGrowth *growthProcess = 0;
    FemMembraneBisect *subdivideProcess = 0;
    MeshSave *meshSave = 0;

    QString savingTimes;
    QStringList listSavingTimes;
    QString meshName;

    //double growthTime = 0;
  };

  class FemMembraneSolver : public fem::FemSolver
  {
  public:
    FemMembraneSolver(const Process &proc) : fem::FemSolver(proc) 
    {
      setName("Model/CCF/02 FEM Solver");
      setDesc("FEM Simulation using triangular membrane elements");

      // Update parameters with our own defaults
      setParmDefault("Stress-Strain", "Model/CCF/05 Stress-Strain");

      // Add derivatives processes
      addParm("Element Derivs", "Process for element derivatives", "Model/CCF/03 Triangle Derivs");
      addParm("Pressure Derivs", "Process for pressure derivatives", "Model/CCF/10 Pressure/10a Pressure Derivs"); 
      addParm("Pressure Edge Derivs", "Process for pressure derivatives", "Model/CCF/10 Pressure/10b Pressure Edge Derivs"); 

      addParm("Dirichlet Derivs", "Process for Dirichlet derivatives", "Model/CCF/16 Dirichlet Derivs");
    }
  };

  class FemMembraneDerivs : public fem::ElementDerivs
  {
  public:
    FemMembraneDerivs(const Process &proc) : ElementDerivs(proc) 
    {
      setName("Model/CCF/03 Triangle Derivs");

      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
    }
  };

  class FemMembraneRefCfg : public fem::SetRefCfg
  {
  public:
    FemMembraneRefCfg(const Process &proc) : SetRefCfg(proc) 
    {
      setName("Model/CCF/04 Reference Configuration");

			addParm("Thickness", "Thickness of the membrane element", "1.0");
			setParmDefault("Element Type", "Linear Triangle");
			setParmDefault("Element Attribute", "Triangle Element");
    }
  };

  class FemMembraneStressStrain : public fem::StressStrain
  {
  public:
    FemMembraneStressStrain(const Process &proc) : StressStrain(proc) 
    {
      setName("Model/CCF/05 Stress-Strain");

      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
    }
  };

  class FemMembraneMaterial : public fem::SetTransIsoMaterial
  {
  public:
    FemMembraneMaterial(const Process &proc) : SetTransIsoMaterial(proc) 
    {
      setName("Model/CCF/06 Material Properties/06a Material Properties Uniform");
    }
  };
  class FemMembraneMaterialMorphogenFaces: public fem::SetTransIsoMaterialMorphogensFaces
  {
  public:
    FemMembraneMaterialMorphogenFaces(const Process &proc) : SetTransIsoMaterialMorphogensFaces(proc) 
    {
      setName("Model/CCF/06 Material Properties/06b Material Properties based on Morphogens (for faces)");
    }
  };


  
  class FemMembraneCellFacesMaterial : public fem::CellsToFacesTransIsoMaterial
  {
  public:
    FemMembraneCellFacesMaterial(const Process &proc) : CellsToFacesTransIsoMaterial(proc) 
    {
      setName("Model/CCF/07 Cells to Faces Material");
    }
  };

  class FemMembraneAnisoDir : public fem::SetAnisoDir
  {
  public:
    FemMembraneAnisoDir(const Process &proc) : SetAnisoDir(proc) 
    {
      setName("Model/CCF/08 Set Anisotropy/08a Set Aniso Dir Uniform");

      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
    }
  };

  class FemMembraneAnisoDirMorphogens : public fem::SetAnisoDirMorphogens
  {
  public:
    FemMembraneAnisoDirMorphogens(const Process &proc) : SetAnisoDirMorphogens(proc) 
    {
      setName("Model/CCF/08 Set Anisotropy/08b Set Aniso Dir Morphogens");

      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
      addParm("Morphogen Signal", "Morphogen Signal name",  "E2Morphogen Signal", QStringList() << "E2Morphogen Signal" << "KPArMorphogen Signal" );
    }
  };

  class FemMembranePressure : public fem::SetPressure
  {
  public:
    FemMembranePressure(const Process &proc) : SetPressure(proc) 
    {
      setName("Model/CCF/09 Set Pressure/09a Set Pressure");
    }
  };

  class FemMembranePressureDerivs : public fem::PressureDerivs
  {
  public:
    FemMembranePressureDerivs(const Process &proc) : PressureDerivs(proc) 
    {
      setName("Model/CCF/10 Pressure/10a Pressure Derivs");
    }
  };

  class FemMembraneEdgePressure : public fem::SetEdgePressure
  {
  public:
    FemMembraneEdgePressure(const Process &proc) : SetEdgePressure(proc) 
    {
      setName("Model/CCF/09 Set Pressure/09b Set Edge Pressure");
    }
  };

  class FemMembranePressureEdgeDerivs : public fem::PressureEdgeDerivs
  {
  public:
    FemMembranePressureEdgeDerivs(const Process &proc) : PressureEdgeDerivs(proc) 
    {
      setName("Model/CCF/10 Pressure/10b Pressure Edge Derivs");
    }
  };



  class FemMembraneSetGrowth : public fem::SetGrowth
  {
  public:
    FemMembraneSetGrowth(const Process &proc) : SetGrowth(proc) 
    {
      setName("Model/CCF/11 Set Growth/11a Set Growth Uniform");
    }
  };

  class FemMembraneSetGrowthMorphogens : public fem::SetGrowthMorphogensFaces
  {
  public:
    FemMembraneSetGrowthMorphogens(const Process &proc) : SetGrowthMorphogensFaces(proc) 
    {
      setName("Model/CCF/11 Set Growth/11b Set Growth Morphogens (on Faces)");
    }
  };

  class FemMembraneCellFacesGrowth : public fem::CellsToFacesGrowth
  {
  public:
    FemMembraneCellFacesGrowth(const Process &proc) : CellsToFacesGrowth(proc) 
    {
      setName("Model/CCF/12 Cells to Faces Growth");
    }
  };

  class FemMembraneGrowth : public fem::Grow
  {
  public:
    FemMembraneGrowth(const Process &proc) : Grow(proc) 
    {
      setName("Model/CCF/13 Growth");
    }
  };

  class FemMembraneComputeCurrentDirections : public fem::ComputeCurrentDirections
  {
    public:
      FemMembraneComputeCurrentDirections(const Process &proc) : ComputeCurrentDirections(proc)
      {
        setName("Model/CCF/17 Compute Current Directions") ;
      }
  };

  class FemMembraneVisDirections : public fem::VisDirections
  {
    public:
      FemMembraneVisDirections(const Process &proc) : VisDirections(proc)
      {
        setName("Model/CCF/18 Visualize Directions") ;
      }
  };

   class FemMembraneSetDirichlet : public fem::SetDirichlet
  {
  public:
    FemMembraneSetDirichlet(const Process &proc) : SetDirichlet(proc) 
    {
      setName("Model/CCF/15 Set Dirichlet");
    }
  };

  class FemMembraneDirichletDerivs : public fem::DirichletDerivs
  {
  public:
    FemMembraneDirichletDerivs(const Process &proc) : DirichletDerivs(proc) 
    {
      setName("Model/CCF/16 Dirichlet Derivs");
    }
  };
 
  class FemMembraneDiffusionAnisotropySolver : public fem::FemMorphogenSolver
  {
    public:
    FemMembraneDiffusionAnisotropySolver(const Process &proc) :FemMorphogenSolver(proc) 
    {
      setName("Model/CCF/Fem Diffusion Anisotropy/00 Diffusion Solver");

      addParm("Element Diffusion Derivs", "Process for element derivatives", "Model/CCF/Fem Diffusion Anisotropy/01 Diffusion Derivatives");
      addParm("Dirichlet Derivs", "Process for fixed concentration derivatives", "Model/CCF/Fem Diffusion Anisotropy/02 Dirichlet Derivatives");
      setParmDefault("Visualize Process", "Model/CCF/Fem Diffusion Anisotropy/05 Morphogen Visualize");
    }
  };
  class FemMembraneDiffusionAnisotropyDerivs: public fem::MorphogenDiffusionDerivs 
  {
  public:
    FemMembraneDiffusionAnisotropyDerivs(const Process &proc) : MorphogenDiffusionDerivs(proc) 
    {
       setName("Model/CCF/Fem Diffusion Anisotropy/01 Diffusion Derivatives");
       addParm("Morphogen Attribute", "Morphogen Attribute name", "AnisotropyE2Morphogen", QStringList() << "AnisotropyE2Morphogen" << "AnisotropyKParMorphogen");
    }
  };
 
  class FemMembraneDiffusionAnisotropyDirichletDerivs: public fem::MorphogenDirichletDerivs 
  {
  public:
    FemMembraneDiffusionAnisotropyDirichletDerivs(const Process &proc) : MorphogenDirichletDerivs(proc) 
    {
       setName("Model/CCF/Fem Diffusion Anisotropy/02 Dirichlet Derivatives");
       addParm("Dirichlet Attribute", "Dirichlet Attribute name", "AnisotropyE2Dirichlet", QStringList() << "AnisotropyE2Dirichlet" << "AnisotropyKParDirichlet");
       addParm("Morphogen Attribute", "Morphogen Attribute name", "AnisotropyE2Morphogen", QStringList() << "AnisotropyE2Morphogen" << "AnisotropyKParMorphogen");

    }
  };

  class FemMembraneDiffusionAnisotropySetDirichlet: public fem::SetMorphogenDirichlet 
  {
  public:
    FemMembraneDiffusionAnisotropySetDirichlet(const Process &proc) : SetMorphogenDirichlet(proc) 
    {
       setName("Model/CCF/Fem Diffusion Anisotropy/03 Set Diffusion Dirchlet");
       addParm("Morphogen Attribute", "Morphogen Attribute name", "AnisotropyE2Morphogen", QStringList() << "AnisotropyE2Morphogen" << "AnisotropyKParMorphogen");
       addParm("Dirichlet Attribute", "Dirichlet Attribute name", "AnisotropyE2Dirichlet", QStringList() << "AnisotropyE2Dirichlet" << "AnisotropyKParDirichlet");

    }
  };

  class FemMembraneCreateAnisotropyMorphogenElement : public fem::CreateMorphogenElement
  {
     public:
       FemMembraneCreateAnisotropyMorphogenElement(const Process &proc) : CreateMorphogenElement(proc)
       {
          setName("Model/CCF/Fem Diffusion Anisotropy/04 Create Diffusion Element");
       }
  };

  class FemMembraneDiffusionAnisotropyVisualize : public fem::MorphogenVisualize
  {
     public:
       FemMembraneDiffusionAnisotropyVisualize(const Process &proc) : MorphogenVisualize(proc)
       {
          setName("Model/CCF/Fem Diffusion Anisotropy/05 Morphogen Visualize");
          addParm("Morphogen Attribute", "Morphogen Attribute name", "AnisotropyE2Morphogen", QStringList() << "AnisotropyE2Morphogen" << "AnisotropyKParMorphogen");
          addParm("Morphogen Signal", "Morphogen Signal name", "E2Morphogen Signal", QStringList() << "E2Morphogen Signal" << "KParMorphogen Signal" );

       }
  };


////////////////////

  class FemMembraneDiffusionGrowthSolver : public fem::FemMorphogenSolver
  {
    public:
    FemMembraneDiffusionGrowthSolver(const Process &proc) :FemMorphogenSolver(proc) 
    {
      setName("Model/CCF/Fem Diffusion Growth/00 Diffusion Solver");

      addParm("Element Diffusion Derivs", "Process for element derivatives", "Model/CCF/Fem Diffusion Growth/01 Diffusion Derivatives");
      addParm("Dirichlet Derivs", "Process for fixed concentration derivatives", "Model/CCF/Fem Diffusion Growth/02 Dirichlet Derivatives");
      setParmDefault("Visualize Process", "Model/CCF/Fem Diffusion Growth/05 Morphogen Visualize");
    }
  };


  class FemMembraneDiffusionGrowthDerivs: public fem::MorphogenDiffusionDerivs 
  {
  public:
    FemMembraneDiffusionGrowthDerivs(const Process &proc) : MorphogenDiffusionDerivs(proc) 
    {
       setName("Model/CCF/Fem Diffusion Growth/01 Diffusion Derivatives");
       addParm("Morphogen Attribute", "Morphogen Attribute name", "GrowthMorphogen1", QStringList() << "GrowthMorphogen1" << "GrowthMorphogen2");
    }
  };
 
  class FemMembraneDiffusionGrowthDirichletDerivs: public fem::MorphogenDirichletDerivs 
  {
  public:
    FemMembraneDiffusionGrowthDirichletDerivs(const Process &proc) : MorphogenDirichletDerivs(proc) 
    {
       setName("Model/CCF/Fem Diffusion Growth/02 Dirichlet Derivatives");
       addParm("Dirichlet Attribute", "Dirichlet Attribute name", "GrowthDirichlet1", QStringList() << "GrowthDirichlet1" << "GrowthDirichlet2");
       addParm("Morphogen Attribute", "Morphogen Attribute name", "GrowthMorphogen1", QStringList() << "GrowthMorphogen1" << "GrowthMorphogen2");

    }
  };

  class FemMembraneDiffusionGrowthSetDirichlet: public fem::SetMorphogenDirichlet 
  {
  public:
    FemMembraneDiffusionGrowthSetDirichlet(const Process &proc) : SetMorphogenDirichlet(proc) 
    {
       setName("Model/CCF/Fem Diffusion Growth/03 Set Diffusion Dirchlet");
       addParm("Morphogen Attribute", "Morphogen Attribute name", "GrowthMorphogen1", QStringList() << "GrowthMorphogen1" << "GrowthMorphogen2");
       addParm("Dirichlet Attribute", "Dirichlet Attribute name", "GrowthDirichlet1", QStringList() << "GrowthDirichlet1" << "GrowthDirichlet2");

    }
  };

  class FemMembraneCreateGrowthMorphogenElement : public fem::CreateMorphogenElement
  {
     public:
       FemMembraneCreateGrowthMorphogenElement(const Process &proc) : CreateMorphogenElement(proc)
       {
          setName("Model/CCF/Fem Diffusion Growth/04 Create Diffusion Element");
       }
  };

  class FemMembraneDiffusionGrowthVisualize : public fem::MorphogenVisualize
  {
     public:
       FemMembraneDiffusionGrowthVisualize(const Process &proc) : MorphogenVisualize(proc)
       {
          setName("Model/CCF/Fem Diffusion Growth/05 Morphogen Visualize");
          addParm("Morphogen Attribute", "Morphogen Attribute name", "GrowthMorphogen1", QStringList() << "GrowthMorphogen1" << "GrowthMorphogen2");
          addParm("Morphogen Signal", "Morphogen Signal name", "Growth1 Morph Signal", QStringList() << "Growth1 Morph Signal" << "Growth2 Morph Signal" );

       }
  };



  /// Subdivide object
  class FemMembraneSubdivide : public mdx::Subdivide
  {
  public:
    FemMembraneSubdivide() {}

    FemMembraneSubdivide(Mesh &mesh, fem::ElasticTriangle3Attr &elementAttr, fem::Triangle3Attr<Point1d> &elementMorphoAnisoAttr, fem::Triangle3Attr<Point1d> &elementMorphoGrowthAttr,  
        fem::TransIsoMaterialAttr &materialAttr, fem::PressureAttr &pressureAttr, fem::PressureEdgeAttr &pressureEdgeAttr, fem::DirichletAttr &dirichletAttr, fem::GrowthAttr &growthAttr, fem::MorphogenDataAttr &morphogenE2DataAttr,
        fem::MorphogenDataAttr &morphogenKParDataAttr,  
        fem::MorphogenDataAttr &morphogenGrowth1DataAttr, fem::MorphogenDataAttr &morphogenGrowth2DataAttr, fem::MorphogenDirichletAttr &moprhogenDirichletE2Attr, 
        fem::MorphogenDirichletAttr &moprhogenDirichletKParAttr, fem::MorphogenDirichletAttr &moprhogenDirichletGrowth1Attr, fem::MorphogenDirichletAttr &moprhogenDirichletGrowth2Attr)  : 
           mdxSubdivide(mesh), elementSubdivide(mesh.indexAttr(), elementAttr), elementMorphogAnisoSubdivide(mesh.indexAttr(), elementMorphoAnisoAttr),  elementMorphogGrowthSubdivide(mesh.indexAttr(), elementMorphoGrowthAttr), 
        materialSubdivide(materialAttr), pressureSubdivide(pressureAttr), pressureEdgeSubdivide(pressureEdgeAttr), dirichletSubdivide(dirichletAttr), growthSubdivide(growthAttr), morphogenE2Subdivide(mesh.indexAttr(), morphogenE2DataAttr), morphogenKParSubdivide(mesh.indexAttr(), morphogenKParDataAttr),
        morphogenGrowth1Subdivide(mesh.indexAttr(), morphogenGrowth1DataAttr), morphogenGrowth2Subdivide(mesh.indexAttr(), morphogenGrowth2DataAttr), morphDirE2Subdivide(mesh.indexAttr(), moprhogenDirichletE2Attr), morphDirKParSubdivide(mesh.indexAttr(), moprhogenDirichletKParAttr), morphDirGrowth1Subdivide(mesh.indexAttr(), moprhogenDirichletGrowth1Attr), morphDirGrowth2Subdivide(mesh.indexAttr(), moprhogenDirichletGrowth2Attr) {}

    // Method to split the element data
    void splitCellUpdate(Dimension dim, const CCStructure &cs, const CCStructure::SplitStruct &ss, 
        CCIndex otherP = CCIndex(), CCIndex otherN = CCIndex(), double interpPos = 0.5) 
    {
      mdxSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);

      // Propagate the material parameters
      if(dim == 2 or dim == 3) {
        elementSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        materialSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        pressureSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        pressureEdgeSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);

        growthSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);

        elementMorphogAnisoSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        elementMorphogGrowthSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);

      
        morphogenE2Subdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        morphogenKParSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        morphogenGrowth1Subdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        morphogenGrowth2Subdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        
      }
      else if(dim == 1)
      {
        dirichletSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        morphDirKParSubdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        morphDirE2Subdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);

        morphDirGrowth1Subdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
        morphDirGrowth2Subdivide.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
       
      }
    }
    MDXSubdivide mdxSubdivide;
    fem::ElasticTriangle3::Subdivide elementSubdivide;
    fem::Triangle3<Point1d>::Subdivide elementMorphogAnisoSubdivide;
    fem::Triangle3<Point1d>::Subdivide elementMorphogGrowthSubdivide;

    fem::TransIsoMaterial::Subdivide materialSubdivide;
    fem::Pressure::Subdivide pressureSubdivide;
    fem::PressureEdge::Subdivide pressureEdgeSubdivide;
    fem::Growth::Subdivide growthSubdivide;
    fem::Dirichlet::Subdivide dirichletSubdivide;

    fem::MorphogenData::Subdivide morphogenE2Subdivide;
    fem::MorphogenData::Subdivide morphogenKParSubdivide;
    fem::MorphogenDirichlet::Subdivide morphDirE2Subdivide;
    fem::MorphogenDirichlet::Subdivide morphDirKParSubdivide;

    fem::MorphogenData::Subdivide morphogenGrowth1Subdivide;
    fem::MorphogenDirichlet::Subdivide morphDirGrowth1Subdivide;
    fem::MorphogenData::Subdivide morphogenGrowth2Subdivide;
    fem::MorphogenDirichlet::Subdivide morphDirGrowth2Subdivide;


  };

  class FemMembraneBisect : public SubdivideBisectTriangle
  {
  public:
    FemMembraneBisect(const Process &proc) : SubdivideBisectTriangle(proc) 
    {
      setName("Model/CCF/14 Subdivide");

      addParm("Element Attribute", "Attribute to store nodal values", "Triangle Element");
      
      addParm("Material Attribute", "Name of the attribute that holds material properties", "TransIso Material");
      addParm("Pressure Attribute", "Name of the attribute that holds pressure", "Fem Pressure");
      addParm("Pressure Edge Attribute", "Name of the attribute that holds pressure", "Fem Edge Pressure");

      addParm("Growth Attribute", "Name of the attribute that holds growth", "Fem Growth");
      addParm("Dirichlet Attribute", "Name of the attribute that holds dirichlet", "Fem Dirichlet");
      addParm("Morphogen Anisotropy Element Attribute", "Attribute to store morphogenetic element for anisotropy", "Morphogen Triangle Element Growth");
      addParm("Morphogen Growth Element Attribute", "Attribute to store morphogenetic element for growth", "Morphogen Triangle Element Growth");

      addParm("Morphogen Anisotropy E2 Data", "Attribute to store morphogenetic concentration for E2", "AnisotropyE2Morphogen");
      addParm("Morphogen Anisotropy KPar Data", "Attribute to store morphogenetic concentration for KPar", "AnisotropyKParMorphogen");
      addParm("Morphogen Anisotropy E2 Dirichlet", "Attribute to store morphogenetic Dirichlet for E2", "AnisotropyE2Dirichlet");
      addParm("Morphogen Anisotropy KPar Dirichlet", "Attribute to store morphogenetic Dirichlet for KPar", "AnisotropyKParDirichelt");

      addParm("Morphogen Growth Data 1", "Attribute to store morphogenetic concentration for growth1", "GrowthMorphogen1");
      addParm("Morphogen Growth Data 2", "Attribute to store morphogenetic concentration for growth2", "GrowthMorphogen2");
      addParm("Morphogen Growth Dirichlet 1", "Attribute to store morphogenetic Dirichlet for growth1", "GrowthDirichlet1");
      addParm("Morphogen Growth Dirichlet 2", "Attribute to store morphogenetic Dirichlet for growth2", "GrowthDirichlet2");
  
    }
    using SubdivideBisectTriangle::run;
    bool run();
  };

  class FemMembraneVisMaterial : public fem::VisTransIsoMaterial
  {
  public:
    FemMembraneVisMaterial(const Process &proc) : VisTransIsoMaterial(proc) 
    {
      setName("Model/CCF/20 Visualize Material");
    }
  };

  class FemMembraneVisGrowth : public fem::VisGrowth
  {
  public:
    FemMembraneVisGrowth(const Process &proc) : VisGrowth(proc) 
    {
      setName("Model/CCF/21 Visualize Growth");
    }
  };


  //class FemMembraneRender : public fem::VectorRender
  //{
  //public:
  //  FemMembraneRender(const Process &proc) : VectorRender(proc) 
  //  {
  //    setName("Model/CCF/22 Vector Render");
  //  }
  //};
  class FemAnisotropyPropagationFailure : public fem::DisplayFailedAnisotropyPropagation
  { 
    public:
    FemAnisotropyPropagationFailure(const Process &proc) : DisplayFailedAnisotropyPropagation(proc) 
    {
      setName("Model/CCF/22 Display Anisotropy Propagation Failure");
    }
  };
 
 
  
}
#endif

