//
// This file is part of MorphoDynamX - http://www.MorphoDynamX.org
// Copyright (C) 2012-2015 Richard S. Smith and collaborators.
//
// If you use MorphoDynamX in your work, please cite:
//   http://dx.doi.org/10.7554/eLife.05864
//
// MorphoDynamX is free software, and is licensed under under the terms of the 
// GNU General (GPL) Public License version 2.0, http://www.gnu.org/licenses.
// 
#ifndef MESH_PROCESS_SELECTION_HPP
#define MESH_PROCESS_SELECTION_HPP

#include <Process.hpp>

namespace mdx
{
  ///\addtogroup MeshProcess
  ///@{
  /**
   * \class MeshSelectAll
   *
   * Select all the vertices in the current mesh.
   */
  class mdxBase_EXPORT MeshSelectAll : public Process 
  {
  public:
    MeshSelectAll(const Process& process) : Process(process) 
    {
      setName("Mesh/Selection/Select All");
      setDesc("Select all vertices and faces of the current mesh");
    }
  
    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh) 
        throw(QString("MeshSelectAll::run No current mesh"));
      QString ccName = mesh->ccName();
      if(ccName.isEmpty())
        throw(QString("MeshSelectAll::run No cell complex"));

      CCStructure &cs = mesh->ccStructure(ccName);
      CCIndexDataAttr &indexAttr = mesh->indexAttr();
      mesh->updateProperties();

      return run(cs, indexAttr);
    }
    bool run(const CCStructure &cs, CCIndexDataAttr &indexAttr);
  
  };
  

//   * \class MeshSelectBadNormals
//   *
//   * Select all the vertices with invalid normals in the current mesh
//   */
//  class mdxBase_EXPORT MeshSelectBadNormals : public Process 
//  {
//  public:
//    MeshSelectBadNormals(const Process& process) : Process(process) {}
//  
//    bool run(const QStringList &parms)
//    {
//      Mesh *mesh = currentMesh();
//      if(!mesh) throw(QString("No current mesh"));
//      return run(mesh);
//    }
//  
//    bool run(Mesh* m);
//  
//    QString name() const { return "Mesh/Selection/Select Bad Normals"; }
//    QString description() const { return 
//      "Select all vertices of the current mesh with normals that cannot be calculated"; }
//    QStringList parmNames() const { return QStringList(); }
//    QStringList parmDescs() const { return QStringList(); }
//  };
//  
  /**
   * \class MeshClearSelection
   *
   * Ensure no vertex is selected in the current mesh.
   */
  class mdxBase_EXPORT MeshClearSelection : public Process 
  {
  public:
    MeshClearSelection(const Process& process) : Process(process) {
      setName("Mesh/Selection/Clear Selection");
      setDesc("Unselect all vertices and faces of the current mesh");
    }
  
    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh) 
        throw(QString("MeshClearSelection::run No current mesh"));
      QString ccName = mesh->ccName();
      if(ccName.isEmpty())
        throw(QString("MeshClearSelection::run No cell complex"));

      CCStructure &cs = mesh->ccStructure(ccName);
      CCIndexDataAttr &indexAttr = mesh->indexAttr();
      mesh->updateProperties();

      return run(cs, indexAttr);
    }
    bool run(const CCStructure &cs, CCIndexDataAttr &indexAttr);
  };

  /**
   * \class MeshInvertSelection
   *
   * Invert the selection status of all the vertices in the current mesh.
   */
  class mdxBase_EXPORT MeshInvertSelection : public Process 
  {
  public:
    MeshInvertSelection(const Process& process) : Process(process) 
    {
      setName("Mesh/Selection/Invert Selection");
      setDesc("Invert the selection");

      addParm("Dimension", "Which dimension to invert", "Faces", dimChoice());
    }
  
    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh) 
        throw QString("%1::run No current mesh").arg(name());
      QString ccName = mesh->ccName();
      if(ccName.isEmpty())
        throw QString("%1:run No cell complex").arg(name());

      CCStructure &cs = mesh->ccStructure(ccName);
      CCIndexDataAttr &indexAttr = mesh->indexAttr();
      mesh->updateProperties();

      IntVec dimension;
      QString dimName = parm("Dimension");
      if(dimName == "Vertices")
        dimension.push_back(0);
      else if(dimName == "Edges")
        throw QString("%1:run Selection not stored on edges").arg(name());
      else if(dimName == "Faces")
        dimension.push_back(2);
      else if(dimName == "Volumes")
        dimension.push_back(3);
      else if(dimName == "All") {
        dimension.push_back(0);
        dimension.push_back(2);
        dimension.push_back(3);
      } else
        throw QString("%1:run Unknown dimension (%2)").arg(name()).arg(dimName);

      return run(cs, indexAttr, dimension);
    }
    bool run(const CCStructure &cs, CCIndexDataAttr &indexAttr, const IntVec &dimensions = {2}); // Default to faces
  };
 
  /**
   * \class MeshCopySelection
   *
   * Copy the selected cells to the clipboard
   */
  class mdxBase_EXPORT MeshCopySelection : public Process 
  {
  public:
    MeshCopySelection(const Process& process) : Process(process) 
    {
      setName("Mesh/Selection/Copy Selection");
      setDesc("Copy the selection to the clipboard");
    }
  
    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh) 
        throw QString("%1::run No current mesh").arg(name());

      CCIndexDataAttr &indexAttr = mesh->indexAttr();

      return run(indexAttr);
    }
    bool run(CCIndexDataAttr &indexAttr);
  }; 

  /**
   * \class MeshPasteSelection
   *
   * Select the cells that are in the clipboard
   */
  class mdxBase_EXPORT MeshPasteSelection : public Process 
  {
  public:
    MeshPasteSelection(const Process& process) : Process(process) 
    {
      setName("Mesh/Selection/Paste Selection");
      setDesc("Select the cells that are in the clipboard");
    }
  
    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh) 
        throw QString("%1::run No current mesh").arg(name());

      CCIndexDataAttr &indexAttr = mesh->indexAttr();
      mesh->updateProperties();

      return run(indexAttr);
    }
    bool run(CCIndexDataAttr &indexAttr);
  }; 

  /**
   * \class MeshSaveSelection
   *
   * Save the selected cells to a file
   */
  class mdxBase_EXPORT MeshSaveSelection : public Process 
  {
  public:
    MeshSaveSelection(const Process& process) : Process(process) 
    {
      setName("Mesh/Selection/Save Selection");
      setDesc("Save the selection to a file");

      addParm("File Name", "File to write selection to", "Selection.txt");
    }
  
    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh) 
        throw QString("%1::run No current mesh").arg(name());

      const QString fileName = parm("File Name");
      if(fileName.isEmpty())
        throw QString("%1::run File name is empty").arg(name());

      CCIndexDataAttr &indexAttr = mesh->indexAttr();

      return run(indexAttr, fileName);
    }
    bool run(CCIndexDataAttr &indexAttr, const QString &fileName);
  }; 

  /**
   * \class MeshLoadSelection
   *
   * Save the selected cells to a file
   */
  class mdxBase_EXPORT MeshLoadSelection : public Process 
  {
  public:
    MeshLoadSelection(const Process& process) : Process(process) 
    {
      setName("Mesh/Selection/Load Selection");
      setDesc("Save the selection to a file");

      addParm("File Name", "File to load selection from", "Selection.txt");
    }
  
    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh) 
        throw QString("%1::run No current mesh").arg(name());

      const QString fileName = parm("File Name");
      if(fileName.isEmpty())
        throw QString("%1::run File name is empty").arg(name());

      CCIndexDataAttr &indexAttr = mesh->indexAttr();
      mesh->updateProperties();

      return run(indexAttr, fileName);
    }
    bool run(CCIndexDataAttr &indexAttr, const QString &fileName);
  }; 

  /**
   * \class MeshSelectVerticesOfFaces
   *
   * Select incident vertices of selected faces
   */
  class mdxBase_EXPORT MeshSelectVerticesOfFaces : public Process 
  {
  public:
    MeshSelectVerticesOfFaces(const Process& process) : Process(process) 
    {
      setName("Mesh/Selection/Select Vertices of Faces");
      setDesc("Select incident vertices of selected faces");
    }
  
    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh) 
        throw QString("%1::run No current mesh").arg(name());
      QString ccName = mesh->ccName();
      if(ccName.isEmpty())
        throw QString("%1::run No cell complex").arg(name());


      auto &cdp = mesh->drawParms(ccName);
      auto &cs = mesh->ccStructure(ccName);
      auto &indexAttr = mesh->indexAttr();
      mesh->updateVertexSelect();

      return run(cs, cdp, indexAttr);
    }
    bool run(const CCStructure &cs, CCDrawParms &cdp, CCIndexDataAttr &indexAttr);
  };

  /**
   * \class MeshSelectIncidentCells
   *
   * Select incident cell of selected cells
   */
  class mdxBase_EXPORT MeshSelectIncidentCells : public Process 
  {
  public:
    MeshSelectIncidentCells(const Process& process) : Process(process) 
    {
      setName("Mesh/Selection/Select Incident Cells");
      setDesc("Select incident vertices of selected faces");
    }
  
    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh) 
        throw QString("%1::run No current mesh").arg(name());
      QString ccName = mesh->ccName();
      if(ccName.isEmpty())
        throw QString("%1::run No cell complex").arg(name());

      auto &cs = mesh->ccStructure(ccName);
      auto &indexAttr = mesh->indexAttr();
      mesh->updateProperties(ccName);

      return run(cs, indexAttr);
    }
    bool run(const CCStructure &cs, CCIndexDataAttr &indexAttr);
  };

  /**
   * \class MeshSelectFacesOfVertices
   *
   * Select incident faces of selected vertices
   */
  class mdxBase_EXPORT MeshSelectFacesOfVertices : public Process 
  {
  public:
    MeshSelectFacesOfVertices(const Process& process) : Process(process) 
    {
      setName("Mesh/Selection/Select Faces of Vertices");
      setDesc("Select incident faces of selected vertices");
    }
  
    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh) 
        throw QString("%1::run No current mesh").arg(name());
      QString ccName = mesh->ccName();
      if(ccName.isEmpty())
        throw QString("%1::run No cell complex").arg(name());

      auto &cdp = mesh->drawParms(ccName);
      auto &cs = mesh->ccStructure(ccName);
      auto &indexAttr = mesh->indexAttr();
      mesh->updateFaceSelect();

      return run(cs, cdp, indexAttr);
    }
    bool run(const CCStructure &cs, CCDrawParms &cdp, CCIndexDataAttr &indexAttr);
  };

  /**
   * \class MeshSelectFacesByLabel
   *
   * Select faces by label
   */
  class mdxBase_EXPORT MeshSelectFacesByLabel : public Process 
  {
  public:
    MeshSelectFacesByLabel(const Process& process) : Process(process) 
    {
      setName("Mesh/Selection/Select Faces by Label");
      setDesc("Select incident faces of selected vertices");

      addParm("Label", "Labels of faces to select, empty for current", "");
    }
  
    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh) 
        throw QString("%1::run No current mesh").arg(name());
      QString ccName = mesh->ccName();
      if(ccName.isEmpty())
        throw QString("%1::run No cell complex").arg(name());

      QStringList labels = parm("Label").split(QRegExp("\\s"));
      if(labels.size() == 0)
        throw QString("%1::run No labels specified").arg(name());

      IntSet labelSet;
      for(const QString &label : labels)
        if(!label.isEmpty())
          labelSet.insert(label.toInt());
      if(labelSet.size() == 0)
        labelSet.insert(selectedLabel());

      auto &cdp = mesh->drawParms(ccName);
      auto &cs = mesh->ccStructure(ccName);
      auto &indexAttr = mesh->indexAttr();
      mesh->updateFaceSelect();

      QString labeling = mesh->labeling();
      const auto &labelMap = labeling == "Labels" ? IntIntAttr() : mesh->labelMap(labeling);
      if(labeling != "Labels" and labelMap.size() == 0)
        throw QString("%1::run Labeling (%2) is empty").arg(name()).arg(labeling);

      return run(cs, cdp, indexAttr, labelMap, labelSet);
    }
    bool run(const CCStructure &cs, CCDrawParms &cdp, CCIndexDataAttr &indexAttr, const IntIntAttr &labelMap, const IntSet &labelSet);
  };

  /**
   * \class MeshSelecFacesByAngle
   *
   * Select faces by angle of the normal
   */
  class mdxBase_EXPORT MeshSelectFacesByAngle : public Process 
  {
  public:
    MeshSelectFacesByAngle(const Process& process) : Process(process) 
    {
      setName("Mesh/Selection/Select Faces by Angle");
      setDesc("Select faces within a given angle of a direction");

      addParm("Direction", "Direction for selection", "0.0 0.0 1.0");
      addParm("Angle", "Select faces within this angle from the direction", "15.0");
    }
  
    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh) 
        throw QString("%1::run No current mesh").arg(name());

      QString ccName = mesh->ccName();
      if(ccName.isEmpty())
        throw QString("%1::run No cell complex").arg(name());
      auto &cs = mesh->ccStructure(ccName);

      Point3d dir = stringToPoint3d(parm("Direction"));
      if(norm(dir) == 0)
        throw QString("%1::run Bad direction").arg(name());

      double angle = parm("Angle").toDouble();
      if(angle < 0 or angle > 90)
        throw QString("%1::run Angle must be between 0-90 degrees").arg(name());

      auto &indexAttr = mesh->indexAttr();
      mesh->updateProperties();
      return run(cs, indexAttr, dir, angle);
    }
    bool run(const CCStructure &cs, CCIndexDataAttr &indexAttr, const Point3d &dir, double angle);
  };

  /**
   * \class MeshSelectCellsByIndex
   *
   * Select cell by CCIndex
   */
  class mdxBase_EXPORT MeshSelectCellsByIndex : public Process 
  {
  public:
    MeshSelectCellsByIndex(const Process& process) : Process(process) 
    {
      setName("Mesh/Selection/Select Cells by Index");
      setDesc("Select cells by index (CCIndex)");

      addParm("Index", "Index(s) of cells to select", "");
    }
  
    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh) 
        throw QString("%1::run No current mesh").arg(name());
      QString ccName = mesh->ccName();
      if(ccName.isEmpty())
        throw QString("%1::run No cell complex").arg(name());

      QStringList indices = parm("Index").split(QRegExp("\\s"));
      if(indices.size() == 0)
        throw QString("%1::run No indices specified").arg(name());

      auto &indexAttr = mesh->indexAttr();

      for(const QString &index : indices)
        if(!index.isEmpty()) {
          CCIndex c(index.toInt());
          if(!c.isPseudocell())
            indexAttr[c].selected = true;
        }
      mesh->updateProperties();

      return true;
    }
    bool run(const CCStructure &cs, CCDrawParms &cdp, CCIndexDataAttr &indexAttr, const IntIntAttr &labelMap, const IntSet &labelSet);
  };

  class mdxBase_EXPORT MeshSelectBorder : public Process
  {
  public:
    MeshSelectBorder(const Process& process) : Process(process)
    {
      setName("Mesh/Selection/Select Border");
      setDesc("Select cells on the border of the cell complex");

      addParm("Dimension", "Dimension of cells to select", "Vertices", dimChoice());
    }

    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh)
        throw QString("%1::run No current mesh").arg(name());
      QString ccName = mesh->ccName();
      if(ccName.isEmpty())
        throw QString("%1:run No cell complex").arg(name());

      CCStructure &cs = mesh->ccStructure(ccName);
      CCIndexDataAttr &indexAttr = mesh->indexAttr();
      mesh->updateProperties();

      IntVec dimension;
      QString dimName = parm("Dimension");
      if(dimName == "Vertices")
        dimension.push_back(0);
      else if(dimName == "Edges")
        throw QString("%1:run Selection not stored on edges").arg(name());
      else if(dimName == "Faces")
        dimension.push_back(2);
      else if(dimName == "Volumes")
        dimension.push_back(3);
      else if(dimName == "All") {
        dimension.push_back(0);
        dimension.push_back(2);
        dimension.push_back(3);
      } else
        throw QString("%1:run Unknown dimension (%2)").arg(name()).arg(dimName);

      return run(cs, indexAttr, dimension);
    }
    bool run(const CCStructure &cs, CCIndexDataAttr &indexAttr, const IntVec &dimensions = {0});
  };

  class mdxBase_EXPORT MeshSelectBorderVertices : public Process
  {
  public:
    MeshSelectBorderVertices(const Process& process) : Process(process)
    {
      setName("Mesh/Selection/Select Border Vertices");
      setDesc("Select vertices on the border of the cell complex");

      addParm("Distance", "Distance between points to select (um), 0 for all", "1.0");
    }

    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh)
        throw QString("%1::run No current mesh").arg(name());
      QString ccName = mesh->ccName();
      if(ccName.isEmpty())
        throw QString("%1:run No cell complex").arg(name());

      CCStructure &cs = mesh->ccStructure(ccName);
      CCIndexDataAttr &indexAttr = mesh->indexAttr();
      mesh->updateProperties();
      return run(cs, indexAttr, parm("Distance").toDouble());
    }
    bool run(const CCStructure &cs, CCIndexDataAttr &indexAttr, double distancee);
  };


//  /**
//   * \class MeshSelectUnlabeled
//   *
//   * Select all the unlabel vertices of the current mesh. Other vertices may either
//   * be unselected, or stay in their current selection state.
//   */
//  class mdxBase_EXPORT MeshSelectUnlabeled : public Process {
//  public:
//    MeshSelectUnlabeled(const Process& process) : Process(process) {}
//  
//    bool run(const QStringList &parms)
//    {
//      Mesh *mesh = currentMesh();
//      if(!mesh) throw(QString("No current mesh"));
//      bool replace = stringToBool(parms[0]);
//      return run(mesh, replace);
//    }
//  
//    bool run(Mesh* m, bool replace);
//  
//    QString name() const { return "Mesh/Selection/Select Unlabeled"; }
//    QString description() const { return 
//      "Add to or replace the selection with the unlabeled vertices."; }
//    QStringList parmNames() const { return QStringList() << "Replace selection"; }
//    QStringList parmDescs() const { return QStringList() << "Replace selection"; }
//    QStringList parmDefaults() const { return QStringList() << "No"; }
////    ParmChoiceMap parmChoice() const
////    {
////      ParmChoiceMap map;
////      map[0] = booleanChoice();
////      return map;
////    }
//  };
//  
//  /**
//   * \class MeshSelectLabeled
//   *
//   * Select all labeled vertices in the current mesh.
//   */
//  class mdxBase_EXPORT MeshSelectLabeled : public Process 
//  {
//  public:
//    MeshSelectLabeled(const Process& process) : Process(process) {}
//  
//    bool run(const QStringList &parms)
//    {
//      Mesh *mesh = currentMesh();
//      if(!mesh) throw(QString("No current mesh"));
//      bool replace = stringToBool(parms[0]);
//      return run(mesh, replace);
//    }
//  
//    bool run(Mesh* m, bool replace);
//  
//    QString name() const { return "Mesh/Selection/Select Labeled"; }
//    QString description() const { return 
//      "Add to or replace the selection with the labeled vertices."; }
//    QStringList parmNames() const { return QStringList() << "Replace selection"; }
//    QStringList parmDescs() const { return QStringList() << "Replace selection"; }
//    QStringList parmDefaults() const { return QStringList() << "No"; }
////    ParmChoiceMap parmChoice() const
////    {
////      ParmChoiceMap map;
////      map[0] = booleanChoice();
////      return map;
////    }  
//  };
//  
//  /**
//   * \class MeshSelectLabel
//   */
//  class mdxBase_EXPORT MeshSelectLabel : public Process 
//  {
//  public:
//    MeshSelectLabel(const Process& process) : Process(process) {}
//  
//    bool run(const QStringList &parms)
//    {
//      Mesh *mesh = currentMesh();
//      if(!mesh) throw(QString("No current mesh"));
//      bool replace = stringToBool(parms[0]);
//      return run(mesh, replace, parms[1].toInt());
//    }
//  
//    bool run(Mesh* m, bool replace, int label);
//  
//    QString name() const { return "Mesh/Selection/Select Label"; }
//    QString description() const { return 
//      "Add to or replace the selection with the vertices of a given label "
//      "(0 for current label)."; }
//    QStringList parmNames() const { return QStringList() 
//      << "Replace selection" << "Label (0 for current)"; }
//    QStringList parmDefaults() const { return QStringList() << "No" << "0"; }
////    ParmChoiceMap parmChoice() const
////    {
////      ParmChoiceMap map;
////      map[0] = booleanChoice();
////      return map;
////    }
//  };
//
//   /**
//   * \class MeshSelectValence
//   */
//  class mdxBase_EXPORT MeshSelectValence : public Process 
//  {
//  public:
//    MeshSelectValence(const Process& process) : Process(process) {}
//  
//    bool run(const QStringList &parms)
//    {
//      Mesh *mesh = currentMesh();
//      if(!mesh) throw(QString("No current mesh"));
//      return run(mesh, parms[0].toInt(), parms[1].toInt());
//    }
//  
//    bool run(Mesh* m, int start, int end);
//  
//    QString name() const { return "Mesh/Selection/Select By Valence"; }
//    QString description() const { return "Select vertices by Valence"; }
//    QStringList parmNames() const { return QStringList() 
//      << "Begin Valence" << "End Valence"; }
//    QStringList parmDefaults() const { return QStringList() << "1" << "5"; }
//  }; 
//
//  /**
//   * \class MeshUnselectLabel ProcessSelection.hpp <MeshProcessSelection.hpp>
//   *
//   * Unselect all the vertices having a given label.
//   */
//  class mdxBase_EXPORT MeshUnselectLabel : public Process 
//  {
//  public:
//    MeshUnselectLabel(const Process& process) : Process(process) {}
//  
//    bool run(const QStringList &parms)
//    {
//      Mesh *mesh = currentMesh();
//      if(!mesh) throw(QString("No current mesh"));
//      return run(mesh, parms[0].toInt());
//    }
//  
//    bool run(Mesh* m, int label);
//  
//    QString name() const { return "Mesh/Selection/Unselect Label"; }
//    QString description() const { return 
//      "Remove the vertices of a given label (0 for current label) from the selection."; }
//    QStringList parmNames() const { return QStringList() << "Label (0 for current)"; }
//    QStringList parmDescs() const { return QStringList() << "Label (0 for current)"; }
//    QStringList parmDefaults() const { return QStringList() << "0"; }
//  };
//  
//  /**
//   * \class MeshSelectClip ProcessSelection.hpp <MeshProcessSelection.hpp>
//   *
//   * Select all vertices within the clipped region.
//   */
//  class mdxBase_EXPORT MeshSelectClip : public Process 
//  {
//  public:
//    MeshSelectClip(const Process& process) : Process(process) {}
//  
//    bool run()
//    {
//      Mesh *mesh = currentMesh();
//      if(!mesh) throw(QString("No current mesh"));
//      return run(mesh);
//    }
//  
//    bool run(Mesh* mesh);
//  
//    QString name() const { return "Mesh/Selection/Select Clip Region"; }
//    QString description() const { return "Add vertices in clip region to selection."; }
//  };
//  
//  /**
//   * \class MeshSelectWholeLabelExtend ProcessSelection.hpp <MeshProcessSelection.hpp>
//   *
//   * Extent the current selection so each label having at least one vertex
//   * selected will be fully selected.
//   */
//  class mdxBase_EXPORT MeshSelectWholeLabelExtend : public Process 
//  {
//  public:
//    MeshSelectWholeLabelExtend(const Process& process) : Process(process) {}
//  
//    bool run()
//    {
//      Mesh *mesh = currentMesh();
//      if(!mesh) throw(QString("No current mesh"));
//      return run(mesh);
//    }
//  
//    bool run(Mesh* mesh);
//  
//    QString name() const { return "Mesh/Selection/Extend to Whole Cells"; }
//    QString description() const { return "Extend Selection to Whole Cells"; }
//  };
//  
//  /**
//   * \class MeshSelectDuplicateCells ProcessSelection.hpp <MeshProcessSelection.hpp>
//   *
//   * Select vertices if the region with their label is not contiguous 
//   * (e.g. will be seen as more than one cell).
//   */
//  class mdxBase_EXPORT MeshSelectDuplicateCells : public Process 
//  {
//  public:
//    MeshSelectDuplicateCells(const Process& process) : Process(process) {}
//  
//    bool run()
//    {
//      Mesh *mesh = currentMesh();
//      if(!mesh) throw(QString("No current mesh"));
//      return run(mesh);
//    }
//  
//    bool run(Mesh* mesh);
//  
//    QString name() const { return "Mesh/Selection/Select Duplicate Cells"; }
//    QString description() const { return "Select cells with duplicate labels."; }
//  };
//  
//  /**
//   * \class ExtendByConnectivity ProcessSelection.hpp <MeshProcessSelection.hpp>
//   *
//   * Extend the current selection to all vertices that are connected to a
//   * currently selected vertex.
//   */
//  class ExtendByConnectivity : public Process 
//  {
//  public:
//    ExtendByConnectivity(const Process& process) : Process(process) {}
//  
//    bool run()
//    {
//      Mesh *mesh = currentMesh();
//      if(!mesh) throw(QString("No current mesh"));
//      return run(mesh);
//    }
//  
//    bool run(Mesh* mesh);
//  
//    QString name() const { return "Mesh/Selection/Extend by Connectivity"; }
//    QString description() const { return "Extend the selection to connected regions"; }
//    QIcon icon() const { return QIcon(":/images/SelectConnected.png"); }
//  };
//  ///@}
//
//  /**
//   * \class SelectByNormal ProcessSelection.hpp <MeshProcessSelection.hpp>
//   *
//   * Selects vertices according to how close their normal is compared to pre-selected reference vertices.
//   * (originally from Cell Maker)
//   */
//class SelectByNormal : public Process
//  {
//    public:
//      SelectByNormal(const Process& process) : Process(process) {}
//
//      bool run(const QStringList &parms)
//      {
//      Mesh *mesh = currentMesh();
//      if(!mesh) throw(QString("No current mesh"));
//        return run(mesh, parms[0].toDouble());
//      }
//
//      bool run(Mesh *mesh, double tolerance);
//
//      QString name() const { return "Mesh/Selection/Extend Selection by Normal"; }
//      QString description() const { return 
//        "Selects vertices according to how close their normal is compared to pre-selected reference vertices."; }
//      QStringList parmNames() const { return QStringList() << "Tolerance"; }
//      QStringList parmDescs() const { return QStringList() << "Tolerance threshold of the scalar of the averaged normal of selected vertices and other vertices"; }
//      QStringList parmDefaults() const { return QStringList() << "0.01"; }
////      ParmChoiceMap parmChoice() const
////      {
////        ParmChoiceMap map;
////        return map;
////      }
//      QIcon icon() const { return QIcon(":/images/SelectConnected.png"); }
//  };
//
//	 /**
//   * \class SelectSharedTriangles ProcessSelection.hpp <MeshProcessSelection.hpp>
//   *
//   * Extend the current selection to all vertices that are connected to a
//   * currently selected vertex.
//   */
//  class SelectSharedTriangles : public Process 
//  {
//  public:
//    SelectSharedTriangles(const Process& process) : Process(process) {}
//  
//    bool run(const QStringList &parms)
//    {
//      Mesh *mesh = currentMesh();
//      if(!mesh) throw(QString("No current mesh"));
//      return run(mesh);
//    }
//  
//    bool run(Mesh* mesh);
//  
//    QString name() const { return "Mesh/Selection/Select Shared Triangles"; }
//    QString description() const { return "Select triangles shared by several cells"; }
//    QStringList parmNames() const { return QStringList(); }
//    QStringList parmDescs() const { return QStringList(); }
//    QIcon icon() const { return QIcon(":/images/Triangle.png"); }
//  };
//  ///@}
//
// /**
//   * \class AreaSelectedTris MeshProcessSelection.hpp <MeshProcessSelection.hpp>
//   *
//   * Calculate the area of selected triangles
//   */
//  class AreaSelectedTris : public Process 
//  {
//  public:
//    AreaSelectedTris(const Process& process) : Process(process) {}
//  
//    bool run(const QStringList &parms)
//    {
//      Mesh *mesh = currentMesh();
//      if(!mesh) throw(QString("No current mesh"));
//      return run(mesh, parms[0]);
//    }
//  
//    bool run(Mesh* mesh, QString mode);
//  
//    QString name() const { return "Mesh/Selection/Area of Selected Triangles"; }
//    QString description() const { return "Returns the area of all selected triangles. Either triangles with all 3 vertices selected or triangles with at least one selected vertex"; }
//    QStringList parmNames() const { return QStringList() << "Mode"; }
//    QStringList parmDescs() const { return QStringList() << "Mode"; }
//    QStringList parmDefaults() const { return QStringList() << "Tris inside Vtxs"; }
//    QIcon icon() const { return QIcon(":/images/Hex.png"); }
//
////    ParmChoiceMap parmChoice() const
////    {
////      ParmChoiceMap map;
////      map[0] = QStringList() << "Tris inside Vtxs" << "Tris neighboring Vtxs";
////      return map;
////    }
//
//  };
//

}

#endif
