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
#include <MeshProcessSelection.hpp>
#include <MeshUtils.hpp>
#include <QApplication>

#include <algorithm>    // std::find

namespace mdx 
{
  bool MeshSelectAll::run(const CCStructure &cs, CCIndexDataAttr &indexAttr)
  {
    for(Dimension dim = 0; dim <= cs.maxDimension; dim++) {
      // Skip edges for now
      if(dim == 1)
        continue;
      const CCIndexVec &cells = cs.cellsOfDimension(dim);
      #pragma omp parallel for
      for(uint i = 0; i < cells.size(); i++)
        indexAttr[cells[i]].selected = true;
    }
    return true;
  }
  REGISTER_PROCESS(MeshSelectAll);

  bool MeshClearSelection::run(const CCStructure &cs, CCIndexDataAttr &indexAttr)
  {
    for(Dimension dim = 0; dim <= cs.maxDimension; dim++) {
      // Skip edges for now
      if(dim == 1)
        continue;
      const CCIndexVec &cells = cs.cellsOfDimension(dim);
      #pragma omp parallel for
      for(uint i = 0; i < cells.size(); i++)
        indexAttr[cells[i]].selected = false;
    }
    return true;
  }
  REGISTER_PROCESS(MeshClearSelection);
  
//  bool MeshSelectBadNormals::run(Mesh* m)
//  {
//    int count = 0;
//    vvGraph &S = m->graph();
//    forall(const vertex& v, S) {
//      v->selected = !setNormal(S, v);
//      if(v->selected)
//        count++;
//    }
//    setStatus("Selected " << count << " vertices with bad normals");
//    m->updateSelection();
//    return true;0!8!-
//  }
//  REGISTER_PROCESS(MeshSelectBadNormals);
//  
//  bool MeshClearSelection::run(Mesh* m)
//  {
//    forall(const vertex& v, m->graph())
//      v->selected = false;
//    m->updateSelection();
//    return true;
//  }
//  REGISTER_PROCESS(MeshClearSelection);
  
  bool MeshInvertSelection::run(const CCStructure &cs, CCIndexDataAttr &indexAttr, const IntVec &dimensions)
  {
    int maxDim = cs.maxDimension;
    for(int dim : dimensions) {
      // Skip edges for now
      if(dim < 0 or dim == 1 or dim > maxDim)
        continue;
      const CCIndexVec &cells = cs.cellsOfDimension(dim);
      #pragma omp parallel for
      for(uint i = 0; i < cells.size(); i++) {
        CCIndexData &cIdx = indexAttr[cells[i]];
        cIdx.selected = not cIdx.selected;
      }
    }
    return true;
  }
  REGISTER_PROCESS(MeshInvertSelection);
   
  bool MeshCopySelection::run(CCIndexDataAttr &indexAttr)
  {
    QString indexList;
    for(auto &pr : indexAttr)
      if(pr.second.selected)
        indexList += QString(" %1").arg(pr.first.value);

    QClipboard *clipboard = QGuiApplication::clipboard();
    if(!indexList.isEmpty())
      clipboard->setText(indexList.remove(0,1));
    return true;
  }
  REGISTER_PROCESS(MeshCopySelection);

  bool MeshPasteSelection::run(CCIndexDataAttr &indexAttr)
  {
    QClipboard *clipboard = QGuiApplication::clipboard();
    QString indices = clipboard->text();
    if(indices.isEmpty())
      return false;
    QStringList indexList = indices.split(QRegExp("\\s+"), QString::SkipEmptyParts);

    for(QString &s : indexList) {
      bool ok = false;
      int idx = s.toInt(&ok);
      if(!ok or idx <= 0)
        continue;
      CCIndex c(idx);
      if(c.isPseudocell())
        continue;
      auto iter = indexAttr.find(c);
      if(iter == indexAttr.end())
        continue;
      iter->second.selected = true;
    }
    return true;
  }
  REGISTER_PROCESS(MeshPasteSelection); 

  bool MeshSaveSelection::run(CCIndexDataAttr &indexAttr, const QString &fileName)
  {
    QFile file(fileName);
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text))
      throw QString("%1::run Cannot open file: %2").arg(name()).arg(fileName);

    QTextStream out(&file);
    for(const auto &pr : indexAttr)
      if(pr.second.selected)
        out << pr.first.value << endl;

    file.close();

    return true;
  }
  REGISTER_PROCESS(MeshSaveSelection);

  bool MeshLoadSelection::run(CCIndexDataAttr &indexAttr, const QString &fileName)
  {
    QFile file(fileName);
    if(!file.open(QIODevice::ReadOnly))
      throw QString("%1::run Cannot open file: %2").arg(name()).arg(fileName);

    QString indices = QString::fromUtf8(file.readAll());
    QStringList indexList = indices.split(QRegExp("\\s+"), QString::SkipEmptyParts);

    for(QString &s : indexList) {
      bool ok = false;
      int idx = s.toInt(&ok);
      if(!ok or idx <= 0)
        continue;
      CCIndex c(idx);
      if(c.isPseudocell())
        continue;
      auto iter = indexAttr.find(c);
      if(iter == indexAttr.end())
        continue;
      iter->second.selected = true;
    }
    file.close();

    return true;
  }
  REGISTER_PROCESS(MeshLoadSelection);

  bool MeshSelectVerticesOfFaces::run(const CCStructure &cs, CCDrawParms &cdp, CCIndexDataAttr &indexAttr)
  {
    for(const CCIndex f : cs.faces()) {
      if(!indexAttr[f].selected) continue;

      std::set<CCIndex> incVtx = cs.incidentCells(f, 0);
      for(const CCIndex v : incVtx) {
        auto &vIdx = indexAttr[v];
        if(!vIdx.selected) {
          vIdx.selected = true;
          cdp.vertexChanged.push_back(v);
        }
      }
    }
    return true;
  }
  REGISTER_PROCESS(MeshSelectVerticesOfFaces);

  bool MeshSelectIncidentCells::run(const CCStructure &cs, CCIndexDataAttr &indexAttr)
  {
    for(const CCIndex f : selectedFaces(cs, indexAttr))
      for(const CCIndex v : cs.incidentCells(f, 0))
        indexAttr[v].selected = true;

    if(cs.maxDimension > 2) {
      for(const CCIndex l : selectedVolumes(cs, indexAttr)) {
        for(const CCIndex v : cs.incidentCells(l, 0))
          indexAttr[v].selected = true;
        for(const CCIndex f : cs.incidentCells(l, 2))
          indexAttr[f].selected = true;
      }
    }
    return true;
  }
  REGISTER_PROCESS(MeshSelectIncidentCells);

  bool MeshSelectFacesOfVertices::run(const CCStructure &cs, CCDrawParms &cdp, CCIndexDataAttr &indexAttr)
  {
    for(const CCIndex v : cs.vertices()) {
      if(!indexAttr[v].selected) 
        continue;

      std::set<CCIndex> incFaces = cs.incidentCells(v, 2);
      for(const CCIndex f : incFaces) {
        auto &fIdx = indexAttr[f];
        if(!fIdx.selected) {
          fIdx.selected = true;
          cdp.faceChanged.push_back(f);
        }
      }
    }

    return true;
  }
  REGISTER_PROCESS(MeshSelectFacesOfVertices);

  bool MeshSelectFacesByLabel::run(const CCStructure &cs, CCDrawParms &cdp, 
                                       CCIndexDataAttr &indexAttr, const IntIntAttr &labelMap, const IntSet &labelSet)
  {
    bool doLabel = false;
    if(labelMap.size() > 0)
      doLabel = true;
    for(const CCIndex f : cs.faces()) {
      auto &fIdx = indexAttr[f];
      int label = fIdx.label;
      if(doLabel)
        label = labelMap[label];
      if(labelSet.count(label) == 0)
        continue;
      if(!fIdx.selected) {
        fIdx.selected = true;
        cdp.faceChanged.push_back(f);
      }
    }

    return true;
  }
  REGISTER_PROCESS(MeshSelectFacesByLabel);

  bool MeshSelectFacesByAngle::run(const CCStructure &cs, CCIndexDataAttr &indexAttr, const Point3d &dir, double angle)
  {
    for(const CCIndex f : cs.faces()) {
      auto &fIdx = indexAttr[f];
      double c = normalized(dir) * normalized(fIdx.nrml);
      double a = acos(c) * 180 / M_PI;
      if(a < angle)
        fIdx.selected = true;
    }
    return true;
  }
  REGISTER_PROCESS(MeshSelectFacesByAngle);

  REGISTER_PROCESS(MeshSelectCellsByIndex);

  bool MeshSelectBorder::run(const CCStructure &cs, CCIndexDataAttr &indexAttr, const IntVec &dimensions)
  {
    for(uint dim : dimensions) {
      const std::vector<CCIndex> &cells = cs.cellsOfDimension(dim);
      #pragma omp parallel for
      for(uint i = 0 ; i < cells.size() ; i++) {
        CCIndex cell = cells[i];
        indexAttr[cell].selected = cs.onBorder(cell);
      }
    }
    return true;
  }
  REGISTER_PROCESS(MeshSelectBorder);

  bool MeshSelectBorderVertices::run(const CCStructure &cs, CCIndexDataAttr &indexAttr, double distance)
  {
    // Find vertex on border
    CCIndex v0;
    for(CCIndex v : cs.vertices())
      if(cs.onBorder(v)) {
        v0 = v;
        break;
      }
    if(v0 == CCIndex::UNDEF)
      return false;

    if(distance < 0)
      distance = 0;

    CellTuple tuple(cs, v0);
    Point3d lastPos = indexAttr[v0].pos;
    indexAttr[v0].selected = true;

    while(true) {
      while(!cs.onBorder(tuple[1]))
        tuple.flip(2,1);
        
      tuple.flip(0);
      CCIndex v = tuple[0];
      if(v == v0)
        break;
      Point3d pos = indexAttr[v].pos;
mdxInfo << "Pos:" << pos << endl;
      if(norm(pos - lastPos) >= distance) {
        indexAttr[v].selected = true;
        lastPos = pos;
mdxInfo << "Selecting pos:" << pos << endl;
      }
      tuple.flip(1);
    }
    
    return true;
  }
  REGISTER_PROCESS(MeshSelectBorderVertices);

//  bool MeshSelectUnlabeled::run(Mesh* m, bool replace)
//  {
//    forall(const vertex& v, m->graph())
//      if(v->label == 0)
//        v->selected = true;
//      else if(replace and v->label != 0)
//        v->selected = false;
//
//    m->correctSelection(true);
//    m->updateSelection();
//
//    setStatus(QString("%1 vertices selected").arg(m->selectedCount()));
//    return true;
//  }
//  REGISTER_PROCESS(MeshSelectUnlabeled);
//  
//  bool MeshSelectLabeled::run(Mesh* m, bool replace)
//  {
//    forall(const vertex& v, m->graph()) {
//      if(v->label > 0)
//        v->selected = true;
//      else if(replace and v->label == 0)
//        v->selected = false;
//    }
//    m->correctSelection(true);
//    m->updateSelection();
//
//    setStatus(QString("%1 vertices selected").arg(m->selectedCount()));
//    return true;
//  }
//  REGISTER_PROCESS(MeshSelectLabeled);
//  
//  bool MeshSelectLabel::run(Mesh* m, bool replace, int label)
//  {
//    if(label <= 0)
//      label = selectedLabel();
//    if(label <= 0)
//      throw(QString("Cannot select label, no current label is defined"));
//  
//    forall(const vertex& v, m->graph()) {
//      if(v->label == label)
//        v->selected = true;
//      else if(replace and v->label != label)
//        v->selected = false;
//    }
//
//    m->correctSelection(true);
//    m->updateSelection();
//
//    setStatus(QString("Selected label %1, vertices selected: %2").arg(label).arg(m->selectedCount()));
//    return true;
//  }
//  REGISTER_PROCESS(MeshSelectLabel);
//
//  bool MeshSelectValence::run(Mesh* m, int start, int end)
//  {
//    vvGraph &S = m->graph();
//    forall(const vertex& v, m->graph()) {
//      if(S.valence(v) >= start and S.valence(v) <= end)
//        v->selected = true;
//    }
//  
//    m->correctSelection(true);
//    m->updateSelection();
//
//    setStatus(QString("%1 vertices selected").arg(m->selectedCount()));
//    return true;
//  }
//  REGISTER_PROCESS(MeshSelectValence); 
//
//  bool MeshUnselectLabel::run(Mesh* m, int label)
//  {
//    if(label <= 0)
//      label = selectedLabel();
//    if(label <= 0)
//      throw(QString("Cannot unselect label, no current label is defined"));
//  
//    forall(const vertex& v, m->graph()) {
//      if(v->label == label)
//        v->selected = false;
//    }
//  
//    m->correctSelection(true);
//  
//    m->updateSelection();
//
//    setStatus(QString("Unselected label %1, vertices selected: %2").arg(label).arg(m->selectedCount()));
//    return true;
//  }
//  REGISTER_PROCESS(MeshUnselectLabel);
//  
//  bool MeshSelectClip::run(Mesh* m)
//  {
//    forall(const vertex& v, m->graph()) {
//      const Stack* s = m->stack();
//      bool clipped = false;
//      Point3f p = Point3f(s->frame().inverseCoordinatesOf(qglviewer::Vec(v->pos)));
//      if(clip1()->isClipped(p))
//        clipped = true;
//      if(clip2()->isClipped(p))
//        clipped = true;
//      if(clip3()->isClipped(p))
//        clipped = true;
//      if(!clipped)
//        v->selected = true;
//    }
//    m->updateSelection();
//
//    setStatus(QString("%1 vertices selected").arg(m->selectedCount()));
//    return true;
//  }
//  REGISTER_PROCESS(MeshSelectClip);
//  
//  bool MeshSelectWholeLabelExtend::run(Mesh* m)
//  {
//    std::set<int> labels;
//  
//    forall(const vertex& v, m->graph())
//      if(v->selected and v->label > 0)
//        labels.insert(v->label);
//  
//    forall(const vertex& v, m->graph())
//      if(labels.find(v->label) != labels.end())
//        v->selected = true;
//  
//    m->correctSelection(true);
//    m->updateSelection();
//
//    setStatus(QString("%1 vertices selected").arg(m->selectedCount()));
//    return true;
//  }
//  REGISTER_PROCESS(MeshSelectWholeLabelExtend);
//  
//  bool MeshSelectDuplicateCells::run(Mesh* m)
//  {
//    // Labels to select
//    std::map<int, int> LabCount;
//    std::map<int, vertex> Labels;
//    std::set<vertex> Vertices;
//  
//    // Grab one vertex for each label
//    vvGraph& S = m->graph();
//    forall(const vertex& v, S) {
//      Vertices.insert(v);
//      Labels[v->label] = v;
//    }
//  
//    // Find contigs
//    while(!Vertices.empty()) {
//      // Grab any vertex and save the label
//      vertex v = *Vertices.begin();
//      Vertices.erase(v);
//      int label = v->label;
//  
//      // Start growing neighbor set
//      std::set<vertex> Nbs;
//      std::set<vertex> NewNbs;
//      Nbs.insert(v);
//      do {
//        forall(const vertex& u, Nbs)
//          forall(const vertex& n, S.neighbors(u))
//            if(Vertices.count(n) > 0 and n->label == label) {
//              Vertices.erase(n);
//              NewNbs.insert(n);
//            }
//        Nbs = NewNbs;
//        NewNbs.clear();
//      } while(!Nbs.empty());
//  
//      // One more region for this label
//      LabCount[label]++;
//    }
//  
//    // Mark labels with more than one region selected
//    forall(const vertex& v, S)
//      if(LabCount[v->label] > 1)
//        v->selected = true;
//      else
//        v->selected = false;
//  
//    m->correctSelection(true);
//    m->updateSelection();
//
//    setStatus(QString("%1 vertices selected").arg(m->selectedCount()));
//    return true;
//  }
//  REGISTER_PROCESS(MeshSelectDuplicateCells);
//  
//  bool ExtendByConnectivity::run(Mesh* m)
//  {
//    vvGraph& S = m->graph();
//    std::vector<vertex> vs = m->activeVertices();
//    // Either all vertices are selected, or none
//    if(vs.size() == S.size())
//      return true;
//    std::set<vertex> selected(vs.begin(), vs.end());
//    // Tabular approach
//    for(size_t i = 0; i < vs.size(); ++i) {
//      vertex v = vs[i];
//      forall(const vertex& n, S.neighbors(v)) {
//        if(selected.count(n))
//          continue;
//        n->selected = true;
//        selected.insert(n);
//        vs.push_back(n);
//      }
//    }
//    m->updateSelection();
//
//    setStatus(QString("%1 vertices selected").arg(m->selectedCount()));
//    return true;
//  }
//  REGISTER_PROCESS(ExtendByConnectivity);
//
//  bool SelectByNormal::run(Mesh *m, double tolerance)
//  {
//    vvGraph& S = m->graph(); 
//
//    Point3d nrml (0,0,0);
//    int counter = 0;
//
//    // take average normal selected vertices
//    forall(const vertex& v, S){
//      if(v->selected){
//        nrml += v->nrml;
//        counter++;
//      }
//    }
//
//    nrml /= counter;
//    nrml /= norm(nrml);
//
//    forall(const vertex& v, S){
//      double dis = v->nrml * nrml; // 1 for identical, -1 for opposite
//      if(dis > 1-tolerance) v->selected = true;
//    }
//
//    m->updateAll();
//
//    setStatus(QString("%1 vertices selected").arg(m->selectedCount()));
//    return true;
//  }
//  REGISTER_PROCESS(SelectByNormal);
//
//  bool SelectSharedTriangles::run(Mesh* m)
//  {
//	  if(m->meshType() != "MDX3D") 
//		  throw(QString("Mesh type (%1) doesn't have shared triangles, mesh type must be (MDX3D)").arg(m->meshType()));
//	  const cellGraph &Cells = m->cells();
//
//    // first unselect everything
//    forall(const vertex& v, m->graph())
//      v->selected = false;
//    m->updateSelection();
//
//    // list of triangle, to find duplicates 
//    std::set<Triangle> triList; 
//	  std::pair<std::set<Triangle>::iterator,bool> newTri, newTriRot;
//
//    // look which triangles belong to more than one cell
//	  forall(const cell &c, Cells)
//      forall(const vertex &v, c->S)
//		    forall(const vertex &n, c->S.neighbors(v)) {
//				  vertex m = c->S.nextTo(v, n);
//			    if(!c->S.uniqueTri(v, n, m))
//				    continue;
//          
//				  newTri = triList.insert(Triangle(v,n,m));
//					// triangle(v,n,m) in one cells will be oriented (v,m,n) in the other cell  
//				  newTriRot = triList.insert(Triangle(v,m,n));
//					// if triangle already exist in another cell, select its vertices
//					if(newTri.second == false or newTriRot.second == false) {
//					  n->selected = true;
//					  v->selected = true;
//					  m->selected = true;
//					}
//				}
//    m->updateSelection();
//
//    setStatus(QString("%1 vertices selected").arg(m->selectedCount()));
//    return true;
//  }
//  REGISTER_PROCESS(SelectSharedTriangles);
//
//  bool AreaSelectedTris::run(Mesh* m, QString mode)
//  {
//
//    double area = 0;
//
//    const std::vector<vertex>& vs = m->selectedVertices();
//    vvGraph& S = m->graph();
//
//    if(mode == "Tris inside Vtxs"){
//      forall(const vertex& v, vs){
//        forall(const vertex& n, S.neighbors(v)){
//          if(!n->selected) continue;
//          vertex m = S.nextTo(v,n);
//          if(!S.uniqueTri(v,n,m) or !m->selected) continue;
//
//          area += triangleArea(v->pos, n->pos, m->pos);
//        }
//      }
//
//    } else if(mode == "Tris neighboring Vtxs"){
//      forall(const vertex& v, S){
//        forall(const vertex& n, S.neighbors(v)){
//          vertex m = S.nextTo(v,n);
//          if(!v->selected and !n->selected and !m->selected) continue;
//          if(!S.uniqueTri(v,n,m)) continue;
//
//          area += triangleArea(v->pos, n->pos, m->pos);
//        }
//      }
//    }
//
//
//    setStatus(QString("Area of selected triangles: %1").arg(area));
//    return true;
//  }
//  REGISTER_PROCESS(AreaSelectedTris);

}
