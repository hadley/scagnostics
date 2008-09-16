#include "GraphMeasures.h"
#include <cstdlib>
#include <cstdio>
#include <algorithm>

using namespace std;

// for using qsort
int dcompare (const void * a, const void * b)
{
  if (*(double*)a < *(double*)b) return -1;
  if (*(double*)a > *(double*)b) return 1;
  
  return 0;
}

inline Edge *Node::shortestEdge(bool mst) {
  Edge *emin = NULL;
  double wmin = DBL_MAX;
  if (neighbors.size() == 0) return NULL;

  list<Edge *>::iterator it = neighbors.begin();
  for (int k = 0; k < neighbors.size() ; ++k, ++it) {
    Edge *e = *it;
    if (mst || !e->otherNode(this)->onMST) {
      double wt = e->getWeight();
      if (wt < wmin) {
        wmin = wt;
        emin = e;
      }
    }
  }

//cout << "node: " << *emin << " min: " << wmin << endl;

  return emin;
}


inline int Node::getMSTChildren(double cutoff, double *maxLength) {
  int count = 0;
  if (isVisited)
    return count;
  isVisited = true;

  list<Edge *>::iterator it = neighbors.begin();
  for (int k =0; k < neighbors.size(); ++k, ++it) {
    Edge *e = *it;
    if (e->onMST) {
      if (e->getWeight() < cutoff) {
        if (!e->otherNode(this)->isVisited) {
          count += e->otherNode(this)->getMSTChildren(cutoff, maxLength);
          double el = e->getWeight();
          if (el > maxLength[0])
            maxLength[0] = el;
        }
      }
    }
  }
  count += this->count; // add count for this node
  return count;
}


inline void Node::setNeighbor(Edge *neighbor) {
  neighbors.push_back(neighbor);

/*cout << "node: " << pointID << "(" << this << ") nb: " << *neighbor << endl;
list<Edge*>::iterator it = neighbors.begin();
for (int k=0; k < neighbors.size();++k, ++it) {
  cout << **it << " ";
}
cout << endl;
*/
}


inline bool Edge::isNewEdge(Node *n) {
  list<Edge *>::iterator it = n->neighbors.begin();
  
  for (int k = 0; k < n->neighbors.size() ; ++k, ++it) {
    if (isEquivalent(**it))
      return false;
  }
  return true;
}

inline int Edge::getRunts(double *maxLength) {
  double cutoff = weight;
  double maxLength1;
  double maxLength2;
  int count1 = p1->getMSTChildren(cutoff, &maxLength1);
  int count2 = p2->getMSTChildren(cutoff, &maxLength2);
  if (count1 < count2) {
    *maxLength = maxLength1;
    return count1;
  } else if (count1 == count2) {        // take more tightly clustered child
    if (maxLength1 < maxLength2)
      *maxLength = maxLength1;
    else
      *maxLength = maxLength2;
    return count1;
  } else {
    *maxLength = maxLength2;
    return count2;
  }
}



const char *Triangulation::scagnosticsLabels[9] = {"Outlying", "Skewed", "Clumpy", "Sparse", "Striated", "Convex", "Skinny", "Stringy", "Monotonic"};

double *Triangulation::compute(BinnedData bdata, bool quick) {
  this->bdata = &bdata;
  this->px = bdata.getXData();
  this->py = bdata.getYData();
  this->np = bdata.getNBins();
  if (this->np < 2)
    return NULL;

  findOutliers(bdata);
  
  if (ALPHA)
    computeAlphaGraph();
/*  if (GABRIEL)
    computeGabrielGraph();
  if (RNG)
    computeRelativeNeighborhoodGraph();
  if (NNG)
    computeNearestNeighborGraph();
  if (DIST)
    computeDistanceGraph(100.);
  if (SPHERE)
    computeSphereOfInfluenceGraph();*/

  if (quick)
    return NULL; 

  computeTotalCount();
  computeAlphaArea();
  computeAlphaPerimeter();
  computeHullArea();
  computeHullPerimeter();
  
  double *measures = computeMeasures();

  // for (int k=0; k < numScagnostics; ++k) 
  //   printf("%3.5f ", measures[k]);
  // printf("\n");

  return measures;
}

inline double *Triangulation::computeMeasures() {
  double *results = new double[numScagnostics];
  // Do not change order of these calls!

  results[OUTLYING] = computeOutlierMeasure();
  results[CLUMPY] = computeClusterMeasure();
  results[SKEWED] = computeMSTEdgeLengthSkewnessMeasure();
  results[CONVEX] = computeConvexityMeasure();
  results[SKINNY] = computeSkinnyMeasure();
  results[STRINGY] = computeStringyMeasure();
  results[STRIATED] = computeStriationMeasure();
  results[SPARSE] = computeSparsenessMeasure();
  results[MONOTONIC] = computeMonotonicityMeasure();

  return results;
}


inline double Triangulation::computeMSTEdgeLengthSkewnessMeasure() {
  if (mstEdges.size() == 0)
    return 0;
  int n = mstEdges.size();
  int n50 = n / 2;
  int n10 = n / 10;
  int n90 = (9 * n) / 10;
  double skewness = (sortedOriginalMSTLengths[n90] - sortedOriginalMSTLengths[n50]) /
      (sortedOriginalMSTLengths[n90] - sortedOriginalMSTLengths[n10]);
  double t = (double) totalCount / 500;
  double correction = .7 + .3 / (1 + t * t);
  return 1 - correction * (1 - skewness);
}

inline double Triangulation::computeStringyMeasure() {
  int count1 = 0;
  int count2 = 0;

  list<Node *>::iterator it = nodes.begin();
  for (int k=0; k < nodes.size(); k++,it++) {
    Node *n = *it;
    if (n->getMSTDegree() == 1)
      count1++;
    if (n->getMSTDegree() == 2)
      count2++;
  }
  double result = (double) count2 / (double) (nodes.size() - count1);
  return result * result * result;
}

inline double Triangulation::computeClusterMeasure() {
  double maxLength;
  double maxValue = 0;
  list<Edge *>::iterator it = mstEdges.begin();
  for (int k=0; k < mstEdges.size(); ++k,++it) {
    Edge *e = *it;
    clearVisits();
    e->onMST = false;  // break MST at this edge
    int runts = e->getRunts(&maxLength);
    e->onMST = true;   // restore this edge to MST
    if (maxLength > 0) {
      double value = runts * (1 - maxLength / e->getWeight());
      if (value > maxValue)
        maxValue = value;
    }
  }
  return 2 * maxValue / totalPeeledCount;
}

inline void Triangulation::clearVisits() {
  list<Node *>::iterator it = nodes.begin();
  for (int k=0; k < nodes.size(); ++k,++it) {
    Node *n = *it;
    n->isVisited = false;
  }
}

inline double Triangulation::computeMonotonicityMeasure() {
  int n = np;
  double *ax = new double[n];
  double *ay = new double[n];
  double *weights = new double[n];
  for (int i = 0; i < n; ++i) {
    ax[i] = px[i];
    ay[i] = py[i];
    weights[i] = counts[i];
  }
//  double *rx = Sort.rank(ax, weights);
//  double *ry = Sort.rank(ay, weights);
  double s = computePearson(n, ax, ay, weights);
  return s * s;
}

inline double Triangulation::computePearson(int n, double *x, double *y, double *weights) {
  double xmean = 0;
  double ymean = 0;
  double xx = 0;
  double yy = 0;
  double xy = 0;
  double wt = 0;
  double sumwt = 0;
  for (int i = 0; i < n; ++i) {
    wt = weights[i];
    if (wt > 0 && !isOutlier[i]) {
      sumwt += wt;
      xx += (x[i] - xmean) * wt * (x[i] - xmean);
      yy += (y[i] - ymean) * wt * (y[i] - ymean);
      xy += (x[i] - xmean) * wt * (y[i] - ymean);
      xmean += (x[i] - xmean) * wt / sumwt;
      ymean += (y[i] - ymean) * wt / sumwt;
    }
  }
  xy /= sqrt(xx * yy);
  return xy;
}

inline double Triangulation::computeSparsenessMeasure() {
  int n = mstEdges.size();
  int n90 = (9 * n) / 10;
  double sparse = min(sortedOriginalMSTLengths[n90] / 1000, 1.);
  double t = (double) totalCount / 500;
  double correction = .7 + .3 / (1 + t * t);
  return correction * sparse;
}

inline double Triangulation::computeStriationMeasure() {
  double numEdges = 0;

  list<Edge *>::iterator it = mstEdges.begin();
  for (int k=0; k < mstEdges.size(); ++k,++it) {
    Edge *e = *it;
    Node *n1 = e->p1;
    Node *n2 = e->p2;
    if (n1->getMSTDegree() == 2 && n2->getMSTDegree() == 2) {
      Edge *e1 = getAdjacentMSTEdge(n1, e);
      Edge *e2 = getAdjacentMSTEdge(n2, e);
      if (cosineOfAdjacentEdges(e, e1, n1) < -.7 && cosineOfAdjacentEdges(e, e2, n2) < -.7)
        ++numEdges;
    }
  }
  return numEdges / (double) mstEdges.size();
}

inline Edge *Triangulation::getAdjacentMSTEdge(Node *n, Edge *e) {
  list<Edge *>::iterator nt = n->neighbors.begin();

  for (int k=0; k < n->neighbors.size(); ++k, ++nt) {
    Edge *et = *nt;
    if (et->onMST && (e != et)) {
      return et;
    }
  }
  return NULL;
}

inline double Triangulation::cosineOfAdjacentEdges(Edge *e1, Edge *e2, Node *n) {
  double v1x = e1->otherNode(n)->x - n->x;
  double v1y = e1->otherNode(n)->y - n->y;
  double v2x = e2->otherNode(n)->x - n->x;
  double v2y = e2->otherNode(n)->y - n->y;
  double v1 = sqrt(v1x * v1x + v1y * v1y);
  double v2 = sqrt(v2x * v2x + v2y * v2y);
  v1x /=  v1;
  v1y /=  v1;
  v2x /=  v2;
  v2y /=  v2;
  return v1x * v2x + v1y * v2y;
}

inline double Triangulation::computeConvexityMeasure() {
  if (hullArea == 0) // points in general position
    return 1;
  else {
    double t = (double) totalCount / 500;
    double correction = .7 + .3 / (1 + t * t);
    double convexity = alphaArea / hullArea;
    return correction * convexity;
  }
}

inline double Triangulation::computeSkinnyMeasure() {
  if (alphaPerimeter > 0)
    return 1 - sqrt(4 * M_PI * alphaArea) / alphaPerimeter;
  else
    return 1;
}

inline void Triangulation::computeAlphaArea() {
  double area = 0;

  list<Triangle *>::iterator tri = triangles.begin();
  for (int k = 0; k < triangles.size() ; ++k, ++tri) {
    Triangle *t=*tri;
    if (t->onComplex) {
      Node *p1 = t->anEdge->p1;
      Node *p2 = t->anEdge->p2;
      Node *p3 = t->anEdge->nextE->p2;
      area += abs(p1->x * p2->y + p1->y * p3->x + p2->x * p3->y
          - p3->x * p2->y - p3->y * p1->x - p1->y * p2->x);
    }
  }
  alphaArea = area / 2;
}

inline void Triangulation::computeHullArea() {
  double area = 0;

  list<Triangle *>::iterator tri = triangles.begin();
  for (int k = 0; k < triangles.size() ; ++k, ++tri) {
    Triangle *t=*tri;
    Node *p1 = t->anEdge->p1;
    Node *p2 = t->anEdge->p2;
    Node *p3 = t->anEdge->nextE->p2;
    area += abs(p1->x * p2->y + p1->y * p3->x + p2->x * p3->y
        - p3->x * p2->y - p3->y * p1->x - p1->y * p2->x);
  }
  hullArea = area / 2.;
}

inline void Triangulation::computeAlphaPerimeter() {
  double sum = 0;

  list<Edge *>::iterator it=edges.begin();
  for (int k = 0; k < edges.size() ; ++k, ++it) {
    Edge *e = *it;
    if (e->onShape) {
      sum += e->getWeight();
    }
  }
  alphaPerimeter = sum;
}

inline void Triangulation::computeHullPerimeter() {
  double sum = 0;
  Edge *e = hullStart;
  do {
    sum += e->p1->distToNode(e->p2->x, e->p2->y);
    e = e->nextH;
  } while (e != hullStart);
  hullPerimeter = sum;
}

inline void Triangulation::computeAlphaGraph() { // requires initializing Edge.onShape = false
  bool deleted;
  double alpha = computeAlphaValue();

  do {
    list<Edge *>::iterator i = edges.begin();
    deleted=false;
    for (int k=0; k < edges.size(); ++i, ++k) {
      Edge *e=*i;
      if (e->inT->onComplex) {
        if (alpha < e->getWeight() / 2) {
          e->inT->onComplex = false;
          deleted = true;
        } else {
          if (e->invE != NULL)
            if (e->invE->inT->onComplex)
              continue;
          if (!edgeIsExposed(alpha, e)) {
            e->inT->onComplex = false;
            deleted = true;
          }
        }
      }

    }
  } while (deleted);  
  markShape();
}

inline bool Triangulation::edgeIsExposed(double alpha, Edge *e) {
  double x1 = e->p1->x;
  double x2 = e->p2->x;
  double y1 = e->p1->y;
  double y2 = e->p2->y;
  double xe = (x1 + x2) / 2;
  double ye = (y1 + y2) / 2;
  double d = sqrt(alpha * alpha - e->getWeight() * e->getWeight() / 4);
  double xt = d * (y2 - y1) / e->getWeight();
  double yt = d * (x2 - x1) / e->getWeight();
  double xc1 = xe + xt;
  double yc1 = ye - yt;
  double xc2 = xe - xt;
  double yc2 = ye + yt;
  bool pointsInCircle1 = pointsInCircle(e->p1, xc1, yc1, alpha) ||
      pointsInCircle(e->p2, xc1, yc1, alpha);
  bool pointsInCircle2 = pointsInCircle(e->p1, xc2, yc2, alpha) ||
      pointsInCircle(e->p2, xc2, yc2, alpha);
  if (pointsInCircle1 && pointsInCircle2)
    return false;
  return true;
}

inline bool Triangulation::pointsInCircle(Node *n, double xc, double yc, double radius) {
  double r = .999 * radius;

  list<Edge *>::iterator i = n->neighbors.begin();
  for (int k=0; k < n->neighbors.size(); ++i, ++k) {
    Edge *e=*i;
    Node *no = e->otherNode(n);
    double dist = no->distToNode(xc, yc);
    if (dist < r)
      return true;
  }
  return false;
}

inline void Triangulation::markShape() {
  list<Edge *>::iterator i = edges.begin();
  for (int k=0; k < edges.size(); ++i, ++k) {
    Edge *e = *i;
    e->onShape = false;
    if (e->inT->onComplex) {
      if (e->invE == NULL) {
        e->onShape = true;
      } else if (!e->invE->inT->onComplex)
        e->onShape = true;    
    }
  }
}

inline double Triangulation::computeAlphaValue() {
  int length = mstEdges.size();
  int n90 = (9 * length) / 10;
  double alpha = sortedOriginalMSTLengths[n90];
  return min(alpha, 100.);
}

inline double *Triangulation::getSortedMSTEdgeLengths() {
  
  double *result = new double[mstEdges.size()];

  list<Edge *>::iterator it = mstEdges.begin();
  for (int k = 0; k < mstEdges.size(); ++k, ++it) {
    result[k] = (*it)->getWeight();
  }
  
  qsort(result, mstEdges.size(), sizeof(double), dcompare);

  return result;
}

inline double Triangulation::computeCutoff(int n, double* lengths) {
  if (n == 0) return 0;
  int n50 = n / 2;
  int n25 = n50 / 2;
  int n75 = n50 + n25;
  return lengths[n75] + 1.5 * (lengths[n75] - lengths[n25]);
}

inline void Triangulation::computeTotalOriginalMSTLengths(int n) {
  for (int i = 0; i < n; ++i)
    totalOriginalMSTLengths += sortedOriginalMSTLengths[i];
}

inline bool Triangulation::computeMSTOutliers(double omega) {
  bool found = false;

  list<Node *>::iterator it = nodes.begin();
  for (int k = 0; k < nodes.size(); ++k, ++it) {
    Node *n = *it;
//cout << endl;
//cout << "node " << *n << " (" << n << ")";
    bool deleteNode = true;
    double sumlength = 0;
    list<Edge *>::iterator ie = n->neighbors.begin();
    for (int j=0; j < n->neighbors.size(); ++j, ++ie) {
      Edge *e = *ie;
//cout << " " << *e;

      if (e->onMST) {
//cout << "mst ";
        double length = e->getWeight();
        sumlength += length;
        if (length < omega)
          deleteNode = false;
      }
    }
    if (deleteNode) {
      totalMSTOutlierLengths += sumlength;
      isOutlier[n->pointID] = true;
      found = true;
    }

  }
  return found;
}

inline void Triangulation::findOutliers(BinnedData bdata) {
  this->counts = bdata.getCounts();
  isOutlier = new bool[np];
  for(int i = 0; i < np; ++i) isOutlier[i]=false;

  computeDT(px, py);

/*cout << "triangles: size " << triangles.size() <<endl;  
  list<Triangle *>::iterator it = triangles.begin();
  for (int k = 0; k < triangles.size() ; ++k, ++it) {
    cout << **it << " ";
  }
cout << endl;

cout << "edges: size " << edges.size()<<endl;  
  list<Edge *>::iterator ig = edges.begin();
  for (int k = 0; k < edges.size() ; ++k, ++ig) {
    cout << **ig << " ";
  }
cout << endl;
*/
  computeMST();

/*cout << "MST edges: size " << mstEdges.size() <<endl;  
  list<Edge *>::iterator ie = mstEdges.begin();
  for (int k = 0; k < mstEdges.size() ; ++k, ++ie) {
    cout << **ie << " ";
  }
cout << endl;
*/
  sortedOriginalMSTLengths = getSortedMSTEdgeLengths();

  double cutoff = computeCutoff(mstEdges.size(), sortedOriginalMSTLengths);

  computeTotalOriginalMSTLengths(mstEdges.size());
//cout << " totalOriginalMSTLengths: " << totalOriginalMSTLengths << endl;

  bool foundNewOutliers = computeMSTOutliers(cutoff);

  double* sortedPeeledMSTLengths;
  while (foundNewOutliers) {
    clear();
    computeDT(px, py);
    computeMST();
    sortedPeeledMSTLengths = getSortedMSTEdgeLengths();
    cutoff = computeCutoff(mstEdges.size(),sortedPeeledMSTLengths);
    foundNewOutliers = computeMSTOutliers(cutoff);
  } 

}

inline void Triangulation::updateMSTEdges(Edge *addEdge, list<Edge *> *mstEdges) {
  mstEdges->push_back(addEdge);
  addEdge->onMST = true;
  addEdge->p1->incDegree(1);
  addEdge->p2->incDegree(1);
}

inline void Triangulation::updateMSTNodes(Node *addNode, list<Node *> *mstNodes) {
  mstNodes->push_back(addNode);
  addNode->onMST = true;
}

inline void Triangulation::computeMST() {
// Prim's algorithm
  if (nodes.size() > 1) {
    list<Node *> mstNodes;
    Node *mstNode = nodes.front();  // pick arbitrary node
    updateMSTNodes(mstNode, &mstNodes);
    
    int count = 1;
    while (count < nodes.size()) {
      Edge *addEdge = NULL;
      double wmin = DBL_MAX;
      Node *nmin = NULL;
  
      list<Node *>::iterator  it = mstNodes.begin();
      int k = 0;
      while (k < mstNodes.size()) {
        mstNode = *it;
        Edge *candidateEdge = mstNode->shortestEdge(false);
        if (candidateEdge != NULL) {
          double wt = candidateEdge->getWeight();
          if (wt < wmin) {
            wmin = wt;
            nmin = mstNode;
            addEdge = candidateEdge;
          }
        }
        ++k;
        ++it;
      }

      if (addEdge != NULL) {
        Node *addNode = addEdge->otherNode(nmin);
        updateMSTNodes(addNode, &mstNodes);
        updateMSTEdges(addEdge, &mstEdges);

        ++count;
      }
    }
  }
}

inline void Triangulation::computeDT(int *px, int *py) {
  totalPeeledCount = 0;
  srand(13579);  // set seed
  for (int i = 0; i < np; ++i) {
    int x = px[i] + (int) (8 * ((double)rand()/RAND_MAX - .5)); // perturb to prevent singularities
    int y = py[i] + (int) (8 * ((double)rand()/RAND_MAX - .5));
    int count = counts[i];
    if (!isOutlier[i]) {
      insert(x, y, count, i);
      totalPeeledCount += count;
    }
  }
  setNeighbors();
  markHull();
}

inline void Triangulation::setNeighbors() {  
  list<Edge *>::iterator it = edges.begin();
  
  for (int k = 0; k < edges.size() ; ++k, ++it) {
    Edge *e = *it;
    if (e->isNewEdge(e->p2)) e->p2->setNeighbor(e);
    if (e->isNewEdge(e->p1)) e->p1->setNeighbor(e);
  }
}

inline void Triangulation::markHull() {
  Edge *e = hullStart;
  if (e != NULL)
    do {
      e->onHull = true;
      e->p1->onHull = true;
      e->p2->onHull = true;
      e = e->nextH;
    } while (*e != *hullStart);
}

inline void Triangulation::insert(int px, int py, int count, int id) {
  int eid;
  Node *nd = new Node(px, py, count, id);
  nodes.push_back(nd);

//cout << *nd << "(" << nd << ")" << endl;

  list<Node *>::iterator p = nodes.begin();
  if (nodes.size() < 3) return;
  if (nodes.size() == 3) {   // create the first triangle
    Node *p1 = *p;  // p points to element 0
    ++p; Node *p2 = *p;  // p points to element 1 
    ++p; Node *p3 = *p;  // p points to element 2
    Edge *e1 = new Edge(p1, p2);
    if (e1->onSide(p3) == 0) {
      nodes.remove(nd);
      return;
    }
    if (e1->onSide(p3) == -1) { // right side
      p = nodes.begin();
      p2 = *p; ++p;  // swap first two elements
      p1 = *p;
      e1->update(p1, p2);
    }
    Edge *e2 = new Edge(p2, p3);
    Edge *e3 = new Edge(p3, p1);
    e1->nextH = e2;
    e2->nextH = e3;
    e3->nextH = e1;
    hullStart = e1;
    triangles.push_back(new Triangle(&edges, e1, e2, e3));

/*cout << "edges: size " << edges.size()<<endl;  
copy(edges.begin(), edges.end(),
          ostream_iterator<Edge,char>(cout," "));
cout << endl;
cout << "insert return 1" << endl;
*/    return;
  }
  actE = edges.front();    
  if (actE->onSide(nd) == -1) {
    if (actE->invE == NULL)
      eid = -1;
    else
      eid = searchEdge(actE->invE, nd);
  } else
    eid = searchEdge(actE, nd);

//cout << "node: " << *nd << " search: " << eid << endl;

  if (eid == 0) {
    nodes.remove(nd);
    return;
  }
  if (eid > 0)
    expandTri(actE, nd, eid);   // nd is inside or on a triangle
  else
    expandHull(nd);                // nd is outside convex hull

}

inline void Triangulation::swapTest(Edge *e11) {
/*
Idea: SwapTest(ab):
  if (ab is an edge on the exterior face) return
  Let d be the vertex to the right of edge ab;
    if (inCirc(p, a, b, d) )                // d violates the incircle test
       Flip edge ab for pd;
       SwapTest(ad); // Fix the new suspect edges
       SwapTest(db);
*/
  Edge *e21 = e11->invE;
  if (e21 == NULL || e21->inT == NULL) return;

  Edge *e12 = e11->nextE;
  Edge *e13 = e12->nextE;
  Edge *e22 = e21->nextE;
  Edge *e23 = e22->nextE;

  if (e11->inT->inCircle(e22->p2) || e21->inT->inCircle(e12->p2)) {

    e11->update(e22->p2, e12->p2);
    e21->update(e12->p2, e22->p2);
    e11->linkSymm(e21);
//cout << "e1 " << e13->p1->y << e22->p1->y << e11->p1->y << endl;
//cout << "e2 " << e23->p1->y << e12->p1->y << e21->p1->y << endl;
    e13->inT->update(e13, e22, e11);
    e23->inT->update(e23, e12, e21);
//cout << "t1 " << *e13->inT << endl;
//cout << "t2 " << *e23->inT << endl;
    e12->asIndex();
    e22->asIndex();
//cout << "t1 " << *e12->inT << endl;
//cout << "t2 " << *e22->inT << endl;
    swapTest(e12);
    swapTest(e22);
    swapTest(e13);
    swapTest(e23);
  }
}

inline void Triangulation::expandTri(Edge *e, Node *nd, int type) {
//cout << "expandTri" << endl;
  Edge *e1 = e;
  Edge *e2 = e1->nextE;
  Edge *e3 = e2->nextE;
  Node *p1=e1->p1;
  Node *p2=e2->p1;
  Node *p3=e3->p1;
  if (type == 2) {   // nd is inside of the triangle
    Edge *e10 = new Edge(p1, nd);
    Edge *e20 = new Edge(p2, nd);
    Edge *e30 = new Edge(p3, nd);
//    cout << "e: " << *e << endl;
    triangles.remove(e->inT);     // remove old triangle
    e->inT->removeEdges(&edges);

    Edge *e100 = e10->makeSymm();
    Edge *e200 = e20->makeSymm();
    Edge *e300 = e30->makeSymm();
    triangles.push_back(new Triangle(&edges, e1, e20, e100));
    triangles.push_back(new Triangle(&edges, e2, e30, e200));
    triangles.push_back(new Triangle(&edges, e3, e10, e300));
    swapTest(e1);   // swap test for the three new triangles
    swapTest(e2);
    swapTest(e3);
  } else  {         // nd is on the edge e
    Edge *e4 = e1->invE;
    if (e4 == NULL || e4->inT == NULL) {          // one triangle involved
      Edge *e30 = new Edge(p3, nd);
      Edge *e02 = new Edge(nd, p2);
      Edge *e10 = new Edge(p1, nd);
      Edge *e03 = e30->makeSymm();
//                shareEdges(e03,e30);
      e10->asIndex();
      e1->mostLeft()->nextH = e10;
      e10->nextH = e02;
      e02->nextH = e1->nextH;
      hullStart = e02;
      triangles.remove(e1->inT);  // remove oldtriangle and add two new triangles
      e->inT->removeEdges(&edges);
      edges.remove(e1);
      edges.push_back(e10);
      edges.push_back(e02);
      edges.push_back(e30);
      edges.push_back(e03);
      triangles.push_back(new Triangle(e2, e30, e02));
      triangles.push_back(new Triangle(e3, e10, e03));
      swapTest(e2);   // swap test for the two new triangles
      swapTest(e3);
      swapTest(e30);
    } else  {      // two triangle involved
      Edge *e5 = e4->nextE;
      Edge *e6 = e5->nextE;
      Node *p4 = e6->p1;
      Edge *e10 = new Edge(p1, nd);
      Edge *e20 = new Edge(p2, nd);
      Edge *e30 = new Edge(p3, nd);
      Edge *e40 = new Edge(p4, nd);
      triangles.remove(e->inT);                   // remove oldtriangle
      e->inT->removeEdges(&edges);
      triangles.remove(e4->inT);               // remove old triangle
      e4->inT->removeEdges(&edges);
      e5->asIndex();   // because e, e4 removed, reset edge sortOrder of node p1 and p2
      e2->asIndex();
      triangles.push_back(new Triangle(&edges, e2, e30, e20->makeSymm()));
      triangles.push_back(new Triangle(&edges, e3, e10, e30->makeSymm()));
      triangles.push_back(new Triangle(&edges, e5, e40, e10->makeSymm()));
      triangles.push_back(new Triangle(&edges, e6, e20, e40->makeSymm()));
      swapTest(e2);   // swap test for the three new triangles
      swapTest(e3);
      swapTest(e5);
      swapTest(e6);
      swapTest(e10);
      swapTest(e20);
      swapTest(e30);
      swapTest(e40);
    }
  }
}

inline int Triangulation::searchEdge(Edge *e, Node *nd) {
  int f2,f3;
  Edge *e0 = NULL;
  f2 = e->nextE->onSide(nd);
  if (f2 == -1) {
    if (e->nextE->invE != NULL)
      return searchEdge(e->nextE->invE, nd);
    else {
      actE = e;
      return -1;
    }
  }
  if (f2 == 0) e0 = e->nextE;

  Edge *ee = e->nextE;
  f3 = ee->nextE->onSide(nd);
  if (f3 == -1) {
    if (ee->nextE->invE != NULL)
      return searchEdge(ee->nextE->invE, nd);
    else {
      actE = ee->nextE;
      return -1;
    }
  }
  if (f3 == 0) e0 = ee->nextE;
  if (e->onSide(nd) == 0) e0 = e;
  if (e0 != NULL) {
    actE = e0;
    if (e0->nextE->onSide(nd) == 0) {
      actE = e0->nextE;
      return 0;
    }
    if (e0->nextE->nextE->onSide(nd) == 0) return 0;
    return 1;
  }
  actE = ee;
  return 2;
}



inline void Triangulation::expandHull(Node *nd) {
//cout << "expandHull" << endl;
/*cout << "triangles: size " << triangles.size() <<endl;  
  list<Triangle *>::iterator it = triangles.begin();
  for (int k = 0; k < triangles.size() ; ++k, ++it) {
    cout << **it << " ";
  }
cout << endl;
*/
  Edge *e1=NULL,*e2=NULL,*e3=NULL,*enext=NULL;
  Edge *e = hullStart;
  Edge *comedge = NULL,*lastbe = NULL;

  while (true) {
    enext = e->nextH;
    if (e->onSide(nd) == -1) {  // right side
      if (lastbe != NULL) {
        e1 = e->makeSymm();
        e2 = new Edge(e->p1, nd);
        e3 = new Edge(nd, e->p2);
        if (comedge == NULL) {
          hullStart = lastbe;
          lastbe->nextH = e2;
          lastbe = e2;
        } else{
          comedge->linkSymm(e2);
        }
        comedge = e3;
        Triangle *trn = new Triangle(&edges, e1, e2, e3);
        triangles.push_back(trn);
        swapTest(e);
      }
    } else {
      if (comedge != NULL) break;
      lastbe = e;
    }
    e = enext;
  }

  lastbe->nextH = e3;
  e3->nextH = e;
  
  
/*cout << "hull is" << endl;
Edge *p = hullStart;
while (p->nextH != hullStart) {
  cout << *p << " ";
  p = p->nextH;
}
cout << *p << " " << *p->nextH << endl;
*/
}

inline void Triangle::update(Edge *e1, Edge *e2, Edge *e3) {
  onComplex = true;
  anEdge = e1;
  e1->nextE = e2;
  e2->nextE = e3;
  e3->nextE = e1;
  e1->inT = this;
  e2->inT = this;
  e3->inT = this;
//  cout <<"setup edges now finding circle" << endl;
  findCircle();
}


inline void Triangle::removeEdges(list<Edge *> *edges) {
//cout << "rm " << *anEdge;
//cout << " " << *anEdge->nextE;
//cout << " " << *anEdge->nextE->nextE << endl;

  edges->remove(anEdge);
  edges->remove(anEdge->nextE);
  edges->remove(anEdge->nextE->nextE);
}

inline void Triangle::findCircle() {
  double x1 = (double) anEdge->p1->x;
  double y1 = (double) anEdge->p1->y;
  double x2 = (double) anEdge->p2->x;
  double y2 = (double) anEdge->p2->y;
  double x3 = (double) anEdge->nextE->p2->x;
  double y3 = (double) anEdge->nextE->p2->y;
  double a = (y2 - y3) * (x2 - x1) - (y2 - y1) * (x2 - x3);
  double a1 = (x1 + x2) * (x2 - x1) + (y2 - y1) * (y1 + y2);
  double a2 = (x2 + x3) * (x2 - x3) + (y2 - y3) * (y2 + y3);
  c_cx = (a1 * (y2 - y3) - a2 * (y2 - y1)) / a / 2.;
  c_cy = (a2 * (x2 - x1) - a1 * (x2 - x3)) / a / 2.;
  c_r = anEdge->p1->distToNode(c_cx, c_cy);
}


