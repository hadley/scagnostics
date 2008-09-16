/*
 * Graph-theoretic Scagnostic Measures
 *
 * Leland Wilkinson (SPSS, Inc.) and Anushka Anand (University of Illinois at Chicago)
 * This program accompanies the paper by Leland Wilkinson, Anushka Anand, and Robert Grossman
 * called Graph-Theoretic Scagnostics
 * Proceedings of the IEEE Symposium on Information Visualization
 * Minneapolis, MN October 23-25, 2005.
 *
 * Delaunay triangulation adapted from Delaunay.java by Marcus Apel (www.geo.tu-freiberg.de/~apelm/)
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software. Supporting documentation must also include a citation of
 * the abovementioned article, Graph-Theoretic Scagnostics
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, THE AUTHORS MAKE NO
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

#ifndef GRAPHMEASURES_H
#define GRAPHMEASURES_H
 
#include <list>
#include <string>
#include <vector>
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <iterator>

#include "Binner.h"
using namespace std;
using std::list;

class Edge;
class Triangle;

class Node {
    int count;     // number of points aggregated at this node
    Edge *anEdge;     // an edge which starts from this node
    char type;      // 1: inner nodes, 2: on convex hull, 3: ...
    int mstDegree;
    int nodeID;

public:
    bool isVisited; // = false;
    bool onMST;
    bool onHull; // = false;
    list<Edge *> neighbors; // nearest Delaunay neighbors list
    int pointID;
    int x,y;       // coordinate X,Y

    Node(int x, int y, int count, int pointID) {
    onHull = false;
    onMST = false;
    isVisited = false;
    mstDegree = 0;

        this->x = x;
        this->y = y;
        this->count = count;
        anEdge = NULL;
        this->pointID = pointID;
    }

  bool operator == (Node p) {
    if (x != p.x) return false;
    if (y != p.y) return false;
    return true;
  }

  bool operator != (Node p) { return(!(*this == p)); }

  void incDegree(int by) { mstDegree += by; }
  int getMSTDegree() { return mstDegree; }
  
    double distToNode(double px, double py) {
        double dx = px - x;
        double dy = py - y;
        return sqrt(dx * dx + dy * dy);
    }

    void setNeighbor(Edge *neighbor);
  
  Edge *shortestEdge(bool mst);

  int getMSTChildren(double cutoff, double *maxLength);

  friend class Edge;
  friend class Triangle;
};

class Edge {
    double a,b,c;               // line equation parameters. aX+bY+c=0
    double weight;

  public:
    bool onHull;// = false;
    bool onMST;// = false;
    bool onShape;// = false;

    Node *p1,*p2;         // start and end point of the edge
    Triangle *inT; // = NULL;      // triangle containing this edge
    Edge *nextH; // = NULL;    // convex hull link
    Edge *invE; // = NULL;     // inverse edge (p2->p1)
    Edge *nextE; // = NULL;    // next edge in the triangle in counterclockwise

    Edge(Node *p1, Node *p2) {
    invE = NULL;     // inverse edge (p2->p1)
    nextE = NULL;    // next edge in the triangle in counterclockwise
    nextH = NULL;    // convex hull link
    inT = NULL;      // triangle containing this edge

    onHull = false;
    onMST = false;
    onShape = false;
    
        update(p1, p2);
    }

  double getWeight() { return weight; }
//**************** --anushka -- *******************************
//commented the following and translated directly from java - makes code 1-2 secs faster on larger data sets (bank, boston)
/*    bool isEquivalent(Edge e) {
    if (e == *this) return true;  // same edges
    // or inverted edge
    return ((*e.p1 == *p2) && (*e.p2 == *p1));
  }

  bool operator == (Edge e) {
    if (*p1 != *e.p1) return false;
    if (*p2 != *e.p2) return false;

    return true;
  }
*/
//translated from java
  bool operator==(Edge e) {
        if ((e.p1->x == this->p1->x) && (e.p2->x == this->p2->x) && (e.p1->y == this->p1->y) && (e.p2->y == this->p2->y))
            return true;
        else
            return false;
    }

    bool isEquivalent(Edge e) {
        if (((e.p1->x == this->p1->x) && (e.p2->x == this->p2->x) && (e.p1->y == this->p1->y) && (e.p2->y == this->p2->y)) ||
                ((e.p1->x == this->p2->x) && (e.p1->y == this->p2->y) && (e.p2->x == this->p1->x) && (e.p2->y == this->p1->y)))
            return true;
        else
            return false;
    }
/*********************************************************/
    
  bool operator != (Edge e) {
    return (!(*this == e));
  }

  void update(Node *p1, Node *p2) {
        this->p1 = p1;
        this->p2 = p2;
        a = p2->y - p1->y;
        b = p1->x - p2->x;
        c = p2->x * p1->y - p1->x * p2->y;
        weight = sqrt(a * a + b * b);
        asIndex();
    }

    void asIndex() {
        p1->anEdge = this;
    }

    int onSide(Node *nd) {
        double s = a * nd->x + b * nd->y + c;
        if (s > 0.0) return 1;
        if (s < 0.0) return -1;
        return 0;
    }

    Edge *makeSymm() {
        Edge *e = new Edge(p2, p1);
        linkSymm(e);
        return e;
    }

    void linkSymm(Edge *e) {
        this->invE = e;
        if (e != NULL) e->invE = this;
    }

    Edge *mostLeft() {
        Edge *ee, *e = this;
        while ((ee = e->nextE->nextE->invE) != NULL && ee != this) e = ee;
        return e->nextE->nextE;
    }

    Edge *mostRight() {
        Edge *ee,*e = this;
        while (e->invE != NULL && (ee = e->invE->nextE) != this) e = ee;
        return e;
    }

    bool isNewEdge(Node *n);

    Node *otherNode(Node *n) {
        if (*n == *p1)
            return p2;
        else
            return p1;
    }


  int getRunts(double *maxLength);

/*

    void deleteSimplex() {
        onShape = false;
        inT.onComplex = false;
        if (invE != NULL) {
            invE.onShape = false;
            invE.inT.onComplex = false;
        }
    }
*/


  friend class Triangle;
};

class Triangle {
    double c_cx;        // center of circle: X
    double c_cy;        // center of circle: Y
    double c_r;         // radius of circle


public:
    bool onComplex;
    Edge *anEdge;        // an edge of this triangle
    Triangle(Edge *e1, Edge *e2, Edge *e3) {
    onComplex  = true;

        update(e1, e2, e3);
    }

    Triangle(list<Edge *> *edges, Edge *e1, Edge *e2, Edge *e3) {
    onComplex  = true;
    list<Edge *> es = *edges;
        update(e1, e2, e3);
        edges->push_back(e1);
        edges->push_back(e2);
        edges->push_back(e3);
    }

  bool operator == (Triangle t) {
    Edge *e1 = anEdge;
    Edge *e2 = e1->nextE;
    Edge *e3 = e2->nextE;

    Edge *te = t.anEdge->nextE;
    Edge *tee = te->nextE;
    if (*e1 == *t.anEdge) {
      return ((*e2 == *te) && (*e3 == *tee));
    }
    if (*e2 == *t.anEdge) {
      return ((*e3 == *te) && (*e1 == *tee));
    }
    if (*e3 == *t.anEdge) {
      return ((*e1 == *te) && (*e2 == *tee));
    }

    return false;
  }

    void update(Edge *e1, Edge *e2, Edge *e3);

    bool inCircle(Node *nd) {
        return nd->distToNode(c_cx, c_cy) < c_r;
    }

    void removeEdges(list<Edge *> *edges);
    void findCircle();
};


class Triangulation {
    static const bool ALPHA = true;      // set Edge.onShape=false before using
    static const bool GABRIEL = false;   // set Edge.onShape=true before using
    static const bool RNG = false;       // set Edge.onShape=true before using
    static const bool NNG = false;       // set Edge.onShape=false before using
    static const bool DIST = false;      // set Edge.onShape=false before using
    static const bool SPHERE = false;    // set Edge.onShape=false before using

    BinnedData *bdata;
    list<Node *> nodes;        // nodes set
    list<Edge *> edges;        // edges set
    list<Triangle *> triangles;    // triangles set
    list<Edge *> mstEdges;     // minimum spanning tree set
    list<Edge> nearestEdges;     // nearest-neighbor edges
    Edge *hullStart;    // entering edge of convex hull
    Edge *actE;
    Node *start, *end;    // for graph diameter calculations
    int totalPeeledCount;
    int totalCount;
    double alphaArea, alphaPerimeter, hullArea, hullPerimeter;  // initialized in constructor
    double totalOriginalMSTLengths;
    double totalMSTOutlierLengths;
    double *sortedOriginalMSTLengths;

  int np;
    int *px, *py, *counts;
    bool *isOutlier;
    static const double FUZZ = .999;

    void findOutliers(BinnedData bdata);
  void computeDT(int *px, int *py);

  void updateMSTEdges(Edge *addEdge, list<Edge *> *mstEdges);
  void updateMSTNodes(Node *addNode, list<Node *> *mstNodes);
  void computeMST();
  double *getSortedMSTEdgeLengths();
  double computeCutoff(int n, double* lengths);
  void computeTotalOriginalMSTLengths(int n);
  bool computeMSTOutliers(double omega);

  double computeAlphaValue();
  void markShape();
  void computeAlphaGraph(); // requires initializing Edge.onShape = false
  bool edgeIsExposed(double alpha, Edge *e);
  bool pointsInCircle(Node *n, double xc, double yc, double radius);
  
  void insert(int px, int py, int count, int id);

    int searchEdge(Edge *e, Node *nd);

  void expandTri(Edge *e, Node *nd, int type);
  void swapTest(Edge *e11);
  void expandHull(Node *nd);
  void setNeighbors();
  void markHull();

  void computeAlphaArea();
  void computeHullArea();
  void computeAlphaPerimeter();
  void computeHullPerimeter();
    void computeTotalCount() {
    totalCount = 0;
        for (int i = 0; i < np; i++) {
            totalCount += counts[i];
        }
    }

  double *computeMeasures();
  double computeMSTEdgeLengthSkewnessMeasure();
  double computeStringyMeasure();
  double computeClusterMeasure();
  double computeMonotonicityMeasure();
  double computePearson(int n, double *x, double *y, double *weights);
  double computeSparsenessMeasure();
  double computeStriationMeasure();
  double computeConvexityMeasure();
  double computeSkinnyMeasure();
    double computeOutlierMeasure() {
        return totalMSTOutlierLengths / totalOriginalMSTLengths;
    }


  Edge *getAdjacentMSTEdge(Node *n, Edge *e);
  double cosineOfAdjacentEdges(Edge *e1, Edge *e2, Node *n);
  void clearVisits();
  
  public:
    static const int numScagnostics = 9;
    static const int OUTLYING = 0, SKEWED = 1, CLUMPY = 2, SPARSE = 3,
    STRIATED = 4, CONVEX = 5, SKINNY = 6, STRINGY = 7;
  static const int MONOTONIC = 8;

  static const char* scagnosticsLabels[9];
  
    Triangulation() {
    alphaArea = 1, alphaPerimeter = 1, hullArea = 1, hullPerimeter = 1;
    clear();
    }

    void clear() {
        nodes.clear();
        edges.clear();
        triangles.clear();
        mstEdges.clear();
    }

    double *compute(BinnedData bdata, bool quick); 
  

  

/*   

    private void computeGabrielGraph() {  // requires initializing Edge.onShape = true
        Iterator i = edges.iterator();
        while (i.hasNext()) {
            Edge e = (Edge) i.next();
            if (e.onShape) {
                double xc = (double) (e.p1.x + e.p2.x) / 2.;
                double yc = (double) (e.p1.y + e.p2.y) / 2.;
                double radius = e.weight / 2.;
                if (pointsInCircle(e.p1, xc, yc, radius) ||
                        pointsInCircle(e.p2, xc, yc, radius)) {
                    e.deleteSimplex();

                }
            }
        }
    }

    private void computeNearestNeighborGraph() {  // requires initializing Edge.onShape = false
        nearestEdges = new Arraylist();
        Iterator i = nodes.iterator();
        while (i.hasNext()) {
            Node n = (Node) i.next();
            Edge e = n.shortestEdge(true);
            if (e != NULL)
                nearestEdges.add(e);
        }
    }

    private void computeDistanceGraph(double distance) {  // requires initializing Edge.onShape = false
        nearestEdges = new Arraylist();
        Iterator i = nodes.iterator();
        while (i.hasNext()) {
            Node n = (Node) i.next();
            Iterator j = n.neighbors.iterator();
            while (j.hasNext()) {
                Edge e = (Edge) j.next();
                if (e.weight < distance) {
                    nearestEdges.add(e);
                }
            }
        }
    }

    private void computeSphereOfInfluenceGraph() {  // requires initializing Edge.onShape = false
        nearestEdges = new Arraylist();
        Iterator i = nodes.iterator();
        while (i.hasNext()) {
            Node ni = (Node) i.next();
            Edge ei = ni.shortestEdge(true);
            if (ei == NULL)
                continue;
            Iterator j = nodes.iterator();
            while (j.hasNext()) {
                Node nj = (Node) j.next();
                Edge ej = nj.shortestEdge(true);
                if (ej == NULL)
                    continue;
                Edge e = new Edge(ni, nj);
                if (e.weight <= ei.weight + ej.weight)
                    nearestEdges.add(e);
            }
        }
    }

    private void computeRelativeNeighborhoodGraph() {  // requires initializing Edge.onShape = true
        Iterator i = edges.iterator();
        while (i.hasNext()) {
            Edge e = (Edge) i.next();
            if (e.onShape) {
                double radius = e.weight;
                if (pointsInLune(e, radius)) {
                    e.deleteSimplex();
                }
            }
        }
    }

    private bool pointsInCircle(Node n, double xc, double yc, double radius) {
        double r = FUZZ * radius;
        Iterator i = n.neighbors.iterator();
        while (i.hasNext()) {
            Edge e = (Edge) i.next();
            Node no = e.otherNode(n);
            double dist = no.distToNode(xc, yc);
            if (dist < r)
                return true;
        }
        return false;
    }

    private bool pointsInLune(Edge e, double radius) {
        double r = FUZZ * radius;
        Iterator i = e.p1.neighbors.iterator();
        while (i.hasNext()) {
            Edge ne = (Edge) i.next();
            Node no = ne.otherNode(e.p1);
            if (no.distToNode(e.p1.x, e.p1.y) < r &&
                    no.distToNode(e.p2.x, e.p2.y) < r) {
                return true;
            }
        }
        Iterator j = e.p2.neighbors.iterator();
        while (j.hasNext()) {
            Edge ne = (Edge) j.next();
            Node no = ne.otherNode(e.p2);
            if (no.distToNode(e.p1.x, e.p1.y) < r &&
                    no.distToNode(e.p2.x, e.p2.y) < r) {
                return true;
            }
        }
        return false;
    }

    private void computeAlphaGraph() { // requires initializing Edge.onShape = false
        bool deleted;
        double alpha = computeAlphaValue();
        do {
            Iterator i = edges.iterator();
            deleted = false;
            while (i.hasNext()) {
                Edge e = (Edge) i.next();
                if (e.inT.onComplex) {
                    if (alpha < e.weight / 2) {
                        e.inT.onComplex = false;
                        deleted = true;
                    } else {
                        if (e.invE != NULL)
                            if (e.invE.inT.onComplex)
                                continue;
                        if (!edgeIsExposed(alpha, e)) {
                            e.inT.onComplex = false;
                            deleted = true;
                        }
                    }
                }
            }
        } while (deleted);
        markShape();
    }

    private void markShape() {
        Iterator i = edges.iterator();
        while (i.hasNext()) {
            Edge e = (Edge) i.next();
            e.onShape = false;
            if (e.inT.onComplex) {
                if (e.invE == NULL) {
                    e.onShape = true;
                } else if (!e.invE.inT.onComplex)
                    e.onShape = true;
            }
        }
    }

    private bool edgeIsExposed(double alpha, Edge e) {
        double x1 = e.p1.x;
        double x2 = e.p2.x;
        double y1 = e.p1.y;
        double y2 = e.p2.y;
        double xe = (x1 + x2) / 2;
        double ye = (y1 + y2) / 2;
        double d = Math.sqrt(alpha * alpha - e.weight * e.weight / 4);
        double xt = d * (y2 - y1) / e.weight;
        double yt = d * (x2 - x1) / e.weight;
        double xc1 = xe + xt;
        double yc1 = ye - yt;
        double xc2 = xe - xt;
        double yc2 = ye + yt;
        bool pointsInCircle1 = pointsInCircle(e.p1, xc1, yc1, alpha) ||
                pointsInCircle(e.p2, xc1, yc1, alpha);
        bool pointsInCircle2 = pointsInCircle(e.p1, xc2, yc2, alpha) ||
                pointsInCircle(e.p2, xc2, yc2, alpha);
        if (pointsInCircle1 && pointsInCircle2)
            return false;
        return true;
    }


    private Edge getAdjacentMSTEdge(Node n, Edge e) {
        Iterator nt = n.neighbors.iterator();
        while (nt.hasNext()) {
            Edge et = (Edge) nt.next();
            if (et.onMST && !e.equals(et)) {
                return et;
            }
        }
        return NULL;
    }

    private double cosineOfAdjacentEdges(Edge e1, Edge e2, Node n) {
        double v1x = e1.otherNode(n).x - n.x;
        double v1y = e1.otherNode(n).y - n.y;
        double v2x = e2.otherNode(n).x - n.x;
        double v2y = e2.otherNode(n).y - n.y;
        double v1 = Math.sqrt(v1x * v1x + v1y * v1y);
        double v2 = Math.sqrt(v2x * v2x + v2y * v2y);
        v1x = v1x / v1;
        v1y = v1y / v1;
        v2x = v2x / v2;
        v2y = v2y / v2;
        return v1x * v2x + v1y * v2y;
    }

    private double computeConvexityMeasure() {
        if (hullArea == 0) // points in general position
            return 1;
        else {
            double t = (double) totalCount / 500;
            double correction = .7 + .3 / (1 + t * t);
            double convexity = alphaArea / hullArea;
            return correction * convexity;
        }
    }

    private double computeSkinnyMeasure() {
        if (alphaPerimeter > 0)
            return 1 - Math.sqrt(4 * Math.PI * alphaArea) / alphaPerimeter;
        else
            return 1;
    }

*/
  friend class Edge;
  friend class Node;
  friend class Triangle;
};

/*
ostream& operator << (ostream& os, const Node& p) {
//  return os << '(' << p.x << ',' << p.y << ',' << p.pointID << ')';
  return os << p.pointID;
}
ostream& operator << (ostream& os, const Edge& e) {
//  return os << *(e.p1) << "(" << e.p1->x <<","<<e.p1->y << ")-" << *(e.p2) << "(" << e.p2->x <<","<<e.p2->y << ")";
  return os << *(e.p1) << "-" << *(e.p2);
}
ostream& operator << (ostream& os, const Triangle& t) {
//  return os << *(t.anEdge) << "-" << *(t.anEdge->nextE) << "-" << *(t.anEdge->nextE->nextE);
  return os << t.anEdge->p1->pointID << "-" << t.anEdge->nextE->p1->pointID << "-" << t.anEdge->nextE->nextE->p1->pointID;
}
*/
#endif 
