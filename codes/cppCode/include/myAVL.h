#include "distContainer.h"
// #include "poissonPoint2d.h"  // for the typedef definition of vec_3DCoord
#include "sphericalPoint.h"
// #include "sphereTwoCorrelation.h" // for the definition of spherePointsContainer

// #include <iostream>


#ifndef __MYAVL_H
#define __MYAVL_H
class Node {
    public:
    // number of subnodes (including self)
    // note that in this case it is possible that many distances have
    // the same value, so that each node represent a bunch of distances
    // of the same value.
    unsigned long count;
    // number of all values in the left subnodes (multiplicity will be treated
    // indepenedently)
    unsigned long leftCount;
    // number of all values in the right subnodes.
    unsigned long rightCount;
    Node * leftNode;
    Node * rightNode;
    Node * prevNode;
    unsigned long height;
    int balance;
    myDouble distance;
    Node() {}
    Node(myDouble dist) : count(1),leftCount(0),rightCount(0),
    leftNode(nullptr),rightNode(nullptr),height(1),distance(dist){}
    void setDist(myDouble dist) { distance = dist;}
    long getCount() const { return count;}

    // // get the height of the node nd
    long getHeight();
    int getBalance();

    bool operator<(const Node & n1) const;
    bool operator>(const Node & n2) const;
    bool operator==(const Node & n2) const;

    // friend void displayNode(Node * nd);
};

typedef void (*nodeFunc)(Node *);
class MyAVL : DistContainer {
    private:
    Node * root;

    MyAVL(const MyAVL &);
    public:
    MyAVL() : root(nullptr){}
    MyAVL(myDouble distance){root = new Node(distance);}
    ~MyAVL();
    // get the height of the node nd
    long height(Node * nd);
    // get the counts (including the node itself)
    unsigned long count(Node * nd);
    // get the total values within the tree
    unsigned long totalCount();
    // update height, balance and count of the node nd
    void update(Node * nd);
    // update the height of the node nd
    void updateHeight(Node * nd);
    // update the counts
    void updateCount(Node * nd);
    // get the balance of the node
    int balance(Node * nd);
    // basic AVL-tree action
    Node * rightRotate(Node * y);
    Node * leftRotate(Node * y);
    // update the balance of the node y
    void updateBalance(Node * y);
    Node * insert(Node * nd, myDouble dist);
    Node * insert(Node * nd, Node * newNode);

    // specific for the distance graph
    // here use the interval [theta, theta + deltaTheta)
    unsigned long getPair(double theta, double deltaTheta);
    void putValue(double val) { exit(TODO_FUNCTION);}
    // help function
    // delete Node
    void deleteNode(Node * nd);
    // find the number of distances greater than or equal with theta starting from
    // node nd
    unsigned long geqDist(Node * nd, myDouble theta);
    // find the number of distances less than or equal with theta, starting 
    // from node nd
    unsigned long leqDist(Node * nd, myDouble theta);
    // traversing the tree
    // 0 pre-order 1 in-order 2 post-order
    void traverse(Node * nd, nodeFunc f, int order = 1);
    void displayTree();

    // read in file
    friend istream & operator>>(istream & is, MyAVL & myAVL);
    // read in vectors
    MyAVL & readVector(const std::vector<SphericalPoint> & vec);
    
};

template <typename T>
const T& max(const T &v1, const T &v2);

// template <>
// const long & max<long>(const long &v1, const long &v2);

void displayNode(Node * nd);
#endif