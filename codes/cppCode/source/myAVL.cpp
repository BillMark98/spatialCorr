#include "myAVL.h"

MyAVL::~MyAVL() {
    deleteNode(root);
}
long MyAVL::height(Node * nd) {
    if(nd == nullptr) {
        return 0;
    }
    else {
        return nd -> height;
    }
}

unsigned long MyAVL::count(Node * nd) {
    if (nd == nullptr) {
        return 0;
    }
    return nd -> count + nd -> leftCount + nd -> rightCount;
}

unsigned long MyAVL::totalCount() {
    if (root == nullptr) {
        return 0;
    }
    return root -> count + root -> leftCount + root -> rightCount;
}

void MyAVL::update(Node * nd) {
    updateHeight(nd);
    updateCount(nd);
}

void MyAVL::updateHeight(Node * nd) {
    nd -> height = 1 + max(height(nd -> leftNode), height(nd -> rightNode));
    // since the updateBalance function also call the height function
    // can actually directly use these two variables in the function without
    // calling a sub function, dunno if it will make the insertion way faster.
    updateBalance(nd);
}
void MyAVL::updateCount(Node * nd) {
    nd -> leftCount = count(nd -> leftNode);
    nd -> rightCount = count(nd -> rightNode);
}
int MyAVL::balance(Node * nd) {
    int bal = height(nd -> leftNode) - height(nd -> rightNode);
    return bal;
}
Node * MyAVL::rightRotate(Node * y) {
    /*
    //           y                   x
    //         /   \                /  \
    //       x       y2   ->       x1    y 
    //      /  \                        /  \
    //     x1   z                      z     y2

    */
    Node * x = y -> leftNode;
    Node * z = x -> rightNode;
    
    x -> rightNode = y;
    y -> leftNode = z;
    
    y -> height = max(height(y -> leftNode), height(y -> rightNode)) + 1;
    x -> height  = max(height(x -> leftNode), height(x -> rightNode)) + 1;

    y -> leftCount = count(z);
    x -> rightCount = count(y);
    return x;
}

Node * MyAVL::leftRotate(Node * y) {
    /*
    //           y                   x
    //         /   \                /  \
    //       y2      x   ->        y    x2 
    //              /  \          / \      
    //             z   x2        y2  z      
    */   
    Node * x = y -> rightNode;
    Node * z = x -> leftNode;
    
    x -> leftNode = y;
    y -> rightNode = z;
    y -> height = max<long>(height(y -> leftNode), height(y -> rightNode)) + 1;
    x -> height  = max<long>(height(x -> leftNode), height(x -> rightNode)) + 1;
    
    y -> rightCount = count(z);
    x -> leftCount = count(y);
    return x;
}

void MyAVL::updateBalance(Node * y) {
    y -> balance = height(y -> leftNode) - height(y -> rightNode);
}

Node * MyAVL::insert(Node * nd, myDouble dist) {
    if(nd == nullptr) {
        return new Node(dist);
    }

    if(dist < nd -> distance) {
        nd -> leftNode =  insert(nd -> leftNode,dist);
    }
    else if(dist > nd -> distance) {
        nd -> rightNode = insert(nd -> rightNode, dist);
    }
    else {
        nd -> count++;
        return nd;
    }

    // update the height,balance, count of the subtree with root at nd
    update(nd);

    // get the balance of the tree
    int balance = nd -> balance;

    // four cases in AVL-tree
    // note that if abs(balance) > 1 means that the nd -> left or rightNode in 
    // the following cases is not nullptr
    // 1. left-left case
    if(balance > 1 && dist < nd -> leftNode -> distance) {
        return rightRotate(nd);
    }
    
    // 2. right-right case
    if(balance < -1 && dist > nd -> rightNode -> distance) {
        return leftRotate(nd);
    }

    /*
    // 3. left-right case
    //    nd              nd              dist
    //   /                /               /  \  
    //  x       ->       dist      ->    x   nd
    //   \              /
    //     dist        x
    */
    if(balance > 1 && dist > nd -> leftNode -> distance) {
        nd -> leftNode = leftRotate(nd -> leftNode);
        return rightRotate(nd);
    }

    /*
    // 4. right-left case
    //     nd               nd                      dist
    //       \                 \                    /  \
   //         x    ->          dist      ->        nd   x
   //        /                   \
   //       dist                  x
   */
   if(balance < -1 && dist < nd -> rightNode -> distance) {
       nd -> rightNode = rightRotate(nd -> rightNode);
       return leftRotate(nd);
   }

   // if nd unchanged
   return nd;
}
Node * insert(Node * nd, Node * newNode) {
    exit(TODO_FUNCTION);
    return nullptr;
}

unsigned long MyAVL::getPair(double theta, double deltaTheta) {
    // since we consider the interval [theta, theta + deltaTheta)
    // which equals #{p \in tree | p's key >= theta} - #{p \in tree | p's key >= theta + deltaTheta}
    unsigned long left = geqDist(root, theta);
    unsigned long right = geqDist(root, theta + deltaTheta);
    return left - right;
}
template <typename T>
const T& max(const T &v1, const T &v2){
    return v1 > v2 ? v1 : v2;
}

void MyAVL::deleteNode(Node * nd) {
    if(nd == nullptr) {
        return;
    }
    deleteNode(nd -> leftNode);
    deleteNode(nd -> rightNode);
    delete nd;
}

unsigned long MyAVL::geqDist(Node * nd, myDouble theta) {
    if (nd == nullptr) {
        return 0;
    }
    /*
                x
               /
              y 
            /  \
           u    z
               / \ 
              w   p

        i) at y if y's key geq theta, then all the points in the right subtree
        (starting from z) and y itself will be counted
        just need to look at the lefttree

        ii) if y's key less than theta, then all the points in the left and y itself
        does not satisfy >= theta, so will not be counted,
        only need to consider the righttree
    */
    if (theta <= nd -> distance) {
        return geqDist(nd -> leftNode, theta) + nd -> count + nd -> rightCount;
    }
    else {
        return geqDist(nd -> rightNode, theta);
    }
}

unsigned long MyAVL::leqDist(Node * nd, myDouble theta) {
    if (nd == nullptr) {
        return 0;
    }
    /*
                x
               /
              y 
            /  \
           u    z
               / \ 
              w   p

        i) at y if y's key leq theta, then all the points in the left subtree
        (starting from z) and y itself will be counted
        just need to look at the righttree

        ii) if y's key greater than theta, then all the points in the right and y itself
        does not satisfy <= theta, so will not be counted,
        only need to consider the lefttree
    */
    if (theta >= nd -> distance) {
        return leqDist(nd -> rightNode, theta) + nd -> count + nd -> leftCount;
    }
    else {
        return leqDist(nd -> leftNode, theta);
    }
}

void MyAVL::traverse(Node * nd, nodeFunc f, int order) {
    if (nd == nullptr) {
        return;
    }
    if (order == 1) {
       traverse(nd -> leftNode, f,1);
       f(nd);
       traverse(nd -> rightNode,f,1);
    }
    else if (order == 0) {
        f(nd);
        traverse(nd -> leftNode,f,0);
        traverse(nd -> rightNode,f,0);
    }
    else {
        traverse(nd -> leftNode,f,2);
        traverse(nd -> rightNode,f,2);
        f(nd);
    }
}

// TO DO. can MyAVL use the function displayNode?
// i guess yes since the function does not need explicitly be called
// like   node.displayNode. The function acts as a friend function
void MyAVL::displayTree() {
    traverse(root, displayNode);
}

void displayNode(Node * nd) {
    cout << "distance: " <<  nd -> distance << "\t count: " << nd -> getCount()
        << "\t height: " << nd -> getHeight()  << "\t balance: " << nd -> getBalance()<< endl;
}

istream & operator>>(istream & is, MyAVL & myAVL) {
    if(myAVL.root != nullptr) {
        myAVL.deleteNode(myAVL.root);
    }
    myDouble dist;
    is >> dist;
    myAVL.root = new Node(dist);
    while(is.good()) {
        is >> dist;
        myAVL.root = myAVL.insert(myAVL.root,dist);
    }
    return is;
}

MyAVL & MyAVL::readVector(const std::vector<SphericalPoint> & vec) {
    typedef std::vector<SphericalPoint>::const_iterator vec3DConstIter;
    // unsigned long len = vec.size();
    if (root != nullptr) {
        deleteNode(root);
    }
    for(vec3DConstIter iter = vec.begin(); iter != vec.end(); iter++) {
        for (vec3DConstIter jter = std::next(iter,1); jter != vec.end(); jter++) {
            root = insert(root,(*iter).sphericalAngleDist((*jter)));
        }
    }
    return (*this);
}


long Node::getHeight() {
    long left,right;
    if(leftNode == nullptr) {
        left = 0;
    }
    else {
        left = leftNode -> getHeight();
    }

    if (rightNode == nullptr) {
        right = 0;
    }
    else {
        right = rightNode -> getHeight();
    }
    return max(left,right) + 1;
}
int Node::getBalance() {
    long left,right;
    if(leftNode == nullptr) {
        left = 0;
    }
    else {
        left = leftNode -> getHeight();
    }

    if (rightNode == nullptr) {
        right = 0;
    }
    else {
        right = rightNode -> getHeight();
    }
    return left - right;
}