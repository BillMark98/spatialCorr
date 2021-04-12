#include "myDefineConst.h"
#include "sphericalPoint.h"
#include <vector>

// typedef std::vector<SphericalPoint> spherePointsContainer;
typedef std::vector<const SphericalPoint &> spherePointsContainer;
class TreeNode {
    public:
        spherePointsContainer pointContainer;
        TreeNode * leftNode;
        TreeNode * rightNode;
        char turns;

        TreeNode(SphericalPoint& p, short turn);
        bool operator<(TreeNode& tn);
        bool operator==(TreeNode& tn);
        bool operator>(TreeNode&);

        bool isEmpty() const;
};

class TwoPointKdTree {
    private:
        TreeNode rootNode;
    public:
        void putPoint(SphericalPoint & p);
        
};