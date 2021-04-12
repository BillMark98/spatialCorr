#include "twoPointkdTree.h"

TreeNode::TreeNode(SphericalPoint& p, short turn) {
    pointContainer.resize(1);
    pointContainer.push_back(p);
    turns = turn;
}

// still needd to modify
bool TreeNode::operator<(TreeNode& tn) {
    if(tn.isEmpty()) {
        cerr << "cant compare treenode < with an empty node. " << endl;
        exit(TREENODE_EMPTY);
    }
    // compare phi
    if(turns == 0) {
        double tnPhi = tn.pointContainer[0].getPhi();
        double phi = pointContainer[0].getPhi();
        double diff = tnPhi - phi;
        if(diff > 0 && diff < M_PI) {
            return true;
        }
        else {
            return false;
        }
    }
    // compare theta
    else if (turns == 1) {
        double tnTheta = tn.pointContainer[0].getTheta();
        double theta = pointContainer [0].getTheta();
        double diff = tnTheta - theta;
        if(diff > 0) {
            return true;
        }
        else {
            return false;
        }
    }
}
bool TreeNode::operator==(TreeNode& tn) {
    if (tn.isEmpty()) {
        cerr << "cant compare treenode == with an empty node. " << endl;
        exit(TREENODE_EMPTY);
    }
    // compare phi
    if (turns == 0) {
        double tnPhi = tn.pointContainer[0].getPhi();
        double phi = pointContainer[0].getPhi();
        double diff = tnPhi - phi;
        if (fabs(diff) < MY_EPSILON) {
            return true;
        }
        else {
            return false;
        }
    }
    // compare theta
    else if (turns == 1) {
        double tnTheta = tn.pointContainer[0].getTheta();
        double theta = pointContainer [0].getTheta();
        double diff = tnTheta - theta;
        if (fabs(diff) < MY_EPSILON) {
            return true;
        }
        else {
            return false;
        }
    }
}
bool TreeNode::operator>(TreeNode& tn) {
    if (tn.isEmpty()) {
        cerr << "cant compare treenode > with an empty node. " << endl;
        exit(TREENODE_EMPTY);
    }
    // compare phi
    if (turns == 0) {
        double tnPhi = tn.pointContainer[0].getPhi();
        double phi = pointContainer[0].getPhi();
        double diff = tnPhi - phi;
        if (fabs(diff) < MY_EPSILON) {
            return true;
        }
        else {
            return false;
        }
    }
    // compare theta
    else if (turns == 1) {
        double tnTheta = tn.pointContainer[0].getTheta();
        double theta = pointContainer [0].getTheta();
        double diff = tnTheta - theta;
        if (fabs(diff) < MY_EPSILON) {
            return true;
        }
        else {
            return false;
        }
    }
}

bool TreeNode::isEmpty() const {
    return pointContainer.size() == 0;
}