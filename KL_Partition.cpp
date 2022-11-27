/*
    Copyright Tiffani Shilts Portland, OR November 2022

    name: KL_Partition.cpp

    abstract: a cpp program that implements a Kernighan-Lin partitioning algorithm on data provided in the chaco input file format

    inputs:

    outputs:

    description:
*/


#include <bits/stdc++.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

// structure for each node which holds it's adjacent nodes
struct Node {
    int node;
    vector<int> adjacent;
    int numAdjNodes = 0;
    bool swapped = 0;
    int internalCost = 0;
    int externalCost = 0;
    int d = 0;
};

// structure to contain gain data for each pair of possible node exchanges
struct NodePair {
    int p1Node = 0;
    int p2Node = 0;
    int gain = 0;
};

// input:
// read file and return a vector of node structures
// output:
vector<Node> readFile(ifstream &file, int *n, int *e, int *totalNumNodesRead, int *totalNumAdjNodesRead)
{
    // declare vector of nodes
    static vector<Node> nodes;
    // declare string variables for read operation
    string line;
    string word;
    // declare int variable to hold next adjacent node
    int nextAdjNode;
    // initialize node number. node "0" is the graph information
    int nodeNumber = 0;

    // read file until eof
    while (getline(file, line))
    {
        // create Node instance for current node
        Node currentNode;
        // bind istringstream to current line in order to use input operator
        istringstream adjacentNodes(line);

        // use input operator to push each adjacent node for the current node into the structure vector member
        while (adjacentNodes >> word)
        {
            nextAdjNode = stoi(word);
            currentNode.adjacent.push_back(nextAdjNode);
        }

        // assign node number and number of adjacent nodes to current node
        currentNode.node = nodeNumber;
        currentNode.numAdjNodes = currentNode.adjacent.size();

        // update running total of adjacent nodes in file to ensure the data is valid
        *totalNumAdjNodesRead += currentNode.numAdjNodes;

        // push current node structure into node vector
        nodes.push_back(currentNode);
        
        // increment node number
        nodeNumber++;
    }

    // assign number of nodes and number of edges from the first line of the file to their respective variables
    *n = nodes[0].adjacent[0];
    *e = nodes[0].adjacent[1];

    // update total nodes read from file to ensure data is valid
    *totalNumNodesRead = nodes.size() - 1;

    // correct for first line
    *totalNumAdjNodesRead -= 2;

    // if the number of nodes is not even, create a dummy node
    if ((*totalNumNodesRead % 2) != 0)
    {
        Node dummyNode;

        dummyNode.node = *totalNumNodesRead + 1;
        
        nodes.push_back(dummyNode);
    }

    // return vector of node structures
    return nodes;
}

// input:
// function to verify the data read from the file is valid
// there must be nodes equal to n
// there must be adjacent nodes equal to 2e
// output:
bool verifyData(vector<Node>& nodes, int n, int e, int totalNumNodesRead, int totalNumAdjNodesRead)
{
    // inform user of number of nodes stated/read and number of edges stated/read
    cout << "the specified number of nodes (n) is " << n << ". the specified number of edges (e) is " << e << "." << endl;
    cout << "the total number of nodes read is " << totalNumNodesRead << ". the total number of adjacent nodes read is " << totalNumAdjNodesRead << "." << endl;
    cout << "\n";

    // error check for: incorrect first line data, incorrect number of adjacent nodes read, and incorrect number of nodes read
    if (nodes[0].numAdjNodes != 2)
    {
        cout << "this file does not have valid graph data. the file must include two numbers on the first line, number of nodes, and number of edges, respectively. please try another file." << endl;
        cout << "\n";
        return 0;
    }
    else if (totalNumAdjNodesRead != (2*e))
    {
        cout << "this file does not have a valid number of adjacent nodes. the graph must include 2e adjacent nodes. please try another file." << endl;
        cout << "\n";
        return 0;
    }
    else if (totalNumNodesRead != n)
    {
        cout << "this file does not have a valid number of nodes. the graph must include n nodes. please try another file." << endl;
        cout << "\n";
        return 0;
    }

    return 1;
}

// input:
// functionality:
// swap the node designators in the vectors representing the partitions
// reorder partitions in ascending order
// remove all node pairs containing either node from the list of possible nodes to swap
// output: 
vector<int>* swapNodes(vector<int> *partitions, vector<NodePair>& nodePairs, int p1Node, int p2Node)
{
    // replace node designator in partition 1 with node to be swapped from partition 2 and sort in ascending order
    replace(partitions[0].begin(), partitions[0].end(), p1Node, p2Node);
    sort(partitions[0].begin(), partitions[0].end());

    // replace node designator in partition 2 with node to be swapped from partition 1 and sort in ascending order
    replace(partitions[1].begin(), partitions[1].end(), p2Node, p1Node);
    sort(partitions[1].begin(), partitions[1].end());

    for (int i = 0; i < nodePairs.size(); i++)
    {
        if ((p1Node == nodePairs[i].p1Node) || (p2Node == nodePairs[i].p2Node))
        {
            nodePairs.erase(nodePairs.begin() + i);
        }
    }

    return partitions;
}

// input: reference to vector of node structures and pointer to int for cost
// functionality:
// create initial partition of (1,...,n/2) and ((n/2)+1,...,n) nodes
// calculate internal cost, external cost, and d for each node
// calculate cost of initial partition
// output: array of two vectors containing the node designators of the nodes in the initial p1 and p2 partitions
vector<int>* initialPartition(vector<Node>& nodes, vector<NodePair>& nodePairs, int *cost, int *maxGain, int *p1NodeSwap, int *p2NodeSwap, int n)
{
    // array of vectors to contain data about which nodes are in each partition
    static vector<int> partitions[2];
    // total number of nodes + dummy node (if there is one)
    int totalNumNodes = 0;
    // n/2
    int partitionOneSize = 0;
    // total number of edges counted while calculating cost -- for parity
    int edgeCount = 0;
    // p1 internal edge count
    int internalEdgeCount1 = 0;
    // p2 internal edge count
    int internalEdgeCount2 = 0;
    // external edge count (p1 and p2 should be equal)
    int externalEdgeCount = 0;
    // variable for the number of connections between nodes a & b
    int Cab = 0;

    // set total number of nodes (including dummy) and n/2
    totalNumNodes = nodes.size() - 1;
    partitionOneSize = totalNumNodes/2;

    // create partitions
    for (int i = 1; i <= partitionOneSize; i++)
    {
        partitions[0].push_back(nodes[i].node);
    }

    for (int i = partitionOneSize + 1; i <= totalNumNodes; i++)
    {
        partitions[1].push_back(nodes[i].node);
    }

    // sort partitions
    sort(partitions[0].begin(), partitions[0].end());
    sort(partitions[1].begin(), partitions[1].end());

    // calculate external cost for nodes in partition 1
    for (int i = 1; i <= partitionOneSize; i++)
    {
        for (int j = partitionOneSize + 1; j <= totalNumNodes; j++)
        {
            for (int k = 0; k < nodes[j].adjacent.size(); k++)
            {
                // cout << "(external) adjacent node is " << nodes[j].adjacent[k] << endl;

                if (nodes[i].node == nodes[j].adjacent[k])
                {
                    nodes[i].externalCost += 1;

                    // cout << "i is " << i << " j is " << j << endl;

                    // cout << "external cost of node " << nodes[i].node << " is now " << nodes[i].externalCost << endl;

                    externalEdgeCount++;
                }
            }
        }
    }

    // calculate internal cost for nodes in partition 1
    for (int i = 1; i <= partitionOneSize; i++)
    {
        for (int j = 1; j <= partitionOneSize; j++)
        {
            for (int k = 0; k < nodes[j].adjacent.size(); k++)
            {
                // cout << "(internal) adjacent node is " << nodes[j].adjacent[k] << endl;

                if ((nodes[i].node == nodes[j].adjacent[k]) && (j != i))
                {
                    nodes[i].internalCost += 1;

                    // cout << "i is " << i << " j is " << j << endl;

                    // cout << "internal cost of node " << nodes[i].node << " is now " << nodes[i].internalCost << endl;

                    internalEdgeCount1++;
                }
            }
        }
    }

    // one edge per two nodes counted
    internalEdgeCount1 = internalEdgeCount1 / 2;

    // calculate external cost for nodes in partition 2
    for (int i = partitionOneSize + 1; i <= totalNumNodes; i++)
    {
        for (int j = 1; j <= partitionOneSize; j++)
        {
            for (int k = 0; k < nodes[j].adjacent.size(); k++)
            {
                // cout << "(external) adjacent node is " << nodes[j].adjacent[k] << endl;
                
                if (nodes[i].node == nodes[j].adjacent[k])
                {
                    nodes[i].externalCost += 1;

                    // cout << "i is " << i << " j is " << j << endl;

                    // cout << "external cost of node " << nodes[i].node << " is now " << nodes[i].externalCost << endl;
                }
            }
        }
    }

    // calculate internal cost for nodes in partition 2
    for (int i = partitionOneSize + 1; i <= totalNumNodes; i++)
    {
        for (int j = partitionOneSize + 1; j <= totalNumNodes; j++)
        {
            for (int k = 0; k < nodes[j].adjacent.size(); k++)
            {
                // cout << "(internal) adjacent node is " << nodes[j].adjacent[k] << endl;

                if ((nodes[i].node == nodes[j].adjacent[k]) && (j != i))
                {
                    nodes[i].internalCost += 1;

                    // cout << "i is " << i << " j is " << j << endl;

                    // cout << "internal cost of node " << nodes[i].node << " is now " << nodes[i].internalCost << endl;

                    internalEdgeCount2++;
                }
            }
        }
    }

    // one edge per two nodes counted
    internalEdgeCount2 = internalEdgeCount2 / 2;

    // update d value for all nodes
    for (int i = 1; i <= totalNumNodes; i++)
    {
        nodes[i].d = nodes[i].externalCost - nodes[i].internalCost;
    }

    // calculate cost of first partition
    for(int i = 1; i <= partitionOneSize; i++)
    {
        *cost += nodes[i].externalCost;
    }

    // calculate cost of second partition for parity
    int costP2 = 0;
    for (int i = partitionOneSize + 1; i <= totalNumNodes; i++)
    {
        costP2 += nodes[i].externalCost;
    }

    // print for error checking
    // cout << "\n";
    cout << "external cost of p1 is " << *cost << endl;
    cout << "\n";

    cout << "external cost of p2 is " << costP2 << endl;
    cout << "\n";

    edgeCount = externalEdgeCount + internalEdgeCount1 + internalEdgeCount2;
    cout << "total number of edges counted is " << edgeCount << endl;
    cout << "\n";

    
    //calculate initial gain for each node pair
    for (int i = 1; i <= partitionOneSize; i++)
    {
        for (int j = partitionOneSize + 1; j <= totalNumNodes; j++)
        {
            NodePair currentNodePair;

            currentNodePair.p1Node = nodes[i].node;
            currentNodePair.p2Node = nodes[j].node;

            for (int k = 0; k < nodes[i].adjacent.size(); k++)
            {
                if (nodes[j].node == nodes[i].adjacent[k])
                {
                    Cab = 1;
                }
            }
            
            currentNodePair.gain = nodes[i].d + nodes[j].d - (2 * Cab);

            // initialize max gain
            if ((i == 1) && (j == (partitionOneSize + 1)))
            {
                *maxGain = currentNodePair.gain;
                *p1NodeSwap = nodes[i].node;
                *p2NodeSwap = nodes[j].node;
            }
            // check if the current gain is larger than the max gain, and update if so
            else
            {
                if (currentNodePair.gain > *maxGain)
                {
                    *maxGain = currentNodePair.gain;
                    *p1NodeSwap = nodes[i].node;
                    *p2NodeSwap = nodes[j].node;
                }
            }

            nodePairs.push_back(currentNodePair);

            // add code to keep a running total of the max gain, and node pair that corresponds to it
        }
    }

    // parity check
    cout << "there should be " << ((n/2) * (n/2)) << " node pairs." << endl;
    cout << "\n";
    cout << "there are " << nodePairs.size() << " node pairs." << endl;
    cout << "\n";

    // temp max gain code to finish checkpoint. move in to original loop which calculates gain

    vector<int> gains;
    int maxGainCheck;

    cout << "the gains are ";

    for(int i = 0; i < nodePairs.size(); i++)
    {
        cout << (i + 1) << ": " << nodePairs[i].gain << " ";
        gains.push_back(nodePairs[i].gain);
    }
    
    cout << endl;
    cout << "\n";

    maxGainCheck = *max_element(gains.begin(), gains.end());

    cout << "the maximum gain is " << maxGainCheck << endl;
    cout << "\n";

    return partitions;
}

/*
vector<Node> KLalgorithm()
{
    
}
*/

// input:
// print out iteration summary in the following format:
// iteration #
// partition 1: list of nodes in increasing order
// partition 2: list of nodes in increasing order
// cost of the partition: calculated by kl algorithm
// output:
void printResults(int iteration, vector<int> &p1, vector<int> &p2, int cost, bool isFirstIteration = 0)
{
    cout << "iteration " << iteration << endl;

    cout << "partition 1: ";
    for (int i: p1)
    {
        cout << i << " ";
    }
    cout << endl;

    cout << "partition 2: ";
    for (int i : p2)
    {
        cout << i << " ";
    }
    cout << endl;

    if(isFirstIteration == 0)
    {
        cout << "cost of the partition: " << cost << endl;
        cout << "\n";
    }
    else if (isFirstIteration == 1)
    {
        cout << "cost of the partition: unknown" << endl;
        cout << "\n";
    }   
}


int main(int argc, char **argv)
{
    vector<Node> nodes;
    vector<NodePair> exchangeableNodes;
    vector<int> *partitions;
    ifstream verifiedFileName;
    string newFileName;
    int cost = 0;
    int maxGain;
    int p1NodeSwap;
    int p2NodeSwap;
    int iteration = 1;
    // declare variables for total number of nodes and number of edges
    int n = 0;
    int e = 0;
    int totalNumNodesRead = 0;
    int totalNumAdjNodesRead = 0;

    ifstream fileName(argv[1]);

    cout << "\n";

    if (fileName.is_open())
    {
        cout << "file was opened successfully!" << endl;
        cout << "\n";
        nodes = readFile(fileName, &n, &e, &totalNumNodesRead, &totalNumAdjNodesRead);
    }
    else
    {
        while (fileName.is_open() != 1)
        {
            cout << "specified file could not be opened. please provide a valid file name." << endl;
            cin >> newFileName;
            fileName.open(newFileName.c_str());
            cout << "\n";
        }

        cout << "file was opened successfully!" << endl;
        cout << "\n";
        nodes = readFile(fileName, &n, &e, &totalNumNodesRead, &totalNumAdjNodesRead);
    }

    while (verifyData(nodes, n, e, totalNumNodesRead, totalNumAdjNodesRead) != 1)
    {
        while (verifiedFileName.is_open() != 1)
        {
            cout << "specified file did not provide valid data. please provide a valid file." << endl;
            cin >> newFileName;
            verifiedFileName.open(newFileName.c_str());
            cout << "\n";
        }

        cout << "file was opened successfully!" << endl;
        cout << "\n";
        nodes = readFile(verifiedFileName, &n, &e, &totalNumNodesRead, &totalNumAdjNodesRead);
    }

    cout << "data was verified successfully!" << endl;
    cout << "\n";
    
    cout << "a valid graph has been obtained from the file. data processing will begin now." << endl;
    cout << "\n";

    partitions = initialPartition(nodes, exchangeableNodes, &cost, &maxGain, &p1NodeSwap, &p2NodeSwap, n);

    printResults(iteration, partitions[0], partitions[1], cost);

    cout << "the maximum gain is " << maxGain << endl;
    cout << "\n";

    cout << "the partition one node to be swapped is " << p1NodeSwap << endl;
    cout << "\n";

    cout << "the partition two node to be swapped is " << p2NodeSwap << endl;
    cout << "\n";

    iteration++;

    partitions = swapNodes(partitions, exchangeableNodes, p1NodeSwap, p2NodeSwap);
    
    cost = maxGain;

    printResults(iteration, partitions[0], partitions[1], cost);

    return 0;
}