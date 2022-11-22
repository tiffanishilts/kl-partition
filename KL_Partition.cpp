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

// declare variables for total number of nodes and number of edges
int n = 0;
int e = 0;
int totalNumNodesRead = 0;
int totalNumAdjNodesRead = 0;

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

struct NodePair {
    int p1Node = 0;
    int p2Node = 0;
    int gain = 0;
};

// function to read file and return a vector of node structures
vector<Node> readFile(ifstream &file)
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
        totalNumAdjNodesRead += currentNode.numAdjNodes;

        // push current node structure into node vector
        nodes.push_back(currentNode);
        
        // increment node number
        nodeNumber++;
    }

    // assign number of nodes and number of edges from the first line of the file to their respective variables
    n = nodes[0].adjacent[0];
    e = nodes[0].adjacent[1];

    // update total nodes read from file to ensure data is valid
    totalNumNodesRead = nodes.size() - 1;

    // correct for first line
    totalNumAdjNodesRead -= 2;

    if ((totalNumNodesRead % 2) != 0)
    {
        Node dummyNode;

        dummyNode.node = totalNumNodesRead + 1;
        dummyNode.internalCost = 0;
        dummyNode.externalCost = 0;
        
        nodes.push_back(dummyNode);
    }

    // return vector of node structures
    return nodes;
}

// function to verify the data read from the file is valid
// there must be nodes equal to n
// there must be adjacent nodes equal to 2e
bool verifyData(vector<Node>& allNodes)
{
    cout << "the specified number of nodes (n) is " << n << ". the specified number of edges (e) is " << e << "." << endl;
    cout << "the total number of nodes read is " << totalNumNodesRead << ". the total number of adjacent nodes read is " << totalNumAdjNodesRead << "." << endl;
    cout << "\n";

    if (allNodes[0].numAdjNodes != 2)
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

vector<int>* initialPartition(vector<Node>& allNodes, int *cost)
{
    static vector<int> partitions[2];
    int totalNumNodes = 0;
    int partitionOneSize = 0;
    int edgeCount = 0;
    int internalEdgeCount1 = 0;
    int internalEdgeCount2 = 0;
    int externalEdgeCount = 0;

    totalNumNodes = allNodes.size() - 1;
    partitionOneSize = totalNumNodes/2;

    for (int i = 1; i <= partitionOneSize; i++)
    {
        partitions[0].push_back(allNodes[i].node);
    }

    for (int i = partitionOneSize + 1; i <= totalNumNodes; i++)
    {
        partitions[1].push_back(allNodes[i].node);
    }

    sort(partitions[0].begin(), partitions[0].end());
    sort(partitions[1].begin(), partitions[1].end());

    // calculate external cost for nodes in partition 1
    for (int i = 1; i <= partitionOneSize; i++)
    {
        for (int j = partitionOneSize + 1; j <= totalNumNodes; j++)
        {
            for (int k = 0; k < allNodes[j].adjacent.size(); k++)
            {
                cout << "(external) adjacent node is " << allNodes[j].adjacent[k] << endl;

                if (allNodes[i].node == allNodes[j].adjacent[k])
                {
                    allNodes[i].externalCost += 1;

                    cout << "i is " << i << " j is " << j << endl;

                    cout << "external cost of node " << allNodes[i].node << " is now " << allNodes[i].externalCost << endl;

                    externalEdgeCount++;
                }
            }
        }
    }

    // calculate internal cost for nodes in partition 2
    for (int i = 1; i <= partitionOneSize; i++)
    {
        for (int j = 1; j <= partitionOneSize; j++)
        {
            for (int k = 0; k < allNodes[j].adjacent.size(); k++)
            {
                cout << "(internal) adjacent node is " << allNodes[j].adjacent[k] << endl;

                if ((allNodes[i].node == allNodes[j].adjacent[k]) && (j != i))
                {
                    allNodes[i].internalCost += 1;

                    cout << "i is " << i << " j is " << j << endl;

                    cout << "internal cost of node " << allNodes[i].node << " is now " << allNodes[i].internalCost << endl;

                    internalEdgeCount1++;
                }
            }
        }
    }

    internalEdgeCount1 = internalEdgeCount1 / 2;

    // calculate external cost for nodes in partition 2
    for (int i = partitionOneSize + 1; i <= totalNumNodes; i++)
    {
        for (int j = 1; j <= partitionOneSize; j++)
        {
            for (int k = 0; k < allNodes[j].adjacent.size(); k++)
            {
                cout << "(external) adjacent node is " << allNodes[j].adjacent[k] << endl;
                
                if (allNodes[i].node == allNodes[j].adjacent[k])
                {
                    allNodes[i].externalCost += 1;

                    cout << "i is " << i << " j is " << j << endl;

                    cout << "external cost of node " << allNodes[i].node << " is now " << allNodes[i].externalCost << endl;
                }
            }
        }
    }

    // calculate internal cost for nodes in partition 2
    for (int i = partitionOneSize + 1; i <= totalNumNodes; i++)
    {
        for (int j = partitionOneSize + 1; j <= totalNumNodes; j++)
        {
            for (int k = 0; k < allNodes[j].adjacent.size(); k++)
            {
                cout << "(internal) adjacent node is " << allNodes[j].adjacent[k] << endl;

                if ((allNodes[i].node == allNodes[j].adjacent[k]) && (j != i))
                {
                    allNodes[i].internalCost += 1;

                    cout << "i is " << i << " j is " << j << endl;

                    cout << "internal cost of node " << allNodes[i].node << " is now " << allNodes[i].internalCost << endl;

                    internalEdgeCount2++;
                }
            }
        }
    }

    internalEdgeCount2 = internalEdgeCount2 / 2;

    // update d value for all nodes
    for (int i = 1; i <= totalNumNodes; i++)
    {
        allNodes[i].d = allNodes[i].externalCost - allNodes[i].internalCost;
        //cost += allNodes[i].externalCost;
    }


    // calculate cost of first partition
    for(int i = 1; i <= partitionOneSize; i++)
    {
        *cost += allNodes[i].externalCost;
    }

    int costP2 = 0;
    for (int i = partitionOneSize + 1; i <= totalNumNodes; i++)
    {
        costP2 += allNodes[i].externalCost;
    }

    cout << "\n";
    cout << "external cost of p1 is " << *cost << endl;
    cout << "\n";

    cout << "external cost of p2 is " << costP2 << endl;
    cout << "\n";

    edgeCount = externalEdgeCount + internalEdgeCount1 + internalEdgeCount2;
    cout << "total number of edges counted is " << edgeCount << endl;
    cout << "\n";

    return partitions;
}

/*
vector<Node> KLalgorithm()
{
    // stub code
}
*/

// print out iteration summary in the following format:
// iteration #
// partition 1: list of nodes in increasing order
// partition 2: list of nodes in increasing order
// cost of the partition = cost calculated by kl algorithm
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
    vector<int> *partitions;

    ifstream fileName(argv[1]);

    cout << "\n";

    if (fileName.is_open())
    {
        cout << "file was opened successfully!" << endl;
        cout << "\n";
        nodes = readFile(fileName);
    }
    else
    {
        string newFileName;

        while (fileName.is_open() != 1)
        {
            cout << "specified file could not be opened. please provide a valid file name." << endl;
            cout << "\n";
            cin >> newFileName;
            fileName.open(newFileName.c_str());
        }

        cout << "file was opened successfully!" << endl;
        cout << "\n";
        nodes = readFile(fileName);
    }

    ifstream verifiedFileName;
    string newFileName2;

    while (verifyData(nodes) != 1)
    {
        while (verifiedFileName.is_open() != 1)
        {
            cout << "specified file did not provide valid data. please provide a valid file." << endl;
            cout << "\n";
            cin >> newFileName2;
            verifiedFileName.open(newFileName2.c_str());
        }

        cout << "file was opened successfully!" << endl;
        cout << "\n";
        nodes = readFile(verifiedFileName);
    }

    cout << "data was verified successfully!" << endl;
    cout << "\n";
    
    cout << "a valid graph has been obtained from the file. data processing will begin now." << endl;
    cout << "\n";

    int cost = 0;

    partitions = initialPartition(nodes, &cost);

    printResults(1, partitions[0], partitions[1], cost);


    return 0;
}
