#include <bits/stdc++.h>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

class kl_partition
{
    private:
        struct Node
        {
            // node designator
            int node;
            // vector of adjacent nodes
            vector<int> adjacent;
            // total number of adjacent nodes
            int numAdjNodes;
            // initial costs
            int internalCost;
            int externalCost;
            // d-value
            int d;
            union
            {
                int Cxa;
                int Cya;
            };
            union
            {
                int Cxb;
                int Cyb;
            };
            int currentPartition;
        };

        // structure to contain gain data for each pair of possible node exchanges
        struct NodePair
        {
            // node from each partition being compared
            int p1Node;
            int p2Node;
            // connection between p1 & p2 node pair
            int Cab;
            // current gain of swapping pair
            int gain;
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

        bool verifyData(vector<Node> &nodes, int n, int e, int totalNumNodesRead, int totalNumAdjNodesRead)
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
            else if (totalNumAdjNodesRead != (2 * e))
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

            cout << "data was verified successfully!" << endl;
            cout << "\n";

            return 1;
        }

        vector<int> *initialPartition(vector<Node> &nodes, vector<NodePair> &nodePairs, int *initialCost, int *maxGain, int *p1NodeSwap, int *p2NodeSwap, int n)
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

            // set total number of nodes (including dummy) and n/2
            totalNumNodes = nodes.size() - 1;
            partitionOneSize = totalNumNodes / 2;

            // create partitions
            for (int i = 1; i <= partitionOneSize; i++)
            {
                partitions[0].push_back(nodes[i].node);
                nodes[i].currentPartition = 1;
            }

            for (int i = partitionOneSize + 1; i <= totalNumNodes; i++)
            {
                partitions[1].push_back(nodes[i].node);
                nodes[i].currentPartition = 2;
            }

            // sort partitions
            sort(partitions[0].begin(), partitions[0].end());
            sort(partitions[1].begin(), partitions[1].end());

            // calculate external cost for nodes in partition 1
            for (int i = 1; i <= partitionOneSize; i++)
            {
                nodes[i].externalCost = 0;

                for (int j = partitionOneSize + 1; j <= totalNumNodes; j++)
                {
                    for (int k = 0; k < nodes[j].adjacent.size(); k++)
                    {
                        if (nodes[i].node == nodes[j].adjacent[k])
                        {
                            nodes[i].externalCost += 1;

                            externalEdgeCount++;
                        }
                    }
                }
            }

            // calculate internal cost for nodes in partition 1
            for (int i = 1; i <= partitionOneSize; i++)
            {
                nodes[i].internalCost = 0;

                for (int j = 1; j <= partitionOneSize; j++)
                {
                    for (int k = 0; k < nodes[j].adjacent.size(); k++)
                    {
                        if ((nodes[i].node == nodes[j].adjacent[k]) && (j != i))
                        {
                            nodes[i].internalCost += 1;

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
                nodes[i].externalCost = 0;

                for (int j = 1; j <= partitionOneSize; j++)
                {
                    for (int k = 0; k < nodes[j].adjacent.size(); k++)
                    {
                        if (nodes[i].node == nodes[j].adjacent[k])
                        {
                            nodes[i].externalCost += 1;
                        }
                    }
                }
            }

            // calculate internal cost for nodes in partition 2
            for (int i = partitionOneSize + 1; i <= totalNumNodes; i++)
            {
                nodes[i].internalCost = 0;

                for (int j = partitionOneSize + 1; j <= totalNumNodes; j++)
                {
                    for (int k = 0; k < nodes[j].adjacent.size(); k++)
                    {
                        if ((nodes[i].node == nodes[j].adjacent[k]) && (j != i))
                        {
                            nodes[i].internalCost += 1;

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
            for (int i = 1; i <= partitionOneSize; i++)
            {
                *initialCost += nodes[i].externalCost;
            }

            // calculate cost of second partition for parity
            int costP2 = 0;
            for (int i = partitionOneSize + 1; i <= totalNumNodes; i++)
            {
                costP2 += nodes[i].externalCost;
            }

            // calculate initial gain for each node pair
            for (int i = 1; i <= partitionOneSize; i++)
            {
                for (int j = partitionOneSize + 1; j <= totalNumNodes; j++)
                {
                    NodePair currentNodePair;

                    currentNodePair.p1Node = nodes[i].node;
                    currentNodePair.p2Node = nodes[j].node;
                    currentNodePair.Cab = 0;

                    for (int k = 0; k < nodes[i].adjacent.size(); k++)
                    {
                        if (nodes[j].node == nodes[i].adjacent[k])
                        {
                            currentNodePair.Cab = 1;
                        }
                    }

                    currentNodePair.gain = nodes[i].d + nodes[j].d - (2 * currentNodePair.Cab);

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
                }
            }

            return partitions;
        }

        vector<int> *swapNodes(vector<int> *partitions, vector<NodePair> &nodePairs, vector<Node> &nodes, vector<Node> &usedNodes, int p1Node, int p2Node)
        {
            // replace node designator in partition 1 with node to be swapped from partition 2 and sort in ascending order
            replace(partitions[0].begin(), partitions[0].end(), p1Node, p2Node);
            sort(partitions[0].begin(), partitions[0].end());

            // replace node designator in partition 2 with node to be swapped from partition 1 and sort in ascending order
            replace(partitions[1].begin(), partitions[1].end(), p2Node, p1Node);
            sort(partitions[1].begin(), partitions[1].end());

            // delete all node pairs to rebuild in gain function
            nodePairs.clear();

            // move nodes to used node vector for second kl iteration
            Node tempNode1;
            int usedNode1Index;
            Node tempNode2;
            int usedNode2Index;

            for (int i = 1; i <= (nodes.size() - 1); i++)
            {
                if (nodes[i].node == p1Node)
                {
                    tempNode1 = nodes[i];
                    tempNode1.currentPartition = 2;
                    usedNode1Index = i;
                }
                else if (nodes[i].node == p2Node)
                {
                    tempNode2 = nodes[i];
                    tempNode2.currentPartition = 1;
                    usedNode2Index = i;
                }
            }

            usedNodes.push_back(tempNode1);
            nodes.erase(nodes.begin() + usedNode1Index);

            usedNodes.push_back(tempNode2);
            nodes.erase(nodes.begin() + (usedNode2Index - 1));

            return partitions;
        }

        void updateD(vector<Node> &nodes, int p1Node, int p2Node, vector<int> &p1, vector<int> &p2, vector<int> usedNodes)
        {
            vector<int> abAdjNodes;
            int p1NodeIndex;
            int p2NodeIndex;
            int adjIndex;

            for (int i = 1; i <= (nodes.size() - 1); i++)
            {
                if (nodes[i].node == p1Node)
                {
                    p1NodeIndex = i;
                }
                else if (nodes[i].node == p2Node)
                {
                    p2NodeIndex = i;
                }
            }

            for (int i = 0; i < nodes[p1NodeIndex].adjacent.size(); i++)
            {
                abAdjNodes.push_back(nodes[p1NodeIndex].adjacent[i]);
            }

            for (int i = 0; i < nodes[p2NodeIndex].adjacent.size(); i++)
            {
                abAdjNodes.push_back(nodes[p2NodeIndex].adjacent[i]);
            }

            for (int i = 0; i < abAdjNodes.size(); i++)
            {
                if (find(usedNodes.begin(), usedNodes.end(), abAdjNodes[i]) != usedNodes.end())
                {
                    continue;
                }

                for (int j = 1; j <= (nodes.size() - 1); j++)
                {
                    if (nodes[j].node == abAdjNodes[i])
                    {
                        adjIndex = j;
                        break;
                    }
                }

                if (nodes[adjIndex].currentPartition == 1)
                {
                    nodes[adjIndex].Cxa = 0;
                    nodes[adjIndex].Cxb = 0;

                    for (int k = 0; k < nodes[p1NodeIndex].adjacent.size(); k++)
                    {
                        if (nodes[p1NodeIndex].adjacent[k] == abAdjNodes[i])
                        {
                            nodes[adjIndex].Cxa = 1;
                            break;
                        }
                    }

                    for (int k = 0; k < nodes[p2NodeIndex].adjacent.size(); k++)
                    {
                        if (nodes[p2NodeIndex].adjacent[k] == abAdjNodes[i])
                        {
                            nodes[adjIndex].Cxb = 1;
                            break;
                        }
                    }

                    nodes[adjIndex].d = nodes[adjIndex].d + (2 * nodes[adjIndex].Cxa) - (2 * nodes[adjIndex].Cxb);
                }

                else if (nodes[adjIndex].currentPartition == 2)
                {
                    nodes[adjIndex].Cya = 0;
                    nodes[adjIndex].Cyb = 0;

                    for (int k = 0; k < nodes[p1NodeIndex].adjacent.size(); k++)
                    {
                        if (nodes[p1NodeIndex].adjacent[k] == abAdjNodes[i])
                        {
                            nodes[adjIndex].Cya = 1;
                            break;
                        }
                    }

                    for (int k = 0; k < nodes[p2NodeIndex].adjacent.size(); k++)
                    {
                        if (nodes[p2NodeIndex].adjacent[k] == abAdjNodes[i])
                        {
                            nodes[adjIndex].Cyb = 1;
                            break;
                        }
                    }

                    nodes[adjIndex].d = nodes[adjIndex].d + (2 * nodes[adjIndex].Cyb) - (2 * nodes[adjIndex].Cya);
                }
            }
        }

        void updateGain(vector<NodePair> &nodePairs, vector<Node> &nodes, int *maxGain, int *p1NodeSwap, int *p2NodeSwap, int initialCost)
        {
            *maxGain = initialCost + 1;

            for (int i = 1; i <= (nodes.size() - 1); i++)
            {
                if (nodes[i].currentPartition == 1)
                {
                    for (int j = 1; j <= (nodes.size() - 1); j++)
                    {
                        if (nodes[j].currentPartition == 2)
                        {
                            NodePair currentNodePair;

                            currentNodePair.p1Node = nodes[i].node;
                            currentNodePair.p2Node = nodes[j].node;
                            currentNodePair.Cab = 0;

                            for (int k = 0; k < nodes[i].adjacent.size(); k++)
                            {
                                if (nodes[j].node == nodes[i].adjacent[k])
                                {
                                    currentNodePair.Cab = 1;
                                }
                            }

                            currentNodePair.gain = nodes[i].d + nodes[j].d - (2 * currentNodePair.Cab);

                            // initialize max gain
                            if (*maxGain == (initialCost + 1))
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
                        }
                    }
                }
            }
        }

        void printResults(int iteration, vector<int> &p1, vector<int> &p2, int cost, bool isFirstIteration = 0)
        {

            if (isFirstIteration == 0)
            {
                cout << "iteration " << iteration << endl;
            }
            else if (isFirstIteration == 1)
            {
                cout << "initial partitions " << endl;
            }

            cout << "partition 1: ";
            for (int i : p1)
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

            cout << "cost of the partition: " << cost << endl;
            cout << "\n";
        }

    public:
        
        vector<Node> partition(int argc, char **argv)
        {
            vector<Node> nodes;
            vector<Node> usedNodes;
            vector<int> usedNodeDesignators;
            vector<NodePair> exchangeableNodes;
            vector<int> *partitions;
            ifstream verifiedFileName;
            string newFileName;
            int initialCost = 0;
            int maxGain;
            int p1NodeSwap;
            int p2NodeSwap;
            int iteration = 1;
            // declare variables for total number of nodes and number of edges
            int n = 0;
            int e = 0;
            int totalNumNodesRead = 0;
            int totalNumAdjNodesRead = 0;
            // total gain
            int Gi;
            // current iteration gain
            vector<int> gi;
            int Gk;
            int k;
            int cost;
            // nodes + dummy (1f there is one)
            int actualNodes;
            vector<vector<int>> partitionHistory;

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

            if ((n % 2) != 0)
            {
                actualNodes = n + 1;
            }
            else
            {
                actualNodes = n;
            }

            partitions = initialPartition(nodes, exchangeableNodes, &initialCost, &maxGain, &p1NodeSwap, &p2NodeSwap, n);
            gi.push_back(maxGain);
            Gi = maxGain;
            cost = initialCost - Gi;
            usedNodeDesignators.push_back(p1NodeSwap);
            usedNodeDesignators.push_back(p2NodeSwap);

            printResults(iteration, partitions[0], partitions[1], initialCost, 1);

            updateD(nodes, p1NodeSwap, p2NodeSwap, partitions[0], partitions[1], usedNodeDesignators);

            partitions = swapNodes(partitions, exchangeableNodes, nodes, usedNodes, p1NodeSwap, p2NodeSwap);

            vector<int> tempVec;
            for (int i = 0; i < partitions[0].size(); i++)
            {
                tempVec.push_back(partitions[0][i]);
            }
            partitionHistory.push_back(tempVec);
            tempVec.clear();

            for (int i = 0; i < partitions[1].size(); i++)
            {
                tempVec.push_back(partitions[1][i]);
            }
            partitionHistory.push_back(tempVec);
            tempVec.clear();

            printResults(iteration, partitions[0], partitions[1], cost);

            for (int i = 0; i < ((actualNodes / 2) - 1); i++)
            {
                iteration++;

                updateGain(exchangeableNodes, nodes, &maxGain, &p1NodeSwap, &p2NodeSwap, initialCost);
                gi.push_back(maxGain);
                Gi += maxGain;
                cost = initialCost - Gi;
                usedNodeDesignators.push_back(p1NodeSwap);
                usedNodeDesignators.push_back(p2NodeSwap);

                updateD(nodes, p1NodeSwap, p2NodeSwap, partitions[0], partitions[1], usedNodeDesignators);

                partitions = swapNodes(partitions, exchangeableNodes, nodes, usedNodes, p1NodeSwap, p2NodeSwap);

                for (int i = 0; i < partitions[0].size(); i++)
                {
                    tempVec.push_back(partitions[0][i]);
                }
                partitionHistory.push_back(tempVec);
                tempVec.clear();

                for (int i = 0; i < partitions[1].size(); i++)
                {
                    tempVec.push_back(partitions[1][i]);
                }
                partitionHistory.push_back(tempVec);
                tempVec.clear();

                printResults(iteration, partitions[0], partitions[1], cost);
            }

            int maxK = 0;
            int gRunningSum = 0;
            int gMaxSum = 0;
            for (int i = 0; i < gi.size(); i++)
            {
                gRunningSum += gi[i];

                if (gRunningSum > gMaxSum)
                {
                    gMaxSum = gRunningSum;
                    maxK = i;
                }
            }
    
            maxK = maxK * 2;

            if ((n % 2) != 0)
            {
                int p1MaxElement = *max_element(partitionHistory[maxK].begin(), partitionHistory[maxK].end());

                int p2MaxElement = *max_element(partitionHistory[maxK + 1].begin(), partitionHistory[maxK + 1].end());

                if (p1MaxElement > p2MaxElement)
                {
                    vector<int>::iterator it = find(partitionHistory[maxK].begin(), partitionHistory[maxK].end(), p1MaxElement);
                    partitionHistory[maxK].erase(partitionHistory[maxK].begin() + (it - partitionHistory[maxK].begin()));
                }
                else
                {
                    vector<int>::iterator it = find(partitionHistory[maxK + 1].begin(), partitionHistory[maxK + 1].end(), p2MaxElement);
                    partitionHistory[maxK + 1].erase(partitionHistory[maxK + 1].begin() + (it - partitionHistory[maxK + 1].begin()));
                }
            }

            cout << "the maximum gain is: " << gMaxSum << endl;
            cout << endl;
            cout << "final partition" << endl;
            cout << "partition one: ";
            for (int i = 0; i < partitionHistory[maxK].size(); i++)
            {
                cout << partitionHistory[maxK][i] << " ";
            }
            cout << endl;
            cout << "partition two: ";
            for (int i = 0; i < partitionHistory[maxK + 1].size(); i++)
            {
                cout << partitionHistory[maxK + 1][i] << " ";
            }
            cout << endl;
            cout << "cost: " << (initialCost - gMaxSum) << endl;
        }
};