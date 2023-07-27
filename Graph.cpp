#include "Graph.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <queue>
#include <limits>
using namespace std;

int getIndex(vector<string> v, string K)
{
    auto it = find(v.begin(), v.end(), K);
    // If element was found
    if (it != v.end()) 
    {
        int index = it - v.begin();
        return index;
    }
    else {
        return -1;
    }
}

Graph::Graph(const char* const & edgelist_csv_fn) {
    // TODO: Initialize a Graph object from a given edge list CSV, where each line `u,v,w` represents an edge between nodes `u` and `v` with weight `w`.

    // 1. create node for each node
    // 2. create pairs for (node, weight) for each node in the total node list
    // 3. add these pairs to vector<pair> adjList
    // 4. num_nodes is adjList.size()
    // 5. num_edges is add all the edges for each node in adjList

    ifstream my_file(edgelist_csv_fn);      // open the file
    string line;                    
    while(getline(my_file, line)) {  // read one line from the file
        stringstream ss(line);      
        string u, v, w; 
        getline(ss, u, ',');     // store first column in "first"
        getline(ss, v, ',');    // store second column in "second"
        getline(ss, w, '\n');    // store third column column in "third"
        
        int i = getIndex(all_nodes, u);
        int j = getIndex(all_nodes, v);
        if(i == -1){
            all_nodes.push_back(u);
            adjList.push_back({});
            i = all_nodes.size() - 1;
        }
        if(j == -1){
            all_nodes.push_back(v);
            adjList.push_back({});
            j = all_nodes.size() - 1;
        }
        edgeList.push_back(make_tuple(stod(w), u, v ));
        adjList[i].push_back(make_pair(v, stod(w)));
        adjList[j].push_back(make_pair(u, stod(w)));
    }
    my_file.close(); 

}

unsigned int Graph::num_nodes() {
    // TODO: Return the number of nodes in this graph.
    return all_nodes.size();
}

vector<string> Graph::nodes() {
    // TODO: Return a `vector` of node labels of all nodes in this graph, in any order.
    return all_nodes;
}

unsigned int Graph::num_edges() {
    // TODO: Return the number of (undirected) edges in this graph.
    unsigned int count = 0;
    for (unsigned int i = 0; i < num_nodes(); i++) {
        count += adjList[i].size();
    }
    return count / 2;
}

unsigned int Graph::num_neighbors(string const & node_label) {
    // TODO: Return the number of neighbors of a given node.
    int count = 0;
    int i = getIndex(all_nodes, node_label);

    if (i != -1) {
        count = adjList[i].size();
    }
    return count;
     
}

double Graph::edge_weight(string const & u_label, string const & v_label) {
    // TODO: Return the weight of the edge between a given pair of nodes, or -1 if there does not exist an edge between the pair of nodes.
    int u_index = getIndex(all_nodes, u_label);
    int v_index = getIndex(all_nodes, v_label);

    if (u_index != -1 && v_index != -1) {
        for (pair<string, double> p : adjList[u_index]) {
            if (p.first == v_label) {
                return p.second;
            }
        }
    }
    return -1;
}

vector<string> Graph::neighbors(string const & node_label) {
    // TODO: Return a `vector` containing the labels of the neighbors of a given node. The neighbors can be in any order within the `vector`.
    vector<string> labels;
    int index = getIndex(all_nodes, node_label);

    if (index != -1) {
        for (pair<string, double> p : adjList[index]) {
            labels.push_back(p.first);
        }
    }
    return labels;
}

vector<string> Graph::shortest_path_unweighted(string const & start_label, string const & end_label) {
    // TODO
    // 1. find the neighbors of the start_label and add to queue
    // 2. for loop through each neighbor until reach end_label and add to queue
    // 3. keep track of the count of the nodes and the previous node for each add 
    // 4. once find the end_label
    queue<string> q;
    // Create a set to mark visited nodes
    unordered_set<string> visited;
    // Create a map to keep track of the previous node for each node
    unordered_map<string, string> previous;
    // Enqueue the start node and mark it as visited
    q.push(start_label);
    visited.insert(start_label);

    // Perform breadth-first search until the end node is found or the queue is empty
    while (!q.empty()) {
        string current = q.front();
        q.pop();

        // Check if the current node is the end node
        if (current == end_label) {
            break;  // Stop the search if the end node is found
        }

        // Get the neighbors of the current node
        vector<string> neighbors = Graph::neighbors(current);
        for (string n : neighbors) {
            if (visited.find(n) == visited.end()) {
                // Enqueue the neighbor if it hasn't been visited
                q.push(n);

                // Mark the neighbor as visited
                visited.insert(n);

                // Set the previous node for the neighbor
                previous[n] = current;
            }
        }
    }
    // If the end node was found, construct the shortest path
    vector<string> shortest_path;
    if (visited.find(end_label) != visited.end()) {
        string current = end_label;
        while (current != start_label) {
            shortest_path.push_back(current);
            current = previous[current];
        }
        shortest_path.push_back(start_label);

        // Reverse the shortest path to get the correct order
        reverse(shortest_path.begin(), shortest_path.end());
    }

    return shortest_path;

}

vector<tuple<string,string,double>> Graph::shortest_path_weighted(string const & start_label, string const & end_label) {
    // TODO
    priority_queue<pair<double, string>, vector<pair<double, string>>, greater<pair<double, string>>> pq;

    // Create a map to store the distance from the start node to each node
    unordered_map<string, double> distance;

    // Create a map to store the previous node in the shortest path
    unordered_map<string, string> previous;

    // Initialize distances to infinity for all nodes except the start node
    for (string node : all_nodes) {
        if (node == start_label) {
            distance[node] = 0.0;
        } else {
            distance[node] = numeric_limits<double>::infinity();
        }
    }

    // Push the start node to the priority queue
    pq.push(make_pair(0.0, start_label));

    while (!pq.empty()) {
        // Get the node with the smallest distance from the priority queue
        string u = pq.top().second;
        pq.pop();

        // If we reach the end node, we have found the shortest path
        if (u == end_label) {
            break;
        }

        // Visit each neighbor of the current node
        for (Pair neighbor : adjList[getIndex(all_nodes, u)]) {
            string v = neighbor.first;
            double weight = neighbor.second;

            // Calculate the new distance from the start node to the neighbor
            double new_distance = distance[u] + weight;

            // If the new distance is smaller than the current distance, update the distance and previous node
            if (new_distance < distance[v]) {
                distance[v] = new_distance;
                previous[v] = u;
                pq.push(make_pair(new_distance, v));
            }
        }
    }

    vector<tuple<string, string, double>> shortest_path;

    // Reconstruct the shortest path from the end node to the start node
    string current_node = end_label;
    while (!current_node.empty() && previous.find(current_node) != previous.end()) {
        string previous_node = previous[current_node];
        double edge_weight = this->edge_weight(previous_node, current_node);
        shortest_path.insert(shortest_path.begin(), make_tuple(previous_node, current_node, edge_weight));
        current_node = previous_node;
    }

    // If there is no shortest path, return an empty vector
    if (shortest_path.empty() && start_label != end_label) {
        return shortest_path;
    }

    // If the start and end nodes are the same, return a single-element vector with the node itself
    if (start_label == end_label) {
        shortest_path.push_back(make_tuple(start_label, start_label, -1.0));
    }

    return shortest_path;
}

vector<vector<string>> Graph::connected_components(double const & threshold) {
    // TODO
    vector<vector<string>> connected;
    unordered_set<string> visited;

    for (string node_label : all_nodes) {
        if (visited.find(node_label) == visited.end()) {
            vector<string> component;
            dfs(node_label, threshold, visited, component);
            connected.push_back(component);
        }
    }

    return connected;
}

void Graph::dfs(const string& node_label, const double& threshold, unordered_set<string>& visited, vector<string>& component) {
    visited.insert(node_label);
    component.push_back(node_label);

    vector<string> neighbors = this->neighbors(node_label);
    for (string neighbor : neighbors) {
        double weight = this->edge_weight(node_label, neighbor);
        if (weight <= threshold && visited.find(neighbor) == visited.end()) {
            dfs(neighbor, threshold, visited, component);
        }
    }
}

double Graph::smallest_connecting_threshold(string const &start_label, string const &end_label) {
    if(start_label == end_label){
        return 0.0;
    }
    unordered_map<string, int> rank;
    // creating a disjoint set for all nodes of the graph
    for (string sen : all_nodes) {
        parent[sen] = sen;
        rank[sen] = 0;
    }

    sort(edgeList.begin(), edgeList.end());
    
    // Perform Union-Find algorithm to find the smallest connecting threshold
    for (tuple<double, string, string> edge: edgeList) {
        double weight = get<0>(edge);
        string u = get<1>(edge);
        string v = get<2>(edge);

        // Find the parent sets of nodes u and v
        string u_set = find(u);
        string v_set = find(v);

        // Perform union operation
        union_string(rank, u_set, v_set);

        // Check if the start and end nodes are now in the same set
        if (find(start_label) == find(end_label)) {
            return weight;
        }
    }
    return -1.0;  // No path exists
}

string Graph::find(string x) {
    if (parent[x] == x)
        return x;
    return find(parent[x]);
}

void Graph::union_string(unordered_map<string, int> &rank, string x, string y) {
    string xset = find(x);
    string yset = find(y);
    if(xset == yset){
        return;
    } 

    int xrank = rank[xset];
    int yrank = rank[yset];

    if(xrank < yrank){
        parent[xset] = yset;
    }
    else if(yrank < xrank){
        parent[yset] = xset;
    }
    else{
        parent[xset] = yset;
        rank[yset] ++;
    }
    
}
