#ifndef _SIMPLEGRAPH_HPP
#define _SIMPLEGRAPH_HPP
#endif

#include <vector>
#include <map> 
#include <iostream>

using namespace std;

template <class T>
class Node {
public:
  Node(T d) { this->data = d; marked = false; }
  ~Node() {}

  vector<Node*> edges;
  bool marked;
  T data;  
private:

 
};
 
template <class T>
class simpleGraph {

public:
  simpleGraph() {}
  
  ~simpleGraph() {
    for (typename vector<Node<T> *>::iterator it = this->nodes.begin(); it != this->nodes.end(); it++) {
      delete *it;
    }

  }
  void addEdge(T x, T y);

  // printGraph();

  vector<Node<T> *> nodes;

  struct dereference_compare {
    template <class I>
    bool operator()(const I& a, const I& b) {
        return *a < *b;
    }
  };

  map<T, Node<T> *, dereference_compare> nodeMap;
  void connectedComponents(int minResultSize, long *clusterNum);
  list<Node<T> *> * bfs(Node<T> *n);
};


template <class T>
void simpleGraph<T>::addEdge(T x, T y) {
  
  if (this->nodeMap.find(x) == this->nodeMap.end()) {
    this->nodeMap[x] = new Node<T>(x);
    this->nodes.push_back(this->nodeMap[x]);
    }
  
  if (this->nodeMap.find(y) == this->nodeMap.end()) {
    this->nodeMap[y] = new Node<T>(y);
    this->nodes.push_back(this->nodeMap[y]);
  }

  this->nodeMap[x]->edges.push_back(this->nodeMap[y]);
  this->nodeMap[y]->edges.push_back(this->nodeMap[x]);
  
}

template <class T>
void simpleGraph<T>::connectedComponents(int minResultSize, long *clusterNum) {

  for (typename vector<Node<T> *>::iterator nodeIt = this->nodes.begin(); nodeIt != this->nodes.end(); nodeIt++ ) {
    if (!(*nodeIt)->marked) {
	list<Node<T> *> * result = this->bfs(*nodeIt);
	if (result->size() >= minResultSize) {
	(*clusterNum)++;
	for (typename list<Node<T> *>::iterator it = result->begin();  it != result->end(); it++) {
	  cout << *clusterNum << "\t" << result->size() << "\t";
	  (*(*it)->data)->print(); 
	}
	cout << endl;
	}
	delete result;
  }
}
    
  //  results.push_back(dfs())
  
}

// breadth-first search to get connected component containing the specified node
// returns a list of Node<T> pointers
template <class T>
list<Node<T> *> * simpleGraph<T>::bfs(Node<T> *n) {
  list<Node<T> *> * results = new list<Node<T> *>;
  results->push_back(n);
  n->marked = true;  
  for(typename vector<Node<T> *>::iterator destIt =  n->edges.begin(); destIt != n->edges.end(); destIt++) {
    if (!(*destIt)->marked) {
      list<Node<T> *> * childResults = bfs(*destIt);
      results->splice(results->end(), *childResults);
      delete childResults;
    }
  }

  return results;
}
