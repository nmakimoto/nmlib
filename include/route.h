// Shortest path finder


#ifndef ROUTE_H
#define ROUTE_H


#include <vector>
#include <map>
#include <queue>
#include <algorithm>
#include <stdexcept>
namespace nmlib{


/******** Prototype ********/

class Route {
public:
  typedef double               Dist;
  typedef int                  Node;
  typedef std::vector<Node>    Path;
  typedef std::map<Node,Node>  N2N;
  typedef std::map<Node,Dist>  N2D;
  typedef std::map<Node,N2D>   NN2D;

  NN2D dist;  // dist[n0][n1] := length of edge (n0,n1)

  Dist length    (const Path& path) const;  // path length (-1 if invalid)
  Dist length    (Node n0, Node n1) const;  // edge length (-1 if invalid)
  Path operator()(Node n0, Node n1) const;  // shortest path from n0 to n1
  Path operator()(Node n0, Node n1, const N2N& prev) const;
  N2N  dijkstra  (Node n0) const;           // shortest-path tree of n0
};


/******** Implementation ********/

inline Route::Dist Route::length(const Path& p) const{
  Dist len=0;
  for(size_t k=0; k+1<p.size(); k++){
    Dist len1=length(p[k],p[k+1]);
    if(len1<0) return Dist(-1);
    len+=len1;
  }
  return len;
}
inline Route::Dist Route::length(Node n0, Node n1) const{
  if(n0==n1) return Dist(0);
  auto i0=dist.find(n0);        if(i0==dist.end())       return Dist(-1);
  auto i1=i0->second.find(n1);  if(i1==i0->second.end()) return Dist(-1);
  return i1->second;
}


inline Route::Path Route::operator()(Node n0, Node n1) const{
  return operator()(n0,n1,dijkstra(n0));
}
inline Route::Path Route::operator()(Node n0, Node n1, const N2N& prev) const{
  if(prev.find(n1)==prev.end()) return Path();
  Path p;
  while(true){
    p.push_back(n1);
    if(n1==n0) break;
    if(n1==prev.at(n1)) return Path();
    n1=prev.at(n1);
  }
  std::reverse(p.begin(),p.end());
  return p;
}


inline Route::N2N Route::dijkstra(Node n_src) const{
  N2N  prev;  // shortest-path tree of node n_src
  N2D  cost;
  auto cmp = [&](Node n1, Node n2){ return cost[n1]>cost[n2]; };
  std::priority_queue<Node,std::vector<Dist>,decltype(cmp)> bdry(cmp);

  prev[n_src]=n_src;
  cost[n_src]=0;
  bdry.push(n_src);
  if(dist.find(n_src)==dist.end()) return prev;

  while(!bdry.empty()){
    Node k1=bdry.top(), k0=prev.at(k1);  // (k0,k1) = argmin_{j0:visited, j1:boundary} cost(j0)+length(j0,j1)
    bdry.pop();
    prev[k1]=k0;

    if(dist.find(k1)==dist.end()) continue;
    const N2D& n2d=dist.at(k1);
    for(const auto& it: n2d){
      Node k2=it.first;
      Dist d =it.second+cost.at(k1);
      if(prev.find(k2)!=prev.end()) continue;
      if(cost.find(k2)!=cost.end() && cost[k2]<=d) continue;
      prev[k2]=k1;
      cost[k2]=d;  // set/update priority
      bdry.push(k2);
    }
  }
  return prev;
}


}  //end of namespace nmlib
#endif //ROUTE_H
