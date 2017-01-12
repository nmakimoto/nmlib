// Shortest path finder


#ifndef ROUTE_H
#define ROUTE_H


#include <vector>
#include <map>
#include <set>
#include <algorithm>
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
    auto it=prev.find(n1);
    if(it==prev.end() || it->first==it->second) return Path();
    n1=it->second;
  }
  std::reverse(p.begin(),p.end());
  return p;
}


inline Route::N2N Route::dijkstra(Node n0) const{
  N2N  prev;
  N2D  cost;
  auto cmp = [&](Node a, Node b){ return cost[a]<cost[b]; };
  std::set<Node,decltype(cmp)> next(cmp);

  cost[n0]=0;
  next.insert(n0);

  while(!next.empty()){
    Node k0=*next.begin();
    Dist d0=cost[k0];
    next.erase(k0);
    //if(k0==n1) break;  // shortest path to goal

    if(dist.find(k0)==dist.end()) continue;
    const N2D& n2d=dist.at(k0);
    for(const auto& it: n2d){
      Node k1=it.first;
      Dist d1=it.second+d0;
      if(cost.find(k1)!=cost.end() && cost[k1]<=d1) continue;
      next.erase(k1);  // for consitency: erase k1 and update cost before re-inserting k1
      cost[k1]=d1;
      prev[k1]=k0;
      next.insert(k1);
    }
  }
  return prev;
}


}  //end of namespace nmlib
#endif //ROUTE_H
