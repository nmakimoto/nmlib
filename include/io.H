// I/O utilities


#ifndef IO_H
#define IO_H


#include <iostream>
#include <map>
#include <vector>


/******** prototype ********/

// map<K,V> I/O (format: n K1 V1 ... Kn Vn)
template<class K,class V> std::ostream& operator<<(std::ostream& str, const std::map<K,V>& k2v);
template<class K,class V> std::istream& operator>>(std::istream& str,       std::map<K,V>& k2v);

// vector<T> I/O (format: n X1 ... Xn)
template<class T> std::ostream& operator<<(std::ostream& str, const std::vector<T>& xx);
template<class T> std::istream& operator>>(std::istream& str,       std::vector<T>& xx);


/******** implementation ********/

template<class K,class V> std::ostream& operator<<(std::ostream& str, const std::map<K,V>& k2v){
  str << k2v.size() << "\n";
  for(typename std::map<K,V>::const_iterator it=k2v.begin(); it!=k2v.end(); it++)
    str << it->first << "\t" << it->second << "\n";
  return str;
}
template<class K,class V> std::istream& operator>>(std::istream& str, std::map<K,V>& k2v){
  int n,i;
  k2v.clear();
  str >> n;
  for(i=0; i<n; i++){ K key; V val; str>>key>>val;  k2v[key]=val; }
  return str;
}

template<class T> std::ostream& operator<<(std::ostream& str, const std::vector<T>& xx){
  str << xx.size() << '\n';
  for(int i=0; i<xx.size(); i++) str << xx[i] << '\n';
  return str;
}
template<class T> std::istream& operator>>(std::istream& str, std::vector<T>& xx){
  int n,i;
  str >> n;
  xx=std::vector<T>(n);
  for(i=0; i<n; i++) str >> xx[i];
  return str;
}


#endif //IO_H
