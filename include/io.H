// I/O utilities


#ifndef IO_H
#define IO_H


#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
namespace nmlib{


/******** Prototype ********/

// map<K,V> I/O (format: n K1 V1 ... Kn Vn)
template<class K,class V> std::ostream& operator<<(std::ostream& str, const std::map<K,V>& k2v);
template<class K,class V> std::istream& operator>>(std::istream& str,       std::map<K,V>& k2v);

// vector<T> I/O (format: n X1 ... Xn)
template<class T> std::ostream& operator<<(std::ostream& str, const std::vector<T>& xx);
template<class T> std::istream& operator>>(std::istream& str,       std::vector<T>& xx);

// string <--> T conversion (via ostream<<T and istream>>T operators)
template<class T> void        str2any(const std::string& str, T& val);
template<class T> T           str2any(const std::string& str);
template<class T> std::string any2str(const T& val,int prec=0);


/******** Implementation ********/

template<class K,class V> std::ostream& operator<<(std::ostream& str, const std::map<K,V>& k2v){
  str << k2v.size() << "\n";
  for(typename std::map<K,V>::const_iterator it=k2v.begin(); it!=k2v.end(); it++)
    str << it->first << "\t" << it->second << "\n";
  return str;
}
template<class K,class V> std::istream& operator>>(std::istream& str, std::map<K,V>& k2v){
  size_t n,i;
  k2v.clear();
  str >> n;
  for(i=0; i<n; i++){ K key; V val; str>>key>>val;  k2v[key]=val; }
  return str;
}

template<class T> std::ostream& operator<<(std::ostream& str, const std::vector<T>& xx){
  str << xx.size() << '\n';
  for(size_t i=0; i<xx.size(); i++) str << xx[i] << '\n';
  return str;
}
template<class T> std::istream& operator>>(std::istream& str, std::vector<T>& xx){
  size_t n,i;
  str >> n;
  xx=std::vector<T>(n);
  for(i=0; i<n; i++) str >> xx[i];
  return str;
}

template<class T> void str2any(const std::string& str, T& val){
  std::stringstream ss;
  ss<<str;
  ss>>val;
  if(ss.fail()) throw std::runtime_error("str2any: bad string");
}
template<class T> T    str2any(const std::string& str){
  T val;
  str2any(str,val);
  return val;
}
template<class T> std::string any2str(const T& val, int prec){
  std::stringstream ss;
  if(prec>0){ ss.precision(prec); ss<<std::scientific; }
  ss<<val;
  return ss.str();
}


}  //namespace nmlib
#endif //IO_H
