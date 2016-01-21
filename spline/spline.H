// スプライン補間


#ifndef SPLINE_H
#define SPLINE_H


#include <vector>
#include <map>
#include <algorithm>


/******** クラス一覧 ********/

template<class T> class spline;
typedef spline<double> Spline;


/******** クラスI/F ********/

template<class T> class spline{
public:
  spline  (void);
  spline  (const std::map<T,T>& x2y);
  void set(const std::map<T,T>& x2y);  // 制御点指定
  std::map<T,T> get(void) const;       // 制御点取得

  T operator()(T x, int deg=3) const;  // 補間値計算

private:
  std::vector<T> xx,yy,aa,bb,cc,dd;  // 制御点, スプライン係数
  void init_abcd(void);              // スプライン係数の計算
};


/******** インライン実装 ********/

template<class T> spline<T>::spline(void){}
template<class T> spline<T>::spline(const std::map<T,T>& x2y){ set(x2y); }
template<class T> void spline<T>::set(const std::map<T,T>& x2y){
  xx=yy=aa=bb=cc=dd=std::vector<T>();
  for(typename std::map<T,T>::const_iterator p=x2y.begin(); p!=x2y.end(); p++){
    xx.push_back(p->first);
    yy.push_back(p->second);
  }
  if(x2y.size()>=3) init_abcd();
}
template<class T> std::map<T,T> spline<T>::get(void) const{
  std::map<T,T> x2y;
  for(int i=0; i<xx.size(); i++) x2y[xx[i]]=yy[i];
  return x2y;
}

// 補間値計算
template<class T> T spline<T>::operator()(T x, int deg) const{
  int n=xx.size();
  int k=std::upper_bound(xx.begin(),xx.end(),x)-xx.begin()-1;  // x∈[xx[k],xx[k+1])

  if(n<3) return (n==0 ? 0 : n==1 ? yy[0] : yy[0]+(yy[1]-yy[0])*(x-xx[0])/(xx[1]-xx[0]));  // 不足
  k=(k<0 ? 0 : n-2<k ? n-2 : k);  // 範囲外
  x-=xx[k];

  if(deg<=0) return (xx[k]+x<xx[k+1]-x ? yy[k] : yy[k+1]);  // 0次(最近傍)
  if(deg<=1) return yy[k]+(yy[k+1]-yy[k])*x/(xx[k+1]-xx[k]);  // 1次(線形)
  return ((aa[k]*x+bb[k])*x+cc[k])*x+dd[k];  // スプライン
}

// スプライン係数の計算
// 制御点 (x0,y0)...(xn,yn) を補間する区分多項式 Sk(x)=ak(x-xk)^3+bk(x-xk)^2+ck(x-xk)+dk (k=0..n-1)の係数
template<class T> void spline<T>::init_abcd(void){
  int n=xx.size()-1;  // 区間数 = 点数-1
  std::vector<T> vv(n+1),pp(n),qq(n),rr(n);

  // 後述の(3)"を解いてb[1]..b[n-1]を得る
  for(int j=1; j<n; j++){
    pp[j] = 2 * (xx[j+1] - xx[j-1]);  // 係数行列Aの対角成分
    qq[j] = xx[j+1]-xx[j];  // Aの下側副対角成分
    rr[j] = xx[j]-xx[j-1];  // Aの上川副対角成分
    vv[j] = 3 * ((yy[j+1]-yy[j])/(xx[j+1]-xx[j]) - (yy[j]-yy[j-1])/(xx[j]-xx[j-1]));  // 初期値=(3)"右辺
  }
  for(int j=1; j<n-1; j++){ T s=qq[j]/pp[j];  vv[j+1]-=s*vv[j];  qq[j]=0;  pp[j+1]-=s*rr[j+1]; }  // Aの下側をsweep
  for(int j=n-1; j>1; j--){ T s=rr[j]/pp[j];  vv[j-1]-=s*vv[j];  rr[j]=0; }  // Aの上側をsweep
  for(int j=1; j<n  ; j++){ vv[j]/=pp[j];  pp[j]=1; }  // Aの対角を正規化
  vv[0]=vv[n]=0;  // (5)'

  // (1)"(2)"(4)"でabdを得る
  aa=bb=cc=dd=std::vector<T>(n);
  for(int k=0; k<n; k++){
    T b0=vv[k], b1=vv[k+1], h=xx[k+1]-xx[k], g=yy[k+1]-yy[k];
    aa[k] = (b1-b0)/(3*h);
    bb[k] = b0;
    cc[k] = g/h-(b1+2*b0)*h/3;
    dd[k] = yy[k];
  }
}
// 係数ak,bk,ck,dkの算出要領は下記。
//
// 制御点を通る(以下x[k+1]-x[k]をh[k]と略記):
//   (1) S[k](0)   =y[k]         for  k=0..n-1
//   (2) S[k](h[k])=y[k+1]       for  k=〃
// 制御点で滑らか:
//   (3) S[k]'(h[k])=S[k+1]'(0)  for  k=0..n-2
//   (4) S[k]"(h[k])=S[k+1]"(0)  for  k=〃
// 両端点で上下に凸でない:
//   (5) S[0]"(0)=S[n-1]"(h[n-1])=0
//
// 展開すると  (以下○[k]を○,○[k+1]を○'などと略記)
//   (1)'                       d =  y   for  k=0..n-1
//   (2)'   ah^3 +  bh^2 + ch + d =  y'  for  k=〃
//   (3)' 3 ah^2 + 2bh   + c      =  c'  for  k=0..n-2
//   (4)' 6 ah   + 2b             = 2b'  for  k=〃
//   (5)' b0=bn=0  (bn:=S[n-1]"(h[n-1])は便宜的に導入)
//
// bについて解くと (以下x'-xをh,y'-yをgなどと略記)
//   (1)"  d = y                 for  k=0..n-1
//   (4)"  a = (b'-b)/3h         for  k=0..n-1  (bn=0としてk=n-1でもOK)
//   (2)"  c = g/h - (b'+2b)h/3  for  k=0..n-1
//   (3)"  bh + 2b'(h'+h) + b"h' = 3{g'/h' - g/h}  for  k=0..n-2
//
// (3)"は3重対角行列を係数とする未知数b[1]..b[n-1] の連立一次方程式であり容易に解ける。
// 残りのa,c,dは(4)"(2)"(1)"で求まる。


#endif //SPLINE_H
