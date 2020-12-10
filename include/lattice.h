#pragma once

#include <enums.h>

template<int d>
class MDLatticeParam {
  private :
    static constexpr int Ndim = d;
    std::array<int, d> nd;

    BCType bc;
    OrderType order;

    mutable bool parity;//false for the full lattice

    int vol;
    int volCB;

  public:
    explicit MDLatticeParam() : nd{0}, bc{BCType::Undefined}, order{OrderType::Undefined}, parity(false), vol(0), volCB(0){}

    MDLatticeParam(const std::array<int, d> nd_, BCType bc_ = BCType::Dirichlet, OrderType order_ = OrderType::Lexicographic, bool parity = false) : nd{nd_}, bc{bc_}, order(order_), parity(parity) {
      vol = 1;
      for(auto i : nd) vol *= i;
      volCB = parity ? vol : vol / 2;
    }
    MDLatticeParam(const MDLatticeParam<d> &param) :
	   nd{param.Extents()},
     bc(param.GetBC()),
     order(param.GetOrder()),
     parity(param.GetParity()),
     vol(param.GetVol()),
     volCB(param.GetVolCB()) {}

    ~MDLatticeParam() {}

    constexpr int GetNdim(){return Ndim;}
    auto Extents() const {return nd;   }
    auto Extent(const int i) const {return nd[i];}

    BCType    GetBC()   {return bc;   }
    OrderType GetOrder(){return order;}

    inline bool GetParity() const {return parity;}
    inline int  GetVol()    const {return vol;}
    inline int  GetVolCB()  const {return volCB;}
    inline int  GetSurface(const int i, const int j)  const {return nd[i]*nd[j];}

    inline void SetParity(bool parity_){parity = parity_;}

};

template<int d>
class Grid {
  private :
    const MDLatticeParam<d> &param;
    int nSites;

  public:
    Grid(const MDLatticeParam<d> &param_) : param(param_){
      //number of regs:
      nSites = !param.GetParity() ? param.GetVol() : param.GetVolCB();
    }

    int  NSites() const {return nSites;}
    int  Extent(const int i) const {return param.Extent(i);}
    auto Extents() const {return param.Extents();}

    ~Grid(){}
};

template<typename T, int d, typename Arg>
class MDLattice : public Grid<d> {

  using MDLattParam = MDLatticeParam<d>;

  private :

  std::vector<T> v;
  //Aux fields:
  std::vector<T> tmp1;
  std::vector<T> tmp2;
  //Model specific args
  Arg &arg;//arg includes MDLatticeParam.

  public:
    MDLattice(const MDLattParam &param, Arg &arg) : Grid<d>(param), v(param.GetVol()), tmp1(param.GetVol()), tmp2(param.GetVol()), arg(arg) {}
    MDLattice(const Grid<d> &grid, Arg &arg) : Grid<d>(grid), v(grid.NSites()), tmp1(grid.NSites()), tmp2(grid.NSites()), arg(arg) {}

    ~MDLattice() {}

    Arg &GetArgs() const {return arg;}

    template<typename Policy>
    void SetColdLattice(const Policy &p, const T &&val){
      std::fill(p, v.begin(), v.end(), val);
      std::fill(p, tmp1.begin(), tmp1.end(), val);
      std::fill(p, tmp2.begin(), tmp2.end(), val);
    }

    template<typename Policy>
    void SetHotLattice(const Policy &p, const T &&val){
      std::fill(p, v.begin(), v.end(), val);
      std::fill(p, tmp1.begin(), tmp1.end(), val);
      std::fill(p, tmp2.begin(), tmp2.end(), val);
    }

    //Collection of default accessors:
    inline T& GetLattPoint(const int i){ return v[i];}
    inline T& GetTmp1Point(const int i){ return tmp1[i];}
    inline T& GetTmp2Point(const int i){ return tmp2[i];}

    std::vector<T>& V()   { return v;   }
    std::vector<T>& Tmp1(){ return tmp1;}
    std::vector<T>& Tmp2(){ return tmp2;}
};
