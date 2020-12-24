#pragma once

#include <enums.h>
#include <iterators.h>
#include <stencil_params.h>
#include <stencil_impl.h>

//requires c++20 flag
//template <typename Float, typename Policy, StencilType ST, StencilPolicy SP, int ND, template<typename data_t, StencilType stencil_t> class Args, std::enable_if<std::is_floating_point<Float>::value, bool> = true>
template <typename Float, typename Policy, StencilType ST, StencilPolicy SP, int ND, template<typename data_t, StencilType stencil_t> class Args>
class FwdEulerIters{
  private:
    Policy &policy;

    const Args<Float, ST> &args;

    std::vector<Float> &v1;
    std::vector<Float> &v2;

    const int outer_range;
    const int inner_range;

  public:
    FwdEulerIters(Policy &policy, const Args<Float, ST> &args, std::vector<Float> &f1, std::vector<Float> &f2, const int outer_range, const int inner_range) :
      policy(policy),
      args(args),
      v1(f1),
      v2(f2),
      outer_range(outer_range),
      inner_range(inner_range) {

      }

    void apply(const int nsteps){
      //Create stencil functor instances:
      std::unique_ptr<GenericNDStencil<Float, ST, ND>> even_t_func_ptr(new GenericNDStencil<Float, ST, ND>(v1, v2, args.c0, args.c1, args.lattice.Extents(), inner_range));
      std::unique_ptr<GenericNDStencil<Float, ST, ND>> odd_t_func_ptr (new GenericNDStencil<Float, ST, ND>(v2, v1, args.c0, args.c1, args.lattice.Extents(), inner_range));
      //launch iterations
      for(int i = 0; i < nsteps; i++) {
        auto &func = (i & 1) == 0 ? *even_t_func_ptr : *odd_t_func_ptr;

        auto in = (i & 1) == 0 ? v2.data() : v1.data();
	auto out= (i & 1) == 0 ? v1.data() : v2.data();
        std::transform(policy,
                      in,
		      in+outer_range,
                      out,
		      out,
                      [&func, &in] (const auto &in_, auto &out_) {
		        const auto i = &in_ - in; //deduce current index
		        return func.template operator()<SP>(i);});
      }
      return;
    }

    ~FwdEulerIters(){}
};
