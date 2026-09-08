module DyadShip

using BlockComponents
using RotationalComponents
using TranslationalComponents
using ThermalComponents
using MultibodyComponents
using ElectricalComponents
using DiscreteComponents
using HVACComponents
using IncompressibleFlowComponents
using ModelingToolkit

include("FlettnerCFDLive.jl")

"""
    flettner_live_fx(t)

Body-frame x-component of the rotor force, as set by the most recent
`FlettnerCFDLive.step!` call. `t` is the Dyad time and is intentionally
unused — taking it as an argument is what keeps MTK from constant-folding
this call away during simplification.
"""
flettner_live_fx(t::Real) = FlettnerCFDLive.Fx()
flettner_live_fy(t::Real) = FlettnerCFDLive.Fy()

@register_symbolic flettner_live_fx(t::Real)::Real
@register_symbolic flettner_live_fy(t::Real)::Real

include("../generated/module.jl")

end # module DyadShip