# Consider renaming FiberDeseign to SIFiber
export FiberDesign

struct FiberDesign
    r::VecFloat  # r[i]: outer radius of Region i (1≤i≤N); core is Region 1
    εmdl::Vector{SellmeierModel}  # relative permittivities of Regions 1, ..., N, plus Region N+1 outside r = r[N]
    μmdl::Vector{SellmeierModel}  # relative permeabilities of Regions 1, ..., N, plus Region N+1 outside r = r[N]

    function FiberDesign(r::AbsVecReal, εmdl::AbsVec{<:SellmeierModel}, μmdl::AbsVec{<:SellmeierModel}=SellmeierModel.(fill(1.0,length(r)+1)))
        @argcheck length(εmdl) == length(μmdl) == length(r)+1

        return new(r,εmdl,μmdl)
    end
end
FiberDesign(r::AbsVecReal, ε::AbsVecReal; μ::AbsVecReal=fill(1.0,length(r)+1)) =
    FiberDesign(r, SellmeierModel.(ε), SellmeierModel.(μ))
# FiberDesign(r::AbsVecReal, εmdl::AbsVec{<:SellmeierModel}; 
#             μmdl::AbsVec{<:SellmeierModel}=SellmeierModel.(fill(1.0,length(r)+1))) =
#     FiberDesign(r, μmdl, εmdl)

Base.length(fd::FiberDesign) = length(fd.r)
