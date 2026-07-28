module PowerFlowStudy

using ..StudyCore

export POWER_FLOW_STUDY, run_power_flow

const POWER_FLOW_STUDY = StudyDescriptor(
    :power_flow,
    "Power flow",
    :static,
    :planned,
    "src/julia/studies/power_flow.jl",
)

run_power_flow(args...; kwargs...) = study_not_implemented(POWER_FLOW_STUDY)

end
