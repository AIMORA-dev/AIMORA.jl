module ReportArtifacts
include(joinpath(@__DIR__, "reports", "artifact_types_and_binary_plot.jl"))
include(joinpath(@__DIR__, "reports", "report_utils.jl"))
include(joinpath(@__DIR__, "reports", "runtime_output_bundle.jl"))
include(joinpath(@__DIR__, "reports", "text_output_suite.jl"))
include(joinpath(@__DIR__, "reports", "finalization_and_results.jl"))
include(joinpath(@__DIR__, "reports", "cable_constants_report.jl"))

end
