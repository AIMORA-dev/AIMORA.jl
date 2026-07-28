module ShortCircuitStudy

using ..StudyCore

export SHORT_CIRCUIT_STUDY, run_short_circuit

const SHORT_CIRCUIT_STUDY = StudyDescriptor(
    :short_circuit,
    "Short circuit",
    :static,
    :planned,
    "src/studies/short_circuit.jl",
)

run_short_circuit(args...; kwargs...) = study_not_implemented(SHORT_CIRCUIT_STUDY)

end
