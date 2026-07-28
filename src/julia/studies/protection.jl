module ProtectionStudy

using ..StudyCore
using ..StudyInputProfiles

export PROTECTION_STUDY, protection_input_profile, run_protection

const PROTECTION_STUDY = StudyDescriptor(
    :protection,
    "Protection coordination",
    :protection,
    :planned,
    "src/julia/studies/protection.jl",
)

protection_input_profile() = input_profile(:protection)

function run_protection(data; kwargs...)
    validate_study_inputs(protection_input_profile(), data)
    study_not_implemented(PROTECTION_STUDY)
end

end
