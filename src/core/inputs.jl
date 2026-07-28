module StudyInputs

using ..ValidationCore

export InputSpec,
       StudyInputProfile,
       required_keys,
       optional_keys,
       missing_required_inputs,
       validate_study_inputs_issues,
       validate_study_inputs

struct InputSpec
    key::Symbol
    label::String
    unit::Union{Nothing,String}
    description::String
end

function InputSpec(
    key::Symbol;
    label::AbstractString = String(key),
    unit::Union{Nothing,AbstractString} = nothing,
    description::AbstractString = "",
)
    return InputSpec(
        key,
        String(label),
        unit === nothing ? nothing : String(unit),
        String(description),
    )
end

struct StudyInputProfile
    study::Symbol
    required::Vector{InputSpec}
    optional::Vector{InputSpec}
    notes::String
end

function StudyInputProfile(
    study::Symbol,
    required::Vector{InputSpec},
    optional::Vector{InputSpec} = InputSpec[];
    notes::AbstractString = "",
)
    return StudyInputProfile(study, required, optional, String(notes))
end

required_keys(profile::StudyInputProfile) = [spec.key for spec in profile.required]
optional_keys(profile::StudyInputProfile) = [spec.key for spec in profile.optional]

function has_input(data, key::Symbol)
    if data isa AbstractDict
        return haskey(data, key) || haskey(data, String(key))
    elseif data isa NamedTuple
        return key in keys(data)
    else
        return key in propertynames(data)
    end
end

function missing_required_inputs(profile::StudyInputProfile, data)
    return [spec for spec in profile.required if !has_input(data, spec.key)]
end

function validate_study_inputs_issues(profile::StudyInputProfile, data)
    result = validation_result(source = "study $(profile.study) inputs")
    missing = missing_required_inputs(profile, data)
    for spec in missing
        add_issue!(
            result,
            missing_data(
                "study $(profile.study)",
                "Missing required input $(spec.key).";
                context = Dict(:study => profile.study, :input => spec.key),
            ),
        )
    end
    return result
end

function validate_study_inputs(profile::StudyInputProfile, data)
    assert_valid!(validate_study_inputs_issues(profile, data))
    return true
end

end
