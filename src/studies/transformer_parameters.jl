module TransformerParameterStudy

using ..TransformerParameterInput:
    OVER41TransformerParameterDeck,
    parse_over41_transformer_parameter_file,
    parse_over41_transformer_parameter_lines
using ..TransformerParameters:
    AbstractTransformerParameterCase,
    MultiphaseTransformerParameterCase,
    MultiphaseTransformerParameterResult,
    SaturableTransformerParameterCase,
    SaturableTransformerParameterResult,
    TransformerShortCircuitCase,
    TransformerShortCircuitResult,
    multiphase_transformer_parameters,
    saturable_transformer_parameters,
    transformer_short_circuit_parameters

export TransformerParameterStudyResult,
       run_transformer_parameter_study,
       run_transformer_parameter_study_file,
       run_transformer_parameter_study_lines

const TransformerParameterCaseResult =
    Union{
        TransformerShortCircuitResult,
        SaturableTransformerParameterResult,
        MultiphaseTransformerParameterResult,
    }

struct TransformerParameterStudyResult
    source::String
    inputs::Vector{AbstractTransformerParameterCase}
    case_results::Vector{TransformerParameterCaseResult}
    generated_branch_count::Int
    physical_checks_passed::Bool
end

_run_case(case::TransformerShortCircuitCase) =
    transformer_short_circuit_parameters(case)
_run_case(case::SaturableTransformerParameterCase) =
    saturable_transformer_parameters(case)
_run_case(case::MultiphaseTransformerParameterCase) =
    multiphase_transformer_parameters(case)

function run_transformer_parameter_study(
    deck::OVER41TransformerParameterDeck,
)
    results = TransformerParameterCaseResult[_run_case(case) for case in deck.cases]
    branch_count = sum(
        length(getproperty(result, :generated_branches))
        for result in results
    )
    checks =
        !isempty(results) &&
        all(result -> getproperty(result, :physical_checks_passed), results) &&
        branch_count == sum(
            case isa TransformerShortCircuitCase ?
                length(case.winding_voltages_kv) :
                case isa SaturableTransformerParameterCase ?
                length(case.windings) :
                case.phase_count * length(case.windings)
            for case in deck.cases
        )
    return TransformerParameterStudyResult(
        deck.source,
        copy(deck.cases),
        results,
        branch_count,
        checks,
    )
end

run_transformer_parameter_study_lines(lines; source::AbstractString = "deck") =
    run_transformer_parameter_study(
        parse_over41_transformer_parameter_lines(lines; source),
    )

run_transformer_parameter_study_file(path::AbstractString) =
    run_transformer_parameter_study(
        parse_over41_transformer_parameter_file(path),
    )

end
