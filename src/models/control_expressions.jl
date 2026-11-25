struct ControlExpressionInstruction
    operation::Symbol
    value::Float64
    name::Symbol
end

struct ControlExpressionProgram
    output_name::Symbol
    instructions::Vector{ControlExpressionInstruction}
    input_names::Vector{Symbol}
    source_text::String
    max_stack_depth::Int
end

mutable struct ControlExpressionRuntime
    program::ControlExpressionProgram
    stack::Vector{Float64}
    random_draw_count::Int
end

function ControlExpressionRuntime(program::ControlExpressionProgram, stack::Vector{Float64})
    return ControlExpressionRuntime(program, stack,
        count(instruction -> instruction.operation == :bpa_uniform_random, program.instructions))
end

ControlExpressionRuntime(program::ControlExpressionProgram) =
    ControlExpressionRuntime(program, zeros(Float64, program.max_stack_depth))

struct _ControlExpressionToken
    kind::Symbol
    text::String
    value::Float64
end

const _CONTROL_EXPRESSION_DOTTED_OPERATORS = Dict(
    ".NOT." => :logical_not,
    ".OR." => :logical_or,
    ".NOR." => :logical_nor,
    ".AND." => :logical_and,
    ".NAND." => :logical_nand,
    ".NE." => :not_equal,
    ".EQ." => :equal,
    ".LT." => :less_than,
    ".LE." => :less_equal,
    ".GE." => :greater_equal,
    ".GT." => :greater_than,
)

const _CONTROL_EXPRESSION_WORD_OPERATORS = Dict(
    "NOT" => :logical_not,
    "OR" => :logical_or,
    "NOR" => :logical_nor,
    "AND" => :logical_and,
    "NAND" => :logical_nand,
)

const _CONTROL_EXPRESSION_FUNCTIONS = Dict(
    "SIN" => :sine,
    "COS" => :cosine,
    "TAN" => :tangent,
    "COTAN" => :cotangent,
    "SINH" => :hyperbolic_sine,
    "COSH" => :hyperbolic_cosine,
    "TANH" => :hyperbolic_tangent,
    "ASIN" => :arc_sine,
    "ACOS" => :arc_cosine,
    "ATAN" => :arc_tangent,
    "EXP" => :exponential,
    "LOG" => :natural_logarithm,
    "LOG10" => :common_logarithm,
    "SQRT" => :square_root,
    "ABS" => :absolute_value,
    "TRUNC" => :truncate,
    "MINUS" => :negate,
    "INVRS" => :reciprocal,
    "RAD" => :degrees_to_radians,
    "DEG" => :radians_to_degrees,
    "SEQ6" => :sequence_six,
    "SIGN" => :sign,
    "NOT" => :logical_not,
)

const _BPA_CONTROL_EXPRESSION_INTRINSICS = Dict(
    :reciprocal => :bpa_reciprocal,
    :exponential => :bpa_exponential,
    :cotangent => :bpa_cotangent,
    :degrees_to_radians => :bpa_degrees_to_radians,
    :radians_to_degrees => :bpa_radians_to_degrees,
    :sequence_six => :bpa_sequence_six,
)

function _control_expression_tokens(source_text::AbstractString)
    chars = collect(uppercase(strip(String(source_text))))
    tokens = _ControlExpressionToken[]
    index = 1
    while index <= length(chars)
        character = chars[index]
        if isspace(character)
            index += 1
            continue
        elseif isdigit(character) ||
               (character == '.' && index < length(chars) && isdigit(chars[index + 1]))
            first_index = index
            while index <= length(chars) && (isdigit(chars[index]) || chars[index] == '.')
                index += 1
            end
            if index <= length(chars) && chars[index] in ('E', 'D')
                index += 1
                if index <= length(chars) && chars[index] in ('+', '-')
                    index += 1
                end
                exponent_start = index
                while index <= length(chars) && isdigit(chars[index])
                    index += 1
                end
                exponent_start < index ||
                    throw(ArgumentError("control expression has an incomplete exponent"))
            end
            numeric_text = String(chars[first_index:(index - 1)])
            numeric_value = tryparse(Float64, replace(numeric_text, 'D' => 'E'))
            numeric_value === nothing &&
                throw(ArgumentError("invalid control-expression number $numeric_text"))
            isfinite(numeric_value) ||
                throw(ArgumentError("control-expression numbers must be finite"))
            push!(tokens, _ControlExpressionToken(:number, numeric_text, numeric_value))
            continue
        elseif character == '.'
            closing_offset = findnext(==('.'), chars, index + 1)
            closing_offset === nothing &&
                throw(ArgumentError("unterminated dotted control-expression operator"))
            operator_text = String(chars[index:closing_offset])
            operation = get(_CONTROL_EXPRESSION_DOTTED_OPERATORS, operator_text, nothing)
            operation === nothing &&
                throw(ArgumentError("unknown control-expression operator $operator_text"))
            push!(tokens, _ControlExpressionToken(operation, operator_text, 0.0))
            index = closing_offset + 1
            continue
        elseif isletter(character) || character == '_'
            first_index = index
            while index <= length(chars) &&
                  (isletter(chars[index]) || isdigit(chars[index]) || chars[index] == '_')
                index += 1
            end
            word = String(chars[first_index:(index - 1)])
            kind = get(_CONTROL_EXPRESSION_WORD_OPERATORS, word, :identifier)
            push!(tokens, _ControlExpressionToken(kind, word, 0.0))
            continue
        elseif character == '*'
            if index < length(chars) && chars[index + 1] == '*'
                push!(tokens, _ControlExpressionToken(:power, "**", 0.0))
                index += 2
            else
                push!(tokens, _ControlExpressionToken(:multiply, "*", 0.0))
                index += 1
            end
            continue
        end
        operation =
            character == '+' ? :add :
            character == '-' ? :subtract :
            character == '/' ? :divide :
            character == '^' ? :power :
            character == '(' ? :open :
            character == ')' ? :close :
            nothing
        operation === nothing &&
            throw(ArgumentError("unknown control-expression character $character"))
        push!(tokens, _ControlExpressionToken(operation, string(character), 0.0))
        index += 1
    end
    push!(tokens, _ControlExpressionToken(:end_of_expression, "", 0.0))
    return tokens
end

mutable struct _ControlExpressionCompiler
    tokens::Vector{_ControlExpressionToken}
    index::Int
    instructions::Vector{ControlExpressionInstruction}
    input_names::Vector{Symbol}
    bpa_factor_arithmetic::Bool
end

_control_expression_current(compiler::_ControlExpressionCompiler) =
    compiler.tokens[compiler.index]

function _control_expression_accept!(compiler::_ControlExpressionCompiler, kind::Symbol)
    _control_expression_current(compiler).kind == kind || return false
    compiler.index += 1
    return true
end

function _control_expression_emit!(
    compiler::_ControlExpressionCompiler,
    operation::Symbol;
    value::Float64=0.0,
    name::Symbol=:none,
)
    push!(compiler.instructions, ControlExpressionInstruction(operation, value, name))
    return compiler
end

function _parse_control_expression_or!(compiler::_ControlExpressionCompiler)
    _parse_control_expression_and!(compiler)
    while _control_expression_current(compiler).kind in (:logical_or, :logical_nor)
        operation = _control_expression_current(compiler).kind
        compiler.index += 1
        _parse_control_expression_and!(compiler)
        _control_expression_emit!(compiler, operation)
    end
end

function _parse_control_expression_and!(compiler::_ControlExpressionCompiler)
    _parse_control_expression_not!(compiler)
    while _control_expression_current(compiler).kind in (:logical_and, :logical_nand)
        operation = _control_expression_current(compiler).kind
        compiler.index += 1
        _parse_control_expression_not!(compiler)
        _control_expression_emit!(compiler, operation)
    end
end

function _parse_control_expression_not!(compiler::_ControlExpressionCompiler)
    if _control_expression_accept!(compiler, :logical_not)
        _parse_control_expression_not!(compiler)
        _control_expression_emit!(compiler, :logical_not)
    else
        _parse_control_expression_comparison!(compiler)
    end
end

function _parse_control_expression_comparison!(compiler::_ControlExpressionCompiler)
    _parse_control_expression_additive!(compiler)
    comparisons = (
        :equal,
        :not_equal,
        :less_than,
        :less_equal,
        :greater_than,
        :greater_equal,
    )
    if _control_expression_current(compiler).kind in comparisons
        operation = _control_expression_current(compiler).kind
        compiler.index += 1
        _parse_control_expression_additive!(compiler)
        _control_expression_emit!(compiler, operation)
    end
end

function _emit_bpa_control_expression_terms!(
    compiler::_ControlExpressionCompiler,
    first_term_end::Int,
    operation_offsets::Vector{Int},
    start_value::Float64,
    combine_operation::Symbol,
)
    # CSUP accumulates the terms from right to left and combines the first term last.
    suffix = compiler.instructions[(first_term_end + 1):end]
    resize!(compiler.instructions, first_term_end)
    _control_expression_emit!(compiler, :literal; value = start_value)
    for index in reverse(eachindex(operation_offsets))
        term_start = index == 1 ? 1 : operation_offsets[index - 1] + 1
        append!(compiler.instructions, @view suffix[term_start:operation_offsets[index]])
    end
    _control_expression_emit!(compiler, combine_operation)
    return
end

function _parse_control_expression_additive!(compiler::_ControlExpressionCompiler)
    _parse_control_expression_multiplicative!(compiler)
    _control_expression_current(compiler).kind in (:add, :subtract) || return

    first_term_end = length(compiler.instructions)
    operation_offsets = compiler.bpa_factor_arithmetic ? Int[] : nothing
    while _control_expression_current(compiler).kind in (:add, :subtract)
        operation = _control_expression_current(compiler).kind
        compiler.index += 1
        _parse_control_expression_multiplicative!(compiler)
        _control_expression_emit!(compiler, operation)
        operation_offsets === nothing ||
            push!(operation_offsets, length(compiler.instructions) - first_term_end)
    end
    operation_offsets === nothing && return
    _emit_bpa_control_expression_terms!(compiler, first_term_end, operation_offsets, 0.0, :add)
end

function _parse_control_expression_multiplicative!(compiler::_ControlExpressionCompiler)
    _parse_control_expression_unary!(compiler)
    _control_expression_current(compiler).kind in (:multiply, :divide) || return

    first_factor_end = length(compiler.instructions)
    operation_offsets = compiler.bpa_factor_arithmetic ? Int[] : nothing
    while _control_expression_current(compiler).kind in (:multiply, :divide)
        operation = _control_expression_current(compiler).kind
        compiler.index += 1
        _parse_control_expression_unary!(compiler)
        if compiler.bpa_factor_arithmetic && operation == :divide
            operation = :bpa_divide
        end
        _control_expression_emit!(compiler, operation)
        operation_offsets === nothing ||
            push!(operation_offsets, length(compiler.instructions) - first_factor_end)
    end
    operation_offsets === nothing && return
    _emit_bpa_control_expression_terms!(
        compiler, first_factor_end, operation_offsets, 1.0, :multiply,
    )
end

function _parse_control_expression_unary!(compiler::_ControlExpressionCompiler)
    if _control_expression_accept!(compiler, :add)
        _parse_control_expression_unary!(compiler)
    elseif _control_expression_accept!(compiler, :subtract)
        _parse_control_expression_unary!(compiler)
        _control_expression_emit!(compiler, :negate)
    else
        _parse_control_expression_power!(compiler)
    end
end

function _parse_control_expression_power!(compiler::_ControlExpressionCompiler)
    _parse_control_expression_primary!(compiler)
    if _control_expression_accept!(compiler, :power)
        _parse_control_expression_unary!(compiler)
        _control_expression_emit!(compiler, :power)
    end
end

function _parse_control_expression_primary!(compiler::_ControlExpressionCompiler)
    token = _control_expression_current(compiler)
    if token.kind == :number
        compiler.index += 1
        _control_expression_emit!(compiler, :literal; value = token.value)
    elseif token.kind == :identifier
        compiler.index += 1
        if _control_expression_accept!(compiler, :open)
            operation = compiler.bpa_factor_arithmetic && token.text == "RAN" ?
                :bpa_uniform_random : get(_CONTROL_EXPRESSION_FUNCTIONS, token.text, nothing)
            operation === nothing &&
                throw(ArgumentError("unknown control-expression function $(token.text)"))
            _parse_control_expression_or!(compiler)
            _control_expression_accept!(compiler, :close) ||
                throw(ArgumentError("control-expression function $(token.text) lacks a closing parenthesis"))
            if compiler.bpa_factor_arithmetic
                operation = get(_BPA_CONTROL_EXPRESSION_INTRINSICS, operation, operation)
            end
            _control_expression_emit!(compiler, operation)
        else
            name = Symbol(token.text)
            name in compiler.input_names || push!(compiler.input_names, name)
            _control_expression_emit!(compiler, :signal; name = name)
        end
    elseif _control_expression_accept!(compiler, :open)
        _parse_control_expression_or!(compiler)
        _control_expression_accept!(compiler, :close) ||
            throw(ArgumentError("control expression lacks a closing parenthesis"))
    else
        throw(ArgumentError("expected a control-expression value, got $(token.text)"))
    end
end

const _CONTROL_EXPRESSION_BINARY_OPERATIONS = (
    :add, :subtract, :multiply, :divide, :bpa_divide, :power,
    :logical_or, :logical_nor, :logical_and, :logical_nand,
    :equal, :not_equal, :less_than, :less_equal, :greater_than, :greater_equal,
)

function _control_expression_max_stack_depth(instructions)
    depth = 0
    maximum_depth = 0
    for instruction in instructions
        if instruction.operation in (:literal, :signal)
            depth += 1
            maximum_depth = max(maximum_depth, depth)
        elseif instruction.operation in _CONTROL_EXPRESSION_BINARY_OPERATIONS
            depth -= 1
        end
        depth > 0 || throw(ArgumentError("invalid compiled control-expression stack"))
    end
    depth == 1 || throw(ArgumentError("control expression did not compile to one result"))
    return maximum_depth
end

"""
    compile_control_expression(output_name, source_text; bpa_factor_arithmetic=false)

Compile a control expression. Fixed-card BPA deck intake selects legacy evaluation order
and guarded division and intrinsics. The retained Linux glibc profile uses system
`log10` and `exp`. Other builds use Julia's functions and may round cutoff or pole
arguments differently.
"""
function compile_control_expression(
    output_name::Symbol,
    source_text::AbstractString;
    bpa_factor_arithmetic::Bool = false,
)
    text = strip(String(source_text))
    isempty(text) && throw(ArgumentError("control expression must not be empty"))
    compiler = _ControlExpressionCompiler(
        _control_expression_tokens(text),
        1,
        ControlExpressionInstruction[],
        Symbol[],
        bpa_factor_arithmetic,
    )
    _parse_control_expression_or!(compiler)
    trailing = _control_expression_current(compiler)
    trailing.kind == :end_of_expression ||
        throw(ArgumentError("unexpected trailing control-expression token $(trailing.text)"))
    depth = _control_expression_max_stack_depth(compiler.instructions)
    return ControlExpressionProgram(
        output_name,
        compiler.instructions,
        compiler.input_names,
        text,
        depth,
    )
end

_control_expression_truth(value::Float64) = value > CONTROL_SYSTEM_ZERO_TOLERANCE

function _control_expression_divide(left::Float64, right::Float64)
    left == 0.0 && return 0.0
    right == 0.0 && return signbit(left) == signbit(right) ? 1.0e300 : -1.0e300
    return left / right
end

# The retained gfortran reference uses glibc for DLOG10 and DEXP.
@static if Sys.islinux() && Base.BinaryPlatforms.libc(Base.BinaryPlatforms.HostPlatform()) == "glibc"
    _bpa_control_log10(value::Float64) = ccall((:log10, "libm.so.6"), Cdouble, (Cdouble,), value)
    _bpa_control_exp(value::Float64) = ccall((:exp, "libm.so.6"), Cdouble, (Cdouble,), value)
else
    _bpa_control_log10(value::Float64) = log10(value)
    _bpa_control_exp(value::Float64) = exp(value)
end

function _bpa_control_expression_divide(left::Float64, right::Float64)
    left == 0.0 && return 0.0
    negative = (left < 0.0) != (right < 0.0)
    saturation = negative ? -CONTROL_SYSTEM_INFINITY : CONTROL_SYSTEM_INFINITY
    right == 0.0 && return saturation
    isfinite(left) && isfinite(right) || return left / right

    # CSUP truncates the decimal exponent difference before comparing it to the bound.
    exponent = trunc(_bpa_control_log10(abs(left)) - _bpa_control_log10(abs(right)))
    limit = log10(CONTROL_SYSTEM_INFINITY)
    exponent < -limit && return 0.0
    exponent > limit && return saturation
    return left / right
end

function _control_expression_unary(operation::Symbol, value::Float64)
    if operation == :negate
        return -value
    elseif operation == :logical_not
        return _control_expression_truth(value) ? 0.0 : 1.0
    elseif operation == :sine
        return sin(value)
    elseif operation == :cosine
        return cos(value)
    elseif operation == :tangent
        return tan(value)
    elseif operation == :cotangent
        return inv(tan(value))
    elseif operation == :bpa_cotangent
        sine = sin(value)
        abs(sine) * CONTROL_SYSTEM_INFINITY > 1.0 ||
            throw(DomainError(value, "BPA TACS cotangent argument is too close to a sine zero"))
        return cos(value) / sine
    elseif operation == :hyperbolic_sine
        return sinh(value)
    elseif operation == :hyperbolic_cosine
        return cosh(value)
    elseif operation == :hyperbolic_tangent
        return tanh(value)
    elseif operation == :arc_sine
        return asin(value)
    elseif operation == :arc_cosine
        return acos(value)
    elseif operation == :arc_tangent
        return atan(value)
    elseif operation == :exponential
        return exp(value)
    elseif operation == :bpa_exponential
        return value >= -87.0 ? _bpa_control_exp(value) : 0.0
    elseif operation == :natural_logarithm
        return value == 0.0 ? -1.0e300 : log(value)
    elseif operation == :common_logarithm
        return value == 0.0 ? -1.0e300 : log10(value)
    elseif operation == :square_root
        return sqrt(value)
    elseif operation == :absolute_value
        return abs(value)
    elseif operation == :truncate
        return trunc(value)
    elseif operation == :reciprocal
        return _control_expression_divide(1.0, value)
    elseif operation == :bpa_reciprocal
        return _bpa_control_expression_divide(1.0, value)
    elseif operation == :degrees_to_radians
        return value * pi / 180.0
    elseif operation == :bpa_degrees_to_radians
        # Keep the intermediate rounding and overflow of CSUP's TWOPI product.
        return value * (2.0 * pi) / 360.0
    elseif operation == :radians_to_degrees
        return value * 180.0 / pi
    elseif operation == :bpa_radians_to_degrees
        return value / (2.0 * pi) * 360.0
    elseif operation == :sequence_six
        integer_value = trunc(Int, value)
        return Float64(mod(integer_value - 1, 6) + 1)
    elseif operation == :bpa_sequence_six
        integer_value = trunc(Int32, value)
        # CSUP's signed quotient cancels in the negative branch, leaving the input plus six.
        sequence = integer_value < 0 ? integer_value + 6 : rem(integer_value, 6)
        return iszero(sequence) ? 6.0 : Float64(sequence)
    elseif operation == :sign
        return value < 0.0 ? -1.0 : 1.0
    end
    throw(ArgumentError("unknown unary control-expression operation $operation"))
end

function _control_expression_binary(operation::Symbol, left::Float64, right::Float64)
    operation == :add && return left + right
    operation == :subtract && return left - right
    operation == :multiply && return left * right
    operation == :divide && return _control_expression_divide(left, right)
    operation == :bpa_divide && return _bpa_control_expression_divide(left, right)
    operation == :power && return left ^ right
    operation == :equal && return left == right ? 1.0 : 0.0
    operation == :not_equal && return left != right ? 1.0 : 0.0
    operation == :less_than && return left < right ? 1.0 : 0.0
    operation == :less_equal && return left <= right ? 1.0 : 0.0
    operation == :greater_than && return left > right ? 1.0 : 0.0
    operation == :greater_equal && return left >= right ? 1.0 : 0.0
    left_truth = _control_expression_truth(left)
    right_truth = _control_expression_truth(right)
    operation == :logical_or && return left_truth || right_truth ? 1.0 : 0.0
    operation == :logical_nor && return left_truth || right_truth ? 0.0 : 1.0
    operation == :logical_and && return left_truth && right_truth ? 1.0 : 0.0
    operation == :logical_nand && return left_truth && right_truth ? 0.0 : 1.0
    throw(ArgumentError("unknown binary control-expression operation $operation"))
end

const _BPA_CONTROL_RANDOM_MODULUS = Int64(2_147_483_647)
const _BPA_CONTROL_RANDOM_MULTIPLIER = Int64(16_807)

function _bpa_control_random_stream(step::Int, draws_per_step::Int)
    period = _BPA_CONTROL_RANDOM_MODULUS - 1
    previous_draws = Int64(mod(widemul(step - 1, draws_per_step), period))
    return Ref(powermod(_BPA_CONTROL_RANDOM_MULTIPLIER, previous_draws,
        _BPA_CONTROL_RANDOM_MODULUS))
end

function _next_bpa_control_random!(state::Ref{Int64})
    next = mod(_BPA_CONTROL_RANDOM_MULTIPLIER * state[], _BPA_CONTROL_RANDOM_MODULUS)
    state[] = next
    # The compiled GNU RAN profile discards the low nine bits before scaling.
    return ldexp(Float64((next - 1) & ~Int64(511)), -31)
end

function evaluate_control_expression!(
    runtime::ControlExpressionRuntime,
    values::AbstractDict{Symbol,<:Real};
    bpa_random_stream::Union{Nothing,Ref{Int64}} = nothing,
)
    stack_index = 0
    for instruction in runtime.program.instructions
        operation = instruction.operation
        if operation == :literal
            stack_index += 1
            runtime.stack[stack_index] = instruction.value
        elseif operation == :signal
            haskey(values, instruction.name) ||
                throw(ArgumentError("missing control-expression signal $(instruction.name)"))
            stack_index += 1
            runtime.stack[stack_index] = Float64(values[instruction.name])
        elseif operation in _CONTROL_EXPRESSION_BINARY_OPERATIONS
            stack_index >= 2 || throw(ArgumentError("control-expression stack underflow"))
            right = runtime.stack[stack_index]
            left = runtime.stack[stack_index - 1]
            stack_index -= 1
            runtime.stack[stack_index] = _control_expression_binary(operation, left, right)
        else
            stack_index >= 1 || throw(ArgumentError("control-expression stack underflow"))
            if operation == :bpa_uniform_random
                bpa_random_stream === nothing && throw(ArgumentError(
                    "RAN requires a fixed-card BPA timestep random stream",
                ))
                runtime.stack[stack_index] = _next_bpa_control_random!(bpa_random_stream)
            else
                runtime.stack[stack_index] =
                    _control_expression_unary(operation, runtime.stack[stack_index])
            end
        end
    end
    stack_index == 1 || throw(ArgumentError("control expression did not produce one result"))
    result = runtime.stack[1]
    isfinite(result) ||
        throw(ArgumentError("control expression $(runtime.program.output_name) produced a nonfinite result"))
    return result
end
