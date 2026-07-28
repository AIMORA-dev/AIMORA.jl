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

function _parse_control_expression_additive!(compiler::_ControlExpressionCompiler)
    _parse_control_expression_multiplicative!(compiler)
    while _control_expression_current(compiler).kind in (:add, :subtract)
        operation = _control_expression_current(compiler).kind
        compiler.index += 1
        _parse_control_expression_multiplicative!(compiler)
        _control_expression_emit!(compiler, operation)
    end
end

function _parse_control_expression_multiplicative!(compiler::_ControlExpressionCompiler)
    _parse_control_expression_unary!(compiler)
    while _control_expression_current(compiler).kind in (:multiply, :divide)
        operation = _control_expression_current(compiler).kind
        compiler.index += 1
        _parse_control_expression_unary!(compiler)
        _control_expression_emit!(compiler, operation)
    end
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
            operation = get(_CONTROL_EXPRESSION_FUNCTIONS, token.text, nothing)
            operation === nothing &&
                throw(ArgumentError("unknown control-expression function $(token.text)"))
            _parse_control_expression_or!(compiler)
            _control_expression_accept!(compiler, :close) ||
                throw(ArgumentError("control-expression function $(token.text) lacks a closing parenthesis"))
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

function _control_expression_max_stack_depth(instructions)
    depth = 0
    maximum_depth = 0
    binary_operations = (
        :add, :subtract, :multiply, :divide, :power,
        :logical_or, :logical_nor, :logical_and, :logical_nand,
        :equal, :not_equal, :less_than, :less_equal, :greater_than, :greater_equal,
    )
    for instruction in instructions
        if instruction.operation in (:literal, :signal)
            depth += 1
            maximum_depth = max(maximum_depth, depth)
        elseif instruction.operation in binary_operations
            depth -= 1
        end
        depth > 0 || throw(ArgumentError("invalid compiled control-expression stack"))
    end
    depth == 1 || throw(ArgumentError("control expression did not compile to one result"))
    return maximum_depth
end

function compile_control_expression(output_name::Symbol, source_text::AbstractString)
    text = strip(String(source_text))
    isempty(text) && throw(ArgumentError("control expression must not be empty"))
    compiler = _ControlExpressionCompiler(
        _control_expression_tokens(text),
        1,
        ControlExpressionInstruction[],
        Symbol[],
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
    elseif operation == :degrees_to_radians
        return value * pi / 180.0
    elseif operation == :radians_to_degrees
        return value * 180.0 / pi
    elseif operation == :sequence_six
        integer_value = trunc(Int, value)
        return Float64(mod(integer_value - 1, 6) + 1)
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

function evaluate_control_expression!(
    runtime::ControlExpressionRuntime,
    values::AbstractDict{Symbol,<:Real},
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
        elseif operation in (
            :add, :subtract, :multiply, :divide, :power,
            :logical_or, :logical_nor, :logical_and, :logical_nand,
            :equal, :not_equal, :less_than, :less_equal, :greater_than, :greater_equal,
        )
            stack_index >= 2 || throw(ArgumentError("control-expression stack underflow"))
            right = runtime.stack[stack_index]
            left = runtime.stack[stack_index - 1]
            stack_index -= 1
            runtime.stack[stack_index] = _control_expression_binary(operation, left, right)
        else
            stack_index >= 1 || throw(ArgumentError("control-expression stack underflow"))
            runtime.stack[stack_index] =
                _control_expression_unary(operation, runtime.stack[stack_index])
        end
    end
    stack_index == 1 || throw(ArgumentError("control expression did not produce one result"))
    result = runtime.stack[1]
    isfinite(result) ||
        throw(ArgumentError("control expression $(runtime.program.output_name) produced a nonfinite result"))
    return result
end
