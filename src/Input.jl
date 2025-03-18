
"""
    Input{T}
Structure in which all the input information is contained. It is used repetively during the whole calculation.
`T` is a parametric Julia type.
"""
struct Input{T <: Number, S <: Number}
    model::Dict{String,AbstractString}
    physical_params::Dict{String,Union{Nothing,T}}
    numerical_params::Dict{String,Union{Nothing,S,T}}
    calculation_details::Dict{String, Union{Nothing,Bool}}
end

# Overloading of the print function for union dictionary
for S in [Number, AbstractString, Bool]
    function Base.show(io::IO, dict::Dict{String,Union{Nothing, T}}) where T <: S
        if isempty(dict)
            println(io, "$(typeof(dict))()")
            return
        end
        # Convert keys and values to strings and collect as tuples.
        items = [(string(k), string(v)) for (k, v) in dict]
        # Sort items by the key's string representation.
        sorted_items = sort(items, by = x -> x[1])
        
        # Determine the maximum widths of the keys and values.
        key_width   = maximum(length.(getindex.(sorted_items, 1)))
        value_width = maximum(length.(getindex.(sorted_items, 2)))
        separator   = " : "
        inner_width = key_width + length(separator) + value_width
    
        # Build the box borders using Unicode box-drawing characters.
        top_border    = "┌" * repeat("─", inner_width) * "┐"
        bottom_border = "└" * repeat("─", inner_width) * "┘"
    
        # Print the box.
        println(io, top_border)
        for (k, v) in sorted_items
            line = "│" * rpad(k, key_width) * separator * rpad(v, value_width) * "│"
            println(io, line)
        end
        println(io, bottom_border)
    end
end

function Base.show(io::IO, dict::Dict{String,AbstractString})
    if isempty(dict)
        println(io, "$(typeof(dict))()")
        return
    end
    # Convert keys and values to strings and collect as tuples.
    items = [(string(k), string(v)) for (k, v) in dict]
    # Sort items by the key's string representation.
    sorted_items = sort(items, by = x -> x[1])
    
    # Determine the maximum widths of the keys and values.
    key_width   = maximum(length.(getindex.(sorted_items, 1)))
    value_width = maximum(length.(getindex.(sorted_items, 2)))
    separator   = " : "
    inner_width = key_width + length(separator) + value_width

    # Build the box borders using Unicode box-drawing characters.
    top_border    = "┌" * repeat("─", inner_width) * "┐"
    bottom_border = "└" * repeat("─", inner_width) * "┘"

    # Print the box.
    println(io, top_border)
    for (k, v) in sorted_items
        line = "│" * rpad(k, key_width) * separator * rpad(v, value_width) * "│"
        println(io, line)
    end
    println(io, bottom_border)
end

# Overloading of the print function for the struct
function Base.show(io::IO, input::Input)
    println("Model:")
    println(input.model)
    println("Physical Parameters:")
    println(input.physical_params)
    println("Numerical Parameters:")
    println(input.numerical_params)
    println("Calculation Details:")
    println(input.calculation_details)
end



"""
    GoodParser(str)
Parse any string with format _number unit_ into a float doing the proper unit change.
Everything will be converter to automatically to atomic units, i.e. ħ = e = m_e = 1.
"""
function GoodParser(str::AbstractString)
    energy_units = ["eV", "meV", "K", "au"]
    time_units = ["fs", "au"]
    laser_units = ["Vm", "Wcm2", "nm", "au"]
    energy_functions = [x * "_to_au" for x in energy_units]
    time_functions = [x * "_to_au" for x in time_units]
    laser_functions = [x * "_to_au" for x in laser_units]

    str_array = split(str, " ")
    if str_array[2] in energy_units
        func = energy_functions[findfirst(x -> x == str_array[2], energy_units)]
        # Call the function to do the automatic unit conversion.
        # Using this method we avoid coding an if statement for each different unit.
        return getfield(Hydrogen, Symbol(func))(parse(Float64, str_array[1]))
    elseif str_array[2] in time_units
        func = time_functions[findfirst(x -> x == str_array[2], time_units)]
        return getfield(Hydrogen, Symbol(func))(parse(Float64, str_array[1]))

    elseif str_array[2] in laser_units
        func = laser_functions[findfirst(x -> x == str_array[2], laser_units)]
        return getfield(Hydrogen, Symbol(func))(parse(Float64, str_array[1]))
    end
    # Note: "au" is repeated into all different arrays. This won't affect the parsing
    # because it will enter the first if statement and do the conversion au_to_au() which
    # does nothing but return the same number. Therefore, it doesn't matter if the code is doing
    # a conversion of a time in atomic units as it was an energy parameter because either way it won't
    # do anything to the number.
end


function IniParserUnits(ini_file, section::String, field::String)
    # Reads the given field of the given section of the INI file knowing it has physical unit.
    # Therefore, it calls GoodParser after reading the field to do the proper unit converison.
    try
        str = IniFile.get(ini_file, section, field)
        return GoodParser(str)
    catch 
        return nothing
    end
end

function IniParserType(ini_file, section::String, field::String, type)
    # Reads the given field of the given section of the INI file knowing its type.
    # This is useful for bool variables and numerical parameters without units.
    try
        param = tryparse(type, IniFile.get(ini_file, section, field))
        return param
    catch
        return nothing
    end
end

function PhysParams(ini_file)

    #Read the section [Physical Parameters] from an INI file and get its fields as a dictionary.
    #Currently, it only works in atomic units. 
    #To not include some variable, just type not included in the INI file.

    #Mandatory params.
    ω = IniParserUnits(ini_file, "Physical Parameters", "Laser frequency")
    E0 = IniParserUnits(ini_file, "Physical Parameters", "Laser amplitude")
    τ = IniParserUnits(ini_file, "Physical Parameters", "Laser Tau")
    if any(x -> x === nothing, [ω, E0, τ])
        error("Mandatory para  meters must be specified in input file!!")
    end

    #Optional params.
    Gap = IniParserUnits(ini_file, "Physical Parameters", "Gap")
    E1_energy = IniParserUnits(ini_file, "Physical Parameters", "Energy E1")
    if any(x -> x === nothing, [Gap, E1_energy])
        @info "Some optional variables are missing. This may cause an error afterwards."
    end
    
    return Dict("ω" => ω, "E0" => E0, "τ" => τ, "Gap" => Gap, "E1_energy" => E1_energy) 
end

function BooleanVars(ini_file)
    #Reads the section [Calculation Details] from the INI file. 
    projected_bool = IniParserType(ini_file, "Calculation Details", "Rydberg calculation", Bool)
    timeprop_bool = IniParserType(ini_file, "Calculation Details", "Time propagation", Bool)
    optim_bool_no_ryd = IniParserType(ini_file, "Calculation Details", "Optimization No Rydbergs", Bool)
    optim_bool_ryd = IniParserType(ini_file, "Calculation Details", "Optimization Rydbergs", Bool)

    if any(x -> x === nothing, [projected_bool, timeprop_bool, optim_bool_no_ryd, optim_bool_ryd])
        error("Some specifications on which calculations you want to make are missing!")
    end
    return Dict("projected_bool" => projected_bool, "timeprop_bool" => timeprop_bool, 
                    "optim_bool_no_ryd" => optim_bool_no_ryd, "optim_bool_ryd" => optim_bool_ryd,
                    )
end

function NumParams(ini_file)
    #Read the section [Numerical Parameters] from an INI file and get its fields as a dictionary.   
    dx = IniParserType(ini_file, "Numerical Parameters", "Box spacing", Float64)
    xmax = IniParserType(ini_file, "Numerical Parameters", "Box length", Int)
    dt = IniParserType(ini_file, "Numerical Parameters", "Time spacing", Float64)

    a_no_rydbergs = IniParserType(ini_file, "Numerical Parameters", "a No Rydbergs", Float64)
    η_no_rydbergs = IniParserType(ini_file, "Numerical Parameters", "eta No Rydbergs", Float64)
    a_rydbergs = IniParserType(ini_file, "Numerical Parameters", "a Rydbergs", Float64)
    η_rydbergs = IniParserType(ini_file, "Numerical Parameters", "eta Rydbergs", Float64)

    if any(x -> x === nothing, [dx, xmax, dt])
        @info "Some optional numerical variables are missing. This may cause an error afterwards."
    end
    return Dict("dx" => dx, "xmax" => xmax,"time_spacing" => dt, "a_no_rydbergs" => a_no_rydbergs,
                     "η_no_rydbergs" => η_no_rydbergs, "a_rydbergs" => a_rydbergs, "η_rydbergs" => η_rydbergs)
end

@inline function ModelInfo(ini_file)
    #Read the section [Model] from an INI file. Furthermore, it loads information specific to that model.
    laser_type  = IniFile.get(ini_file, "Model", "Laser type")
    return Dict("laser type" => laser_type)
end


"""
    ReadInput("filename")

Reads the input parameters from an INI file divided in various sections. These are:
    · [Model]                --> Loads information specific to the model being used.
    · [Physical Parameters]  --> Gets the fields of this section as a dictionary.
    · [Numerical Parameters] --> Gets the fields of this section as a dictionary.
    · [Calculation Details]  --> Specific details of the calculation we want to do.
"""
function ReadInput(filename::String)
    ini_file = IniFile.read(IniFile.Inifile(),filename)
    #Obtain the model
    model = ModelInfo(ini_file)
    #Load everything else
    phys_param = PhysParams(ini_file)
    num_param = NumParams(ini_file)
    calculation_details = BooleanVars(ini_file)
    #Create the struct
    input = Input{Float64,Int}(model, phys_param, num_param, calculation_details)
    println("The given input is:              (Remainder: all physical quantities are given in atomic units)")
    print(input) 
    return input
end
