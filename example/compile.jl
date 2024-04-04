function GVcompileC_so(partition, datatype::DataType=Float64; c_source=joinpath(@__DIR__, "source_codeGV", "func_sigmaGV.c"),
    lib_path=joinpath(@__DIR__, "source_codeGV"), lib_name="sigmaGV", compiler::String="gcc", isnative::Bool=false)
    lib = joinpath(lib_path, "$lib_name.so")
    ### compile C source file to *.so library
    if isnative
        command = `$compiler -shared -fPIC -o $lib $c_source -O3 -march=native`
    else
        command = `$compiler -shared -fPIC -o $lib $c_source -O3`
    end
    run(command)

    ### generate the C wrapper Julia functions and save them in the same path of C library
    GV_Cwrapper(partition, datatype, lib_path=lib_path, lib_name=lib_name)
end