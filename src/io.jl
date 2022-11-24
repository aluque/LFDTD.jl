function save(fname, mesh::Mesh, fields::Fields)
    jldsave(fname, compress=true; mesh, fields)
end

function saveobs(fname, observer)
    jldsave(fname, compress=true; observer)
end

function saveprobes(fname, locations, time, values)
    jldsave(fname, compress=true; locations, time, values)
end

function meshvar(fields, mesh, var)
    vars = Dict(:en => (f -> en(f), true, true, "Reduced field (Td)"),
                :eabs => (f -> f.eabs, true, true, "Field magnitude (V/m)"),
                :ne => (f -> f.ne, true, true, "Electron density \$\\mathdefault{m^{-3}}\$"))

    func, rinface, zinface, label = vars[var]

    r = rinface ? rface(mesh) : rcenter(mesh)
    z = zinface ? zface(mesh) : zcenter(mesh)
    v = valid(func(fields), mesh, rinface, zinface)
    
    return (r ./ co.kilo, z ./ co.kilo, v')
end
