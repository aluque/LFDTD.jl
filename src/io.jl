function save(fname, mesh::Mesh, fields::Fields)
    jldsave(fname, compress=true; mesh, fields)
    @info "\u21b3 $fname"
end
