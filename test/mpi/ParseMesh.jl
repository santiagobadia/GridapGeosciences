module ParseMesh

  using Printf
  using LinearMaps
  using SparseMatricesCSR
  using PartitionedArrays
  using Test
  using FillArrays
  using Gridap
  using GridapPETSc
  using GridapGeosciences
  using GridapDistributed
  using GridapSolvers
  using GridapP4est

  num_refinements = 0

  function write_connectivity(filename,cell_cells,cell_edges,cell_verts,edge_verts,coords)

    io_buf_cv = IOBuffer()
    io_buf_ce = IOBuffer()
    io_buf_cc = IOBuffer()
    io_buf_ev = IOBuffer()
    io_buf_nx = IOBuffer()
    io_buf_ny = IOBuffer()

    for c in 1:length(cell_verts)
      @printf(io_buf_cv, "  %d, %d, %d, %d\n", cell_verts[c][1], cell_verts[c][2], cell_verts[c][3], cell_verts[c][4])
      @printf(io_buf_ce, "  %d, %d, %d, %d\n", cell_edges[c][1], cell_edges[c][2], cell_edges[c][3], cell_edges[c][4])
      @printf(io_buf_cc, "  %d, %d, %d, %d\n", cell_cells[c][1], cell_cells[c][2], cell_cells[c][3], cell_cells[c][4])
    end
    for e in 1:length(edge_verts)
      @printf(io_buf_ev, "  %d, %d\n", edge_verts[e][1], edge_verts[e][2])
    end
    for v in 1:length(coords)
      @printf(io_buf_nx, "%f\n", coords[v][1])
      @printf(io_buf_ny, "%f\n", coords[v][2])
    end
    open(filename,"w") do file
      write(file," dynamics_face_nodes =\n")
      write(file,take!(io_buf_cv))
      write(file," dynamics_face_edges =\n")
      write(file,take!(io_buf_ce))
      write(file," dynamics_face_links =\n")
      write(file,take!(io_buf_cc))
      write(file," dynamics_edge_nodes =\n")
      write(file,take!(io_buf_ev))
      write(file," dynamics_node_x =\n")
      write(file,take!(io_buf_nx))
      write(file," dynamics_node_y =\n")
      write(file,take!(io_buf_ny))
    end

    # TODO: write the node geometry
    # TODO: write the cell links between MG meshes
  end

  function main(distribute,parts)
    ranks = distribute(LinearIndices((prod(parts),)))

    # Change directory to the location of the script, where the mesh data files are located 
    cd(@__DIR__)

    coarse_model, cell_panels, coarse_cell_wise_vertex_coordinates = parse_cubed_sphere_coarse_model("williamson-5-C12/connectivity-gridapgeo.txt",
                                                                                                     "williamson-5-C12/geometry-gridapgeo.txt")

    model = CubedSphereDiscreteModel(ranks,
                                     coarse_model,
                                     coarse_cell_wise_vertex_coordinates,
                                     cell_panels,
                                     num_refinements;
                                     radius=6371220.0,
                                     adaptive=false,
				     order=1)

    topo = Gridap.Geometry.get_grid_topology(coarse_model)
    #topo = Gridap.Geometry.get_grid_topology(model.cubed_sphere_linear_model)
    cell_cells = Gridap.Geometry.get_faces(topo,2,1)
    cell_edges = Gridap.Geometry.get_faces(topo,2,1)
    cell_verts = Gridap.Geometry.get_faces(topo,2,0)
    edge_verts = Gridap.Geometry.get_faces(topo,1,0)

    coords = Gridap.Geometry.get_node_coordinates(coarse_model)
    #coords = Gridap.Geometry.get_node_coordinates(model)
    #coords = Gridap.Geometry.get_node_coordinates(get_grid(model.cubed_sphere_linear_model))

    write_connectivity("parsed_geometry.txt",cell_cells,cell_edges,cell_verts,edge_verts,coords)
  end

  with_mpi() do distribute 
    main(distribute,1)
  end
end # module
