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

  function adapt_model(ranks,model)
    cell_partition=get_cell_gids(model.octree_model.dmodel)
    ref_coarse_flags=map(ranks,partition(cell_partition)) do rank,indices
        flags=zeros(Cint,length(indices))
        flags.=refine_flag
    end
    GridapP4est.adapt(model,ref_coarse_flags)
  end

  function write_connectivity(filename,topo,model)
    cell_cells = Gridap.Geometry.get_faces(topo,2,2)
    cell_edges = Gridap.Geometry.get_faces(topo,2,1)
    cell_verts = Gridap.Geometry.get_faces(topo,2,0)
    edge_verts = Gridap.Geometry.get_faces(topo,1,0)
    edge_cells = Gridap.Geometry.get_faces(topo,1,2)
    coords     = Gridap.Geometry.get_node_coordinates(get_grid(model))

    io_buf_cv = IOBuffer()
    io_buf_ce = IOBuffer()
    io_buf_cc = IOBuffer()
    io_buf_ev = IOBuffer()
    io_buf_nx = IOBuffer()
    io_buf_ny = IOBuffer()

    _cell_cells = Vector{Int32}(undef, 4)

    for c in 1:length(cell_verts)
      @printf(io_buf_cv, "  %d, %d, %d, %d\n", cell_verts[c][1], cell_verts[c][2], cell_verts[c][3], cell_verts[c][4])
      @printf(io_buf_ce, "  %d, %d, %d, %d\n", cell_edges[c][1], cell_edges[c][2], cell_edges[c][3], cell_edges[c][4])
      i_cell = cell_cells[c][1]
      for d in 1:4
	i_edge = cell_edges[c][d]
	if i_cell == edge_cells[i_edge][1]
	  _cell_cells[d] = edge_cells[i_edge][2]
	else
	  _cell_cells[d] = edge_cells[i_edge][1]
	end
      end
      @printf(io_buf_cc, "  %d, %d, %d, %d\n", _cell_cells[1], _cell_cells[2], _cell_cells[3], _cell_cells[4])
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

    # TODO: write the cell links between MG meshes
  end

  function main(distribute,parts)
    ranks = distribute(LinearIndices((prod(parts),)))

    # Change directory to the location of the script, where the mesh data files are located 
    cd(@__DIR__)

    coarse_model, cell_panels, coarse_vertex_coords = parse_cubed_sphere_coarse_model("C12-regular/connectivity-gridapgeo.txt",
                                                                                      "C12-regular/geometry-gridapgeo.txt")

    num_refinements = 3

    model = CubedSphereDiscreteModel(ranks,
                                     coarse_model,
                                     coarse_vertex_coords,
                                     cell_panels,
                                     0;
                                     radius=6371220.0,
                                     adaptive=true,
                                     order=1)

    for refinement_level in 1:num_refinements
      topo = Gridap.Geometry.get_grid_topology(model)
      write_connectivity("parsed_geometry_$(refinement_level).txt",local_views(topo).item_ref[],local_views(model).item_ref[])
      model,_ = adapt_model(ranks,model)
    end
    topo = Gridap.Geometry.get_grid_topology(model)
    refinement_level = num_refinements + 1
    write_connectivity("parsed_geometry_$(refinement_level).txt",local_views(topo).item_ref[],local_views(model).item_ref[])
  end

  with_mpi() do distribute 
    main(distribute,1)
  end
end # module
