use bevy_mesh::Mesh as BevyMesh;
use bevy_mesh::Indices;
use crate::core::mesh::MeshBackend;
use std::collections::HashMap;

/// A backend for the Mesh trait that wraps a Bevy Mesh.
/// 
/// It treats groups of vertices sharing the same geometric position as a single topological vertex,
/// allowing it to handle UV-seams and other redundant vertex data correctly as a manifold.
pub struct BevyBackend {
    pub mesh: BevyMesh,
    pub twin: Vec<Option<usize>>,
    pub canonical_vertex_halfedge: Vec<Option<usize>>,
    /// Maps a raw Bevy vertex index to its canonical topological index.
    pub raw_to_canonical: Vec<usize>,
    /// Maps a canonical topological index to its first encountered raw Bevy vertex index.
    pub canonical_to_raw: Vec<usize>,
    pub edge_halfedge: Vec<Option<usize>>,
    pub face_halfedge: Vec<Option<usize>>,
    pub boundary_halfedge: Vec<Option<usize>>,
    pub halfedge_edge: Vec<usize>,
    pub num_edges: usize,
    pub num_canonical_vertices: usize,
}

impl BevyBackend {
    pub fn new(mesh: BevyMesh) -> Self {
        let num_vertices_raw = mesh.count_vertices();
        let indices = mesh.indices().expect("Mesh must have indices");
        let num_indices = match indices {
            Indices::U16(v) => v.len(),
            Indices::U32(v) => v.len(),
        };
        let num_faces = num_indices / 3;
        let num_halfedges = num_indices;

        // Map vertex indices to canonical indices based on position to handle UV seams
        let mut raw_to_canonical = Vec::with_capacity(num_vertices_raw);
        let mut canonical_to_raw = Vec::new();
        if let Some(bevy_mesh::VertexAttributeValues::Float32x3(pos_slice)) = mesh.attribute(BevyMesh::ATTRIBUTE_POSITION) {
            let mut pos_to_idx = HashMap::new();
            for (v_idx, pos) in pos_slice.iter().enumerate() {
                // Use rounded coordinates to handle floating point noise at seams
                let key = ( (pos[0] * 1e4) as i32, (pos[1] * 1e4) as i32, (pos[2] * 1e4) as i32 );
                if let Some(&canonical_idx) = pos_to_idx.get(&key) {
                    raw_to_canonical.push(canonical_idx);
                } else {
                    let new_idx = canonical_to_raw.len();
                    pos_to_idx.insert(key, new_idx);
                    raw_to_canonical.push(new_idx);
                    canonical_to_raw.push(v_idx);
                }
            }
        } else {
            // Fallback: no positions, just use identity mapping
            raw_to_canonical = (0..num_vertices_raw).collect();
            canonical_to_raw = (0..num_vertices_raw).collect();
        }
        let num_canonical_vertices = canonical_to_raw.len();

        let mut twin = vec![None; num_halfedges];
        let mut canonical_vertex_halfedge = vec![None; num_canonical_vertices];
        let mut face_halfedge = vec![None; num_faces];
        
        let mut edge_map = HashMap::new();
        let mut next_edge_idx = 0;
        let mut halfedge_edge_indices = vec![0; num_halfedges];

        for f_idx in 0..num_faces {
            face_halfedge[f_idx] = Some(f_idx * 3);
            for i in 0..3 {
                let h_idx = f_idx * 3 + i;
                let v_curr_raw = match indices {
                    Indices::U16(v) => v[h_idx] as usize,
                    Indices::U32(v) => v[h_idx] as usize,
                };
                let v_next_raw = match indices {
                    Indices::U16(v) => v[f_idx * 3 + (i + 1) % 3] as usize,
                    Indices::U32(v) => v[f_idx * 3 + (i + 1) % 3] as usize,
                };

                let v_curr = raw_to_canonical[v_curr_raw];
                let v_next = raw_to_canonical[v_next_raw];

                canonical_vertex_halfedge[v_curr] = Some(h_idx);

                let key = if v_curr < v_next { (v_curr, v_next) } else { (v_next, v_curr) };
                if let Some(&e_idx) = edge_map.get(&key) {
                    halfedge_edge_indices[h_idx] = e_idx;
                } else {
                    edge_map.insert(key, next_edge_idx);
                    halfedge_edge_indices[h_idx] = next_edge_idx;
                    next_edge_idx += 1;
                }
            }
        }

        // Second pass for twins (using canonical vertices)
        let mut edge_to_halfedges = HashMap::new();
        for h_idx in 0..num_halfedges {
            let v_curr_raw = match indices {
                Indices::U16(v) => v[h_idx] as usize,
                Indices::U32(v) => v[h_idx] as usize,
            };
            let h_next = if h_idx % 3 == 2 { h_idx - 2 } else { h_idx + 1 };
            let v_next_raw = match indices {
                Indices::U16(v) => v[h_next] as usize,
                Indices::U32(v) => v[h_next] as usize,
            };
            
            let v_curr = raw_to_canonical[v_curr_raw];
            let v_next = raw_to_canonical[v_next_raw];

            let key = (v_curr, v_next);
            edge_to_halfedges.insert(key, h_idx);
        }

        for h_idx in 0..num_halfedges {
            let v_curr_raw = match indices {
                Indices::U16(v) => v[h_idx] as usize,
                Indices::U32(v) => v[h_idx] as usize,
            };
            let h_next = if h_idx % 3 == 2 { h_idx - 2 } else { h_idx + 1 };
            let v_next_raw = match indices {
                Indices::U16(v) => v[h_next] as usize,
                Indices::U32(v) => v[h_next] as usize,
            };

            let v_curr = raw_to_canonical[v_curr_raw];
            let v_next = raw_to_canonical[v_next_raw];

            if let Some(&twin_idx) = edge_to_halfedges.get(&(v_next, v_curr)) {
                twin[h_idx] = Some(twin_idx);
            }
        }

        // Edge halfedge
        let mut edge_halfedge_indices = vec![None; next_edge_idx];
        for h_idx in 0..num_halfedges {
            let e_idx = halfedge_edge_indices[h_idx];
            if edge_halfedge_indices[e_idx].is_none() {
                edge_halfedge_indices[e_idx] = Some(h_idx);
            }
        }

        Self {
            mesh,
            twin,
            canonical_vertex_halfedge,
            raw_to_canonical,
            canonical_to_raw,
            edge_halfedge: edge_halfedge_indices,
            face_halfedge,
            boundary_halfedge: Vec::new(),
            halfedge_edge: halfedge_edge_indices,
            num_edges: next_edge_idx,
            num_canonical_vertices,
        }
    }
}

impl MeshBackend for BevyBackend {
    fn num_vertices(&self) -> usize { self.num_canonical_vertices }
    fn num_edges(&self) -> usize { self.num_edges }
    fn num_faces(&self) -> usize {
        let indices = self.mesh.indices().expect("Mesh must have indices");
        match indices {
            Indices::U16(v) => v.len() / 3,
            Indices::U32(v) => v.len() / 3,
        }
    }
    fn num_halfedges(&self) -> usize {
        let indices = self.mesh.indices().expect("Mesh must have indices");
        match indices {
            Indices::U16(v) => v.len(),
            Indices::U32(v) => v.len(),
        }
    }
    fn num_corners(&self) -> usize { self.num_halfedges() }
    fn num_boundaries(&self) -> usize { self.boundary_halfedge.len() }
    fn num_generators(&self) -> usize { 0 }

    fn vertex_halfedge(&self, v_idx: usize) -> Option<usize> { self.canonical_vertex_halfedge.get(v_idx).copied().flatten() }
    fn edge_halfedge(&self, e_idx: usize) -> Option<usize> { self.edge_halfedge.get(e_idx).copied().flatten() }
    fn face_halfedge(&self, f_idx: usize) -> Option<usize> { self.face_halfedge.get(f_idx).copied().flatten() }
    fn corner_halfedge(&self, c_idx: usize) -> Option<usize> { Some(c_idx) }
    fn boundary_halfedge(&self, b_idx: usize) -> Option<usize> { self.boundary_halfedge.get(b_idx).copied().flatten() }
    fn generator(&self, _g_idx: usize) -> &[usize] { &[] }

    fn halfedge_vertex(&self, h_idx: usize) -> Option<usize> {
        let indices = self.mesh.indices().expect("Mesh must have indices");
        let raw_idx = match indices {
            Indices::U16(v) => v[h_idx] as usize,
            Indices::U32(v) => v[h_idx] as usize,
        };
        Some(self.raw_to_canonical[raw_idx])
    }
    fn halfedge_edge(&self, h_idx: usize) -> Option<usize> {
        self.halfedge_edge.get(h_idx).copied()
    }
    fn halfedge_face(&self, h_idx: usize) -> Option<usize> { Some(h_idx / 3) }
    fn halfedge_corner(&self, h_idx: usize) -> Option<usize> { Some(h_idx) }
    fn halfedge_next(&self, h_idx: usize) -> Option<usize> {
        Some(if h_idx % 3 == 2 { h_idx - 2 } else { h_idx + 1 })
    }
    fn halfedge_prev(&self, h_idx: usize) -> Option<usize> {
        Some(if h_idx % 3 == 0 { h_idx + 2 } else { h_idx - 1 })
    }
    fn halfedge_twin(&self, h_idx: usize) -> Option<usize> { self.twin.get(h_idx).copied().flatten() }
    fn halfedge_on_boundary(&self, h_idx: usize) -> bool { self.twin[h_idx].is_none() }

    fn add_generator(&mut self, _generator: Vec<usize>) {
    }
}
