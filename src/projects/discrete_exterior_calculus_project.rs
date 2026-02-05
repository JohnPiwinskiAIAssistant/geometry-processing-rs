use crate::core::geometry::Geometry;
use crate::core::dec::DEC;
use crate::linear_algebra::SparseMatrix;

pub struct DiscreteExteriorCalculusProject<'a> {
    pub geometry: &'a Geometry<'a>,
}

impl<'a> DiscreteExteriorCalculusProject<'a> {
    pub fn new(geometry: &'a Geometry<'a>) -> Self {
        Self { geometry }
    }

    pub fn build_hodge_star_0_form(&self) -> SparseMatrix {
        DEC::build_hodge_star_0_form(self.geometry)
    }

    pub fn build_hodge_star_1_form(&self) -> SparseMatrix {
        DEC::build_hodge_star_1_form(self.geometry)
    }

    pub fn build_hodge_star_2_form(&self) -> SparseMatrix {
        DEC::build_hodge_star_2_form(self.geometry)
    }

    pub fn build_exterior_derivative_0_form(&self) -> SparseMatrix {
        DEC::build_exterior_derivative_0_form(self.geometry)
    }

    pub fn build_exterior_derivative_1_form(&self) -> SparseMatrix {
        DEC::build_exterior_derivative_1_form(self.geometry)
    }
}
