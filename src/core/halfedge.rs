#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Halfedge {
    pub index: usize,
}

impl Halfedge {
    pub fn new(index: usize) -> Self {
        Self { index }
    }
}
