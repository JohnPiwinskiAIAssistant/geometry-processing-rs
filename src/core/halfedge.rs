#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Halfedge {
    pub vertex: Option<usize>,
    pub edge: Option<usize>,
    pub face: Option<usize>,
    pub corner: Option<usize>,
    pub next: Option<usize>,
    pub prev: Option<usize>,
    pub twin: Option<usize>,
    pub on_boundary: bool,
    pub index: usize,
}

impl Halfedge {
    pub fn new(index: usize) -> Self {
        Self {
            vertex: None,
            edge: None,
            face: None,
            corner: None,
            next: None,
            prev: None,
            twin: None,
            on_boundary: false,
            index,
        }
    }
}
