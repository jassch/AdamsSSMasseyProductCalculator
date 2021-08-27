
mod element;
mod generator;

pub use element::AdamsElement;
pub use generator::AdamsGenerator;

/// type synonym for (s,t) bidegrees
pub type Bidegree = (u32, i32);