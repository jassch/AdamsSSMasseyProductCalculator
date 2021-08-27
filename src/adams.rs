
mod element;
mod generator;
mod multiplication;

pub use element::AdamsElement;
pub use generator::AdamsGenerator;
pub use multiplication::AdamsMultiplication;

/// type synonym for (s,t) bidegrees
pub type Bidegree = (u32, i32);