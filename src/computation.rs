
use std::fmt;
use std::fmt::{Display, Formatter};


#[derive(Copy, Clone, PartialEq, PartialOrd, Eq, Ord, Debug)]
pub enum ComputationResult<T, E, PartialT=T> {
    FullyComputed(T),
    PartiallyComputed(PartialT),
    Uncomputed,
    Error(E),
}

impl <T, E, PartialT> ComputationResult<T, E, PartialT> {
    /// converter function takes self by move
    pub fn full_value(self) -> Option<T> {
        match self {
            ComputationResult::FullyComputed(t) => Some(t),
            _ => None
        }
    }
    /// get value by reference
    pub fn full_value_ref(&self) -> Option<&T> {
        match self {
            Self::FullyComputed(t) => Some(t),
            _ => None
        }
    }
    /// converter function takes self by move
    pub fn fully_computed(&self) -> bool {
        match self {
            Self::FullyComputed(_) => true,
            _ => false
        }
    }
    /// convert object to result, consuming the object
    pub fn to_result<F>(self, uncomputed: E, handle_partial: F) -> Result<T,E> where
        F: Fn(PartialT) -> Result<T,E> {
        match self {
            Self::FullyComputed(t) => Ok(t),
            Self::PartiallyComputed(p) => handle_partial(p),
            Self::Uncomputed => Err(uncomputed),
            Self::Error(e) => Err(e)
        }
    }
}

impl <T,E> ComputationResult<T,E> {
    pub fn to_result_or_partial(self, uncomputed: E) -> Result<T,E> {
        match self {
            Self::FullyComputed(t) => Ok(t),
            Self::PartiallyComputed(p) => Ok(p),
            Self::Uncomputed => Err(uncomputed),
            Self::Error(e) => Err(e)
        }
    }
}

impl <T,E,PartialT> Display for ComputationResult<T,E,PartialT> where
    T: Display, E: Display, PartialT: Display 
{
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        match self {
            Self::FullyComputed(t) => write!(f, "full({})", t),
            Self::PartiallyComputed(p) => write!(f, "partial({})", p),
            Self::Uncomputed => write!(f, "uncomputed"),
            Self::Error(e) => write!(f, "error({})", e),
        }
    }
}