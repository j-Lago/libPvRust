use std::fmt;
use crate::pvcell::PvCell;
use crate::series::{Series};

// pub static mut SOLVER_CALLS: usize = 0;

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ParallelSolver {
    pub max_iter: usize,
    pub tol_i: f64,
}
impl Default for ParallelSolver {
    fn default() -> Self {
        ParallelSolver {
            max_iter: 1000,
            tol_i: 0.1,
        }
    }
}

#[derive(Clone)]
#[allow(dead_code)]
pub struct Parallel {
    pub elements: Vec<Series>,
    pub solver: ParallelSolver,
}

impl fmt::Debug for Parallel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // write!(f, "Parallel{:?}", &self.elements)
        write!(f, "Parallel[\n")?;
        for e in self.elements.iter() {
            write!(f, "  {:?}, \n", e)?;
        }
        write!(f, "]")
    }
}

#[allow(dead_code)]
impl Parallel {
    pub fn new(elements: Vec<Series>) -> Self {
        Parallel{ elements, solver: ParallelSolver::default() }
    }

    pub fn empty() -> Self {
        Parallel{ elements: Vec::new(), solver: ParallelSolver::default() }
    }

    pub fn with_solver(settings: ParallelSolver) -> Self {
        println!("{:?}", &settings);
        todo!()
    }
    
    pub fn len(&self) -> usize {
        self.elements.len()
    }

    pub fn find_parallel_equivalent(&self, other: &Series) -> Option<usize> {
        for (k, s) in self.elements.iter().enumerate(){
            if s.is_parallel_equivalent(other) {
                return Some(k);
            }
        }
        
        return None;
    }

    pub fn reduce(&self) -> (Parallel, Vec<u32>, Vec<Vec<u32>>) {
        let mut reduced: Parallel = Parallel::empty();
        let mut origin_to_reduced: Vec<u32> = vec![0; self.len()];
        let mut reduced_to_origin: Vec<Vec<u32>> = vec![];

        for i in 0..self.elements.len() {
            
            let (copy, _, _) = self.elements[i].reduce();
            match reduced.find_parallel_equivalent(&copy) {
                Some(j) => {
                    for (r, c) in reduced.elements[j].elements.iter_mut().zip(copy.iter()) {
                        r.np += c.np;
                    }
                    // reduced_to_origin[j].push(i as u32);
                    // origin_to_reduced[i] = j as u32;
                }
                None => {
                    // origin_to_reduced[i] = reduced.elements.len() as u32;
                    reduced.elements.push(copy.clone());
                    // reduced_to_origin.push(vec![i as u32]);
                }
            }
        }
        return (reduced, origin_to_reduced, reduced_to_origin);
    }

}