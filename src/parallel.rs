use crate::series::{Series};

pub static mut SOLVER_CALLS: usize = 0;

#[derive(Debug, Clone)]
struct ParallelSolver {
    max_iter: usize,
    tol_i: f64,
}
impl Default for ParallelSolver {
    fn default() -> Self {
        ParallelSolver {
            max_iter: 1000,
            tol_i: 0.1,
        }
    }
}

#[derive(Debug, Clone)]
struct Parallel {
    elements: Vec<Series>,
    solver: ParallelSolver,
}

#[allow(dead_code)]
impl Parallel {
    fn new(elements: Vec<Series>) -> Self {
        Parallel{ elements, solver: ParallelSolver::default() }
    }

    fn empty() -> Self {
        Parallel{ elements: Vec::new(), solver: ParallelSolver::default() }
    }

    fn with_solver(settings: ParallelSolver) -> Self {
        todo!()
    }

}