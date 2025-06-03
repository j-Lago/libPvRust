
use crate::pvcell::{CellState, PvCell};
use std::ops::Index;
use std::iter::IntoIterator;
use tracing::{warn, error};
use std::fmt;

pub static mut SOLVER_CALLS: usize = 0;

#[derive(Debug, Clone)]
pub struct SeriesSolver {
    pub max_iter: usize,
    pub tol_v: f64,
    pub min_g: f64,
}

impl Default for SeriesSolver {
    fn default() -> Self {
        SeriesSolver {
            max_iter: 1000,
            tol_v: 0.1,
            min_g: 0.00001,
        }
    }
}

#[derive(Clone)]
pub struct Series {
    elements: Vec<PvCell>,
    solver: SeriesSolver,
}

#[allow(dead_code)]
impl Series {
    pub fn empty() -> Series {
        return Series {
            elements: Vec::new(),
            solver: SeriesSolver::default()
        };
    }

    pub fn new(panels: Vec<PvCell>) -> Series {
        return Series {
            elements: panels,
            solver: SeriesSolver::default()
        };
    }

    /// builder: let s: Series = Series::new(...).with_solver(...);
    pub fn with_solver(mut self, settings: SeriesSolver) -> Series { self.solver = settings; return self; }

    pub fn len(&self) -> usize {
        return self.elements.len();
    }

    pub fn iter(&self) -> impl Iterator<Item = &PvCell> {
        self.elements.iter()
    }

    pub fn push(&mut self, element: PvCell){
        self.elements.push(element);
    }

    pub fn states_uniform_conditions(&self, irrad_ef: f64, cell_temp: f64) -> Vec<CellState> {
        let mut states: Vec<CellState> = Vec::with_capacity(self.len());
        for pnl in self.elements.iter(){
            states.push(pnl.compute_state(irrad_ef, cell_temp));
        }
        return states;
    }

    pub fn vs_from_i(&self, states: &Vec<CellState>, i: f64) -> Vec<f64> {
        let mut voltages: Vec<f64> = Vec::with_capacity(self.len());
        for (k, pnl) in self.iter().enumerate(){
            voltages.push(pnl.v_from_i(&states[k], i))
        }
        voltages
    }

    pub fn v_from_i(&self, states: &Vec<CellState>, i: f64) -> f64 {
        let voltages = self.vs_from_i(states, i);
        voltages.iter().sum()
    }

    pub fn i_from_v(&self, states: &Vec<CellState>, v_str: f64) -> f64 {
        unsafe{ SOLVER_CALLS += 1; }

        let mut sum_voc: f64 = 0.0;
        let mut il: f64 = f64::INFINITY;
        let mut i0: f64 = f64::INFINITY;
        for (k, pnl) in self.elements.iter().enumerate() {
            sum_voc += pnl.v_oc_ref * (pnl.ns as f64);
            il = il.min(pnl.i_l_ref * (pnl.np as f64));
            i0 = i0.min(pnl.solve_i(&states[k], 0.0));
        }
        let mut g: f64 = il / sum_voc;
        // let mut g: f64 = 10.0 / sum_voc;

        let mut dv1: f64 = 0.0;
        for _ in 0..self.solver.max_iter {
            let v0 = self.v_from_i(states, i0);
            let dv = v0 - v_str;
            if dv.abs() < self.solver.tol_v {
                return i0;
            }
            if dv * dv1 < 0.0 {  // mudança de sinal
                g /= 2.0;
            }
            g = g.max(self.solver.min_g);
            i0 += dv * g;
            dv1 = dv;
        }

        if i0.is_normal() {
            warn!("({:p}) Series::i_from_v(v_str={:e}) nao convergiu (tol={:e}, max_iter={}) -> (i={})",
                &self, v_str, self.solver.tol_v, self.solver.max_iter, i0);
        } else {
            error!("({:p}) Series::i_from_v(v_str={:e}) nao convergiu (tol={:e}, max_iter={}) -> (i={})",
                &self, v_str, self.solver.tol_v, self.solver.max_iter, i0);
        }
        return i0;
    }

    pub fn find(&self, other: &PvCell) -> Option<usize> {
        for (k, pnl) in self.elements.iter().enumerate(){
            if pnl.is_series_equivalent(other) {
                return Some(k);
            }
        }
        return None;
    }

    pub fn reduce(&self) -> (Series, Vec<u32>, Vec<Vec<u32>>) {
        let mut reduced: Series = Series::empty();
        let mut origin_to_reduced: Vec<u32> = vec![0; self.len()];
        let mut reduced_to_origin: Vec<Vec<u32>> = vec![];

        for i in 0..self.elements.len() {
            match reduced.find(&self.elements[i]) {
                Some(j) => {
                    reduced.elements[j].ns += self.elements[i].ns;
                    reduced_to_origin[j].push(i as u32);
                    origin_to_reduced[i] = j as u32;
                }
                None => {
                    origin_to_reduced[i] = reduced.elements.len() as u32;
                    reduced.elements.push(self.elements[i].clone());
                    reduced_to_origin.push(vec![i as u32]);
                }
            }
        }
        return (reduced, origin_to_reduced, reduced_to_origin);
    }
}

impl fmt::Debug for Series {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "PvString {:?}", &self.elements)
    }
}

impl Index<usize> for Series {
    type Output = PvCell;
    fn index(&self, index: usize) -> &Self::Output {
        &self.elements[index]
    }
}

impl IntoIterator for Series {
    type Item = PvCell;
    type IntoIter = std::vec::IntoIter<PvCell>;
    fn into_iter(self) -> Self::IntoIter {
        return self.elements.into_iter();
    }
}

// impl StringStates{
//     #[allow(dead_code)]
//     pub fn iter(&self) -> impl Iterator<Item = &CellState> {
//         self.elements.iter()
//     }
//
//     #[allow(dead_code)]
//     fn len(&self) -> usize {
//         self.elements.len()
//     }
// }

// impl IntoIterator for StringStates {
//     type Item = CellState;
//     type IntoIter = std::vec::IntoIter<CellState>;
//     fn into_iter(self) -> Self::IntoIter {
//         return self.elements.into_iter();
//     }
// }

// impl Index<usize> for StringStates {
//     type Output = CellState;
//     fn index(&self, i: usize) -> &Self::Output {
//         &self.elements[i]
//     }
// }
